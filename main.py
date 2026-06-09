import argparse
import contextlib
import csv
import io
import random
import re
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

from LFG import algo
import diffusionModel
import gpu_diffusion
from RRSet import mf_rr
import numpy as np
import tools


BASE_DIR = Path(__file__).resolve().parent


DATASETS: Dict[str, Dict[str, object]] = {
    "NetScience": {
        "path": BASE_DIR / "LFG" / "ca-netscience.txt",
        "reader": tools.readGraph_undirect,
    },
    "Wiki": {
        "path": BASE_DIR / "LFG" / "soc-wiki-Vote.txt",
        "reader": tools.readGraph_direct,
    },
    "Gnutella": {
        "path": BASE_DIR / "LFG" / "p2p-Gnutella08.txt",
        "reader": tools.readGraph_direct,
    },
}


GRAPH_BUILDERS: Dict[str, Dict[int, Callable]] = {
    "MFIC": {
        2: tools.generateGraph_factor_IC_2,
        3: tools.generateGraph_factor_IC_3,
        5: tools.generateGraph_factor_IC_5,
    },
    "MFLT": {
        2: tools.generateGraph_factor_LT_2,
        3: tools.generateGraph_factor_LT_3,
        5: tools.generateGraph_factor_LT_5,
    },
    "MFTRIGGER": {
        2: tools.generateGraph_factor_TRIGGER_2,
        3: tools.generateGraph_factor_TRIGGER_3,
        5: tools.generateGraph_factor_TRIGGER_5,
    },
}


REALIZATION_BUILDERS: Dict[str, Callable] = {
    "MFIC": tools.generateCollectionOfRealizationIC,
    "MFLT": tools.generateCollectionOfRealizationLT,
    "MFTRIGGER": tools.generateCollectionOfRealizationTrigger,
}


PROGRESS_PATTERN = re.compile(r"k =\s*(\d+)\s*,?\s*Influence =\s*([0-9.]+)")


def _progress(message: str, enabled: bool = True):
    if not enabled:
        return
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", file=sys.stderr, flush=True)


class _TeeCapture(io.StringIO):
    """Capture legacy algorithm output while still showing it in the log."""

    def __init__(self, stream):
        super().__init__()
        self._stream = stream

    def write(self, text):
        self._stream.write(text)
        self._stream.flush()
        return super().write(text)

    def flush(self):
        self._stream.flush()
        return super().flush()


def _safe_size(value) -> str:
    try:
        return str(len(value))
    except TypeError:
        return "unknown"


def _normalize_algorithm_name(name: str) -> str:
    return name.strip().lower()


def _algorithm_order(algorithms: Sequence[str]) -> List[str]:
    return [_normalize_algorithm_name(name) for name in algorithms]


def _run_baseline(name: str, graph, generated_graph, realization_collection, k: int, theta: int, seed: int):
    if name == "lfg":
        return algo.LFGreedy(graph, realization_collection, k)
    if name == "greedy":
        return algo.greedy(graph, realization_collection, k)
    if name == "maxdegree":
        return algo.maxDegree(graph, realization_collection, k)
    if name == "random":
        random.seed(seed)
        return algo.randomm(graph, realization_collection, k)
    if name == "gwdm":
        return algo.greedyWDM(graph, realization_collection, k)
    raise ValueError(f"unsupported algorithm: {name}")


def _run_incremental_baseline(
    name: str,
    graph,
    generated_graph,
    realization_collection,
    max_k: int,
    theta: int,
    seed: int,
    progress: bool = True,
):
    if max_k <= 0:
        return {}

    start_time = time.time()
    _progress(f"Algorithm {_canonical_label(name)} started; target max_k={max_k}", progress)
    buffer = _TeeCapture(sys.stderr) if progress else io.StringIO()
    with contextlib.redirect_stdout(buffer):
        _run_baseline(name, graph, generated_graph, realization_collection, max_k, theta, seed)

    results = {}
    for line in buffer.getvalue().splitlines():
        match = PROGRESS_PATTERN.search(line)
        if match:
            k = int(match.group(1))
            influence = float(match.group(2))
            results[k] = influence
    elapsed = time.time() - start_time
    _progress(
        f"Algorithm {_canonical_label(name)} finished; "
        f"captured {len(results)} budget points in {elapsed:.1f}s",
        progress,
    )
    return results


def _canonical_label(name: str) -> str:
    labels = {
        "lfg": "LFG",
        "greedy": "Greedy",
        "maxdegree": "MaxDegree",
        "random": "Random",
        "gwdm": "GWDM",
        "mf-rr": "MF-RR",
        "mf_rr": "MF-RR",
        "mfrr": "MF-RR",
        "mf-imm": "MF-IMM",
        "mf_imm": "MF-IMM",
        "mfimm": "MF-IMM",
    }
    return labels[_normalize_algorithm_name(name)]


def _rr_model_name(model: str) -> str:
    if model == "MFIC":
        return "IC"
    if model == "MFLT":
        return "LT"
    if model == "MFTRIGGER":
        return "TRIGGER"
    raise ValueError(f"unsupported model: {model}")


def run_experiment(
    dataset: str,
    model: str,
    factor_count: int,
    budgets: Sequence[int],
    theta: int,
    realization_count: int,
    seed: int,
    algorithms: Sequence[str],
    imm_epsilon: float = 0.5,
    imm_ell: float = 1.0,
    imm_max_theta: int = 0,
    greedy_from_lfg: bool = False,
    progress: bool = True,
    checkpoint_path: Optional[Path] = None,
    use_gpu: bool = True,
    gpu_device: str = "cuda",
    require_gpu: bool = False,
):
    if dataset not in DATASETS:
        raise ValueError(f"unsupported dataset: {dataset}")
    if model not in GRAPH_BUILDERS:
        raise ValueError(f"unsupported model: {model}")
    if factor_count not in GRAPH_BUILDERS[model]:
        raise ValueError(f"unsupported factor count for {model}: {factor_count}")

    _progress(
        f"Experiment started: dataset={dataset}, model={model}, factors={factor_count}, "
        f"budgets={list(budgets)}, theta={theta}, realizations={realization_count}, seed={seed}, "
        f"use_gpu={use_gpu}, gpu_device={gpu_device}, "
        f"algorithms={list(algorithms)}",
        progress,
    )

    dataset_info = DATASETS[dataset]
    _progress(f"Loading graph from {dataset_info['path']}", progress)
    graph = dataset_info["reader"](dataset_info["path"])
    _progress(
        f"Original graph loaded: nodes={_safe_size(graph.nodes)}, edges={_safe_size(graph.edges)}",
        progress,
    )
    random.seed(seed)
    np.random.seed(seed)
    _progress("Generating auxiliary graph G'", progress)
    generated_graph = GRAPH_BUILDERS[model][factor_count](graph)
    _progress(
        f"Auxiliary graph generated: nodes={_safe_size(generated_graph.nodes)}, "
        f"edges={_safe_size(generated_graph.edges)}",
        progress,
    )
    random.seed(seed + 1)
    np.random.seed(seed + 1)
    _progress(f"Generating {realization_count} realizations", progress)
    realization_collection = REALIZATION_BUILDERS[model](generated_graph, number=realization_count)
    _progress(f"Realization collection generated: count={len(realization_collection)}", progress)
    if use_gpu and model == "MFLT":
        if gpu_diffusion.cuda_available() or not gpu_device.startswith("cuda"):
            _progress(f"Converting MFLT realizations to GPU evaluator on {gpu_device}", progress)
            realization_collection = gpu_diffusion.GPURealizationCollection(
                generated_graph,
                realization_collection,
                device=gpu_device,
                progress=lambda msg: _progress(msg, progress),
            )
            _progress("GPU evaluator enabled for diffusionModel.compute(...)", progress)
        elif require_gpu:
            raise RuntimeError("GPU was required, but CUDA is not available")
        else:
            _progress("CUDA is not available; falling back to CPU evaluator", progress)
    elif use_gpu and model != "MFLT":
        _progress("GPU evaluator is currently enabled only for MFLT; using CPU evaluator for this model", progress)

    normalized_algorithms = _algorithm_order(algorithms)
    if greedy_from_lfg and "greedy" in normalized_algorithms and "lfg" not in normalized_algorithms:
        raise ValueError("--greedy-from-lfg requires LFG to be included in --algorithms")
    max_budget = max(budgets) if budgets else 0
    baseline_progress = {}
    for algorithm in normalized_algorithms:
        if algorithm in {"mf-rr", "mf_rr", "mfrr", "mf-imm", "mf_imm", "mfimm"}:
            continue
        if algorithm == "greedy" and greedy_from_lfg:
            continue
        baseline_progress[algorithm] = _run_incremental_baseline(
            algorithm,
            graph,
            generated_graph,
            realization_collection,
            max_budget,
            theta,
            seed,
            progress=progress,
        )

    rows = []
    completed_algorithm_labels = set()
    for k in budgets:
        _progress(f"Budget k={k} started", progress)
        row = {"size": k}
        if k == 0:
            for algorithm in normalized_algorithms:
                row[_canonical_label(algorithm)] = 0.0
            rows.append(row)
            if checkpoint_path:
                _write_csv(rows, checkpoint_path)
                _progress(f"Checkpoint saved to {checkpoint_path}", progress)
            continue

        for algorithm in normalized_algorithms:
            label = _canonical_label(algorithm)
            if algorithm in {"mf-rr", "mf_rr", "mfrr"}:
                rr_model = _rr_model_name(model)
                _progress(f"{label} k={k}: seed selection started", progress)
                seeds = mf_rr.mf_rr_greedy(
                    graph,
                    generated_graph,
                    rr_model,
                    k,
                    theta,
                    seed=seed,
                    progress=lambda msg: _progress(f"{label} k={k}: {msg}", progress),
                )
                _progress(f"{label} k={k}: evaluating {len(seeds)} seeds", progress)
                influence = diffusionModel.compute(realization_collection, seeds)
            elif algorithm in {"mf-imm", "mf_imm", "mfimm"}:
                rr_model = _rr_model_name(model)
                _progress(f"{label} k={k}: seed selection started", progress)
                seeds = mf_rr.mf_imm_greedy(
                    graph,
                    generated_graph,
                    rr_model,
                    k,
                    epsilon=imm_epsilon,
                    ell=imm_ell,
                    seed=seed,
                    max_theta=imm_max_theta or None,
                    progress=lambda msg: _progress(f"{label} k={k}: {msg}", progress),
                )
                _progress(f"{label} k={k}: evaluating {len(seeds)} seeds", progress)
                influence = diffusionModel.compute(realization_collection, seeds)
            elif algorithm == "greedy" and greedy_from_lfg:
                influence = baseline_progress["lfg"][k]
            else:
                influence = baseline_progress[algorithm][k]
            row[label] = influence
            _progress(f"{label} k={k}: influence={influence}", progress)
            if k == max_budget and label not in completed_algorithm_labels:
                completed_algorithm_labels.add(label)
                _progress(f"Algorithm {label} completed through k={max_budget}", progress)
        rows.append(row)
        _progress(f"Budget k={k} completed: {row}", progress)
        if checkpoint_path:
            _write_csv(rows, checkpoint_path)
            _progress(f"Checkpoint saved to {checkpoint_path}", progress)
    _progress("Experiment finished", progress)
    return rows


def _write_csv(rows, path: Path):
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Run multi-factor experiments with RR-set baselines.")
    parser.add_argument("--dataset", choices=sorted(DATASETS.keys()), default="NetScience")
    parser.add_argument("--model", choices=sorted(GRAPH_BUILDERS.keys()), default="MFIC")
    parser.add_argument("--factors", type=int, choices=[2, 3, 5], default=2)
    parser.add_argument("--theta", type=int, default=1000)
    parser.add_argument("--imm-epsilon", type=float, default=0.5)
    parser.add_argument("--imm-ell", type=float, default=1.0)
    parser.add_argument("--imm-max-theta", type=int, default=0)
    parser.add_argument("--realizations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-k", type=int, default=50)
    parser.add_argument("--step", type=int, default=5)
    parser.add_argument(
        "--algorithms",
        nargs="+",
        default=["LFG", "Greedy", "MaxDegree", "Random", "GWDM", "MF-RR"],
    )
    parser.add_argument(
        "--greedy-from-lfg",
        action="store_true",
        help="Fill the Greedy column with LFG values instead of running the slow Greedy baseline.",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--quiet-progress",
        action="store_true",
        help="Disable progress messages. By default progress is printed to stderr.",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Disable the GPU evaluator and use the original CPU evaluator.",
    )
    parser.add_argument("--gpu-device", default="cuda")
    parser.add_argument(
        "--require-gpu",
        action="store_true",
        help="Fail immediately if --no-gpu is not set and CUDA is unavailable.",
    )
    args = parser.parse_args()

    budgets = list(range(0, args.max_k + 1, args.step))
    rows = run_experiment(
        dataset=args.dataset,
        model=args.model,
        factor_count=args.factors,
        budgets=budgets,
        theta=args.theta,
        realization_count=args.realizations,
        seed=args.seed,
        algorithms=args.algorithms,
        imm_epsilon=args.imm_epsilon,
        imm_ell=args.imm_ell,
        imm_max_theta=args.imm_max_theta,
        greedy_from_lfg=args.greedy_from_lfg,
        progress=not args.quiet_progress,
        checkpoint_path=args.output,
        use_gpu=not args.no_gpu,
        gpu_device=args.gpu_device,
        require_gpu=args.require_gpu,
    )

    if args.output:
        _write_csv(rows, args.output)
    else:
        writer = csv.DictWriter(
            __import__("sys").stdout,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
