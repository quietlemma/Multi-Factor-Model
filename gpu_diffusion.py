import time
from typing import Callable, Dict, Iterable, Optional, Sequence, Tuple

import numpy as np


OriginNode = Tuple[int, str]


def cuda_available() -> bool:
    try:
        import torch

        return torch.cuda.is_available()
    except Exception:
        return False


def _node_sort_key(node):
    return node[0], node[1]


def _safe_progress(progress: Optional[Callable[[str], None]], message: str):
    if progress:
        progress(message)


class GPURealizationCollection:
    """GPU-backed realization collection with the same influence semantics as CPU BFS."""

    def __init__(
        self,
        generated_graph,
        realization_collection: Sequence,
        device: str = "cuda",
        progress: Optional[Callable[[str], None]] = None,
    ):
        import torch

        if device.startswith("cuda") and not torch.cuda.is_available():
            raise RuntimeError("CUDA is not available, cannot build GPU realization collection")

        self.torch = torch
        self.device = torch.device(device)
        self.realization_count = len(realization_collection)
        self.nodes = sorted(generated_graph.nodes, key=_node_sort_key)
        self.node_to_idx: Dict[OriginNode, int] = {node: idx for idx, node in enumerate(self.nodes)}
        self.node_count = len(self.nodes)
        self.total_node_count = self.realization_count * self.node_count
        self.origin_indices = torch.tensor(
            [idx for idx, node in enumerate(self.nodes) if node[1] == "origin"],
            dtype=torch.long,
            device=self.device,
        )
        self.realization_offsets = (
            torch.arange(self.realization_count, dtype=torch.long, device=self.device) * self.node_count
        )

        _safe_progress(
            progress,
            "GPU evaluator: indexing realization edges "
            f"(R={self.realization_count}, nodes={self.node_count})",
        )
        started_at = time.time()
        edge_counts = np.fromiter(
            (sum(len(children) for children in realization.children.values()) for realization in realization_collection),
            dtype=np.int64,
            count=self.realization_count,
        )
        total_edges = int(edge_counts.sum())
        rows = np.empty(total_edges, dtype=np.int64)
        cols = np.empty(total_edges, dtype=np.int64)

        cursor = 0
        progress_interval = max(1, self.realization_count // 10)
        for ridx, realization in enumerate(realization_collection):
            base = ridx * self.node_count
            for src, children in realization.children.items():
                src_idx = self.node_to_idx[src]
                for dst in children:
                    rows[cursor] = base + self.node_to_idx[dst]
                    cols[cursor] = base + src_idx
                    cursor += 1
            if ridx == 0 or ridx + 1 == self.realization_count or (ridx + 1) % progress_interval == 0:
                _safe_progress(
                    progress,
                    f"GPU evaluator: indexed {ridx + 1}/{self.realization_count} realizations",
                )

        if cursor != total_edges:
            rows = rows[:cursor]
            cols = cols[:cursor]
            total_edges = cursor

        _safe_progress(
            progress,
            f"GPU evaluator: moving sparse adjacency to {self.device} with {total_edges} live edges",
        )
        indices = torch.from_numpy(np.vstack((rows, cols))).to(self.device)
        values = torch.ones(total_edges, dtype=torch.float32, device=self.device)
        self.adj_t = torch.sparse_coo_tensor(
            indices,
            values,
            size=(self.total_node_count, self.total_node_count),
            device=self.device,
        ).coalesce()

        del rows, cols, indices, values
        if self.device.type == "cuda":
            torch.cuda.empty_cache()

        elapsed = time.time() - started_at
        _safe_progress(progress, f"GPU evaluator: ready in {elapsed:.1f}s")

    def __len__(self):
        return self.realization_count

    def _seed_flat_indices(self, seeds: Iterable[OriginNode]):
        seed_node_indices = [
            self.node_to_idx[seed]
            for seed in seeds
            if seed in self.node_to_idx
        ]
        if not seed_node_indices:
            return None
        seed_node_indices = self.torch.tensor(seed_node_indices, dtype=self.torch.long, device=self.device)
        return (self.realization_offsets[:, None] + seed_node_indices[None, :]).reshape(-1)

    def compute_influence(self, seeds) -> float:
        torch = self.torch
        seed_flat_indices = self._seed_flat_indices(seeds)
        if seed_flat_indices is None:
            return 0.0

        visited = torch.zeros(self.total_node_count, dtype=torch.bool, device=self.device)
        visited[seed_flat_indices] = True
        frontier = torch.zeros((self.total_node_count, 1), dtype=torch.float32, device=self.device)
        frontier[seed_flat_indices, 0] = 1.0

        for _ in range(self.node_count):
            reached = torch.sparse.mm(self.adj_t, frontier).squeeze(1) > 0
            new_nodes = reached & ~visited
            if not bool(new_nodes.any().item()):
                break
            visited |= new_nodes
            frontier = new_nodes.to(dtype=torch.float32).unsqueeze(1)

        origin_counts = visited.view(self.realization_count, self.node_count).index_select(
            1,
            self.origin_indices,
        ).sum(dim=1)
        result = origin_counts.to(dtype=torch.float32).mean().item()

        del visited, frontier, origin_counts
        return result
