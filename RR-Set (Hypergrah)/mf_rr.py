import math
import random
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple


OriginNode = Tuple[int, str]
RRSet = Set[Tuple[int, str]]


def _origin_nodes(original_graph) -> List[OriginNode]:
    return sorted((node, "origin") for node in original_graph.nodes)


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _sample_rr_ic(generated_graph, root: OriginNode, rng: random.Random) -> RRSet:
    rr_set: RRSet = {root}
    queue: List[OriginNode] = [root]
    while queue:
        current = queue.pop(0)
        for parent in generated_graph.get_parentss(current):
            if parent in rr_set:
                continue
            if rng.random() <= generated_graph.edges[(parent, current)]:
                rr_set.add(parent)
                queue.append(parent)
    return rr_set


def _weighted_parent_choice(generated_graph, node, rng: random.Random):
    parents = list(generated_graph.get_parentss(node))
    if not parents:
        return None
    weights = [generated_graph.edges[(parent, node)] for parent in parents]
    total = sum(weights)
    if total <= 0:
        return None
    target = rng.random() * total
    cumulative = 0.0
    for parent, weight in zip(parents, weights):
        cumulative += weight
        if target <= cumulative:
            return parent
    return parents[-1]


def _sample_rr_lt(generated_graph, root: OriginNode, rng: random.Random) -> RRSet:
    rr_set: RRSet = {root}
    queue: List[OriginNode] = [root]
    while queue:
        current = queue.pop(0)
        parent = _weighted_parent_choice(generated_graph, current, rng)
        if parent is None or parent in rr_set:
            continue
        rr_set.add(parent)
        queue.append(parent)
    return rr_set


def _trigger_probability(generated_graph, node) -> float:
    probabilities = getattr(generated_graph, "trigger_probabilities", None)
    if probabilities is not None and node in probabilities:
        return probabilities[node]

    parents = list(generated_graph.get_parentss(node))
    if not parents:
        return 0.0
    return max(generated_graph.edges.get((parent, node), 0.0) for parent in parents)


def _sample_rr_trigger(generated_graph, root: OriginNode, rng: random.Random) -> RRSet:
    rr_set: RRSet = {root}
    queue: List[OriginNode] = [root]
    while queue:
        current = queue.pop(0)
        parents = list(generated_graph.get_parentss(current))
        if not parents:
            continue
        if rng.random() > _trigger_probability(generated_graph, current):
            continue
        for parent in parents:
            if parent in rr_set:
                continue
            rr_set.add(parent)
            queue.append(parent)
    return rr_set


def _sampler_for_model(model_name: str):
    if model_name == "IC":
        return _sample_rr_ic
    if model_name == "LT":
        return _sample_rr_lt
    if model_name == "TRIGGER":
        return _sample_rr_trigger
    raise ValueError("model must be 'IC', 'LT', or 'TRIGGER'")


def _sample_rr_sets_with_rng(
    origin_nodes: Sequence[OriginNode],
    generated_graph,
    model_name: str,
    theta: int,
    rng: random.Random,
    progress: Optional[Callable[[str], None]] = None,
    progress_label: str = "RR sets",
) -> List[RRSet]:
    sampler = _sampler_for_model(model_name)
    rr_sets: List[RRSet] = []
    progress_interval = max(1, theta // 10) if progress and theta > 0 else 0
    for idx in range(theta):
        root = rng.choice(origin_nodes)
        rr_sets.append(sampler(generated_graph, root, rng))
        sampled = idx + 1
        if progress and (sampled == 1 or sampled == theta or sampled % progress_interval == 0):
            progress(f"{progress_label}: sampled {sampled}/{theta} RR sets")
    return rr_sets


def sample_rr_sets(original_graph, generated_graph, model: str, theta: int, seed: Optional[int] = None) -> List[RRSet]:
    if theta < 0:
        raise ValueError("theta must be non-negative")
    model_name = model.upper()
    _sampler_for_model(model_name)

    origin_nodes = _origin_nodes(original_graph)
    if not origin_nodes or theta == 0:
        return []

    rng = random.Random(seed)
    return _sample_rr_sets_with_rng(origin_nodes, generated_graph, model_name, theta, rng)


def _build_origin_cover_index(origin_nodes: Sequence[OriginNode], rr_sets: Sequence[RRSet]) -> Dict[OriginNode, Set[int]]:
    cover_index: Dict[OriginNode, Set[int]] = {node: set() for node in origin_nodes}
    for idx, rr_set in enumerate(rr_sets):
        for node in rr_set:
            if node in cover_index:
                cover_index[node].add(idx)
    return cover_index


def _max_cover_greedy(origin_nodes: Sequence[OriginNode], rr_sets: Sequence[RRSet], k: int) -> Tuple[Set[OriginNode], int]:
    cover_index = _build_origin_cover_index(origin_nodes, rr_sets)
    covered_rr: Set[int] = set()
    selected: Set[OriginNode] = set()
    target_k = min(k, len(origin_nodes))

    for _ in range(target_k):
        best_node = None
        best_gain = -1
        for node in origin_nodes:
            if node in selected:
                continue
            gain = len(cover_index[node] - covered_rr)
            if gain > best_gain:
                best_gain = gain
                best_node = node
        if best_node is None:
            break
        selected.add(best_node)
        covered_rr.update(cover_index[best_node])
    return selected, len(covered_rr)


def mf_rr_greedy(
    original_graph,
    generated_graph,
    model: str,
    k: int,
    theta: int,
    seed: Optional[int] = None,
    progress: Optional[Callable[[str], None]] = None,
) -> Set[OriginNode]:
    if k <= 0:
        return set()

    origin_nodes = _origin_nodes(original_graph)
    if not origin_nodes:
        return set()

    model_name = model.upper()
    _sampler_for_model(model_name)

    rng = random.Random(seed)
    rr_sets = _sample_rr_sets_with_rng(
        origin_nodes,
        generated_graph,
        model_name,
        theta,
        rng,
        progress=progress,
        progress_label="MF-RR",
    )
    if not rr_sets:
        return set()

    if progress:
        progress(f"MF-RR: running max-cover greedy with k={k}")
    selected, _ = _max_cover_greedy(origin_nodes, rr_sets, k)
    if progress:
        progress(f"MF-RR: selected {len(selected)} seeds")
    return selected


def _estimate_imm_lower_bound(
    origin_nodes: Sequence[OriginNode],
    generated_graph,
    model_name: str,
    k: int,
    epsilon: float,
    ell: float,
    rng: random.Random,
    progress: Optional[Callable[[str], None]] = None,
) -> float:
    n = len(origin_nodes)
    if n <= 1:
        return 1.0

    epsilon_prime = math.sqrt(2.0) * epsilon
    log_n = math.log(n)
    log_log_n = math.log(max(math.log(n, 2), 1.000000001))
    log_choose = _log_comb(n, min(k, n))
    lambda_prime = (2.0 + 2.0 * epsilon_prime / 3.0) * (
        log_choose + ell * log_n + log_log_n
    ) * n / (epsilon_prime ** 2)

    rr_sets: List[RRSet] = []
    lower_bound = 1.0
    max_round = max(1, math.ceil(math.log(n, 2)))
    for i in range(1, max_round):
        x = n / (2 ** i)
        theta_i = math.ceil(lambda_prime / x)
        needed = theta_i - len(rr_sets)
        if progress:
            progress(
                f"MF-IMM lower bound round {i}/{max_round - 1}: "
                f"x={x:.4f}, target_theta={theta_i}, extra={max(needed, 0)}"
            )
        if needed > 0:
            rr_sets.extend(
                _sample_rr_sets_with_rng(
                    origin_nodes,
                    generated_graph,
                    model_name,
                    needed,
                    rng,
                    progress=progress,
                    progress_label=f"MF-IMM lower bound round {i}",
                )
            )

        _, covered_count = _max_cover_greedy(origin_nodes, rr_sets, k)
        estimated_spread = n * covered_count / len(rr_sets)
        if progress:
            progress(
                f"MF-IMM lower bound round {i}: "
                f"covered={covered_count}/{len(rr_sets)}, estimate={estimated_spread:.4f}"
            )
        if estimated_spread >= (1.0 + epsilon_prime) * x:
            lower_bound = estimated_spread / (1.0 + epsilon_prime)
            if progress:
                progress(f"MF-IMM lower bound accepted: {lower_bound:.4f}")
            break
    return max(lower_bound, 1.0)


def _imm_theta(n: int, k: int, epsilon: float, ell: float, lower_bound: float) -> int:
    log_n = math.log(max(n, 2))
    log_choose = _log_comb(n, min(k, n))
    alpha = math.sqrt(ell * log_n + math.log(2.0))
    beta = math.sqrt((1.0 - 1.0 / math.e) * (log_choose + ell * log_n + math.log(2.0)))
    lambda_star = 2.0 * n * (((1.0 - 1.0 / math.e) * alpha + beta) ** 2) / (epsilon ** 2)
    return max(1, math.ceil(lambda_star / max(lower_bound, 1.0)))


def mf_imm_greedy(
    original_graph,
    generated_graph,
    model: str,
    k: int,
    epsilon: float = 0.5,
    ell: float = 1.0,
    seed: Optional[int] = None,
    max_theta: Optional[int] = None,
    progress: Optional[Callable[[str], None]] = None,
) -> Set[OriginNode]:
    if k <= 0:
        return set()
    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    if ell <= 0:
        raise ValueError("ell must be positive")

    model_name = model.upper()
    _sampler_for_model(model_name)

    origin_nodes = _origin_nodes(original_graph)
    if not origin_nodes:
        return set()

    if len(origin_nodes) > 1:
        ell = ell * (1.0 + math.log(2.0) / math.log(len(origin_nodes)))
        if progress:
            progress(f"MF-IMM: adjusted ell={ell:.6f}")

    rng = random.Random(seed)
    target_k = min(k, len(origin_nodes))
    lower_bound = _estimate_imm_lower_bound(
        origin_nodes,
        generated_graph,
        model_name,
        target_k,
        epsilon,
        ell,
        rng,
        progress=progress,
    )
    if progress:
        progress(f"MF-IMM: lower_bound={lower_bound:.4f}")
    theta = _imm_theta(len(origin_nodes), target_k, epsilon, ell, lower_bound)
    if max_theta is not None:
        if max_theta <= 0:
            raise ValueError("max_theta must be positive when provided")
        original_theta = theta
        theta = min(theta, max_theta)
        if progress and theta != original_theta:
            progress(f"MF-IMM: theta capped from {original_theta} to {theta}")
    if progress:
        progress(f"MF-IMM: final theta={theta}")

    rr_sets = _sample_rr_sets_with_rng(
        origin_nodes,
        generated_graph,
        model_name,
        theta,
        rng,
        progress=progress,
        progress_label="MF-IMM final sampling",
    )
    if progress:
        progress(f"MF-IMM: running max-cover greedy with k={target_k}")
    selected, _ = _max_cover_greedy(origin_nodes, rr_sets, target_k)
    if progress:
        progress(f"MF-IMM: selected {len(selected)} seeds")
    return selected
