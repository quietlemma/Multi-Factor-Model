import networkx as nx
import numpy as np
import ndlib
import ndlib.models.epidemics as ep
import ndlib.models.ModelConfig as mc
import statistics as s
import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from random import uniform, seed

from collections import Counter
import operator
import copy
from xflow.diffusion.SI import SI
from xflow.diffusion.IC import IC
from xflow.diffusion.LT import LT
from xflow.diffusion.SI import SI_multifactor
from xflow.diffusion.IC import IC_multifactor
from xflow.diffusion.LT import LT_multifactor

# random

# baselines: simulation based

def _user_candidates(g):
    candidates = [
        node for node, data in g.nodes(data=True)
        if data.get("node_type") == "user"
    ]
    if candidates:
        return candidates
    return [node for node in g.nodes() if not isinstance(node, tuple)]

def _weighted_graph_copy(g, config):
    g_copy = g.__class__()
    g_copy.add_nodes_from(g.nodes(data=True))
    g_copy.add_edges_from(g.edges)
    for a, b in g_copy.edges():
        weight = config.config["edges"]['threshold'][(a, b)]
        g_copy[a][b]['weight'] = weight
    return g_copy

def _select_top_user(score, candidates):
    user_score = {node: score[node] for node in candidates if node in score}
    if not user_score:
        return None
    return sorted(user_score, key=user_score.get, reverse=True)[0]

def _diffusion_user_spread(g, config, seeds, rounds=1, model='SI', beta=0.1):
    if model == "IC":
        result = IC_multifactor(g, config, seeds, rounds)
    elif model == "LT":
        result = LT_multifactor(g, config, seeds, rounds)
    elif model == "SI":
        result = SI_multifactor(g, config, seeds, rounds, beta)
    else:
        raise ValueError("Unsupported diffusion model: %s" % model)
    return s.mean(result)

def greedy_multifactor(g, config, budget, rounds=1, model='IC', beta=0.1):
    selected = []
    candidates = _user_candidates(g)
    budget = min(budget, len(candidates))

    for _ in range(budget):
        max_spread = 0
        index = -1
        for node in candidates:
            seeds = selected + [node]
            spread = _diffusion_user_spread(g, config, seeds, rounds, model, beta)
            if spread > max_spread:
                max_spread = spread
                index = node

        if index == -1:
            break
        selected.append(index)
        candidates.remove(index)

    print(selected)
    return selected

def celf_multifactor(g, config, budget, rounds=1, model='IC', beta=1):
    candidates = _user_candidates(g)
    budget = min(budget, len(candidates))
    if budget <= 0:
        print([])
        return []

    marg_gain = [
        _diffusion_user_spread(g, config, [node], rounds, model, beta)
        for node in candidates
    ]
    Q = sorted(zip(candidates, marg_gain), key=lambda x: x[1], reverse=True)

    selected = [Q[0][0]]
    spread = Q[0][1]
    Q = Q[1:]

    for _ in range(budget - 1):
        check = False
        while Q and not check:
            current = Q[0][0]
            current_spread = _diffusion_user_spread(
                g, config, selected + [current], rounds, model, beta
            )
            Q[0] = (current, current_spread - spread)
            Q = sorted(Q, key=lambda x: x[1], reverse=True)
            check = Q[0][0] == current

        if not Q:
            break
        selected.append(Q[0][0])
        spread = _diffusion_user_spread(g, config, selected, rounds, model, beta)
        Q = Q[1:]

    print(selected)
    return selected

def celfpp_multifactor(g, config, budget, rounds=1, model='IC', beta=0.1):
    candidates = _user_candidates(g)
    budget = min(budget, len(candidates))
    if budget <= 0:
        print([])
        return []

    marg_gain = [
        _diffusion_user_spread(g, config, [node], rounds, model, beta)
        for node in candidates
    ]
    Q = sorted(zip(candidates, marg_gain), key=lambda x: x[1], reverse=True)

    selected = [Q[0][0]]
    spread = Q[0][1]
    Q = Q[1:]
    last_seed = selected[0]

    for _ in range(budget - 1):
        check = False
        while Q and not check:
            current, old_gain = Q[0][0], Q[0][1]

            if current != last_seed:
                current_spread = _diffusion_user_spread(
                    g, config, selected + [current], rounds, model, beta
                )
                new_gain = current_spread - spread
            else:
                new_gain = old_gain

            Q[0] = (current, new_gain)
            Q = sorted(Q, key=lambda x: x[1], reverse=True)
            check = Q[0][0] == current

        if not Q:
            break
        selected.append(Q[0][0])
        spread = _diffusion_user_spread(g, config, selected, rounds, model, beta)
        last_seed = Q[0][0]
        Q = Q[1:]

    print(selected)
    return selected

def eigen_multifactor(g, config, budget):
    g_eig = _weighted_graph_copy(g, config)
    candidates = _user_candidates(g_eig)
    budget = min(budget, len(candidates))
    eig = []

    for _ in range(budget):
        try:
            eigen = nx.eigenvector_centrality_numpy(g_eig, weight='weight')
        except nx.NetworkXException:
            eigen = nx.eigenvector_centrality(g_eig, weight='weight', max_iter=1000)
        selected = _select_top_user(eigen, candidates)
        if selected is None:
            break
        eig.append(selected)
        candidates.remove(selected)
        g_eig.remove_node(selected)

    print(eig)
    return eig

def degree_multifactor(g, config, budget):
    g_deg = _weighted_graph_copy(g, config)
    candidates = _user_candidates(g_deg)
    budget = min(budget, len(candidates))
    deg = []

    for _ in range(budget):
        degree = nx.centrality.degree_centrality(g_deg)
        selected = _select_top_user(degree, candidates)
        if selected is None:
            break
        deg.append(selected)
        candidates.remove(selected)
        g_deg.remove_node(selected)

    print(deg)
    return deg

def pi_multifactor(g, config, budget):
    g_greedy = _weighted_graph_copy(g, config)
    candidates = _user_candidates(g_greedy)
    budget = min(budget, len(candidates))
    result = []

    for _ in range(budget):
        n = g_greedy.number_of_nodes()
        if n == 0:
            break

        nodes = list(g_greedy.nodes())
        I = np.ones((n, 1))
        C = np.ones((n, n))
        N = np.ones((n, n))
        A = nx.convert_matrix.to_numpy_array(g_greedy, nodelist=nodes, weight='weight')

        for i in range(5):
            B = np.power(A, i + 1)
            D = C - B
            N = np.multiply(N, D)

        P = C - N
        pi = np.matmul(P, I)
        value = {nodes[i]: pi[i, 0] for i in range(n)}

        selected = _select_top_user(value, candidates)
        if selected is None:
            break
        result.append(selected)
        candidates.remove(selected)
        g_greedy.remove_node(selected)

    print(result)
    return result

def sigma_multifactor(g, config, budget):
    g_greedy = _weighted_graph_copy(g, config)
    candidates = _user_candidates(g_greedy)
    budget = min(budget, len(candidates))
    result = []

    for _ in range(budget):
        n = g_greedy.number_of_nodes()
        if n == 0:
            break

        nodes = list(g_greedy.nodes())
        I = np.ones((n, 1))
        A = nx.convert_matrix.to_numpy_array(g_greedy, nodelist=nodes, weight='weight')

        sigma = I
        for i in range(5):
            B = np.power(A, i + 1)
            C = np.matmul(B, I)
            sigma += C

        value = {nodes[i]: sigma[i, 0] for i in range(n)}

        selected = _select_top_user(value, candidates)
        if selected is None:
            break
        result.append(selected)
        candidates.remove(selected)
        g_greedy.remove_node(selected)

    print(result)
    return result

def RIS_multifactor(
        g,
        config,
        budget,
        rounds=100,
        model='IC',
        beta=0.1):

    user_nodes = _user_candidates(g)
    if len(user_nodes) == 0:
        print([])
        return []
    budget = min(
        budget,
        len(user_nodes)
    )

    # Generate RR sets
    R = [
        get_RRS(
            g,
            config,
            model
        )
        for _ in range(rounds)
    ]

    selected = []

    for _ in range(budget):

        flat_map = [
            node
            for subset in R
            for node in subset
        ]

        if len(flat_map) == 0:
            break

        seed = Counter(
            flat_map
        ).most_common(1)[0][0]

        selected.append(seed)

        removed_num = 0

        new_R = []

        for rrs in R:

            if seed in rrs:
                removed_num += 1
            else:
                new_R.append(rrs)

        R = new_R

        for _ in range(removed_num):

            R.append(
                get_RRS(
                    g,
                    config,
                    model
                )
            )

    print(selected)

    return selected


def get_RRS(g, config, model='IC'):
    """
    Generate Reverse Reachable Set for multifactor graph.

    Supports:
        IC_multifactor
        LT_multifactor

    Only user nodes are returned in the final RRS.
    """

    # user node
    user_nodes = [
        node for node, data in g.nodes(data=True)
        if data.get("node_type") == "user"
    ]
    if len(user_nodes) == 0:
        return []
    source = random.choice(user_nodes)

    # ========================== IC ==========================
    if model == "IC":

        active_edges = []

        for u, v in g.edges():

            p = config.config["edges"]['threshold'][(u, v)]

            if uniform(0, 1) < p:
                active_edges.append((u, v))

        g_sub = nx.DiGraph()
        g_sub.add_nodes_from(g.nodes(data=True))
        g_sub.add_edges_from(active_edges)

    # ========================== LT ==========================
    elif model == "LT":

        active_edges = []

        for v in g.nodes():
            in_edges = list(g.in_edges(v))

            if len(in_edges) == 0:
                continue
            
            weights = [
                config.config["edges"]['threshold'][e]
                for e in in_edges
                ]
                
            total = sum(weights)
            
            if total == 0:
                continue
            
            probs = [w / total for w in weights]
            
            chosen_edge = random.choices(
                in_edges,
                weights=probs,
                k=1
                )[0]
                
            active_edges.append(chosen_edge)

        g_sub = nx.DiGraph()
        g_sub.add_nodes_from(g.nodes(data=True))
        g_sub.add_edges_from(active_edges)

    else:
        raise ValueError(
            f"Unsupported diffusion model: {model}"
        )

    if source not in g_sub:
        return [source]

    # Reverse Reachable
    g_rev = g_sub.reverse()

    reachable = set(
        nx.dfs_preorder_nodes(
            g_rev,
            source
        )
    )

    rrs = [
        node for node in reachable
        if g.nodes[node].get("node_type") == "user"
    ]

    return rrs