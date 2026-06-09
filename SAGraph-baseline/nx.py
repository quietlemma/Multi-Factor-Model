import networkx as nx
import random
from collections import defaultdict

import ndlib.models.epidemics as ep
import ndlib.models.ModelConfig as mc
import json

def _selected_weibo_interactions(static_data, dynamic_data, seed_set):
    selected = []
    for user_id in dynamic_data:
        if user_id not in seed_set:
            continue
        interact_list = list(set([uttr["interact_id"] for uttr in dynamic_data[user_id]]))
        reserved_list = []
        backup_list = []
        for interact_id in interact_list:
            interact_id = str(interact_id)
            if interact_id not in static_data:
                continue
            if static_data[interact_id]["user_followers"] > 100000:
                reserved_list.append(interact_id)
            else:
                backup_list.append(interact_id)
        backup_list = random.sample(backup_list, min(len(backup_list), 100))
        if len(reserved_list) > 0:
            reserved_list.extend(backup_list)
        else:
            reserved_list = backup_list
        for item in dynamic_data[user_id]:
            interact_id = str(item["interact_id"])
            if (user_id in seed_set and interact_id in reserved_list) or interact_id in seed_set:
                selected.append((user_id, interact_id))
    return selected

def _add_factor_edge(g, config, source, target, weight=1):
    g.add_edge(source, target, weight=weight)
    config.add_edge_configuration("threshold", (source, target), weight)

def weibo_multifactor(static_data, dynamic_data, seed_set, vertex_dict, vertex2userdict):
    g = nx.DiGraph()
    config = mc.Configuration()
    selected_interactions = _selected_weibo_interactions(static_data, dynamic_data, seed_set)
    selected_pairs = set(selected_interactions)
    user_vertices = set()
    interaction_sources = defaultdict(set)

    for user_id, interact_id in selected_interactions:
        from_vertex = vertex_dict[user_id]
        to_vertex = vertex_dict[interact_id]
        user_vertices.add(from_vertex)
        user_vertices.add(to_vertex)
        g.add_node(from_vertex, node_type="user", user_id=user_id)
        g.add_node(to_vertex, node_type="user", user_id=interact_id)
        interaction_sources[to_vertex].add(from_vertex)        

    for to_vertex, from_vertices in interaction_sources.items():
        factor_node = ("factor", "interaction", to_vertex)
        g.add_node(
            factor_node,
            node_type="factor",
            factor_type="interaction",
            target=to_vertex,
            sources=sorted(from_vertices),
        )
        _add_factor_edge(g, config, factor_node, to_vertex, 1)
        for from_vertex in from_vertices:
            _add_factor_edge(g, config, from_vertex, factor_node, 1)

    user_vertices = list(user_vertices)
    interest_sources = defaultdict(set)
    for i in range(len(user_vertices)):
        for j in range(i + 1, len(user_vertices)):
            from_vertex = user_vertices[i]
            to_vertex = user_vertices[j]
            user_from = vertex2userdict[from_vertex]
            user_to = vertex2userdict[to_vertex]
            from_interests = set(static_data[user_from].get("user_interests", []))
            to_interests = set(static_data[user_to].get("user_interests", []))
            common_interests = from_interests & to_interests
            if not common_interests:
                continue
            if (user_from, user_to) not in selected_pairs and (user_to, user_from) not in selected_pairs:
                continue
            interest_sources[to_vertex].add(from_vertex)
            interest_sources[from_vertex].add(to_vertex)
    
    for to_vertex, from_vertices in interest_sources.items():
        factor_node = ("factor", "interest", to_vertex)
        g.add_node(
        factor_node,
        node_type="factor",
        factor_type="interest",
        target=to_vertex,
        sources=sorted(from_vertices),
    )    
        _add_factor_edge(g, config, factor_node, to_vertex, 1)
        for from_vertex in from_vertices:
            _add_factor_edge(g, config, from_vertex, factor_node, 1)

    return g, config

connSW(1000, 0.1)
