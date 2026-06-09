import tools
import random
import copy
from collections import deque

# def computeIC(graph, seeds, R=200):
#     influence = 0
#     for i in range(R):
#         queue = []
#         queue.extend(seeds)
#         checked = copy.deepcopy(seeds)
#         while len(queue) != 0:
#             current_node = queue.pop(0)
#             children = graph.get_children(current_node)
#             for child in children:
#                 if child not in checked:
#                     rate = graph.edges[(current_node, child)]
#                     if tools.isHappened(rate):
#                         checked.add(child)
#                         queue.append(child)
#         influence += len(checked)
#     influence = influence/R
#     return influence

# def computeLT(graph, seeds, R=200):
#     influence = 0
#     for i in range(R):
#         threshold = {}
#         for node in graph.nodes:
#             threshold[node] = random.random()
#         sum_weight = {}
#         for node in graph.nodes:
#             sum_weight[node] = 0
#         influence_node = copy.deepcopy(seeds)
#         new_activate = copy.deepcopy(seeds)
#         while len(new_activate):
#             Snew = set()
#             for u in new_activate:
#                 for v in graph.get_children(u):
#                     if v not in influence_node:
#                         sum_weight[v] += graph.edges[(u, v)]
#                         if sum_weight[v] >= threshold[v]:
#                             Snew.add(v)
#                             influence_node.add(v)
#             new_activate = copy.deepcopy(Snew)
#         influence += len(influence_node)
#     influence = influence/R
#     return influence

def compute(realizationCollection, seeds):
    if hasattr(realizationCollection, "compute_influence"):
        return realizationCollection.compute_influence(seeds)

    influence = 0
    for generatedGraph in realizationCollection:
        queue = deque(seeds)
        checked = set(seeds)
        while len(queue) != 0:
            current_node = queue.popleft()
            children = generatedGraph.get_children(current_node)
            for child in children:
                if child not in checked:
                    checked.add(child)
                    queue.append(child)
        influence += countOriginNode(checked)
    influence = influence/len(realizationCollection)
    return influence

def computeWithoutDM(graph, seeds):
    queue = deque(seeds)
    checked = set(seeds)
    while len(queue) != 0:
        current_node = queue.popleft()
        children = graph.get_children(current_node)
        for child in children:
            if child not in checked:
                checked.add(child)
                queue.append(child)
    return len(checked)

def countOriginNode(nodeset):
    count = 0
    for node in nodeset:
        if node[1] == "origin":
            count += 1
    return count
