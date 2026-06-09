import copy
import diffusionModel
import random

def LFGreedy(graph, realizationCollection, k):
    seeds = set()
    influence = 0
    nodeList = []
    for node in graph.nodes:
        seed = set()
        seed.add((node, "origin"))
        seedInfluence = diffusionModel.compute(realizationCollection, seed)
        nodeList.append((node, seedInfluence))
    nodeList.sort(key=lambda nodeList:nodeList[1], reverse=True)
    seed = nodeList[0][0]
    influence = nodeList[0][1]
    seeds.add((seed, "origin"))
    del(nodeList[0])
    print("k =  1, Influence = " + str(diffusionModel.compute(realizationCollection, seeds)))

    for i in range(2, k + 1):
        check = False
        while not check:
            seed = nodeList[0][0]
            currentSeeds = copy.deepcopy(seeds)
            currentSeeds.add((seed, "origin"))
            nodeList[0] = (seed, diffusionModel.compute(realizationCollection, currentSeeds) - influence)
            nodeList.sort(key=lambda nodeList: nodeList[1], reverse=True)
            check = (seed == nodeList[0][0])
        seeds.add((nodeList[0][0], "origin"))
        influence += nodeList[0][1]
        del (nodeList[0])
        print("k = ", str(i), ", Influence = " + str(diffusionModel.compute(realizationCollection, seeds)))
    return seeds

def greedy(graph, realizationCollection, k):
    seeds = set()
    originNodes = set()
    for node in graph.nodes:
            originNodes.add((node, "origin"))
    for i in range(1, k + 1):
        candidate = originNodes - seeds
        candidate = list(candidate)
        max_seed = 0
        max_seed_result = 0
        for node in candidate:
            current_seeds = copy.deepcopy(seeds)
            current_seeds.add(node)
            current_seeds_result = diffusionModel.compute(realizationCollection, current_seeds)
            if current_seeds_result > max_seed_result:
                max_seed = node
                max_seed_result = current_seeds_result
        seeds.add(max_seed)
        print("k = ", str(i), ", Influence = " + str(diffusionModel.compute(realizationCollection, seeds)))
    return seeds

def greedyWDM(graph, realizationCollection, k):
    seeds_origin = set()
    seeds = set()
    for i in range(1, k + 1):
        candidate = graph.nodes - seeds_origin
        candidate = list(candidate)
        max_seed = 0
        max_seed_result = 0
        for node in candidate:
            current_seeds = copy.deepcopy(seeds_origin)
            current_seeds.add(node)
            current_seeds_result = diffusionModel.computeWithoutDM(graph, current_seeds)
            if current_seeds_result > max_seed_result:
                max_seed = node
                max_seed_result = current_seeds_result
        seeds_origin.add(max_seed)
        seeds.add((max_seed, "origin"))
        print("k = ", str(i), ", Influence = " + str(diffusionModel.compute(realizationCollection, seeds)))
    return seeds

def maxDegree(graph, realizationCollection, k):
    node_degree = {}
    for node in graph.nodes:
        node_degree[node] = len(graph.get_children(node))
    node_degree = sorted(node_degree.items(), key=lambda item: item[1], reverse=True)
    seeds = set()
    for i in range(1, k + 1):
        node = node_degree[i-1][0]
        seeds.add((node, "origin"))
        print("k = ", str(i), ", Influence = " + str(diffusionModel.compute(realizationCollection, seeds)))
    return seeds

def randomm(graph, realizationCollection, k):
    candidate = copy.deepcopy(graph.nodes)
    candidate = list(candidate)
    seeds = set()
    for i in range(1, k + 1):
        node = random.choice(candidate)
        candidate.remove(node)
        seeds.add((node, "origin"))
        print("k = ", str(i), ", Influence = " + str(diffusionModel.compute(realizationCollection, seeds)))
    return seeds