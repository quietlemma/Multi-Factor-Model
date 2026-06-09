import random
import pickle
import numpy as np
import copy
import numpy as np

# read a directed graph
def readGraph_direct(path):
    parentss = {}
    children = {}
    edges = set()
    nodes = set()
    f = open(path, 'r')
    for line in f.readlines():
        line = line.strip()
        if not len(line) or line.startswith('#'):
            continue
        row = line.split()
        src = int(row[0])
        dst = int(row[1])
        nodes.add(src)
        nodes.add(dst)
        if children.get(src) is None:
            children[src] = set()
        if parentss.get(dst) is None:
            parentss[dst] = set()
        edges.add((src, dst))
        children[src].add(dst)
        parentss[dst].add(src)
    return Graph(nodes, edges, children, parentss)

# read an undirected graph
def readGraph_undirect(path):
    parentss = {}
    children = {}
    edges = set()
    nodes = set()
    f = open(path, 'r')
    for line in f.readlines():
        line = line.strip()
        if not len(line) or line.startswith('#'):
            continue
        row = line.split()
        src = int(row[0])
        dst = int(row[1])
        nodes.add(src)
        nodes.add(dst)
        if children.get(src) is None:
            children[src] = set()
        if children.get(dst) is None:
            children[dst] = set()
        if parentss.get(src) is None:
            parentss[src] = set()
        if parentss.get(dst) is None:
            parentss[dst] = set()

        edges.add((src, dst))
        edges.add((dst, src))
        children[src].add(dst)
        children[dst].add(src)
        parentss[src].add(dst)
        parentss[dst].add(src)
    return Graph(nodes, edges, children, parentss)

class Graph:
    nodes = None
    edges = None
    children = None
    parentss = None
    def __init__(self, nodes, edges, children, parentss):
        self.nodes = nodes
        self.edges = edges
        self.children = children
        self.parentss = parentss
    def get_children(self, node):
        itsChildren = self.children.get(node)
        if itsChildren is None:
            return set()
        return self.children[node]
    def get_parentss(self, node):
        itsParentss = self.parentss.get(node)
        if itsParentss is None:
            return set()
        return self.parentss[node]

def isHappened(prob):
    if prob == 1:
        return True
    if prob == 0:
        return False
    rand = random.random()
    if rand <= prob:
        return True
    else:
        return False

# Generate auxillary graph for the Multi-Factor IC model
# Number of factors is 2
def generateGraph_factor_IC_2(graph):
    parentss = {}
    children = {}
    edges = {}
    nodes = set()
    for node in graph.nodes:
        nodes.add((node, "origin"))
        if len(graph.get_parentss(node)) > 0:
            parentss[(node, "origin")] = set()
            # factor 1
            nodes.add((node, "factor_1"))
            parentss[(node, "origin")].add((node, "factor_1"))
            children[(node, "factor_1")] = set()
            children[(node, "factor_1")].add((node, "origin"))
            edges[((node, "factor_1"),(node, "origin"))] = random.uniform(0, 0.2)
            # factor 2
            nodes.add((node, "factor_2"))
            parentss[(node, "origin")].add((node, "factor_2"))
            children[(node, "factor_2")] = set()
            children[(node, "factor_2")].add((node, "origin"))
            edges[((node, "factor_2"), (node, "origin"))] = random.uniform(0.2, 0.3)

            parentss[(node, "factor_1")] = set()
            parentss[(node, "factor_2")] = set()
            for parent in graph.get_parentss(node):
                parentss[(node, "factor_1")].add((parent, "origin"))
                parentss[(node, "factor_2")].add((parent, "origin"))
                children[(parent, "origin")] = set()
                children[(parent, "origin")].add((node, "factor_1"))
                children[(parent, "origin")].add((node, "factor_2"))
                edges[((parent, "origin"), (node, "factor_1"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_2"))] = 1 / len(graph.get_parentss(node))
    return Graph(nodes, edges, children, parentss)

# Number of factors is 3
def generateGraph_factor_IC_3(graph):
    parentss = {}
    children = {}
    edges = {}
    nodes = set()
    for node in graph.nodes:
        nodes.add((node, "origin"))
        if len(graph.get_parentss(node)) > 0:
            parentss[(node, "origin")] = set()
            # factor 1
            nodes.add((node, "factor_1"))
            parentss[(node, "origin")].add((node, "factor_1"))
            children[(node, "factor_1")] = set()
            children[(node, "factor_1")].add((node, "origin"))
            edges[((node, "factor_1"),(node, "origin"))] = random.uniform(0, 0.2)
            # factor 2
            nodes.add((node, "factor_2"))
            parentss[(node, "origin")].add((node, "factor_2"))
            children[(node, "factor_2")] = set()
            children[(node, "factor_2")].add((node, "origin"))
            edges[((node, "factor_2"), (node, "origin"))] = random.uniform(0.2, 0.3)
            # factor 3
            nodes.add((node, "factor_3"))
            parentss[(node, "origin")].add((node, "factor_3"))
            children[(node, "factor_3")] = set()
            children[(node, "factor_3")].add((node, "origin"))
            edges[((node, "factor_3"), (node, "origin"))] = random.uniform(0, 0.2)

            parentss[(node, "factor_1")] = set()
            parentss[(node, "factor_2")] = set()
            parentss[(node, "factor_3")] = set()
            for parent in graph.get_parentss(node):
                parentss[(node, "factor_1")].add((parent, "origin"))
                parentss[(node, "factor_2")].add((parent, "origin"))
                parentss[(node, "factor_3")].add((parent, "origin"))
                children[(parent, "origin")] = set()
                children[(parent, "origin")].add((node, "factor_1"))
                children[(parent, "origin")].add((node, "factor_2"))
                children[(parent, "origin")].add((node, "factor_3"))
                edges[((parent, "origin"), (node, "factor_1"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_2"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_3"))] = 1 / len(graph.get_parentss(node))
    return Graph(nodes, edges, children, parentss)

# Number of factors is 5
def generateGraph_factor_IC_5(graph):
    parentss = {}
    children = {}
    edges = {}
    nodes = set()
    for node in graph.nodes:
        nodes.add((node, "origin"))
        if len(graph.get_parentss(node)) > 0:
            parentss[(node, "origin")] = set()
            # factor 1
            nodes.add((node, "factor_1"))
            parentss[(node, "origin")].add((node, "factor_1"))
            children[(node, "factor_1")] = set()
            children[(node, "factor_1")].add((node, "origin"))
            edges[((node, "factor_1"),(node, "origin"))] = random.uniform(0, 0.2)
            # factor 2
            nodes.add((node, "factor_2"))
            parentss[(node, "origin")].add((node, "factor_2"))
            children[(node, "factor_2")] = set()
            children[(node, "factor_2")].add((node, "origin"))
            edges[((node, "factor_2"), (node, "origin"))] = random.uniform(0.2, 0.3)
            # factor 3
            nodes.add((node, "factor_3"))
            parentss[(node, "origin")].add((node, "factor_3"))
            children[(node, "factor_3")] = set()
            children[(node, "factor_3")].add((node, "origin"))
            edges[((node, "factor_3"), (node, "origin"))] = random.uniform(0, 0.2)
            # factor 4
            nodes.add((node, "factor_4"))
            parentss[(node, "origin")].add((node, "factor_4"))
            children[(node, "factor_4")] = set()
            children[(node, "factor_4")].add((node, "origin"))
            edges[((node, "factor_4"), (node, "origin"))] = random.uniform(0, 0.2)
            # factor 5
            nodes.add((node, "factor_5"))
            parentss[(node, "origin")].add((node, "factor_5"))
            children[(node, "factor_5")] = set()
            children[(node, "factor_5")].add((node, "origin"))
            edges[((node, "factor_5"), (node, "origin"))] = random.uniform(0, 0.2)

            parentss[(node, "factor_1")] = set()
            parentss[(node, "factor_2")] = set()
            parentss[(node, "factor_3")] = set()
            parentss[(node, "factor_4")] = set()
            parentss[(node, "factor_5")] = set()
            for parent in graph.get_parentss(node):
                parentss[(node, "factor_1")].add((parent, "origin"))
                parentss[(node, "factor_2")].add((parent, "origin"))
                parentss[(node, "factor_3")].add((parent, "origin"))
                parentss[(node, "factor_4")].add((parent, "origin"))
                parentss[(node, "factor_5")].add((parent, "origin"))
                children[(parent, "origin")] = set()
                children[(parent, "origin")].add((node, "factor_1"))
                children[(parent, "origin")].add((node, "factor_2"))
                children[(parent, "origin")].add((node, "factor_3"))
                children[(parent, "origin")].add((node, "factor_4"))
                children[(parent, "origin")].add((node, "factor_5"))
                edges[((parent, "origin"), (node, "factor_1"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_2"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_3"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_4"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_5"))] = 1 / len(graph.get_parentss(node))
    return Graph(nodes, edges, children, parentss)

def generateRealizationIC(generatedGraph):
    nodes = copy.deepcopy(generatedGraph.nodes)
    edges = set()
    parentss = {}
    children = {}
    for edge in generatedGraph.edges:
        if random.random() <= generatedGraph.edges[edge]:
            edges.add(edge)
            if children.get(edge[0]) is None:
                children[edge[0]] = set()
            if parentss.get(edge[1]) is None:
                parentss[edge[1]] = set()
            children[edge[0]].add(edge[1])
            parentss[edge[1]].add(edge[1])
    return Graph(nodes, edges, children, parentss)

def generateCollectionOfRealizationIC(generatedGraph, number=1000):
    realizationCollection = []
    for i in range(number):
        realization = generateRealizationIC(generatedGraph)
        realizationCollection.append(realization)
    return realizationCollection

# Generate auxillary graph for the Multi-Factor LT model
# Number of factors is 2
def generateGraph_factor_LT_2(graph):
    parentss = {}
    children = {}
    edges = {}
    nodes = set()
    for node in graph.nodes:
        nodes.add((node, "origin"))
        if len(graph.get_parentss(node)) > 0:
            parentss[(node, "origin")] = set()
            # factor 1
            nodes.add((node, "factor_1"))
            parentss[(node, "origin")].add((node, "factor_1"))
            children[(node, "factor_1")] = set()
            children[(node, "factor_1")].add((node, "origin"))
            edges[((node, "factor_1"),(node, "origin"))] = random.uniform(0, 0.3)
            # factor 2
            nodes.add((node, "factor_2"))
            parentss[(node, "origin")].add((node, "factor_2"))
            children[(node, "factor_2")] = set()
            children[(node, "factor_2")].add((node, "origin"))
            edges[((node, "factor_2"), (node, "origin"))] = 1 - edges[((node, "factor_1"),(node, "origin"))]

            parentss[(node, "factor_1")] = set()
            parentss[(node, "factor_2")] = set()
            for parent in graph.get_parentss(node):
                parentss[(node, "factor_1")].add((parent, "origin"))
                parentss[(node, "factor_2")].add((parent, "origin"))
                children[(parent, "origin")] = set()
                children[(parent, "origin")].add((node, "factor_1"))
                children[(parent, "origin")].add((node, "factor_2"))
                edges[((parent, "origin"), (node, "factor_1"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_2"))] = 1 / len(graph.get_parentss(node))
    return Graph(nodes, edges, children, parentss)

# Number of factors is 3
def generateGraph_factor_LT_3(graph):
    parentss = {}
    children = {}
    edges = {}
    nodes = set()
    for node in graph.nodes:
        nodes.add((node, "origin"))
        if len(graph.get_parentss(node)) > 0:
            parentss[(node, "origin")] = set()
            # factor 1
            nodes.add((node, "factor_1"))
            parentss[(node, "origin")].add((node, "factor_1"))
            children[(node, "factor_1")] = set()
            children[(node, "factor_1")].add((node, "origin"))
            edges[((node, "factor_1"),(node, "origin"))] = random.uniform(0, 0.3)
            # factor 3
            nodes.add((node, "factor_3"))
            parentss[(node, "origin")].add((node, "factor_3"))
            children[(node, "factor_3")] = set()
            children[(node, "factor_3")].add((node, "origin"))
            edges[((node, "factor_3"), (node, "origin"))] = random.uniform(0, 0.3)
            # factor 2
            nodes.add((node, "factor_2"))
            parentss[(node, "origin")].add((node, "factor_2"))
            children[(node, "factor_2")] = set()
            children[(node, "factor_2")].add((node, "origin"))
            edges[((node, "factor_2"), (node, "origin"))] = 1 - edges[((node, "factor_1"),(node, "origin"))] - edges[((node, "factor_3"), (node, "origin"))]

            parentss[(node, "factor_1")] = set()
            parentss[(node, "factor_2")] = set()
            parentss[(node, "factor_3")] = set()
            for parent in graph.get_parentss(node):
                parentss[(node, "factor_1")].add((parent, "origin"))
                parentss[(node, "factor_2")].add((parent, "origin"))
                parentss[(node, "factor_3")].add((parent, "origin"))
                children[(parent, "origin")] = set()
                children[(parent, "origin")].add((node, "factor_1"))
                children[(parent, "origin")].add((node, "factor_2"))
                children[(parent, "origin")].add((node, "factor_3"))
                edges[((parent, "origin"), (node, "factor_1"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_2"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_3"))] = 1 / len(graph.get_parentss(node))
    return Graph(nodes, edges, children, parentss)

# Number of factors is 5
def generateGraph_factor_LT_5(graph):
    parentss = {}
    children = {}
    edges = {}
    nodes = set()
    for node in graph.nodes:
        nodes.add((node, "origin"))
        if len(graph.get_parentss(node)) > 0:
            parentss[(node, "origin")] = set()
            # factor 1
            nodes.add((node, "factor_1"))
            parentss[(node, "origin")].add((node, "factor_1"))
            children[(node, "factor_1")] = set()
            children[(node, "factor_1")].add((node, "origin"))
            edges[((node, "factor_1"),(node, "origin"))] = random.uniform(0, 0.2)
            # factor 3
            nodes.add((node, "factor_3"))
            parentss[(node, "origin")].add((node, "factor_3"))
            children[(node, "factor_3")] = set()
            children[(node, "factor_3")].add((node, "origin"))
            edges[((node, "factor_3"), (node, "origin"))] = random.uniform(0, 0.2)
            # factor 4
            nodes.add((node, "factor_4"))
            parentss[(node, "origin")].add((node, "factor_4"))
            children[(node, "factor_4")] = set()
            children[(node, "factor_4")].add((node, "origin"))
            edges[((node, "factor_4"), (node, "origin"))] = random.uniform(0, 0.2)
            # factor 5
            nodes.add((node, "factor_5"))
            parentss[(node, "origin")].add((node, "factor_5"))
            children[(node, "factor_5")] = set()
            children[(node, "factor_5")].add((node, "origin"))
            edges[((node, "factor_5"), (node, "origin"))] = random.uniform(0, 0.2)
            # factor 2
            nodes.add((node, "factor_2"))
            parentss[(node, "origin")].add((node, "factor_2"))
            children[(node, "factor_2")] = set()
            children[(node, "factor_2")].add((node, "origin"))
            edges[((node, "factor_2"), (node, "origin"))] = 1 - edges[((node, "factor_1"),(node, "origin"))] - edges[((node, "factor_3"), (node, "origin"))] - edges[((node, "factor_4"), (node, "origin"))] - edges[((node, "factor_5"), (node, "origin"))]

            parentss[(node, "factor_1")] = set()
            parentss[(node, "factor_2")] = set()
            parentss[(node, "factor_3")] = set()
            parentss[(node, "factor_4")] = set()
            parentss[(node, "factor_5")] = set()
            for parent in graph.get_parentss(node):
                parentss[(node, "factor_1")].add((parent, "origin"))
                parentss[(node, "factor_2")].add((parent, "origin"))
                parentss[(node, "factor_3")].add((parent, "origin"))
                parentss[(node, "factor_4")].add((parent, "origin"))
                parentss[(node, "factor_5")].add((parent, "origin"))
                children[(parent, "origin")] = set()
                children[(parent, "origin")].add((node, "factor_1"))
                children[(parent, "origin")].add((node, "factor_2"))
                children[(parent, "origin")].add((node, "factor_3"))
                children[(parent, "origin")].add((node, "factor_4"))
                children[(parent, "origin")].add((node, "factor_5"))
                edges[((parent, "origin"), (node, "factor_1"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_2"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_3"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_4"))] = 1 / len(graph.get_parentss(node))
                edges[((parent, "origin"), (node, "factor_5"))] = 1 / len(graph.get_parentss(node))
    return Graph(nodes, edges, children, parentss)

def generateRealizationLT(generatedGraph):
    nodes = copy.deepcopy(generatedGraph.nodes)
    edges = set()
    parentss = {}
    children = {}
    for node in generatedGraph.nodes:
        if len(generatedGraph.get_parentss(node)) == 0:
            continue
        plist = []
        plistNumber = []
        number = 0
        plistDistribution = []
        for parent in generatedGraph.get_parentss(node):
            plist.append(parent)
            plistNumber.append(number)
            number += 1
            plistDistribution.append(generatedGraph.edges[(parent, node)])
        plistDistribution = np.array(plistDistribution)
        selectedParentNumber = np.random.choice(plistNumber, p=plistDistribution.ravel())
        src = plist[selectedParentNumber]
        dst = node
        if children.get(src) is None:
            children[src] = set()
        if parentss.get(dst) is None:
            parentss[dst] = set()
        children[src].add(dst)
        parentss[dst].add(src)
    return Graph(nodes, edges, children, parentss)

def generateCollectionOfRealizationLT(generatedGraph, number=1000):
    realizationCollection = []
    for i in range(number):
        realization = generateRealizationLT(generatedGraph)
        realizationCollection.append(realization)
    return realizationCollection

# Generate auxillary graph for the Multi-Factor Triggering model
def _trigger_origin_probability(factor_index, factor_count):
    if factor_count == 2:
        return random.uniform(0, 0.3) if factor_index == 1 else random.uniform(0.2, 0.3)
    if factor_count == 3:
        return random.uniform(0, 0.3) if factor_index in {1, 3} else random.uniform(0.2, 0.4)
    return random.uniform(0, 0.2) if factor_index != 2 else random.uniform(0.2, 0.4)

def _generateGraph_factor_TRIGGER(graph, factor_count):
    parentss = {}
    children = {}
    edges = {}
    nodes = set()
    trigger_probabilities = {}
    for node in graph.nodes:
        origin_node = (node, "origin")
        nodes.add(origin_node)
        if len(graph.get_parentss(node)) > 0:
            parentss[origin_node] = set()
            origin_trigger_probability = 0
            factor_nodes = []
            for factor_index in range(1, factor_count + 1):
                factor_node = (node, "factor_" + str(factor_index))
                factor_nodes.append(factor_node)
                nodes.add(factor_node)
                parentss[origin_node].add(factor_node)
                children[factor_node] = set()
                children[factor_node].add(origin_node)
                origin_trigger_probability = max(
                    origin_trigger_probability,
                    _trigger_origin_probability(factor_index, factor_count),
                )

            trigger_probabilities[origin_node] = origin_trigger_probability
            for factor_node in factor_nodes:
                edges[(factor_node, origin_node)] = origin_trigger_probability
                parentss[factor_node] = set()
                factor_trigger_probability = 1 / len(graph.get_parentss(node))
                trigger_probabilities[factor_node] = factor_trigger_probability
                for parent in graph.get_parentss(node):
                    parent_node = (parent, "origin")
                    parentss[factor_node].add(parent_node)
                    if children.get(parent_node) is None:
                        children[parent_node] = set()
                    children[parent_node].add(factor_node)
                    edges[(parent_node, factor_node)] = factor_trigger_probability
    generatedGraph = Graph(nodes, edges, children, parentss)
    generatedGraph.trigger_probabilities = trigger_probabilities
    return generatedGraph

def generateGraph_factor_TRIGGER_2(graph):
    return _generateGraph_factor_TRIGGER(graph, 2)

def generateGraph_factor_TRIGGER_3(graph):
    return _generateGraph_factor_TRIGGER(graph, 3)

def generateGraph_factor_TRIGGER_5(graph):
    return _generateGraph_factor_TRIGGER(graph, 5)

def _getTriggerProbability(generatedGraph, node):
    trigger_probabilities = getattr(generatedGraph, "trigger_probabilities", None)
    if trigger_probabilities is not None and node in trigger_probabilities:
        return trigger_probabilities[node]
    parents = generatedGraph.get_parentss(node)
    if len(parents) == 0:
        return 0
    return max(generatedGraph.edges.get((parent, node), 0) for parent in parents)

def generateRealizationTrigger(generatedGraph):
    nodes = copy.deepcopy(generatedGraph.nodes)
    edges = set()
    parentss = {}
    children = {}
    for node in generatedGraph.nodes:
        parents = generatedGraph.get_parentss(node)
        if len(parents) == 0:
            continue
        if random.random() > _getTriggerProbability(generatedGraph, node):
            continue
        for parent in parents:
            edge = (parent, node)
            edges.add(edge)
            if children.get(parent) is None:
                children[parent] = set()
            if parentss.get(node) is None:
                parentss[node] = set()
            children[parent].add(node)
            parentss[node].add(parent)
    return Graph(nodes, edges, children, parentss)

def generateCollectionOfRealizationTrigger(generatedGraph, number=1000):
    realizationCollection = []
    for i in range(number):
        realization = generateRealizationTrigger(generatedGraph)
        realizationCollection.append(realization)
    return realizationCollection
