#pragma once

/// Node list
typedef std::vector<uint32_t> Nodelist;
/// Edge structure, neighbor id and the edge weight
typedef std::pair<uint32_t, float> Edge;
/// Edgelist structure from one source/target node
typedef std::vector<Edge> Edgelist;
/// Graph structure
typedef std::vector<Edgelist> Graph;
/// One forward reachable set
typedef std::vector<size_t> FRset;
/// A set of forward reachable sets
typedef std::vector<FRset> FRsets;
/// One reverse reachable set
typedef std::vector<uint32_t> RRset;
/// A set of reverse reachable sets
typedef std::vector<RRset> RRsets;

/// Cascade models: IC, IC2, IC2Mix, LT, OLO
enum CascadeModel { IC, IC2, IC2MIX, LT, OLO };

/// One k-list prefix result recorded from the same greedy pass.
typedef struct PrefixRecord
{
	int k = 0;
	size_t coveredRR = 0;
	double selfEstimatedInfluence = 0.0;
	double validationInfluence = 0.0;
	double approximationForKMax = -1.0;
	double runningTimeSec = 0.0;
	Nodelist seeds;
} PrefixRecord;

/// Node element with id and a property value
typedef struct NodeElement
{
	int id;
	double value;
} NodeEleType;

/// Smaller operation for node element
struct smaller
{
	bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const
	{
		return (Ele1.value < Ele2.value);
	}
};

/// Greater operation for node element
struct greater
{
	bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const
	{
		return (Ele1.value > Ele2.value);
	}
};
