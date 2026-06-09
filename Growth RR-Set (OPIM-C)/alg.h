#pragma once

class Alg
{
private:
	/// __numV: number of nodes in the graph.
	uint32_t __numV;
	/// __numE: number of edges in the graph.
	size_t __numE;
	/// __numRRsets: number of RR sets.
	size_t __numRRsets = 0;
	/// Fixed total RR budget for OPIM-C. 0 means use the original adaptive schedule.
	size_t __fixedRRsetsTotal = 0;
	/// Cap the greedy candidate pool by initial RR coverage. 0 means disabled.
	size_t __statefulCandidatePoolCap = 0;
	/// Upper bound in the last round for __mode=1.
	double __boundLast = DBL_MAX;
	/// The minimum upper bound among all rounds for __model=2.
	double __boundMin = DBL_MAX;
	/// Two hyper-graphs, one is used for selecting seeds and the other is used for validating influence. 
	THyperGraph __hyperG, __hyperGVldt;
	/// Result object.
	TResult& __tRes;
	/// Seed set.
	Nodelist __vecSeed;
	/// Maximum coverage by lazy updating.
	double max_cover_lazy(const int targetSize, const int mode = 2);
	/// Maximum coverage by maintaining the top-k marginal coverage.
	double max_cover_topk(const int targetSize);
	/// Maximum coverage.
	double max_cover(const int targetSize, const int mode = 2);
	/// For stateful models, build a fresh validation RR sample with a floor budget
	/// so prefix validation is less noisy than the seed-selection sample.
	size_t choose_stateful_validation_rrsets_per_graph(const size_t rrsetsPerGraph);
	double stateful_validation_influence_fresh(const Nodelist& vecSeed, const size_t rrsetsPerGraph);
	void fill_stateful_validation_prefixes(std::vector<PrefixRecord>& prefixRecords, const size_t rrsetsPerGraph);
	/// Maximum coverage with prefix checkpoints from one greedy pass.
	double max_cover_lazy_prefix(const int targetSize, const int mode, const std::vector<int>& kList,
	                             std::vector<PrefixRecord>& prefixRecords);
public:
	Alg(const Graph& graph, TResult& tRes) : __hyperG(graph), __hyperGVldt(graph), __tRes(tRes)
	{
		__numV = __hyperG.get_nodes();
		__numE = __hyperG.get_edges();
	}
	~Alg()
	{
	}
	/// Set cascade model.
	void set_cascade_model(const CascadeModel model);
	/// Restrict RR roots, candidate seeds, and influence estimates to origin nodes.
	void set_origin_node_count(const uint32_t originN);
	/// Fix the total RR budget for OPIM-C/OPIM-C k-list. 0 keeps the adaptive schedule.
	void set_fixed_rrsets_total(const size_t totalRRsets);
	/// Cap greedy candidates to the top initial-RR-coverage pool. 0 disables capping.
	void set_stateful_candidate_pool_cap(const size_t cap);
	/// Evaluate influence spread for the seed set constructed
	double effic_inf_valid_algo(const double evalDelta = 1e-3, const double evalEps = 0.01);
	/// Evaluate influence spread for a given seed set
	double effic_inf_valid_algo(const Nodelist vecSeed, const double evalDelta = 1e-3, const double evalEps = 0.01);
	/// OPIM: return the approximation guarantee alpha for the greedy algorithm when a given number of RR sets are generated.
	double opim(const int targetSize, const size_t numRRsets, const double delta, const int mode = 2,
	            const double evalDelta = 1e-3, const double evalEps = 0.01);
	/// OPIM-C: return (epsilon, delta)-approximate solution for influence maximization.  
	double opimc(const int targetSize, const double epsilon, const double delta, const int mode = 2,
	             const double evalDelta = 1e-3, const double evalEps = 0.01);
	/// OPIM prefix mode: sample once at k_max and report k-list prefixes from one greedy pass.
	double opim_klist(const std::vector<int>& kList, const size_t numRRsets, const double delta, const int mode,
	                  const double evalDelta, const double evalEps, const std::string& outFileName,
	                  const std::string& outFolder);
	/// OPIM-C prefix mode: run OPIM-C once at k_max and report k-list prefixes from one greedy pass.
	double opimc_klist(const std::vector<int>& kList, const double epsilon, const double delta, const int mode,
	                   const double evalDelta, const double evalEps, const std::string& outFileName,
	                   const std::string& outFolder);
};

using TAlg = Alg;
using PAlg = std::shared_ptr<TAlg>;
