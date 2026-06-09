#include "stdafx.h"

static size_t checked_size_from_double(const double value, const std::string& label)
{
	if (!std::isfinite(value) || value < 0.0
		|| value > static_cast<double>(std::numeric_limits<size_t>::max()))
	{
		std::cout << "Error: " << label << " is out of size_t range: " << value << std::endl;
		exit(1);
	}
	return static_cast<size_t>(value);
}

static void write_prefix_outputs(const std::string& outFolder, const std::string& outFileName,
                                 const std::vector<PrefixRecord>& records,
                                 const size_t rrsetsPerGraph,
                                 const size_t rrsetsTotalBudget)
{
	TIO::mkdir_absence(outFolder.c_str());
	const auto seedDir = outFolder + "/seed";
	TIO::mkdir_absence(seedDir.c_str());
	const auto csvPath = outFolder + "/" + outFileName + "_prefix.csv";
	std::ofstream csv(csvPath);
	csv << "k,approximation_for_kmax,time_sec,self_estimated_influence,validation_influence,covered_rr,rr_sets,rr_sets_total_budget,seed_file,seeds\n";
	for (const auto& record : records)
	{
		const auto seedFile = seedDir + "/seed_" + outFileName + "_k" + std::to_string(record.k);
		std::ofstream seedOut(seedFile);
		for (auto seed : record.seeds)
		{
			seedOut << seed << '\n';
		}
		seedOut.close();

		csv << record.k << ','
		    << record.approximationForKMax << ','
		    << record.runningTimeSec << ','
		    << record.selfEstimatedInfluence << ','
		    << record.validationInfluence << ','
			    << record.coveredRR << ','
			    << rrsetsPerGraph << ','
			    << rrsetsTotalBudget << ','
			    << seedFile << ','
		    << '"';
		for (size_t i = 0; i < record.seeds.size(); ++i)
		{
			if (i) csv << ' ';
			csv << record.seeds[i];
		}
		csv << '"' << '\n';
	}
	LogProgress("klist_outputs_written csv=" + csvPath);
}

typedef std::pair<size_t, uint32_t> CandidateCoveragePair;

struct CandidatePoolSelection
{
	bool enabled = false;
	size_t cap = 0;
	size_t nonzeroCandidates = 0;
	size_t selectedCandidates = 0;
	std::vector<uint8_t> allowed;
};

struct CandidatePoolMinHeapCompare
{
	bool operator()(const CandidateCoveragePair& lhs, const CandidateCoveragePair& rhs) const
	{
		if (lhs.first != rhs.first) return lhs.first > rhs.first;
		return lhs.second < rhs.second;
	}
};

static CandidatePoolSelection make_candidate_pool_selection(const THyperGraph& hyperG, const size_t candidatePoolCap)
{
	CandidatePoolSelection selection;
	const auto candidateN = hyperG.get_origin_nodes();
	if (candidatePoolCap == 0 || candidatePoolCap >= candidateN)
	{
		return selection;
	}

	selection.enabled = true;
	selection.cap = candidatePoolCap;
	std::priority_queue<CandidateCoveragePair, std::vector<CandidateCoveragePair>, CandidatePoolMinHeapCompare> heap;
	for (uint32_t candidate = 0; candidate < candidateN; ++candidate)
	{
		const auto coverage = hyperG._FRsets[candidate].size();
		if (coverage == 0) continue;
		selection.nonzeroCandidates++;
		CandidateCoveragePair item(coverage, candidate);
		if (heap.size() < candidatePoolCap)
		{
			heap.push(item);
			continue;
		}
		const auto& worst = heap.top();
		if (item.first > worst.first || (item.first == worst.first && item.second < worst.second))
		{
			heap.pop();
			heap.push(item);
		}
	}

	selection.allowed.assign(candidateN, 0);
	while (!heap.empty())
	{
		selection.allowed[heap.top().second] = 1;
		selection.selectedCandidates++;
		heap.pop();
	}
	return selection;
}

static void log_candidate_pool_selection(const std::string& label, const CandidatePoolSelection& selection)
{
	if (!selection.enabled) return;
	LogProgress(label
	            + " candidate_pool_cap=" + std::to_string(selection.cap)
	            + " nonzero_candidates=" + std::to_string(selection.nonzeroCandidates)
	            + " selected_candidates=" + std::to_string(selection.selectedCandidates));
}

double Alg::max_cover_lazy(const int targetSize, const int mode)
{
	LogProgress("max_cover_lazy_start k=" + std::to_string(targetSize)
	            + " rrsets=" + std::to_string(__numRRsets));
	if (__hyperG.uses_stateful_seed_selection_model())
	{
		(void)mode;
		const std::vector<int> wanted{targetSize};
		const auto prefixRecords = __hyperG.build_seedset_prefix_stateful(wanted, __statefulCandidatePoolCap);
		__vecSeed.clear();
		if (!prefixRecords.empty())
		{
			__vecSeed = prefixRecords.back().seeds;
			const auto finalInf = prefixRecords.back().selfEstimatedInfluence;
			__boundLast = finalInf / (1.0 - 1.0 / exp(1.0));
			__boundMin = __boundLast;
			LogProgress("max_cover_stateful_done k=" + std::to_string(targetSize)
			            + " selected=" + std::to_string(__vecSeed.size())
			            + " influence=" + std::to_string(finalInf));
			return finalInf;
		}
		__boundLast = DBL_MAX;
		__boundMin = DBL_MAX;
		return 0.0;
	}
	// mode: optimization mode.
	// 0->no optimization, 
	// 1->optimization with the upper bound in the last round,
	// 2->optimization with minimum upper bound among all rounds [Default].
	__boundLast = DBL_MAX, __boundMin = DBL_MAX;
	const auto candidateN = __hyperG.get_origin_nodes();
	const auto candidatePool = make_candidate_pool_selection(__hyperG, __statefulCandidatePoolCap);
	log_candidate_pool_selection("candidate_pool_lazy", candidatePool);
	FRset coverage(candidateN, 0);
	size_t maxDeg = 0;
	for (auto i = candidateN; i--;)
	{
		if (candidatePool.enabled && !candidatePool.allowed[i]) continue;
		const auto deg = __hyperG._FRsets[i].size();
		coverage[i] = deg;
		if (deg > maxDeg) maxDeg = deg;
	}
	RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
	for (auto i = candidateN; i--;)
	{
		if (candidatePool.enabled && !candidatePool.allowed[i]) continue;
		if (coverage[i] == 0) continue;
		degMap[coverage[i]].push_back(i);
	}
	size_t sumInf = 0;

	// check if an edge is removed
	std::vector<bool> edgeMark(__numRRsets, false);

	__vecSeed.clear();
	for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
	{
		auto& vecNode = degMap[deg];
		for (auto idx = vecNode.size(); idx--;)
		{
			auto argmaxIdx = vecNode[idx];
			const auto currDeg = coverage[argmaxIdx];
			if (deg > currDeg)
			{
				degMap[currDeg].push_back(argmaxIdx);
				continue;
			}
			if (mode == 2 || (mode == 1 && __vecSeed.size() == targetSize))
			{
				// Find upper bound
				auto topk = targetSize;
				auto degBound = deg;
				FRset vecBound(targetSize);
				// Initialize vecBound
				auto idxBound = idx + 1;
				while (topk && idxBound--)
				{
					vecBound[--topk] = coverage[degMap[degBound][idxBound]];
				}
				while (topk && --degBound)
				{
					idxBound = degMap[degBound].size();
					while (topk && idxBound--)
					{
						vecBound[--topk] = coverage[degMap[degBound][idxBound]];
					}
				}
				make_min_heap(vecBound);

				// Find the top-k marginal coverage
				auto flag = topk == 0;
				while (flag && idxBound--)
				{
					const auto currDegBound = coverage[degMap[degBound][idxBound]];
					if (vecBound[0] >= degBound)
					{
						flag = false;
					}
					else if (vecBound[0] < currDegBound)
					{
						min_heap_replace_min_value(vecBound, currDegBound);
					}
				}
				while (flag && --degBound)
				{
					idxBound = degMap[degBound].size();
					while (flag && idxBound--)
					{
						const auto currDegBound = coverage[degMap[degBound][idxBound]];
						if (vecBound[0] >= degBound)
						{
							flag = false;
						}
						else if (vecBound[0] < currDegBound)
						{
							min_heap_replace_min_value(vecBound, currDegBound);
						}
					}
				}
				__boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * candidateN / __numRRsets;
				if (__boundMin > __boundLast) __boundMin = __boundLast;
			}
			if (__vecSeed.size() >= targetSize)
			{
				// Top-k influential nodes constructed
				const auto finalInf = 1.0 * sumInf * candidateN / __numRRsets;
				std::cout << "  >>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << __boundMin <<
					", last-bound: " << __boundLast << '\n';
				LogProgress("max_cover_lazy_done k=" + std::to_string(targetSize)
				            + " selected=" + std::to_string(__vecSeed.size())
				            + " influence=" + std::to_string(finalInf));
				return finalInf;
			}
			sumInf += currDeg;
			__vecSeed.push_back(argmaxIdx);
			coverage[argmaxIdx] = 0;
			for (auto edgeIdx : __hyperG._FRsets[argmaxIdx])
			{
				if (edgeMark[edgeIdx]) continue;
				edgeMark[edgeIdx] = true;
				for (auto nodeIdx : __hyperG._RRsets[edgeIdx])
				{
					if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
					coverage[nodeIdx]--;
				}
			}
		}
		degMap.pop_back();
	}
	if (candidatePool.enabled && __vecSeed.size() < static_cast<size_t>(targetSize))
	{
		LogProgress("candidate_pool_exhausted k=" + std::to_string(targetSize)
		            + " selected=" + std::to_string(__vecSeed.size()));
	}
	return __numRRsets > 0 ? 1.0 * sumInf * candidateN / __numRRsets : 0.0;
}

double Alg::max_cover_lazy_prefix(const int targetSize, const int mode, const std::vector<int>& kList,
                                  std::vector<PrefixRecord>& prefixRecords)
{
	LogProgress("max_cover_lazy_prefix_start kmax=" + std::to_string(targetSize)
	            + " rrsets=" + std::to_string(__numRRsets));
	if (__hyperG.uses_stateful_seed_selection_model())
	{
		(void)mode;
		prefixRecords = __hyperG.build_seedset_prefix_stateful(kList, __statefulCandidatePoolCap);
		__vecSeed.clear();
		if (!prefixRecords.empty())
		{
			__vecSeed = prefixRecords.back().seeds;
			const auto finalInf = prefixRecords.back().selfEstimatedInfluence;
			__boundLast = finalInf / (1.0 - 1.0 / exp(1.0));
			__boundMin = __boundLast;
			LogProgress("max_cover_stateful_prefix_done kmax=" + std::to_string(targetSize)
			            + " selected=" + std::to_string(__vecSeed.size())
			            + " influence=" + std::to_string(finalInf));
			return finalInf;
		}
		__boundLast = DBL_MAX;
		__boundMin = DBL_MAX;
		return 0.0;
	}
	__boundLast = DBL_MAX, __boundMin = DBL_MAX;
	prefixRecords.clear();
	const auto candidateN = __hyperG.get_origin_nodes();
	const auto candidatePool = make_candidate_pool_selection(__hyperG, __statefulCandidatePoolCap);
	log_candidate_pool_selection("candidate_pool_prefix", candidatePool);
	FRset coverage(candidateN, 0);
	size_t maxDeg = 0;
	for (auto i = candidateN; i--;)
	{
		if (candidatePool.enabled && !candidatePool.allowed[i]) continue;
		const auto deg = __hyperG._FRsets[i].size();
		coverage[i] = deg;
		if (deg > maxDeg) maxDeg = deg;
	}
	RRsets degMap(maxDeg + 1);
	for (auto i = candidateN; i--;)
	{
		if (candidatePool.enabled && !candidatePool.allowed[i]) continue;
		if (coverage[i] == 0) continue;
		degMap[coverage[i]].push_back(i);
	}

	size_t nextCheckpoint = 0;
	while (nextCheckpoint < kList.size() && kList[nextCheckpoint] == 0)
	{
		PrefixRecord record;
		record.k = 0;
		prefixRecords.push_back(std::move(record));
		nextCheckpoint++;
	}

	size_t sumInf = 0;
	std::vector<bool> edgeMark(__numRRsets, false);
	__vecSeed.clear();

	auto record_prefix = [&]()
	{
		while (nextCheckpoint < kList.size()
			&& kList[nextCheckpoint] <= static_cast<int>(__vecSeed.size()))
		{
			PrefixRecord record;
			record.k = kList[nextCheckpoint];
			record.coveredRR = sumInf;
			record.selfEstimatedInfluence = 1.0 * sumInf * candidateN / __numRRsets;
			record.seeds.assign(__vecSeed.begin(), __vecSeed.begin() + record.k);
			prefixRecords.push_back(std::move(record));
			LogProgress("klist_prefix_recorded k=" + std::to_string(kList[nextCheckpoint])
			            + " self_estimated_influence="
			            + std::to_string(prefixRecords.back().selfEstimatedInfluence));
			nextCheckpoint++;
		}
	};

	for (auto deg = maxDeg; deg > 0; deg--)
	{
		auto& vecNode = degMap[deg];
		for (auto idx = vecNode.size(); idx--;)
		{
			auto argmaxIdx = vecNode[idx];
			const auto currDeg = coverage[argmaxIdx];
			if (deg > currDeg)
			{
				degMap[currDeg].push_back(argmaxIdx);
				continue;
			}
			if (mode == 2 || (mode == 1 && __vecSeed.size() == static_cast<size_t>(targetSize)))
			{
				auto topk = targetSize;
				auto degBound = deg;
				FRset vecBound(targetSize);
				auto idxBound = idx + 1;
				while (topk && idxBound--)
				{
					vecBound[--topk] = coverage[degMap[degBound][idxBound]];
				}
				while (topk && --degBound)
				{
					idxBound = degMap[degBound].size();
					while (topk && idxBound--)
					{
						vecBound[--topk] = coverage[degMap[degBound][idxBound]];
					}
				}
				make_min_heap(vecBound);

				auto flag = topk == 0;
				while (flag && idxBound--)
				{
					const auto currDegBound = coverage[degMap[degBound][idxBound]];
					if (vecBound[0] >= degBound)
					{
						flag = false;
					}
					else if (vecBound[0] < currDegBound)
					{
						min_heap_replace_min_value(vecBound, currDegBound);
					}
				}
				while (flag && --degBound)
				{
					idxBound = degMap[degBound].size();
					while (flag && idxBound--)
					{
						const auto currDegBound = coverage[degMap[degBound][idxBound]];
						if (vecBound[0] >= degBound)
						{
							flag = false;
						}
						else if (vecBound[0] < currDegBound)
						{
							min_heap_replace_min_value(vecBound, currDegBound);
						}
					}
				}
				__boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * candidateN / __numRRsets;
				if (__boundMin > __boundLast) __boundMin = __boundLast;
			}

			sumInf += currDeg;
			__vecSeed.push_back(argmaxIdx);
			coverage[argmaxIdx] = 0;
			for (auto edgeIdx : __hyperG._FRsets[argmaxIdx])
			{
				if (edgeMark[edgeIdx]) continue;
				edgeMark[edgeIdx] = true;
				for (auto nodeIdx : __hyperG._RRsets[edgeIdx])
				{
					if (coverage[nodeIdx] == 0) continue;
					coverage[nodeIdx]--;
				}
			}
			record_prefix();
			if (__vecSeed.size() >= static_cast<size_t>(targetSize))
			{
				const auto finalInf = 1.0 * sumInf * candidateN / __numRRsets;
				LogProgress("max_cover_lazy_prefix_done kmax=" + std::to_string(targetSize)
				            + " selected=" + std::to_string(__vecSeed.size())
				            + " influence=" + std::to_string(finalInf));
				return finalInf;
			}
		}
		degMap.pop_back();
	}
	if (candidatePool.enabled && __vecSeed.size() < static_cast<size_t>(targetSize))
	{
		LogProgress("candidate_pool_exhausted_prefix kmax=" + std::to_string(targetSize)
		            + " selected=" + std::to_string(__vecSeed.size()));
	}
	while (nextCheckpoint < kList.size())
	{
		PrefixRecord record;
		record.k = kList[nextCheckpoint];
		record.coveredRR = sumInf;
		record.selfEstimatedInfluence = 1.0 * sumInf * candidateN / __numRRsets;
		record.seeds = __vecSeed;
		prefixRecords.push_back(std::move(record));
		nextCheckpoint++;
	}
	return 1.0 * sumInf * candidateN / __numRRsets;
}

double Alg::max_cover_topk(const int targetSize)
{
	LogProgress("max_cover_topk_start k=" + std::to_string(targetSize)
	            + " rrsets=" + std::to_string(__numRRsets));
	const auto candidateN = __hyperG.get_origin_nodes();
	FRset coverage(candidateN, 0);
	size_t maxDeg = 0;
	for (auto i = candidateN; i--;)
	{
		const auto deg = __hyperG._FRsets[i].size();
		coverage[i] = deg;
		if (deg > maxDeg) maxDeg = deg;
	}
	RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
	for (auto i = candidateN; i--;)
	{
		//if (coverage[i] == 0) continue;
		degMap[coverage[i]].push_back(i);
	}
	Nodelist sortedNode(candidateN); // sortedNode: record the sorted nodes in ascending order of degree
	Nodelist nodePosition(candidateN); // nodePosition: record the position of each origin node in the sortedNode
	Nodelist degreePosition(maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
	uint32_t idxSort = 0;
	size_t idxDegree = 0;
	for (auto& nodes : degMap)
	{
		degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
		idxDegree++;
		for (auto& node : nodes)
		{
			nodePosition[node] = idxSort;
			sortedNode[idxSort++] = node;
		}
	}
	// check if an edge is removed
	std::vector<bool> edgeMark(__numRRsets, false);
	// record the total of top-k marginal gains
	size_t sumTopk = 0;
	for (auto deg = maxDeg + 1; deg--;)
	{
		if (degreePosition[deg] <= candidateN - targetSize)
		{
			sumTopk += deg * (degreePosition[deg + 1] - (candidateN - targetSize));
			break;
		}
		sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
	}
	__boundMin = 1.0 * sumTopk;
	__vecSeed.clear();
	size_t sumInf = 0;
	/*
	* sortedNode: position -> node
	* nodePosition: node -> position
	* degreePosition: degree -> position (start position of this degree)
	* coverage: node -> degree
	* e.g., swap the position of a node with the start position of its degree
	* swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
	*/
	for (auto k = targetSize; k--;)
	{
		const auto seed = sortedNode.back();
		sortedNode.pop_back();
		const auto newNumV = sortedNode.size();
		sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
		sumInf += coverage[seed];
		__vecSeed.push_back(seed);
		coverage[seed] = 0;
		for (auto edgeIdx : __hyperG._FRsets[seed])
		{
			if (edgeMark[edgeIdx]) continue;
			edgeMark[edgeIdx] = true;
			for (auto nodeIdx : __hyperG._RRsets[edgeIdx])
			{
				if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
				const auto currPos = nodePosition[nodeIdx]; // The current position
				const auto currDeg = coverage[nodeIdx]; // The current degree
				const auto startPos = degreePosition[currDeg]; // The start position of this degree
				const auto startNode = sortedNode[startPos]; // The node with the start position
				// Swap this node to the start position with the same degree, and update their positions in nodePosition
				std::swap(sortedNode[currPos], sortedNode[startPos]);
				nodePosition[nodeIdx] = startPos;
				nodePosition[startNode] = currPos;
				// Increase the start position of this degree by 1, and decrease the degree of this node by 1
				degreePosition[currDeg]++;
				coverage[nodeIdx]--;
				// If the start position of this degree is in top-k, reduce topk by 1
				if (startPos >= newNumV - targetSize) sumTopk--;
			}
		}
		__boundLast = 1.0 * (sumInf + sumTopk);
		if (__boundMin > __boundLast) __boundMin = __boundLast;
	}
	__boundMin *= 1.0 * candidateN / __numRRsets;
	__boundLast *= 1.0 * candidateN / __numRRsets;
	const auto finalInf = 1.0 * sumInf * candidateN / __numRRsets;
	std::cout << "  >>>[greedy-topk] influence: " << finalInf << ", min-bound: " << __boundMin <<
		", last-bound: " << __boundLast << '\n';
	LogProgress("max_cover_topk_done k=" + std::to_string(targetSize)
	            + " selected=" + std::to_string(__vecSeed.size())
	            + " influence=" + std::to_string(finalInf));
	return finalInf;
}

double Alg::max_cover(const int targetSize, const int mode)
{
	if (__hyperG.uses_stateful_seed_selection_model()) return max_cover_lazy(targetSize, mode);
	if (__statefulCandidatePoolCap > 0) return max_cover_lazy(targetSize, mode);
	if (targetSize >= 1000) return max_cover_topk(targetSize);
	return max_cover_lazy(targetSize, mode);
}

size_t Alg::choose_stateful_validation_rrsets_per_graph(const size_t rrsetsPerGraph)
{
	if (!__hyperG.uses_stateful_seed_selection_model())
		return rrsetsPerGraph;
	// Small stateful RR samples are extremely noisy and can collapse the validation
	// column to zero even when the prefix curve itself is reasonable. Keep selection
	// on the user-requested budget, but validate on a less noisy fresh RR sample.
	return std::max<size_t>(rrsetsPerGraph, 1000);
}

double Alg::stateful_validation_influence_fresh(const Nodelist& vecSeed, const size_t rrsetsPerGraph)
{
	if (!__hyperG.uses_stateful_seed_selection_model())
		return __hyperGVldt.self_inf_cal(vecSeed);

	const auto validationRRsets = choose_stateful_validation_rrsets_per_graph(rrsetsPerGraph);
	THyperGraph validator(__hyperG._graph);
	validator.set_cascade_model(__hyperG._cascadeModel);
	validator.set_origin_node_count(__hyperG.get_origin_nodes());
	validator.build_n_RRsets(validationRRsets);
	LogProgress("stateful_validation_fresh_rr"
	            " select_rr=" + std::to_string(rrsetsPerGraph)
	            + " validate_rr=" + std::to_string(validator.get_RR_sets_size()));
	return validator.self_inf_cal(vecSeed);
}

void Alg::fill_stateful_validation_prefixes(std::vector<PrefixRecord>& prefixRecords, const size_t rrsetsPerGraph)
{
	if (!__hyperG.uses_stateful_seed_selection_model())
	{
		for (auto& record : prefixRecords)
			record.validationInfluence = __hyperGVldt.self_inf_cal(record.seeds);
		return;
	}

	const auto validationRRsets = choose_stateful_validation_rrsets_per_graph(rrsetsPerGraph);
	THyperGraph validator(__hyperG._graph);
	validator.set_cascade_model(__hyperG._cascadeModel);
	validator.set_origin_node_count(__hyperG.get_origin_nodes());
	validator.build_n_RRsets(validationRRsets);
	LogProgress("stateful_validation_prefix_rr"
	            " select_rr=" + std::to_string(rrsetsPerGraph)
	            + " validate_rr=" + std::to_string(validator.get_RR_sets_size())
	            + " prefixes=" + std::to_string(prefixRecords.size()));
	for (auto& record : prefixRecords)
		record.validationInfluence = validator.self_inf_cal(record.seeds);
}

void Alg::set_cascade_model(const CascadeModel model)
{
	__hyperG.set_cascade_model(model);
	__hyperGVldt.set_cascade_model(model);
}

void Alg::set_origin_node_count(const uint32_t originN)
{
	__hyperG.set_origin_node_count(originN);
	__hyperGVldt.set_origin_node_count(originN);
}

void Alg::set_fixed_rrsets_total(const size_t totalRRsets)
{
	__fixedRRsetsTotal = totalRRsets;
}

void Alg::set_stateful_candidate_pool_cap(const size_t cap)
{
	__statefulCandidatePoolCap = cap;
}

double Alg::effic_inf_valid_algo(const double evalDelta, const double evalEps)
{
	return effic_inf_valid_algo(__vecSeed, evalDelta, evalEps);
}

double Alg::effic_inf_valid_algo(const Nodelist vecSeed, const double evalDelta, const double evalEps)
{
	Timer EvalTimer("Inf. Eval.");
	std::cout << "  >>>Evaluating influence with evaleps=" << evalEps
	          << ", evaldelta=" << evalDelta << " ...\n";
	const auto inf = __hyperG.effic_inf_valid_algo(vecSeed, evalDelta, evalEps);
	//const auto inf = __hyperGVldt.effic_inf_valid_algo(vecSeed);
	std::cout << "  >>>Down! influence: " << inf << ", time used (sec): " << EvalTimer.get_total_time() << '\n';
	return inf;
}

double Alg::opim(const int targetSize, const size_t numRRsets, const double delta, const int mode,
                 const double evalDelta, const double evalEps)
{
	LogProgress("opim_start k=" + std::to_string(targetSize)
	            + " samplesize=" + std::to_string(numRRsets)
	            + " delta=" + std::to_string(delta));
	Timer timerOPIM("OPIM");
	const double e = exp(1);
	const double approx = 1 - 1.0 / e;
	const double a1 = log(2.0 / delta);
	const double a2 = log(2.0 / delta);
	__hyperG.build_n_RRsets(numRRsets / 2); // R1
	__hyperGVldt.build_n_RRsets(numRRsets / 2); // R2
	__numRRsets = __hyperG.get_RR_sets_size();
	const auto candidateN = __hyperG.get_origin_nodes();
	const auto time1 = timerOPIM.get_operation_time();
	const auto infSelf = max_cover(targetSize, mode);
	const auto time2 = timerOPIM.get_operation_time();
	const auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
	const auto degVldt = infVldt * __numRRsets / candidateN;
	auto upperBound = infSelf / approx;
	if (mode == 1) upperBound = __boundLast;
	else if (mode == 2) upperBound = __boundMin;
	const auto upperDegOPT = upperBound * __numRRsets / candidateN;
	const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
	const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
	const auto approxOPIM = lowerSelect / upperOPT;
	//guaranteeOpt = opt_error(degVldt, degSelfUpper, delta);
	__tRes.set_approximation(approxOPIM);
	__tRes.set_running_time(timerOPIM.get_total_time());
	__tRes.set_influence(infVldt);
	__tRes.set_influence_original(infSelf);
	__tRes.set_seed_vec(__vecSeed);
	__tRes.set_RR_sets_size(__numRRsets);
	std::cout << "==>OPIM approx. (max-cover): " << approxOPIM << " (" << infSelf / upperBound << ")\n";
	std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
	print_memory_snapshot("opim_before_final_eval_release_rr");
	__hyperGVldt.release_memory();
	__hyperG.release_rr_index(false);
	print_memory_snapshot("opim_after_final_eval_release_rr");
	if (!__hyperG.uses_stateful_seed_selection_model())
		__tRes.set_influence(effic_inf_valid_algo(evalDelta, evalEps));
	LogProgress("opim_done k=" + std::to_string(targetSize)
	            + " influence=" + std::to_string(__tRes.get_influence()));
	return approxOPIM;
}

double Alg::opim_klist(const std::vector<int>& kList, const size_t numRRsets, const double delta, const int mode,
                       const double evalDelta, const double evalEps, const std::string& outFileName,
                       const std::string& outFolder)
{
	(void)evalDelta;
	(void)evalEps;
	if (kList.empty() || kList.back() <= 0)
	{
		std::cout << "Error: OPIM k-list needs at least one positive k.\n";
		exit(1);
	}
	const int targetSize = kList.back();
	LogProgress("opim_klist_start kmax=" + std::to_string(targetSize)
	            + " samplesize=" + std::to_string(numRRsets));
	Timer timerOPIM("OPIM-klist");
	const double e = exp(1);
	const double approx = 1 - 1.0 / e;
	const double a1 = log(2.0 / delta);
	const double a2 = log(2.0 / delta);
	__hyperG.build_n_RRsets(numRRsets / 2);
	__hyperGVldt.build_n_RRsets(numRRsets / 2);
	__numRRsets = __hyperG.get_RR_sets_size();
	const auto candidateN = __hyperG.get_origin_nodes();
	const auto time1 = timerOPIM.get_operation_time();
	std::vector<PrefixRecord> prefixRecords;
	const auto infSelf = max_cover_lazy_prefix(targetSize, mode, kList, prefixRecords);
	const auto time2 = timerOPIM.get_operation_time();
	double infVldt = 0.0;
	if (__hyperG.uses_stateful_seed_selection_model())
	{
		const auto validationRRsets = choose_stateful_validation_rrsets_per_graph(__numRRsets);
		THyperGraph validator(__hyperG._graph);
		validator.set_cascade_model(__hyperG._cascadeModel);
		validator.set_origin_node_count(__hyperG.get_origin_nodes());
		validator.build_n_RRsets(validationRRsets);
		LogProgress("stateful_validation_prefix_rr"
		            " select_rr=" + std::to_string(__numRRsets)
		            + " validate_rr=" + std::to_string(validator.get_RR_sets_size())
		            + " prefixes=" + std::to_string(prefixRecords.size()));
		infVldt = validator.self_inf_cal(__vecSeed);
		for (auto& record : prefixRecords)
			record.validationInfluence = validator.self_inf_cal(record.seeds);
		if (!prefixRecords.empty() && prefixRecords.back().k == targetSize)
			prefixRecords.back().validationInfluence = infVldt;
	}
	else
	{
		infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
		for (auto& record : prefixRecords)
			record.validationInfluence = __hyperGVldt.self_inf_cal(record.seeds);
	}
	const auto degVldt = infVldt * __numRRsets / candidateN;
	auto upperBound = infSelf / approx;
	if (mode == 1) upperBound = __boundLast;
	else if (mode == 2) upperBound = __boundMin;
	const auto upperDegOPT = upperBound * __numRRsets / candidateN;
	const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
	const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
	const auto approxOPIM = lowerSelect / upperOPT;

	for (auto& record : prefixRecords)
	{
		record.approximationForKMax = approxOPIM;
		record.runningTimeSec = timerOPIM.get_total_time();
	}
	write_prefix_outputs(outFolder, outFileName, prefixRecords, __numRRsets, __numRRsets * 2);

	__tRes.set_approximation(approxOPIM);
	__tRes.set_running_time(timerOPIM.get_total_time());
	__tRes.set_influence(infVldt);
	__tRes.set_influence_original(infSelf);
	__tRes.set_seed_vec(__vecSeed);
	__tRes.set_RR_sets_size(__numRRsets);
	std::cout << "==>OPIM k-list approx. (kmax): " << approxOPIM << " (" << infSelf / upperBound << ")\n";
	std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
	print_memory_snapshot("opim_klist_before_release_rr");
	__hyperGVldt.release_memory();
	__hyperG.release_rr_index(false);
	print_memory_snapshot("opim_klist_after_release_rr");
	LogProgress("opim_klist_done kmax=" + std::to_string(targetSize)
	            + " csv=" + outFolder + "/" + outFileName + "_prefix.csv");
	return approxOPIM;
}

double Alg::opimc(const int targetSize, const double epsilon, const double delta, const int mode,
                  const double evalDelta, const double evalEps)
{
	LogProgress("opimc_start k=" + std::to_string(targetSize)
	            + " eps=" + std::to_string(epsilon)
	            + " delta=" + std::to_string(delta)
	            + " fixed_total_rr=" + std::to_string(__fixedRRsetsTotal)
	            + " candidate_pool_cap=" + std::to_string(__statefulCandidatePoolCap));
	Timer timerOPIMC("OPIM-C");
	if (targetSize <= 0)
	{
		std::cout << "Error: OPIM-C targetSize must be positive. Handle k=0 as influence 0 outside OPIM.\n";
		exit(1);
	}
	const double e = exp(1);
	const double approx = 1 - 1.0 / e;
	const double alpha = sqrt(log(6.0 / delta));
	const auto candidateN = __hyperG.get_origin_nodes();
	if (static_cast<size_t>(targetSize) > candidateN)
	{
		std::cout << "Error: targetSize exceeds origin candidate count: "
		          << targetSize << " > " << candidateN << std::endl;
		exit(1);
	}
	const double beta = sqrt((1 - 1 / e) * (logcnk(candidateN, targetSize) + log(6.0 / delta)));
	const double rrCoeff = pow2((1 - 1 / e) * alpha + beta);
	const auto numRbase = checked_size_from_double(2.0 * rrCoeff, "numRbase");
	const auto maxNumR = checked_size_from_double(
		2.0 * static_cast<double>(candidateN) * rrCoeff / targetSize / pow2(epsilon) + 1.0,
		"maxNumR");
	if (numRbase == 0)
	{
		std::cout << "Error: numRbase is zero" << std::endl;
		exit(1);
	}
	const auto numIter = maxNumR <= numRbase ? size_t(1) : static_cast<size_t>(log2(maxNumR / numRbase)) + 1;
	const double a1 = log(numIter * 3.0 / delta);
	const double a2 = log(numIter * 3.0 / delta);
	double time1 = 0.0, time2 = 0.0;
	if (__fixedRRsetsTotal > 0)
	{
		const auto requestedTotalRR = __fixedRRsetsTotal;
		const auto numR = std::max<size_t>(1, requestedTotalRR / 2);
		LogProgress("opimc_fixed_rr_override requested_total_rr=" + std::to_string(requestedTotalRR)
		            + " per_graph_rr=" + std::to_string(numR)
		            + " actual_total_rr=" + std::to_string(numR * 2));
		timerOPIMC.get_operation_time();
		__hyperG.build_n_RRsets(numR);
		__hyperGVldt.build_n_RRsets(numR);
		__numRRsets = __hyperG.get_RR_sets_size();
		time1 += timerOPIMC.get_operation_time();
		const auto infSelf = max_cover(targetSize, mode);
		time2 += timerOPIMC.get_operation_time();
		const auto infVldt = __hyperG.uses_stateful_seed_selection_model()
			? stateful_validation_influence_fresh(__vecSeed, __numRRsets)
			: __hyperGVldt.self_inf_cal(__vecSeed);
		const auto degVldt = infVldt * __numRRsets / candidateN;
		auto upperBound = infSelf / approx;
		if (mode == 1) upperBound = __boundLast;
		else if (mode == 2) upperBound = __boundMin;
		const auto upperDegOPT = upperBound * __numRRsets / candidateN;
		const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
		const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
		const auto approxOPIMC = lowerSelect / upperOPT;
			std::cout << "==>OPIM-C fixed RR approx. (max-cover): " << approxOPIMC
			          << " (" << infSelf / upperBound << "), #RR sets (per graph): " << __numRRsets
			          << ", total budget: " << __numRRsets * 2 << '\n';
		std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
		__tRes.set_approximation(approxOPIMC);
		__tRes.set_running_time(timerOPIMC.get_total_time());
		__tRes.set_influence(infVldt);
		__tRes.set_influence_original(infSelf);
		__tRes.set_seed_vec(__vecSeed);
		__tRes.set_RR_sets_size(__numRRsets);
		print_memory_snapshot("opimc_fixed_before_final_eval_release_rr");
		__hyperGVldt.release_memory();
		__hyperG.release_rr_index(false);
		print_memory_snapshot("opimc_fixed_after_final_eval_release_rr");
		if (!__hyperG.uses_stateful_seed_selection_model())
			__tRes.set_influence(effic_inf_valid_algo(evalDelta, evalEps));
		LogProgress("opimc_fixed_rr_done k=" + std::to_string(targetSize)
		            + " influence=" + std::to_string(__tRes.get_influence())
		            + " rrsets=" + std::to_string(__tRes.get_RRsets_size())
		            + " approx=" + std::to_string(approxOPIMC));
		return __tRes.get_influence();
	}
	size_t numR = numRbase;
	for (size_t idx = 0; idx < numIter; idx++)
	{
		LogProgress("opimc_iter_start iter=" + std::to_string(idx + 1)
		            + "/" + std::to_string(numIter)
		            + " numR=" + std::to_string(numR)
		            + " maxNumR=" + std::to_string(maxNumR));
		timerOPIMC.get_operation_time();
		//dsfmt_gv_init_gen_rand(idx);
		__hyperG.build_n_RRsets(numR); // R1
		//dsfmt_gv_init_gen_rand(idx + 1000000);
		__hyperGVldt.build_n_RRsets(numR); // R2
		__numRRsets = __hyperG.get_RR_sets_size();
		time1 += timerOPIMC.get_operation_time();
		const auto infSelf = max_cover(targetSize, mode);
		time2 += timerOPIMC.get_operation_time();
			const auto infVldt = __hyperG.uses_stateful_seed_selection_model()
				? stateful_validation_influence_fresh(__vecSeed, __numRRsets)
				: __hyperGVldt.self_inf_cal(__vecSeed);
		const auto degVldt = infVldt * __numRRsets / candidateN;
		auto upperBound = infSelf / approx;
		if (mode == 1) upperBound = __boundLast;
		else if (mode == 2) upperBound = __boundMin;
		const auto upperDegOPT = upperBound * __numRRsets / candidateN;
		const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
		const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
		const auto approxOPIMC = lowerSelect / upperOPT;
		//guaranteeOpt = opt_error(degVldt, degSelfUpper, delta);
		std::cout << " -->OPIM-C (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << approxOPIMC <<
			" (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
		LogProgress("opimc_iter_result iter=" + std::to_string(idx + 1)
		            + "/" + std::to_string(numIter)
		            + " approx=" + std::to_string(approxOPIMC)
		            + " inf_self=" + std::to_string(infSelf)
		            + " inf_r2=" + std::to_string(infVldt)
		            + " rrsets=" + std::to_string(__numRRsets));
		// Check whether the requirement is satisfied
		if (approxOPIMC >= approx - epsilon)
		{
			__tRes.set_approximation(approxOPIMC);
			__tRes.set_running_time(timerOPIMC.get_total_time());
			__tRes.set_influence(infVldt);
			__tRes.set_influence_original(infSelf);
			__tRes.set_seed_vec(__vecSeed);
			__tRes.set_RR_sets_size(__numRRsets);
			std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
			std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
				print_memory_snapshot("opimc_before_final_eval_release_rr");
				__hyperGVldt.release_memory();
				__hyperG.release_rr_index(false);
				print_memory_snapshot("opimc_after_final_eval_release_rr");
				if (!__hyperG.uses_stateful_seed_selection_model())
					__tRes.set_influence(effic_inf_valid_algo(evalDelta, evalEps));
				LogProgress("opimc_done k=" + std::to_string(targetSize)
				            + " influence=" + std::to_string(__tRes.get_influence())
			            + " rrsets=" + std::to_string(__tRes.get_RRsets_size()));
			return __tRes.get_influence();
		}
		if (idx + 1 < numIter)
		{
			if (numR > std::numeric_limits<size_t>::max() / 2)
			{
				std::cout << "Error: RR-set count overflow while doubling: " << numR << std::endl;
				exit(1);
			}
			numR *= 2;
			if (numR > maxNumR) numR = maxNumR;
		}
	}
	return 0.0;
}

double Alg::opimc_klist(const std::vector<int>& kList, const double epsilon, const double delta, const int mode,
                        const double evalDelta, const double evalEps, const std::string& outFileName,
                        const std::string& outFolder)
{
	(void)evalDelta;
	(void)evalEps;
	if (kList.empty() || kList.back() <= 0)
	{
		std::cout << "Error: OPIM-C k-list needs at least one positive k.\n";
		exit(1);
	}
	const int targetSize = kList.back();
	LogProgress("opimc_klist_start kmax=" + std::to_string(targetSize)
	            + " eps=" + std::to_string(epsilon)
	            + " delta=" + std::to_string(delta)
	            + " fixed_total_rr=" + std::to_string(__fixedRRsetsTotal)
	            + " candidate_pool_cap=" + std::to_string(__statefulCandidatePoolCap));
	Timer timerOPIMC("OPIM-C-klist");
	const double e = exp(1);
	const double approx = 1 - 1.0 / e;
	const double alpha = sqrt(log(6.0 / delta));
	const auto candidateN = __hyperG.get_origin_nodes();
	if (static_cast<size_t>(targetSize) > candidateN)
	{
		std::cout << "Error: targetSize exceeds origin candidate count: "
		          << targetSize << " > " << candidateN << std::endl;
		exit(1);
	}
	const double beta = sqrt((1 - 1 / e) * (logcnk(candidateN, targetSize) + log(6.0 / delta)));
	const double rrCoeff = pow2((1 - 1 / e) * alpha + beta);
	const auto numRbase = checked_size_from_double(2.0 * rrCoeff, "numRbase");
	const auto maxNumR = checked_size_from_double(
		2.0 * static_cast<double>(candidateN) * rrCoeff / targetSize / pow2(epsilon) + 1.0,
		"maxNumR");
	if (numRbase == 0)
	{
		std::cout << "Error: numRbase is zero" << std::endl;
		exit(1);
	}
	const auto numIter = maxNumR <= numRbase ? size_t(1) : static_cast<size_t>(log2(maxNumR / numRbase)) + 1;
	const double a1 = log(numIter * 3.0 / delta);
	const double a2 = log(numIter * 3.0 / delta);
	double time1 = 0.0, time2 = 0.0;
	if (__fixedRRsetsTotal > 0)
	{
		const auto requestedTotalRR = __fixedRRsetsTotal;
		const auto numR = std::max<size_t>(1, requestedTotalRR / 2);
		LogProgress("opimc_klist_fixed_rr_override requested_total_rr=" + std::to_string(requestedTotalRR)
		            + " per_graph_rr=" + std::to_string(numR)
		            + " actual_total_rr=" + std::to_string(numR * 2));
			timerOPIMC.get_operation_time();
			__hyperG.build_n_RRsets(numR);
			__hyperGVldt.build_n_RRsets(numR);
			__numRRsets = __hyperG.get_RR_sets_size();
			time1 += timerOPIMC.get_operation_time();
			std::vector<PrefixRecord> prefixRecords;
			const auto infSelf = max_cover_lazy_prefix(targetSize, mode, kList, prefixRecords);
			time2 += timerOPIMC.get_operation_time();
			double infVldt = 0.0;
			if (__hyperG.uses_stateful_seed_selection_model())
			{
				const auto validationRRsets = choose_stateful_validation_rrsets_per_graph(__numRRsets);
				THyperGraph validator(__hyperG._graph);
				validator.set_cascade_model(__hyperG._cascadeModel);
				validator.set_origin_node_count(__hyperG.get_origin_nodes());
				validator.build_n_RRsets(validationRRsets);
				LogProgress("stateful_validation_prefix_rr"
				            " select_rr=" + std::to_string(__numRRsets)
				            + " validate_rr=" + std::to_string(validator.get_RR_sets_size())
				            + " prefixes=" + std::to_string(prefixRecords.size()));
				infVldt = validator.self_inf_cal(__vecSeed);
				for (auto& record : prefixRecords)
					record.validationInfluence = validator.self_inf_cal(record.seeds);
				if (!prefixRecords.empty() && prefixRecords.back().k == targetSize)
					prefixRecords.back().validationInfluence = infVldt;
			}
			else
			{
				infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
				for (auto& record : prefixRecords)
					record.validationInfluence = __hyperGVldt.self_inf_cal(record.seeds);
			}
			const auto degVldt = infVldt * __numRRsets / candidateN;
		auto upperBound = infSelf / approx;
		if (mode == 1) upperBound = __boundLast;
		else if (mode == 2) upperBound = __boundMin;
		const auto upperDegOPT = upperBound * __numRRsets / candidateN;
		const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
		const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
		const auto approxOPIMC = lowerSelect / upperOPT;
			for (auto& record : prefixRecords)
			{
				record.approximationForKMax = approxOPIMC;
				record.runningTimeSec = timerOPIMC.get_total_time();
			}
		write_prefix_outputs(outFolder, outFileName, prefixRecords, __numRRsets, __numRRsets * 2);
		__tRes.set_approximation(approxOPIMC);
		__tRes.set_running_time(timerOPIMC.get_total_time());
		__tRes.set_influence(infVldt);
		__tRes.set_influence_original(infSelf);
		__tRes.set_seed_vec(__vecSeed);
		__tRes.set_RR_sets_size(__numRRsets);
		std::cout << "==>Influence via R2 at kmax (fixed RR): " << infVldt
		          << ", time: " << __tRes.get_running_time() << '\n';
		std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
		print_memory_snapshot("opimc_klist_fixed_before_release_rr");
		__hyperGVldt.release_memory();
		__hyperG.release_rr_index(false);
		print_memory_snapshot("opimc_klist_fixed_after_release_rr");
		LogProgress("opimc_klist_fixed_rr_done kmax=" + std::to_string(targetSize)
		            + " csv=" + outFolder + "/" + outFileName + "_prefix.csv"
		            + " approx=" + std::to_string(approxOPIMC));
		return __tRes.get_influence();
	}
	size_t numR = numRbase;
	for (size_t idx = 0; idx < numIter; idx++)
	{
		LogProgress("opimc_klist_iter_start iter=" + std::to_string(idx + 1)
		            + "/" + std::to_string(numIter)
		            + " numR=" + std::to_string(numR)
		            + " maxNumR=" + std::to_string(maxNumR));
		timerOPIMC.get_operation_time();
		__hyperG.build_n_RRsets(numR);
		__hyperGVldt.build_n_RRsets(numR);
		__numRRsets = __hyperG.get_RR_sets_size();
		time1 += timerOPIMC.get_operation_time();
		std::vector<PrefixRecord> prefixRecords;
		const auto infSelf = max_cover_lazy_prefix(targetSize, mode, kList, prefixRecords);
		time2 += timerOPIMC.get_operation_time();
			const auto infVldt = __hyperG.uses_stateful_seed_selection_model()
				? stateful_validation_influence_fresh(__vecSeed, __numRRsets)
				: __hyperGVldt.self_inf_cal(__vecSeed);
		const auto degVldt = infVldt * __numRRsets / candidateN;
		auto upperBound = infSelf / approx;
		if (mode == 1) upperBound = __boundLast;
		else if (mode == 2) upperBound = __boundMin;
		const auto upperDegOPT = upperBound * __numRRsets / candidateN;
		const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
		const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
		const auto approxOPIMC = lowerSelect / upperOPT;
		std::cout << " -->OPIM-C-klist (" << idx + 1 << "/" << numIter << ") approx. (kmax): " << approxOPIMC <<
			" (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
		LogProgress("opimc_klist_iter_result iter=" + std::to_string(idx + 1)
		            + "/" + std::to_string(numIter)
		            + " approx=" + std::to_string(approxOPIMC)
		            + " inf_self_kmax=" + std::to_string(infSelf)
		            + " inf_r2_kmax=" + std::to_string(infVldt)
		            + " rrsets=" + std::to_string(__numRRsets));
		if (approxOPIMC >= approx - epsilon)
		{
				for (auto& record : prefixRecords)
				{
					record.approximationForKMax = approxOPIMC;
					record.runningTimeSec = timerOPIMC.get_total_time();
				}
			write_prefix_outputs(outFolder, outFileName, prefixRecords, __numRRsets, __numRRsets * 2);

			__tRes.set_approximation(approxOPIMC);
			__tRes.set_running_time(timerOPIMC.get_total_time());
			__tRes.set_influence(infVldt);
			__tRes.set_influence_original(infSelf);
			__tRes.set_seed_vec(__vecSeed);
			__tRes.set_RR_sets_size(__numRRsets);
			std::cout << "==>Influence via R2 at kmax: " << infVldt
			          << ", time: " << __tRes.get_running_time() << '\n';
			std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
			print_memory_snapshot("opimc_klist_before_release_rr");
			__hyperGVldt.release_memory();
			__hyperG.release_rr_index(false);
			print_memory_snapshot("opimc_klist_after_release_rr");
			LogProgress("opimc_klist_done kmax=" + std::to_string(targetSize)
			            + " csv=" + outFolder + "/" + outFileName + "_prefix.csv");
			return __tRes.get_influence();
		}
		if (idx + 1 < numIter)
		{
			if (numR > std::numeric_limits<size_t>::max() / 2)
			{
				std::cout << "Error: RR-set count overflow while doubling: " << numR << std::endl;
				exit(1);
			}
			numR *= 2;
			if (numR > maxNumR) numR = maxNumR;
		}
	}
	return 0.0;
}
