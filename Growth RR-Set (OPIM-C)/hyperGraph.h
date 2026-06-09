#pragma once

class HyperGraph
{
private:
	/// __numV: number of nodes in the graph.
	uint32_t __numV = 0;
	/// __numE: number of edges in the graph.
	size_t __numE = 0;
	/// __numRRsets: number of RR sets.
	size_t __numRRsets = 0;
		/// __originN: number of origin nodes used for root sampling and influence counting.
		uint32_t __originN = 0;
		std::vector<bool> __vecVisitBool;
		Nodelist __vecVisitNode;
		std::vector<int> __statefulFactorOwner;
		std::vector<int> __statefulFactorIndex;
		std::vector<int> __statefulOriginLocal;
		int __statefulFactorCount = 0;

		struct StatefulUse
		{
			int originLocal = -1;
			int factorIndex = -1;
		};

		struct StatefulSampleRef
		{
			int sampleId = -1;
			int localIndex = -1;
		};

		struct StatefulSample
		{
			bool valid = false;
			int rootLocal = -1;
			std::vector<int> origins;
			std::vector<std::vector<StatefulUse>> uses;
		};

		struct StatefulSampleState
		{
			bool covered = false;
			std::vector<uint8_t> activeOrigin;
			std::vector<uint8_t> factorSatisfied;
		};

		struct StatefulRollbackLog
		{
			bool coveredChanged = false;
			bool coveredBefore = false;
			std::vector<int> activeOrigins;
			std::vector<int> factorSatisfiedIndices;
		};

		struct StatefulGain
		{
			size_t covered = 0;
			size_t partial = 0;
		};

		struct StatefulCandidateDiagnostic
		{
			int candidate = -1;
			StatefulGain gain;
		};

		std::vector<std::vector<StatefulSampleRef>> _statefulHyperG;
		std::vector<StatefulSample> _statefulSamples;

	/// Initialization
	void init_hypergraph()
	{
		if (_graph.size() > std::numeric_limits<uint32_t>::max())
		{
			std::cout << "Error: graph has too many nodes for OPIM uint32_t node ids: "
			          << _graph.size() << std::endl;
			exit(1);
		}
		__numV = static_cast<uint32_t>(_graph.size());
		__originN = __numV;
		for (const auto& nbrs : _graph) __numE += nbrs.size();
		_FRsets = FRsets(__originN);
		__vecVisitBool = std::vector<bool>(__numV);
		__vecVisitNode = Nodelist(__numV);
	}

	bool is_origin_node(const uint32_t node) const
	{
		return node < __originN;
	}

	bool is_seed_node(const std::vector<bool>& vecBoolSeed, const uint32_t node) const
	{
		return is_origin_node(node) && vecBoolSeed[node];
	}

		bool uses_stateful_factor_gate() const
		{
			return _cascadeModel == IC2 || _cascadeModel == IC2MIX;
		}

		bool uses_stateful_seed_selection() const
		{
			return _cascadeModel == IC2MIX;
		}

		bool has_only_factor_parents(const uint32_t node) const
		{
		if (!is_origin_node(node))
			return false;
		const auto& nbrs = _graph[node];
		if (nbrs.empty())
			return false;
		for (const auto& nbr : nbrs)
		{
			if (is_origin_node(nbr.first))
				return false;
			}
			return true;
		}

		void clear_stateful_metadata()
		{
			std::vector<int>().swap(__statefulFactorOwner);
			std::vector<int>().swap(__statefulFactorIndex);
			std::vector<int>().swap(__statefulOriginLocal);
			__statefulFactorCount = 0;
		}

		void ensure_stateful_factor_metadata()
		{
			if (!uses_stateful_seed_selection())
				return;
			if (__statefulFactorOwner.size() == __numV && __statefulFactorIndex.size() == __numV && __statefulFactorCount > 0)
				return;

			__statefulFactorOwner.assign(__numV, -1);
			__statefulFactorIndex.assign(__numV, -1);
			__statefulFactorCount = 0;
			for (uint32_t origin = 0; origin < __originN; ++origin)
			{
				int factorIndex = 0;
				for (const auto& nbr : _graph[origin])
				{
					if (nbr.first < __originN || nbr.first >= __numV)
						continue;
					__statefulFactorOwner[nbr.first] = static_cast<int>(origin);
					__statefulFactorIndex[nbr.first] = factorIndex++;
				}
				if (factorIndex > __statefulFactorCount)
					__statefulFactorCount = factorIndex;
			}
			if (__statefulFactorCount <= 0)
				__statefulFactorCount = 1;
		}

		void ensure_stateful_workspace()
		{
			if (_statefulHyperG.size() != __originN)
				_statefulHyperG.assign(__originN, std::vector<StatefulSampleRef>());
			if (__statefulOriginLocal.size() != __originN)
				__statefulOriginLocal.assign(__originN, -1);
		}

	void ensure_origin_buckets()
	{
		if (_FRsets.size() != __originN)
		{
			_FRsets = FRsets(__originN);
		}
	}

	void add_to_cover_index_if_origin(const uint32_t node, const size_t hyperIdx)
	{
		if (is_origin_node(node))
		{
			_FRsets[node].push_back(hyperIdx);
		}
	}

	void append_origin_nodes_to_rrset(const size_t numVisitNode)
	{
		RRset rrset;
		rrset.reserve(numVisitNode);
		for (size_t i = 0; i < numVisitNode; i++)
		{
			const auto node = __vecVisitNode[i];
			if (is_origin_node(node))
			{
				rrset.push_back(node);
			}
			__vecVisitBool[node] = false;
		}
		_RRsets.push_back(std::move(rrset));
	}

		void clear_visit_marks(const size_t numVisitNode)
		{
			for (size_t i = 0; i < numVisitNode; i++)
			{
				__vecVisitBool[__vecVisitNode[i]] = false;
			}
		}

		int add_stateful_origin(StatefulSample& sample, const uint32_t origin)
		{
			int& local = __statefulOriginLocal[origin];
			if (local >= 0)
				return local;
			local = static_cast<int>(sample.origins.size());
			sample.origins.push_back(static_cast<int>(origin));
			sample.uses.emplace_back();
			return local;
		}

		void append_stateful_use(StatefulSample& sample, std::vector<uint8_t>& expanded,
		                         std::deque<int>& originQueue, const uint32_t parentOrigin,
		                         const int currentLocal, const int factorIndex)
		{
			const auto parentLocal = add_stateful_origin(sample, parentOrigin);
			StatefulUse use;
			use.originLocal = currentLocal;
			use.factorIndex = factorIndex;
			sample.uses[static_cast<size_t>(parentLocal)].push_back(use);
			if (parentLocal >= static_cast<int>(expanded.size()))
				expanded.resize(static_cast<size_t>(parentLocal) + 1, 0);
			if (!expanded[static_cast<size_t>(parentLocal)])
				originQueue.push_back(parentLocal);
		}

		void reset_stateful_origin_local(const StatefulSample& sample)
		{
			for (const auto origin : sample.origins)
			{
				if (origin >= 0 && static_cast<size_t>(origin) < __statefulOriginLocal.size())
					__statefulOriginLocal[static_cast<size_t>(origin)] = -1;
			}
		}

		void build_one_stateful_rrset(const uint32_t uStart, const size_t hyperIdx)
		{
			ensure_stateful_factor_metadata();
			ensure_stateful_workspace();
			if (_statefulSamples.size() <= hyperIdx)
				_statefulSamples.resize(hyperIdx + 1);

			StatefulSample sample;
			sample.valid = true;
			sample.rootLocal = add_stateful_origin(sample, uStart);

			std::deque<int> originQueue;
			std::vector<uint8_t> expanded(1, 0);
			originQueue.push_back(sample.rootLocal);

			while (!originQueue.empty())
			{
				const auto currentLocal = originQueue.front();
				originQueue.pop_front();
				if (currentLocal < 0 || currentLocal >= static_cast<int>(sample.origins.size()))
					continue;
				if (currentLocal >= static_cast<int>(expanded.size()))
					expanded.resize(static_cast<size_t>(currentLocal) + 1, 0);
				if (expanded[static_cast<size_t>(currentLocal)])
					continue;
				expanded[static_cast<size_t>(currentLocal)] = 1;

				const auto origin = static_cast<uint32_t>(sample.origins[static_cast<size_t>(currentLocal)]);
				if (_graph[origin].empty())
					continue;

				bool allFactorParents = true;
				for (const auto& nbr : _graph[origin])
				{
					if (nbr.first < __originN || nbr.first >= __numV)
					{
						allFactorParents = false;
						break;
					}
				}

				if (!allFactorParents)
				{
					for (const auto& nbr : _graph[origin])
					{
						const auto parentOrigin = nbr.first;
						if (parentOrigin >= __originN)
							continue;
						if (dsfmt_gv_genrand_open_close() > nbr.second)
							continue;
						append_stateful_use(sample, expanded, originQueue, parentOrigin, currentLocal, -1);
					}
					continue;
				}

				bool allSuccess = true;
				for (const auto& nbr : _graph[origin])
				{
					if (dsfmt_gv_genrand_open_close() > nbr.second)
					{
						allSuccess = false;
						break;
					}
				}
				if (!allSuccess)
					continue;

				for (const auto& factorNbr : _graph[origin])
				{
					const auto factorNode = factorNbr.first;
					if (factorNode < __originN || factorNode >= __numV)
						continue;
					const auto factorIndex = __statefulFactorIndex[factorNode];
					if (factorIndex < 0 || factorIndex >= __statefulFactorCount)
						continue;
					for (const auto& parentNbr : _graph[factorNode])
					{
						const auto parentOrigin = parentNbr.first;
						if (parentOrigin >= __originN)
							continue;
						if (dsfmt_gv_genrand_open_close() > parentNbr.second)
							continue;
						append_stateful_use(sample, expanded, originQueue, parentOrigin, currentLocal, factorIndex);
					}
				}
			}

			reset_stateful_origin_local(sample);
			_statefulSamples[hyperIdx] = std::move(sample);
			for (size_t local = 0; local < _statefulSamples[hyperIdx].origins.size(); ++local)
			{
				const auto origin = _statefulSamples[hyperIdx].origins[local];
				if (origin < 0 || static_cast<uint32_t>(origin) >= __originN)
					continue;
				_statefulHyperG[static_cast<size_t>(origin)].push_back(
					StatefulSampleRef{static_cast<int>(hyperIdx), static_cast<int>(local)});
			}
		}

		bool all_stateful_factors_satisfied(const StatefulSampleState& state, const int local) const
		{
			for (int factor = 0; factor < __statefulFactorCount; ++factor)
			{
				const auto index = static_cast<size_t>(local) * static_cast<size_t>(__statefulFactorCount)
				                 + static_cast<size_t>(factor);
				if (!state.factorSatisfied[index])
					return false;
			}
			return true;
		}

		bool stateful_root_is_covered(const StatefulSample& sample, const StatefulSampleState& state) const
		{
			return sample.rootLocal >= 0
			    && sample.rootLocal < static_cast<int>(state.activeOrigin.size())
			    && state.activeOrigin[static_cast<size_t>(sample.rootLocal)];
		}

		void rollback_stateful_trial(StatefulSampleState& state, const StatefulRollbackLog& rollback) const
		{
			for (int i = static_cast<int>(rollback.factorSatisfiedIndices.size()) - 1; i >= 0; --i)
			{
				const auto index = rollback.factorSatisfiedIndices[static_cast<size_t>(i)];
				if (index >= 0 && index < static_cast<int>(state.factorSatisfied.size()))
					state.factorSatisfied[static_cast<size_t>(index)] = 0;
			}
			for (int i = static_cast<int>(rollback.activeOrigins.size()) - 1; i >= 0; --i)
			{
				const auto local = rollback.activeOrigins[static_cast<size_t>(i)];
				if (local >= 0 && local < static_cast<int>(state.activeOrigin.size()))
					state.activeOrigin[static_cast<size_t>(local)] = 0;
			}
			if (rollback.coveredChanged)
				state.covered = rollback.coveredBefore;
		}

		StatefulGain propagate_stateful_seed(const StatefulSample& sample, StatefulSampleState& state, const int seedLocal,
		                                     StatefulRollbackLog* rollback = nullptr) const
		{
			StatefulGain gain;
			if (state.covered || !sample.valid)
				return gain;
			if (seedLocal < 0 || seedLocal >= static_cast<int>(sample.origins.size()))
				return gain;

			const auto wasCovered = state.covered;
			std::deque<int> work;
			if (!state.activeOrigin[static_cast<size_t>(seedLocal)])
			{
				if (rollback)
					rollback->activeOrigins.push_back(seedLocal);
				state.activeOrigin[static_cast<size_t>(seedLocal)] = 1;
				gain.partial++;
				work.push_back(seedLocal);
			}

			while (!work.empty())
			{
				const auto activeLocal = work.front();
				work.pop_front();
				for (const auto& use : sample.uses[static_cast<size_t>(activeLocal)])
				{
					const auto dst = use.originLocal;
					if (dst < 0 || dst >= static_cast<int>(sample.origins.size()))
						continue;
					if (use.factorIndex < 0)
					{
						if (!state.activeOrigin[static_cast<size_t>(dst)])
						{
							if (rollback)
								rollback->activeOrigins.push_back(dst);
							state.activeOrigin[static_cast<size_t>(dst)] = 1;
							gain.partial++;
							work.push_back(dst);
						}
						continue;
					}

					const auto satIndex = static_cast<size_t>(dst) * static_cast<size_t>(__statefulFactorCount)
					                    + static_cast<size_t>(use.factorIndex);
					if (state.factorSatisfied[satIndex])
						continue;
					if (rollback)
						rollback->factorSatisfiedIndices.push_back(static_cast<int>(satIndex));
					state.factorSatisfied[satIndex] = 1;
					gain.partial++;
					if (!state.activeOrigin[static_cast<size_t>(dst)] && all_stateful_factors_satisfied(state, dst))
					{
						if (rollback)
							rollback->activeOrigins.push_back(dst);
						state.activeOrigin[static_cast<size_t>(dst)] = 1;
						gain.partial++;
						work.push_back(dst);
					}
				}
			}

			if (!wasCovered && stateful_root_is_covered(sample, state))
			{
				if (rollback)
				{
					rollback->coveredChanged = true;
					rollback->coveredBefore = wasCovered;
				}
				state.covered = true;
				gain.covered = 1;
			}
			return gain;
		}

		std::vector<StatefulSampleState> make_stateful_states() const
		{
			std::vector<StatefulSampleState> states(_statefulSamples.size());
			for (size_t i = 0; i < _statefulSamples.size(); ++i)
			{
				if (!_statefulSamples[i].valid)
					continue;
				const auto localCount = _statefulSamples[i].origins.size();
				states[i].activeOrigin.assign(localCount, 0);
				states[i].factorSatisfied.assign(localCount * static_cast<size_t>(__statefulFactorCount), 0);
			}
			return states;
		}

		size_t stateful_sample_count() const
		{
			size_t total = 0;
			for (const auto& sample : _statefulSamples)
			{
				if (sample.valid)
					++total;
			}
			return total;
		}

		StatefulCandidateDiagnostic evaluate_stateful_candidate(const int candidate,
		                                                        std::vector<StatefulSampleState>& states) const
		{
			StatefulCandidateDiagnostic total;
			total.candidate = candidate;
			if (candidate < 0 || static_cast<uint32_t>(candidate) >= __originN)
				return total;
			if (static_cast<size_t>(candidate) >= _statefulHyperG.size())
				return total;
			for (const auto& ref : _statefulHyperG[static_cast<size_t>(candidate)])
			{
				if (ref.sampleId < 0 || static_cast<size_t>(ref.sampleId) >= _statefulSamples.size())
					continue;
				const auto& sample = _statefulSamples[static_cast<size_t>(ref.sampleId)];
				auto& state = states[static_cast<size_t>(ref.sampleId)];
				if (!sample.valid || state.covered)
					continue;
				StatefulRollbackLog rollback;
				const auto gain = propagate_stateful_seed(sample, state, ref.localIndex, &rollback);
				rollback_stateful_trial(state, rollback);
				total.gain.covered += gain.covered;
				total.gain.partial += gain.partial;
			}
			return total;
		}

		StatefulGain commit_stateful_candidate(const int candidate, std::vector<StatefulSampleState>& states) const
		{
			StatefulGain total;
			if (candidate < 0 || static_cast<uint32_t>(candidate) >= __originN)
				return total;
			if (static_cast<size_t>(candidate) >= _statefulHyperG.size())
				return total;
			for (const auto& ref : _statefulHyperG[static_cast<size_t>(candidate)])
			{
				if (ref.sampleId < 0 || static_cast<size_t>(ref.sampleId) >= _statefulSamples.size())
					continue;
				const auto& sample = _statefulSamples[static_cast<size_t>(ref.sampleId)];
				auto& state = states[static_cast<size_t>(ref.sampleId)];
				if (!sample.valid || state.covered)
					continue;
				const auto gain = propagate_stateful_seed(sample, state, ref.localIndex);
				total.covered += gain.covered;
				total.partial += gain.partial;
			}
			return total;
		}

		bool better_stateful_candidate(const StatefulCandidateDiagnostic& lhs,
		                               const StatefulCandidateDiagnostic& rhs) const
		{
			if (rhs.candidate < 0)
				return true;
			if (lhs.gain.covered != rhs.gain.covered)
				return lhs.gain.covered > rhs.gain.covered;
			if (lhs.gain.partial != rhs.gain.partial)
				return lhs.gain.partial > rhs.gain.partial;
			return lhs.candidate < rhs.candidate;
		}

public:
	/// _graph: reverse graph
	const Graph& _graph;
	/// _FRsets: origin-only forward cover sets. _FRsets[i] contains RR-set ids covered by origin node i.
	FRsets _FRsets;
	/// _RRsets: origin-only reverse reachable sets. Factor nodes are traversed but not stored here.
	RRsets _RRsets;
	/// _cascadeModel: the cascade model, default is IC
	CascadeModel _cascadeModel = IC;

	explicit HyperGraph(const Graph& graph) : _graph(graph)
	{
		init_hypergraph();
	}

	/// Set cascade model
		void set_cascade_model(const CascadeModel model)
		{
			if (_cascadeModel != model)
			{
				clear_stateful_metadata();
				release_rr_index(true);
			}
			_cascadeModel = model;
		}

	/// Restrict RR roots and influence estimation to origin nodes.
	void set_origin_node_count(const uint32_t originN)
	{
		const auto nextOriginN = (originN == 0 || originN > __numV) ? __numV : originN;
		if (nextOriginN != __originN)
		{
			__originN = nextOriginN;
			refresh_hypergraph();
		}
		else
		{
			ensure_origin_buckets();
		}
	}

	/// Returns the number of origin nodes used by origin-only mode.
	uint32_t get_origin_nodes() const
	{
		return __originN;
	}

	/// Returns the number of nodes in the graph.
	uint32_t get_nodes() const
	{
		return __numV;
	}

	/// Returns the number of edges in the graph.
	size_t get_edges() const
	{
		return __numE;
	}

	/// Returns the number of RR sets in the graph.
		size_t get_RR_sets_size() const
		{
			return __numRRsets;
		}

		bool uses_stateful_seed_selection_model() const
		{
			return uses_stateful_seed_selection();
		}

	/// Get out degree
	std::vector<size_t> get_out_degree() const
	{
		std::vector<size_t> outDeg(__numV);
		for (const auto& nbrs : _graph)
		{
			for (const auto& nbr : nbrs)
			{
				outDeg[nbr.first]++;
			}
		}
		return outDeg;
	}

	/// Generate a set of n RR sets
	void build_n_RRsets(const size_t numSamples)
	{
		if (__originN == 0)
		{
			std::cout << "Error: no origin nodes available for RR-set sampling" << std::endl;
			exit(1);
		}
		ensure_origin_buckets();
		const auto prevSize = __numRRsets;
		__numRRsets = __numRRsets > numSamples ? __numRRsets : numSamples;
		const auto todo = numSamples > prevSize ? numSamples - prevSize : size_t(0);
		const auto progressStep = todo > 10 ? todo / 10 : size_t(1);
		LogProgress("rr_build_start target_rr=" + std::to_string(numSamples)
		            + " existing_rr=" + std::to_string(prevSize)
		            + " origin_n=" + std::to_string(__originN)
		            + " graph_nodes=" + std::to_string(__numV)
		            + " model=" + std::to_string(static_cast<int>(_cascadeModel)));
			for (auto i = prevSize; i < numSamples; i++)
			{
				const auto root = dsfmt_gv_genrand_uint32_range(__originN);
				if (uses_stateful_seed_selection())
					build_one_stateful_rrset(root, i);
				else
					build_one_RRset(root, i);
				const auto built = i + 1 - prevSize;
				if (built == todo || built % progressStep == 0)
				{
				LogProgress("rr_progress target_rr=" + std::to_string(numSamples)
				            + " built_this_call=" + std::to_string(built)
				            + "/" + std::to_string(todo)
				            + " total_rr=" + std::to_string(i + 1));
			}
		}
		LogProgress("rr_build_done target_rr=" + std::to_string(numSamples));
	}

	/// Generate one RR set
	void build_one_RRset(const uint32_t uStart, const size_t hyperIdx)
	{
		size_t numVisitNode = 0, currIdx = 0;
		add_to_cover_index_if_origin(uStart, hyperIdx);
		__vecVisitNode[numVisitNode++] = uStart;
		__vecVisitBool[uStart] = true;
		while (currIdx < numVisitNode)
		{
			const auto expand = __vecVisitNode[currIdx++];
			if (_cascadeModel == IC || _cascadeModel == IC2 || _cascadeModel == IC2MIX)
			{
				if (uses_stateful_factor_gate() && has_only_factor_parents(expand))
				{
					bool allSuccess = true;
					for (const auto& nbr : _graph[expand])
					{
						const auto randDouble = dsfmt_gv_genrand_open_close();
						if (randDouble > nbr.second)
						{
							allSuccess = false;
							break;
						}
					}
					if (!allSuccess)
						continue;
					for (const auto& nbr : _graph[expand])
					{
						const auto nbrId = nbr.first;
						if (__vecVisitBool[nbrId])
							continue;
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
						add_to_cover_index_if_origin(nbrId, hyperIdx);
					}
				}
				else
				{
					for (const auto& nbr : _graph[expand])
					{
						const auto nbrId = nbr.first;
						if (__vecVisitBool[nbrId])
							continue;
						const auto randDouble = dsfmt_gv_genrand_open_close();
						if (randDouble > nbr.second)
							continue;
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
						add_to_cover_index_if_origin(nbrId, hyperIdx);
					}
				}
			}
			else if (_cascadeModel == LT)
			{
				if (_graph[expand].empty())
					continue;
				const auto nextNbrIdx = gen_random_node_by_weight_LT(_graph[expand]);
				if (nextNbrIdx >= _graph[expand].size()) break; // No element activated
				const auto nbrId = _graph[expand][nextNbrIdx].first;
				if (__vecVisitBool[nbrId]) break; // Stop, no further node activated
				__vecVisitNode[numVisitNode++] = nbrId;
				__vecVisitBool[nbrId] = true;
				add_to_cover_index_if_origin(nbrId, hyperIdx);
			}
			else if (_cascadeModel == OLO)
			{
				if (_graph[expand].empty())
					continue;
				// Strict OLO: sample once for the current node, then keep all parents or none.
				const auto randDouble = dsfmt_gv_genrand_open_close();
				const auto listenProb = _graph[expand][0].second;
				if (randDouble > listenProb)
					continue;
				for (const auto& nbr : _graph[expand])
				{
					const auto nbrId = nbr.first;
					if (__vecVisitBool[nbrId])
						continue;
					__vecVisitNode[numVisitNode++] = nbrId;
					__vecVisitBool[nbrId] = true;
					add_to_cover_index_if_origin(nbrId, hyperIdx);
				}
			}
		}
		append_origin_nodes_to_rrset(numVisitNode);
	}

	/// Evaluate the influence spread of a seed set on current generated RR sets
		double self_inf_cal(const Nodelist& vecSeed)
		{
			if (__numRRsets == 0) return 0.0;
			if (uses_stateful_seed_selection())
			{
				const auto sampleCount = stateful_sample_count();
				if (sampleCount == 0)
					return 0.0;
				auto states = make_stateful_states();
				size_t influence = 0;
				for (auto seed : vecSeed)
				{
					if (!is_origin_node(seed))
						continue;
					influence += commit_stateful_candidate(static_cast<int>(seed), states).covered;
				}
				return 1.0 * influence * __originN / sampleCount;
			}
			std::vector<bool> vecBoolVst = std::vector<bool>(__numRRsets);
			for (auto seed : vecSeed)
			{
			if (!is_origin_node(seed)) continue;
			for (auto node : _FRsets[seed])
			{
				vecBoolVst[node] = true;
			}
			}
			return 1.0 * std::count(vecBoolVst.begin(), vecBoolVst.end(), true) * __originN / __numRRsets;
		}

		std::vector<PrefixRecord> build_seedset_prefix_stateful(const std::vector<int>& requestedK,
		                                                        const size_t candidatePoolCap) const
		{
			std::vector<PrefixRecord> result;
			if (!uses_stateful_seed_selection())
				return result;
			if (requestedK.empty())
				return result;

			std::vector<int> wanted = requestedK;
			std::sort(wanted.begin(), wanted.end());
			wanted.erase(std::unique(wanted.begin(), wanted.end()), wanted.end());
			for (auto& value : wanted)
				value = std::min<int>(value, static_cast<int>(__originN));
			wanted.erase(std::unique(wanted.begin(), wanted.end()), wanted.end());
			if (wanted.empty())
				return result;

			size_t nextWantedIndex = 0;
			while (nextWantedIndex < wanted.size() && wanted[nextWantedIndex] <= 0)
			{
				PrefixRecord record;
				record.k = wanted[nextWantedIndex++];
				result.push_back(std::move(record));
			}

			const auto sampleCount = stateful_sample_count();
			auto states = make_stateful_states();
			std::vector<uint8_t> selected(__originN, 0);
			std::vector<int> candidates;
			candidates.reserve(__originN);

			if (candidatePoolCap == 0 || candidatePoolCap >= __originN)
			{
				for (uint32_t candidate = 0; candidate < __originN; ++candidate)
				{
					if (candidate < _statefulHyperG.size() && !_statefulHyperG[candidate].empty())
						candidates.push_back(static_cast<int>(candidate));
				}
			}
			else
			{
				typedef std::pair<size_t, uint32_t> CandidateCoveragePair;
				struct CandidateMinHeapCompare
				{
					bool operator()(const CandidateCoveragePair& lhs, const CandidateCoveragePair& rhs) const
					{
						if (lhs.first != rhs.first) return lhs.first > rhs.first;
						return lhs.second < rhs.second;
					}
				};
				std::priority_queue<CandidateCoveragePair, std::vector<CandidateCoveragePair>, CandidateMinHeapCompare> heap;
				for (uint32_t candidate = 0; candidate < __originN; ++candidate)
				{
					if (candidate >= _statefulHyperG.size() || _statefulHyperG[candidate].empty())
						continue;
					CandidateCoveragePair item(_statefulHyperG[candidate].size(), candidate);
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
				std::vector<uint8_t> inPool(__originN, 0);
				while (!heap.empty())
				{
					inPool[heap.top().second] = 1;
					heap.pop();
				}
				const auto rootDirectCap = std::max<size_t>(5000, candidatePoolCap * 10);
				size_t rootDirectAdded = 0;
				for (const auto& sample : _statefulSamples)
				{
					if (!sample.valid || sample.rootLocal < 0 || static_cast<size_t>(sample.rootLocal) >= sample.uses.size())
						continue;
					for (size_t parentLocal = 0; parentLocal < sample.uses.size() && rootDirectAdded < rootDirectCap; ++parentLocal)
					{
						for (const auto& use : sample.uses[parentLocal])
						{
							if (use.originLocal != sample.rootLocal || use.factorIndex < 0)
								continue;
							const auto candidate = sample.origins[parentLocal];
							if (candidate < 0 || static_cast<uint32_t>(candidate) >= __originN)
								continue;
							if (candidate >= static_cast<int>(inPool.size()) || inPool[static_cast<size_t>(candidate)])
								continue;
							inPool[static_cast<size_t>(candidate)] = 1;
							++rootDirectAdded;
						}
					}
					if (rootDirectAdded >= rootDirectCap)
						break;
				}
				for (uint32_t candidate = 0; candidate < __originN; ++candidate)
				{
					if (inPool[candidate] && candidate < _statefulHyperG.size() && !_statefulHyperG[candidate].empty())
						candidates.push_back(static_cast<int>(candidate));
				}
			}

			std::sort(candidates.begin(), candidates.end());
			Nodelist seeds;
			seeds.reserve(static_cast<size_t>(wanted.back()));
			size_t influence = 0;
			while (seeds.size() < static_cast<size_t>(wanted.back()))
			{
				StatefulCandidateDiagnostic best;
				for (const auto candidate : candidates)
				{
					if (selected[static_cast<size_t>(candidate)])
						continue;
					const auto diag = evaluate_stateful_candidate(candidate, states);
					if (better_stateful_candidate(diag, best))
						best = diag;
				}
				if (best.candidate < 0)
					break;
				selected[static_cast<size_t>(best.candidate)] = 1;
				const auto committed = commit_stateful_candidate(best.candidate, states);
				influence += committed.covered;
				seeds.push_back(static_cast<uint32_t>(best.candidate));
				while (nextWantedIndex < wanted.size()
				    && wanted[nextWantedIndex] <= static_cast<int>(seeds.size()))
				{
					PrefixRecord record;
					record.k = wanted[nextWantedIndex];
					record.coveredRR = influence;
					record.selfEstimatedInfluence = sampleCount > 0 ? 1.0 * influence * __originN / sampleCount : 0.0;
					record.seeds.assign(seeds.begin(), seeds.begin() + record.k);
					result.push_back(std::move(record));
					++nextWantedIndex;
				}
			}

			while (nextWantedIndex < wanted.size())
			{
				PrefixRecord record;
				record.k = wanted[nextWantedIndex++];
				record.coveredRR = influence;
				record.selfEstimatedInfluence = sampleCount > 0 ? 1.0 * influence * __originN / sampleCount : 0.0;
				record.seeds = seeds;
				result.push_back(std::move(record));
			}
			return result;
		}

	/// Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
		double effic_inf_valid_algo(const Nodelist& vecSeed, const double delta = 1e-3, const double eps = 0.01)
		{
			if (uses_stateful_seed_selection())
			{
				return self_inf_cal(vecSeed);
			}
			const double c = 2.0 * (exp(1.0) - 2.0);
			const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
		size_t numHyperEdge = 0;
		size_t numCoverd = 0;
		std::vector<bool> vecBoolSeed(__originN);
		for (auto seed : vecSeed)
		{
			if (is_origin_node(seed)) vecBoolSeed[seed] = true;
		}

		while (numCoverd < LambdaL)
		{
			numHyperEdge++;
			size_t numVisitNode = 0, currIdx = 0;
			const auto uStart = dsfmt_gv_genrand_uint32_range(__originN);
			if (vecBoolSeed[uStart])
			{
				// Stop, this sample is covered
				numCoverd++;
				continue;
			}
			__vecVisitNode[numVisitNode++] = uStart;
			__vecVisitBool[uStart] = true;
			while (currIdx < numVisitNode)
			{
				const auto expand = __vecVisitNode[currIdx++];
				if (_cascadeModel == IC || _cascadeModel == IC2 || _cascadeModel == IC2MIX)
				{
					if (uses_stateful_factor_gate() && has_only_factor_parents(expand))
					{
						bool allSuccess = true;
						for (const auto& nbr : _graph[expand])
						{
							const auto randDouble = dsfmt_gv_genrand_open_close();
							if (randDouble > nbr.second)
							{
								allSuccess = false;
								break;
							}
						}
						if (!allSuccess)
							continue;
						for (const auto& nbr : _graph[expand])
						{
							if (is_seed_node(vecBoolSeed, nbr.first))
							{
								// Stop, this sample is covered
								numCoverd++;
								goto postProcess;
							}
						}
						for (const auto& nbr : _graph[expand])
						{
							const auto nbrId = nbr.first;
							if (__vecVisitBool[nbrId])
								continue;
							__vecVisitNode[numVisitNode++] = nbrId;
							__vecVisitBool[nbrId] = true;
						}
					}
					else
					{
						for (const auto& nbr : _graph[expand])
						{
							const auto nbrId = nbr.first;
							if (__vecVisitBool[nbrId])
								continue;
							const auto randDouble = dsfmt_gv_genrand_open_close();
							if (randDouble > nbr.second)
								continue;
							if (is_seed_node(vecBoolSeed, nbrId))
							{
								// Stop, this sample is covered
								numCoverd++;
								goto postProcess;
							}
							__vecVisitNode[numVisitNode++] = nbrId;
							__vecVisitBool[nbrId] = true;
						}
					}
				}
				else if (_cascadeModel == LT)
				{
					if (_graph[expand].empty())
						continue;
					const auto nextNbrIdx = gen_random_node_by_weight_LT(_graph[expand]);
					if (nextNbrIdx >= _graph[expand].size()) break; // No element activated
					const auto nbrId = _graph[expand][nextNbrIdx].first;
					if (__vecVisitBool[nbrId]) break; // Stop, no further node activated
					if (is_seed_node(vecBoolSeed, nbrId))
					{
						// Stop, this sample is covered
						numCoverd++;
						goto postProcess;
					}
					__vecVisitNode[numVisitNode++] = nbrId;
					__vecVisitBool[nbrId] = true;
				}
				else if (_cascadeModel == OLO)
				{
					if (_graph[expand].empty())
						continue;
					const auto randDouble = dsfmt_gv_genrand_open_close();
					const auto listenProb = _graph[expand][0].second;
					if (randDouble > listenProb)
						continue;
					for (const auto& nbr : _graph[expand])
					{
						if (is_seed_node(vecBoolSeed, nbr.first))
						{
							// Stop, this sample is covered
							numCoverd++;
							goto postProcess;
						}
					}
					for (const auto& nbr : _graph[expand])
					{
						const auto nbrId = nbr.first;
						if (__vecVisitBool[nbrId])
							continue;
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
					}
				}
			}
		postProcess:
			clear_visit_marks(numVisitNode);
		}
		return 1.0 * numCoverd * __originN / numHyperEdge;
	}

	/// Efficiently evaluate the influence spread of a seed set with a given number of RR sets to test
	double effic_inf_valid_algo_with_samplesize(const std::vector<uint32_t>& vecSeed, const size_t numSamples)
	{
		size_t numHyperEdge = 0;
		size_t numCoverd = 0;
		std::vector<bool> vecBoolSeed(__originN);
		for (auto seed : vecSeed)
		{
			if (is_origin_node(seed)) vecBoolSeed[seed] = true;
		}
		while (++numHyperEdge < numSamples)
		{
			size_t numVisitNode = 0, currIdx = 0;
			const auto uStart = dsfmt_gv_genrand_uint32_range(__originN);
			if (vecBoolSeed[uStart])
			{
				// Stop, this sample is covered
				numCoverd++;
				continue;
			}
			__vecVisitNode[numVisitNode++] = uStart;
			__vecVisitBool[uStart] = true;
			while (currIdx < numVisitNode)
			{
				const auto expand = __vecVisitNode[currIdx++];
				if (_cascadeModel == IC || _cascadeModel == IC2 || _cascadeModel == IC2MIX)
				{
					if (uses_stateful_factor_gate() && has_only_factor_parents(expand))
					{
						bool allSuccess = true;
						for (const auto& nbr : _graph[expand])
						{
							const auto randDouble = dsfmt_gv_genrand_open_close();
							if (randDouble > nbr.second)
							{
								allSuccess = false;
								break;
							}
						}
						if (!allSuccess)
							continue;
						for (const auto& nbr : _graph[expand])
						{
							if (is_seed_node(vecBoolSeed, nbr.first))
							{
								// Stop, this sample is covered
								numCoverd++;
								goto postProcess;
							}
						}
						for (const auto& nbr : _graph[expand])
						{
							const auto nbrId = nbr.first;
							if (__vecVisitBool[nbrId])
								continue;
							__vecVisitNode[numVisitNode++] = nbrId;
							__vecVisitBool[nbrId] = true;
						}
					}
					else
					{
						for (const auto& nbr : _graph[expand])
						{
							const auto nbrId = nbr.first;
							if (__vecVisitBool[nbrId])
								continue;
							const auto randDouble = dsfmt_gv_genrand_open_close();
							if (randDouble > nbr.second)
								continue;
							if (is_seed_node(vecBoolSeed, nbrId))
							{
								// Stop, this sample is covered
								numCoverd++;
								goto postProcess;
							}
							__vecVisitNode[numVisitNode++] = nbrId;
							__vecVisitBool[nbrId] = true;
						}
					}
				}
				else if (_cascadeModel == LT)
				{
					if (_graph[expand].empty())
						continue;
					const auto nextNbrIdx = gen_random_node_by_weight_LT(_graph[expand]);
					if (nextNbrIdx >= _graph[expand].size()) break; // No element activated
					const auto nbrId = _graph[expand][nextNbrIdx].first;
					if (__vecVisitBool[nbrId]) break; // Stop, no further node activated
					if (is_seed_node(vecBoolSeed, nbrId))
					{
						// Stop, this sample is covered
						numCoverd++;
						goto postProcess;
					}
					__vecVisitNode[numVisitNode++] = nbrId;
					__vecVisitBool[nbrId] = true;
				}
				else if (_cascadeModel == OLO)
				{
					if (_graph[expand].empty())
						continue;
					const auto randDouble = dsfmt_gv_genrand_open_close();
					const auto listenProb = _graph[expand][0].second;
					if (randDouble > listenProb)
						continue;
					for (const auto& nbr : _graph[expand])
					{
						if (is_seed_node(vecBoolSeed, nbr.first))
						{
							// Stop, this sample is covered
							numCoverd++;
							goto postProcess;
						}
					}
					for (const auto& nbr : _graph[expand])
					{
						const auto nbrId = nbr.first;
						if (__vecVisitBool[nbrId])
							continue;
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
					}
				}
			}
		postProcess:
			clear_visit_marks(numVisitNode);
		}
		return 1.0 * numCoverd * __originN / numHyperEdge;
	}

	/// Release RR sets and forward cover index. Keep origin buckets when more RR sampling will follow.
		void release_rr_index(const bool keepOriginBuckets = true)
		{
			LogProgress("release_rr_index_start rrsets=" + std::to_string(_RRsets.size())
			            + " fr_buckets=" + std::to_string(_FRsets.size())
			            + " keep_origin_buckets=" + std::to_string(keepOriginBuckets ? 1 : 0));
		for (auto& rrset : _RRsets)
		{
			RRset().swap(rrset);
		}
		RRsets().swap(_RRsets);
		for (auto& frset : _FRsets)
		{
			FRset().swap(frset);
		}
		if (keepOriginBuckets)
		{
			_FRsets = FRsets(__originN);
		}
			else
			{
				FRsets().swap(_FRsets);
			}
			for (auto& refs : _statefulHyperG)
			{
				std::vector<StatefulSampleRef>().swap(refs);
			}
			std::vector<std::vector<StatefulSampleRef>>().swap(_statefulHyperG);
			for (auto& sample : _statefulSamples)
			{
				std::vector<int>().swap(sample.origins);
				for (auto& uses : sample.uses)
				{
					std::vector<StatefulUse>().swap(uses);
				}
				std::vector<std::vector<StatefulUse>>().swap(sample.uses);
			}
			std::vector<StatefulSample>().swap(_statefulSamples);
			std::vector<int>().swap(__statefulOriginLocal);
			__numRRsets = 0;
			trim_allocator_memory();
			LogProgress("release_rr_index_done");
		}

	/// Refresh the hypergraph
		void refresh_hypergraph()
		{
			clear_stateful_metadata();
			release_rr_index(true);
		}

	/// Release memory
	void release_memory()
	{
		LogProgress("release_hypergraph_memory_start");
		release_rr_index(false);
		std::vector<bool>().swap(__vecVisitBool);
		Nodelist().swap(__vecVisitNode);
		trim_allocator_memory();
		LogProgress("release_hypergraph_memory_done");
	}
};

using THyperGraph = HyperGraph;
using PHyperGraph = std::shared_ptr<THyperGraph>;
