#include <sys/stat.h>
#include <sys/types.h>

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

using BinaryEdge = pair<uint32_t, float>;
using BinaryGraph = vector<vector<BinaryEdge>>;

struct Options
{
	string input;
	string output;
	string graphname;
	string model = "MFIC";
	string mfltParentGate = "original";
	string outputFormat = "text";
	double icvariantsP = 0.778;
	double ic2MixRatio = 0.10;
	int factors = 2;
	uint64_t seed = 0;
	uint64_t progressEdges = 10000000ULL;
	uint64_t compactN = 0;
	bool undirected = false;
	bool skipHeader = false;
};

static void usage(const char* prog)
{
	cerr
		<< "Usage: " << prog << " -input EDGE_LIST -output OUT_DIR "
		<< "-model MFIC|MFIC2|MFIC2MIX|MFLT|MFOLO -factors 2|3|5 [-seed 0] [-undirected]\n"
		<< "       [-compact-n ORIGINAL_NODE_COUNT] [-progress-edges 10000000]\n"
		<< "       [-skip-header]\n"
		<< "       [-gname graph_lt.inf]\n"
		<< "       [-format text|bin|both] [-bin] [-both]\n"
		<< "       [-icvariants-p 0.778] [-ic2-mix-ratio 0.10]\n"
		<< "       [-mflt-parent-gate original|factor-aware] [-mflt-gate] [-mflt-original]\n\n"
		<< "Text output is OPIM weighted format: first line is \"n m\", followed by src dst prob.\n"
		<< "Bin output writes OPIM's serialized reverse graph directly; run OPIM with -func=1, no formatting step.\n"
		<< "Default MFLT parent gate is original, i.e. origin->factor prob = 1 / indeg(dst).\n"
		<< "Use -mflt-parent-gate factor-aware or -mflt-gate to enable the later gate / indeg(dst) variant.\n";
}

static bool parse_args(int argc, char** argv, Options& opt)
{
	for (int i = 1; i < argc; ++i)
	{
		string arg = argv[i];
		if (arg == "-input" && i + 1 < argc) opt.input = argv[++i];
		else if (arg == "-output" && i + 1 < argc) opt.output = argv[++i];
		else if (arg == "-gname" && i + 1 < argc) opt.graphname = argv[++i];
		else if (arg == "-model" && i + 1 < argc) opt.model = argv[++i];
		else if (arg == "-factors" && i + 1 < argc) opt.factors = static_cast<int>(strtol(argv[++i], nullptr, 10));
		else if (arg == "-seed" && i + 1 < argc) opt.seed = strtoull(argv[++i], nullptr, 10);
		else if (arg == "-progress-edges" && i + 1 < argc) opt.progressEdges = strtoull(argv[++i], nullptr, 10);
		else if (arg == "-compact-n" && i + 1 < argc) opt.compactN = strtoull(argv[++i], nullptr, 10);
		else if (arg == "-icvariants-p" && i + 1 < argc) opt.icvariantsP = atof(argv[++i]);
		else if (arg == "-ic2-mix-ratio" && i + 1 < argc) opt.ic2MixRatio = atof(argv[++i]);
		else if ((arg == "-format" || arg == "-output-format") && i + 1 < argc) opt.outputFormat = argv[++i];
		else if (arg == "-mflt-parent-gate" && i + 1 < argc) opt.mfltParentGate = argv[++i];
		else if (arg == "-mflt-gate") opt.mfltParentGate = "factor-aware";
		else if (arg == "-mflt-original") opt.mfltParentGate = "original";
		else if (arg == "-bin" || arg == "-binary") opt.outputFormat = "bin";
		else if (arg == "-both") opt.outputFormat = "both";
		else if (arg == "-skip-header" || arg == "-input-has-header") opt.skipHeader = true;
		else if (arg == "-undirected") opt.undirected = true;
		else if (arg == "-help" || arg == "--help" || arg == "-h") return false;
		else
		{
			cerr << "Unknown or incomplete argument: " << arg << "\n";
			return false;
		}
	}

	if (opt.input.empty() || opt.output.empty()) return false;
	if (!(opt.model == "MFIC" || opt.model == "MFIC2" || opt.model == "MFIC2MIX" || opt.model == "MFLT" || opt.model == "MFOLO")) return false;
	if (!(opt.factors == 2 || opt.factors == 3 || opt.factors == 5)) return false;
	if (!(opt.mfltParentGate == "original" || opt.mfltParentGate == "factor-aware")) return false;
	if (!(opt.outputFormat == "text" || opt.outputFormat == "bin" || opt.outputFormat == "both")) return false;
	if (opt.compactN > numeric_limits<uint32_t>::max()) return false;
	if (opt.ic2MixRatio < 0.0 || opt.ic2MixRatio > 1.0) return false;
	if (opt.graphname.empty())
	{
		const string suffix = opt.outputFormat == "bin" ? ".bin" : ".inf";
		if (opt.model == "MFIC") opt.graphname = "graph_ic" + suffix;
		else if (opt.model == "MFIC2") opt.graphname = "graph_ic2" + suffix;
		else if (opt.model == "MFIC2MIX") opt.graphname = "graph_ic2mix" + suffix;
		else if (opt.model == "MFLT") opt.graphname = "graph_lt" + suffix;
		else opt.graphname = "graph_olo" + suffix;
	}
	return true;
}

static bool create_directories(const string& path)
{
	if (path.empty()) return false;
	for (size_t i = 1; i <= path.size(); ++i)
	{
		if (i < path.size() && path[i] != '/') continue;
		string current = path.substr(0, i);
		while (current.size() > 1 && current[current.size() - 1] == '/') current.erase(current.size() - 1);
		if (current.empty() || current == "/") continue;
		if (mkdir(current.c_str(), 0755) != 0 && errno != EEXIST)
		{
			cerr << "Failed to create directory " << current << ": " << strerror(errno) << "\n";
			return false;
		}
	}
	return true;
}

static string join_path(string dir, const string& file)
{
	if (!dir.empty() && dir[dir.size() - 1] != '/') dir.push_back('/');
	return dir + file;
}

static bool ends_with(const string& value, const string& suffix)
{
	return value.size() >= suffix.size()
		&& value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static string binary_graph_path(const string& outputDir, const string& graphname)
{
	const string graphPath = join_path(outputDir, graphname);
	if (ends_with(graphname, ".bin") || ends_with(graphname, ".vec.rvs.graph")) return graphPath;
	return graphPath + ".vec.rvs.graph";
}

template <typename T>
static void write_binary_value(ofstream& out, const T& value)
{
	out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

static bool save_opim_reverse_graph_binary(const string& path, const BinaryGraph& graph)
{
	static_assert(sizeof(BinaryEdge) == sizeof(uint32_t) + sizeof(float),
	              "OPIM binary edge layout must match pair<uint32_t,float> serialization.");
	ofstream out(path.c_str(), ios::binary);
	if (!out)
	{
		cerr << "Cannot write binary graph file: " << path << "\n";
		return false;
	}
	const size_t outerSize = graph.size();
	write_binary_value(out, outerSize);
	for (const auto& nbrs : graph)
	{
		const size_t innerSize = nbrs.size();
		write_binary_value(out, innerSize);
		if (!nbrs.empty())
		{
			out.write(reinterpret_cast<const char*>(nbrs.data()),
			          static_cast<streamsize>(sizeof(BinaryEdge) * nbrs.size()));
		}
	}
	return static_cast<bool>(out);
}

static bool parse_edge_line(const string& line, long long& src, long long& dst)
{
	if (line.empty() || line[0] == '#') return false;
	stringstream ss(line);
	return static_cast<bool>(ss >> src >> dst);
}

static uint32_t get_or_add_id(unordered_map<long long, uint32_t>& idMap, vector<uint64_t>& indeg, long long rawId)
{
	auto it = idMap.find(rawId);
	if (it != idMap.end()) return it->second;
	if (idMap.size() >= numeric_limits<uint32_t>::max())
	{
		cerr << "Too many origin nodes for OPIM uint32_t node ids: " << idMap.size() << "\n";
		exit(2);
	}
	uint32_t id = static_cast<uint32_t>(idMap.size());
	idMap.emplace(rawId, id);
	indeg.push_back(0);
	return id;
}

static uint32_t find_id(const unordered_map<long long, uint32_t>& idMap, long long rawId)
{
	auto it = idMap.find(rawId);
	if (it == idMap.end())
	{
		cerr << "Internal error: raw id not found in second pass: " << rawId << "\n";
		exit(2);
	}
	return it->second;
}

static uint32_t compact_id(long long rawId, uint64_t compactN)
{
	if (rawId < 0 || static_cast<uint64_t>(rawId) >= compactN)
	{
		cerr << "Compact input id out of range: " << rawId << " not in [0, " << compactN << ")\n";
		exit(2);
	}
	return static_cast<uint32_t>(rawId);
}

static vector<double> factor_to_origin_weights(const string& model, int factors, mt19937_64& rng)
{
	vector<double> weights(factors, 0.0);
	if (model == "MFOLO")
	{
		for (int i = 0; i < factors; ++i)
		{
			const int factorIndex = i + 1;
			if (factors == 2)
			{
				if (factorIndex == 1) weights[i] = uniform_real_distribution<double>(0.0, 0.3)(rng);
				else weights[i] = uniform_real_distribution<double>(0.2, 0.3)(rng);
			}
			else if (factors == 3)
			{
				if (factorIndex == 2) weights[i] = uniform_real_distribution<double>(0.2, 0.4)(rng);
				else weights[i] = uniform_real_distribution<double>(0.0, 0.3)(rng);
			}
			else
			{
				if (factorIndex == 2) weights[i] = uniform_real_distribution<double>(0.2, 0.4)(rng);
				else weights[i] = uniform_real_distribution<double>(0.0, 0.2)(rng);
			}
		}
		return weights;
	}
	if (model == "MFIC")
	{
		weights[0] = uniform_real_distribution<double>(0.0, 0.2)(rng);
		if (factors >= 2) weights[1] = uniform_real_distribution<double>(0.2, 0.3)(rng);
		for (int i = 2; i < factors; ++i) weights[i] = uniform_real_distribution<double>(0.0, 0.2)(rng);
		return weights;
	}
	if (factors == 2)
	{
		weights[0] = uniform_real_distribution<double>(0.0, 0.3)(rng);
		weights[1] = 1.0 - weights[0];
	}
	else if (factors == 3)
	{
		weights[0] = uniform_real_distribution<double>(0.0, 0.3)(rng);
		weights[2] = uniform_real_distribution<double>(0.0, 0.3)(rng);
		weights[1] = 1.0 - weights[0] - weights[2];
	}
	else
	{
		weights[0] = uniform_real_distribution<double>(0.0, 0.2)(rng);
		weights[2] = uniform_real_distribution<double>(0.0, 0.2)(rng);
		weights[3] = uniform_real_distribution<double>(0.0, 0.2)(rng);
		weights[4] = uniform_real_distribution<double>(0.0, 0.2)(rng);
		weights[1] = 1.0 - weights[0] - weights[2] - weights[3] - weights[4];
	}
	return weights;
}

static uint64_t splitmix64_stateless(uint64_t x)
{
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return x ^ (x >> 31);
}

static bool choose_ic2_mix_target(const int node, const uint64_t seed, const double ratio)
{
	if (ratio <= 0.0) return false;
	if (ratio >= 1.0) return true;
	uint64_t x = splitmix64_stateless(seed ^ 0x9e3779b97f4a7c15ULL);
	x = splitmix64_stateless(x ^ static_cast<uint64_t>(node));
	const double scale = 1.0 / 9007199254740992.0;
	const double u = static_cast<double>(x >> 11) * scale;
	return u < ratio;
}

static double uniform01_from_keys(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
{
	uint64_t x = splitmix64_stateless(a ^ 0x243f6a8885a308d3ULL);
	x = splitmix64_stateless(x ^ b);
	x = splitmix64_stateless(x ^ c);
	x = splitmix64_stateless(x ^ d);
	const double scale = 1.0 / 9007199254740992.0;
	return static_cast<double>(x >> 11) * scale;
}

static pair<double, double> mflt_parent_gate_range(int factors, int factorIndex)
{
	if (factors == 2)
	{
		if (factorIndex == 1) return make_pair(0.10, 0.45);
		return make_pair(0.55, 0.95);
	}
	if (factors == 3)
	{
		if (factorIndex == 1) return make_pair(0.05, 0.30);
		if (factorIndex == 2) return make_pair(0.35, 0.90);
		return make_pair(0.15, 0.55);
	}
	if (factorIndex == 1) return make_pair(0.02, 0.18);
	if (factorIndex == 2) return make_pair(0.45, 0.95);
	if (factorIndex == 3) return make_pair(0.12, 0.45);
	if (factorIndex == 4) return make_pair(0.05, 0.22);
	return make_pair(0.02, 0.18);
}

static double mflt_parent_gate(int factors, int factorIndex, long long srcRaw, long long dstRaw, uint64_t seed)
{
	const auto range = mflt_parent_gate_range(factors, factorIndex);
	const double u = uniform01_from_keys(
		seed,
		static_cast<uint64_t>(srcRaw),
		static_cast<uint64_t>(dstRaw),
		static_cast<uint64_t>(factorIndex));
	return range.first + (range.second - range.first) * u;
}

static double aggregate_origin_listen_probability(const vector<double>& weights)
{
	double blockedProb = 1.0;
	for (double w : weights) blockedProb *= 1.0 - w;
	double listenProb = 1.0 - blockedProb;
	if (listenProb < 0.0) listenProb = 0.0;
	if (listenProb > 1.0) listenProb = 1.0;
	return listenProb;
}

int main(int argc, char** argv)
{
	ios::sync_with_stdio(false);

	Options opt;
	if (!parse_args(argc, argv, opt))
	{
		usage(argv[0]);
		return 1;
	}
	if (!create_directories(opt.output)) return 1;

	const bool compactInput = opt.compactN > 0;
	unordered_map<long long, uint32_t> idMap;
	if (!compactInput) idMap.reserve(1 << 20);
	vector<uint64_t> indeg;
	if (compactInput) indeg.assign(static_cast<size_t>(opt.compactN), 0);
	uint64_t originalEdges = 0;

	cerr << "[opim_mf_aux_graph] pass1: scan input and count indegrees\n";
	{
		ifstream in(opt.input.c_str());
		if (!in)
		{
			cerr << "Cannot open input file: " << opt.input << "\n";
			return 1;
		}
		string line;
		long long srcRaw, dstRaw;
		uint64_t lineNo = 0;
		auto add_directed = [&](long long sRaw, long long dRaw)
		{
			const uint32_t src = compactInput ? compact_id(sRaw, opt.compactN) : get_or_add_id(idMap, indeg, sRaw);
			const uint32_t dst = compactInput ? compact_id(dRaw, opt.compactN) : get_or_add_id(idMap, indeg, dRaw);
			(void)src;
			indeg[dst]++;
			originalEdges++;
		};
		while (getline(in, line))
		{
			lineNo++;
			if (opt.skipHeader && lineNo == 1) continue;
			if (!parse_edge_line(line, srcRaw, dstRaw)) continue;
			add_directed(srcRaw, dstRaw);
			if (opt.undirected) add_directed(dstRaw, srcRaw);
			if (opt.progressEdges > 0 && originalEdges % opt.progressEdges == 0)
			{
				cerr << "[opim_mf_aux_graph] pass1 edges=" << originalEdges
				     << " nodes=" << (compactInput ? static_cast<size_t>(opt.compactN) : idMap.size()) << "\n";
			}
		}
	}

	const uint64_t originN64 = compactInput ? opt.compactN : static_cast<uint64_t>(idMap.size());
	if (originN64 >= numeric_limits<uint32_t>::max())
	{
		cerr << "originN exceeds OPIM uint32_t node id capacity: " << originN64 << "\n";
		return 2;
	}
	const uint32_t originN = static_cast<uint32_t>(originN64);
	const uint32_t noFactor = numeric_limits<uint32_t>::max();
	vector<uint32_t> factorBase(static_cast<size_t>(originN), noFactor);
	vector<char> usesFactorTarget(static_cast<size_t>(originN), 0);
	uint64_t factorOwnerCount = 0;
	uint64_t mixedDirectTargetCount = 0;
	uint64_t mixedStatefulTargetCount = 0;
	uint64_t auxNodes = originN;
	for (uint32_t node = 0; node < originN; ++node)
	{
		if (indeg[node] == 0) continue;
		bool useFactor = true;
		if (opt.model == "MFIC2MIX")
			useFactor = choose_ic2_mix_target(static_cast<int>(node), opt.seed, opt.ic2MixRatio);
		usesFactorTarget[node] = useFactor ? 1 : 0;
		if (!useFactor)
		{
			mixedDirectTargetCount++;
			continue;
		}
		if (auxNodes + static_cast<uint64_t>(opt.factors) >= numeric_limits<uint32_t>::max())
		{
			cerr << "Auxiliary graph node ids exceed OPIM uint32_t capacity. "
			     << "Current auxNodes=" << auxNodes << ", factors=" << opt.factors << "\n";
			return 2;
		}
		factorBase[node] = static_cast<uint32_t>(auxNodes);
		auxNodes += opt.factors;
		factorOwnerCount++;
		if (opt.model == "MFIC2MIX")
			mixedStatefulTargetCount++;
	}
	uint64_t auxEdges = 0;
	if (opt.model == "MFIC2MIX")
	{
		for (uint32_t node = 0; node < originN; ++node)
		{
			if (indeg[node] == 0) continue;
			if (usesFactorTarget[node])
				auxEdges += indeg[node] * static_cast<uint64_t>(opt.factors)
				          + static_cast<uint64_t>(opt.factors);
			else
				auxEdges += indeg[node];
		}
	}
	else
	{
		auxEdges = originalEdges * static_cast<uint64_t>(opt.factors)
		         + factorOwnerCount * static_cast<uint64_t>(opt.factors);
	}

	const string graphPath = join_path(opt.output, opt.graphname);
	const string binGraphPath = binary_graph_path(opt.output, opt.graphname);
	const string attrPath = join_path(opt.output, "attribute.txt");
	const bool writeText = opt.outputFormat == "text" || opt.outputFormat == "both";
	const bool writeBin = opt.outputFormat == "bin" || opt.outputFormat == "both";

	ofstream out;
	if (writeText)
	{
		out.open(graphPath.c_str());
		if (!out)
		{
			cerr << "Cannot write graph file: " << graphPath << "\n";
			return 1;
		}
		out << setprecision(12);
		out << auxNodes << ' ' << auxEdges << '\n';
	}

	BinaryGraph reverseGraph;
	if (writeBin) reverseGraph.resize(static_cast<size_t>(auxNodes));

	ofstream attr(attrPath.c_str());
	if (attr)
	{
		attr << "n=" << auxNodes << "\n";
		attr << "m=" << auxEdges << "\n";
		attr << "origin_n=" << originN << "\n";
		attr << "format=" << opt.outputFormat << "\n";
		attr << "icvariants_p=" << opt.icvariantsP << "\n";
		attr << "ic2_mix_ratio=" << opt.ic2MixRatio << "\n";
		if (writeText) attr << "text_graph=" << opt.graphname << "\n";
		if (writeBin) attr << "bin_graph=" << binGraphPath << "\n";
	}

	cerr << "[opim_mf_aux_graph] original_nodes=" << originN
	     << " original_edges=" << originalEdges
	     << " aux_nodes=" << auxNodes
	     << " aux_edges=" << auxEdges
	     << " model=" << opt.model
	     << " factors=" << opt.factors
	     << " icvariants_p=" << opt.icvariantsP
	     << " ic2_mix_ratio=" << opt.ic2MixRatio
	     << " mixed_stateful_targets=" << mixedStatefulTargetCount
	     << " mixed_direct_targets=" << mixedDirectTargetCount
	     << " mflt_parent_gate=" << opt.mfltParentGate
	     << " output_format=" << opt.outputFormat
	     << "\n";

	mt19937_64 rng(opt.seed);
	uint64_t written = 0;
	auto write_edge = [&](uint32_t src, uint32_t dst, double prob)
	{
		if (writeText) out << src << ' ' << dst << ' ' << prob << '\n';
		if (writeBin)
		{
			reverseGraph[static_cast<size_t>(dst)].push_back(
				BinaryEdge(static_cast<uint32_t>(src), static_cast<float>(prob)));
		}
		written++;
	};

	cerr << "[opim_mf_aux_graph] write factor -> origin edges\n";
	for (uint32_t node = 0; node < originN; ++node)
	{
		if (factorBase[node] == noFactor) continue;
		const vector<double> weights = factor_to_origin_weights(opt.model, opt.factors, rng);
		const double originListenProb = opt.model == "MFOLO" ? aggregate_origin_listen_probability(weights) : 0.0;
		for (int factor = 0; factor < opt.factors; ++factor)
		{
			double prob = opt.model == "MFOLO" ? originListenProb : weights[factor];
			if (opt.model == "MFIC2" || opt.model == "MFIC2MIX")
				prob = opt.icvariantsP;
			write_edge(factorBase[node] + factor, node, prob);
		}
	}

	cerr << "[opim_mf_aux_graph] pass2: write origin -> factor edges\n";
	{
		ifstream in(opt.input.c_str());
		if (!in)
		{
			cerr << "Cannot reopen input file: " << opt.input << "\n";
			return 1;
		}
		string line;
		long long srcRaw, dstRaw;
		uint64_t pass2Edges = 0;
		uint64_t lineNo = 0;
		auto write_directed = [&](long long sRaw, long long dRaw)
		{
			const uint32_t src = compactInput ? compact_id(sRaw, originN) : find_id(idMap, sRaw);
			const uint32_t dst = compactInput ? compact_id(dRaw, originN) : find_id(idMap, dRaw);
			pass2Edges++;
			if (indeg[dst] == 0) return;
			if (factorBase[dst] == noFactor)
			{
				write_edge(src, dst, 1.0 / static_cast<double>(indeg[dst]));
				if (opt.progressEdges > 0 && pass2Edges % opt.progressEdges == 0)
				{
					cerr << "[opim_mf_aux_graph] pass2 input_edges=" << pass2Edges
					     << "/" << originalEdges
					     << " written_aux_edges=" << written
					     << "/" << auxEdges << "\n";
				}
				return;
			}
			for (int factor = 0; factor < opt.factors; ++factor)
			{
				double prob = 1.0 / static_cast<double>(indeg[dst]);
				if (opt.model == "MFLT" && opt.mfltParentGate == "factor-aware")
				{
					const double gate = mflt_parent_gate(opt.factors, factor + 1, sRaw, dRaw, opt.seed);
					prob = gate / static_cast<double>(indeg[dst]);
				}
				write_edge(src, factorBase[dst] + factor, prob);
			}
			if (opt.progressEdges > 0 && pass2Edges % opt.progressEdges == 0)
			{
				cerr << "[opim_mf_aux_graph] pass2 input_edges=" << pass2Edges
				     << "/" << originalEdges
				     << " written_aux_edges=" << written
				     << "/" << auxEdges << "\n";
			}
		};
		while (getline(in, line))
		{
			lineNo++;
			if (opt.skipHeader && lineNo == 1) continue;
			if (!parse_edge_line(line, srcRaw, dstRaw)) continue;
			write_directed(srcRaw, dstRaw);
			if (opt.undirected) write_directed(dstRaw, srcRaw);
		}
	}

	if (writeText) out.close();
	if (written != auxEdges)
	{
		cerr << "Written edge count mismatch: expected " << auxEdges << ", got " << written << "\n";
		return 2;
	}
	if (writeBin)
	{
		cerr << "[opim_mf_aux_graph] save OPIM binary reverse graph: " << binGraphPath << "\n";
		if (!save_opim_reverse_graph_binary(binGraphPath, reverseGraph)) return 3;
	}

	if (writeText) cerr << "[opim_mf_aux_graph] text graph: " << graphPath << "\n";
	if (writeBin) cerr << "[opim_mf_aux_graph] binary graph: " << binGraphPath << "\n";
	cerr << "[opim_mf_aux_graph] done\n";
	if (writeText && !writeBin)
	{
		cerr << "[opim_mf_aux_graph] format with:\n";
		cerr << "  ./OPIM1.1.o -func=0 -dir=" << opt.output << " -gname=" << opt.graphname << " -mode=w\n";
	}
	else
	{
		cerr << "[opim_mf_aux_graph] run OPIM directly with:\n";
		cerr << "  ./OPIM1.1.o -func=1 -dir=" << opt.output << " -gname=" << opt.graphname << " -mode=2 -pdist=load\n";
	}
	return 0;
}
