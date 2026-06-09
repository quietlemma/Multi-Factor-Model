#include <sys/stat.h>
#include <sys/types.h>

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct Options {
    string input;
    string output;
    string model = "MFIC";
    int factors = 2;
    double olo_p = 0.8;
    uint64_t seed = 0;
    bool undirected = false;
    int compact_n = -1;
    uint64_t progress_edges = 10000000ULL;
    bool binary = false;
    bool both = false;
};

static void usage(const char *prog) {
    cerr
        << "Usage: " << prog << " -input EDGE_LIST -output OUT_DIR "
        << "-model MFIC|MFLT|MFOLO -factors 2|3|5 [-seed 0] [-undirected]\n"
        << "       [-compact-n ORIGINAL_NODE_COUNT] [-progress-edges 10000000]\n"
        << "       [-olo-p 0.8]  # compatibility option; current MFOLO uses factor-aware random probabilities\n"
        << "       [-binary | -both]\n\n"
        << "Output files:\n"
        << "  OUT_DIR/attribute.txt with n, m, origin_n\n"
        << "  OUT_DIR/graph_ic.inf for MFIC, OUT_DIR/graph_lt.inf for MFLT,\n"
        << "  or OUT_DIR/graph_olo.inf for MFOLO\n"
        << "  OUT_DIR/graph_ic.bin / graph_lt.bin / graph_olo.bin when -binary or -both is used\n\n"
        << "Node ids in the output are remapped as:\n"
        << "  0 .. origin_n-1            original origin nodes\n"
        << "  origin_n .. n-1            generated factor nodes\n";
}

static bool parse_args(int argc, char **argv, Options &opt) {
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-input" && i + 1 < argc) {
            opt.input = argv[++i];
        } else if (arg == "-output" && i + 1 < argc) {
            opt.output = argv[++i];
        } else if (arg == "-model" && i + 1 < argc) {
            opt.model = argv[++i];
        } else if (arg == "-factors" && i + 1 < argc) {
            opt.factors = atoi(argv[++i]);
        } else if (arg == "-olo-p" && i + 1 < argc) {
            opt.olo_p = atof(argv[++i]);
        } else if (arg == "-seed" && i + 1 < argc) {
            opt.seed = strtoull(argv[++i], nullptr, 10);
        } else if (arg == "-undirected") {
            opt.undirected = true;
        } else if (arg == "-compact-n" && i + 1 < argc) {
            opt.compact_n = atoi(argv[++i]);
        } else if (arg == "-progress-edges" && i + 1 < argc) {
            opt.progress_edges = strtoull(argv[++i], nullptr, 10);
        } else if (arg == "-binary") {
            opt.binary = true;
        } else if (arg == "-both") {
            opt.both = true;
        } else if (arg == "-help" || arg == "--help" || arg == "-h") {
            return false;
        } else {
            cerr << "Unknown or incomplete argument: " << arg << "\n";
            return false;
        }
    }

    if (opt.input.empty() || opt.output.empty()) return false;
    if (!(opt.model == "MFIC" || opt.model == "MFLT" || opt.model == "MFOLO")) return false;
    if (!(opt.factors == 2 || opt.factors == 3 || opt.factors == 5)) return false;
    if (opt.olo_p < 0.0 || opt.olo_p > 1.0) return false;
    if (opt.compact_n == 0 || opt.compact_n < -1) return false;
    if (opt.binary && opt.both) return false;
    return true;
}

static bool create_directories(const string &path) {
    if (path.empty()) return false;
    for (size_t i = 1; i <= path.size(); ++i) {
        if (i < path.size() && path[i] != '/') continue;
        string current = path.substr(0, i);
        while (current.size() > 1 && current[current.size() - 1] == '/') {
            current.erase(current.size() - 1);
        }
        if (current.empty() || current == "/") continue;
        if (mkdir(current.c_str(), 0755) != 0 && errno != EEXIST) {
            cerr << "Failed to create directory " << current << ": " << strerror(errno) << "\n";
            return false;
        }
    }
    return true;
}

static string join_path(string dir, const string &file) {
    if (!dir.empty() && dir[dir.size() - 1] != '/') dir.push_back('/');
    return dir + file;
}

static string binary_path_from_text_path(const string &text_path) {
    if (text_path.size() >= 4 && text_path.substr(text_path.size() - 4) == ".inf") {
        return text_path.substr(0, text_path.size() - 4) + ".bin";
    }
    return text_path + ".bin";
}

static bool parse_edge_line(const string &line, long long &src, long long &dst) {
    if (line.empty() || line[0] == '#') return false;
    stringstream ss(line);
    return static_cast<bool>(ss >> src >> dst);
}

static int get_or_add_id(unordered_map<long long, int> &id_map, vector<uint64_t> &indeg, long long raw_id) {
    auto it = id_map.find(raw_id);
    if (it != id_map.end()) return it->second;
    int id = static_cast<int>(id_map.size());
    id_map.emplace(raw_id, id);
    indeg.push_back(0);
    return id;
}

static int find_id(const unordered_map<long long, int> &id_map, long long raw_id) {
    auto it = id_map.find(raw_id);
    if (it == id_map.end()) {
        cerr << "Internal error: raw id not found in second pass: " << raw_id << "\n";
        exit(2);
    }
    return it->second;
}

static int compact_id(long long raw_id, int compact_n) {
    if (raw_id < 0 || raw_id >= compact_n) {
        cerr << "Compact input id out of range: " << raw_id
             << " not in [0, " << compact_n << ")\n";
        exit(2);
    }
    return static_cast<int>(raw_id);
}

static vector<double> factor_to_origin_weights(const string &model, int factors, mt19937_64 &rng) {
    vector<double> weights(factors, 0.0);
    if (model == "MFOLO") {
        for (int i = 0; i < factors; ++i) {
            int factor_index = i + 1;
            if (factors == 2) {
                if (factor_index == 1) {
                    uniform_real_distribution<double> d(0.0, 0.3);
                    weights[i] = d(rng);
                } else {
                    uniform_real_distribution<double> d(0.2, 0.3);
                    weights[i] = d(rng);
                }
            } else if (factors == 3) {
                if (factor_index == 1 || factor_index == 3) {
                    uniform_real_distribution<double> d(0.0, 0.3);
                    weights[i] = d(rng);
                } else {
                    uniform_real_distribution<double> d(0.2, 0.4);
                    weights[i] = d(rng);
                }
            } else {
                if (factor_index != 2) {
                    uniform_real_distribution<double> d(0.0, 0.2);
                    weights[i] = d(rng);
                } else {
                    uniform_real_distribution<double> d(0.2, 0.4);
                    weights[i] = d(rng);
                }
            }
        }
        return weights;
    }
    if (model == "MFIC") {
        uniform_real_distribution<double> low(0.0, 0.2);
        uniform_real_distribution<double> mid(0.2, 0.3);
        weights[0] = low(rng);
        if (factors >= 2) weights[1] = mid(rng);
        for (int i = 2; i < factors; ++i) weights[i] = low(rng);
        return weights;
    }

    if (factors == 2) {
        uniform_real_distribution<double> d(0.0, 0.3);
        weights[0] = d(rng);
        weights[1] = 1.0 - weights[0];
    } else if (factors == 3) {
        uniform_real_distribution<double> d(0.0, 0.3);
        weights[0] = d(rng);
        weights[2] = d(rng);
        weights[1] = 1.0 - weights[0] - weights[2];
    } else {
        uniform_real_distribution<double> d(0.0, 0.2);
        weights[0] = d(rng);
        weights[2] = d(rng);
        weights[3] = d(rng);
        weights[4] = d(rng);
        weights[1] = 1.0 - weights[0] - weights[2] - weights[3] - weights[4];
    }
    return weights;
}

static uint64_t splitmix64_stateless(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static double uniform01_from_keys(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
    uint64_t x = splitmix64_stateless(a ^ 0x243f6a8885a308d3ULL);
    x = splitmix64_stateless(x ^ b);
    x = splitmix64_stateless(x ^ c);
    x = splitmix64_stateless(x ^ d);
    const double scale = 1.0 / 9007199254740992.0;  // 2^53
    return static_cast<double>(x >> 11) * scale;
}

static pair<double, double> mflt_parent_gate_range(int factors, int factor_index) {
    if (factors == 2) {
        if (factor_index == 1) return make_pair(0.10, 0.45);
        return make_pair(0.55, 0.95);
    }
    if (factors == 3) {
        if (factor_index == 1) return make_pair(0.05, 0.30);
        if (factor_index == 2) return make_pair(0.35, 0.90);
        return make_pair(0.15, 0.55);
    }
    if (factor_index == 1) return make_pair(0.02, 0.18);
    if (factor_index == 2) return make_pair(0.45, 0.95);
    if (factor_index == 3) return make_pair(0.12, 0.45);
    if (factor_index == 4) return make_pair(0.05, 0.22);
    return make_pair(0.02, 0.18);
}

static double mflt_parent_gate(int factors, int factor_index, long long src_raw, long long dst_raw, uint64_t seed) {
    pair<double, double> range = mflt_parent_gate_range(factors, factor_index);
    double u = uniform01_from_keys(
        seed,
        static_cast<uint64_t>(src_raw),
        static_cast<uint64_t>(dst_raw),
        static_cast<uint64_t>(factor_index));
    return range.first + (range.second - range.first) * u;
}

static double aggregate_origin_listen_probability(const vector<double> &weights) {
    // Strict OLO needs one node-level listen probability for each destination node.
    double blocked_prob = 1.0;
    for (double w : weights) {
        blocked_prob *= (1.0 - w);
    }
    double listen_prob = 1.0 - blocked_prob;
    if (listen_prob < 0.0) listen_prob = 0.0;
    if (listen_prob > 1.0) listen_prob = 1.0;
    return listen_prob;
}

struct BinaryEdgeRecord {
    int src;
    int dst;
    float prob;
};

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);

    Options opt;
    if (!parse_args(argc, argv, opt)) {
        usage(argv[0]);
        return 1;
    }
    if (!create_directories(opt.output)) return 1;

    const bool compact_input = opt.compact_n > 0;
    unordered_map<long long, int> id_map;
    if (!compact_input) id_map.reserve(1 << 20);
    vector<uint64_t> indeg;
    if (compact_input) indeg.assign(opt.compact_n, 0);
    uint64_t original_edges = 0;

    cerr << "[mf_aux_graph] pass1: scan input and count indegrees\n";
    if (compact_input) {
        cerr << "[mf_aux_graph] compact input mode: origin_n=" << opt.compact_n
             << ", skip raw-id remapping\n";
    }
    {
        ifstream in(opt.input.c_str());
        if (!in) {
            cerr << "Cannot open input file: " << opt.input << "\n";
            return 1;
        }
        string line;
        long long src_raw, dst_raw;
        auto add_directed = [&](long long s_raw, long long d_raw) {
            int s = compact_input ? compact_id(s_raw, opt.compact_n)
                                  : get_or_add_id(id_map, indeg, s_raw);
            int d = compact_input ? compact_id(d_raw, opt.compact_n)
                                  : get_or_add_id(id_map, indeg, d_raw);
            (void)s;
            indeg[d]++;
            original_edges++;
        };
        while (getline(in, line)) {
            if (!parse_edge_line(line, src_raw, dst_raw)) continue;
            add_directed(src_raw, dst_raw);
            if (opt.undirected) add_directed(dst_raw, src_raw);
            if (opt.progress_edges > 0 && original_edges % opt.progress_edges == 0) {
                cerr << "[mf_aux_graph] pass1 edges=" << original_edges
                     << " nodes=" << (compact_input ? static_cast<size_t>(opt.compact_n) : id_map.size())
                     << "\n";
            }
        }
    }

    const int origin_n = compact_input ? opt.compact_n : static_cast<int>(id_map.size());
    vector<int> factor_base(origin_n, -1);
    int64_t aux_nodes = origin_n;
    uint64_t factor_owner_count = 0;
    for (int node = 0; node < origin_n; ++node) {
        if (indeg[node] == 0) continue;
        factor_base[node] = static_cast<int>(aux_nodes);
        aux_nodes += opt.factors;
        factor_owner_count++;
    }

    uint64_t aux_edges = original_edges * static_cast<uint64_t>(opt.factors)
                       + factor_owner_count * static_cast<uint64_t>(opt.factors);

    string graph_name;
    if (opt.model == "MFIC")
        graph_name = "graph_ic.inf";
    else if (opt.model == "MFLT")
        graph_name = "graph_lt.inf";
    else
        graph_name = "graph_olo.inf";
    const string graph_path = join_path(opt.output, graph_name);
    const string graph_bin_path = binary_path_from_text_path(graph_path);
    const string attr_path = join_path(opt.output, "attribute.txt");

    {
        ofstream attr(attr_path.c_str());
        if (!attr) {
            cerr << "Cannot write attribute file: " << attr_path << "\n";
            return 1;
        }
        attr << "n=" << aux_nodes << "\n";
        attr << "m=" << aux_edges << "\n";
        attr << "origin_n=" << origin_n << "\n";
    }

        cerr << "[mf_aux_graph] original_nodes=" << origin_n
             << " original_edges=" << original_edges
             << " aux_nodes=" << aux_nodes
             << " aux_edges=" << aux_edges
             << " model=" << opt.model
             << " factors=" << opt.factors
             << "\n";

    const bool write_text = !opt.binary;
    const bool write_bin = opt.binary || opt.both;

    ofstream out_text;
    if (write_text) {
        out_text.open(graph_path.c_str());
        if (!out_text) {
            cerr << "Cannot write graph file: " << graph_path << "\n";
            return 1;
        }
        out_text << setprecision(12);
    }

    ofstream out_bin;
    if (write_bin) {
        out_bin.open(graph_bin_path.c_str(), ios::binary);
        if (!out_bin) {
            cerr << "Cannot write binary graph file: " << graph_bin_path << "\n";
            return 1;
        }
    }

    mt19937_64 rng(opt.seed);
    uint64_t written = 0;

    auto write_edge = [&](int src, int dst, double prob) {
        if (write_text) {
            out_text << src << ' ' << dst << ' ' << prob << '\n';
        }
        if (write_bin) {
            BinaryEdgeRecord rec{src, dst, static_cast<float>(prob)};
            out_bin.write(reinterpret_cast<const char *>(&rec), sizeof(rec));
        }
        written++;
    };

    cerr << "[mf_aux_graph] write factor -> origin edges\n";
    for (int node = 0; node < origin_n; ++node) {
        if (factor_base[node] < 0) continue;
        vector<double> weights = factor_to_origin_weights(opt.model, opt.factors, rng);
        double origin_listen_prob = (opt.model == "MFOLO")
            ? aggregate_origin_listen_probability(weights)
            : 0.0;
        for (int factor = 0; factor < opt.factors; ++factor) {
            double prob = (opt.model == "MFOLO") ? origin_listen_prob : weights[factor];
            write_edge(factor_base[node] + factor, node, prob);
        }
    }

    cerr << "[mf_aux_graph] pass2: write origin -> factor edges\n";
    {
        ifstream in(opt.input.c_str());
        if (!in) {
            cerr << "Cannot reopen input file: " << opt.input << "\n";
            return 1;
        }
        string line;
        long long src_raw, dst_raw;
        uint64_t pass2_edges = 0;
        auto write_directed = [&](long long s_raw, long long d_raw) {
            int src = compact_input ? compact_id(s_raw, origin_n) : find_id(id_map, s_raw);
            int dst = compact_input ? compact_id(d_raw, origin_n) : find_id(id_map, d_raw);
            pass2_edges++;
            if (factor_base[dst] < 0 || indeg[dst] == 0) return;
            for (int factor = 0; factor < opt.factors; ++factor) {
                double prob = 1.0 / static_cast<double>(indeg[dst]);
                if (opt.model == "MFLT") {
                    // Goal-B variant: each factor keeps its own parent gate profile.
                    // This makes the effective origin-level LT distribution truly factor-aware.
                    const double gate = mflt_parent_gate(opt.factors, factor + 1, s_raw, d_raw, opt.seed);
                    prob = gate / static_cast<double>(indeg[dst]);
                }
                write_edge(src, factor_base[dst] + factor, prob);
            }
            if (opt.progress_edges > 0 && pass2_edges % opt.progress_edges == 0) {
                cerr << "[mf_aux_graph] pass2 input_edges=" << pass2_edges
                     << "/" << original_edges
                     << " written_aux_edges=" << written
                     << "/" << aux_edges << "\n";
            }
        };
        while (getline(in, line)) {
            if (!parse_edge_line(line, src_raw, dst_raw)) continue;
            write_directed(src_raw, dst_raw);
            if (opt.undirected) write_directed(dst_raw, src_raw);
        }
    }

    if (write_text) out_text.close();
    if (write_bin) out_bin.close();
    if (written != aux_edges) {
        cerr << "Written edge count mismatch: expected " << aux_edges << ", got " << written << "\n";
        return 2;
    }

    if (write_text) cerr << "[mf_aux_graph] text_done: " << graph_path << "\n";
    if (write_bin) cerr << "[mf_aux_graph] binary_done: " << graph_bin_path << "\n";
    cerr << "[mf_aux_graph] run IMM with:\n";
    cerr << "  ./imm_discrete -dataset " << opt.output << "/ -k 50 -model "
         << ((opt.model == "MFIC") ? "IC" : (opt.model == "MFLT" ? "LT" : "OLO")) << " -epsilon 0.5\n";
    return 0;
}
