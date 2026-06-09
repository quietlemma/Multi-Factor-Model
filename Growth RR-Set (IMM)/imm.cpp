#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"

class Argument{
public:
    int k = 0;
    vector<int> k_list;
    string dataset;
    double epsilon = 0.0;
    string model;
    double T = 0.0;
};

#include "graph.h"
#include "infgraph.h"
#include "imm.h"



void run_with_parameter(InfGraph &g, const Argument & arg)
{
        Timer::clearAll();
        g.init_hyper_graph();
        sfmt_init_gen_rand(&g.sfmtSeed, 95082);

        cout << "--------------------------------------------------------------------------------" << endl;
        cout << arg.dataset << " k=" << arg.k << " epsilon=" << arg.epsilon <<   " " << arg.model << endl;

        if (arg.k == 0)
        {
            g.seedSet.clear();
            cout << "estimated_influence=0" << endl;
            INFO(g.seedSet);
            Timer::show();
            return;
        }

        double estimated_influence = Imm::InfluenceMaximize(g, arg);
        cout << "estimated_influence=" << estimated_influence << endl;

        INFO(g.seedSet);
       
    Timer::show();
}

static vector<int> parse_k_list(const string &text)
{
    vector<int> result;
    string token;
    stringstream ss(text);
    while (getline(ss, token, ','))
    {
        if (token.empty()) continue;
        result.push_back(atoi(token.c_str()));
    }
    return result;
}

static vector<int> normalize_k_list(const vector<int> &input)
{
    vector<int> result;
    for (int value : input)
    {
        if (value < 0) continue;
        result.push_back(value);
    }
    sort(result.begin(), result.end());
    result.erase(unique(result.begin(), result.end()), result.end());
    return result;
}

static void print_result_block(InfGraph &g, const Argument &arg, int k, double estimated_influence,
                               const vector<int> &seed_set, bool show_total_time, double total_seconds)
{
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << arg.dataset << " k=" << k << " epsilon=" << arg.epsilon << " " << arg.model << endl;
    cout << "estimated_influence=" << estimated_influence << endl;
    g.seedSet = seed_set;
    INFO(g.seedSet);
    if (!show_total_time)
        return;
    cout << "########## Timer ##########" << endl;
    char str[100];
    sprintf(str, "%.6lf", total_seconds);
    string s = str;
    if ((int)s.size() < 15) s = " " + s;
    char t[200];
    memset(t, 0, sizeof t);
    sprintf(t, "Spend %s seconds on %s", s.c_str(), "InfluenceMaximize(Total Time)");
    cout << t << endl;
}

static void run_with_reused_k_list(InfGraph &g, const Argument &arg)
{
    vector<int> k_list = normalize_k_list(arg.k_list);
    ASSERT(!k_list.empty());

    Timer::clearAll();
    g.init_hyper_graph();
    sfmt_init_gen_rand(&g.sfmtSeed, 95082);
    cout << "reuse_k_list_mode=";
    for (size_t i = 0; i < k_list.size(); ++i)
    {
        if (i) cout << ",";
        cout << k_list[i];
    }
    cout << endl;

    int kmax = k_list.back();
    if (kmax == 0)
    {
        print_result_block(g, arg, 0, 0.0, vector<int>(), false, 0.0);
        return;
    }

    vector<int> positive_k;
    for (int value : k_list)
        if (value > 0)
            positive_k.push_back(value);

    Argument run_arg = arg;
    run_arg.k = kmax;
    vector<PrefixInfluenceEntry> entries = Imm::InfluenceMaximizeWithReuse(g, run_arg, positive_k);

    size_t next_entry = 0;
    for (int current_k : k_list)
    {
        if (current_k == 0)
        {
            print_result_block(g, arg, 0, 0.0, vector<int>(), false, 0.0);
            continue;
        }
        ASSERT(next_entry < entries.size());
        ASSERT(entries[next_entry].k == current_k);
        print_result_block(
            g,
            arg,
            current_k,
            entries[next_entry].estimated_influence,
            entries[next_entry].seed_set,
            true,
            entries[next_entry].total_seconds);
        next_entry++;
    }
}

void Run(int argn, char **argv)
{
    Argument arg;


    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-help") || argv[i] == string("--help") || argn == 1)
        {
            cout << "./tim -dataset *** -epsilon *** -k ***  -model IC|LT|OLO|TR|CONT " << endl;
            return ;
        }
        if (argv[i] == string("-dataset")) 
            arg.dataset = argv[i + 1];
        if (argv[i] == string("-epsilon")) 
            arg.epsilon = atof(argv[i + 1]);
        if (argv[i] == string("-T")) 
            arg.T = atof(argv[i + 1]);
        if (argv[i] == string("-k")) 
            arg.k = atoi(argv[i + 1]);
        if (argv[i] == string("-k-list"))
            arg.k_list = parse_k_list(argv[i + 1]);
        if (argv[i] == string("-model"))
            arg.model = argv[i + 1];
    }
    ASSERT(arg.dataset != "");
    ASSERT(arg.model == "IC" || arg.model == "LT" || arg.model == "OLO" || arg.model == "TR" || arg.model=="CONT");

    string graph_file;
    if (arg.model == "IC")
        graph_file = arg.dataset + "graph_ic.inf";
    else if (arg.model == "LT")
        graph_file = arg.dataset + "graph_lt.inf";
    else if (arg.model == "OLO")
        graph_file = arg.dataset + "graph_olo.inf";
    else if (arg.model == "TR")
        graph_file = arg.dataset + "graph_tr.inf";
    else if (arg.model == "CONT")
        graph_file = arg.dataset + "graph_cont.inf";
    else
        ASSERT(false);

    InfGraph g(arg.dataset, graph_file);


    if (arg.model == "IC")
        g.setInfuModel(InfGraph::IC);
    else if (arg.model == "LT")
        g.setInfuModel(InfGraph::LT);
    else if (arg.model == "OLO")
        g.setInfuModel(InfGraph::OLO);
    else if (arg.model == "TR")
        g.setInfuModel(InfGraph::IC);
    else if (arg.model == "CONT")
        g.setInfuModel(InfGraph::CONT);
    else
        ASSERT(false);

    INFO(arg.T);

    arg.k_list = normalize_k_list(arg.k_list);
    if (arg.k_list.empty())
        arg.k_list.push_back(max(arg.k, 0));

    if (arg.k_list.size() > 1)
    {
        run_with_reused_k_list(g, arg);
        return;
    }

    Argument current_arg = arg;
    current_arg.k = arg.k_list[0];
    run_with_parameter(g, current_arg);
}


int main(int argn, char **argv)
{
    __head_version = "v1";
    OutputInfo info(argn, argv);


    Run( argn, argv );
}
