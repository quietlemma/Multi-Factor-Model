#define HEAD_INFO
//#define HEAD_TRACE
#include "sfmt/SFMT.h"
#include "head.h"
using namespace std;
typedef double (*pf)(int, int);
class Graph
{
public:
    int n, k, origin_n;
    int64 m;
    vector<int> inDeg;
    vector<vector<int>> gT;

#ifdef CONTINUOUS
    vector<vector<dpair>> probT;
#endif
#ifdef DISCRETE
    vector<vector<double>> probT;
#endif


    enum InfluModel {IC, LT, OLO, CONT};
    InfluModel influModel;
    void setInfuModel(InfluModel p)
    {
        influModel = p;
        TRACE(influModel == IC);
        TRACE(influModel == LT);
        TRACE(influModel == OLO);
        TRACE(influModel == CONT);
    }

    string folder;
    string graph_file;
    string binaryGraphPath() const
    {
        if (graph_file.size() >= 4 && graph_file.substr(graph_file.size() - 4) == ".inf")
            return graph_file.substr(0, graph_file.size() - 4) + ".bin";
        return graph_file + ".bin";
    }
    bool hasBinaryGraph() const
    {
        return access(binaryGraphPath().c_str(), F_OK) == 0;
    }
    void readNM()
    {
        ifstream cin((folder + "attribute.txt").c_str());
        ASSERT(!cin == false);
        origin_n = -1;
        string s;
        while (cin >> s)
        {
            if (s.substr(0, 2) == "n=")
            {
                n = atoi(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 2) == "m=")
            {
                m = atoll(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 9) == "origin_n=")
            {
                origin_n = atoi(s.substr(9).c_str());
                continue;
            }
            ASSERT(false);
        }
        if (origin_n < 0)
            origin_n = n;
        ASSERT(origin_n > 0);
        ASSERT(origin_n <= n);
        TRACE(n, m );
        cin.close();
    }
#ifdef DISCRETE
    void add_edge(int a, int b, double p)
    {
        probT[b].push_back(p);
        gT[b].push_back(a);
        inDeg[b]++;
    }
#endif
#ifdef CONTINUOUS
    void add_edge(int a, int b, double p1, double p2)
    {
        probT[b].push_back(MP(p1, p2));
        gT[b].push_back(a);
        inDeg[b]++;
    }
#endif
    vector<bool> hasnode;
    void readGraph()
    {
        FILE *fin = fopen((graph_file).c_str(), "r");
        ASSERT(fin != nullptr);
        int64 readCnt = 0;
        for (int64 i = 0; i < m; i++)
        {
            readCnt ++;
            int a, b;
#ifdef DISCRETE
            double p;
            int c = fscanf(fin, "%d%d%lf", &a, &b, &p);
            ASSERTT(c == 3, a, b, p, c);
#endif
#ifdef CONTINUOUS
            double p1, p2;
            int c = fscanf(fin, "%d%d%lf%lf", &a, &b, &p1, &p2);
            ASSERT(c == 4);
#endif

            //TRACE_LINE(a, b);
            ASSERT( a < n );
            ASSERT( b < n );
            hasnode[a] = true;
            hasnode[b] = true;
#ifdef DISCRETE
            add_edge(a, b, p);
#endif
#ifdef CONTINUOUS
            add_edge(a, b, p1, p2);
#endif
        }
        TRACE_LINE_END();
        int s = 0;
        for (int i = 0; i < n; i++)
            if (hasnode[i])
                s++;
        INFO(s);
        ASSERT(readCnt == m);
        fclose(fin);
    }
#ifdef DISCRETE
    void readGraphBin()
    {
        string graph_file_bin = binaryGraphPath();
        INFO(graph_file_bin);
        FILE *fin = fopen(graph_file_bin.c_str(), "rb");
        ASSERT(fin != nullptr);
        struct EdgeRecord
        {
            int a;
            int b;
            float p;
        };
        const size_t chunkRecords = 1 << 20;
        vector<EdgeRecord> buf(chunkRecords);
        int64 readCnt = 0;
        while (readCnt < m)
        {
            size_t need = static_cast<size_t>(min<int64>(chunkRecords, m - readCnt));
            size_t got = fread(buf.data(), sizeof(EdgeRecord), need, fin);
            ASSERTT(got == need, readCnt, need, got);
            for (size_t i = 0; i < got; i++)
            {
                int a = buf[i].a;
                int b = buf[i].b;
                float p = buf[i].p;
                ASSERT(a < n);
                ASSERT(b < n);
                hasnode[a] = true;
                hasnode[b] = true;
                add_edge(a, b, p);
            }
            readCnt += got;
        }
        int s = 0;
        for (int i = 0; i < n; i++)
            if (hasnode[i])
                s++;
        INFO(s);
        ASSERT(readCnt == m);
        fclose(fin);
    }
#endif
    Graph(string folder, string graph_file): folder(folder), graph_file(graph_file)
    {
        readNM();

        //init vector
        FOR(i, n)
        {
            gT.push_back(vector<int>());
            hasnode.push_back(false);
#ifdef DISCRETE
            probT.push_back(vector<double>());
#endif
#ifdef CONTINUOUS
            probT.push_back(vector<dpair>());
#endif
            //hyperGT.push_back(vector<int>());
            inDeg.push_back(0);
        }
        if (hasBinaryGraph())
        {
            INFO("read_binary_graph", binaryGraphPath());
            readGraphBin();
        }
        else
        {
            INFO("read_text_graph", graph_file);
            readGraph();
        }
        //system("sleep 10000");
    }

};
double sqr(double t)
{
    return t * t;
}
