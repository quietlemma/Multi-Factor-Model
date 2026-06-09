
typedef pair<double,int> dipair;


#include "iheap.h"
#include <queue>	
#include <utility>  

struct CompareBySecond {
	bool operator()(pair<int, int> a, pair<int, int> b)
	{
		return a.second < b.second;
	}
};

struct SeedPrefixRecord {
    int k;
    double coverage_ratio;
    vector<int> seed_set;
    double total_seconds;
};


class InfGraph: public Graph
{
private:
    vector<bool> visit;
    vector<int> visit_mark;
public:
    vector<vector<int>> hyperG;
    vector<vector<int>> hyperGT;

    InfGraph(string folder, string graph_file): Graph(folder, graph_file)
    {
        sfmt_init_gen_rand(&sfmtSeed , 95082);
        init_hyper_graph();
        visit = vector<bool> (n);
        visit_mark = vector<int> (n);
    }


    void init_hyper_graph(){
        hyperG.clear();
        for (int i = 0; i < n; i++)
            hyperG.push_back(vector<int>());
        hyperGT.clear();
    }
    void build_hyper_graph_r(int64 R, const Argument & arg)
    {
        if( R > INT_MAX ){
            cout<<"Error:R too large"<<endl;
            exit(1);
        }
        //INFO("build_hyper_graph_r", R);

        int prevSize = hyperGT.size();
        int64 delta = R - prevSize;
        uint64 buildStart = rdtsc();
        INFO("build_hyper_graph_r begin", prevSize, R, delta);
        while ((int)hyperGT.size() <= R)
            hyperGT.push_back( vector<int>() );



        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            random_number.push_back(  sfmt_genrand_uint32(&sfmtSeed) % origin_n);
        }

        //trying BFS start from same node
        int64 rrProgressStep = max<int64>(1, delta / 20);

        for (int i = prevSize; i < R; i++)
        {
#ifdef CONTINUOUS
            BuildHypergraphNode(random_number[i], i, arg );
#endif
#ifdef DISCRETE
            BuildHypergraphNode(random_number[i], i );
#endif
            int64 completed = (int64)i - prevSize + 1;
            if (completed % rrProgressStep == 0 || i == R - 1)
            {
                double elapsedSec = (rdtsc() - buildStart) / TIMES_PER_SEC;
                INFO("build_hyper_graph_r rr_progress", completed, delta, elapsedSec);
            }
        }


        int totAddedElement = 0;
        int64 edgeProgressStep = max<int64>(1, delta / 20);
        for (int i = prevSize; i < R; i++)
        {
            for (int t : hyperGT[i])
            {
                hyperG[t].push_back(i);
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
            int64 completed = (int64)i - prevSize + 1;
            if (completed % edgeProgressStep == 0 || i == R - 1)
            {
                double elapsedSec = (rdtsc() - buildStart) / TIMES_PER_SEC;
                INFO("build_hyper_graph_r hyperG_progress", completed, delta, totAddedElement, elapsedSec);
            }
        }
        double totalElapsedSec = (rdtsc() - buildStart) / TIMES_PER_SEC;
        INFO("build_hyper_graph_r finish", delta, totAddedElement, totalElapsedSec);
    }

#ifdef DISCRETE
#include "discrete_rrset.h"
#endif
#ifdef CONTINUOUS
#include "continuous_rrset.h"
#endif

    //return the number of edges visited
    deque<int> q;
    sfmt_t sfmtSeed;
	vector<int> seedSet;

	//This is build on Mapped Priority Queue
	double build_seedset(int k)
	{

		priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
		vector<int>coverage(n, 0);

		for (int i = 0; i < origin_n; i++)
		{
			pair<int, int>tep(make_pair(i, (int)hyperG[i].size()));
			heap.push(tep);
			coverage[i] = (int)hyperG[i].size();
		}

		int maxInd;

		long long influence = 0;
		long long numEdge = hyperGT.size();

		// check if an edge is removed
		vector<bool> edgeMark(numEdge, false);
		// check if an node is remained in the heap
		vector<bool> nodeMark(n + 1, false);
		for (int i = 0; i < origin_n; i++)
			nodeMark[i] = true;

		seedSet.clear();
		int targetK = min(k, origin_n);
		while ((int)seedSet.size()<targetK)
		{
			pair<int, int>ele = heap.top();
			heap.pop();
			if (ele.second > coverage[ele.first])
			{
				ele.second = coverage[ele.first];
				heap.push(ele);
				continue;
			}

			maxInd = ele.first;
			vector<int>e = hyperG[maxInd];  
			influence += coverage[maxInd];
			seedSet.push_back(maxInd);
			nodeMark[maxInd] = false;

			for (unsigned int j = 0; j < e.size(); ++j){
				if (edgeMark[e[j]])continue;

				vector<int>nList = hyperGT[e[j]];
				for (unsigned int l = 0; l < nList.size(); ++l){
					if (nodeMark[nList[l]])coverage[nList[l]]--;
				}
				edgeMark[e[j]] = true;
			}
		}
		return 1.0*influence / hyperGT.size();
	}

    vector<SeedPrefixRecord> build_seedset_prefix(const vector<int> &requested_k, double base_seconds = 0.0)
    {
        vector<SeedPrefixRecord> result;
        if (requested_k.empty())
            return result;

        vector<int> wanted = requested_k;
        sort(wanted.begin(), wanted.end());
        wanted.erase(unique(wanted.begin(), wanted.end()), wanted.end());
        for (int &value : wanted)
            value = min(value, origin_n);
        wanted.erase(unique(wanted.begin(), wanted.end()), wanted.end());
        if (wanted.empty())
            return result;

        priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond> heap;
        vector<int> coverage(n, 0);

        for (int i = 0; i < origin_n; i++)
        {
            pair<int, int> tep(make_pair(i, (int)hyperG[i].size()));
            heap.push(tep);
            coverage[i] = (int)hyperG[i].size();
        }

        int maxInd;
        long long influence = 0;
        long long numEdge = hyperGT.size();

        vector<bool> edgeMark(numEdge, false);
        vector<bool> nodeMark(n + 1, false);
        for (int i = 0; i < origin_n; i++)
            nodeMark[i] = true;

        seedSet.clear();
        int targetK = min(wanted.back(), origin_n);
        size_t nextWantedIndex = 0;
        while (nextWantedIndex < wanted.size() && wanted[nextWantedIndex] <= 0)
            nextWantedIndex++;

        uint64 greedyStart = rdtsc();
        while ((int)seedSet.size() < targetK)
        {
            pair<int, int> ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])
            {
                ele.second = coverage[ele.first];
                heap.push(ele);
                continue;
            }

            maxInd = ele.first;
            vector<int> e = hyperG[maxInd];
            influence += coverage[maxInd];
            seedSet.push_back(maxInd);
            nodeMark[maxInd] = false;

            for (unsigned int j = 0; j < e.size(); ++j) {
                if (edgeMark[e[j]]) continue;

                vector<int> nList = hyperGT[e[j]];
                for (unsigned int l = 0; l < nList.size(); ++l) {
                    if (nodeMark[nList[l]]) coverage[nList[l]]--;
                }
                edgeMark[e[j]] = true;
            }

            while (nextWantedIndex < wanted.size() && wanted[nextWantedIndex] == (int)seedSet.size())
            {
                SeedPrefixRecord record;
                record.k = wanted[nextWantedIndex];
                record.coverage_ratio = 1.0 * influence / hyperGT.size();
                record.seed_set = seedSet;
                record.total_seconds = base_seconds + (rdtsc() - greedyStart) / TIMES_PER_SEC;
                result.push_back(record);
                nextWantedIndex++;
            }
        }
        return result;
    }

};


