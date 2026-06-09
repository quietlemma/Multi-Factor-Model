class Math{
    public:
        static double log2(int n){
            return log(n) / log(2);
        }
        static double logcnk(int n, int k) {
            double ans = 0;
            for (int i = n - k + 1; i <= n; i++)
            {
                ans += log(i);
            }
            for (int i = 1; i <= k; i++)
            {
                ans -= log(i);
            }
            return ans;
        }
};

struct PrefixInfluenceEntry {
    int k;
    double estimated_influence;
    vector<int> seed_set;
    double total_seconds;
};

class Imm
{
    private:
        //static InfGraph g;
        //static int k;
        //static map<string, string> arg;

        static double step1(InfGraph &g, const Argument & arg)
        {
            double epsilon_prime = arg.epsilon * sqrt(2);
            int effective_n = g.origin_n;
            Timer t(1, "step1");
            INFO("step1 begin", effective_n, arg.k, arg.epsilon, epsilon_prime);
            for (int x = 1; ; x++)
            {
                int64 ci = (2+2/3 * epsilon_prime)* ( log(effective_n) + Math::logcnk(effective_n, arg.k) + log(Math::log2(effective_n))) * pow(2.0, x) / (epsilon_prime* epsilon_prime);
                INFO("step1 round begin", x, ci, (int64)g.hyperGT.size());
                g.build_hyper_graph_r(ci, arg);
                INFO("step1 rrsets ready", x, (int64)g.hyperGT.size());
                
                INFO("step1 seedset begin", x, arg.k);
				double ept = g.build_seedset(arg.k);
                INFO("step1 seedset finish", x, ept);
                //double estimate_influence = ept * g.n;
                //INFO(x, estimate_influence);
                if (ept > 1 / pow(2.0, x))
                {
                    double OPT_prime = ept * effective_n / (1+epsilon_prime);
                    INFO("step1 round accepted", x, ept, OPT_prime);
                    //INFO("step1", OPT_prime);
                    //INFO("step1", OPT_prime * (1+epsilon_prime));
                    return OPT_prime;
                }
                INFO("step1 round continue", x, ept, 1 / pow(2.0, x));
            }
            ASSERT(false);
            return -1;
        }
        static int64 step2_rr_count(int effective_n, const Argument &arg, double OPT_prime)
        {
            double e = exp(1);
            double alpha = sqrt(log(effective_n) + log(2));
            double beta = sqrt((1-1/e) * (Math::logcnk(effective_n, arg.k) + log(effective_n) + log(2)));
            return 2.0 * effective_n *  sqr((1-1/e) * alpha + beta) /  OPT_prime / arg.epsilon / arg.epsilon;
        }
        static double step2(InfGraph &g, const Argument & arg, double OPT_prime)
        {
            Timer t(2, "step2");
            ASSERT(OPT_prime > 0);
            int effective_n = g.origin_n;
            double e = exp(1);
            double alpha = sqrt(log(effective_n) + log(2));
            double beta = sqrt((1-1/e) * (Math::logcnk(effective_n, arg.k) + log(effective_n) + log(2)));

            int64 R = step2_rr_count(effective_n, arg, OPT_prime);

            //R/=100;
            INFO("step2 begin", effective_n, arg.k, OPT_prime, R);
            INFO("step2 begin detail", alpha, beta, (int64)g.hyperGT.size());
            g.build_hyper_graph_r(R, arg);
            INFO("step2 rrsets ready", R, (int64)g.hyperGT.size());
            
            INFO("step2 seedset begin", arg.k);
			double opt = g.build_seedset(arg.k)*effective_n;
            INFO("step2 seedset finish", arg.k, opt);
            return opt;
        }
        static vector<PrefixInfluenceEntry> step2_with_prefixes(
            InfGraph &g, const Argument &arg, double OPT_prime, const vector<int> &requested_k, uint64 total_start)
        {
            Timer t(2, "step2");
            ASSERT(OPT_prime > 0);
            int effective_n = g.origin_n;
            double e = exp(1);
            double alpha = sqrt(log(effective_n) + log(2));
            double beta = sqrt((1-1/e) * (Math::logcnk(effective_n, arg.k) + log(effective_n) + log(2)));

            int64 R = step2_rr_count(effective_n, arg, OPT_prime);
            INFO("step2 begin", effective_n, arg.k, OPT_prime, R);
            INFO("step2 begin detail", alpha, beta, (int64)g.hyperGT.size());
            g.build_hyper_graph_r(R, arg);
            INFO("step2 rrsets ready", R, (int64)g.hyperGT.size());

            INFO("step2 seedset begin", arg.k);
            double base_seconds = (rdtsc() - total_start) / TIMES_PER_SEC;
            vector<SeedPrefixRecord> prefix_records = g.build_seedset_prefix(requested_k, base_seconds);
            vector<PrefixInfluenceEntry> result;
            for (const SeedPrefixRecord &record : prefix_records)
            {
                PrefixInfluenceEntry entry;
                entry.k = record.k;
                entry.estimated_influence = record.coverage_ratio * effective_n;
                entry.seed_set = record.seed_set;
                entry.total_seconds = record.total_seconds;
                result.push_back(entry);
                INFO("step2 prefix result", entry.k, entry.estimated_influence, entry.total_seconds);
            }
            if (!result.empty())
                INFO("step2 seedset finish", arg.k, result.back().estimated_influence);
            return result;
        }
    public:
        static double InfluenceMaximize(InfGraph &g, const Argument &arg)
        {
            Timer t(100, "InfluenceMaximize(Total Time)");

            INFO("########## Step1 ##########");

            // debugging mode lalala
            double OPT_prime;
            OPT_prime = step1(g, arg ); //estimate OPT



            INFO("########## Step2 ##########");


            double opt_lower_bound = OPT_prime;
            INFO(opt_lower_bound);
            double estimated_influence = step2(g, arg, OPT_prime);
            INFO("step2 finish");
            return estimated_influence;

        }
        static vector<PrefixInfluenceEntry> InfluenceMaximizeWithReuse(
            InfGraph &g, const Argument &arg, const vector<int> &requested_k)
        {
            Timer t(100, "InfluenceMaximize(Total Time)");
            uint64 total_start = rdtsc();

            INFO("########## Step1 ##########");
            double OPT_prime = step1(g, arg);

            INFO("########## Step2 ##########");
            double opt_lower_bound = OPT_prime;
            INFO(opt_lower_bound);

            vector<PrefixInfluenceEntry> result = step2_with_prefixes(g, arg, OPT_prime, requested_k, total_start);
            INFO("step2 finish");
            return result;
        }

};
