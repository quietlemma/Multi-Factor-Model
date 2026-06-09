/**
* @file OPIM.cpp
* @brief This project implements the OPIM and OPIM-C for the following paper:
* Jing Tang, Xueyan Tang, Xiaokui Xiao, Junsong Yuan, "Online Processing Algorithms for Influence Maximization," in Proc. ACM SIGMOD, 2018.
*
* @author Jing Tang (Nanyang Technological University)
*
* Copyright (C) 2018 Jing Tang and Nanyang Technological University. All rights reserved.
*
*/

#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"

int main(int argc, char* argv[])
{
	const TArgument Arg(argc, argv);
	// Randomize the seed for generating random numbers. Use -rngseed for reproducible experiments.
	dsfmt_gv_init_gen_rand(Arg._rngSeed);
	const std::string infilename = Arg._dir + "/" + Arg._graphname;
	if (Arg._func == 0 || Arg._func == 2)
	{
		// Format the graph
		GraphBase::format_graph(infilename, Arg._mode);
		if (Arg._func == 0) return 1;
	}

	std::cout << "---The Begin of " << Arg._outFileName << "---\n";
	LogProgress("opim_process_start out=" + Arg._outFileName
	            + " graph=" + Arg._graphname
	            + " alg=" + Arg._algName
	            + " k=" + std::to_string(Arg._seedsize));
	Timer mainTimer("main");
	
	// Load the reverse graph
	Graph graph = GraphBase::load_graph(infilename, true, Arg._probDist, Arg._probEdge);
	if (Arg._model == LT)
	{
		// Normalize the propagation probabilities in accumulation format for LT cascade model for quickly generating RR sets
		to_normal_accum_prob(graph);
	}
	print_memory_snapshot("after_graph_load");
	LogProgress("graph_loaded nodes=" + std::to_string(graph.size()));
	// Initialize a result object to record the results
	TResult tRes;
	TAlg tAlg(graph, tRes);
	tAlg.set_cascade_model(Arg._model); // Set propagation model
	tAlg.set_origin_node_count(Arg._originN); // Origin-only mode for auxiliary graphs when provided.
	tAlg.set_fixed_rrsets_total(Arg._maxRRsets);
	tAlg.set_stateful_candidate_pool_cap(Arg._statefulCandidatePoolCap);

	std::cout << "  ==>Graph loaded for RIS! total time used (sec): " << mainTimer.get_total_time() << '\n';
	int mode = 2; // Default is to use the minimum upper bound among all the rounds
	if (Arg._mode == "0" || Arg._mode == "vanilla")
	{
		mode = 0;
	}
	else if (Arg._mode == "1" || Arg._mode == "last")
	{
		mode = 1;
	}
	auto delta = Arg._delta;
	if (delta < 0) delta = 1.0 / graph.size();
	if (Arg._algName == "opim-c" || Arg._algName == "OPIM-C")
	{
		if (!Arg._kList.empty())
		{
			tAlg.opimc_klist(Arg._kList, Arg._eps, delta, mode, Arg._evalDelta, Arg._evalEps,
			                 Arg._outFileName, Arg._resultFolder);
		}
		else
		{
			tAlg.opimc(Arg._seedsize, Arg._eps, delta, mode, Arg._evalDelta, Arg._evalEps);
		}
	}
	else if (Arg._algName == "opim" || Arg._algName == "OPIM")
	{
		if (!Arg._kList.empty())
		{
			tAlg.opim_klist(Arg._kList, Arg._samplesize, delta, mode, Arg._evalDelta, Arg._evalEps,
			                Arg._outFileName, Arg._resultFolder);
		}
		else
		{
			tAlg.opim(Arg._seedsize, Arg._samplesize, delta, mode, Arg._evalDelta, Arg._evalEps);
		}
	}
	TIO::write_result(Arg._outFileName, tRes, Arg._resultFolder);
	TIO::write_order_seeds(Arg._outFileName, tRes, Arg._resultFolder);
	LogProgress("opim_process_done out=" + Arg._outFileName
	            + " result_dir=" + Arg._resultFolder);
	std::cout << "---The End of " << Arg._outFileName << "---\n";
	return 0;
}
