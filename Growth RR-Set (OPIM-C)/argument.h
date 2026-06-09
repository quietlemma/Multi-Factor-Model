#pragma once


class Argument
{
public:
	int _func = 1; // Function parameter. 0: format graph, 1: maximize profit, 2: format graph and then maximize profit.
	int _seedsize = 50; // The number of nodes to be selected. Default is 50.
	size_t _samplesize = 1000; // The number of RR sets to be generated.
	float _probEdge = float(0.1); // For the UNI setting, every edge has the same diffusion probability.
	double _eps = 0.1; // Error threshold 1-1/e-epsilon. 
	double _delta = -1.0; // Failure probability delta. Default is 1/#nodes.
	double _evalEps = 0.01; // Error threshold for the final RR-set influence validation.
	double _evalDelta = 1e-3; // Failure probability for the final RR-set influence validation.
	uint32_t _originN = 0; // Number of origin nodes. 0 means use all graph nodes.
	uint32_t _rngSeed = static_cast<uint32_t>(time(nullptr)); // Random seed for RR-set sampling.
	size_t _maxRRsets = 0; // Fixed total RR budget for OPIM-C. 0 means adaptive doubling.
	size_t _statefulCandidatePoolCap = 0; // 0 means disable candidate pool capping.
	std::vector<int> _kList; // Optional prefix budgets, e.g., 0,5,10,...,50.
	CascadeModel _model = IC; // Cascade models: IC, IC2, IC2Mix, LT, OLO. Default is IC.
	std::string _graphname = "facebook"; // Graph name. Default is "facebook".
	std::string _mode = "g"; // Format graph --> g: graph only [default, using WC], w: with edge property.
							 // OPIM or OPIM-C --> 0: vanilla, 1: last-bound, 2: min-bound [default].  
	std::string _dir = "graphInfo"; // Directory
	std::string _resultFolder = "result"; // Result folder. Default is "test".
	std::string _algName = "OPIM"; // Algorithm. Default is oneHop.
	std::string _probDist = "load"; // Probability distribution for diffusion model. Option: load, WC, TR, UNI. Default is loaded from the file.
	std::string _outFileName; // File name of the result
	
	Argument(int argc, char* argv[])
	{
		std::string param, value;
		for (int ind = 1; ind < argc; ind++)
		{
			if (argv[ind][0] != '-') break;
			std::stringstream sstr(argv[ind]);
			getline(sstr, param, '=');
			getline(sstr, value, '=');
			if (!param.compare("-func")) _func = stoi(value);
			else if (!param.compare("-seedsize")) _seedsize = stoi(value);
			else if (!param.compare("-samplesize")) _samplesize = stoull(value);
			else if (!param.compare("-pedge")) _probEdge = stof(value);
			else if (!param.compare("-eps")) _eps = stod(value);
			else if (!param.compare("-delta")) _delta = stod(value);
			else if (!param.compare("-evaleps")) _evalEps = stod(value);
			else if (!param.compare("-evaldelta")) _evalDelta = stod(value);
			else if (!param.compare("-origin-n")) _originN = stoul(value);
			else if (!param.compare("-rngseed")) _rngSeed = stoul(value);
			else if (!param.compare("-max-rrsets")) _maxRRsets = stoull(value);
			else if (!param.compare("-stateful-candidate-pool-cap")) _statefulCandidatePoolCap = stoull(value);
			else if (!param.compare("-klist"))
			{
				std::stringstream kstr(value);
				std::string item;
				while (std::getline(kstr, item, ','))
				{
					if (!item.empty()) _kList.push_back(stoi(item));
				}
			}
			else if (!param.compare("-model"))
			{
				std::string modelName = value;
				for (auto& ch : modelName) ch = static_cast<char>(toupper(ch));
				if (modelName == "LT") _model = LT;
				else if (modelName == "IC2") _model = IC2;
				else if (modelName == "IC2MIX") _model = IC2MIX;
				else if (modelName == "OLO") _model = OLO;
				else _model = IC;
			}
			else if (!param.compare("-gname")) _graphname = value;
			else if (!param.compare("-mode")) _mode = value;
			else if (!param.compare("-dir")) _dir = value;
			else if (!param.compare("-outpath")) _resultFolder = value;
			else if (!param.compare("-alg")) _algName = value;
			else if (!param.compare("-pdist")) _probDist = value;
		}
		if (!_kList.empty())
		{
			std::sort(_kList.begin(), _kList.end());
			_kList.erase(std::remove(_kList.begin(), _kList.end(), -1), _kList.end());
			_kList.erase(std::unique(_kList.begin(), _kList.end()), _kList.end());
			_kList.erase(std::remove_if(_kList.begin(), _kList.end(), [](int k) { return k < 0; }), _kList.end());
			if (!_kList.empty()) _seedsize = _kList.back();
		}
		std::string postfix = "_minBound"; // Default is to use the minimum upper bound among all the rounds
		if (_mode == "0" || _mode == "vanilla")
		{
			postfix = "_vanilla";
		}
		else if (_mode == "1" || _mode == "last")
		{
			postfix = "_lastBound";
		}
		_outFileName = TIO::get_out_file_name(_graphname, _algName + postfix, _seedsize, _probDist, _probEdge);
		if (_model == LT) _outFileName = "LT_" + _outFileName;
		else if (_model == IC2) _outFileName = "IC2_" + _outFileName;
		else if (_model == IC2MIX) _outFileName = "IC2Mix_" + _outFileName;
		else if (_model == OLO) _outFileName = "OLO_" + _outFileName;
		if (_algName == "OPIM" || _algName == "opim") _outFileName += "_s" + std::to_string(_samplesize);
		if (_maxRRsets > 0) _outFileName += "_rr" + std::to_string(_maxRRsets);
		if (_statefulCandidatePoolCap > 0) _outFileName += "_cap" + std::to_string(_statefulCandidatePoolCap);
	}

	std::string get_outfilename_with_alg(const std::string& algName) const
	{
		return TIO::get_out_file_name(_graphname, algName, _seedsize, _probDist, _probEdge);
	}
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;
