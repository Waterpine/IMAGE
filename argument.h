//
// Created by sbian on 2019/9/3.
//

#ifndef BIM_ARGUMENT_H
#define BIM_ARGUMENT_H

class Argument
{
public:
    int _func = 2; // Function parameter. 0: format graph, 1: maximize profit, 2: format graph and then maximize profit.
    int _seedsize = 50; // The number of nodes to be selected. Default is 50.
    size_t _samplesize = 1000; // The number of RR sets to be generated.
    float _probEdge = float(0.1); // For the UNI setting, every edge has the same diffusion probability.
    double _budget = 5000;  // The budget is assigned for influence maximization
    double _eps = 0.1; // Error threshold 1-1/e-epsilon.
    double _feps = 0.1; // fast greedy algorithm epsilon
    double _delta = -1.0; // Failure probability delta. Default is 1/#nodes.
    double _percent = 0.99999; // The _percent nodes which have least degree's node cost will be set as 1 * degree
    double _budpercent = 0.0001; // budget percent 0.0001 is default, if larger than 0.0001 too many rr-sets
    double _normalpha = 0.2; // 0.99999 people node cost alpha
    CascadeModel _model = IC; // Cascade models: IC, LT. Default is IC.
    std::string _graphname = "facebook"; // Graph name. Default is "facebook".
    std::string _mode = "g"; // Format graph --> g: graph only [default, using WC], w: with edge property.
    // OPIM or OPIM-C --> 0: vanilla, 1: last-bound, 2: min-bound [default].
    std::string _dir = "graphInfo"; // Directory
    std::string _resultFolder = "result"; // Result folder. Default is "test".
    std::string _algName = "opim-b"; // Algorithm. Default is oneHop.
    std::string _probDist = "load"; // Probability distribution for diffusion model. Option: load, WC, TR, UNI.
    // Default is loaded from the file.
    std::string _outFileName; // File name of the result

    Argument(int argc, char* argv[])
    {
        std::string param, value;
        for (int ind = 1; ind < argc; ind++)
        {
            if (argv[ind][0] != '-') {
                break;
            }
            std::stringstream sstr(argv[ind]);
            getline(sstr, param, '=');
            getline(sstr, value, '=');
            if (!param.compare("-func")) {
                _func = stoi(value);
            }
            else if (!param.compare("-seedsize")) {
                _seedsize = stoi(value);
            }
            else if (!param.compare("-samplesize")) {
                _samplesize = stoull(value);
            }
            else if (!param.compare("-pedge")) {
                _probEdge = stof(value);
            }
            else if (!param.compare("-eps")) {
                _eps = stod(value);
            }
            else if (!param.compare("-feps")) {
                _feps = stod(value);
            }
            else if (!param.compare("-delta")) {
                _delta = stod(value);
            }
            else if (!param.compare("-budget")) {
                _budget = stod(value);
            }
            else if (!param.compare("-budper")) {
                _budpercent = stod(value);
            }
            else if (!param.compare("-nalpha")) {
                _normalpha = stod(value);
            }
            else if (!param.compare("-model")) {
                _model = value == "LT" ? LT : IC;
            }
            else if (!param.compare("-gname")) {
                _graphname = value;
            }
            else if (!param.compare("-mode")) {
                _mode = value;
            }
            else if (!param.compare("-dir")) {
                _dir = value;
            }
            else if (!param.compare("-outpath")) {
                _resultFolder = value;
            }
            else if (!param.compare("-alg")) {
                _algName = value;
            }
            else if (!param.compare("-pdist")) {
                _probDist = value;
            }
        }
        std::string postfix = "_minBound"; // Default is to use the minimum upper bound among all the rounds
        if (_mode == "0" || _mode == "vanilla") {
            postfix = "_vanilla";
        }
        else if (_mode == "1" || _mode == "last") {
            postfix = "_lastBound";
        }
        if (_algName == "opim-b" || _algName == "OPIM-B" ||
            _algName == "opim-b-n" || _algName == "OPIM-B-N" ||
            _algName == "opim-b-fe" || _algName == "OPIM-B-FE" ||
            _algName == "opim-ba" || _algName == "OPIM-BA" ||
            _algName == "opim-b-fa" || _algName == "OPIM-B-FA" ||
            _algName == "opim-a" || _algName == "OPIM-A" ||
            _algName == "opim-a-f" || _algName == "OPIM-A-F") {
            _outFileName = TIO::get_out_file_name_budget(_graphname, _algName + postfix, _budpercent,
                    _normalpha, _percent, _eps, _feps, _probDist, _probEdge);
        }
        else {
            _outFileName = TIO::get_out_file_name(_graphname, _algName + postfix, _seedsize, _probDist, _probEdge);
        }
        if (_model == LT) {
            _outFileName = "LT_" + _outFileName;
        }
        if (_algName == "OPIM" || _algName == "opim") {
            _outFileName += "_s" + std::to_string(_samplesize);
        }
    }

    std::string get_outfilename_with_alg(const std::string& algName) const
    {
        return TIO::get_out_file_name(_graphname, algName, _seedsize, _probDist, _probEdge);
    }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;

#endif //BIM_ARGUMENT_H
