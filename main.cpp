/**
* @file main.cpp
* @brief This project implements IMAGE-BR IMAGE:
* @author Song Bian The Chinese University of Hong Kong
*
*/

#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"

using namespace std;

int main(int argc, char* argv[])
{
    // Randomize the seed for generating random numbers
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
    const TArgument Arg(argc, argv);
    const string infilename = Arg._dir + "/" + Arg._graphname;
    if (Arg._func == 0 || Arg._func == 2) {
        // Format the graph
        GraphBase::format_graph(infilename, Arg._mode);
        if (Arg._func == 0) {
            return 1;
        }
    }

    cout << "---The Begin of " << Arg._outFileName << "---\n";
    Timer mainTimer("main");
    cout << "Format Graph Successfully" << endl;

    // Load the reverse graph
    Graph graph = GraphBase::load_graph(infilename, false, Arg._probDist, Arg._probEdge);
    Graph resgraph = GraphBase::load_graph(infilename, true, Arg._probDist, Arg._probEdge);
    cout << "Load Graph Successfully" << endl;

    vector<size_t> vec_degree(graph.size(), 0);
    std::ofstream outputFile("degree_" + Arg._graphname);
    size_t count_degree = 0;
    for (size_t i = 0; i < graph.size(); i++) {
        vec_degree[i] = graph[i].size();
        count_degree = count_degree + graph[i].size();
    }
    std::cout << "# Node: " << graph.size() << std::endl;
    std::cout << "# Edge: " << count_degree << std::endl;
    std::cout << "Ave degree: " << (double) count_degree / graph.size() << std::endl;
    outputFile.close();
    sort(vec_degree.begin(), vec_degree.end(), [](auto & left, auto & right){
        return left < right;
    });

    // Assign budget for each node
    vector<double> budget_list(graph.size(), 0);
    default_random_engine engine;
    uniform_int_distribution<unsigned> u(5, 10);
    for (size_t i = 0; i < graph.size(); i++) {
        if (graph[i].size() > vec_degree[int(Arg._percent * graph.size())]) {
            budget_list[i] = (double) 10 * 0.001 * graph[i].size();
//            budget_list[i] = (double) u(engine) * 0.001 * graph[i].size();
        }
        else {
            budget_list[i] = (double) 10 * Arg._normalpha * graph[i].size();
//            budget_list[i] = (double) u(engine) * Arg._normalpha * graph[i].size();
        }
    }
    double real_budget = (double) Arg._budpercent * int(Arg._percent * graph.size());
//    double real_budget = (double) Arg._budpercent * int(graph.size());
    std::cout << "0.3 percent degree: " << vec_degree[int(0.3 * graph.size())] << std::endl;
    std::cout << "0.5 percent degree: " << vec_degree[int(0.5 * graph.size())] << std::endl;
    std::cout << "0.7 percent degree: " << vec_degree[int(0.7 * graph.size())] << std::endl;
    std::cout << "0.9 percent degree: " << vec_degree[int(0.9 * graph.size())] << std::endl;
    std::cout << "0.95 percent degree: " << vec_degree[int(0.95 * graph.size())] << std::endl;
    std::cout << "0.99 percent degree: " << vec_degree[int(0.99 * graph.size())] << std::endl;
    std::cout << "0.999 percent degree: " << vec_degree[int(0.999 * graph.size())] << std::endl;
    std::cout << "0.9999 percent degree: " << vec_degree[int(0.9999 * graph.size())] << std::endl;
    std::cout << "0.99995 percent degree: " << vec_degree[int(0.99995 * graph.size())] << std::endl;
    std::cout << "0.99999 percent degree: " << vec_degree[int(0.99999 * graph.size())] << std::endl;
    std::cout << "0.999999 percent degree: " << vec_degree[int(0.999999 * graph.size())] << std::endl;
    std::cout << "0.9999999 percent degree: " << vec_degree[int(0.9999999 * graph.size())] << std::endl;
    std::cout << "0.99999999 percent degree: " << vec_degree[int(0.99999999 * graph.size())] << std::endl;
    std::cout << "# nodes: " << int(Arg._percent * graph.size()) << std::endl;
    std::cout << "real budget: " << real_budget << std::endl;

    if (Arg._algName == "CELF" || Arg._algName == "celf" ||
        Arg._algName == "CELF++" || Arg._algName == "celf++") {
        if (Arg._algName == "CELF" || Arg._algName == "celf") {
            cout << "start CELF algorithm" << endl;
            Timer EvalTimer;
            vector<size_t> res = CELF::celfalgo(graph, Arg._budget, budget_list, Arg._model);
            cout << "Time used (sec): " << EvalTimer.get_total_time() << endl;
            for(unsigned int re : res) {
                cout << re << " ";
            }
            cout << endl;
            cout << "Seed Size: " << res.size() << endl;
        }
        else if (Arg._algName == "CELF++" || Arg._algName == "celf++") {
            cout << "start CELF++ algorithm" << endl;
            Timer EvalTimer;
            vector<size_t> res = CELFPP::celfppalgo(graph, Arg._budget, budget_list, Arg._model);
            cout << "Time used (sec): " << EvalTimer.get_total_time() << endl;
            for(unsigned int re : res) {
                cout << re << " ";
            }
            cout << endl;
            cout << "Seed Size: " << res.size() << endl;
        }
    }
    // IMAGE == opim-b-fe IMAGE-BR = opim-b
    else if (Arg._algName == "opim-c" || Arg._algName == "OPIM-C" ||
             Arg._algName == "opim" || Arg._algName == "OPIM" ||
             Arg._algName == "opim-b" || Arg._algName == "OPIM-B" ||
             Arg._algName == "opim-b-n" || Arg._algName == "OPIM-B-N" ||
             Arg._algName == "opim-b-fe" || Arg._algName == "OPIM-B-FE" ||
             Arg._algName == "opim-ba" || Arg._algName == "OPIM-BA" ||
             Arg._algName == "opim-b-fa" || Arg._algName == "OPIM-B-FA" ||
             Arg._algName == "opim-a" || Arg._algName == "OPIM-A" ||
             Arg._algName == "opim-a-f" || Arg._algName == "OPIM-A-F") {
        if (Arg._model == LT) {
            // Normalize the propagation probabilities in accumulation format for
            // LT cascade model for quickly generating RR sets
            to_normal_accum_prob(resgraph);
        }
        // Initialize a result object to record the results
        TResult tRes;
        TAlg tAlg(resgraph, tRes);
        tAlg.set_cascade_model(Arg._model); // Set propagation model

        cout << "  ==>Graph loaded for RIS! total time used (sec): " << mainTimer.get_total_time() << '\n';
        int mode = 2; // Default is to use the minimum upper bound among all the rounds
        if (Arg._mode == "0" || Arg._mode == "vanilla") {
            mode = 0;
        }
        else if (Arg._mode == "1" || Arg._mode == "last") {
            mode = 1;
        }
        auto delta = Arg._delta;
        if (delta < 0) {
            delta = 1.0 / resgraph.size();
        }
        if (Arg._algName == "opim-c" || Arg._algName == "OPIM-C") {
            tAlg.opimc(Arg._seedsize, Arg._eps, delta, mode);
        }
        else if (Arg._algName == "opim" || Arg._algName == "OPIM") {
            tAlg.opim(Arg._seedsize, Arg._samplesize, delta, mode);
        }
        else if (Arg._algName == "opim-b" || Arg._algName == "OPIM-B") {
            cout << "OPIM-B" << endl;
            tAlg.opimb(graph, real_budget, budget_list, Arg._eps, delta, mode);
        }
        else if (Arg._algName == "opim-b-n" || Arg._algName == "OPIM-B-N") {
            cout << "OPIM-B-NAIVE" << endl;
            tAlg.opimb_naive(graph, real_budget, budget_list, Arg._eps, delta, 0);
        }
        else if (Arg._algName == "opim-b-fe" || Arg._algName == "OPIM-B-FE") {
            cout << "OPIM-B-FAST-END" << endl;
            if (Arg._model == LT && Arg._graphname == "twitter") {
                tAlg.opimb_fast_end(graph, real_budget, budget_list, Arg._eps, Arg._feps, delta, "LT",
                        Arg._graphname, 1);
            }
            else {
                if (Arg._model == LT) {
                    tAlg.opimb_fast_end(graph, real_budget, budget_list, Arg._eps, Arg._feps, delta, "LT",
                            Arg._graphname, mode);
                }
                else {
                    tAlg.opimb_fast_end(graph, real_budget, budget_list, Arg._eps, Arg._feps, delta, "IC",
                            Arg._graphname, mode);
                }
            }
        }
        else if (Arg._algName == "opim-ba" || Arg._algName == "OPIM-BA") {
            cout << "OPIM-B-APPROX" << endl;
            tAlg.opimb_approx(graph, real_budget, budget_list, Arg._eps, delta, mode);
        }
        else if (Arg._algName == "opim-b-fa" || Arg._algName == "OPIM-B-FA") {
            cout << "OPIM-B-FAST-APPROX" << endl;
            tAlg.opimb_fast_approx(graph, real_budget, budget_list, Arg._eps, Arg._feps, delta, mode);
        }
        else if (Arg._algName == "opim-a" || Arg._algName == "OPIM-A") {
            cout << "OPIM-APPROX" << endl;
            if (Arg._model == LT) {
                tAlg.opim_approx(graph, real_budget, budget_list, Arg._eps, Arg._feps, delta, "LT",
                                 Arg._graphname, mode);
            }
            else {
                tAlg.opim_approx(graph, real_budget, budget_list, Arg._eps, Arg._feps, delta, "IC",
                                 Arg._graphname, mode);
            }

        }
        else if (Arg._algName == "opim-a-f" || Arg._algName == "OPIM-A-F") {
            cout << "OPIM-APPROX-FAST" << endl;
            if (Arg._model == LT) {
                tAlg.opim_approx_fast(graph, real_budget, budget_list, Arg._eps, Arg._feps, delta, "LT",
                        Arg._graphname, mode);
            }
            else {
                tAlg.opim_approx_fast(graph, real_budget, budget_list, Arg._eps, Arg._feps, delta, "IC",
                        Arg._graphname, mode);
            }
        }
        TIO::write_result(Arg._outFileName, tRes, Arg._resultFolder);
        TIO::write_order_seeds(Arg._outFileName, tRes, Arg._resultFolder, graph);
        if (Arg._algName == "opim-b" || Arg._algName == "OPIM-B" ||
            Arg._algName == "opim-b-n" || Arg._algName == "OPIM-B-N" ||
            Arg._algName == "opim-b-fe" || Arg._algName == "OPIM-B-FE" ||
            Arg._algName == "opim-ba" || Arg._algName == "OPIM-BA" ||
            Arg._algName == "opim-b-fa" || Arg._algName == "OPIM-B-FA" ||
            Arg._algName == "opim-a" || Arg._algName == "OPIM-A" ||
            Arg._algName == "opim-a-f" || Arg._algName == "OPIM-A-F") {
            TIO::write_budget_order_seeds(Arg._outFileName, tRes, Arg._resultFolder, graph);
        }
        cout << "---The End of " << Arg._outFileName << "---\n";
    }

    return 0;
}
