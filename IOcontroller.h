//
// Created by sbian on 2019/9/5.
//

#ifndef BIM_IOCONTROLLER_H
#define BIM_IOCONTROLLER_H

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

class IOcontroller
{
public:
    static void mkdir_absence(const char* outFolder)
    {
#if defined(_WIN32)
        CreateDirectoryA(outFolder, nullptr); // can be used on Windows
#else
        mkdir(outFolder, 0733); // can be used on non-Windows
#endif
    }

    /// Save a serialized file
    template <class T>
    static void save_file(const std::string filename, const T& output)
    {
        std::ofstream outfile(filename, std::ios::binary);
        if (!outfile.eof() && !outfile.fail())
        {
            StreamType res;
            serialize(output, res);
            outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
            outfile.close();
            res.clear();
            std::cout << "Save file successfully: " << filename << '\n';
        }
        else
        {
            std::cout << "Save file failed: " + filename << '\n';
            exit(1);
        }
    }

    /// Load a serialized file
    template <class T>
    static void load_file(const std::string filename, T& input)
    {
        std::ifstream infile(filename, std::ios::binary);
        if (!infile.eof() && !infile.fail())
        {
            infile.seekg(0, std::ios_base::end);
            const std::streampos fileSize = infile.tellg();
            infile.seekg(0, std::ios_base::beg);
            std::vector<uint8_t> res(fileSize);
            infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
            infile.close();
            input.clear();
            auto it = res.cbegin();
            input = deserialize<T>(it, res.cend());
            res.clear();
        }
        else
        {
            std::cout << "Cannot open file: " + filename << '\n';
            exit(1);
        }
    }

    /// Save graph structure to a file
    static void save_graph_struct(const std::string graphName, const Graph& vecGraph, const bool isReverse)
    {
        std::string postfix = ".vec.graph";
        if (isReverse) postfix = ".vec.rvs.graph";
        const std::string filename = graphName + postfix;
        save_file(filename, vecGraph);
    }

    /// Load graph structure from a file
    static void load_graph_struct(const std::string graphName, Graph& vecGraph, const bool isReverse)
    {
        std::string postfix = ".vec.graph";
        if (isReverse) postfix = ".vec.rvs.graph";
        const std::string filename = graphName + postfix;
        load_file(filename, vecGraph);
    }

    /// Get out-file name
    static std::string get_out_file_name(const std::string graphName, const std::string algName, const int seedsize,
                                         const std::string probDist, const float probEdge)
    {
        if (probDist == "UNI")
        {
            return graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist + std::
            to_string(probEdge);
        }
        return graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist;
    }

    /// New methods get out-file name
    static std::string get_out_file_name_budget(const std::string graphName, const std::string algName,
            const double budpercent, const double normalpha, const double percent, const double eps,
            const double feps, const std::string probDist, const float probEdge) {
        if (probDist == "UNI") {
            return graphName + "_" + algName + "_" + std::to_string(budpercent) +
                   "_" + std::to_string(normalpha) + "_" + std::to_string(percent) +
                   "_" + std::to_string(eps) + "_" + std::to_string(feps) +
                   "_" + probDist + std::to_string(probEdge);
        }
        return graphName + "_" + algName + "_" + std::to_string(budpercent) +
               "_" + std::to_string(normalpha) + "_" + std::to_string(percent) +
               "_" + std::to_string(eps) + "_" + std::to_string(feps) + "_" + probDist;
    }

    /// Print the results
    static void write_result(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder)
    {
        const auto approx = resultObj.get_approximation();
        const auto runTime = resultObj.get_running_time();
        const auto influence = resultObj.get_influence();
        const auto influenceOriginal = resultObj.get_influence_original();
        const auto seedSize = resultObj.get_seed_size();
        const auto RRsetsSize = resultObj.get_RRsets_size();
        const auto Boundmin_influence = resultObj.get_boundmin_inf();
        const auto Boundmin_cost_influence = resultObj.get_boundmin_inf_cost();

        std::cout << "   --------------------" << std::endl;
        std::cout << "  |Approx.: " << approx << std::endl;
        std::cout << "  |Time (sec): " << runTime << std::endl;
        std::cout << "  |Influence: " << influence << std::endl;
        std::cout << "  |Self-estimated influence: " << influenceOriginal << std::endl;
        std::cout << "  |#Seeds: " << seedSize << std::endl;
        std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
        std::cout << "  |BoundMin Influence: " << Boundmin_influence << std::endl;
        std::cout << "  |BoundMinCost Influence: " << Boundmin_cost_influence << std::endl;
        std::cout << "   --------------------" << std::endl;
        mkdir_absence(outFolder.c_str());
        std::ofstream outFileNew(outFolder + "/" + outFileName);
        if (outFileNew.is_open())
        {
            outFileNew << "Approx.: " << approx << std::endl;
            outFileNew << "Time (sec): " << runTime << std::endl;
            outFileNew << "Influence: " << influence << std::endl;
            outFileNew << "Self-estimated influence: " << influenceOriginal << std::endl;
            outFileNew << "#Seeds: " << seedSize << std::endl;
            outFileNew << "#RR sets: " << RRsetsSize << std::endl;
            outFileNew << "BoundMin Influence: " << Boundmin_influence << std::endl;
            outFileNew << "BoundMinCost Influence: " << Boundmin_cost_influence << std::endl;
            outFileNew.close();
        }
    }

    /// Print the seeds
    static void write_order_seeds(const std::string& outFileName, const TResult& resultObj,
            const std::string& outFolder, const Graph & graph)
    {
        auto vecSeed = resultObj.get_seed_vec();
        mkdir_absence(outFolder.c_str());
        const auto outpath = outFolder + "/seed";
        mkdir_absence(outpath.c_str());
        std::ofstream outFile(outpath + "/seed_" + outFileName);
        for (auto i = 0; i < vecSeed.size(); i++)
        {
            outFile << vecSeed[i] << " " << graph[i].size() << '\n';
        }
        outFile.close();
    }

    /// Print budget the budget seeds
    static void write_budget_order_seeds(const std::string & outFileName, const TResult & resultObj,
            const std::string & outFolder, const Graph & graph) {
        auto vecSeed_inf = resultObj.get_seed_inf_vec();
        auto vecSeed_inf_cost = resultObj.get_seed_inf_cost_vec();
        mkdir_absence(outFolder.c_str());
        const auto outpath = outFolder + "/seed";
        mkdir_absence(outpath.c_str());
        std::ofstream outFile_inf(outpath + "/seed_" + outFileName + "_inf");
        for (auto i = 0; i < vecSeed_inf.size(); i++)
        {
            outFile_inf << vecSeed_inf[i] << " " << graph[i].size() << '\n';
        }
        outFile_inf.close();
        std::ofstream outFile_inf_cost(outpath + "/seed_" + outFileName + "_inf_cost");
        for (auto i = 0; i < vecSeed_inf_cost.size(); i++)
        {
            outFile_inf_cost << vecSeed_inf_cost[i] << " " << graph[i].size() << '\n';
        }
        outFile_inf_cost.close();
    }
};

using TIO = IOcontroller;
using PIO = std::shared_ptr<IOcontroller>;

#endif //BIM_IOCONTROLLER_H
