//
// Created by sbian on 2019/9/9.
//

#ifndef BIM_ALG_H
#define BIM_ALG_H

class Alg
{
private:
    /// __numV: number of nodes in the graph.
    uint32_t __numV;
    /// __numE: number of edges in the graph.
    size_t __numE;
    /// __numRRsets: number of RR sets.
    size_t __numRRsets = 0;
    /// flag choose upper bound
    bool __is_inf_cost = false;
    /// Upper bound in the last round for __mode=1.
    double __boundLast = DBL_MAX;
    /// The minimum upper bound among all rounds for __model=2.
    double __boundMin = DBL_MAX;
    /// Upper bound in the last round for __mode=1 influence first.
    double __boundLast_inf = DBL_MAX;
    /// The minimum upper bound among all rounds for __model=2 influence first.
    double __boundMin_inf = DBL_MAX;
    /// Upper bound in the last round for __mode=1 influence/cost first.
    double __boundLast_inf_cost = DBL_MAX;
    /// The minimum upper bound among all rounds for __model=2 influence/cost first.
    double __boundMin_inf_cost = DBL_MAX;
    /// Two hyper-graphs, one is used for selecting seeds and the other is used for validating influence.
    THyperGraph __hyperG, __hyperGVldt;
    /// Result object.
    TResult& __tRes;
    /// Seed set.
    Nodelist __vecSeed;
    /// Seed set influence first
    Nodelist __vecSeed_inf;
    /// Seed set influence/cost first
    Nodelist __vecSeed_inf_cost;
    /// Maximum coverage by lazy updating.
    double max_cover_lazy(const int targetSize, const int mode = 2);
    /// Maximum coverage by maintaining the top-k marginal coverage.
    double max_cover_topk(const int targetSize);
    /// Maximum coverage by lazy updating budget (reverse index).
    double max_cover_lazy_budget(const Graph & graph, const double budget,
            const std::vector<double> & budget_list, double & timeb, const int mode = 2);
    /// Maximum coverage by lazy updating budget fast (reverse index). End version
    double max_cover_lazy_budget_fast_end(const Graph & graph, const double budget,
            const std::vector<double> & budget_list, const double epsilon, double & timeb, const int mode = 2);
    /// Maximum coverage.
    double max_cover(const int targetSize, const int mode = 2);
    /// Maximum coverage budget (reverse index).
    double max_cover_budget(const Graph & graph, const double budget,
            const std::vector<double> & budget_list, double & timeb, const int mode = 2);
    /// Maximum coverage budget fast end (reverse index).
    double max_cover_budget_fast_end(const Graph & graph, const double budget,
            const std::vector<double> & budget_list, const double epsilon, double & timeb, const int mode);
    /// Maximum element of vector.
    static double max_vector(const std::vector<double> & budget_list);
    /// Maximum element of vector.
    static double max_vector_nonzero(const std::vector<double> & budget_list);
    /// Mimimum element of vector.
    static double min_vector(const std::vector<double> & budget_list);
    /// Miminum upper bound
    static double min_bound(const double a, const double b);
    /// Mimimum element of vector not zero.
    static double min_vector_nonzero(const std::vector<double> & budget_list);
    /// Minimum element of vector not zero with mask.
    static double min_vector_nonzero_mask(const std::vector<double> & budget_list,
            const std::vector<bool> & mask);
    /// Max element
    static double max_element(double a, double b);
    /// return k_max
    static int cal_k_max(std::vector<double> budget_list, double budget);
    /// return k_min
    static int cal_k_min(std::vector<double> budget_list, double budget);
    /// ProductLog Function
    static double productLog(double realz, double imagez);
public:
    Alg(const Graph& graph, TResult& tRes) : __hyperG(graph), __hyperGVldt(graph), __tRes(tRes)
    {
        __numV = __hyperG.get_nodes();
        __numE = __hyperG.get_edges();
    }
    ~Alg()
    {
    }
    /// Set cascade model.
    void set_cascade_model(const CascadeModel model);
    /// Evaluate influence spread for the seed set constructed
    double effic_inf_valid_algo();
    /// Evaluate influence spread for a given seed set
    double effic_inf_valid_algo(const Nodelist vecSeed);
    /// OPIM: return the approximation guarantee alpha for the greedy algorithm when a given number of RR sets are generated.
    double opim(const int targetSize, const size_t numRRsets, const double delta, const int mode = 2);
    /// OPIM-C: return (epsilon, delta)-approximate solution for influence maximization.
    double opimc(const int targetSize, const double epsilon, const double delta, const int mode = 2);
    /// OPIM-B (reverse index): return (epsilon, delta)-approximate solution for budget influence maximization.
    double opimb(const Graph & graph, const double budget, const std::vector<double> & budget_list,
            const double epsilon, const double delta, const int mode = 2);
    /// OPIM-B-APPROX (reverse index): return (epsilon, delta)-approximate solution for budget influence maximization.
    double opimb_approx(const Graph & graph, const double budget, const std::vector<double> & budget_list,
            const double epsilon, const double delta, const int mode = 2);
    /// OPIM-APPROX (reverse index): return approximation under same # RR sets
    double opim_approx(const Graph & graph, const double budget, const std::vector<double> & budget_list,
            const double epsilon, const double feps, const double delta, const std::string model,
            const std::string graphname, const int mode = 2);
    /// OPIM-APPROX-FAST (reverse index): return approximation under same # RR sets however feps 0.01 0.05 0.10
    double opim_approx_fast(const Graph & graph, const double budget, const std::vector<double> & budget_list,
            const double epsilon, const double feps, const double delta, const std::string model,
            const std::string graphname, const int mode = 2);
    /// OPIM-B (reverse index without proposed upper bound)
    double opimb_naive(const Graph & graph, const double budget, const std::vector<double> & budget_list,
            const double epsilon, const double delta, const int mode = 0);
    /// OPIM-B-FAST-END (reverse index): return (epsilon, delta)-approximate solution for budget influence maximization.
    double opimb_fast_end(const Graph & graph, const double budget, const std::vector<double> & budget_list,
            const double epsilon, const double feps, const double delta, const std::string model,
            std::string graphname, const int mode = 2);
    /// OPIM-B-FAST-APPROX (reverse index): return (epsilon, delta)-approximate solution for budget influence maximization.
    double opimb_fast_approx(const Graph & graph, const double budget, const std::vector<double> & budget_list,
            const double epsilon, const double feps, const double delta, const int mode = 2);
};

using TAlg = Alg;
using PAlg = std::shared_ptr<TAlg>;

#endif //BIM_ALG_H
