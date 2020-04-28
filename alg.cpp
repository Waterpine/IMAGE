//
// Created by sbian on 2019/9/9.
//

#include "stdafx.h"

double Alg::max_cover_lazy(const int targetSize, const int mode)
{
    // mode: optimization mode.
    // 0->no optimization,
    // 1->optimization with the upper bound in the last round,
    // 2->optimization with minimum upper bound among all rounds [Default].
    __boundLast = DBL_MAX, __boundMin = DBL_MAX;
    FRset coverage(__numV, 0);
    size_t maxDeg = 0;
    for (auto i = __numV; i--;)
    {
        const auto deg = __hyperG._FRsets[i].size();
        coverage[i] = deg;
        if (deg > maxDeg) maxDeg = deg;
    }
    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
    for (auto i = __numV; i--;)
    {
        if (coverage[i] == 0) continue;
        degMap[coverage[i]].push_back(i);
    }
    size_t sumInf = 0;

    // check if an edge is removed
    std::vector<bool> edgeMark(__numRRsets, false);

    __vecSeed.clear();
    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& vecNode = degMap[deg];
        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];
            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }
            if (mode == 2 || (mode == 1 && __vecSeed.size() == targetSize))
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                Nodelist vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;
                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }
                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();
                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }
                make_min_heap(vecBound);

                // Find the top-k marginal coverage
                auto flag = topk == 0;
                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];
                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        min_heap_replace_min_value(vecBound, currDegBound);
                    }
                }
                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();
                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];
                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            min_heap_replace_min_value(vecBound, currDegBound);
                        }
                    }
                }
                __boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) *
                        __numV / __numRRsets;
                if (__boundMin > __boundLast) __boundMin = __boundLast;
            }
            if (__vecSeed.size() >= targetSize)
            {
                // Top-k influential nodes constructed
                const auto finalInf = 1.0 * sumInf * __numV / __numRRsets;
                std::cout << "  >>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << __boundMin <<
                          ", last-bound: " << __boundLast << '\n';
                return finalInf;
            }
            sumInf = sumInf + currDeg;
            __vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;
            for (auto edgeIdx : __hyperG._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;
                edgeMark[edgeIdx] = true;
                for (auto nodeIdx : __hyperG._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
                    coverage[nodeIdx]--;
                }
            }
        }
        degMap.pop_back();
    }
    return 1.0 * __numV; // All RR sets are covered.
}

double Alg::max_cover_topk(const int targetSize)
{
    FRset coverage(__numV, 0);
    size_t maxDeg = 0;
    for (auto i = __numV; i--;)
    {
        const auto deg = __hyperG._FRsets[i].size();
        coverage[i] = deg;
        if (deg > maxDeg) maxDeg = deg;
    }
    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
    for (auto i = __numV; i--;)
    {
        //if (coverage[i] == 0) continue;
        degMap[coverage[i]].push_back(i);
    }
    Nodelist sortedNode(__numV); // sortedNode: record the sorted nodes in ascending order of degree
    Nodelist nodePosition(__numV); // nodePosition: record the position of each node in the sortedNode
    Nodelist degreePosition(maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
    uint32_t idxSort = 0;
    size_t idxDegree = 0;
    for (auto& nodes : degMap)
    {
        degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
        idxDegree++;
        for (auto& node : nodes)
        {
            nodePosition[node] = idxSort;
            sortedNode[idxSort++] = node;
        }
    }
    // check if an edge is removed
    std::vector<bool> edgeMark(__numRRsets, false);
    // record the total of top-k marginal gains
    size_t sumTopk = 0;
    for (auto deg = maxDeg + 1; deg--;)
    {
        if (degreePosition[deg] <= __numV - targetSize)
        {
            sumTopk += deg * (degreePosition[deg + 1] - (__numV - targetSize));
            break;
        }
        sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
    }
    __boundMin = 1.0 * sumTopk;
    __vecSeed.clear();
    size_t sumInf = 0;
    /*
    * sortedNode: position -> node
    * nodePosition: node -> position
    * degreePosition: degree -> position (start position of this degree)
    * coverage: node -> degree
    * e.g., swap the position of a node with the start position of its degree
    * swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
    */
    for (auto k = targetSize; k--;)
    {
        const auto seed = sortedNode.back();
        sortedNode.pop_back();
        const auto newNumV = sortedNode.size();
        sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
        sumInf += coverage[seed];
        __vecSeed.push_back(seed);
        coverage[seed] = 0;
        for (auto edgeIdx : __hyperG._FRsets[seed])
        {
            if (edgeMark[edgeIdx]) continue;
            edgeMark[edgeIdx] = true;
            for (auto nodeIdx : __hyperG._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
                const auto currPos = nodePosition[nodeIdx]; // The current position
                const auto currDeg = coverage[nodeIdx]; // The current degree
                const auto startPos = degreePosition[currDeg]; // The start position of this degree
                const auto startNode = sortedNode[startPos]; // The node with the start position
                // Swap this node to the start position with the same degree, and update their positions in nodePosition
                std::swap(sortedNode[currPos], sortedNode[startPos]);
                nodePosition[nodeIdx] = startPos;
                nodePosition[startNode] = currPos;
                // Increase the start position of this degree by 1, and decrease the degree of this node by 1
                degreePosition[currDeg]++;
                coverage[nodeIdx]--;
                // If the start position of this degree is in top-k, reduce topk by 1
                if (startPos >= newNumV - targetSize) sumTopk--;
            }
        }
        __boundLast = 1.0 * (sumInf + sumTopk);
        if (__boundMin > __boundLast) __boundMin = __boundLast;
    }
    __boundMin *= 1.0 * __numV / __numRRsets;
    __boundLast *= 1.0 * __numV / __numRRsets;
    const auto finalInf = 1.0 * sumInf * __numV / __numRRsets;
    std::cout << "  >>>[greedy-topk] influence: " << finalInf << ", min-bound: " << __boundMin <<
              ", last-bound: " << __boundLast << '\n';
    return finalInf;
}

double Alg::max_cover_lazy_budget(const Graph & graph, const double budget,
        const std::vector<double> & budget_list, double & timeb, const int mode) {
    // mode: optimization mode.
    // 0->no optimization,
    // 1->optimization with the upper bound in the last round,
    // 2->optimization with minimum upper bound among all rounds [Default].
    const double EPS = 0.0000000001;
    __boundLast = DBL_MAX, __boundMin = DBL_MAX;
    __boundLast_inf = DBL_MAX, __boundMin_inf = DBL_MAX;
    __boundLast_inf_cost = DBL_MAX, __boundMin_inf_cost = DBL_MAX;
    FRset coverage(__numV, 0);
    size_t maxDeg = 0;
    for (auto i = __numV; i--;) {
        const size_t deg = __hyperG._FRsets[i].size();
        coverage[i] = deg;
        if (deg > maxDeg) {
            maxDeg = deg;
        }
    }
    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
    for (auto i = __numV; i--;) {
        if (coverage[i] == 0) {
            continue;
        }
        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    int flag = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(__numRRsets, false);
    // check if an edge is removed
    std::vector<bool> node_sel(__numV, false);
    double min_cost = min_vector_nonzero_mask(budget_list, node_sel);
    std::cout << "min_cost: " << min_cost << std::endl;

    __vecSeed_inf.clear();
    for (size_t deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        if (min_cost > budget) {
            break;
        }
        auto & vecNode = degMap[deg];
        for (auto idx = vecNode.size(); idx--;) {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];
            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }
            if (budget_list[argmaxIdx] <= budget) {
                if (mode == 2) {
                    // Find upper bound
                    auto degBound = deg;
                    Nodelist vecBound;
                    // Initialize vecBound
                    auto idxBound = idx;
                    vecBound.push_back(coverage[degMap[degBound][idxBound]]);
                    __boundLast_inf = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0))
                                             + sumInf) * __numV /__numRRsets;
                    if (__boundMin_inf > __boundLast_inf) {
                        __boundMin_inf = __boundLast_inf;
                    }
                }
                if (flag == 0) {
                    sumInf = sumInf + coverage[argmaxIdx];
                    __vecSeed_inf.push_back(argmaxIdx);
                    coverage[argmaxIdx] = 0;
                    for (auto edgeIdx : __hyperG._FRsets[argmaxIdx])
                    {
                        if (edgeMark[edgeIdx]) continue;
                        edgeMark[edgeIdx] = true;
                        for (auto nodeIdx : __hyperG._RRsets[edgeIdx])
                        {
                            if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
                            coverage[nodeIdx]--;
                        }
                    }
                }
                flag = flag + 1;
            }
            if (flag == 2) {
                break;
            }
        }
        if (flag == 2) {
            break;
        }
        degMap.pop_back();
    }
    // Top-k influential nodes constructed
    const auto finalInf = 1.0 * sumInf * __numV / __numRRsets;
    std::cout << std::endl;
    std::cout << "  >>>[greedy-lazy] influence: " << sumInf << " " << finalInf <<
              ", min-bound: " << __boundMin_inf <<
              ", last-bound: " << __boundLast_inf << '\n';

    // cost influence most
    std::vector<double> coverage_cost(__numV, 0);
    size_t maxDeg_cost = 0;
    for (auto i = __numV; i--;) {
        const size_t deg_cost = __hyperG._FRsets[i].size();
        if (fabs(budget_list[i]) <= EPS || deg_cost == 0) {
            coverage_cost[i] = 0;
        }
        else {
            coverage_cost[i] = (double) deg_cost / budget_list[i];
        }
    }

    // Normalization
    for (auto i = __numV; i--;) {
        if (int(coverage_cost[i]) + 1 > maxDeg_cost) {
            maxDeg_cost = int(coverage_cost[i]) + 1;
        }
    }
    std::cout << "maxDeg: " << maxDeg_cost << std::endl;

    // int(u.inf/u.cost) -> sort(pair<u.id, u.inf/u.cost>)
    std::vector<std::vector<std::pair<size_t, double>>> degMap_cost(maxDeg_cost + 1);
    for (auto i = __numV; i--;) {
        if (fabs(coverage_cost[i]) <= EPS) {
            continue;
        }
        degMap_cost[int(coverage_cost[i]) + 1].push_back(std::make_pair(i, coverage_cost[i]));
    }

    size_t sumInf_cost = 0;
    double sumBudget_cost = 0;
    // label miss nodes
    std::vector<std::pair<size_t, double>> miss_nodes_vec;
    // check if an edge is removed
    std::vector<bool> edgeMark_cost(__numRRsets, false);
    std::vector<bool> node_sel_cost(__numV, false);
    min_cost = min_vector_nonzero_mask(budget_list, node_sel_cost);

    Timer timerBound("Bound");
    __vecSeed_inf_cost.clear();
    for (size_t deg = maxDeg_cost; deg > 0; deg--) // Enusre deg > 0
    {
        if (sumBudget_cost + min_cost > budget) {
            break;
        }
        std::sort(degMap_cost[deg].begin(), degMap_cost[deg].end(), [](auto & left, auto & right) {
            return std::get<1>(left) < std::get<1>(right);
        });
        auto & vecNode = degMap_cost[deg];
        for (auto idx = vecNode.size(); idx--;) {
            auto argmaxIdx = std::get<0>(vecNode[idx]);
            const size_t currDeg = int(coverage_cost[argmaxIdx]) + 1;
            if (deg > currDeg) {
                degMap_cost[currDeg].push_back(std::make_pair(argmaxIdx, coverage_cost[argmaxIdx]));
                degMap_cost[deg].pop_back();
                continue;
            }
            if (sumBudget_cost + budget_list[argmaxIdx] <= budget) {
                if (mode == 2) {
                    timerBound.get_operation_time();
                    // Find upper bound
                    auto degBound = deg;
                    Nodelist vecBound;
                    double vecBound_Budget = 0;
                    // Initialize vecBound
                    auto idxBound = idx + 1;
                    int miss_idxBound = miss_nodes_vec.size() - 1;
                    bool is_stop = false;
                    for (size_t miss_idx = miss_nodes_vec.size(); miss_idx--; ) {
                        std::get<1>(miss_nodes_vec[miss_idx]) = coverage_cost[std::get<0>(miss_nodes_vec[miss_idx])];
                    }
                    std::sort(miss_nodes_vec.begin(), miss_nodes_vec.end(), [](auto & left, auto & right) {
                        return std::get<1>(left) < std::get<1>(right);
                    });
                    if (is_stop == false) {
                        while (vecBound_Budget + min_cost <= budget && idxBound--) {
                            if (miss_idxBound >= 0) {
                                if (std::get<1>(degMap_cost[degBound][idxBound]) >= std::get<1>(miss_nodes_vec[miss_idxBound])) {
                                    if (vecBound_Budget + budget_list[std::get<0>(degMap_cost[degBound][idxBound])] <= budget) {
                                        vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                                           budget_list[std::get<0>(degMap_cost[degBound][idxBound])]);
                                        vecBound_Budget += budget_list[std::get<0>(degMap_cost[degBound][idxBound])];
                                    }
                                    else {
                                        vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                                           (budget - vecBound_Budget));
                                        vecBound_Budget = budget;
                                        is_stop = true;
                                        break;
                                    }
                                }
                                else {
                                    if (vecBound_Budget + budget_list[std::get<0>(miss_nodes_vec[miss_idxBound])] <= budget) {
                                        vecBound.push_back(std::get<1>(miss_nodes_vec[miss_idxBound]) *
                                                           budget_list[std::get<0>(miss_nodes_vec[miss_idxBound])]);
                                        vecBound_Budget += budget_list[std::get<0>(miss_nodes_vec[miss_idxBound])];
                                        miss_idxBound = miss_idxBound - 1;
                                    }
                                    else {
                                        vecBound.push_back(std::get<1>(miss_nodes_vec[miss_idxBound]) *
                                                           (budget - vecBound_Budget));
                                        vecBound_Budget = budget;
                                        is_stop = true;
                                        break;
                                    }
                                }
                            }
                            else {
                                if (vecBound_Budget + budget_list[std::get<0>(degMap_cost[degBound][idxBound])] <= budget) {
                                    vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                                       budget_list[std::get<0>(degMap_cost[degBound][idxBound])]);
                                    vecBound_Budget += budget_list[std::get<0>(degMap_cost[degBound][idxBound])];
                                }
                                else {
                                    vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                                       (budget - vecBound_Budget));
                                    vecBound_Budget = budget;
                                    is_stop = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (is_stop == false) {
                        while (vecBound_Budget + min_cost <= budget && --degBound) {
                            std::sort(degMap_cost[degBound].begin(), degMap_cost[degBound].end(), [](auto & left, auto & right) {
                                return std::get<1>(left) < std::get<1>(right);
                            });
                            idxBound = degMap_cost[degBound].size();
                            while (vecBound_Budget + min_cost <= budget && idxBound--) {
                                if (miss_idxBound >= 0) {
                                    if (std::get<1>(degMap_cost[degBound][idxBound]) >= std::get<1>(miss_nodes_vec[miss_idxBound])) {
                                        if (vecBound_Budget + budget_list[std::get<0>(degMap_cost[degBound][idxBound])] <= budget) {
                                            vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                                               budget_list[std::get<0>(degMap_cost[degBound][idxBound])]);
                                            vecBound_Budget += budget_list[std::get<0>(degMap_cost[degBound][idxBound])];
                                        }
                                        else {
                                            vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                                               (budget - vecBound_Budget));
                                            vecBound_Budget = budget;
                                            is_stop = true;
                                            break;
                                        }
                                    }
                                    else {
                                        if (vecBound_Budget + budget_list[std::get<0>(miss_nodes_vec[miss_idxBound])] <= budget) {
                                            vecBound.push_back(std::get<1>(miss_nodes_vec[miss_idxBound]) *
                                                               budget_list[std::get<0>(miss_nodes_vec[miss_idxBound])]);
                                            vecBound_Budget += budget_list[std::get<0>(miss_nodes_vec[miss_idxBound])];
                                            miss_idxBound = miss_idxBound - 1;
                                        }
                                        else {
                                            vecBound.push_back(std::get<1>(miss_nodes_vec[miss_idxBound]) *
                                                               (budget - vecBound_Budget));
                                            vecBound_Budget = budget;
                                            is_stop = true;
                                            break;
                                        }
                                    }
                                }
                                else {
                                    if (vecBound_Budget + budget_list[std::get<0>(degMap_cost[degBound][idxBound])] <= budget) {
                                        vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                                           budget_list[std::get<0>(degMap_cost[degBound][idxBound])]);
                                        vecBound_Budget += budget_list[std::get<0>(degMap_cost[degBound][idxBound])];
                                    }
                                    else {
                                        vecBound.push_back(std::get<1>(degMap_cost[degBound][idxBound]) *
                                                           (budget - vecBound_Budget));
                                        vecBound_Budget = budget;
                                        is_stop = true;
                                        break;
                                    }
                                }
                            }
                            if (is_stop == true) {
                                break;
                            }
                        }
                    }
                    std::sort(vecBound.begin(), vecBound.end());
                    make_min_heap(vecBound);

                    __boundLast_inf_cost = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0))
                                                  + sumInf_cost) * __numV /__numRRsets;
                    if (__boundMin_inf_cost > __boundLast_inf_cost) {
                        __boundMin_inf_cost = __boundLast_inf_cost;
                    }
                    timeb += timerBound.get_operation_time();
                }
                sumInf_cost = sumInf_cost + coverage_cost[argmaxIdx] * budget_list[argmaxIdx];
                __vecSeed_inf_cost.push_back(argmaxIdx);
                sumBudget_cost = sumBudget_cost + budget_list[argmaxIdx];
                coverage_cost[argmaxIdx] = 0;
                node_sel_cost[argmaxIdx] = true;
                min_cost = min_vector_nonzero_mask(budget_list, node_sel_cost);
                for (auto edgeIdx : __hyperG._FRsets[argmaxIdx]) {
                    if (edgeMark_cost[edgeIdx]) {
                        continue;
                    }
                    edgeMark_cost[edgeIdx] = true;
                    for (auto nodeIdx : __hyperG._RRsets[edgeIdx]) {
                        if (coverage_cost[nodeIdx] == 0) {
                            continue; // This node is seed, skip
                        }
                        coverage_cost[nodeIdx] = (double) (coverage_cost[nodeIdx] * budget_list[nodeIdx] - 1) /
                                                 budget_list[nodeIdx];
                    }
                }
                degMap_cost[deg].pop_back();
                for (auto degup = deg; degup > 0; degup--) {
                    auto & vecNodeup = degMap_cost[degup];
                    for(auto idxup = vecNodeup.size(); idxup--;) {
                        std::get<1>(degMap_cost[degup][idxup]) = coverage_cost[std::get<0>(degMap_cost[degup][idxup])];
                    }
                }
                std::sort(degMap_cost[deg].begin(), degMap_cost[deg].end(), [](auto & left, auto & right) {
                    return std::get<1>(left) < std::get<1>(right);
                });
            }
            else{
                miss_nodes_vec.push_back(std::make_pair(argmaxIdx, coverage_cost[argmaxIdx]));
                degMap_cost[deg].pop_back();
            }
            if (sumBudget_cost + min_cost > budget) {
                break;
            }
        }
        degMap_cost.pop_back();
    }
    const auto finalInf_cost = 1.0 * sumInf_cost * __numV / __numRRsets;
    std::cout << "  >>>[greedy-lazy] influence_cost: " << sumInf_cost << " " <<
              finalInf_cost << ", min-bound_cost: " <<
              __boundMin_inf_cost << ", last-bound_cost: " << __boundLast_inf_cost << '\n';

    if (finalInf > finalInf_cost) {
        __vecSeed.assign(__vecSeed_inf.begin(), __vecSeed_inf.end());
        std::cout << "Inf: " <<  __vecSeed.size() << std::endl;
        __is_inf_cost = false;
        return finalInf;
    }
    else {
        __vecSeed.assign(__vecSeed_inf_cost.begin(), __vecSeed_inf_cost.end());
        std::cout << "Cost Inf: " << __vecSeed.size() << std::endl;
        __is_inf_cost = true;
        return finalInf_cost;
    }
}

double Alg::max_cover_lazy_budget_fast_end(const Graph & graph, const double budget,
        const std::vector<double> & budget_list, const double epsilon, double & timeb, const int mode) {
    // mode: optimization mode.
    // 0->no optimization,
    // 1->optimization with the upper bound in the last round,
    // 2->optimization with minimum upper bound among all rounds [Default].
    const double EPS = 0.0000000001;
    __boundLast = DBL_MAX, __boundMin = DBL_MAX;
    __boundLast_inf = DBL_MAX, __boundMin_inf = DBL_MAX;
    __boundLast_inf_cost = DBL_MAX, __boundMin_inf_cost = DBL_MAX;
    FRset coverage(__numV, 0);
    size_t maxDeg = 0;
    for (auto i = __numV; i--;) {
        const size_t deg = __hyperG._FRsets[i].size();
        coverage[i] = deg;
        if (deg > maxDeg) {
            maxDeg = deg;
        }
    }
    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
    for (auto i = __numV; i--;) {
        if (coverage[i] == 0) {
            continue;
        }
        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    int flag = 0;
    std::vector<bool> edgeMark(__numRRsets, false);
    // check if an edge is removed
    std::vector<bool> node_sel(__numV, false);
    double min_cost = min_vector_nonzero_mask(budget_list, node_sel);
    std::cout << "min_cost: " << min_cost << std::endl;

    __vecSeed_inf.clear();
    for (size_t deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        if (min_cost > budget) {
            break;
        }
        auto & vecNode = degMap[deg];
        for (auto idx = vecNode.size(); idx--;) {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];
            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }
            if (budget_list[argmaxIdx] <= budget) {
                if (mode == 2) {
                    // Find upper bound
                    auto degBound = deg;
                    Nodelist vecBound;
                    // Initialize vecBound
                    auto idxBound = idx;
                    vecBound.push_back(coverage[degMap[degBound][idxBound]]);
                    __boundLast_inf = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0))
                                             + sumInf) * __numV /__numRRsets;
                    if (__boundMin_inf > __boundLast_inf) {
                        __boundMin_inf = __boundLast_inf;
                    }
                }
                if (flag == 0) {
                    sumInf = sumInf + coverage[argmaxIdx];
                    __vecSeed_inf.push_back(argmaxIdx);
                    coverage[argmaxIdx] = 0;
                    for (auto edgeIdx : __hyperG._FRsets[argmaxIdx])
                    {
                        if (edgeMark[edgeIdx]) continue;
                        edgeMark[edgeIdx] = true;
                        for (auto nodeIdx : __hyperG._RRsets[edgeIdx])
                        {
                            if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
                            coverage[nodeIdx]--;
                        }
                    }
                }
                flag = flag + 1;
            }
            if (flag == 2) {
                break;
            }
        }
        if (flag == 2) {
            break;
        }
        degMap.pop_back();
    }
    // Top-k influential nodes constructed
    const auto finalInf = 1.0 * sumInf * __numV / __numRRsets;
    std::cout << std::endl;
    std::cout << "  >>>[greedy-lazy] influence: " << sumInf << " " << finalInf <<
              ", min-bound: " << __boundMin_inf <<
              ", last-bound: " << __boundLast_inf << '\n';

    // cost influence most
    std::vector<double> coverage_cost(__numV, 0);
    double maxDeg_cost = 0;
    double minDeg_cost = DBL_MAX;
    for (auto i = __numV; i--;) {
        const size_t deg_cost = __hyperG._FRsets[i].size();
        if (fabs(budget_list[i]) <= EPS || deg_cost == 0) {
            coverage_cost[i] = 0;
        }
        else {
            coverage_cost[i] = (double) deg_cost / budget_list[i];
        }
        if (coverage_cost[i] > maxDeg_cost) {
            maxDeg_cost = coverage_cost[i];
        }
        if (coverage_cost[i] < minDeg_cost) {
            if (fabs(coverage_cost[i]) > EPS) {
                minDeg_cost = coverage_cost[i];
            }
        }
    }
    std::cout << "minDeg_cost: " <<  minDeg_cost << std::endl;

    size_t maxIdx_cost = 0;
    for (double deg_thr = maxDeg_cost; deg_thr >= minDeg_cost * (1 - epsilon); deg_thr = deg_thr * (1 - epsilon)) {
        maxIdx_cost++;
    }

    // initial reverse index
    RRsets degMap_cost(maxIdx_cost); // degMap: map degree threshold to the nodes over this threshold
    std::cout << "maxIdx_cost: " << maxIdx_cost << std::endl;
    size_t tmpidx_cost = maxIdx_cost - 1;
    for (double deg_thr = maxDeg_cost; deg_thr >= minDeg_cost * (1 - epsilon); deg_thr = deg_thr * (1 - epsilon)) {
        for (auto i = __numV; i--; ) {
            if (coverage_cost[i] >= deg_thr && coverage_cost[i] <= deg_thr / (1 - epsilon)) {
                degMap_cost[tmpidx_cost].push_back(i);
            }
        }
        if (tmpidx_cost == 0) {
            break;
        }
        tmpidx_cost--;
    }

    size_t sumInf_cost = 0;
    double sumBudget_cost = 0;
    // miss nodes
    std::vector<size_t> miss_nodes_vec;
    RRsets missMap_cost(maxIdx_cost);
    // check if an edge is removed
    std::vector<bool> edgeMark_cost(__numRRsets, false);
    std::vector<bool> node_sel_cost(__numV, false);
    min_cost = min_vector_nonzero_mask(budget_list, node_sel_cost);

    Timer timerBound("Bound");
    __vecSeed_inf_cost.clear();
    for (size_t degidx = degMap_cost.size(); degidx--; ) {
        if (sumBudget_cost + min_cost > budget) {
            break;
        }
        auto & vecNode = degMap_cost[degidx];
        for (size_t idx = vecNode.size(); idx--;) {
            auto argmaxIdx = vecNode[idx];
            const double currDeg = coverage_cost[argmaxIdx];
            if (sumBudget_cost + budget_list[argmaxIdx] <= budget) {
                if (mode == 2 || (mode == 1 && sumBudget_cost + budget_list[argmaxIdx] + min_cost > budget)) {
                    timerBound.get_operation_time();
                    // Find upper bound
                    auto degBound = degidx;
                    Nodelist vecBound;
                    double vecBound_Budget = 0;
                    // Initialize vecBound
                    auto idxBound = idx + 1;
                    size_t lastNum = 0;
                    double vecBound_Budget_Last = 0;
                    bool is_stop = false;
                    size_t tmp_miss_idx_cost = 0;
                    for (auto miss_idx = miss_nodes_vec.size(); miss_idx--; ) {
                        tmp_miss_idx_cost = maxIdx_cost - 1;
                        for (double miss_deg_thr = maxDeg_cost;
                                miss_deg_thr >= minDeg_cost * (1 - epsilon); miss_deg_thr = miss_deg_thr * (1 - epsilon)) {
                            if (coverage_cost[miss_nodes_vec[miss_idx]] >= miss_deg_thr &&
                                    coverage_cost[miss_nodes_vec[miss_idx]] <= miss_deg_thr / (1 - epsilon)) {
                                missMap_cost[tmp_miss_idx_cost].push_back(miss_nodes_vec[miss_idx]);
                                break;
                            }
                            if (tmp_miss_idx_cost == 0) {
                                break;
                            }
                            tmp_miss_idx_cost--;
                        }
                    }
                    auto miss_idxBound = missMap_cost[degBound].size();
                    if (is_stop == false) {
                        while (vecBound_Budget + min_cost <= budget && idxBound--) {
                            if (vecBound_Budget + budget_list[degMap_cost[degBound][idxBound]] <= budget) {
                                vecBound.push_back(coverage_cost[degMap_cost[degBound][idxBound]] *
                                                   budget_list[degMap_cost[degBound][idxBound]]);
                                vecBound_Budget = vecBound_Budget + budget_list[degMap_cost[degBound][idxBound]];
                                lastNum = lastNum + 1;
                            }
                            else {
                                vecBound.push_back(coverage_cost[degMap_cost[degBound][idxBound]] *
                                                   (budget - vecBound_Budget));
                                vecBound_Budget = budget;
                                is_stop = true;
                                break;
                            }
                        }
                        while (vecBound_Budget + min_cost <= budget && miss_idxBound--) {
                            if (vecBound_Budget + budget_list[missMap_cost[degBound][miss_idxBound]] <= budget) {
                                vecBound.push_back(coverage_cost[missMap_cost[degBound][miss_idxBound]] *
                                                   budget_list[missMap_cost[degBound][miss_idxBound]]);
                                vecBound_Budget = vecBound_Budget + budget_list[missMap_cost[degBound][miss_idxBound]];
                                lastNum = lastNum + 1;
                            }
                            else {
                                vecBound.push_back(coverage_cost[missMap_cost[degBound][miss_idxBound]] *
                                                   (budget - vecBound_Budget));
                                vecBound_Budget = budget;
                                is_stop = true;
                                break;
                            }
                        }
                    }
                    if (is_stop == false) {
                        while (vecBound_Budget + min_cost <= budget && --degBound) {
                            lastNum = 0;
                            vecBound_Budget_Last = vecBound_Budget;
                            idxBound = degMap_cost[degBound].size();
                            miss_idxBound = missMap_cost[degBound].size();
                            while (vecBound_Budget + min_cost <= budget && idxBound--) {
                                if (vecBound_Budget + budget_list[degMap_cost[degBound][idxBound]] <= budget) {
                                    vecBound.push_back(coverage_cost[degMap_cost[degBound][idxBound]] *
                                                       budget_list[degMap_cost[degBound][idxBound]]);
                                    vecBound_Budget = vecBound_Budget + budget_list[degMap_cost[degBound][idxBound]];
                                    lastNum = lastNum + 1;
                                }
                                else {
                                    vecBound.push_back(coverage_cost[degMap_cost[degBound][idxBound]] *
                                                       (budget - vecBound_Budget));
                                    vecBound_Budget = budget;
                                    is_stop = true;
                                    break;
                                }
                            }
                            while (vecBound_Budget + min_cost <= budget && miss_idxBound--) {
                                if (vecBound_Budget + budget_list[missMap_cost[degBound][miss_idxBound]] <= budget) {
                                    vecBound.push_back(coverage_cost[missMap_cost[degBound][miss_idxBound]] *
                                                       budget_list[missMap_cost[degBound][miss_idxBound]]);
                                    vecBound_Budget = vecBound_Budget + budget_list[missMap_cost[degBound][miss_idxBound]];
                                    lastNum = lastNum + 1;
                                }
                                else {
                                    vecBound.push_back(coverage_cost[missMap_cost[degBound][miss_idxBound]] *
                                                       (budget - vecBound_Budget));
                                    vecBound_Budget = budget;
                                    is_stop = true;
                                    break;
                                }
                            }
                            if (is_stop == true) {
                                break;
                            }
                        }
                    }
                    // roll back
                    if (lastNum != 0 && idxBound != 0) {
                        while (lastNum) {
                            vecBound.pop_back();
                            lastNum = lastNum - 1;
                        }
                        vecBound_Budget = vecBound_Budget_Last;
                        std::vector<std::pair<size_t, double>> vecNodeLastPair(degMap_cost[degBound].size()
                            + missMap_cost[degBound].size());
                        for (size_t idxpair = degMap_cost[degBound].size(); idxpair--; ) {
                            vecNodeLastPair[idxpair] = std::make_pair(degMap_cost[degBound][idxpair],
                                                                      coverage_cost[degMap_cost[degBound][idxpair]]);
                        }
                        for (size_t miss_idxpair = missMap_cost[degBound].size(); miss_idxpair--; ){
                            vecNodeLastPair[miss_idxpair + degMap_cost[degBound].size()] =
                                    std::make_pair(missMap_cost[degBound][miss_idxpair],
                                            coverage_cost[missMap_cost[degBound][miss_idxpair]]);
                        }
                        std::sort(vecNodeLastPair.begin(), vecNodeLastPair.end(), [](auto & left, auto & right){
                            return std::get<1>(left) < std::get<1>(right);
                        });
                        // choose max
                        idxBound = degMap_cost[degBound].size() + missMap_cost[degBound].size();
                        while (vecBound_Budget + min_cost <= budget && idxBound--) {
                            if (vecBound_Budget + budget_list[std::get<0>(vecNodeLastPair[idxBound])] <= budget) {
                                vecBound.push_back(coverage_cost[std::get<0>(vecNodeLastPair[idxBound])] *
                                                   budget_list[std::get<0>(vecNodeLastPair[idxBound])]);
                                vecBound_Budget = vecBound_Budget + budget_list[std::get<0>(vecNodeLastPair[idxBound])];
                            }
                            else {
                                vecBound.push_back(coverage_cost[std::get<0>(vecNodeLastPair[idxBound])] *
                                                   (budget- vecBound_Budget));
                                vecBound_Budget = budget;
                            }
                        }
                    }
                    else if (idxBound == 0 && vecBound_Budget < budget) {
                        double max_coverage_ratio = 0.0;
                        degBound = degBound - 1;
                        for (size_t idxpair = degMap_cost[degBound].size(); idxpair--; ) {
                            if (coverage_cost[degMap_cost[degBound][idxpair]] > max_coverage_ratio) {
                                max_coverage_ratio = coverage_cost[degMap_cost[degBound][idxpair]];
                            }
                        }
                        for (size_t miss_idxpair = missMap_cost[degBound].size(); miss_idxpair--; ) {
                            if (coverage_cost[missMap_cost[degBound][miss_idxpair]] > max_coverage_ratio) {
                                max_coverage_ratio = coverage_cost[missMap_cost[degBound][miss_idxpair]];
                            }
                        }
                        vecBound.push_back(max_coverage_ratio * (budget- vecBound_Budget));
                    }
                    std::sort(vecBound.begin(), vecBound.end());
                    make_min_heap(vecBound);
                    __boundLast_inf_cost = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0))
                                                  + sumInf_cost) * __numV /__numRRsets;
                    if (__boundMin_inf_cost > __boundLast_inf_cost) {
                        __boundMin_inf_cost = __boundLast_inf_cost;
                    }
                    timeb += timerBound.get_operation_time();
                }
                sumInf_cost = sumInf_cost + currDeg * budget_list[argmaxIdx];
                __vecSeed_inf_cost.push_back(argmaxIdx);
                sumBudget_cost = sumBudget_cost + budget_list[argmaxIdx];
                coverage_cost[argmaxIdx] = 0;
                node_sel_cost[argmaxIdx] = true;
                min_cost = min_vector_nonzero_mask(budget_list, node_sel_cost);
                for (auto edgeIdx : __hyperG._FRsets[argmaxIdx]) {
                    if (edgeMark_cost[edgeIdx]) {
                        continue;
                    }
                    edgeMark_cost[edgeIdx] = true;
                    for (auto nodeIdx : __hyperG._RRsets[edgeIdx]) {
                        if (fabs(coverage_cost[nodeIdx]) <= EPS) {
                            continue; // This node is seed, skip
                        }
                        if (fabs(budget_list[nodeIdx]) >= EPS) {
                            coverage_cost[nodeIdx] = (coverage_cost[nodeIdx] * budget_list[nodeIdx] - 1) /
                                                     budget_list[nodeIdx];
                        } else {
                            coverage_cost[nodeIdx] = 0;
                        }
                    }
                }
                degMap_cost[degidx].pop_back();
            }
            else {
                if (fabs(coverage_cost[argmaxIdx]) > EPS) {
                    miss_nodes_vec.push_back(argmaxIdx);
//                    std::cout << "miss node size: " << miss_nodes_vec.size() << std::endl;
                }
                degMap_cost[degidx].pop_back();
            }
            if (sumBudget_cost + min_cost > budget) {
                break;
            }
        }
        degMap_cost.pop_back();
    }
    const auto finalInf_cost = 1.0 * sumInf_cost * __numV / __numRRsets;
    std::cout << "  >>>[greedy-lazy] influence_cost: " << sumInf_cost << " " <<
              finalInf_cost << ", min-bound_cost: " <<
              __boundMin_inf_cost << ", last-bound_cost: " << __boundLast_inf_cost << '\n';

    if (finalInf > finalInf_cost) {
        __vecSeed.clear();
        __vecSeed.assign(__vecSeed_inf.begin(), __vecSeed_inf.end());
        std::cout << "Inf: " <<  __vecSeed.size() << std::endl;
        __is_inf_cost = false;
        return finalInf;
    }
    else {
        __vecSeed.clear();
        __vecSeed.assign(__vecSeed_inf_cost.begin(), __vecSeed_inf_cost.end());
        std::cout << "Cost Inf: " << __vecSeed.size() << std::endl;
        __is_inf_cost = true;
        return finalInf_cost;
    }
}

double Alg::max_cover(const int targetSize, const int mode)
{
    if (targetSize >= 1000) return max_cover_topk(targetSize);
    return max_cover_lazy(targetSize, mode);
}

double Alg::max_cover_budget(const Graph & graph, const double budget,
        const std::vector<double> & budget_list, double & timeb, const int mode) {
    return max_cover_lazy_budget(graph, budget, budget_list, timeb, mode);
}

double Alg::max_cover_budget_fast_end(const Graph & graph, const double budget,
        const std::vector<double> & budget_list, const double epsilon, double & timeb, const int mode) {
    return max_cover_lazy_budget_fast_end(graph, budget, budget_list, epsilon, timeb, mode);
}

double Alg::max_vector(const std::vector<double> & budget_list) {
    double tmp = WINT_MIN;
    for (double i : budget_list) {
        if (i > tmp) {
            tmp = i;
        }
    }
    return tmp;
}

double Alg::max_vector_nonzero(const std::vector<double> & budget_list) {
    double tmp = WINT_MIN;
    const double esp = 0.0001;
    for (double i : budget_list) {
        if (i > tmp && fabs(i) > esp) {
            tmp = i;
        }
    }
    return tmp;
}

double Alg::min_bound(const double a, const double b) {
    if (a > b) {
        return b;
    }
    else {
        return a;
    }
}

double Alg::min_vector(const std::vector<double> & budget_list) {
    double tmp = WINT_MAX;
    for (double i : budget_list) {
        if (i < tmp) {
            tmp = i;
        }
    }
    return tmp;
}

double Alg::min_vector_nonzero(const std::vector<double> & budget_list) {
    double tmp = WINT_MAX;
    const double esp = 0.0001;
    for (double i : budget_list) {
        if (i < tmp && fabs(i) > esp) {
            tmp = i;
        }
    }
    return tmp;
}

double Alg::min_vector_nonzero_mask(const std::vector<double> & budget_list, const std::vector<bool> & mask) {
    double tmp = WINT_MAX;
    const double esp = 0.0001;
    for (size_t i = 0; i < budget_list.size(); i++) {
        if (budget_list[i] < tmp && fabs(budget_list[i]) > esp && mask[i] == false) {
            tmp = budget_list[i];
        }
    }
    return tmp;
}

/// return k_max
int Alg::cal_k_max(std::vector<double> budget_list, double budget) {
    int count = 0;
    std::sort(budget_list.begin(), budget_list.end());
    for(int i = 0; i < budget_list.size(); i++) {
        if (budget > budget_list[i] && fabs(budget_list[i]) > 0.001) {
            budget = budget - budget_list[i];
            count = count + 1;
        }
    }
    return count;
}

/// return k_min
int Alg::cal_k_min(std::vector<double> budget_list, double budget) {
    int count = 0;
    std::sort(budget_list.begin(), budget_list.end(), std::greater<double>());
    for(int i = 0; i < budget_list.size(); i++) {
        if (budget > budget_list[i] && fabs(budget_list[i]) > 0.001) {
            budget = budget - budget_list[i];
            count = count + 1;
        }
    }
    return count;
}

double Alg::max_element(double a, double b) {
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

double Alg::productLog(double realz, double imagez) {
    if (fabs(realz) + fabs(imagez) < 1e-10) {
        return realz;
    }
    double lnx = log(realz);
    double zx = lnx;
    double zy = 0;
    double F[2];
    double Fx[2];
    double Fy[2];
    double temp = 0;
    double x = zx;
    double y = zy;
    double x0 = 0;
    double y0 = 0;
    Fy[0] = log(x * x + y * y) / 2 + x - zx;
    Fy[1] = y - zy + atan2(y, x);
    double error = fabs(Fy[0]) + fabs(Fy[1]);
    int loopn = 1000;
    double w = 1;
    while (loopn > 0 && w > -1.05 && error > 1e-10) {
        w = x * x + y * y;
        F[0] = x / w + 1;
        F[1] = y / w;
        w = F[0] * F[0] + F[1] * F[1];
        Fx[0] = (Fy[0] * F[0] - Fy[1] * F[1]) / w;
        Fx[1] = (Fy[1] * F[0] - Fy[0] * F[1]) / w;
        w = 1;
        while (w >= -1){
            x0 = x - w * Fx[0];
            y0 = y - w * Fx[1];
            Fy[0] = log(x0 * x0 + y0 * y0) / 2 + x0 - zx;
            Fy[1] = y0 - zy + atan2(y0, x0);
            temp = fabs(Fy[0]) + fabs(Fy[1]);
            if (temp < error) {
                error = temp;
                x = x0;
                y = y0;
                w = -1.03;
            }
            else {
                w = w - 0.1;
            }
        }
    }
    return x;
}

void Alg::set_cascade_model(const CascadeModel model)
{
    __hyperG.set_cascade_model(model);
    __hyperGVldt.set_cascade_model(model);
}

double Alg::effic_inf_valid_algo()
{
    return effic_inf_valid_algo(__vecSeed);
}

double Alg::effic_inf_valid_algo(const Nodelist vecSeed)
{
    Timer EvalTimer("Inf. Eval.");
    std::cout << "  >>>Evaluating influence in [0.99,1.01]*EPT with prob. 99.9% ...\n";
    const auto inf = __hyperG.effic_inf_valid_algo(vecSeed);
    //const auto inf = __hyperGVldt.effic_inf_valid_algo(vecSeed);
    std::cout << "  >>>Down! influence: " << inf << ", time used (sec): " << EvalTimer.get_total_time() << '\n';
    return inf;
}

double Alg::opim(const int targetSize, const size_t numRRsets, const double delta, const int mode)
{
    Timer timerOPIM("OPIM");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double a1 = log(2.0 / delta);
    const double a2 = log(2.0 / delta);
    __hyperG.build_n_RRsets(numRRsets / 2); // R1
    __hyperGVldt.build_n_RRsets(numRRsets / 2); // R2
    __numRRsets = __hyperG.get_RR_sets_size();
    const auto time1 = timerOPIM.get_operation_time();
    const auto infSelf = max_cover(targetSize, mode);
    const auto time2 = timerOPIM.get_operation_time();
    const auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
    const auto degVldt = infVldt * __numRRsets / __numV;
    auto upperBound = infSelf / approx;
    if (mode == 1) upperBound = __boundLast;
    else if (mode == 2) upperBound = __boundMin;
    const auto upperDegOPT = upperBound * __numRRsets / __numV;
    const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
    const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
    const auto approxOPIM = lowerSelect / upperOPT;
    //guaranteeOpt = opt_error(degVldt, degSelfUpper, delta);
    __tRes.set_approximation(approxOPIM);
    __tRes.set_running_time(timerOPIM.get_total_time());
    __tRes.set_influence(infVldt);
    __tRes.set_influence_original(infSelf);
    __tRes.set_seed_vec(__vecSeed);
    __tRes.set_RR_sets_size(__numRRsets * 2);
    std::cout << "==>OPIM approx. (max-cover): " << approxOPIM << " (" << infSelf / upperBound << ")\n";
    std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
    __tRes.set_influence(effic_inf_valid_algo());
    return approxOPIM;
}

double Alg::opimc(const int targetSize, const double epsilon, const double delta, const int mode)
{
    Timer timerOPIMC("OPIM-C");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(__numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto maxNumR = size_t(2.0 * __numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0;
    for (auto idx = 0; idx < numIter; idx++)
    {
        const auto numR = numRbase << idx;
        timerOPIMC.get_operation_time();
        //dsfmt_gv_init_gen_rand(idx);
        __hyperG.build_n_RRsets(numR); // R1
        //dsfmt_gv_init_gen_rand(idx + 1000000);
        __hyperGVldt.build_n_RRsets(numR); // R2
        __numRRsets = __hyperG.get_RR_sets_size();
        time1 += timerOPIMC.get_operation_time();
        const auto infSelf = max_cover(targetSize, mode);
        time2 += timerOPIMC.get_operation_time();
        const auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        const auto degVldt = infVldt * __numRRsets / __numV;
        auto upperBound = infSelf / approx;
        if (mode == 1) upperBound = __boundLast;
        else if (mode == 2) upperBound = __boundMin;
        const auto upperDegOPT = upperBound * __numRRsets / __numV;
        const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto approxOPIMC = lowerSelect / upperOPT;
        //guaranteeOpt = opt_error(degVldt, degSelfUpper, delta);
        std::cout << " -->OPIM-C (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << approxOPIMC <<
                  " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        // Check whether the requirement is satisfied
        if (approxOPIMC >= approx - epsilon)
        {
            __tRes.set_approximation(approxOPIMC);
            __tRes.set_running_time(timerOPIMC.get_total_time());
            __tRes.set_influence(infVldt);
            __tRes.set_influence_original(infSelf);
            __tRes.set_seed_vec(__vecSeed);
            __tRes.set_seed_inf_vec(__vecSeed_inf);
            __tRes.set_seed_inf_cost_vec(__vecSeed_inf_cost);
            __tRes.set_boundmin_inf(__boundMin_inf);
            __tRes.set_boundmin_inf_cost(__boundMin_inf_cost);
            __tRes.set_RR_sets_size(__numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
            __tRes.set_influence(effic_inf_valid_algo());
            return __tRes.get_influence();
        }
    }
    return 0.0;
}

double Alg::opimb(const Graph & graph, const double budget, const std::vector<double> & budget_list,
        const double epsilon, const double delta, const int mode) {
    Timer timerOPIMC("OPIM-B");
    const double e_alpha = exp(0.4384521543);
    const double approx = 1 - 1.0 / e_alpha;
    const double alpha = sqrt(log(6.0 / delta));
    const int k_max = cal_k_max(budget_list, budget);
    const int k_min = cal_k_min(budget_list, budget);
    const auto theta_zero = (size_t) ((2.0 * __numV * pow2((1 - 1 / e_alpha) * alpha
            + sqrt((1 - 1 / e_alpha) * (k_max * log(__numV) + log(6.0 / delta))))) / (__numV));
    const auto theta_max = (size_t) ((2.0 * __numV * pow2((1 - 1 / e_alpha) * alpha
            + sqrt((1 - 1 / e_alpha) * (k_max * log(__numV) + log(6.0 / delta))))) / (pow2(epsilon) * k_max));
    const auto numIter = (size_t) log2(theta_max / theta_zero) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0;
    double timeb = 0.0;
    int satisfy_count = 0;
    size_t idx = 0;
    while (true) {
        std::cout << idx << std::endl;
        const auto numR = theta_zero << idx;
        timerOPIMC.get_operation_time();
        __hyperG.build_n_RRsets(numR); // R1
        __hyperGVldt.build_n_RRsets(numR); // R2
        std::cout << "# RR sets: " << numR << std::endl;
        __numRRsets = __hyperG.get_RR_sets_size();
        time1 += timerOPIMC.get_operation_time();
        const auto infSelf = max_cover_budget(graph, budget, budget_list, timeb, mode);
        time2 += timerOPIMC.get_operation_time();
        const auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        const auto degVldt = infVldt * __numRRsets / __numV;
        auto upperBound = infSelf / approx;
        auto upperBound_med = infSelf / approx;
        if (mode == 1) {
            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
                upperBound = min_bound(upperBound, upperBound_med);
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        std::cout << "========upper bound ratio: " << upperBound / upperBound_med << "========" << std::endl;
        const auto upperDegOPT = upperBound * __numRRsets / __numV;
        auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-B (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << approxOPIMC <<
                  " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        // Check whether the requirement is satisfied
        if (approxOPIMC >= approx - epsilon && satisfy_count == 1)
        {
            __tRes.set_approximation(approxOPIMC);
            __tRes.set_running_time(timerOPIMC.get_total_time());
            __tRes.set_influence(infVldt);
            __tRes.set_influence_original(infSelf);
            __tRes.set_seed_vec(__vecSeed);
            __tRes.set_seed_inf_vec(__vecSeed_inf);
            __tRes.set_seed_inf_cost_vec(__vecSeed_inf_cost);
            __tRes.set_boundmin_inf(__boundMin_inf);
            __tRes.set_boundmin_inf_cost(__boundMin_inf_cost);
            __tRes.set_RR_sets_size(__numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 - timeb << '\n';
            std::cout << "==>Time for Bound Computation: " << timeb << '\n';
            return 0;
        }
        else if (approxOPIMC >= approx - epsilon && satisfy_count == 0) {
            satisfy_count = satisfy_count + 1;
        }
        idx = idx + 1;
    }
}

double Alg::opimb_approx(const Graph & graph, const double budget, const std::vector<double> & budget_list,
        const double epsilon, const double delta, const int mode) {
    Timer timerOPIMC("OPIM-B-APPROX");
    const double approx = epsilon;
    const double alpha = sqrt(log(6.0 / delta));
    const int k_max = cal_k_max(budget_list, budget);
    const int k_min = cal_k_min(budget_list, budget);
    const auto theta_zero = (size_t) ((2.0 * __numV * pow2(approx * alpha
            + sqrt(approx * (k_max * log(__numV) + log(6.0 / delta))))) / (__numV));
    const auto theta_max = (size_t) ((2.0 * __numV * pow2(approx * alpha
            + sqrt(approx * (k_max * log(__numV) + log(6.0 / delta))))) / (pow2(epsilon) * k_max));
    const auto numIter = (size_t) log2(theta_max / theta_zero) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0;
    double timeb = 0.0;
    int satisfy_count = 0;
    size_t idx = 0;
    while (true) {
        std::cout << idx << std::endl;
        const auto numR = theta_zero << idx;
        timerOPIMC.get_operation_time();
        __hyperG.build_n_RRsets(numR); // R1
        __hyperGVldt.build_n_RRsets(numR); // R2
        std::cout << "# RR sets: " << numR << std::endl;
        __numRRsets = __hyperG.get_RR_sets_size();
        time1 += timerOPIMC.get_operation_time();
        const auto infSelf = max_cover_budget(graph, budget, budget_list, timeb, mode);
        time2 += timerOPIMC.get_operation_time();
        const auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        const auto degVldt = infVldt * __numRRsets / __numV;
        auto upperBound = infSelf / approx;
        auto upperBound_med = infSelf / approx;
        if (mode == 1) {
            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
                upperBound = min_bound(upperBound, upperBound_med);
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        std::cout << "========upper bound ratio: " << upperBound / upperBound_med << "========" << std::endl;
        const auto upperDegOPT = upperBound * __numRRsets / __numV;
        auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-B-APPROX (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << approxOPIMC <<
                  " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        // Check whether the requirement is satisfied
        if (approxOPIMC >= approx && satisfy_count == 1)
        {
            __tRes.set_approximation(approxOPIMC);
            __tRes.set_running_time(timerOPIMC.get_total_time());
            __tRes.set_influence(infVldt);
            __tRes.set_influence_original(infSelf);
            __tRes.set_seed_vec(__vecSeed);
            __tRes.set_seed_inf_vec(__vecSeed_inf);
            __tRes.set_seed_inf_cost_vec(__vecSeed_inf_cost);
            __tRes.set_boundmin_inf(__boundMin_inf);
            __tRes.set_boundmin_inf_cost(__boundMin_inf_cost);
            __tRes.set_RR_sets_size(__numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 - timeb << '\n';
            std::cout << "==>Time for Bound Computation: " << timeb << '\n';
            return 0;
        }
        else if (approxOPIMC >= approx && satisfy_count == 0) {
            satisfy_count = satisfy_count + 1;
        }
        idx = idx + 1;
    }
}

double Alg::opimb_naive(const Graph & graph, const double budget, const std::vector<double> & budget_list,
        const double epsilon, const double delta, const int mode) {
    Timer timerOPIMC("OPIM-B-NAIVE");
    std::cout << "mode: " << mode << std::endl;
    const double e_alpha = exp(1);
    const double approx = 0.5 * (1 - 1.0 / e_alpha);
    const double alpha = sqrt(log(6.0 / delta));
    const int k_max = cal_k_max(budget_list, budget);
    const int k_min = cal_k_min(budget_list, budget);
    int targetSize_min = int(budget / max_vector(budget_list));
    int targetSize_max = int(budget / min_vector_nonzero(budget_list));
    std::cout << "k_max: " << k_max << std::endl;
    std::cout << "k_min: " << k_min << std::endl;
    std::cout << "targetSize_max: " << targetSize_max << std::endl;
    std::cout << "targetSize_min: " << targetSize_min << std::endl;
    const auto theta_zero = (size_t) ((2.0 * __numV * pow2(0.5 * (1 - 1 / e_alpha) * alpha
            + sqrt(0.5 * (1 - 1 / e_alpha) * (k_max * log(__numV) + log(6.0 / delta))))) / (__numV));
    const auto theta_max = (size_t) ((2.0 * __numV * pow2(0.5 * (1 - 1 / e_alpha) * alpha
            + sqrt(0.5 * (1 - 1 / e_alpha) * (k_max * log(__numV) + log(6.0 / delta))))) / (pow2(epsilon) * k_max));
    const auto numIter = (size_t) log2(theta_max / theta_zero) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0;
    double timeb = 0.0;
    size_t idx = 0;
    int satisfy_count = 0;
    while (true) {
        std::cout << idx << std::endl;
        const auto numR = theta_zero << idx;
        timerOPIMC.get_operation_time();
        __hyperG.build_n_RRsets(numR); // R1
        __hyperGVldt.build_n_RRsets(numR); // R2
        std::cout << "# RR sets: " << numR << std::endl;
        __numRRsets = __hyperG.get_RR_sets_size();
        time1 += timerOPIMC.get_operation_time();
        const auto infSelf = max_cover_budget(graph, budget, budget_list, timeb, mode);
        time2 += timerOPIMC.get_operation_time();
        const auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        const auto degVldt = infVldt * __numRRsets / __numV;
        auto upperBound = infSelf / approx;
        const auto upperDegOPT = upperBound * __numRRsets / __numV;
        auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-B-NAIVE (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << approxOPIMC <<
                  " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        // Check whether the requirement is satisfied
        if (approxOPIMC >= approx - epsilon && satisfy_count == 1)
        {
            __tRes.set_approximation(approxOPIMC);
            __tRes.set_running_time(timerOPIMC.get_total_time());
            __tRes.set_influence(infVldt);
            __tRes.set_influence_original(infSelf);
            __tRes.set_seed_vec(__vecSeed);
            __tRes.set_seed_inf_vec(__vecSeed_inf);
            __tRes.set_seed_inf_cost_vec(__vecSeed_inf_cost);
            __tRes.set_boundmin_inf(__boundMin_inf);
            __tRes.set_boundmin_inf_cost(__boundMin_inf_cost);
            __tRes.set_RR_sets_size(__numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 - timeb << '\n';
            std::cout << "==>Time for Bound Computation: " << timeb << '\n';
            return 0;
        }
        else if (approxOPIMC >= approx - epsilon && satisfy_count == 0) {
            satisfy_count = satisfy_count + 1;
        }
        idx = idx + 1;
    }
}

double Alg::opimb_fast_end(const Graph & graph, const double budget, const std::vector<double> & budget_list,
        const double epsilon, const double feps, const double delta, const std::string model,
        std::string graphname, const int mode) {
    Timer timerOPIMC("OPIM-B-FAST-END");
    std::ofstream outfile;
    double e_alpha;
    if (fabs(feps - 0.01) <= 0.000001) {
        e_alpha = exp(0.437888 * (1 - feps));
    }
    else if (fabs(feps - 0.05) <= 0.000001) {
        e_alpha = exp(0.435641 * (1 - feps));
    }
    else if (fabs(feps - 0.1) <= 0.000001) {
        e_alpha = exp(0.432857 * (1 - feps));
    }
    else if (fabs(feps - 0.15) <= 0.000001) {
        e_alpha = exp(0.430099 * (1 - feps));
    }
    else if (fabs(feps - 0.2) <= 0.000001) {
        e_alpha = exp(0.427369 * (1 - feps));
    }
    else {
        throw "feps error!";
    }
    const double approx = 1 - 1.0 / e_alpha;
    std::cout << "Approx: " << approx << std::endl;
    const double alpha = sqrt(log(6.0 / delta));
    const int k_max = cal_k_max(budget_list, budget);
    const int k_min = cal_k_min(budget_list, budget);
    const auto theta_zero = (size_t) ((2.0 * __numV * pow2((1 - 1 / e_alpha) * alpha
            + sqrt((1 - 1 / e_alpha) * (k_max * log(__numV) + log(6.0 / delta))))) / (__numV));
    const auto theta_max = (size_t) ((2.0 * __numV * pow2((1 - 1 / e_alpha) * alpha
            + sqrt((1 - 1 / e_alpha) * (k_max * log(__numV) + log(6.0 / delta))))) / (pow2(epsilon) * k_max));
    const auto numIter = (size_t) log2(theta_max / theta_zero) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0;
    double timeb = 0.0;
    int satisfy_count = 0;
    size_t idx = 0;
    while (true) {
        std::cout << idx << std::endl;
        const auto numR = theta_zero << idx;
        timerOPIMC.get_operation_time();
        __hyperG.build_n_RRsets(numR); // R1
        __hyperGVldt.build_n_RRsets(numR); // R2
        std::cout << "# RR sets: " << numR << std::endl;
        __numRRsets = __hyperG.get_RR_sets_size();
        time1 += timerOPIMC.get_operation_time();
        const auto infSelf = max_cover_budget_fast_end(graph, budget, budget_list, feps, timeb, mode);
        time2 += timerOPIMC.get_operation_time();
        const auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        const auto degVldt = infVldt * __numRRsets / __numV;
        auto upperBound = infSelf;
        auto upperBound_med = infSelf / approx;
        if (mode == 1) {
//            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
                upperBound = min_bound(upperBound, upperBound_med);
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
                upperBound = min_bound(upperBound, upperBound_med);
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        std::cout << "========upper bound ratio: " << upperBound / upperBound_med << "========" << std::endl;
        const auto upperDegOPT = upperBound * __numRRsets / __numV;
        auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-B-FAST-END (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " <<
            approxOPIMC << " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        // Check whether the requirement is satisfied
        if (approxOPIMC >= approx - epsilon && satisfy_count == 1)
        {
            __tRes.set_approximation(approxOPIMC);
            __tRes.set_running_time(timerOPIMC.get_total_time());
            __tRes.set_influence(infVldt);
            __tRes.set_influence_original(infSelf);
            __tRes.set_seed_vec(__vecSeed);
            __tRes.set_seed_inf_vec(__vecSeed_inf);
            __tRes.set_seed_inf_cost_vec(__vecSeed_inf_cost);
            __tRes.set_boundmin_inf(__boundMin_inf);
            __tRes.set_boundmin_inf_cost(__boundMin_inf_cost);
            __tRes.set_RR_sets_size(__numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 - timeb << '\n';
            std::cout << "==>Time for Bound Computation: " << timeb << '\n';
            std::cout << "Seed size: " << __vecSeed.size() << std::endl;
            time_t now = time(0);
            char* dt = ctime(&now);
            outfile.open("fast/" + model + "_" + graphname + "_" + std::to_string(epsilon) + "_" + std::to_string(feps)
                         + "_" + dt + ".txt");
            outfile << __tRes.get_running_time() << "," << "," << __numRRsets * 2 << "," <<
            infVldt << "," << "," << infSelf << "," << "," << __vecSeed.size() << "," <<
            approxOPIMC << "," << "," << time1 << "," << timeb << "," << time2 - timeb;
            return 0;
        }
        else if (approxOPIMC >= approx - epsilon && satisfy_count == 0) {
            satisfy_count = satisfy_count + 1;
        }
        idx = idx + 1;
    }
}

double Alg::opim_approx(const Graph & graph, const double budget, const std::vector<double> & budget_list,
        const double epsilon, const double feps, const double delta, const std::string model,
        const std::string graphname, const int mode) {
    Timer timerOPIMC("OPIM-APPORX");
    std::ofstream outfile;
    const auto numIter = 12;
    const auto theta_0 = 1000;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1_n = 0.0, time2_n = 0.0;
    double timeb_n = 0.0;
    double time1_b = 0.0, time2_b = 0.0;
    double timeb_b = 0.0;
    double time1_f = 0.0, time2_f = 0.0;
    double timeb_f = 0.0;
    double time_rr = 0.0;
    std::vector<double> inf_n, inf_b, inf_f;
    std::vector<double> approx_n, approx_b, approx_f;
    std::vector<double> time_n, time_b, time_f;
    size_t idx = 0;
    while (true) {
        std::cout << idx << std::endl;
        const auto numR = theta_0 << idx;
        timerOPIMC.get_operation_time();
        __hyperG.build_n_RRsets(numR); // R1
        __hyperGVldt.build_n_RRsets(numR); // R2
        std::cout << "# RR sets: " << numR << std::endl;
        __numRRsets = __hyperG.get_RR_sets_size();
        /// OPIM-B-N
        time_rr += timerOPIMC.get_operation_time();
        auto infSelf = max_cover_budget(graph, budget, budget_list, timeb_n, 0);
        time2_n += timerOPIMC.get_operation_time();
        auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        auto degVldt = infVldt * __numRRsets / __numV;
        auto upperBound = infSelf / (0.5 * (1 - 1.0 / exp(1)));
        auto upperDegOPT = upperBound * __numRRsets / __numV;
        auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        auto approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-APPROX-N (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " <<
                  approxOPIMC << " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
        std::cout << "==>Time for RR sets and greedy: " << time_rr << ", " << time2_n + time_rr << '\n';
        std::cout << "==>Time for Bound Computation: " << timeb_n << '\n';
        std::cout << "Seed size: " << __vecSeed.size() << std::endl;
        approx_n.push_back(approxOPIMC);
        inf_n.push_back(infVldt);
        time_n.push_back(time2_n + time_rr);
        /// OPIM-B
        time1_b += timerOPIMC.get_operation_time();
        infSelf = max_cover_budget(graph, budget, budget_list, timeb_b, mode);
        time2_b += timerOPIMC.get_operation_time();
        infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        degVldt = infVldt * __numRRsets / __numV;
        upperBound = infSelf;
        if (mode == 1) {
            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        upperDegOPT = upperBound * __numRRsets / __numV;
        lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-APPROX-B (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " <<
                  approxOPIMC << " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
        std::cout << "==>Time for RR sets and greedy: " << time_rr << ", " << time2_b + time_rr << '\n';
        std::cout << "==>Time for Bound Computation: " << timeb_b << '\n';
        std::cout << "Seed size: " << __vecSeed.size() << std::endl;
        approx_b.push_back(approxOPIMC);
        inf_b.push_back(infVldt);
        time_b.push_back(time2_b + time_rr);
        /// OPIM-B-F
        time1_f += timerOPIMC.get_operation_time();
        infSelf = max_cover_budget_fast_end(graph, budget, budget_list, feps, timeb_f, mode);
        time2_f += timerOPIMC.get_operation_time();
        infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        degVldt = infVldt * __numRRsets / __numV;
        upperBound = infSelf;
        if (mode == 1) {
            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        upperDegOPT = upperBound * __numRRsets / __numV;
        lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-APPROX-F (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " <<
                  approxOPIMC << " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
        std::cout << "==>Time for RR sets and greedy: " << time_rr << ", " << time2_f + time_rr << '\n';
        std::cout << "==>Time for Bound Computation: " << timeb_f << '\n';
        std::cout << "Seed size: " << __vecSeed.size() << std::endl;
        approx_f.push_back(approxOPIMC);
        inf_f.push_back(infVldt);
        time_f.push_back(time2_f + time_rr);

        idx = idx + 1;
        if (idx == numIter) {
            break;
        }
    }
    time_t now = time(0);
    char* dt = ctime(&now);
    outfile.open("rr/" + model + "_" + graphname + "_" + dt + ".txt");
    for(int i = 0; i < approx_n.size(); i++) {
        outfile << approx_n[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < time_n.size(); i++) {
        outfile << time_n[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < inf_n.size(); i++) {
        outfile << inf_n[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < approx_b.size(); i++) {
        outfile << approx_b[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < time_b.size(); i++) {
        outfile << time_b[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < inf_b.size(); i++) {
        outfile << inf_b[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < approx_f.size(); i++) {
        outfile << approx_f[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < time_f.size(); i++) {
        outfile << time_f[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < inf_f.size(); i++) {
        outfile << inf_f[i] << ",";
    }
    outfile << std::endl;
    outfile.close();
}

double Alg::opim_approx_fast(const Graph & graph, const double budget, const std::vector<double> & budget_list,
        const double epsilon, const double feps, const double delta, const std::string model,
        const std::string graphname, const int mode) {
    Timer timerOPIMC("OPIM-APPORX");
    std::ofstream outfile;
    const auto numIter = 12;
    const auto theta_0 = 1000;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1_n = 0.0, time2_n = 0.0;
    double timeb_n = 0.0;
    double time1_b = 0.0, time2_b = 0.0;
    double timeb_b = 0.0;
    double time1_f = 0.0, time2_f = 0.0;
    double timeb_f = 0.0;
    double time1_fs = 0.0, time2_fs = 0.0;
    double timeb_fs = 0.0;
    double time1_ft = 0.0, time2_ft = 0.0;
    double timeb_ft = 0.0;
    double time_rr = 0.0;
    std::vector<double> inf_n, inf_b, inf_f;
    std::vector<double> approx_n, approx_b, approx_f;
    std::vector<double> time_n, time_b, time_f;
    size_t idx = 0;
    while (true) {
        std::cout << idx << std::endl;
        const auto numR = theta_0 << idx;
        timerOPIMC.get_operation_time();
        __hyperG.build_n_RRsets(numR); // R1
        __hyperGVldt.build_n_RRsets(numR); // R2
        std::cout << "# RR sets: " << numR << std::endl;
        __numRRsets = __hyperG.get_RR_sets_size();
        /// OPIM-B-F 0.01
        time_rr += timerOPIMC.get_operation_time();
        auto infSelf = max_cover_budget_fast_end(graph, budget, budget_list, 0.01, timeb_n, mode);
        time2_n += timerOPIMC.get_operation_time();
        auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        auto degVldt = infVldt * __numRRsets / __numV;
        auto upperBound = infSelf;
        if (mode == 1) {
            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        auto upperDegOPT = upperBound * __numRRsets / __numV;
        auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        auto approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-APPROX-F-0.01 (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " <<
                  approxOPIMC << " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
        std::cout << "==>Time for RR sets and greedy: " << time_rr << ", " << time2_n + time_rr << '\n';
        std::cout << "==>Time for Bound Computation: " << timeb_n << '\n';
        std::cout << "Seed size: " << __vecSeed.size() << std::endl;
        approx_n.push_back(approxOPIMC);
        inf_n.push_back(infVldt);
        time_n.push_back(time2_n + time_rr);

        /// OPIM-B-F 0.05
        time1_b += timerOPIMC.get_operation_time();
        infSelf = max_cover_budget_fast_end(graph, budget, budget_list, 0.05, timeb_b, mode);
        time2_b += timerOPIMC.get_operation_time();
        infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        degVldt = infVldt * __numRRsets / __numV;
        upperBound = infSelf;
        if (mode == 1) {
            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        upperDegOPT = upperBound * __numRRsets / __numV;
        lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-APPROX-F-0.05 (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " <<
                  approxOPIMC << " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
        std::cout << "==>Time for RR sets and greedy: " << time_rr << ", " << time2_b + time_rr << '\n';
        std::cout << "==>Time for Bound Computation: " << timeb_b << '\n';
        std::cout << "Seed size: " << __vecSeed.size() << std::endl;
        approx_b.push_back(approxOPIMC);
        inf_b.push_back(infVldt);
        time_b.push_back(time2_b + time_rr);

        /// OPIM-B-F 0.1
        time1_f += timerOPIMC.get_operation_time();
        infSelf = max_cover_budget_fast_end(graph, budget, budget_list, 0.1, timeb_f, mode);
        time2_f += timerOPIMC.get_operation_time();
        infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        degVldt = infVldt * __numRRsets / __numV;
        upperBound = infSelf;
        if (mode == 1) {
            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        upperDegOPT = upperBound * __numRRsets / __numV;
        lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-APPROX-F-0.1 (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " <<
                  approxOPIMC << " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
        std::cout << "==>Time for RR sets and greedy: " << time_rr << ", " << time2_f + time_rr << '\n';
        std::cout << "==>Time for Bound Computation: " << timeb_f << '\n';
        std::cout << "Seed size: " << __vecSeed.size() << std::endl;
        approx_f.push_back(approxOPIMC);
        inf_f.push_back(infVldt);
        time_f.push_back(time2_f + time_rr);
        idx = idx + 1;
        if (idx == numIter) {
            break;
        }
    }
    time_t now = time(0);
    char* dt = ctime(&now);
    outfile.open("xi/" + model + "_" + graphname + "_" + dt + ".txt");
    for(int i = 0; i < approx_n.size(); i++) {
        outfile << approx_n[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < time_n.size(); i++) {
        outfile << time_n[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < inf_n.size(); i++) {
        outfile << inf_n[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < approx_b.size(); i++) {
        outfile << approx_b[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < time_b.size(); i++) {
        outfile << time_b[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < inf_b.size(); i++) {
        outfile << inf_b[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < approx_f.size(); i++) {
        outfile << approx_f[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < time_f.size(); i++) {
        outfile << time_f[i] << ",";
    }
    outfile << std::endl;
    for(int i = 0; i < inf_f.size(); i++) {
        outfile << inf_f[i] << ",";
    }
    outfile << std::endl;
    outfile.close();
}

double Alg::opimb_fast_approx(const Graph & graph, const double budget, const std::vector<double> & budget_list,
                           const double epsilon, const double feps, const double delta, const int mode) {
    Timer timerOPIMC("OPIM-B-FAST-APPORX");
    const double approx = epsilon;
    std::cout << "Approx: " << approx << std::endl;
    const double alpha = sqrt(log(6.0 / delta));
    const int k_max = cal_k_max(budget_list, budget);
    const int k_min = cal_k_min(budget_list, budget);
    const auto theta_zero = (size_t) ((2.0 * __numV * pow2(approx * alpha
            + sqrt(approx * (k_max * log(__numV) + log(6.0 / delta))))) / (__numV));
    const auto theta_max = (size_t) ((2.0 * __numV * pow2(approx * alpha
            + sqrt(approx * (k_max * log(__numV) + log(6.0 / delta))))) / (pow2(epsilon) * k_max));
    const auto numIter = (size_t) log2(theta_max / theta_zero) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0;
    double timeb = 0.0;
    int satisfy_count = 0;
    size_t idx = 0;
    while (true) {
        std::cout << idx << std::endl;
        const auto numR = theta_zero << idx;
        timerOPIMC.get_operation_time();
        __hyperG.build_n_RRsets(numR); // R1
        __hyperGVldt.build_n_RRsets(numR); // R2
        std::cout << "# RR sets: " << numR << std::endl;
        __numRRsets = __hyperG.get_RR_sets_size();
        time1 += timerOPIMC.get_operation_time();
        const auto infSelf = max_cover_budget_fast_end(graph, budget, budget_list, feps, timeb, mode);
        time2 += timerOPIMC.get_operation_time();
        const auto infVldt = __hyperGVldt.self_inf_cal(__vecSeed);
        const auto degVldt = infVldt * __numRRsets / __numV;
        auto upperBound = infSelf;
        auto upperBound_med = infSelf / approx;
        if (mode == 1) {
            upperBound = max_element(__boundLast_inf, __boundLast_inf_cost);
        }
        else if (mode == 2) {
            if (__is_inf_cost) {
                upperBound = __boundMin_inf_cost;
                upperBound = min_bound(upperBound, upperBound_med);
            }
            else {
                upperBound = __boundMin_inf;
            }
        }
        std::cout << "========upper bound ratio: " << upperBound / upperBound_med << "========" << std::endl;
        const auto upperDegOPT = upperBound * __numRRsets / __numV;
        auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        if (lowerSelect < 0) {
            lowerSelect = 0;
        }
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto approxOPIMC = lowerSelect / upperOPT;
        std::cout << " -->OPIM-B-FAST-APPROX (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " <<
                  approxOPIMC << " (" << infSelf / upperBound << "), #RR sets: " << __numRRsets << '\n';
        // Check whether the requirement is satisfied
        if (approxOPIMC >= approx && satisfy_count == 1)
        {
            __tRes.set_approximation(approxOPIMC);
            __tRes.set_running_time(timerOPIMC.get_total_time());
            __tRes.set_influence(infVldt);
            __tRes.set_influence_original(infSelf);
            __tRes.set_seed_vec(__vecSeed);
            __tRes.set_seed_inf_vec(__vecSeed_inf);
            __tRes.set_seed_inf_cost_vec(__vecSeed_inf_cost);
            __tRes.set_boundmin_inf(__boundMin_inf);
            __tRes.set_boundmin_inf_cost(__boundMin_inf_cost);
            __tRes.set_RR_sets_size(__numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << __tRes.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 - timeb << '\n';
            std::cout << "==>Time for Bound Computation: " << timeb << '\n';
            std::cout << "Seed size: " << __vecSeed.size() << std::endl;
            return 0;
        }
        else if (approxOPIMC >= approx && satisfy_count == 0) {
            satisfy_count = satisfy_count + 1;
        }
        idx = idx + 1;
    }
}
