//
// Created by sbian on 2019/9/9.
//

#ifndef BIM_CELFPP_H
#define BIM_CELFPP_H

#include "graphBase.h"
#include <ctime>

class CELFPP {
public:
    static double sum_budget(const std::vector<double> & budget_list, const std::vector<size_t> & vecSeed) {
        double sum = 0;
        for (size_t i = 0; i < vecSeed.size(); i++) {
            sum += budget_list[vecSeed[i]];
        }
        return sum;
    }
    static std::vector<size_t> celfppalgo (const Graph & graph, const double budget,
            const std::vector<double> & budget_list, const CascadeModel model) {
        std::vector<size_t> S_UC;
        std::vector<size_t> S_CB;
        // u u.mg1 u.prev_best u.mg2 u.flag u.cost
        std::deque<std::tuple<size_t, double, int, double, size_t, double>> Q_UC;
        // u u.mg1/u.cost u.prev_best u.mg2/(u.cost + u.prev_best.cost) u.flag u.cost
        std::deque<std::tuple<size_t, double, int, double, size_t, double>> Q_CB;

        int last_seed = -1;
        int cur_best_UC = -1;
        double cur_best_UC_mg1 = 0;
        int cur_best_CB = -1;
        double cur_best_CB_mg1 = 0;
        double mg1_UC = 0;
        double mg1_CB = 0;
        double mg2_UC = 0;
        double mg2_CB = 0;
        std::cout << "Initial deque start" << std::endl;
        std::vector<size_t> vecSeed;
        for (size_t i = 0; i < graph.size(); i++) {
            vecSeed.push_back(i);

            // initial UC
            mg1_UC = GraphBase::inf_eval(graph, vecSeed, model);
            for (size_t j = 0; j < vecSeed.size(); j++) {
                std::cout << vecSeed[j] << " ";
            }
            std::cout << std::endl;
            if (cur_best_UC != -1) {
                vecSeed.push_back(cur_best_UC);
                mg2_UC = GraphBase::inf_eval(graph, vecSeed, model);
                for (size_t j = 0; j < vecSeed.size(); j++) {
                    std::cout << vecSeed[j] << " ";
                }
                std::cout << std::endl;
                vecSeed.pop_back();
            }
            else {
                mg2_UC = mg1_UC;
            }
            // u u.mg1 u.prev_best u.mg2 u.flag u.cost
            Q_UC.push_back(std::make_tuple(i, mg1_UC, cur_best_UC, mg2_UC, 0, budget_list[i]));
            if (cur_best_UC == -1) {
                cur_best_UC = i;
                cur_best_UC_mg1 = mg1_UC;
            }
            else {
                if(cur_best_UC_mg1 < mg1_UC) {
                    cur_best_UC = i;
                    cur_best_UC_mg1 = mg1_UC;
                }
            }

            // initial CB
            mg1_CB = mg1_UC / budget_list[i];
            if (cur_best_CB != -1) {
                vecSeed.push_back(cur_best_CB);
                mg2_CB = GraphBase::inf_eval(graph, vecSeed, model) / (budget_list[i] + budget_list[cur_best_CB]);
                for (size_t j = 0; j < vecSeed.size(); j++) {
                    std::cout << vecSeed[j] << " ";
                }
                std::cout << std::endl;
                vecSeed.pop_back();
            }
            else {
                mg2_CB = mg1_CB;
            }

            // u u.mg1/u.cost u.prev_best u.mg2/(u.cost + u.prev_best.cost) u.flag u.cost
            Q_CB.push_back(std::make_tuple(i, mg1_CB, cur_best_CB, mg2_CB, 0, budget_list[i]));
            if (cur_best_CB == -1) {
                cur_best_CB = i;
                cur_best_CB_mg1 = mg1_CB;
            }
            else {
                if (cur_best_CB_mg1 < mg1_CB) {
                    cur_best_CB = i;
                    cur_best_CB_mg1 = mg1_CB;
                }
            }
            vecSeed.pop_back();
        }
        std::cout << "Initial deque end" << std::endl;

        // u u.mg1 u.prev_best u.mg2 u.flag u.cost
        std::sort(Q_UC.begin(), Q_UC.end(), [](auto & left, auto & right) {
            if (std::get<1>(left) == std::get<1>(right)) {
                return std::get<5>(left) < std::get<5>(right);
            }
            else {
                return std::get<1>(left) > std::get<1>(right);
            }
        });

        // u u.mg1/u.cost u.prev_best u.mg2/(u.cost + u.prev_best.cost) u.flag u.cost
        std::sort(Q_CB.begin(), Q_CB.end(), [](auto & left, auto & right) {
            if (std::get<1>(left) == std::get<1>(right)) {
                return std::get<5>(left) < std::get<5>(right);
            }
            else {
                return std::get<1>(left) > std::get<1>(right);
            }
        });

        // u u.mg1 u.prev_best u.mg2 u.flag u.cost
        std::cout << "Start UC" << std::endl;
        size_t budget_UC = 0;
        double S_UC_inf = 0;
        while (budget_UC < budget) {
            while (std::get<5>(Q_UC.front()) + budget_UC > budget) {
                Q_UC.pop_front();
                if (Q_UC.empty()) {
                    break;
                }
            }
            if (Q_UC.empty()) {
                break;
            }
            // u u.mg1 u.prev_best u.mg2 u.flag u.cost
            std::tuple<size_t, double, int, double, size_t, double> u = Q_UC.front();
            if (std::get<4>(u) == S_UC.size()) {
                S_UC.push_back(std::get<0>(u));
                Q_UC.pop_front();
                budget_UC = budget_UC + std::get<5>(u);
                S_UC_inf = GraphBase::inf_eval(graph, S_UC, model);
                last_seed = std::get<0>(u);
                std::cout << "UC push new seed: " << std::get<0>(u) << " " << S_UC.size() << " " <<
                          sum_budget(budget_list, S_UC) << " " << S_UC_inf << std::endl;
                continue;
            }
            else if (std::get<2>(u) == last_seed) {
                std::get<1>(Q_UC.front()) = std::get<3>(Q_UC.front());   // u u.mg1 u.prev_best u.mg2 u.flag u.cost
                std::cout << "last seed" << std::endl;
            }
            else {
                S_UC.push_back(std::get<0>(u));
                std::get<1>(Q_UC.front()) = GraphBase::inf_eval(graph, S_UC, model) - S_UC_inf;
                S_UC.pop_back();
                std::get<2>(Q_UC.front()) = cur_best_UC;
                S_UC.push_back(cur_best_UC);
                S_UC.push_back(std::get<0>(u));
                std::get<3>(Q_UC.front()) = GraphBase::inf_eval(graph, S_UC, model) - S_UC_inf;
                S_UC.pop_back();
                S_UC.pop_back();
            }

            if (cur_best_UC != -1){
                if (cur_best_UC_mg1 < std::get<1>(Q_UC.front())) {
                    cur_best_UC = std::get<0>(Q_UC.front());
                }
            }
            else {
                throw "cur_best_UC should not be -1";
            }

            std::get<4>(Q_UC.front()) = S_UC.size();
            std::sort(Q_UC.begin(), Q_UC.end(), [](auto & left, auto & right) {
                if (std::get<1>(left) == std::get<1>(right)) {
                    return std::get<5>(left) < std::get<5>(right);
                }
                else {
                    return std::get<1>(left) > std::get<1>(right);
                }
            });
        }
        std::cout << "End UC" << std::endl;

        // u u.mg1/u.cost u.prev_best u.mg2/(u.cost + u.prev_best.cost) u.flag u.cost
        std::cout << "Start CB" << std::endl;
        size_t budget_CB = 0;
        double S_CB_inf = 0;
        while (budget_CB < budget) {
            while (std::get<5>(Q_CB.front()) + budget_CB > budget) {
                Q_CB.pop_front();
                if (Q_CB.empty()) {
                    break;
                }
            }
            if (Q_CB.empty()) {
                break;
            }
            // u u.mg1/u.cost u.prev_best u.mg2/(u.cost + u.prev_best.cost) u.flag u.cost
            std::tuple<size_t, double, int, double, size_t, double> u = Q_CB.front();
            if (std::get<4>(u) == S_CB.size()) {
                S_CB.push_back(std::get<0>(u));
                Q_CB.pop_front();
                budget_CB = budget_CB + std::get<5>(u);
                S_CB_inf = GraphBase::inf_eval(graph, S_CB, model);
                last_seed = std::get<0>(u);
                std::cout << "UC push new seed: " << std::get<0>(u) << " " << S_CB.size() << " " <<
                          sum_budget(budget_list, S_CB) << " " << S_CB_inf << std::endl;
                continue;
            }
            else if (std::get<2>(u) == last_seed) {
                // u u.mg1/u.cost u.prev_best u.mg2/(u.cost + u.prev_best.cost) u.flag u.cost
                std::get<1>(Q_CB.front()) = std::get<3>(Q_CB.front());
            }
            else {
                S_CB.push_back(std::get<0>(u));
                std::get<1>(Q_CB.front()) = (GraphBase::inf_eval(graph, S_CB, model) - S_CB_inf) / std::get<5>(u);
                S_CB.pop_back();
                std::get<2>(Q_CB.front()) = cur_best_CB;
                S_CB.push_back(cur_best_CB);
                S_CB.push_back(std::get<0>(u));
                std::get<3>(Q_CB.front()) = (GraphBase::inf_eval(graph, S_CB, model) - S_CB_inf) /
                        (std::get<5>(u) + budget_list[cur_best_CB]);
                S_CB.pop_back();
                S_CB.pop_back();
            }
            if (cur_best_CB != -1){
                if(cur_best_CB_mg1 < std::get<1>(Q_CB.front())){
                    cur_best_CB = std::get<0>(Q_CB.front());
                }
            }
            else {
                throw "cur_best_CB should not be -1";
            }
            std::get<4>(Q_CB.front()) = S_CB.size();
            std::sort(Q_CB.begin(), Q_CB.end(), [](auto & left, auto & right) {
                if (std::get<1>(left) == std::get<1>(right)) {
                    return std::get<5>(left) < std::get<5>(right);
                }
                else {
                    return std::get<1>(left) > std::get<1>(right);
                }
            });
        }
        std::cout << "End CB" << std::endl;

        if (S_CB_inf > S_UC_inf) {
            std::cout << "CB_inf: " << S_CB_inf << std::endl;
            return S_CB;
        }
        else {
            std::cout << "UC_inf: " << S_UC_inf << std::endl;
            return S_UC;
        }
    }
};

#endif //BIM_CELFPP_H
