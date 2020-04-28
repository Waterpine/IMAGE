//
// Created by sbian on 2019/9/9.
//

#ifndef BIM_CELF_H
#define BIM_CELF_H

#include "graphBase.h"
#include <ctime>

class CELF {
public:
    static double sum_budget(const std::vector<double> & budget_list, const std::vector<size_t> & vecSeed) {
        double sum = 0;
        for(size_t i = 0; i < vecSeed.size(); i++) {
            sum += budget_list[vecSeed[i]];
        }
        return sum;
    }
    static std::vector<size_t> celfalgo(const Graph& graph, const double budget,
            const std::vector<double> & budget_list, const CascadeModel model) {
        std::vector<size_t> S_UC;
        std::vector<size_t> S_CB;
        std::deque<std::tuple<size_t, double, double, size_t>> Q_UC;  // u u.mg u.cost u.iter
        std::deque<std::tuple<size_t, double, double, size_t>> Q_CB;  // u u.mg/u.cost u.cost u.iter

        std::vector<size_t> vecSeed(1);

        clock_t startTime, endTime;
        startTime = clock();
        std::cout << "Initial deque start" << std::endl;
        for (size_t i = 0; i < graph.size(); i++){
            vecSeed[0] = i;
            double mg = GraphBase::inf_eval(graph, vecSeed, model);
            Q_UC.push_back(std::make_tuple(i, mg, budget_list[i], 0));
            Q_CB.push_back(std::make_tuple(i, mg / budget_list[i], budget_list[i], 0));
        }
        std::cout << "Initial deque end" << std::endl;
        endTime = clock();
        std::cout << double((endTime - startTime) / CLOCKS_PER_SEC) << " s" << std::endl;

        std::sort(Q_UC.begin(), Q_UC.end(), [](auto & left, auto & right) {
            if (std::get<1>(left) == std::get<1>(right)) {
                return std::get<2>(left) < std::get<2>(right);
            }
            else {
                return std::get<1>(left) > std::get<1>(right);
            }
        });

        std::sort(Q_CB.begin(), Q_CB.end(), [](auto & left, auto & right) {
            if (std::get<1>(left) == std::get<1>(right)) {
                return std::get<2>(left) < std::get<2>(right);
            }
            else {
                return std::get<1>(left) > std::get<1>(right);
            }
        });

        std::cout << "Start UC" << std::endl;
        double budget_UC = 0;
        double S_UC_inf = 0;
        while(budget_UC < budget) {
            while(std::get<2>(Q_UC.front()) + budget_UC > budget){  // u u.mg u.cost u.iter
                Q_UC.pop_front();
                if(Q_UC.empty()){
                    break;
                }
            }
            if(Q_UC.empty()){
                break;
            }
            std::tuple<size_t, double, double, size_t> u = Q_UC.front();  // u u.mg u.cost u.iter
            if(std::get<3>(u) == S_UC.size()) {
                S_UC.push_back(std::get<0>(u));
                Q_UC.pop_front();
                budget_UC = budget_UC + std::get<2>(u);
                S_UC_inf = GraphBase::inf_eval(graph, S_UC, model);
                std::cout << "UC push new seed: " << std::get<0>(u) << " " << S_UC.size() << " " <<
                sum_budget(budget_list, S_UC) << " " << S_UC_inf << std::endl;
            }
            else {
                S_UC.push_back(std::get<0>(u));
                std::get<1>(Q_UC.front()) = GraphBase::inf_eval(graph, S_UC, model) - S_UC_inf;
                S_UC.pop_back();
                std::get<3>(Q_UC.front()) = S_UC.size();
                std::sort(Q_UC.begin(), Q_UC.end(), [](auto & left, auto & right) {
                    if (std::get<1>(left) == std::get<1>(right)) {
                        return std::get<2>(left) < std::get<2>(right);
                    }
                    else {
                        return std::get<1>(left) > std::get<1>(right);
                    }
                });
            }
        }
        std::cout << "End UC" << std::endl;

        std::cout << "Start CB" << std::endl;
        double budget_CB = 0;
        double S_CB_inf = 0;
        while(budget_CB < budget) {
            while(std::get<2>(Q_CB.front()) + budget_CB > budget){
                Q_CB.pop_front();
                if(Q_CB.empty()){
                    break;
                }
            }
            if(Q_CB.empty()) {
                break;
            }
            std::tuple<size_t, double, double, size_t> u = Q_CB.front();  // u u.mg/u.cost u.cost u.iter
            if(std::get<3>(u) == S_CB.size()) {
                S_CB.push_back(std::get<0>(u));
                Q_CB.pop_front();
                budget_CB = budget_CB + std::get<2>(u);
                S_CB_inf = GraphBase::inf_eval(graph, S_CB, model);
                std::cout << "CB push new seed: " << std::get<0>(u) << " " << S_CB.size() << " " <<
                sum_budget(budget_list, S_CB) << " " << S_CB_inf << std::endl;
            }
            else {
                S_CB.push_back(std::get<0>(u));
                std::get<1>(Q_CB.front()) = (GraphBase::inf_eval(graph, S_CB, model) - S_CB_inf) / std::get<2>(u);
                S_CB.pop_back();
                std::get<3>(Q_CB.front()) = S_CB.size();
                std::sort(Q_CB.begin(), Q_CB.end(), [](auto & left, auto & right) {
                    if (std::get<1>(left) == std::get<1>(right)) {
                        return std::get<2>(left) < std::get<2>(right);
                    }
                    else {
                        return std::get<1>(left) > std::get<1>(right);
                    }
                });
            }
        }
        std::cout << "CB end" << std::endl;

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

#endif //BIM_CELF_H
