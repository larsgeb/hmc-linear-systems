//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP
#define HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP

namespace HMC {
    AlgebraLib::Matrix std_to_inv_cov(AlgebraLib::Vector stdv) {
        AlgebraLib::Matrix inv_cov(stdv.size(), stdv.size());

        for (auto it = stdv.begin(); it != stdv.end(); it++) {
            inv_cov[it - stdv.begin()][it - stdv.begin()] = 1.0 / (*it * *it);
        }
        return inv_cov;
    }
}

#include "ForwardModel.hpp"
#include "Prior.hpp"
#include "Data.hpp"
#include "Posterior.hpp"
#include "Sampler.hpp"

#endif //HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP
