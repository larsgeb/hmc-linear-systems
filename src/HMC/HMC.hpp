//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP
#define HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP

namespace hmc {
    algebra_lib::matrix std_to_inv_cov(algebra_lib::vector stdv) {
        algebra_lib::matrix inv_cov(stdv.size(), stdv.size());

        for (auto it = stdv.begin(); it != stdv.end(); it++) {
            inv_cov[it - stdv.begin()][it - stdv.begin()] = 1.0 / (*it * *it);
        }
        return inv_cov;
    }
}

#include "forwardmodel.hpp"
#include "prior.hpp"
#include "data.hpp"
#include "posterior.hpp"
#include "sampler.hpp"

#endif //HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP
