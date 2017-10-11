//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP
#define HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP

using namespace algebra_lib;

namespace hmc {
    sparse_matrix std_to_inv_cov(sparse_vector stdv) {
        sparse_matrix inv_cov(stdv.size(), stdv.size());

        for (int it = 0; it < stdv.size(); it++) {
            inv_cov(it)(it) = 1.0 / (stdv[it] * stdv[it]);
        }
        return inv_cov;
    }
}

#include "forward_model.hpp"
#include "prior.hpp"
#include "data.hpp"
#include "posterior.hpp"
#include "sampler.hpp"

#endif //HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP
