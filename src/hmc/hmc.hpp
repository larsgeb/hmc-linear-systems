//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP
#define HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP

#include <armadillo>

namespace hmc {

    arma::mat std_to_inv_cov(const arma::dcolvec &stdv) {
        arma::mat inv_cov = diagmat(arma::pow(arma::square(stdv), -1));
        return inv_cov;
    }
}


#include "linearSampler.hpp"

#endif //HMC_LINEAR_SYSTEM_HAMILTONIANMONTECARLO_HPP_HPP
