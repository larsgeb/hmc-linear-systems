//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_POSTERIOR_HPP
#define HMC_LINEAR_SYSTEM_POSTERIOR_HPP

#include "prior.hpp"
#include "data.hpp"

namespace hmc {
    class Posterior {
        double misfit(arma::vec &parameters, prior &prior, data &data);

        arma::vec gradient_misfit(arma::vec &parameters, prior &prior, data &data);
    };
}

#endif //HMC_LINEAR_SYSTEM_POSTERIOR_HPP
