//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP
#define HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP

#include "armadillo"

using namespace arma;

namespace hmc {
    class forward_model {
    public:
        explicit forward_model(const mat &_model_matrix);

        forward_model();

        arma::mat _g;

    };
}

#endif //HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP
