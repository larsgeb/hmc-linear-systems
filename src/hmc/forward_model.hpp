//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP
#define HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP

#include "AlgebraLib/src/algebra_lib/algebra_lib.hpp"

using namespace algebra_lib;

namespace hmc {
    class forward_model {
    public:
        explicit forward_model(const matrix &_model_matrix);

        const matrix _g;

    };
}

#endif //HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP
