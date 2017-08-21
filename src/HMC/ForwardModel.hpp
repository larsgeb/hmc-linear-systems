//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP
#define HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP

#include <algebra_lib/src/algebra_lib/algebra_lib.hpp>

namespace hmc {
    class ForwardModel {
    public:
        explicit ForwardModel(const algebra_lib::matrix &_model_matrix);

        const algebra_lib::matrix _g;

    };
}

#endif //HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP
