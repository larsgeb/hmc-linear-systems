//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP
#define HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP

#include <SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp>

namespace HMC {
    class ForwardModel {
    public:
        explicit ForwardModel(const AlgebraLib::Matrix &_model_matrix);

        const AlgebraLib::Matrix _g;

    };
}

#endif //HMC_LINEAR_SYSTEM_FORWARDMODEL_HPP
