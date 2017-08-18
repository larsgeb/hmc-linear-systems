//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_POSTERIOR_HPP
#define HMC_LINEAR_SYSTEM_POSTERIOR_HPP

#include <SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp>
#include "Prior.hpp"
#include "Data.hpp"

namespace HMC {
    class Posterior {
        double misfit(AlgebraLib::Vector &parameters, Prior &prior, Data &data);

        AlgebraLib::Vector gradient_misfit(AlgebraLib::Vector &parameters, Prior &prior, Data &data);
    };
}

#endif //HMC_LINEAR_SYSTEM_POSTERIOR_HPP
