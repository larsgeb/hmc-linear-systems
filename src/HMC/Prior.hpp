//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_PRIOR_HPP
#define HMC_LINEAR_SYSTEM_PRIOR_HPP

#include <SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp>

namespace HMC{
    class Prior {
    public:
        // Constructors
        Prior(AlgebraLib::Vector &mean, AlgebraLib::Vector &std);

        // Fields
        const AlgebraLib::Vector _means;
        const AlgebraLib::Matrix _covariance;
        const AlgebraLib::Matrix _inv_cov_m;

        // Member functions
        double misfit(AlgebraLib::Vector &parameters);
        AlgebraLib::Vector gradient_misfit(AlgebraLib::Vector &parameters);
    };
}



#endif //HMC_LINEAR_SYSTEM_PRIOR_HPP
