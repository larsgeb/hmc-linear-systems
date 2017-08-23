//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_PRIOR_HPP
#define HMC_LINEAR_SYSTEM_PRIOR_HPP

#include <AlgebraLib/src/algebra_lib/algebra_lib.hpp>

using namespace algebra_lib;

namespace hmc{
    class prior {
    public:
        // Constructors
        prior(sparse_vector &mean, sparse_vector &std);

        // Fields
        const sparse_vector _means;
        const sparse_matrix _covariance;
        const sparse_matrix _inv_cov_m;

        // Member functions
        double misfit(sparse_vector &parameters);
        sparse_vector gradient_misfit(sparse_vector &parameters);
    };
}



#endif //HMC_LINEAR_SYSTEM_PRIOR_HPP
