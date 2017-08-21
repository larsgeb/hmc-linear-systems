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
        prior(vector &mean, vector &std);

        // Fields
        const vector _means;
        const matrix _covariance;
        const matrix _inv_cov_m;

        // Member functions
        double misfit(vector &parameters);
        vector gradient_misfit(vector &parameters);
    };
}



#endif //HMC_LINEAR_SYSTEM_PRIOR_HPP
