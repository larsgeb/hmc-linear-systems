//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_PRIOR_HPP
#define HMC_LINEAR_SYSTEM_PRIOR_HPP

#include <algebra_lib/src/algebra_lib/algebra_lib.hpp>

namespace hmc{
    class prior {
    public:
        // Constructors
        prior(algebra_lib::vector &mean, algebra_lib::vector &std);

        // Fields
        const algebra_lib::vector _means;
        const algebra_lib::matrix _covariance;
        const algebra_lib::matrix _inv_cov_m;

        // Member functions
        double misfit(algebra_lib::vector &parameters);
        algebra_lib::vector gradient_misfit(algebra_lib::vector &parameters);
    };
}



#endif //HMC_LINEAR_SYSTEM_PRIOR_HPP
