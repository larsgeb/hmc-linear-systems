//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_POSTERIOR_HPP
#define HMC_LINEAR_SYSTEM_POSTERIOR_HPP

#include <algebra_lib/src/algebra_lib/algebra_lib.hpp>
#include "prior.hpp"
#include "data.hpp"

namespace hmc {
    class Posterior {
        double misfit(algebra_lib::vector &parameters, prior &prior, data &data);

        algebra_lib::vector gradient_misfit(algebra_lib::vector &parameters, prior &prior, data &data);
    };
}

#endif //HMC_LINEAR_SYSTEM_POSTERIOR_HPP
