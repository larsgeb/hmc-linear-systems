//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_POSTERIOR_HPP
#define HMC_LINEAR_SYSTEM_POSTERIOR_HPP

#include <AlgebraLib/src/algebra_lib/algebra_lib.hpp>
#include "prior.hpp"
#include "data.hpp"

using namespace algebra_lib;

namespace hmc {
    class Posterior {
        double misfit(vector &parameters, prior &prior, data &data);

        vector gradient_misfit(vector &parameters, prior &prior, data &data);
    };
}

#endif //HMC_LINEAR_SYSTEM_POSTERIOR_HPP
