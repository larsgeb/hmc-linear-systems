//
// Created by Lars Gebraad on 18-8-17.
//

#include "posterior.hpp"

using namespace hmc;
using namespace algebra_lib;

double posterior::misfit(sparse_vector &parameters, prior &prior, data &data) {
    return prior.misfit(parameters) + data.misfit(parameters);
}

sparse_vector posterior::gradient_misfit(sparse_vector &parameters, prior &prior, data &data) {
    return prior.gradient_misfit(parameters) + data.gradient_misfit(parameters);
}
