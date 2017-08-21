//
// Created by Lars Gebraad on 18-8-17.
//

#include "posterior.hpp"

using namespace hmc;
using namespace algebra_lib;

double Posterior::misfit(vector &parameters, prior &prior, data &data) {
    return prior.misfit(parameters) + data.misfit(parameters);
}

vector Posterior::gradient_misfit(vector &parameters, prior &prior, data &data) {
    return prior.gradient_misfit(parameters) + data.gradient_misfit(parameters);
}
