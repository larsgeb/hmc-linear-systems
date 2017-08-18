//
// Created by Lars Gebraad on 18-8-17.
//

#include "Posterior.hpp"

using namespace HMC;
using namespace AlgebraLib;

double Posterior::misfit(Vector &parameters, Prior &prior, Data &data) {
    return prior.misfit(parameters) + data.misfit(parameters);
}

Vector Posterior::gradient_misfit(Vector &parameters, Prior &prior, Data &data) {
    return prior.gradient_misfit(parameters) + data.gradient_misfit(parameters);
}
