//
// Created by Lars Gebraad on 18-8-17.
//

#include "posterior.hpp"

using namespace hmc;
using namespace algebra_lib;

double Posterior::misfit(arma::vec &parameters, prior &prior, data &data) {
    return prior.misfit(parameters) + data.misfit(parameters);
}

arma::vec Posterior::gradient_misfit(arma::vec &parameters, prior &prior, data &data) {
    return prior.gradient_misfit(parameters) + data.gradient_misfit(parameters);
}
