//
// Created by Lars Gebraad on 18-8-17.
//

#include <armadillo>
#include <utility>
#include "prior.hpp"

namespace hmc {
    prior::prior(arma::vec mean, arma::vec stdv) :
            _means(std::move(mean)),
            _covariance(arma::diagmat(arma::square(stdv))),
            _inv_cov_m(arma::diagmat(1.0 / arma::square(stdv))) {}

    double prior::misfit(arma::vec &parameters) {
        return arma::as_scalar(0.5 * ((parameters - _means).t() * _inv_cov_m * (parameters - _means)));
    }

    arma::vec prior::gradient_misfit(arma::vec &parameters) {
        return _inv_cov_m * (parameters - _means);
    }

    prior::prior() = default;
}