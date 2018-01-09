//
// Created by Lars Gebraad on 18-8-17.
//

#include <armadillo/armadillo-8.200.2/include/armadillo>
#include <utility>
#include "prior.hpp"

namespace hmc {
    prior::prior(double mean, double std, arma::uword dim) {
        _means = mean * arma::ones<arma::vec>(dim);
        _covariance =  (std * std) * arma::diagmat(arma::ones<arma::vec>(dim));
        _inv_cov_m = (1.0 / (std * std)) * arma::diagmat(arma::ones<arma::vec>(dim));
    }

    prior::prior(char *mean_vector_file, char *cov_matrix_file) {
        _means.load(mean_vector_file);
        _covariance.load(cov_matrix_file);
        _inv_cov_m = arma::inv_sympd(_covariance);
    }

    prior::prior() = default;
}