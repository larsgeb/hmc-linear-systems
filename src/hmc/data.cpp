//
// Created by Lars Gebraad on 18-8-17.
//

#include "data.hpp"

namespace hmc {
    // TODO implement constructors with predefined data covariance matrices.
    data::data(forward_model &forward_model, arma::vec &data, double cov, bool percentage)
            : _numberData(data.n_elem),
              _observedData(data),
              _inv_cov_d(percentage ?
                         calculate_inverse_data_covariance_percentual(cov) :
                         calculate_inverse_data_covariance(cov)),
              _G(forward_model._g),
              _tG_invCd_d(forward_model._g.t() * _inv_cov_d * _observedData),
              _tG_invCd_G((forward_model._g.t() * _inv_cov_d) * forward_model._g) {}

    arma::mat data::calculate_inverse_data_covariance(double std) {
        // Covariance is interpreted as absolute value of all datapoints.
        return (1.0 / (std * std)) * ones<mat>(_numberData);
    }

    arma::mat data::calculate_inverse_data_covariance_percentual(double percentage) {
        // Covariance is interpreted as percentual value of measured value.
        arma::mat invStd = arma::diagmat(1.0 / arma::square(_observedData * (percentage / 100.0)));
        return conv_to<arma::mat>::from(invStd);
    }

    double data::misfit(arma::vec &in_parameters) {
        double a = conv_to<double>::from(0.5 * _G * (in_parameters - _observedData).t() * (_inv_cov_d) * _G *
                                         (in_parameters - _observedData));
        return a;
    }

    arma::vec data::gradient_misfit(arma::vec &parameters) {
        return _tG_invCd_G * parameters - _tG_invCd_d;
    }
}