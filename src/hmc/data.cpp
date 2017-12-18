//
// Created by Lars Gebraad on 18-8-17.
//

#include "data.hpp"

namespace hmc {
    data::data(hmc::forward_model &forward_model, char *input_data_file, double cov, bool percentage) {
        // Constructor with fixed covariance (percentual or absolute), implies NO correlation
        _observedData.load(input_data_file);
        _numberData = _observedData.n_elem;
        _inv_cov_d = percentage ?
                     calculate_inverse_data_covariance_percentual(cov) :
                     calculate_inverse_data_covariance(cov);
    };

    data::data(hmc::forward_model &forward_model, char *input_data_file, char *cov_matrix_file) {
        // Constructor with supplied matrix
        _observedData.load(input_data_file);
        _numberData = _observedData.n_elem;
        arma::mat cov_matrix;
        cov_matrix.load(cov_matrix_file);
        _inv_cov_d = arma::inv_sympd(cov_matrix);
    };

    arma::mat data::calculate_inverse_data_covariance(double std) {
        // Covariance is interpreted as absolute value of all datapoints.
        return (1.0 / (std * std)) * arma::diagmat(arma::ones<arma::vec>(_numberData));
    }

    arma::mat data::calculate_inverse_data_covariance_percentual(double percentage) {
        // Covariance is interpreted as percentual value of measured value.
        arma::mat invStd = arma::diagmat(1.0 / arma::square(_observedData * (percentage / 100.0)));
        return conv_to<arma::mat>::from(invStd);
    }

    data::data() = default;
}
