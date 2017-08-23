//
// Created by Lars Gebraad on 18-8-17.
//

#include "data.hpp"

namespace hmc {
    // TODO implement constructors with predefined data covariance matrices.
    data::data(forward_model &forward_model, algebra_lib::sparse_vector &data, double cov, bool percentage)
            : _numberData(data.size()),
              _observedData(data),
              _inv_cov_d(percentage ?
                         calculate_inverse_data_covariance_percentual(cov) :
                         calculate_inverse_data_covariance(cov)),
              _G(forward_model._g),
              _tG_invCd_d(forward_model._g.Transpose() * _inv_cov_d * _observedData),
              _tG_invCd_G(forward_model._g.Transpose() * _inv_cov_d * forward_model._g) {}

    algebra_lib::sparse_matrix data::calculate_inverse_data_covariance(double std) {
        // Covariance is interpreted as absolute value of all datapoints.
        return (1.0 / (std * std)) * (algebra_lib::sparse_matrix(_numberData, _numberData).Unit());
    }

    algebra_lib::sparse_matrix data::calculate_inverse_data_covariance_percentual(double percentage) {
        // Covariance is interpreted as percentual value of measured value.
        algebra_lib::sparse_vector invStd = algebra_lib::ElementWiseDivision(1.0, (percentage / 100.0) * _observedData);

        return algebra_lib::VectorToDiagonal(algebra_lib::ElementWiseMultiplication(invStd, invStd));

    }

    double data::misfit(algebra_lib::sparse_vector &in_parameters) {
        return 0.5 * ((_G * in_parameters - _observedData).TransposeSelf()
                      * _inv_cov_d * (_G * in_parameters - _observedData)
        );
    }

    algebra_lib::sparse_vector data::gradient_misfit(algebra_lib::sparse_vector &parameters) {
        return _tG_invCd_G * parameters - _tG_invCd_d;
    }
}