//
// Created by Lars Gebraad on 18-8-17.
//

#include "Data.hpp"

namespace HMC {
    // TODO implement constructors with predefined data covariance matrices.
    Data::Data(ForwardModel &forward_model, AlgebraLib::Vector &data, double cov, bool percentage)
            : _numberData(data.size()),
              _observedData(data),
              _inv_cov_d(percentage ?
                         calculate_inverse_data_covariance_percentual(cov) :
                         calculate_inverse_data_covariance(cov)),
              _G(forward_model._g),
              _tG_invCd_d(forward_model._g.Transpose() * _inv_cov_d * _observedData),
              _tG_invCd_G(forward_model._g.Transpose() * _inv_cov_d * forward_model._g) {}

    AlgebraLib::Matrix Data::calculate_inverse_data_covariance(double std) {
        // Covariance is interpreted as absolute value of all datapoints.
        return (1.0 / (std * std)) * (AlgebraLib::Matrix(_numberData, _numberData).Unit());
    }

    AlgebraLib::Matrix Data::calculate_inverse_data_covariance_percentual(double percentage) {
        // Covariance is interpreted as percentual value of measured value.
        AlgebraLib::Vector invStd = AlgebraLib::ElementWiseDivision(1.0, (percentage / 100.0) * _observedData);

        return AlgebraLib::VectorToDiagonal(AlgebraLib::ElementWiseMultiplication(invStd, invStd));

    }

    double Data::misfit(AlgebraLib::Vector &in_parameters) {
        return 0.5 * ((_G * in_parameters - _observedData).TransposeSelf()
                      * _inv_cov_d * (_G * in_parameters - _observedData)
        );
    }

    AlgebraLib::Vector Data::gradient_misfit(AlgebraLib::Vector &parameters) {
        return _tG_invCd_G * parameters - _tG_invCd_d;
    }
}