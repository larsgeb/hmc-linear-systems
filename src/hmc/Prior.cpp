//
// Created by Lars Gebraad on 18-8-17.
//

#include "prior.hpp"

namespace hmc {
    prior::prior(algebra_lib::vector &mean, algebra_lib::vector &stdv) :
            _means(mean),
            _covariance(algebra_lib::VectorToDiagonal(algebra_lib::ElementWiseMultiplication(stdv, stdv))),
            _inv_cov_m(algebra_lib::ElementWiseDivision(1.0, _covariance)) {}

    double prior::misfit(algebra_lib::vector &parameters) {
        return 0.5 * ((parameters - _means).TransposeSelf() * _inv_cov_m * (parameters - _means));
    }

    algebra_lib::vector prior::gradient_misfit(algebra_lib::vector &parameters) {
        return _inv_cov_m * (parameters - _means);
    }
}