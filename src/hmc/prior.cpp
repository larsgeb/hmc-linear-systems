//
// Created by Lars Gebraad on 18-8-17.
//

#include "prior.hpp"

using namespace algebra_lib;

namespace hmc {
    prior::prior(vector &mean, vector &stdv) :
            _means(mean),
            _covariance(VectorToDiagonal(ElementWiseMultiplication(stdv, stdv))),
            _inv_cov_m(ElementWiseDivision(1.0, _covariance)) {}

    double prior::misfit(vector &parameters) {
        return 0.5 * ((parameters - _means).TransposeSelf() * _inv_cov_m * (parameters - _means));
    }

    vector prior::gradient_misfit(vector &parameters) {
        return _inv_cov_m * (parameters - _means);
    }
}