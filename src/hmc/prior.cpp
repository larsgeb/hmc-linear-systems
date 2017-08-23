//
// Created by Lars Gebraad on 18-8-17.
//

#include "prior.hpp"

using namespace algebra_lib;

namespace hmc {
    prior::prior(sparse_vector &mean, sparse_vector &stdv) :
            _means(mean),
            _covariance(VectorToDiagonal(ElementWiseMultiplication(stdv, stdv))),
            _inv_cov_m(ElementWiseDivision(1.0, _covariance)) {}

    double prior::misfit(sparse_vector &parameters) {
        return 0.5 * ((parameters - _means).TransposeSelf() * _inv_cov_m * (parameters - _means));
    }

    sparse_vector prior::gradient_misfit(sparse_vector &parameters) {
        return _inv_cov_m * (parameters - _means);
    }
}