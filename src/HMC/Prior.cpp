//
// Created by Lars Gebraad on 18-8-17.
//

#include "Prior.hpp"

namespace HMC {
    Prior::Prior(AlgebraLib::Vector &mean, AlgebraLib::Vector &stdv) :
            _means(mean),
            _covariance(AlgebraLib::VectorToDiagonal(AlgebraLib::ElementWiseMultiplication(stdv, stdv))),
            _inv_cov_m(AlgebraLib::ElementWiseDivision(1.0, _covariance)) {}

    double Prior::misfit(AlgebraLib::Vector &parameters) {
        return 0.5 * ((parameters - _means).TransposeSelf() * _inv_cov_m * (parameters - _means));
    }

    AlgebraLib::Vector Prior::gradient_misfit(AlgebraLib::Vector &parameters) {
        return _inv_cov_m * (parameters - _means);
    }
}