//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_DATA_HPP
#define HMC_LINEAR_SYSTEM_DATA_HPP

#include <SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp>
#include "ForwardModel.hpp"

namespace HMC {
    class Data {
    public:
        // Constructors
        Data(ForwardModel &forward_model, AlgebraLib::Vector &data, double cov, bool percentage);

        // Member fields
        const unsigned long _numberData;
        const AlgebraLib::Vector _observedData;
        const AlgebraLib::Matrix _inv_cov_d;
        const AlgebraLib::Matrix _G;
        // Precomputed matrices
        const AlgebraLib::Vector _tG_invCd_d;
        const AlgebraLib::Matrix _tG_invCd_G;

        // Member functions
        double misfit(AlgebraLib::Vector &in_parameters);
        AlgebraLib::Vector gradient_misfit(AlgebraLib::Vector &parameters);

    private:
        // Member functions
        AlgebraLib::Matrix calculate_inverse_data_covariance(double std);
        AlgebraLib::Matrix calculate_inverse_data_covariance_percentual(double percentage);
    };
}

#endif //HMC_LINEAR_SYSTEM_DATA_HPP
