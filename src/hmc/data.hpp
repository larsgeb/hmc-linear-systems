//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_DATA_HPP
#define HMC_LINEAR_SYSTEM_DATA_HPP

#include <AlgebraLib/src/algebra_lib/algebra_lib.hpp>
#include "forward_model.hpp"

using namespace algebra_lib;

namespace hmc {
    class data {
    public:
        // Constructors
        data(forward_model &forward_model, sparse_vector &data, double cov, bool percentage);

        // Member fields
        const unsigned int _numberData;
        const sparse_vector _observedData;
        const sparse_matrix _inv_cov_d;
        const sparse_matrix _G;
        // Precomputed matrices
        const sparse_vector _tG_invCd_d;
        const sparse_matrix _tG_invCd_G;

        // Member functions
        double misfit(sparse_vector &in_parameters);

        sparse_vector gradient_misfit(sparse_vector &parameters);

    private:
        // Member functions
        sparse_matrix calculate_inverse_data_covariance(double std);

        sparse_matrix calculate_inverse_data_covariance_percentual(double percentage);
    };
}

#endif //HMC_LINEAR_SYSTEM_DATA_HPP
