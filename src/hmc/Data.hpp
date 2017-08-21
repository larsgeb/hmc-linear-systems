//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_DATA_HPP
#define HMC_LINEAR_SYSTEM_DATA_HPP

#include <algebra_lib/src/algebra_lib/algebra_lib.hpp>
#include "forwardmodel.hpp"

namespace hmc {
    class data {
    public:
        // Constructors
        data(ForwardModel &forward_model, algebra_lib::vector &data, double cov, bool percentage);

        // Member fields
        const unsigned long _numberData;
        const algebra_lib::vector _observedData;
        const algebra_lib::matrix _inv_cov_d;
        const algebra_lib::matrix _G;
        // Precomputed matrices
        const algebra_lib::vector _tG_invCd_d;
        const algebra_lib::matrix _tG_invCd_G;

        // Member functions
        double misfit(algebra_lib::vector &in_parameters);
        algebra_lib::vector gradient_misfit(algebra_lib::vector &parameters);

    private:
        // Member functions
        algebra_lib::matrix calculate_inverse_data_covariance(double std);
        algebra_lib::matrix calculate_inverse_data_covariance_percentual(double percentage);
    };
}

#endif //HMC_LINEAR_SYSTEM_DATA_HPP
