//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_DATA_HPP
#define HMC_LINEAR_SYSTEM_DATA_HPP

#include <AlgebraLib/src/algebra_lib/algebra_lib.hpp>
#include "forward_model.hpp"

namespace hmc {
    class data {
    public:
        // Constructors
        data(forward_model &forward_model, arma::vec &data, double cov, bool percentage);

        // Member fields
        const unsigned long _numberData;
        const arma::vec _observedData;
        const arma::mat _inv_cov_d;
        const arma::mat _G;
        // Precomputed matrices
        const arma::vec _tG_invCd_d;
        const arma::mat _tG_invCd_G;

        // Member functions
        double misfit(arma::vec &in_parameters);

        arma::vec gradient_misfit(arma::vec &parameters);

    private:
        // Member functions
        arma::mat calculate_inverse_data_covariance(double std);

        arma::mat calculate_inverse_data_covariance_percentual(double percentage);
    };
}

#endif //HMC_LINEAR_SYSTEM_DATA_HPP
