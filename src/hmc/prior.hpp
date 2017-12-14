//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_PRIOR_HPP
#define HMC_LINEAR_SYSTEM_PRIOR_HPP

#include <armadillo/armadillo-8.200.2/include/armadillo>

namespace hmc{
    class prior {
    public:
        // Constructors
        prior(arma::vec mean, arma::vec std);

        prior();


        // Fields
        arma::vec _means;
        arma::mat _covariance;
        arma::mat _inv_cov_m;

        // Member functions
        double misfit(arma::vec &parameters);
        arma::vec gradient_misfit(arma::vec &parameters);
    };
}



#endif //HMC_LINEAR_SYSTEM_PRIOR_HPP
