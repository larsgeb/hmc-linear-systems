//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_PRIOR_HPP
#define HMC_LINEAR_SYSTEM_PRIOR_HPP

#include <armadillo>

namespace hmc{
    class prior {
    public:
        // Constructors
        prior(char * mean_vector_file, char * cov_matrix_file);
        prior(double mean, double std, arma::uword dim);
        prior();


        // Fields
        arma::vec _means;
        arma::mat _covariance;
        arma::mat _inv_cov_m;

    };
}



#endif //HMC_LINEAR_SYSTEM_PRIOR_HPP
