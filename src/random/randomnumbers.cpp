//
// Created by Lars Gebraad on 7/10/17.
//
#include "randomnumbers.hpp"
#include <cmath>
#include <random>
#include <armadillo/armadillo-8.200.2/include/armadillo>

// todo rewrite using C++11 random functions.
// Random number generators
/* Uniformly distributed, double-valued random numbers. ---------------------------*/
double randf(double min, double max) {
    return (max - min) * (double) rand() / RAND_MAX + min;
}

double randn(double mean, double stdv) {
    double x;

    double z1 = (double) rand() / RAND_MAX;
    double z2 = (double) rand() / RAND_MAX;

    x = sqrt(-2.0 * log(z1)) * cos(2.0 * PI * z2);
    x = stdv * x + mean;

    return x;

}

arma::vec randn(arma::vec means, arma::vec cov) {
    return means + randn(std::move(cov));
}

arma::vec randn(arma::vec cov) {
    // Zero mean
    arma::vec samples(cov.size(), true);
    for (int iStd = 0; iStd < cov.size(); ++iStd) {
        samples[iStd] = (randn(0, sqrt(cov[iStd])));
    }
    return samples;
}

arma::vec randn_Cholesky(arma::vec mean, arma::mat CholeskyLower_CovarianceMatrix) {
    return mean + randn_Cholesky(std::move(CholeskyLower_CovarianceMatrix));
}

arma::vec randn_Cholesky(arma::mat CholeskyLower_CovarianceMatrix) {
    // Assumes zero mean
    arma::vec uncorrelated(CholeskyLower_CovarianceMatrix.n_rows, true);

    for (int i = 0; i < CholeskyLower_CovarianceMatrix.n_rows; ++i) {
        uncorrelated[i] = randn(0, 1);
    }

    return CholeskyLower_CovarianceMatrix * uncorrelated;

}

arma::vec randn(arma::mat DiagonalCovarianceMatrix) {
    // Generate uncorrelated samples from diagonal.
    arma::vec samples(DiagonalCovarianceMatrix.n_rows, true);

    for (int i = 0; i < DiagonalCovarianceMatrix.n_rows; i++) {
        samples[i] = randn(0, sqrt(DiagonalCovarianceMatrix(i, i)));
    }
    return samples;
}


