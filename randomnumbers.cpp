//
// Created by Lars Gebraad on 7/10/17.
//
#include <algebra_lib/src/algebra_lib/algebra_lib.hpp>
#include "randomnumbers.hpp"
#include <cmath>

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


algebra_lib::vector randn(algebra_lib::vector means, algebra_lib::vector cov) {
    return means + randn(std::move(cov));
}

algebra_lib::vector randn(algebra_lib::vector cov) {
    // Zero mean
    algebra_lib::vector samples(cov.size(), true);
    for (int iStd = 0; iStd < cov.size(); ++iStd) {
        samples[iStd] = (randn(0, sqrt(cov[iStd])));
    }
    return samples;
}

algebra_lib::vector
randn_Cholesky(algebra_lib::vector mean, algebra_lib::matrix CholeskyLower_CovarianceMatrix) {
    return mean + randn_Cholesky(std::move(CholeskyLower_CovarianceMatrix));
}

algebra_lib::vector randn_Cholesky(algebra_lib::matrix CholeskyLower_CovarianceMatrix) {
    // Assumes zero mean
    algebra_lib::vector uncorrelated(CholeskyLower_CovarianceMatrix.rows(), true);

    for (int i = 0; i < CholeskyLower_CovarianceMatrix.rows(); ++i) {
        uncorrelated[i] = randn(0, 1);
    }

    return CholeskyLower_CovarianceMatrix * uncorrelated;

}

algebra_lib::vector randn(algebra_lib::matrix DiagonalCovarianceMatrix) {
    algebra_lib::vector samples(DiagonalCovarianceMatrix.rows(), true);

    for (int i = 0; i < DiagonalCovarianceMatrix.rows(); i++) {
        samples[i] = randn(0, sqrt(DiagonalCovarianceMatrix[i][i]));
    }
    return samples;
}
