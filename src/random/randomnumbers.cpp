//
// Created by Lars Gebraad on 7/10/17.
//
#include <AlgebraLib/src/algebra_lib/algebra_lib.hpp>
#include "randomnumbers.hpp"
#include <cmath>
#include <random>

using namespace algebra_lib;

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

vector randn(const vector &means, const vector &cov) {
    return means + randn(cov);
}

sparse_vector randn(const sparse_vector &means, const sparse_vector &cov) {
    return means + randn(cov);
}

vector randn(const vector &cov) {
    vector samples(cov.size(), true);
    for (int iStd = 0; iStd < cov.size(); ++iStd) {
        samples[iStd] = (randn(0, sqrt(cov[iStd])));
    }
    return samples;
}

sparse_vector randn(const sparse_vector &cov) {
    sparse_vector samples(cov.size(), cov.isColumn());
    for (auto &&row : cov) {
        samples(row.first) = (randn(0, sqrt(cov[row.first])));
    }
    return samples;
}

vector randn_Cholesky(const vector &mean, const matrix &CholeskyLower_CovarianceMatrix) {
    return mean + randn_Cholesky(CholeskyLower_CovarianceMatrix);
}

sparse_vector randn_Cholesky(const sparse_vector &mean, const sparse_matrix &CholeskyLower_CovarianceMatrix) {
    return mean + randn_Cholesky(CholeskyLower_CovarianceMatrix);
}

vector randn_Cholesky(const matrix &CholeskyLower_CovarianceMatrix) {
    // Assumes zero mean
    vector uncorrelated(CholeskyLower_CovarianceMatrix.rows(), true);

    for (int i = 0; i < CholeskyLower_CovarianceMatrix.rows(); ++i) {
        uncorrelated[i] = randn(0, 1);
    }

    return CholeskyLower_CovarianceMatrix * uncorrelated;

}

sparse_vector randn_Cholesky(const sparse_matrix &CholeskyLower_CovarianceMatrix) {
    // Assumes zero mean
    sparse_vector uncorrelated(CholeskyLower_CovarianceMatrix.rows(), true);

    for (int i = 0; i < CholeskyLower_CovarianceMatrix.rows(); ++i) {
        uncorrelated(i) = randn(0, 1);
    }

    return CholeskyLower_CovarianceMatrix * uncorrelated;
}

vector randn(const matrix &DiagonalCovarianceMatrix) {
    vector samples(DiagonalCovarianceMatrix.rows(), true);

    for (int i = 0; i < DiagonalCovarianceMatrix.rows(); i++) {
        samples[i] = randn(0, sqrt(DiagonalCovarianceMatrix[i][i]));
    }
    return samples;
}

sparse_vector randn(const sparse_matrix &DiagonalCovarianceMatrix) {
    sparse_vector samples(DiagonalCovarianceMatrix.rows(), true);

    for (auto i = 0; i < DiagonalCovarianceMatrix.rows(); i++) {
        samples(i) = randn(0, sqrt(DiagonalCovarianceMatrix[i][i]));
    }
    return samples;
}

