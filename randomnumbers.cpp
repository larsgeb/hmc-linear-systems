//
// Created by Lars Gebraad on 7/10/17.
//
#include <SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp>
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


AlgebraLib::Vector randn(AlgebraLib::Vector means, AlgebraLib::Vector stdv) {
    return means + randn(std::move(stdv));
}

AlgebraLib::Vector randn(AlgebraLib::Vector stdv) {
    // Zero mean
    AlgebraLib::Vector samples(stdv.size(), true);
    for (int iStd = 0; iStd < stdv.size(); ++iStd) {
        samples[iStd] = (randn(0, stdv[iStd]));
    }
    return samples;
}

AlgebraLib::Vector
randn_Cholesky(AlgebraLib::Vector mean, AlgebraLib::Matrix CholeskyLower_CovarianceMatrix) {
    return mean + randn_Cholesky(std::move(CholeskyLower_CovarianceMatrix));
}

AlgebraLib::Vector randn_Cholesky(AlgebraLib::Matrix CholeskyLower_CovarianceMatrix) {
    // Assumes zero mean
    AlgebraLib::Vector uncorrelated(CholeskyLower_CovarianceMatrix.rows(), true);

    for (int i = 0; i < CholeskyLower_CovarianceMatrix.rows(); ++i) {
        uncorrelated[i] = randn(0, 1);
    }

    return CholeskyLower_CovarianceMatrix * uncorrelated;

}

AlgebraLib::Vector randn(AlgebraLib::Matrix DiagonalCovarianceMatrix) {
    AlgebraLib::Vector samples(DiagonalCovarianceMatrix.rows(), true);

    for (int i = 0; i < DiagonalCovarianceMatrix.rows(); i++) {
        samples[i] = randn(0, sqrt(DiagonalCovarianceMatrix[i][i]));
    }
    return samples;
}
