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

vector randn(vector means, vector cov) {
    return means + randn(std::move(cov));
}

vector randn(vector cov) {
    // Zero mean
    vector samples(cov.size(), true);
    for (int iStd = 0; iStd < cov.size(); ++iStd) {
        samples[iStd] = (randn(0, sqrt(cov[iStd])));
    }
    return samples;
}

vector
randn_Cholesky(vector mean, matrix CholeskyLower_CovarianceMatrix) {
    return mean + randn_Cholesky(std::move(CholeskyLower_CovarianceMatrix));
}

vector randn_Cholesky(matrix CholeskyLower_CovarianceMatrix) {
    // Assumes zero mean
    vector uncorrelated(CholeskyLower_CovarianceMatrix.rows(), true);

    for (int i = 0; i < CholeskyLower_CovarianceMatrix.rows(); ++i) {
        uncorrelated[i] = randn(0, 1);
    }

    return CholeskyLower_CovarianceMatrix * uncorrelated;

}

vector randn(matrix DiagonalCovarianceMatrix) {
    vector samples(DiagonalCovarianceMatrix.rows(), true);

    for (int i = 0; i < DiagonalCovarianceMatrix.rows(); i++) {
        samples[i] = randn(0, sqrt(DiagonalCovarianceMatrix[i][i]));
    }
    return samples;
}

/*
// using C++11 rand library
double randrandf(double min, double max) {

    return std::uniform_real_distribution<double>(min, max).operator();
}

// using C++11 rand library
double randrandn(double mean, double stdv) {
    double x;

    double z1 = (double) rand() / RAND_MAX;
    double z2 = (double) rand() / RAND_MAX;

    x = sqrt(-2.0 * log(z1)) * cos(2.0 * PI * z2);
    x = stdv * x + mean;

    return x;

}*/

