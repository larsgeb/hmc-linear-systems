//
// Created by Lars Gebraad on 7/10/17.
//
#include <cstdlib>
#include "randomnumbers.hpp"
#include "linearalgebra.hpp"
#include <cmath>
#include <utility>

// Random number generators
/* Uniformly distributed, double-valued random numbers. ---------------------------*/
double randf(double min, double max) {
    return (max - min) * (double) rand() / RAND_MAX + min;
}

/* Normally distributed, double-valued random numbers. ----------------------------*/
/* This function implements the Box-Muller transform to obtain a pair of
 normally distributed random numbers with a given mean and standard deviation. */
void randn(double mean, double stdv, double *x1, double *x2) {
    double z1 = (double) rand() / RAND_MAX;
    double z2 = (double) rand() / RAND_MAX;

    *x1 = sqrt(-2.0 * log(z1)) * cos(2.0 * PI * z2);
    *x2 = sqrt(-2.0 * log(z1)) * sin(2.0 * PI * z2);

    *x1 = stdv * (*x1) + mean;
    *x2 = stdv * (*x2) + mean;
}

double randn(double mean, double stdv) {
    double x;

    double z1 = (double) rand() / RAND_MAX;
    double z2 = (double) rand() / RAND_MAX;

    x = sqrt(-2.0 * log(z1)) * cos(2.0 * PI * z2);
    x = stdv * x + mean;

    return x;
}

std::vector<double> randn(std::vector<double> mean, std::vector<std::vector<double>> CholeskyLower_CovarianceMatrix) {

    std::vector<double> sample;
    std::vector<double> uncorrelated;
    uncorrelated.reserve(mean.size());

    for (int i = 0; i < mean.size(); ++i) {
        uncorrelated.push_back(randn(0, 1));
    }
    sample = mean + std::move(CholeskyLower_CovarianceMatrix) * uncorrelated;

    double a = sample[0];
    return sample;
}
