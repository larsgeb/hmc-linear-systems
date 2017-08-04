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

double randn(double mean, double stdv) {
    double x;

    double z1 = (double) rand() / RAND_MAX;
    double z2 = (double) rand() / RAND_MAX;

    x = sqrt(-2.0 * log(z1)) * cos(2.0 * PI * z2);
    x = stdv * x + mean;

    return x;
}


std::vector<double> randn(std::vector<double> means, std::vector<double> std) {
    std::vector<double> samples;
    samples.reserve(means.size());
    for (int i = 0; i < means.size(); i++){
        samples.push_back(randn(means[i],std[i]));
    }
    return samples;
}

std::vector<double> randn(std::vector<double> mean, std::vector<std::vector<double>> CholeskyLower_CovarianceMatrix) {
    return mean + randn(std::move(CholeskyLower_CovarianceMatrix));
}

std::vector<double> randn(std::vector<std::vector<double>> CholeskyLower_CovarianceMatrix) {
    // Assumes zero mean
    std::vector<double> uncorrelated;
    uncorrelated.reserve(CholeskyLower_CovarianceMatrix.size());
    for (int i = 0; i < CholeskyLower_CovarianceMatrix.size(); ++i) {
        uncorrelated.push_back(randn(0, 1));
    }
    return std::move(CholeskyLower_CovarianceMatrix) * uncorrelated;

}
