//
// Created by Lars Gebraad on 7/10/17.
//

/*! @file
 * @brief Set of functions to draw from (multivariate, correlated) normal distributions.
 *
 * This set of functions allows one to sample from mutliple types of normal distributions. All are based on a uniform
 * number generator which transforms to normally distributed samples using the Box-Müller transform.
 *
 */

#ifndef HMC_VSP_RANDOMNUMBERS_HPP
#define HMC_VSP_RANDOMNUMBERS_HPP

#include <vector>

const double PI = 3.14159265358979323846264338327;
/*!
 * @brief Draws from Gaussian \f$ \mathcal{N} (\mu,\sigma) \f$ (mean, standard deviation) using Box-Müller transform.
 * @param mean double containing \f$ \mu \f$
 * @param stdv double containing \f$ \sigma \f$
 * @return double, sample from the distribution
 */
double randn(double mean, double stdv);

std::vector<double> randn(std::vector<double> means, std::vector<double> stdv);

std::vector<double> randn(std::vector<double> stdv);
/**
 * @brief Drawing non-zero mean samples from a correlated multivariate Gaussian.
 * Invokes randn_Cholesky(std::vector std::vector double CholeskyLower_CovarianceMatrix) and adds mean.
 * @param mean std::vector of doubles containging n means.
 * @param CholeskyLower_CovarianceMatrix std::vector of std::vector of doubles, containing the n x n Lower Cholesky matrix
 * of the n x n covariance matrix, must be square and lower triangular.
 * @return std::vector of doubles containing the non-zero mean correlated samples.
 */
std::vector<double> randn_Cholesky(std::vector<double> mean, std::vector<std::vector<double>>
CholeskyLower_CovarianceMatrix);

/**
 * @brief Drawing zero mean samples from a correlated multivariate Gaussian. This algorithm uses the lower Cholesky matrix
 * of the (positive definite ermitian) covariance matrix to transform (affine transform) n uncorrelated samples with standard
 * deviation 1 to the right covariances. The mean is assumed zero.
 * @param CholeskyLower_CovarianceMatrix std::vector of std::vector of doubles, containing the n x n Lower Cholesky matrix
 * of the n x n covariance matrix, must be square and lower triangular.
 * @return std::vector of doubles containing the zero mean correlated samples.
 */
std::vector<double> randn_Cholesky(std::vector<std::vector<double>> CholeskyLower_CovarianceMatrix);

/**
 * @brief Drawing n zero mean samples from an n x n diagonal variance matrix. No correlation between the parameters.
 * @param DiagonalCovarianceMatrix std::vector of a std::vector of doubles, containing on the diagonal the variance, or
 * standard deviation squared.
 * @return std::vector of doubles containing n samples.
 */
std::vector<double> randn(std::vector<std::vector<double>> DiagonalCovarianceMatrix);

double randf(double min, double max);

#endif //HMC_VSP_RANDOMNUMBERS_HPP
