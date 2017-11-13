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
#include <armadillo>

const double PI = 3.14159265358979323846264338327;

/*!
 * @brief Draws from Gaussian \f$ \mathcal{N} (\mu,\sigma) \f$ (mean, standard deviation) using Box-Müller transform.
 * @param mean double containing \f$ \mu \f$
 * @param stdv double containing \f$ \sigma \f$
 * @return double, sample from the distribution
 */
double randn(double mean, double stdv);

/*!
 * @brief Draws from uncorrelated Gaussians \f$ \mathcal{N} (\boldsymbol \mu,\boldsymbol{\sigma}) \f$ (vectors of mean, standard deviation) using Box-Müller transform. Loops over both vectors and calls randn(double mean, double stdv) every
 * iteration.
 * @param mean vector containing \f$ \mu_i \f$
 * @param cov vector containing \f$ \sigma_i \f$
 * @return Vector of samples from the distributions.
 */
arma::vec randn(arma::vec means, arma::vec cov);

/*!
 * @brief Draws zero-mean samples from uncorrelated Gaussians \f$ \mathcal{N} (\boldsymbol 0,\boldsymbol{\sigma}) \f$
 * (standard deviation) using Box-Müller transform. Loops over both vectors and calls randn(double mean, double stdv) every
 * iteration.
 * @param cov vector containing \f$ \sigma_i \f$
 * @return Vector of samples from the distributions.
 */
arma::vec randn(arma::vec cov);

/**
 * @brief Drawing non-zero mean samples from an \f$ n \f$ dimensional correlated Gaussian.
 * Invokes randn_Cholesky(std::vector<std::vector<double>> CholeskyLower_CovarianceMatrix) and adds mean.
 * @param mean vector containing \f$ n \f$ means.
 * @param CholeskyLower_CovarianceMatrix Matrix containing the n x n Lower Cholesky matrix
 * of the n x n covariance matrix \f$ \boldsymbol \Sigma \f$, must be square and lower triangular.
 * @return Vector containing the non-zero mean correlated samples.
 */
arma::vec randn_Cholesky(arma::vec mean, arma::mat CholeskyLower_CovarianceMatrix);

/**
 * @brief Drawing non-zero mean samples from an \f$ n \f$ dimensional correlated Gaussian. This algorithm uses the lower
 * Cholesky matrix of the (positive definite Hermitian) covariance matrix to transform (affine transform) n uncorrelated
 * samples with standard deviation 1 to the right covariances. The mean is assumed zero.
 * @param CholeskyLower_CovarianceMatrix Matrix containing the n x n Lower Cholesky matrix
 * of the n x n covariance matrix \f$ \boldsymbol \Sigma \f$, must be square and lower triangular.
 * @return Vector containing the zero mean correlated samples.
 */
arma::vec randn_Cholesky(arma::mat CholeskyLower_CovarianceMatrix);

/**
 * @brief Drawing n zero mean samples from \f$ \mathcal{N} (\boldsymbol 0,\boldsymbol{\sigma}) \f$. No correlation is present between the parameters.
 * @param DiagonalCovarianceMatrix Matrix containing on the diagonal the variance, or standard deviation squared.
 * @return Vector containing n samples.
 */
arma::vec randn(arma::mat DiagonalCovarianceMatrix);

/**
 * @brief Draw uniformly distributed samples between two numbers.
 * @param min Minimum of the distribution.
 * @param max Maximum of the distribution.
 * @return Sample.
 */
double randf(double min, double max);

#endif //HMC_VSP_RANDOMNUMBERS_HPP
