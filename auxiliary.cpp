// Created by Lars Gebraad on 7/10/17.

#include <utility>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "auxiliary.hpp"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"

/* ----------------------------------------------------------------------------------------------------------------------- *
 * Class for Gaussian Distributed prior information about model parameters for a VSP probabilistic inversion.              *
 * ----------------------------------------------------------------------------------------------------------------------- */
prior::~prior() = default;

prior::prior() = default;

// Setting prior data manually
prior::prior(AlgebraLib::Vector mean, AlgebraLib::Vector std) {
    _mean = std::move(mean);
    _std = std::move(std);
    _numberParameters = _mean.size();
    setInverseCovarianceMatrix();
}

double prior::misfit(AlgebraLib::Vector parameters) {
    AlgebraLib::Vector parameterDifference = parameters - _mean;
    return 0.5 * (parameterDifference * (_inverseCovarianceMatrix * parameterDifference));
}

AlgebraLib::Vector prior::gradientMisfit(AlgebraLib::Vector parameters) {
    AlgebraLib::Vector parameters_diff = parameters - _mean;
    AlgebraLib::Vector gradient(parameters.size());

    for (int q = 0; q < parameters.size(); q++) {
        gradient[q] = 0.5 * (_inverseCovarianceMatrix.getColumn(q) * parameters_diff +
                             _inverseCovarianceMatrix[q] * parameters_diff);
    }
    return gradient;
}

// Set prior inverse covariance matrix, or mass matrix. Only diagonal entries are filled, no correlation is described.
void prior::setInverseCovarianceMatrix() {
    _inverseCovarianceMatrix = AlgebraLib::Matrix (_numberParameters, _numberParameters);

    for (int i = 0; i < _numberParameters; i++) {
        _inverseCovarianceMatrix[i][i] = 1.0 / (_std[i] * _std[i]);
    }
}

// Copy constructor
prior::prior(const prior &in_prior) {
    _mean = in_prior._mean;
    _std = in_prior._std;
    _inverseCovarianceMatrix = in_prior._inverseCovarianceMatrix;
    _numberParameters = in_prior._numberParameters;
}

/* -----------------------------------------------------------------------------------------------------------------------
 * Data class, for generating new data or loading previously generated data. Also calculates inverse data covariance
 * matrix and is able to compute data misfits.
 * ----------------------------------------------------------------------------------------------------------------------- */
data::data(const char *filename) {
    _observedData = AlgebraLib::ReadVector(filename);
}

data::data(const char *filename, double percentage) {
    _observedData = AlgebraLib::ReadVector(filename);
    _numberData = _observedData.size();
    setICDMatrix_percentual(percentage);
}

data::data() = default;

double data::misfit(AlgebraLib::Vector q, forwardModel m) {
    AlgebraLib::Vector dataDifference = m.calculateData(std::move(q)) - _observedData;
    return 0.5 * (dataDifference * (_inverseCD * dataDifference));
}

void data::setICDMatrix(double std) {
    AlgebraLib::Matrix _inverseCD(_numberData, _numberData);
    for (int i = 0; i < _inverseCD.rows(); i++) {
        _inverseCD[i][i] = 1 / (std * std); // computing inverse variance
    }
}

void data::setICDMatrix_percentual(double percentage) {
    _inverseCD = AlgebraLib::Matrix(_numberData, _numberData);

    for (int i = 0; i < _inverseCD.rows(); i++) {

        double std = fabs(_observedData[i] * (percentage / 100.0));
        _inverseCD[i][i] = 1.0 / (std * std); // computing inverse variance
    }
}

void data::readData(const char *filename) {
    _observedData = AlgebraLib::ReadVector(filename);
}

void data::writeData(const char *filename) {
    AlgebraLib::WriteVector(_observedData, filename);
}

void data::setMisfitParameterDataMatrix(AlgebraLib::Matrix designMatrix) {
    _misfitParameterDataMatrix = designMatrix.Transpose() * _inverseCD;
    setMisfitParameterMatrix(designMatrix);
}

void data::setMisfitParameterMatrix(AlgebraLib::Matrix designMatrix) {
    _misfitParameterMatrix = _misfitParameterDataMatrix * designMatrix;
}

AlgebraLib::Vector data::gradientMisfit(AlgebraLib::Vector parameters) {
    AlgebraLib::Vector gradient;

    for (int q = 0; q < parameters.size(); q++) {
        // I am -fairly- sure of this matrix equation derivative
        gradient[q] = (0.5 * (_misfitParameterMatrix.getColumn(q) * parameters +
                              _misfitParameterMatrix[q] * parameters +
                              -2 * (_misfitParameterDataMatrix[q] * _observedData))
        );
    }
    return gradient;
};

/* -----------------------------------------------------------------------------------------------------------------------
 * Forward model class.
 * ----------------------------------------------------------------------------------------------------------------------- */
void forwardModel::constructUnitDesignMatrix(unsigned long numberParameters) {
    // Make square zero matrix
    _designMatrix = AlgebraLib::Matrix(numberParameters, numberParameters);

    // Set diagonal entries to 1
    for (int i = 0; i < numberParameters; i++) {
        _designMatrix[i][i] = 1.0;
    }
}

forwardModel::forwardModel(unsigned long numberParameters) {
    _numberParameters = numberParameters;
    constructUnitDesignMatrix(_numberParameters);
}

AlgebraLib::Vector forwardModel::calculateData(AlgebraLib::Vector parameters) {
    return (_designMatrix * parameters);
}

forwardModel::forwardModel(const char *filename) {
    // Read file for observed data
    _designMatrix = AlgebraLib::ReadMatrix(filename);
    _numberParameters = _designMatrix.columns();
}

forwardModel::forwardModel() = default;

double posterior::misfit(AlgebraLib::Vector parameters, prior &in_prior, data &in_data, forwardModel m) {
    return in_prior.misfit(parameters) + in_data.misfit(parameters, std::move(m));
}

AlgebraLib::Vector posterior::gradientMisfit(AlgebraLib::Vector parameters, prior &in_prior, data &in_data) {
    return (in_data.gradientMisfit(parameters) + in_prior.gradientMisfit(parameters));
}

#pragma clang diagnostic pop