// Created by Lars Gebraad on 7/10/17.

#include <utility>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "auxiliary.hpp"
#include "linearalgebra.hpp"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"

/* ----------------------------------------------------------------------------------------------------------------------- *
 * Class for Gaussian Distributed prior information about model parameters for a VSP probabilistic inversion.              *
 * ----------------------------------------------------------------------------------------------------------------------- */
prior::~prior() = default;

prior::prior() = default;

// Setting prior data manually
prior::prior(std::vector<double> mean, std::vector<double> std) {
    _mean.clear();
    _std.clear();
    _mean = std::move(mean);
    _std = std::move(std);
    _numberParameters = _mean.size();
    setInverseCovarianceMatrix();
}

double prior::misfit(std::vector<double> parameters) {
    std::vector<double> parameterDifference = (std::move(parameters) - _mean);
    return 0.5 * (parameterDifference * (_inverseCovarianceMatrix * parameterDifference));
}

std::vector<double> prior::gradientMisfit(std::vector<double> parameters) {
    std::vector<double> parameters_diff = (parameters - _mean);
    std::vector<double> gradient;
    gradient.clear();
    for (int q = 0; q < parameters.size(); q++) {
        gradient.push_back(0.5 * ((GetMatrixColumn(_inverseCovarianceMatrix, q) * parameters_diff) +
                                  (GetMatrixRow(_inverseCovarianceMatrix, q) * parameters_diff)));
    }
    return gradient;
}

// Set prior inverse covariance matrix, or mass matrix. Only diagonal entries are filled, no correlation is described.
void prior::setInverseCovarianceMatrix() {
    std::vector<double> zeroRow(_numberParameters, 0.0);
    for (int i = 0; i < _numberParameters; i++) {
        _inverseCovarianceMatrix.push_back(zeroRow);
        _inverseCovarianceMatrix[i][i] = 1 / (_std[i] * _std[i]);
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
    readData(filename);
}

data::data(const char *filename, double percentage) {
    readData(filename);
    setICDMatrix_percentual(percentage);
}

data::data() = default;

double data::misfit(std::vector<double> q, forwardModel m) {
    std::vector<double> dataDifference = m.calculateData(std::move(q)) - _observedData;
    return 0.5 * (dataDifference * (_inverseCD * dataDifference));
}

void data::setICDMatrix(double std) {
    std::vector<double> zeroRow((unsigned long) _numberData, 0.0);
    _inverseCD.clear();
    for (int i = 0; i < _numberData; i++) {
        _inverseCD.push_back(zeroRow);
        _inverseCD[i][i] = 1 / (std * std); // computing inverse variance
    }
}

void data::setICDMatrix_percentual(double percentage) {
    double std;
    std::vector<double> zeroRow((unsigned long) _numberData, 0.0);
    _inverseCD.clear();
    for (int i = 0; i < _numberData; i++) {
        _inverseCD.push_back(zeroRow);
        std = _observedData[i] * (percentage / 100.0);
        _inverseCD[i][i] = 1 / (std * std); // computing inverse variance
    }
}

void data::readData(const char *filename) {
    // Read file for observed data
    double a;
    _observedData.clear();
    std::ifstream infile(filename);

    // Ignore first two lines
    infile.ignore(500,'\n');
    infile.ignore(500,'\n');
    infile.ignore(500,'\n');

    infile >> _numberData;
    for (int i = 0; i < _numberData; i++) {
        infile >> a;
        _observedData.push_back(a);
    }
    infile.close();
}

void data::writeData(const char *filename) {
    // Write data
    std::ofstream outfile;
    outfile.open(filename);

    outfile << "# Data generated from hardcoded class, which reads from matrix.txt" << std::endl;
    outfile << "# Line 4: number of data points. Line 5: data in sequential order." << std::endl;
    outfile << "# The first three lines are always ignored, no matter the characters (up to 500 characters per line)." << std::endl;

    outfile << _numberData << std::endl;
    for (int i = 0; i < _numberData; i++) {
        outfile << _observedData[i] << " ";
    }
    outfile.close();
}

void data::setMisfitParameterDataMatrix(std::vector<std::vector<double>> designMatrix) {
    _misfitParameterDataMatrix = TransposeMatrix(designMatrix) * _inverseCD;
    setMisfitParameterMatrix(designMatrix);
}

void data::setMisfitParameterMatrix(std::vector<std::vector<double>> designMatrix) {
    _misfitParameterMatrix = _misfitParameterDataMatrix * std::move(designMatrix);
}

std::vector<double> data::gradientMisfit(std::vector<double> parameters) {
    std::vector<double> gradient;
    gradient.clear();
    for (int q = 0; q < parameters.size(); q++) {
        // I am -fairly- sure of this matrix equation derivative
        gradient.push_back(0.5 * (GetMatrixColumn(_misfitParameterMatrix, q) * parameters +
                                  GetMatrixRow(_misfitParameterMatrix, q) * parameters +
                                  -2 * (GetMatrixRow(_misfitParameterDataMatrix, q) * _observedData)
        ));
    }
    return gradient;
};

/* -----------------------------------------------------------------------------------------------------------------------
 * Forward model class.
 * ----------------------------------------------------------------------------------------------------------------------- */
void forwardModel::constructUnitDesignMatrix(int numberParameters) {
    // Make square zero matrix
    _designMatrix.clear();
    std::vector<double> zeroRow((unsigned long) numberParameters, 0);
    _designMatrix.insert(_designMatrix.end(), (unsigned long) numberParameters, zeroRow);

    // Set diagonal entries to 1
    for (int i = 0; i < numberParameters; i++) {
        _designMatrix[i][i] = 1;
    }
}

forwardModel::forwardModel(int numberParameters) {
    _numberParameters = numberParameters;
    constructUnitDesignMatrix(numberParameters);
}

std::vector<double> forwardModel::calculateData(std::vector<double> parameters) {
    return (_designMatrix * std::move(parameters));
}

forwardModel::forwardModel(const char *filename) {
    // Read file for observed data
    double element;
    std::vector<double> row;
    int numberData;
    _designMatrix.clear();

    std::ifstream infile(filename);

    // Ignore first two lines
    infile.ignore(500,'\n');
    infile.ignore(500,'\n');
    infile.ignore(500,'\n');

    infile >> numberData;
    infile >> _numberParameters;
    for (int i = 0; i < numberData; i++) {
        for (int j = 0; j < _numberParameters; j++) {
            infile >> element;
            row.push_back(element);
        }
        _designMatrix.push_back(row);
        row.clear();
    }
    infile.close();
}

forwardModel::forwardModel() = default;

double posterior::misfit(std::vector<double> parameters, prior &in_prior, data &in_data, forwardModel m) {
    return in_prior.misfit(parameters) + in_data.misfit(parameters, std::move(m));
}

std::vector<double> posterior::gradientMisfit(std::vector<double> parameters, prior &in_prior, data &in_data) {
    return (in_data.gradientMisfit(parameters) + in_prior.gradientMisfit(parameters));
}

void printVector(std::vector<double> A) {
    for (int i = 0; i < A.size(); i++) {
        std::cout << "Component " << i + 1 << ": " << A[i] << std::endl;
    }
}

#pragma clang diagnostic pop