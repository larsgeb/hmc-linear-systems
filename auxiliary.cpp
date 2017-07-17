// Created by Lars Gebraad on 7/10/17.

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "auxiliary.hpp"
#include "linearalgebra.hpp"

/* -----------------------------------------------------------------------------------------------------------------------
 * Class for Gaussian Distributed prior information about model parameters for a VSP probabilistic inversion.
 * ----------------------------------------------------------------------------------------------------------------------- */
prior::~prior() {}

prior::prior() {}

// Setting prior data manually
prior::prior(std::vector<double> mean, std::vector<double> std) {
    _mean.clear();
    _std.clear();
    _mean = mean;
    _std = std;
    _numberParameters = _mean.size();
    setMassMatrix();
}

double prior::misfit(std::vector<double> q) {
    std::vector<double> parameterDifference = VectorDifference(q, _mean);
    return 0.5 * VectorVectorProduct(parameterDifference, MatrixVectorProduct
            (_massMatrix, parameterDifference));
}

// Set prior inverse covariance matrix, or mass matrix. Only diagonal entries are filled, no correlation is described.
void prior::setMassMatrix() {
    std::vector<double> zeroRow(_numberParameters, 0.0);
    for (int i = 0; i < _numberParameters; i++) {
        _massMatrix.push_back(zeroRow);
        _massMatrix[i][i] = 1 / (_std[i] * _std[i]);
    }
}

// Read prior information on model parameters from file
/*void prior::readPrior(const char *filename) {
    FILE *fid;
    char str[1000];

    fid = fopen(filename, "r");

    int numberLayers;
    fgets(str, 1000, fid); // Move line
    fscanf(fid, "%i", &numberLayers);
    fgets(str, 1000, fid);
    _numberParameters = (unsigned long) (numberLayers * 2 + 1);

    double temp;
    *//*- Layer Speeds -------------------------------------------------------------------------------------------------*//*
    fgets(str, 1000, fid);
    // Scan all prior mean layer speeds (number of layers + half space)
    for (int i = 0; i < numberLayers + 1; i++) {
        fscanf(fid, "%lg", &temp);
        _mean.push_back(temp);
    }
    fgets(str, 1000, fid);
    // Scan all prior sigma layer speeds (number of layers + half space)
    for (int i = 0; i < numberLayers + 1; i++) {
        fscanf(fid, "%lg", &temp);
        _std.push_back(temp);
    }
    fgets(str, 1000, fid);

    *//*- Layer thicknesses --------------------------------------------------------------------------------------------*//*
    fgets(str, 1000, fid);
    for (int i = numberLayers + 1; i < numberLayers * 2 + 1; i++) {
        fscanf(fid, "%lg", &temp);
        _mean.push_back(temp);
    }
    fgets(str, 1000, fid);
    for (int i = numberLayers + 1; i < numberLayers * 2 + 1; i++) {
        fscanf(fid, "%lg", &temp);
        _std.push_back(temp);
    }
    fgets(str, 1000, fid);

    fclose(fid);
}*/

/* -----------------------------------------------------------------------------------------------------------------------
 * Data class, for generating new data or loading previously generated data. Also calculates inverse data covariance
 * matrix and is able to compute data misfits.
 * ----------------------------------------------------------------------------------------------------------------------- */
data::data(int numberData) {
    _numberData = numberData;
    setICDMatrix(0.1); // Hardcoded measurement error. Update as you see fit.
}

data::data() {

}

double data::misfit(std::vector<double> q, forwardModel m) {
    std::vector<double> dataDifference = VectorDifference(m.calculateData(q), _observedData);
    return 0.5 * VectorVectorProduct(dataDifference, MatrixVectorProduct(_inverseCD, dataDifference));
}

void data::setICDMatrix(double std) {
    std::vector<double> zeroRow((unsigned long) _numberData, 0.0);
    for (int i = 0; i < _numberData; i++) {
        _inverseCD.push_back(zeroRow);
        _inverseCD[i][i] = 1 / (std * std); // computing inverse variance
    }
}

void data::readData(const char *filename) {
    // Read file for observed data
    double a;
    _observedData.clear();
    std::ifstream infile(filename);
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
    outfile << _numberData << " ";
    for (int i = 0; i < _numberData; i++) {
        outfile << _observedData[i] << " ";
    }
    outfile.close();
}

void data::setMisfitMatrix(std::vector<std::vector<double>> designMatrix) {
    _misfitMatrix = MatrixMatrixProduct(TransposeMatrix(designMatrix), MatrixMatrixProduct(_inverseCD, designMatrix));
}

/* -----------------------------------------------------------------------------------------------------------------------
 * Forward model class.
 * ----------------------------------------------------------------------------------------------------------------------- */
void forwardModel::constructDesignMatrix(int numberParameters) {
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
    constructDesignMatrix(numberParameters);
}

std::vector<double> forwardModel::calculateData(std::vector<double> parameters) {
    return MatrixVectorProduct(_designMatrix, parameters);
}

forwardModel::forwardModel() {}

/* -----------------------------------------------------------------------------------------------------------------------
 * Taylor expansion class. Contains probability classes by reference.
 * ----------------------------------------------------------------------------------------------------------------------- */
/*void taylorExpansion::calculate0(std::vector<double> expansionPoint) {
    // Zeroth order Taylor-series term is just the function value at the expansion point
    _functionValue = _posterior.misfit(expansionPoint, _prior, _data);
}

void taylorExpansion::calculate1(std::vector<double> expansionPoint) {
    // First order Taylor-series term is the first derivative at the expansion point, a vector containing the derivatives
    // in all parameter dimensions.
    _firstDerivativeValue.clear();
    for (int i = 0; i < _prior._numberParameters; i++) {
        // A central finite difference scheme of order 2
        std::vector<double> point2 = expansionPoint;
        std::vector<double> point1 = expansionPoint;
        point1[i] = expansionPoint[i] * (1 - _stepRatio);
        point2[i] = expansionPoint[i] * (1 + _stepRatio);
        double deltaMisfit = _posterior.misfit(point2, _prior, _data) - _posterior.misfit(point1, _prior, _data);
        _firstDerivativeValue.push_back(deltaMisfit / (2 * expansionPoint[i] * _stepRatio));
        point1.clear();
        point2.clear();
    }
}

void taylorExpansion::calculate2(std::vector<double> expansionPoint) {
    // Second order Taylor-series term is the derivative of the derivative in two directions at the expansion point, a
    // vector containing the derivatives of the first derivatives is contained in a matrix.

    // Generating a zero matrix of dimension _prior.numberParameters.
    _secondDerivativeValue.clear();
    std::vector<double> zeroRow(_prior._numberParameters, 0);
    _secondDerivativeValue.insert(_secondDerivativeValue.end(), _prior._numberParameters, zeroRow);

    // The matrix should be symmetric, which allows for approximately half the actual computation if every element was
    // computed.
    for (int dimension1 = 0; dimension1 < (int) _prior._numberParameters; dimension1++) {
        for (int dimension2 = 0; dimension2 <= dimension1; dimension2++) {

            if (dimension1 == dimension2) {
                // Finite difference stencil for accuracy 2 for second derivative of one variable
                std::vector<double> point1 = expansionPoint, point2 = expansionPoint;
                point2[dimension1] = expansionPoint[dimension1] * (1 + _stepRatio);
                point1[dimension1] = expansionPoint[dimension1] * (1 - _stepRatio);
                _secondDerivativeValue[dimension1][dimension2] =
                        (1 / pow(_stepRatio * expansionPoint[dimension1], 2)) *
                        (_posterior.misfit(point1, _prior, _data) - 2 *
                                                                    _functionValue +
                         _posterior.misfit(point2, _prior, _data));
            } else {
                // Finite difference stencil for accuracy 2 for second derivative of mixed variables
                std::vector<double> point1 = expansionPoint, point2 = expansionPoint, point3 = expansionPoint, point4 = expansionPoint;
                point1[dimension1] *= (1 - _stepRatio);
                point1[dimension2] *= (1 + _stepRatio);

                point2[dimension1] *= (1 + _stepRatio);
                point2[dimension2] *= (1 + _stepRatio);

                point3[dimension1] *= (1 - _stepRatio);
                point3[dimension2] *= (1 - _stepRatio);

                point4[dimension1] *= (1 + _stepRatio);
                point4[dimension2] *= (1 - _stepRatio);

                _secondDerivativeValue[dimension1][dimension2] =
                        (1 / (4 * pow(_stepRatio, 2) * expansionPoint[dimension1] * expansionPoint[dimension2])) *
                        (_posterior.misfit(point3, _prior, _data) + _posterior.misfit(point2, _prior, _data) -
                         _posterior.misfit(point1, _prior, _data) - _posterior.misfit(point4, _prior, _data));
                _secondDerivativeValue[dimension2][dimension1] = _secondDerivativeValue[dimension1][dimension2];
            }

        }
    }

}

taylorExpansion::taylorExpansion() {

}

taylorExpansion::taylorExpansion(std::vector<double> parameters, double stepRatio, prior &in_prior, data &in_data,
                                 posterior &in_posterior) {
    _stepRatio = stepRatio;
    _expansionPoint = parameters;
    _prior = in_prior;
    _data = in_data;
    _posterior = in_posterior;
    calculate0(_expansionPoint);
    calculate1(_expansionPoint);
    calculate2(_expansionPoint);
}

taylorExpansion::~taylorExpansion() {

}

void taylorExpansion::updateExpansion(std::vector<double> in_expansionPoint) {
    _expansionPoint = in_expansionPoint;
    calculate0(in_expansionPoint);
    calculate1(in_expansionPoint);
    calculate2(in_expansionPoint);
}

std::vector<double> taylorExpansion::gradient(std::vector<double> q) {
    // The gradient of the Taylor expansion up to second order is the first order term derivatives plus 2 * the derivative
    // of the second order term evaluated at the current point
    std::vector<double> gradient;
    gradient =
            VectorSum(_firstDerivativeValue,
                      MatrixVectorProduct(_secondDerivativeValue, VectorDifference(q, _expansionPoint)));
    return gradient;
}*/

double posterior::misfit(std::vector<double> parameters, prior &in_prior, data &in_data, forwardModel m) {
    return in_prior.misfit(parameters) + in_data.misfit(parameters, m);
}

//std::vector<double> posterior::gradientMisfit
