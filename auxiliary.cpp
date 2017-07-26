// Created by Lars Gebraad on 7/10/17.

#include <sstream>
#include "auxiliary.hpp"

/* ----------------------------------------------------------------------------------------------------------------------- *
 * Class for Gaussian Distributed prior information about model parameters for a VSP probabilistic inversion.              *
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
    setInverseCovarianceMatrix();
}

double prior::misfit(std::vector<double> parameters) {
    std::vector<double> parameterDifference = VectorDifference(parameters, _mean);
    return 0.5 * VectorVectorProduct(parameterDifference, MatrixVectorProduct
            (_inverseCovarianceMatrix, parameterDifference));
}

std::vector<double> prior::gradientMisfit(std::vector<double> parameters) {
    std::vector<double> parameters_diff = VectorDifference(parameters, _mean);
    std::vector<double> gradient;
    gradient.clear();
    for (int q = 0; q < parameters.size(); q++) {
        // Mind that as other masses are assigned, this formula should actually use the inverse covariance matrix, which
        // is not explicitly defined. TODO Implement the inverse covariance matrix explicitly
        gradient.push_back(0.5 * (VectorVectorProduct(GetMatrixColumn(_inverseCovarianceMatrix, q), parameters_diff) +
                                  VectorVectorProduct(GetMatrixRow(_inverseCovarianceMatrix, q), parameters_diff)));
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

/* -----------------------------------------------------------------------------------------------------------------------
 * Data class, for generating new data or loading previously generated data. Also calculates inverse data covariance
 * matrix and is able to compute data misfits.
 * ----------------------------------------------------------------------------------------------------------------------- */
data::data(int numberData, double measurementError) {
    _numberData = numberData;
    setICDMatrix(measurementError); // Hardcoded measurement error. Update as you see fit.
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

void data::readData(const char *folder, int numberSources, int numberReceivers) {
    // Read file for observed data
    _observedData.clear();
    for (int source = 0; source < numberSources; source++) {
        double a;
        std::stringstream filenameStream;
        filenameStream << "straight-ray-solver/DATA_synthetic/Recorded_time_source_" << source + 1 << ".txt";
        std::string filename = filenameStream.str();
        std::ifstream infile(filename);
        for (int receiver = 0; receiver < numberReceivers; receiver++) {
            infile >> a;
            infile >> a;
            infile >> a;
            _observedData.push_back(a);
        }
        infile.close();
    }
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

void data::setMisfitParameterDataMatrix(std::vector<std::vector<double>> designMatrix) {
    _misfitParameterDataMatrix = MatrixMatrixProduct(TransposeMatrix(designMatrix), _inverseCD);
    setMisfitParameterMatrix(designMatrix);
}

void data::setMisfitParameterMatrix(std::vector<std::vector<double>> designMatrix) {
    _misfitParameterMatrix = MatrixMatrixProduct(_misfitParameterDataMatrix, designMatrix);
}

std::vector<double> data::gradientMisfit(std::vector<double> parameters) {
    std::vector<double> gradient;
    gradient.clear();
    for (int q = 0; q < parameters.size(); q++) {
        // I am -fairly- sure of this matrix equation derivative
        gradient.push_back(0.5 * (
                                   VectorVectorProduct(GetMatrixColumn(_misfitParameterMatrix, q), parameters) +
                                   VectorVectorProduct(GetMatrixRow(_misfitParameterMatrix, q), parameters) +
                                   - 2 * VectorVectorProduct(GetMatrixRow(_misfitParameterDataMatrix,q),_observedData)
                           )
        );
    }
    return gradient;
};



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

std::vector<double> posterior::gradientMisfit(std::vector<double> parameters, prior &in_prior, data &in_data) {
    return VectorSum(in_data.gradientMisfit(parameters), in_prior.gradientMisfit(parameters));
}

void printVector(std::vector<double> A) {
    for (int i = 0; i < A.size(); i++) {
        std::cout << "Component " << i + 1 << ": " << A[i] << std::endl;
    }
}
