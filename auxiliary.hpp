//
// Created by Lars Gebraad on 7/10/17.
// Remake of the original aux.cpp, with a bit more logic to it.
//

#ifndef HMC_VSP_AUXILIARY_HPP
#define HMC_VSP_AUXILIARY_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "auxiliary.hpp"
#include "linearalgebra.hpp"

// Update this include for the wanted forward model, along with the constructors
#include "tomographyForwardModel.hpp"

class prior;

class data;

class posterior;

//class taylorExpansion;

class prior {
public:
    // Fields
    unsigned long _numberParameters;
    std::vector<double> _mean;
    std::vector<double> _std;
    // Mind that as other masses are assigned, the function prior::misfitGradient should actually use the inverse covariance
    // matrix, which is not explicitly defined.
    std::vector<std::vector<double> > _inverseCovarianceMatrix;

    // Constructors and destructor
    prior();

    prior(std::vector<double> mean, std::vector<double> std);

    ~prior();

    // Member functions
    double misfit(std::vector<double> parameters);

    std::vector<double> gradientMisfit(std::vector<double> parameters);

private:

    // Mind that as other masses are assigned, the function prior::misfitGradient should actually use the inverse covariance
    // matrix, which is not explicitly defined.
    void setMassMatrix();
};

class data {
public:
    data(int numberData, double measurementError);

    data();

    int _numberData;
    std::vector<double> _observedData;
    std::vector<std::vector<double> > _inverseCD;
    std::vector<std::vector<double> > _misfitParameterDataMatrix; // This is the forward model transposed times inverseCD
    std::vector<std::vector<double> > _misfitParameterMatrix; // This is the forward model transposed times inverseCD times forward
    // model. Pre-calculated for speed.

    void setICDMatrix(double std);

    void readData(const char *folder, int numberSources, int numberReceivers);

    void writeData(const char *filename);

    double misfit(std::vector<double> in_parameters, forwardModel m);

    void setMisfitParameterDataMatrix(std::vector<std::vector<double>> designMatrix);

    void setMisfitParameterMatrix(std::vector<std::vector<double>> designMatrix);

    std::vector<double> gradientMisfit(std::vector<double> parameters);
};

class posterior {
public:
    double misfit(std::vector<double> parameters, prior &in_prior, data &in_data, forwardModel m);

    std::vector<double> gradientMisfit(std::vector<double> parameters, prior &in_prior, data &in_data);
};

void printVector(std::vector<double> A);

/*class taylorExpansion {
public:
    taylorExpansion();

    // Constructor and destructor
    taylorExpansion(std::vector<double> parameters, double stepRatio, prior &in_prior, data &in_data, posterior
    &in_posterior);

    ~taylorExpansion();

    // Fields
    double _stepRatio;
    std::vector<double> _expansionPoint;
    double _functionValue;
    std::vector<double> _firstDerivativeValue;
    std::vector<std::vector<double> > _secondDerivativeValue;
    prior _prior;
    data _data;
    posterior _posterior;

    // Member functions
    std::vector<double> gradient(std::vector<double> q);

    void updateExpansion(std::vector<double> in_expansionPoint);

private:
    void calculate0(std::vector<double> expansionPoint);

    void calculate1(std::vector<double> expansionPoint);

    void calculate2(std::vector<double> expansionPoint);
};*/

#endif //HMC_VSP_AUXILIARY_HPP
