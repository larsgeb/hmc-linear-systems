//
// Created by Lars Gebraad on 7/10/17.
// Remake of the original aux.cpp, with a bit more logic to it.
//

#ifndef HMC_VSP_AUXILIARY_HPP
#define HMC_VSP_AUXILIARY_HPP

#include "SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp"

class prior;

class data;

class posterior;

class forwardModel;

/**
 * @brief Prior information in parameter space.
 */
class prior {
public:
    // Fields
    unsigned long _numberParameters;
    AlgebraLib::Vector _mean;
    AlgebraLib::Vector _std;
    AlgebraLib::Matrix _inverseCovarianceMatrix;

    // Constructors and destructor
    prior();

    // Constructor needed for prior information
    explicit prior(AlgebraLib::Vector mean, AlgebraLib::Vector std);

    // Copy constructor
    prior(const prior &);

    ~prior();

    // Member functions
    // Explicit misfit, without pre-computation
    double misfit(AlgebraLib::Vector parameters);

    // Explicit misfit gradient, without pre-computation
    AlgebraLib::Vector gradientMisfit(AlgebraLib::Vector parameters);

private:
    void setInverseCovarianceMatrix();
};

class data {
public:

    data();

    explicit data(const char *filename);

    data(const char *filename, double percentage);

    unsigned long _numberData;
    AlgebraLib::Vector _observedData;
    AlgebraLib::Matrix _inverseCD;
    AlgebraLib::Matrix _misfitParameterDataMatrix;
    AlgebraLib::Matrix _misfitParameterMatrix;

    void setICDMatrix(double std);

    void setICDMatrix_percentual(double percentage);

    void readData(const char *filename);

    void writeData(const char *filename);

    double misfit(AlgebraLib::Vector in_parameters, forwardModel m);

    void setMisfitParameterDataMatrix(AlgebraLib::Matrix designMatrix);

    void setMisfitParameterMatrix(AlgebraLib::Matrix designMatrix);

    AlgebraLib::Vector gradientMisfit(AlgebraLib::Vector parameters);
};

/* ----------------------------------------------------------------------------------------------------------------------- *
 * Class is redundant when pre-computation is used within the monte carlo class.
 * ----------------------------------------------------------------------------------------------------------------------- */
class posterior {
public:
    posterior() = default;

    double misfit(AlgebraLib::Vector parameters, prior &in_prior, data &in_data, forwardModel m);

    AlgebraLib::Vector gradientMisfit(AlgebraLib::Vector parameters, prior &in_prior, data &in_data);
};

class forwardModel {
public:
    // Fields
    unsigned long _numberParameters;
    AlgebraLib::Matrix _designMatrix;

    // Constructors & destructors
    // Constructor which creates a unit forward model of dimensions nP x nP
    explicit forwardModel(unsigned long numberParameters);

    explicit forwardModel(const char *filename);

    forwardModel();

    // Member functions
    void constructUnitDesignMatrix(unsigned long numberParameters);

    AlgebraLib::Vector calculateData(AlgebraLib::Vector parameters);
};


void printVector(AlgebraLib::Vector A);

#endif //HMC_VSP_AUXILIARY_HPP
