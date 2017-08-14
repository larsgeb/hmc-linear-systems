//
// Created by Lars Gebraad on 7/10/17.
// Remake of the original aux.cpp, with a bit more logic to it.
//

#ifndef HMC_VSP_AUXILIARY_HPP
#define HMC_VSP_AUXILIARY_HPP

class prior;

class data;

class posterior;

class forwardModel;

/**
 * @brief Prior information in parameter space.
 * \addtogroup Probabilisticinversion classes
 */
class prior {
public:
    // Fields
    unsigned long _numberParameters;
    std::vector<double> _mean;
    std::vector<double> _std;
    std::vector<std::vector<double> > _inverseCovarianceMatrix;

    // Constructors and destructor
    prior();

    // Constructor needed for prior information
    explicit prior(std::vector<double> mean, std::vector<double> std);

    // Copy constructor
    prior(const prior &);

    ~prior();

    // Member functions
    // Explicit misfit, without pre-computation
    double misfit(std::vector<double> parameters);

    // Explicit misfit gradient, without pre-computation
    std::vector<double> gradientMisfit(std::vector<double> parameters);

private:
    void setInverseCovarianceMatrix();
};

class data {
public:

    data();

    explicit data(const char *filename);

    data(const char *filename, double percentage);

    int _numberData;
    std::vector<double> _observedData;
    std::vector<std::vector<double> > _inverseCD;
    std::vector<std::vector<double> > _misfitParameterDataMatrix;
    std::vector<std::vector<double> > _misfitParameterMatrix;

    void setICDMatrix(double std);

    void setICDMatrix_percentual(double percentage);

    void readData(const char *filename);

    void writeData(const char *filename);

    double misfit(std::vector<double> in_parameters, forwardModel m);

    void setMisfitParameterDataMatrix(std::vector<std::vector<double>> designMatrix);

    void setMisfitParameterMatrix(std::vector<std::vector<double>> designMatrix);

    std::vector<double> gradientMisfit(std::vector<double> parameters);
};

/* ----------------------------------------------------------------------------------------------------------------------- *
 * Class is redundant when pre-computation is used within the monte carlo class.
 * ----------------------------------------------------------------------------------------------------------------------- */
class posterior {
public:
    posterior() = default;

    double misfit(std::vector<double> parameters, prior &in_prior, data &in_data, forwardModel m);

    std::vector<double> gradientMisfit(std::vector<double> parameters, prior &in_prior, data &in_data);
};

class forwardModel {
public:
    // Fields
    int _numberParameters;
    std::vector<std::vector<double>> _designMatrix;

    // Constructors & destructors
    // Constructor which creates a unit forward model of dimensions nP x nP
    explicit forwardModel(int numberParameters);

    explicit forwardModel(const char *filename);

    forwardModel();

    // Member functions
    void constructUnitDesignMatrix(int numberParameters);

    std::vector<double> calculateData(std::vector<double> parameters);
};


void printVector(std::vector<double> A);

#endif //HMC_VSP_AUXILIARY_HPP
