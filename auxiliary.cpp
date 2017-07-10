// Created by Lars Gebraad on 7/10/17.

#include <vector>
#include <iostream>
#include <fstream>
#include "auxiliary.hpp"
#include "linearalgebra.cpp"

/* -----------------------------------------------------------------------------------------------------------------------
 * Class for Gaussian Distributed prior information about model parameters for a VSP probabilistic inversion.
 * ----------------------------------------------------------------------------------------------------------------------- */
prior::~prior() {}

// Setting prior
prior::prior(const char *filename) {
    readPrior(filename);
    setMassMatrix();
}

// Setting prior data manually
prior::prior(std::vector<double> mean, std::vector<double> std,
             std::vector<double>) {
    _mean = mean;
    _std = std;
    _numberParameters = _mean.size();
    mean.clear();
    std.clear();
    setMassMatrix();
}

double prior::misfit(std::vector<double> q) {
    std::vector<double> parameterDifference = VectorDifference(q, _mean);
    return 0.5 * VectorVectorProduct(parameterDifference, MatrixVectorProduct
            (_massMatrix, parameterDifference));
}

// Set prior inverse covariance matrix, or mass matrix
void prior::setMassMatrix() {
    std::vector<double> zeroRow(_numberParameters, 0.0);
    for (int i = 0; i < _numberParameters; i++) {
        _massMatrix.push_back(zeroRow);
        _massMatrix[i][i] = 1 / (_std[i] * _std[i]);
    }
}

// Read prior information on model parameters from file
void prior::readPrior(const char *filename) {
    FILE *fid;
    char str[1000];

    fid = fopen(filename, "r");

    int numberLayers;
    fgets(str, 1000, fid); // Move line
    fscanf(fid, "%i", &numberLayers);
    fgets(str, 1000, fid);
    _numberParameters = (unsigned long) (numberLayers * 2 + 1);

    double temp;
    /*- Layer Speeds -------------------------------------------------------------------------------------------------*/
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

    /*- Layer thicknesses --------------------------------------------------------------------------------------------*/
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
}

/* -----------------------------------------------------------------------------------------------------------------------
 * Data class, for generating new data or loading previously generated data. Also calculates inverse data covariance
 * matrix and is able to compute data misfits.
 * ----------------------------------------------------------------------------------------------------------------------- */
data::data() {
    FILE *fid2;
    char str[1000];

    /* Read VSP setup. */
    fid2 = fopen("INPUT/setup.txt", "r");

    double receiverSpacing, initialLocation, dataError;

    fgets(str, 1000, fid2);
    fscanf(fid2, "%i %lf %lf %lf", &_numberReceivers, &receiverSpacing,
           &initialLocation, &dataError);

    for (int i = 0; i < _numberReceivers; i++) {
        _depthReceivers.push_back(initialLocation + receiverSpacing * i);
    }

    setICDMatrix(dataError); // Hardcoded measurement error. Update as you see fit.
    fclose(fid2);
}

void data::setICDMatrix(double std) {
    std::vector<double> zeroRow((unsigned long) _numberReceivers, 0.0);
    for (int i = 0; i < _numberReceivers; i++) {
        _inverseCD.push_back(zeroRow);
        _inverseCD[i][i] = 1 / (std * std);
    }
}

/* -----------------------------------------------------------------------------------------------------------------------
 * Forward model class. Is only used statically.
 * ----------------------------------------------------------------------------------------------------------------------- */
std::vector<double> forwardVSP::forwardModel(std::vector<double> parameters, std::vector<double> locations) {
    // This is the forward model code for a simple VSP where first arrivals are calculated from layer thicknesses and speeds.
    // Remember, positive Z is downward
    // parameters is filled in the following way: 0 - nLayers: speeds || nLayers+1 - nLayers*2: thicknesses
    std::vector<double> travelTime;
    std::vector<double> top;
    std::vector<double> base;
    std::vector<double> thickness;
    std::vector<double> speed;
    int numberLayers = (int) ((parameters.size() - 1) / 2);

    top.push_back(0.0);
    for (int i = 0; i < numberLayers; i++) {
        thickness.push_back(parameters[numberLayers+1+i]);
        base.push_back(top[i] + thickness[i]);
        top.push_back(base[i]);
    }
    travelTime.clear();
    for (int i = 0; i < locations.size(); i++) {// Looping through receivers
        travelTime.push_back(0.0);
        for (int j = 0; j < numberLayers; j++) {// Looping through layers
            if (locations[i] >= base[j]) {
                travelTime[i] += thickness[j] / parameters[j];
            } else if (locations[i] > top[j]) {
                travelTime[i] += (locations[i] - top[j]) / parameters[j];
            }
        }
        if (locations[i] > base[numberLayers - 1]) {
            travelTime[i] += (locations[i] - base[numberLayers - 1]) /  parameters[numberLayers-1];
        }
    }
    return travelTime;
}

void forwardVSP::writeData(const char *filename, std::vector<double> travelTime, std::vector<double> locations) {
    FILE *fid;

    fid = fopen(filename, "w");

    for (int i = 0; i < locations.size(); i++) {
        fprintf(fid, "%.15lg %.15lg", locations[i], travelTime[i]);
        fprintf(fid, "\n");
    }
    fclose(fid);
}

void forwardVSP::generateSynthetics(const char *filename, std::vector<double> locations) {
    std::vector<double > parameters;

    FILE *fid;
    char str[1000];
    int numberLayers;

    fid = fopen(filename, "r");

    fgets(str, 1000, fid); // Move line
    fscanf(fid, "%d", &numberLayers);
    fgets(str, 1000, fid);

    int numberParameters = numberLayers * 2 + 1;
    double temp;

    /*- Layer Speeds -------------------------------------------------------------------------------------------------*/
    fgets(str, 1000, fid);
    // Scan all priori mean layer speeds (number of layers + half space)
    for (int i = 0; i < numberLayers + 1; i++) {
        fscanf(fid, "%lg", &temp);
        parameters.push_back(temp);
    }
    fgets(str, 1000, fid);

    /*- Layer thicknesses --------------------------------------------------------------------------------------------*/
    fgets(str, 1000, fid);
    for (int i = numberLayers + 1; i < numberLayers * 2 + 1; i++) {
        fscanf(fid, "%lg", &temp);
        parameters.push_back(temp);
    }
    fgets(str, 1000, fid);

    fclose(fid);
    forwardVSP::writeData("DATA/synthetics.txt",forwardVSP::forwardModel(parameters,locations),locations);
}

/* -----------------------------------------------------------------------------------------------------------------------
 * Taylor expansion class. Contains probability classes by reference.
 * ----------------------------------------------------------------------------------------------------------------------- */
void taylorExpansion::calculate0() {

}

void taylorExpansion::calculate1() {

}

void taylorExpansion::calculate2() {

}

double posterior::misfit(std::vector<double> parameters, prior &in_prior, data &in_data) {
    return in_prior.misfit(parameters);
}
