//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"

int main() {

    // Load the observed data
    data observedData;
    observedData.readData("DATA/synthetics.txt");

    // Load the prior information
    prior priorConstraints("INPUT/prior_information.txt");

    // Creating a posterior object (doesn't contain much except of the misfit function)
    posterior posteriorPDF;

    // Do a Taylor expansion of the misfit function to avoid long computations. Based on a simple forward finite
    // difference scheme.
    std::vector<double> startModel = priorConstraints._mean;
    taylorExpansion expansionOfMisfit(startModel, 1.00000001, priorConstraints, observedData, posteriorPDF);


    return 0;
}