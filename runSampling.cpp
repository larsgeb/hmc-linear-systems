//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include <iostream>
#include <math.h>

int main() {

    // Load the observed data
    data observedData;
    observedData.readData("DATA/synthetics.txt");

    // Load the prior information
    prior priorConstraints("INPUT/prior_information.txt");

    // Creating a posterior object (doesn't contain much except of the misfit function)
    posterior posteriorPDF;

    // Make montecarlo object for propagating the model through phase space
    montecarlo mc(priorConstraints._mean, priorConstraints, observedData, posteriorPDF, 100, 0.01, 500);

    std::vector<double> testModel;
    // 2500.0 2600.0 2700.0 2800.0 2900.0 Prior model means (Gaussian distribution)
    // 2500.0 3000.0 3000.0 2500.0 3000.0 Actual model used for the synthetics
    testModel.push_back(1500.0);
    testModel.push_back(1500.0);
    testModel.push_back(1500.0);
    testModel.push_back(1500.0);
    testModel.push_back(1500.0);
    testModel.push_back(500.0);
    testModel.push_back(500.0);
    testModel.push_back(500.0);
    testModel.push_back(500.0);

    // Evaluate the local gradient and local misfit at starting model to check if Taylor Expansion works well. (Done in
    // debugging with breakpoints)
    std::vector<double> localGradient = mc._misfitApproximation.gradient(testModel);
    double localMisfit = mc._posterior.misfit(testModel, priorConstraints, observedData);
    std::cout << "Starting model misfit: " << localMisfit << std::endl;

    // Propose new model
    mc.propose_hamilton();

    // Output new model
    std::cout << "There are " << mc._prior._numberParameters << " parameters. The first "
              << ceil(((float) mc._prior._numberParameters) / 2) << " are layer speeds, the last "
              << floor(((float) mc._prior._numberParameters) / 2) << " are layer thicknesses." << std::endl;
    for (int i = 0; i < mc._proposedModel.size(); i++) {
        std::cout << "Parameter " << i + 1 << ": " << mc._proposedModel[i] << std::endl;
    }
    std::cout << "Proposed model misfit: " << mc._posterior.misfit(mc._proposedModel, priorConstraints, observedData);

    return 0;
}