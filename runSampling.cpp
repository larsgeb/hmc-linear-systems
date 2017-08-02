//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include <iostream>
#include <fstream>
#include <ctime>

int main() {

    int numberParameters = 10;

    // Load the observed data
    data observedData(10, 1); // Define number of observables
    observedData.readData("OUTPUT/synthetics.txt");
    observedData.setICDMatrix_percentual(10.0);

    // Create prior information (Gaussian distribution)
    std::vector<double> means;
    std::vector<double> std;

    // Create design matrix within forwardModel object
    forwardModel model(numberParameters);

    for (int i = 0; i < numberParameters; i++) {
        model._designMatrix[i][i] = (i + 1) * (i + 1);
        means.push_back(50);
        std.push_back(50);
    }

    prior priorInfo(means, std);

    observedData.setMisfitParameterDataMatrix(model._designMatrix);

    // Create posterior object
    posterior posterior1;

    // Maximum step size is approximately smallest constrained dimension of mass
    // matrix, i.e. the square root of the smallest mass.
    montecarlo mc(priorInfo, observedData, posterior1, model, 10, 0.05, 1000);

    /* ---- The actual sampling ---- */
    std::clock_t start;
    start = std::clock();
    mc.sample(true);
    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    return 0;
}