//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include "randomnumbers.hpp"
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>

int main() {

    // Load the observed data

    data observedData(2, 0.5); // Define number of parameters
    observedData.readData("OUTPUT/synthetics.txt");

    // Create prior information (Gaussian distribution)
    std::vector<double> means{2, 2};
    std::vector<double> std{1, 1};
    prior priorInfo(means, std);

    // Create design matrix within forwardModel object
    forwardModel model(2);

    model._designMatrix[1][1] = 2;

    observedData.setMisfitParameterDataMatrix(model._designMatrix); // Check this result with a few awkwardly sized
    // matrices
    // Create posterior object
    posterior posterior1;

    // Maximum stepsize is approximately smallest constrained dimension of mass
    // matrix, i.e. the square root of the smallest mass.
    montecarlo mc(priorInfo, observedData, posterior1, model, 50, 0.05, 20);

    /* ---- The actual sampling ---- */
    std::clock_t start;
    start = std::clock();
    mc.sample(true);
    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;


    return 0;
}