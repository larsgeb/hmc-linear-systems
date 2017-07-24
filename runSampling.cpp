//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "linearalgebra.hpp"
#include "montecarlo.hpp"
#include "randomnumbers.hpp"
#include <math.h>
#include <iostream>
#include <fstream>

int main() {

    // Load the observed data
    data observedData(2, 0.5); // Define number of parameters
    observedData.readData("OUTPUT/synthetics.txt");

    // Create prior information (Gaussian distribution)
    std::vector<double> means{2, 2};
    std::vector<double> std{1, 1};
    prior priorInfo(means, std);

    // Create design matrix within forwardModel object
    forwardModel model(2); // Creates identity matrix, update afterwards using model._designMatrix[i][j]

    model._designMatrix[1][1] = 2;
//    model._designMatrix[2][2] = 3;
//    model._designMatrix[3][3] = 4;
//    model._designMatrix[4][4] = 5;

    observedData.setMisfitParameterDataMatrix(model._designMatrix); // Check this result with a few awkwardly sized
    // matrices
    // Create posterior object
    posterior posterior1;

    montecarlo mc(priorInfo, observedData, posterior1, model, 10, 0.05, 1000);

    /* ---- The actual sampling ---- */
    mc.sample(true);

    return 0;
}