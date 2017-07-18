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
    data observedData(5,1); // Define number of parameters
    observedData.readData("OUTPUT/synthetics.txt");

    // Create prior information (Gaussian distribution)
    std::vector<double> means{2.5, 2.5, 2.5, 2.5, 2.5};
    std::vector<double> std{5,5,5,5,5};
    prior priorInfo(means, std);

    // Create design matrix within forwardModel object
    forwardModel forwardModel1(5); // Define number of parameters
    observedData.setMisfitParameterDataMatrix(forwardModel1._designMatrix); // Check this result with a few awkwardly sized
    // matrices
    // Create posterior object
    posterior posterior1;

    montecarlo mc(priorInfo, observedData, posterior1, forwardModel1, 10, 0.1, 100000);

    /* ---- Works well! ---- */
    // This piece of code investigates the gradient in the neighborhood of expected parameters 1 & 2
    // Use it in conjunction with plot_gradient.py
    /*std::ofstream outfile;
    outfile.open("OUTPUT/gradient.txt");
    for (double q1 = 0; q1 < 5; q1 += 0.25) {
        for (double q2 = 0; q2 < 5; q2 += 0.25) {
            std::vector<double> params{q1, q2, 3, 4, 5};
            std::vector<double> gradient = mc._posterior.gradientMisfit(params, mc._prior, mc._data);
            outfile << q1 << " " << q2 << " " << gradient[0] << " " << gradient[1] << std::endl;
        }
    }
    outfile.close();*/

    /* ---- The actual sampling ---- */
    mc.sample(true);

    return 0;
}