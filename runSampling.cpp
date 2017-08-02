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
    // Load the observed data
    double percentualCovariance = 5.0;
    data observedData("OUTPUT/synthetics.txt", percentualCovariance);

    // Create design matrix within forwardModel object
    forwardModel model("INPUT/matrix.txt");
    std::vector<double> means;
    std::vector<double> std;
    for (int i = 0; i < model._numberParameters; i++) {
        means.push_back(10.0 * double(i + 1));
        std.push_back(10.0 * double(i + 1));
    }
    prior priorInfo(means, std);

    // Is unstable as of yet, momentum samples have to be drawn from n-dimensional correlated Gaussian,
    // for which one needs to LU-decompose the mass matrix. This might be implemented later.
    bool boolGeneralisedMomentum = false;

    montecarlo mc(priorInfo, observedData, model, 10, 0.05, 100000, boolGeneralisedMomentum);

    /* ---- The actual sampling ---- */
    std::clock_t start;
    start = std::clock();
    mc.sample(true);
    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    return 0;
}