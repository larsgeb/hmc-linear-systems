//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include <iostream>
#include <ctime>

int main() {
    // Load the observed data
    double percentualCovariance = 5.0;
    data observedData("INPUT/tomography_synthetics.txt", percentualCovariance);

    // Create design matrix within forwardModel object
    forwardModel model("INPUT/tomography_matrix.txt");
    std::vector<double> means;
    std::vector<double> std;
    for (int i = 0; i < model._numberParameters; i++) {
        means.push_back(1.0/1500.0);
        std.push_back(0.0005);
    }
    prior priorInfo(means, std);

    bool boolGeneralisedMomentumPropose = false;
    bool boolGeneralisedMomentumKinetic = false;

    montecarlo mc(priorInfo, observedData, model, 10, 0.1, 1000000, boolGeneralisedMomentumPropose,
                  boolGeneralisedMomentumKinetic);

    /* ---- The actual sampling ---- */
    std::clock_t start;
    start = std::clock();
    mc.sample(true);
    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    return 0;
}