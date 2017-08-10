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
    data observedData("INPUT/synthetics.txt", percentualCovariance);

    // Create design matrix within forwardModel object
    forwardModel model("INPUT/matrix.txt");
    std::vector<double> means;
    std::vector<double> std;
    for (int i = 0; i < model._numberParameters; i++) {
        means.push_back(15);
        std.push_back(10);
    }
    prior priorInfo(means, std);
    bool boolGeneralisedMomentumPropose = true;
    bool boolGeneralisedMomentumKinetic = true;
    bool boolNormalizeMomentum = false;
    bool evaluateHamiltonianBeforeLeap = true;

    montecarlo mc(priorInfo, observedData, model, 10, 0.1, 100000, boolGeneralisedMomentumPropose,
                  boolGeneralisedMomentumKinetic, boolNormalizeMomentum, evaluateHamiltonianBeforeLeap);

    /* ---- The actual sampling ---- */
    std::clock_t start;
    start = std::clock();
    mc.sample(true);
    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    return 0;
}
