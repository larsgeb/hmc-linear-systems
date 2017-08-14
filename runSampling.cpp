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
    forwardModel model(150);
//    forwardModel model("INPUT/matrix.txt");
    std::vector<double> means;
    std::vector<double> std;
    for (int i = 0; i < model._numberParameters; i++) {
        model._designMatrix[i][i] = (i + 1.0) * (i + 1.0);
        means.push_back(35.0);
        std.push_back(100.0);
    }
    prior priorInfo(means, std);
    bool boolGeneralisedMomentumPropose = true;
    bool boolGeneralisedMomentumKinetic = false;
    bool boolNormalizeMomentum = false;
    bool evaluateHamiltonianBeforeLeap = true;
    bool boolHMC = true;


    int numberProposals = 100000000;

    montecarlo mc(priorInfo, observedData, model, 10, 0.1, numberProposals, boolGeneralisedMomentumPropose,
                  boolGeneralisedMomentumKinetic, boolNormalizeMomentum, evaluateHamiltonianBeforeLeap);
//    "\033[1;31mbold red text\033[0m"
    /* ---- The actual sampling ---- */
    std::cout << "Inversion of linear model using MCMC sampling." << std::endl;
    std::cout << "Selected method; \033[1;34m" << (boolHMC ? "HMC" : "Metropolis-Hastings") << "\033[0m with following options:"
              << std::endl;
    std::cout << "\t parameters:   \033[1;32m" << model._numberParameters << "\033[0m" << std::endl;
    std::cout << "\t proposals:    \033[1;32m" << numberProposals << "\033[0m" << std::endl;

    if (evaluateHamiltonianBeforeLeap)
        std::cout << "\t - Exploiting Hamiltonian invariance by evaluating before propagation" << std::endl;
    std::cout << "\t - Use generalised mass matrix with" << (boolGeneralisedMomentumPropose ? "" : "out")
              << " off diagonal entries" << std::endl;
    if (boolGeneralisedMomentumKinetic) std::cout << "\t - Use generalised momentum for kinetic energy" << std::endl;

    std::clock_t start;
    start = std::clock();
    mc.sample(boolHMC);
    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;
    return 0;
}
