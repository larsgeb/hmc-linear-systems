//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"

int main() {

    // Load the observed data
    forwardModel forwardModel1(5); // Define number of parameters
    data observedData(forwardModel1._numberData, 0.001); // Define number of parameters
    observedData.readData("straight-ray-solver/DATA_synthetics", forwardModel1._numberSources,
                          forwardModel1._numberReceivers);

    // Create prior information (Gaussian distribution)
    std::vector<double> means;
    std::vector<double> std;
    for (int parameter = 0; parameter < forwardModel1._numberParameters; parameter++) {
        means.push_back(0.002);
        std.push_back(0.002);
    }
    prior priorInfo(means, std);

    // Create design matrix within tomographyForwardModel object
    observedData.setMisfitParameterDataMatrix(forwardModel1._designMatrix); // Check this result with a few awkwardly sized
    // matrices
    // Create posterior object
    posterior posterior1;
    montecarlo mc(priorInfo, observedData, posterior1, forwardModel1, 1, 0.0001, 500);

    /* ---- The actual sampling ---- */
    clock_t start;
    start = clock();
    std::cout << "Metropolis Hastings sampling: " << std::endl;
    mc.sample(false);
    std::cout << ", elapsed time: " << (double) (clock() - start) / CLOCKS_PER_SEC << std::endl;

    start = clock();
    std::cout << "Hamiltonian sampling: " << std::endl;
    mc.sample(true);
    std::cout << ", elapsed time: " << (double) (clock() - start) / CLOCKS_PER_SEC << std::endl;

    return 0;
}