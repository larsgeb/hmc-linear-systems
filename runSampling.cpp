//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"

int main() {

    // Load the observed data
    forwardModel forwardModel1(5); // Define number of parameters
    data observedData(forwardModel1._numberData, 1); // Define number of parameters
    observedData.readData("straight-ray-solver/DATA_synthetics", forwardModel1._numberSources,
                          forwardModel1._numberReceivers);

    // Create prior information (Gaussian distribution)
    std::vector<double> means;
    std::vector<double> std;
    for (int parameter = 0; parameter < forwardModel1._numberParameters; parameter++) {
        means.push_back(0.00066666);
        std.push_back(0.00066666*0.1);
    }
    prior priorInfo(means, std);

    // Create design matrix within tomographyForwardModel object
    observedData.setMisfitParameterDataMatrix(forwardModel1._designMatrix); // Check this result with a few awkwardly sized
    // matrices
    // Create posterior object
    posterior posterior1;

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
    montecarlo mc(priorInfo, observedData, posterior1, forwardModel1, 10, 0.01, 100);
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