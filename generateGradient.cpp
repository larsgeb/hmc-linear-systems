//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"

int main() {
    try {



        // Load the observed data
        forwardModel forwardModel1(5); // Define number of parameters
        data observedData(forwardModel1._numberData, 0.001); // Define number of parameters
        observedData.readData("straight-ray-solver/DATA_synthetics", forwardModel1._numberSources,
                              forwardModel1._numberReceivers);

        // Create prior information (Gaussian distribution)
        std::vector<double> means;
        std::vector<double> std;
        for (int parameter = 0; parameter < forwardModel1._numberParameters; parameter++) {
            means.push_back(0.0001);
            std.push_back(0.00005);
        }

        prior priorInfo(means, std);
        // Create design matrix within tomographyForwardModel object
        observedData.setMisfitParameterDataMatrix(
                forwardModel1._designMatrix); // Check this result with a few awkwardly sized
        // matrices
        // Create posterior object
        posterior posterior1;
        montecarlo mc(priorInfo, observedData, posterior1, forwardModel1, 10, 0.001, 100);

        /* ---- Works well! ---- */
        // This piece of code investigates the gradient in the neighborhood of expected parameters 1 & 2
        // Use it in conjunction with plot_gradient.py
        std::ofstream outfile;
        outfile.open("OUTPUT/gradient.txt");
        for (double q1 = 0.0006; q1 < 0.0014; q1 += 0.00005) {
            for (double q2 = 0.0006; q2 < 0.0014; q2 += 0.00005) {
                means[0] = q1;
                means[1] = q2;
                std::vector<double> gradient = mc._posterior.gradientMisfit(means, mc._prior, mc._data);
                outfile << q1 << " " << q2 << " " << gradient[0] << " " << gradient[1] << std::endl;
            }
        }
        outfile.close();

        return EXIT_SUCCESS;
    } catch (std::exception *) {
        return EXIT_FAILURE;
    }
}