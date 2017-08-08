//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include <fstream>
#include "auxiliary.hpp"
#include "montecarlo.hpp"

int main() {
    try {
        // Load the observed data
        forwardModel forwardModel1(2); // Define number of parameters
        forwardModel1._designMatrix[0][0] = 1;
        forwardModel1._designMatrix[1][1] = 2;
        data observedData("OUTPUT/tomography_synthetics.txt", 0.5); // Define number of parameters

        // Create prior information (Gaussian distribution)
        std::vector<double> means;
        std::vector<double> std;
        for (int parameter = 0; parameter < forwardModel1._numberParameters; parameter++) {
            means.push_back(2);
            std.push_back(1);
        }

        prior priorInfo(means, std);
        // Create design matrix within tomographyForwardModel object
        observedData.setMisfitParameterDataMatrix(forwardModel1._designMatrix);
        // Create posterior object
        posterior posterior1;
        montecarlo mc(priorInfo, observedData, posterior1, forwardModel1, 10, 0.001, 100);

        /* ---- Works well! ---- */
        // This piece of code investigates the gradient in the neighborhood of expected parameters 1 & 2
        // Use it in conjunction with plot_gradient.py
        std::ofstream outfile;
        outfile.open("OUTPUT/gradient.txt");
        for (double q1 = 0.4; q1 <= 2.0; q1 += 0.1) {
            for (double q2 = 2.5; q2 <= 3.4; q2 += 0.1) {
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