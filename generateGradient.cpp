//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include <fstream>
#include "auxiliary.hpp"
#include "montecarlo.hpp"

int main() {
    try {
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


        /* ---- Works well! ---- */
        // This piece of code investigates the gradient in the neighborhood of expected parameters 1 & 2
        // Use it in conjunction with plot_gradient.py
        std::ofstream outfile;
        outfile.open("OUTPUT/gradient.txt");
        for (double q1 = 0.4; q1 <= 2.0; q1 += 0.1) {
            for (double q2 = 2.5; q2 <= 3.4; q2 += 0.1) {
                means[0] = q1;
                means[1] = q2;

                std::vector<double> gradient = mc.precomp_misfitGrad(means);
                outfile << q1 << " " << q2 << " " << gradient[0] << " " << gradient[1] << std::endl;
            }
        }
        outfile.close();

        return EXIT_SUCCESS;
    } catch (std::exception *) {
        return EXIT_FAILURE;
    }
}