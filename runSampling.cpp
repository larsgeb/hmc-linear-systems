//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include "randomnumbers.hpp"
#include <iostream>
#include <math.h>

int main() {

    // Load the observed data
    data observedData;
    observedData.readData("DATA/synthetics.txt");

    // Load the prior information
    prior priorConstraints("INPUT/prior_information.txt");

    // Creating a posterior object (doesn't contain much except of the misfit function)
    posterior posteriorPDF;

    // Make montecarlo object for propagating the model through phase space
    montecarlo mc(priorConstraints._mean, priorConstraints, observedData, posteriorPDF, 100, 0.001, 5000);

    bool hamilton = true;
    double x = (hamilton ? mc.chi() : mc.chi());
    double x_new;
    int accepted = 0;

    FILE *pfile;
    pfile = fopen("OUTPUT/samples.txt", "w");
    mc.write_sample(pfile, x, 0);
    for (int it = 1; it < mc._iterations; it++) {
        hamilton ? mc.propose_hamilton() : mc.propose_metropolis();
        x_new = (hamilton ? mc.chi() : mc.chi());

        if ((x_new < x) || (exp(x - x_new) > randf(0.0, 1.0))) {
            accepted++;
            x = x_new;

/*            for (int i = 0; i < mc._proposedModel.size(); i++) {
                std::cout << "Parameter " << i + 1 << ": " << mc._proposedModel[i] << std::endl;
            }
            std::cout << "Proposed model misfit: " << mc._posterior.misfit(mc._proposedModel, priorConstraints,
                                                                           observedData) << std::endl;*/
            mc._currentModel = mc._proposedModel;
        }
        mc.write_sample(pfile, x, it);
    }
    fclose(pfile);
    std::cout << "Number of accepted models: " << accepted;

    return 0;
}