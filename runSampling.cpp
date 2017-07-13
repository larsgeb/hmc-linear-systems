//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include "randomnumbers.hpp"
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
    montecarlo mc(priorConstraints._mean, priorConstraints, observedData, posteriorPDF, 500, 0.00001, 50000);

    bool hamilton = true;
    mc.sample(hamilton);


    return 0;
}