//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include "randomnumbers.hpp"

int main() {

    // Load the observed data
    data observedData;
    observedData.readData("DATA/synthetics.txt");

    // Load the prior information
    prior priorConstraints("INPUT/prior_information.txt");

    // Creating a posterior object (doesn't contain much except of the misfit function)
    posterior posteriorPDF;

    montecarlo mc(priorConstraints._mean, priorConstraints, observedData, posteriorPDF, 100, 5, 500);

    mc.propose_hamilton();

    return 0;
}