//
// Created by Lars Gebraad on 7/10/17.
// Creates synthetic data based on setup.txt (for receivers) and parameters_synthetics_model.txt (for model).
//
#include <vector>
#include "auxiliary.hpp"

int main() {
    // Easy as can be.
    data synthetics(5); // Needed for location of receivers, hardcoded into class
    forwardModel model(5);

    std::vector<double> parameters;
    parameters.push_back(1);
    parameters.push_back(2);
    parameters.push_back(3);
    parameters.push_back(4);
    parameters.push_back(5);

    synthetics._observedData = model.calculateData(parameters);
    synthetics.writeData("OUTPUT/synthetics.txt");

    return 0;
}