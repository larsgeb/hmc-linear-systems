//
// Created by Lars Gebraad on 7/10/17.
// Creates synthetic data based on setup.txt (for receivers) and parameters_synthetics_model.txt (for model).
//
#include <vector>
#include "auxiliary.hpp"

int main() {
    // Easy as can be.
    data synthetics(2,1); // Needed for location of receivers, hardcoded into class
    forwardModel model(2); // Creates identity matrix, update afterwards using model._designMatrix[i][j]

    model._designMatrix[1][1] = 2;

    std::vector<double> parameters;
    parameters.push_back(1);
    parameters.push_back(3);

    synthetics._observedData = model.calculateData(parameters);
    synthetics.writeData("OUTPUT/synthetics.txt");

    return 0;
}