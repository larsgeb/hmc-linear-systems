//
// Created by Lars Gebraad on 7/10/17.
// Creates synthetic data based on setup.txt (for receivers) and parameters_synthetics_model.txt (for model).
//
#include <vector>
#include <iostream>
#include "auxiliary.hpp"

int main() {
    // Easy as can be.
    forwardModel model("INPUT/matrix.txt");

    std::vector<double> parameters;
    for (int i = 0; i < model._numberParameters; i++) {
        parameters.push_back(10.0 * double(i + 1));
    }

    data synthetics;
    synthetics._observedData = model.calculateData(parameters);
    synthetics._numberData = static_cast<int>(synthetics._observedData.size());
    synthetics.writeData("OUTPUT/synthetics.txt");

    return 0;
}