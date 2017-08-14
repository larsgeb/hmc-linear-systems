//
// Created by Lars Gebraad on 7/10/17.
// Creates synthetic data based on setup.txt (for receivers) and parameters_synthetics_model.txt (for model).
//
#include <vector>
#include <iostream>
#include "auxiliary.hpp"

int main() {
    // Easy as can be.
    forwardModel model(150);
//    forwardModel model("INPUT/matrix.txt");

    std::vector<double> parameters;
    for (int i = 0; i < model._numberParameters; i++) {
        parameters.push_back(i + 1.0);
        model._designMatrix[i][i] = (i + 1.0) * (i + 1.0);
    }

    data synthetics;
    synthetics._observedData = model.calculateData(parameters);
    synthetics._numberData = static_cast<int>(synthetics._observedData.size());
    synthetics.writeData("INPUT/synthetics.txt");

    return 0;
}