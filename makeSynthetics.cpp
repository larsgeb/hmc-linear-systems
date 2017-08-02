//
// Created by Lars Gebraad on 7/10/17.
// Creates synthetic data based on setup.txt (for receivers) and parameters_synthetics_model.txt (for model).
//
#include <vector>
#include "auxiliary.hpp"

int main() {
    int numberParameters = 10;

    // Easy as can be.
    data synthetics(numberParameters, 1); // Needed for location of receivers, hardcoded into class
    forwardModel model(numberParameters); // Creates identity matrix, update afterwards using model._designMatrix[i][j]
    
    std::vector<double> parameters;
    parameters.reserve(static_cast<unsigned long>(numberParameters));
    for (int i = 0; i < numberParameters; i++) {
        model._designMatrix[i][i] = (i+1)*(i+1);
        parameters.push_back(double(i+1));
    }

    synthetics._observedData = model.calculateData(parameters);
    synthetics.writeData("OUTPUT/synthetics.txt");

    return 0;
}