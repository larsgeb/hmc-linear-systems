//
// Created by Lars Gebraad on 7/10/17.
// Creates synthetic data based on setup.txt (for receivers) and parameters_synthetics_model.txt (for model).
//
#include <vector>
#include "auxiliary.hpp"

int main() {
    data synthetics(5, 1);
    forwardModel model(5);

    std::vector<double> parameters{10, 10, 10, 10, 10};
    synthetics._observedData = model.calculateData(parameters);
    synthetics.writeData("OUTPUT/synthetics.txt");

    return 0;
}