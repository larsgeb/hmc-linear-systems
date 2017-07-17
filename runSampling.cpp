//
// Created by Lars Gebraad on 7/11/17.
//

#include <vector>
#include "auxiliary.hpp"
#include "linearalgebra.hpp"
#include "montecarlo.hpp"
#include "randomnumbers.hpp"
#include <math.h>

int main() {

    // Load the observed data
    data observedData(5);
    observedData.readData("OUTPUT/synthetics.txt");

    // Create prior information (Gaussian distribution)
    std::vector<double> means{2.5, 2.5, 2.5, 2.5, 2.5};
    std::vector<double> std{2, 2, 2, 2, 2};
    prior priorInfo(means, std);
    means.clear();
    std.clear();

    // Create design matrix within forwardModel object
    forwardModel forwardModel1(5);

    // Create posterior object
    posterior posterior1;

    montecarlo mc(priorInfo, observedData, posterior1, forwardModel1, 20, 0.001, 500);

    std::vector<std::vector<double>> A = {{1,2},{3,4},{5,6}};
    std::vector<std::vector<double>> B = {{1,2,3,4},{5,6,7,8}};

    std::vector<std::vector<double>> P = MatrixMatrixProduct(A,B);

    return 0;
}