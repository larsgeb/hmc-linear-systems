//
// Created by Lars Gebraad on 8/3/17.
//
#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include "linearalgebra.hpp"
#include "randomnumbers.hpp"
#include <iostream>
#include <fstream>
#include <ctime>

int main() {
    std::vector<double> meanC{0, 0};
    std::vector<std::vector<double>> covC{{4,  -4},
                                          {-4, 68}};

    std::vector<std::vector<double>> Cholesky_covC = CholeskyDecompose(covC);

    std::ofstream outfile;
    outfile.open("OUTPUT/multivariate.txt");


    for (int i = 0; i < 10000000; ++i) {
        std::vector<double> C = randn(meanC, Cholesky_covC);
        outfile << C[0] << "  " << C[1] << std::endl;
    }

    outfile.close();

    return 0;
}