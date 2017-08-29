//
// Created by Lars Gebraad on 8/3/17.
//
#include <vector>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include "linearalgebra.hpp"
#include "src/random/randomnumbers.hpp"
#include <iostream>
#include <fstream>
#include <ctime>

int main() {
    std::vector<double> meanC{1,2,3,4,5,6,7,8};
    std::vector<std::vector<double>> covC{{4,  -4, 0,  0,  0,  0,  0,  -1},
                                          {-4, 68, 0,  0,  0,  0,  0,  0},
                                          {0,  0,  4,  -4, 0,  0,  0,  0},
                                          {0,  0,  -4, 68, 0,  0,  0,  0},
                                          {0,  0,  0,  0,  4,  -4, 0,  0},
                                          {0,  0,  0,  0,  -4, 68, 0,  0},
                                          {0,  0,  0,  0,  0,  0,  4,  -4},
                                          {-1,  0,  0,  0,  0,  0,  -4, 10}};

    std::vector<std::vector<double>> Cholesky_covC = CholeskyDecompose(covC);

    std::ofstream outfile;
    outfile.open("OUTPUT/multivariate_samples.txt");
    int samples = 100000;

    outfile << meanC.size() << " " << samples << std::endl;

    for (int i = 0; i < samples; ++i) {
        std::vector<double> C = randn_Cholesky(Cholesky_covC);
        for (double j : C) {
            outfile << j << "  ";
        }
        outfile << std::endl;
    }

    outfile.close();

    return 0;
}