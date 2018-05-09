/*
 * Main executable for running sampling of linear models with a provided Gaussian mean.
 * Created by Lars Gebraad on 7/11/17.
 */

#include <cstdlib>
#include "linearSampler.hpp"

int main(int argc, char *argv[]) {

    // Standard settings
    hmc::InversionSettings settings(argc, argv);

    // Creating the linearSampler
    hmc::linearSampler sampler(settings);



    // Running the actual linearSampler
    sampler.sample();
//
//    arma::mat m = {"0.984699  0.877115  0.585528  0.0173705  -0.918611  1  1  1  1  1  1  1  1  1  1  1  1"};
//    m = m.t();
//
//    std::cout << 0.5 * m.t() * sampler.A * m + m.t() * sampler.B + sampler.C;
//
//    std::cout << std::endl << sampler.massMatrix << std::endl;
//
//    std::cout << std::endl << 2 * sampler.A * m - sampler.B << std::endl;

//    std::cout << sampler.massMatrix;

    return EXIT_SUCCESS;
}
