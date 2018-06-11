/*
 * Main executable for running sampling of linear models q(x) = xt A x + bt x + c
 * Created by Lars Gebraad on 07-11-17.
 * Last modified 29-05-18.
 */

#include <cstdlib>
#include "../hmc/linearSampler.hpp"

int main(int argc, char *argv[]) {

    // Standard settings
    hmc::InversionSettings settings(argc, argv);

    // Creating the linearSampler
    hmc::linearSampler sampler(settings);

    // Running the actual linearSampler
    sampler.sample();

    return EXIT_SUCCESS;
}
