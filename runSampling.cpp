//
// Created by Lars Gebraad on 7/11/17.
//
#include "AlgebraLib/src/algebra_lib/algebra_lib.hpp"
#include "src/hmc/hmc.hpp"
#include <ctime>

int main() {

    using namespace algebra_lib;
    using namespace hmc;

    // Create prior information
    sparse_vector means(10201, true);
    sparse_vector std(means.size(), true);
    for (unsigned int i = 0; i < means.size(); i++) {
        means(i) = 1.0 / 1500.0;
        std(i) = 0.006;
    }

    sparse_matrix forward_matrix = ReadSparseMatrix("INPUT/matrix_checkerboard_lr_100x100.txt");
    forward_model model(forward_matrix);

    // Load the observed data
    sparse_vector synthData = ReadSparseVector("INPUT/Recorded_time_sources_100x100.txt");
    data data(model, synthData, 5.0, true);
    prior prior(means, std);

    GenerateInversionSettings settings;
    settings
            .setSamples(1000)
            .setGravity(0.0000001) // Choose it such that oscillations are around explorative (no slow exploration)
            .setTimeStep(0.0001) // ideally below 2 std_dev for 1d problems for stability
            .setErgodicity(true)
            .setHamiltonianMonteCarlo(true)
            .setGenMomPropose(true)
            .setGenMomKinetic(true)
            .setOutfile(const_cast<char *>("OUTPUT/samples1.txt"));

    sampler sampler1(prior, data, model, settings);

    /* ---- The actual sampling ---- */
    std::clock_t start;
    start = std::clock();

    sampler1.sample();

    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;

    return EXIT_SUCCESS;
}
