//
// Created by Lars Gebraad on 7/11/17.
//
#include "algebra_lib/src/algebra_lib/algebra_lib.hpp"
#include <src/hmc/hmc.hpp>
#include <ctime>

int main() {

    using namespace algebra_lib;
    using namespace hmc;

    vector means(150, true);
    matrix g(means.size(), means.size());
    vector std(means.size(), true);
    for (int i = 0; i < means.size(); i++) {
        g[i][i] = (i + 1.0) * (i + 1.0);
        means[i] = (i + 1.0) ;
        std[i] = (10.0);
    }

    ForwardModel model(g);

    // Load the observed data
    vector synthData = ReadVector("INPUT/synthetics.txt");
    data data(model, synthData, 5.0, true);
    prior prior(means, std);

    GenerateInversionSettings settings;
    settings.setSamples(2000000).setGravity(10);
    sampler sampler(prior, data, model, settings);


    /* ---- The actual sampling ---- */
    std::clock_t start;
    start = std::clock();

    sampler.sample();


    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;

    return EXIT_SUCCESS;
}
