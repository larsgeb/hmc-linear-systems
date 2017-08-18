//
// Created by Lars Gebraad on 7/11/17.
//
#include "SparseLinearAlgebra/src/AlgebraLib/AlgebraLib.hpp"
#include <src/HMC/HMC.hpp>
#include <ctime>

int main() {

    AlgebraLib::Vector means(150, true);
    AlgebraLib::Matrix g(means.size(), means.size());
    AlgebraLib::Vector std(means.size(), true);
    for (int i = 0; i < means.size(); i++) {
        g[i][i] = (i + 1.0) * (i + 1.0);
        means[i] = (35.0);
        std[i] = (100.0);
    }

    HMC::ForwardModel model(g);

    // Load the observed data
    AlgebraLib::Vector synthData = AlgebraLib::ReadVector("INPUT/synthetics.txt");
    HMC::Data data(model, synthData, 5.0, true);

    AlgebraLib::Vector synth = AlgebraLib::ReadVector("INPUT/synthetics.txt");

    HMC::Prior prior(means, std);

    HMC::GenerateInversionSettings settings = HMC::GenerateInversionSettings()
            .setSamples(100000)
            .setGravity(1.0);

    HMC::Sampler sampler(prior, data, model, settings);

    { // Some nicely formatted settings.
        std::cout << "Inversion of linear model using MCMC sampling." << std::endl;
        std::cout << "Selected method; \033[1;34m" << (settings._hamiltonianMonteCarlo ? "HMC" : "Metropolis-Hastings")
                  << "\033[0m with following options:"
                  << std::endl;
        std::cout << "\t parameters:   \033[1;32m" << means.size() << "\033[0m" << std::endl;
        std::cout << "\t proposals:    \033[1;32m" << settings._proposals << "\033[0m" << std::endl;

        if (settings._testBefore) {
            std::cout << "\t - Exploiting conservation of energy by evaluating before propagation" << std::endl;
        }
        std::cout << "\t - Use generalised mass matrix with" << (settings._genMomPropose ? "" : "out")
                  << " off diagonal entries" << std::endl;
        if (settings._genMomKinetic) std::cout << "\t - Use generalised momentum for kinetic energy" << std::endl;
    }

    /* ---- The actual sampling ---- */
    std::clock_t start;
    start = std::clock();
    sampler.sample(settings._hamiltonianMonteCarlo);
    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;

    return EXIT_SUCCESS;
}
