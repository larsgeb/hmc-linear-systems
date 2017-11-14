//
// Created by Lars Gebraad on 7/11/17.
//
#include "src/hmc/hmc.hpp"
#include <ctime>
#include <armadillo/armadillo-8.200.2/include/armadillo>
#include <omp.h>

int main() {

    std::cout << std::endl << "Metropolis Hastings/Hamiltonian Monte Carlo Sampler" << std::endl
        << "Lars Gebraad, 2017" << std::endl << std::endl;

    std::clock_t start;
    start = std::clock();

    // Create prior information
    std::cout << "Loading forward matrix ..." << std::endl;
    arma::dcolvec means(10201);
    arma::dcolvec std(means.size());
    for (int i = 0; i < means.size(); i++) {
        means[i] = i%1 == 1 ? 1.0 / 1000.0 : 1.0 / 2000.0; // i%1 > everything 1/1e3, i%2 > alternating checkerboard
        std[i] = 0.01;
    }

    arma::mat forward_matrix;
    forward_matrix.load("INPUT/matrix_checkerboard_lr_100x100_arma.txt");

    hmc::forward_model model(forward_matrix);

    std::cout << "Forward matrix loading time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl << std::endl;
    start = std::clock();    

    // Load the observed data
    // TODO replace workaround with proper loading
    std::cout << "Loading data ..." << std::endl;
    arma::rowvec synthDataRow;
    synthDataRow.load("INPUT/Recorded_time_sources_checkerboard_lr_100x100_arma.txt");
    arma::vec synthData = synthDataRow.t();
    synthDataRow.clear();

    std::cout << "Data loading time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl << std::endl;
    start = std::clock();

    std::cout << "Generating required matrices and setting up the inversion ..." << std::endl;
    hmc::data data(model, synthData, 5.0, true);
    hmc::prior prior(means, std);
    hmc::GenerateInversionSettings settings;

    settings
            .setSamples(5000)
            .setGravity(1); // Choose it such that oscillations are around explorative (no slow exploration)
    settings
            .setTimeStep(0.1); // ideally below 2 std_dev for 1d problems for stability
    settings
            .setErgodicity(true)
            .setHamiltonianMonteCarlo(true)
            .setGenMomPropose(true)
            .setGenMomKinetic(true)
            .setOutfile(const_cast<char *>("OUTPUT/samples1.txt"));

    settings.setAcceptanceFactor(0.0000000000000001);
    settings.setSamples(10000);	
    
    hmc::sampler sampler1(prior, data, model, settings);
	
    std::cout << "Set up time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl << std::endl;
    start = std::clock();

    sampler1.sample();

    std::cout << std::endl << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s" << std::endl;

    return EXIT_SUCCESS;
}
