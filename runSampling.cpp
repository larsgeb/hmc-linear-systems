//
// Created by Lars Gebraad on 7/11/17.
//
#include "src/hmc/hmc.hpp"
#include <ctime>
#include <armadillo>
#include <omp.h>
#include <time.h>
#include <sys/time.h>

double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        //  Handle error
        return 0;
    }
    return (double) time.tv_sec + (double) time.tv_usec * .000001;
}

double get_cpu_time() {
    return (double) clock() / CLOCKS_PER_SEC;
}

int main() {
    std::cout << std::endl << "Metropolis Hastings/Hamiltonian Monte Carlo Sampler" << std::endl
              << "Lars Gebraad, 2017" << std::endl << std::endl;
    std::clock_t startCPU;
    auto startWall = get_wall_time();

    // Create prior information
    arma::dcolvec means(121);
    arma::dcolvec std(means.size());
    for (int i = 0; i < means.size(); i++) {
        means[i] = 1.0 / 1000.0; // i%1 > everything 1/1e3, i%2 > alternating checkerboard
        std[i] = 0.01;
    }

    // Load forward matrix
    startCPU = std::clock();
    startWall = get_wall_time();
    std::cout << "Loading forward matrix ..." << std::endl;
    arma::mat forward_matrix;
    forward_matrix.load("INPUT/matrix_checkerboard_lr_10x10_arma.txt");
    hmc::forward_model model(forward_matrix);
    std::cout << "Forward matrix loading time CPU: " << (std::clock() - startCPU) / (double)
            (CLOCKS_PER_SEC) << "s, wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;


    // Load the observed data
    // TODO replace workaround with proper loading
    startCPU = std::clock();
    startWall = get_wall_time();
    std::cout << "Loading data ..." << std::endl;
    arma::rowvec synthDataRow;
    synthDataRow.load("INPUT/Recorded_time_sources_checkerboard_lr_10x10_arma.txt");
    arma::vec synthData = synthDataRow.t();
    synthDataRow.clear();
    std::cout << "Data loading time CPU: " << (std::clock() - startCPU) / (double) (CLOCKS_PER_SEC) <<
              "s, wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;

    // Loading data and prior objects
    startCPU = std::clock();
    startWall = get_wall_time();
    std::cout << "Generating required matrices and setting up the inversion ..." << std::endl;
    hmc::data data(model, synthData, 5.0, true);
    hmc::prior prior(means, std);

    // Settings for the sampler
    hmc::GenerateInversionSettings settings;
    settings.setSamples(10000);

    // Choose it such that oscillations are around explorative (no slow exploration)
    settings.setGravity(0.0001);

    // ideally below 2 std_dev for 1d problems for stability, based on PDE analysis
    settings.setTimeStepFromGrav_nSteps();

    settings.setErgodicity(true);
    settings.setAcceptanceFactor(0.0001);
    settings.setHamiltonianMonteCarlo(true);
    settings.setGenMomPropose(true);
    settings.setGenMomKinetic(true);
    settings.setOutfile(const_cast<char *>("OUTPUT/samples1.txt"));;

    // Creating the sampler
    hmc::sampler sampler1(prior, data, model, settings);
    std::cout << "Set up time CPU: " << (std::clock() - startCPU) / (double) (CLOCKS_PER_SEC) << "s, "
            "wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;

    // Running the actual sampler
    startCPU = std::clock();
    startWall = get_wall_time();
    sampler1.sample();
    std::cout << std::endl << "Sampling time CPU: " << (std::clock() - startCPU) / (double) (CLOCKS_PER_SEC) << "s, "
            "wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;

    return EXIT_SUCCESS;
}


