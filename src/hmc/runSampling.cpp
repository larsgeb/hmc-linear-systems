/*
 * Main executable for running sampling of linear models with a provided Gaussian mean.
 * Created by Lars Gebraad on 7/11/17.
 */

#include "src/hmc/hmc.hpp"

double get_wall_time();

double get_cpu_time();

int main(int argc, char *argv[]) {

    // Standard settings
    hmc::InversionSettings settings(argc, argv);

    std::cout << std::endl << "Metropolis Hastings/Hamiltonian Monte Carlo Sampler" << std::endl << "Lars Gebraad, "
            "2017" << std::endl << "Use --help or -h to display the documentation." << std::endl << std::endl;

    // Load forward matrix
    auto startCPU = get_cpu_time();
    auto startWall = get_wall_time();
    std::cout << "Loading forward matrix ..." << std::endl;
    arma::mat forward_matrix;
    forward_matrix.load(settings._inputMatrix);
    hmc::forward_model model(forward_matrix);
    std::cout << "Forward matrix loading time CPU: " << (std::clock() - startCPU) / (double)
            (CLOCKS_PER_SEC) << "s, wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;


    // Load the observed data
    startCPU = std::clock();
    startWall = get_wall_time();
    std::cout << "Loading data ..." << std::endl;
    arma::vec synthData;
    synthData.load(settings._inputData);
    std::cout << "Data loading time CPU: " << (std::clock() - startCPU) / (double) (CLOCKS_PER_SEC) <<
              "s, wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;

    // Loading data and prior objects
    startCPU = std::clock();
    startWall = get_wall_time();
    std::cout << "Generating required matrices and setting up the inversion ..." << std::endl;
    hmc::data data(model, synthData, 0.5, true);
    arma::dcolvec means(forward_matrix.n_cols);
    arma::dcolvec std(means.size());
    for (int i = 0; i < means.size(); i++) {
        means[i] = settings._means;
        std[i] = settings._std_dev;
    }
    hmc::prior prior(means, std);

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

double get_wall_time() {
    struct timeval time{};
    if (gettimeofday(&time, nullptr)) {
        // Handle error
        return 0;
    }
    return (double) time.tv_sec + (double) time.tv_usec * .000001;
}

double get_cpu_time() {
    return (double) clock() / CLOCKS_PER_SEC;
}