/*
 * Main executable for running sampling of linear models with a provided Gaussian mean.
 * Created by Lars Gebraad on 7/11/17.
 */

#include "src/hmc/hmc.hpp"
#include <ctime>
#include <sstream>
#include <armadillo>
#include <omp.h>
//#include <time.h>
#include <sys/time.h>

double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        // Handle error
        return 0;
    }
    return (double) time.tv_sec + (double) time.tv_usec * .000001;
}

double get_cpu_time() {
    return (double) clock() / CLOCKS_PER_SEC;
}

int main(int argc, char *argv[]) {

    // Standard settings
    hmc::InversionSettings settings;

    // Parse command line flags, it is very easy to implement new flags, but only done out of necessity
    for (int i = 1; i < argc; i++) {
        if (i + 1 != argc) {
            if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--temperature") == 0) {
                settings._temperature = atof(argv[i + 1]);
                i++;
            }

            if (strcmp(argv[i], "-nt") == 0 || strcmp(argv[i], "--trajectorysteps") == 0) {
                settings._trajectorySteps = atof(argv[i + 1]);
                i++;
            }

            if (strcmp(argv[i], "-dt") == 0 || strcmp(argv[i], "--timestep") == 0) {
                settings._timeStep = atof(argv[i + 1]);
                i++;
            }

            if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--ergodic") == 0) {
                std::stringstream ss(argv[i + 1]);
                bool b;

                if (!(ss >> std::boolalpha >> b)) {
                    std::cout << "Incorrect value for boolean parameter (ergodicity), settings to default." <<
                              std::endl;
                } else {
                    settings._ergodic = b;
                }
                i++;
            }

            if (strcmp(argv[i], "-gms") == 0 || strcmp(argv[i], "--generalizedmass") == 0) {
                std::stringstream ss(argv[i + 1]);
                bool b;

                if (!(ss >> std::boolalpha >> b)) {
                    std::cout << "Incorrect value for boolean parameter (generalized mass), settings to default." <<
                              std::endl;
                } else {
                    settings._genMomPropose = b;
                }
                i++;
            }

            if (strcmp(argv[i], "-gmm") == 0 || strcmp(argv[i], "--generalizedmomentum") == 0) {
                std::stringstream ss(argv[i + 1]);
                bool b;

                if (!(ss >> std::boolalpha >> b)) {
                    std::cout << "Incorrect value for boolean parameter (generalized momentum), settings to default." <<
                              std::endl;
                } else {
                    settings._genMomKinetic = b;
                }
                i++;
            }

            if (strcmp(argv[i], "-Hb") == 0 || strcmp(argv[i], "--hamiltonianbefore") == 0) {
                std::stringstream ss(argv[i + 1]);
                bool b;

                if (!(ss >> std::boolalpha >> b)) {
                    std::cout << "Incorrect value for boolean parameter (hamiltonian invariance exploit), settings to "
                            "default." << std::endl;
                } else {
                    settings._testBefore = b;
                }
                i++;
            }
        }
    }


    std::cout << std::endl << "Metropolis Hastings/Hamiltonian Monte Carlo Sampler" << std::endl
              << "Lars Gebraad, 2017" << std::endl << std::endl;
    std::clock_t startCPU;
    auto startWall = get_wall_time();

    // Create prior information
    arma::dcolvec means(121);
    arma::dcolvec std(means.size());
    for (int i = 0; i < means.size(); i++) {
        means[i] = 1.0 / 1500.0; // i%1 > everything 1/1e3, i%2 > alternating checkerboard
        std[i] = 0.001;
        // Increasing this value to 0.1 will allow the algorithm to diverge far from relatively
        // expected speeds, and thus for this underdetermined problem come up with very good fitting models, but
        // unexpected values. Settings it to 0.0001 will pull many values towards the value 1.0/1500.0, giving to
        // much importance to the prior information. As always one has to use `logical' values which are relevant to
        // the data.
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
    hmc::data data(model, synthData, 0.5, true);
    hmc::prior prior(means, std);

    // Settings for the sampler

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


