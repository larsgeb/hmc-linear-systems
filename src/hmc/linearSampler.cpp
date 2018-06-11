/*
 * Created by Lars Gebraad on 18-08-17.
 * Last modified 29-05-18.
 */
#include <cmath>
#include <sstream>
#include <iomanip>
#include "linearSampler.hpp"
#include "../random/randomnumbers.hpp"

using namespace arma;

namespace hmc {
    linearSampler::linearSampler(InversionSettings settings) {
        // Window settings
        _window = settings._window;

        // Output files
        _outputSamples = settings._outputSamplesFile;
        _outputTrajectory = settings._outputTrajectoryFile;

        // Forward model
        A_file = settings.A_file;
        B_file = settings.B_file;
        C_file = settings.C_file;

        // Tuning parameters
        dt = settings._timeStep;
        temperature = settings._temperature;
        proposals = settings._proposals;
        nt = settings._trajectorySteps;
        massMatrixType = settings._massMatrixType;

        // Initialise random number generator
        srand((unsigned int) time(nullptr));

        std::cout << std::endl << "Hamiltonian Monte Carlo Sampler" << std::endl << "Lars Gebraad, 2017" << std::endl
                  << "Use --help or -h to display the documentation." << std::endl << std::endl;

        // Load forward matrix
        auto startCPU = get_cpu_time();
        auto startWall = get_wall_time();
        std::cout << "Loading equation ..." << std::endl;

        A.load(A_file);
        B.load(B_file);
        mat C_mat;
        C_mat.load(C_file);
        C = C_mat[0];

        // Start pre-computation
        startCPU = std::clock();
        startWall = get_wall_time();

        massMatrix = 0.5 * (A + A.t());
        CholeskyLowerMassMatrix = arma::chol(massMatrix, "lower");
        invMass = inv(massMatrix);

        // Set starting proposal
        propose_momentum();
        _proposedModel = -pinv((A + A.t()), 1e-1 * datum::eps, "std") * B;

        // Set starting model
        _currentModel = _proposedModel;
//        _currentMomentum = _proposedMomentum;

        // Do analysis of the product _A * massMatrix to determine optimal time step
        if (settings._adaptTimestep) {
            double maxFrequency;
            arma::vec eigval;
            arma::mat eigvec;
            switch (massMatrixType) {
                case 0:
                    dt = (2.0 * PI / nt);
                    break;
                case 1:
                case 2:
                    eig_sym(eigval, eigvec, invMass * A);
                    maxFrequency = arma::max(eigval);
                    eigval.clear();
                    eigvec.clear();
                    dt = (2.0 * PI / nt) * 0.61497 / sqrt(maxFrequency); // Randomization, if 0.61497 ==> 1: oscillatory samples
                    break;

                default:
                    break;
            }
        }

        std::cout << "Set up time CPU: " << (std::clock() - startCPU) / (double) (CLOCKS_PER_SEC)
                  << "s, wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;

        // Output settings
        std::cout << "Inversion of linear model using MCMC sampling." << std::endl;
        std::cout << "\033[1;34m Hamiltonian Monte Carlo\033[0m with following options:" << std::endl;
        std::cout << "\t parameters:        \033[1;32m" << _currentModel.size() << "\033[0m" << std::endl;
        std::cout << "\t proposals:         \033[1;32m" << proposals << "\033[0m" << std::endl;
        std::cout << "\t temperature:       \033[1;32m" << temperature << "\033[0m" << std::endl;
        std::cout << "\t timestep:          \033[1;32m" << dt << "\033[0m" << std::endl;
        std::cout << "\t number of steps:   \033[1;32m" << nt << "\033[0m" << std::endl << std::endl;
        std::cout << "\t Optimal timestep:  \033[1;32m" << (settings._adaptTimestep ? "true" : "false") << "\033[0m" << std::endl;
//        std::cout << "\t mass matrix type:  \033[1;32m" << (massMatrixType == 0 ? "full optimal matrix" :
//                                                            (massMatrixType == 1 ? "diagonal optimal matrix" : "unit matrix"))
//                  << "\033[0m" << std::endl << std::endl;
        std::cout << "\t output samples:    \033[1;32m" << _outputSamples << "\033[0m" << std::endl;
        std::cout << "\t output trajectory: \033[1;32m" << _outputTrajectory << "\033[0m" << std::endl << std::endl;

    };

    void linearSampler::setStarting(arma::vec &model) {
        _currentModel = model;
        _proposedModel = model;
    }

    void linearSampler::propose_momentum() { // todo replace with general covariance matrix
        // Draw random prior momenta according to the distribution defined by the mass matrix.
//        _proposedMomentum = randn(conv_to<mat>::from(massMatrix));
        _proposedMomentum = randn_Cholesky(CholeskyLowerMassMatrix);
    }

    double linearSampler::misfit() {
        return as_scalar(_proposedModel.t() * (A * _proposedModel) + B.t() * _proposedModel + C);
    }

    double linearSampler::kineticEnergy() {
        return as_scalar(0.5 * _proposedMomentum.t() * invMass * _proposedMomentum);

    }

    double linearSampler::chi() {
        return misfit();
    }

    double linearSampler::energy() {
        return chi() + kineticEnergy();
    }

    void linearSampler::sample_neal() {
        // Sample the distribution using the modified algorithm
        double x = energy();
        double x_new;
        int accepted = 1;

        // Open output file
        std::ofstream samplesfile;
        samplesfile.open(_outputSamples);

        write_sample(samplesfile, chi());

        // Write progress to console
        std::cout << "[" << std::setw(3) << (int) (100.0 * double(0) / proposals) << "%] "
                  << std::string(((unsigned long) ((_window.ws_col - 7) * 0 / proposals)), *"=") <<
                  "\r" << std::flush;


        // Start sampling
        for (int it = 1; it < proposals; it++) {

            // Write progress to console
            if (it % 100 == 0) {
                std::cout << "[" << std::setw(3) << (int) (100.0 * double(it) / proposals) << "%] "
                          << std::string(((unsigned long) ((_window.ws_col - 7) * it / proposals)), *"=") <<
                          "\r" << std::flush;
            }

            // Propose
            propose_momentum();
            x = energy();
            leap_frog(it == proposals - 1);

            x_new = energy();

            double result_exponent = exp((x - x_new) / temperature);
            if ((x_new < x) || (result_exponent > randf(0.0, 1.0))) {
                accepted++;
                x = x_new;
                _currentModel = _proposedModel;
                write_sample(samplesfile, chi());
            }
        }

        std::cout << "[" << 100 << "%] " << std::string((unsigned long) (_window.ws_col - 7), *"=") << "\r\n"
                  << std::flush;
        std::cout << "Number of accepted models: " << accepted << std::endl;

        // Write result
//        samplesfile << accepted << std::endl;
        samplesfile.close();
    }

    void linearSampler::leap_frog(bool writeTrajectory) {

        // start proposal at current momentum
        _proposedModel = _currentModel;
        // Acts as starting momentum
//        _currentMomentum = _proposedMomentum;

        arma::vec misfitGrad;

        std::ofstream trajectoryfile;

        unsigned long local_nt = nt;
        double local_dt = dt;

        local_nt = static_cast<unsigned long>(nt * randf(0.5, 1.5));
        local_dt = dt * randf(0.5, 1.5);

        for (int it = 0; it < local_nt; it++) {

            // this allows for non-symmetric matrices, usually same as 2Ax + B
            misfitGrad = (A.t() * _proposedModel + A * _proposedModel + B);
            _proposedMomentum = _proposedMomentum - 0.5 * local_dt * misfitGrad;

            if (writeTrajectory) write_sample(trajectoryfile, chi());

            _proposedModel = _proposedModel + local_dt * (invMass * _proposedMomentum);

            misfitGrad = (A.t() * _proposedModel + A * _proposedModel + B);
            _proposedMomentum = _proposedMomentum - 0.5 * local_dt * misfitGrad;

        }
        if (writeTrajectory) trajectoryfile.close();
    }

    void linearSampler::write_sample(std::ofstream &outfile, double misfit) {
        for (double j : _proposedModel) {
            outfile << std::setprecision(20) << j << "  ";
        }
        outfile << misfit;
        outfile << std::endl;

    }

    void linearSampler::sample() {
        // Start timers
        auto startCPU = std::clock();
        auto startWall = get_wall_time();

        // Allow for other methods, remnant of old structure
        sample_neal();

        // Output sampling time
        std::cout << "Sampling time CPU: " << (std::clock() - startCPU) / (double) (CLOCKS_PER_SEC)
                  << "s, wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;
    }


} // namespace hmc


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