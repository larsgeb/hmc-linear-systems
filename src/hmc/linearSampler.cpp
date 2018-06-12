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
        window = settings._window;

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

        // Show version
        std::cout << std::endl << "Hamiltonian Monte Carlo Sampler" << std::endl << "Lars Gebraad, version 2 - Summer 2018" << std::endl
                  << "Use --help or -h to display the documentation." << std::endl << std::endl;

        // Load quadratic form
        auto startCPU = get_cpu_time();
        auto startWall = get_wall_time();
        std::cout << "Loading equation ..." << std::endl;
        A.load(A_file);
        B.load(B_file);
        mat C_mat;
        C_mat.load(C_file);
        C = C_mat[0];
        std::cout << "Matrices loaded." << std::endl;

        // Start pre-computation
        startCPU = std::clock();
        startWall = get_wall_time();
        // Perform mass pre-computations
        At = A.t();
        massMatrix = 0.5 * (A + At);
        if (arma::approx_equal(massMatrix, A, "rel_tol", 0.01)) {
            symmetricA = true;
            // If matrix is symmetric, we don't need to `remember' its transpose.
            At.clear();
        }
        // Check for mass matrix type and modify accordingly.
        if (massMatrixType == 1) {
            massMatrix = diagvec(massMatrix);
        } else if (massMatrixType == 2) {
            // This is not optimally placed, as when we truly want to use eye, we shouldn't want to compute A+At for
            // the mass matrix. However, this is still necessary to calculate the symmetry condition.
            massMatrix = ones(massMatrix.n_rows, 1);
        }
        // Perform necessary precomputations
        if (massMatrixType == 0) {
            std::cout << "Performing Cholesky decomposition." << std::endl;
            CholeskyLowerMassMatrix = arma::chol(massMatrix, "lower");
            std::cout << "Performed Cholesky decomposition." << std::endl;
            std::cout << "Inverting Cholesky decomposition." << std::endl;
            mat invChol = inv(CholeskyLowerMassMatrix);
            std::cout << "Inverted Cholesky decomposition." << std::endl;
            std::cout << "Inverting mass using Cholesky decomposition." << std::endl;
            invMass = invChol.t() * invChol;
            std::cout << "Inverted mass using Cholesky decomposition." << std::endl;
        } else {
            invMass = 1.0 / massMatrix;
        }

        // Set starting proposal
        propose_momentum();
        _proposedModel = -pinv((A + A.t()), 1e-1 * datum::eps, "std") * B;

        // Set starting model
        _currentModel = _proposedModel;

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
                    eig_sym(eigval, eigvec, diagmat(invMass) * A);
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
        std::cout << "\t mass matrix type:  \033[1;32m" << (massMatrixType == 0 ? "full optimal matrix" :
                                                            (massMatrixType == 1 ? "diagonal optimal matrix" : "unit matrix"))
                  << "\033[0m" << std::endl << std::endl;
        std::cout << "\t output samples:    \033[1;32m" << _outputSamples << "\033[0m" << std::endl;
        std::cout << "\t output trajectory: \033[1;32m" << _outputTrajectory << "\033[0m" << std::endl;
        std::cout << "\t Diagonal matrix:   \033[1;32m" << (symmetricA ? "yes" : "no") << "\033[0m" << std::endl << std::endl;
    };

    void linearSampler::setStarting(arma::vec &model) {
        _currentModel = model;
        _proposedModel = model;
    }

    void linearSampler::propose_momentum() {
        // Draw random prior momenta according to the distribution defined by the mass matrix.
        if (massMatrixType == 0) {
            _proposedMomentum = randn_Cholesky(CholeskyLowerMassMatrix);
        } else {
            _proposedMomentum = randn(conv_to<vec>::from(massMatrix));
        }
    }

    double linearSampler::misfit() {
        return as_scalar(_proposedModel.t() * (A * _proposedModel) + B.t() * _proposedModel + C);
    }

    double linearSampler::kineticEnergy() {
        return 0.5 * ((massMatrixType == 0) ?
                      as_scalar(_proposedMomentum.t() * invMass * _proposedMomentum) :
                      dot(_proposedMomentum % invMass, _proposedMomentum));

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

        // Open output file and write starting model
        std::ofstream samplesfile;
        samplesfile.open(_outputSamples);
        write_sample(samplesfile, chi());

        // Write progress in percentages to console
        std::cout << "[" << std::setw(3) << (int) (100.0 * double(0) / proposals) << "%] "
                  << std::string(((unsigned long) ((window.ws_col - 7) * 0 / proposals)), *"=") <<
                  "\r" << std::flush;

        // Perform sampling
        for (int it = 1; it < proposals; it++) {
            // Write progress to console every 100 steps
            if (it % 100 == 0) {
                std::cout << "[" << std::setw(3) << (int) (100.0 * double(it) / proposals) << "%] "
                          << std::string(((unsigned long) ((window.ws_col - 7) * it / proposals)), *"=") <<
                          "\r" << std::flush;
            }

            // Propose new momentum and propagate
            propose_momentum();
            x = energy();
            leap_frog(it == proposals - 1);

            // Calculate new Hamiltonian
            x_new = energy();

            // Evaluate acceptance criterion
            double result_exponent = exp((x - x_new) / temperature);
            if ((x_new < x) || (result_exponent > randf(0.0, 1.0))) {
                accepted++;
                x = x_new;
                _currentModel = _proposedModel;
                write_sample(samplesfile, chi());
            }
        }

        // Write out 100% at the end
        std::cout << "[" << 100 << "%] " << std::string((unsigned long) (window.ws_col - 7), *"=") << "\r\n"
                  << std::flush;
        std::cout << "Number of accepted models: " << accepted << std::endl;

        // Close output file
        samplesfile.close();
    }

    void linearSampler::leap_frog(bool writeTrajectory) {

        // start proposal at current momentum
        _proposedModel = _currentModel;

        arma::vec misfitGrad;
        std::ofstream trajectoryfile;

        // Randomize settings as to ensure ergodicity
        unsigned long local_nt = nt;
        double local_dt = dt;
        local_nt = static_cast<unsigned long>(nt * randf(0.5, 1.5));
        local_dt = dt * randf(0.5, 1.5);

        // Time integrate Hamiltons equations
        for (int it = 0; it < local_nt; it++) {

            // this allows for non-symmetric matrices, usually same as 2Ax + B
            misfitGrad = (A.t() * _proposedModel + A * _proposedModel + B);
            _proposedMomentum = _proposedMomentum - 0.5 * local_dt * misfitGrad;
            if (writeTrajectory) write_sample(trajectoryfile, chi());
            _proposedModel = _proposedModel + local_dt *
                                              (massMatrixType == 0 ?
                                               conv_to<vec>::from(invMass * _proposedMomentum) :
                                               conv_to<vec>::from(invMass % _proposedMomentum));

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