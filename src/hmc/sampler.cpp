//
// Created by Lars Gebraad on 18-8-17.
//
#include <src/random/randomnumbers.hpp>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "forward_model.hpp"
#include "data.hpp"
#include "prior.hpp"
#include "sampler.hpp"

namespace hmc {
    sampler::sampler(InversionSettings settings) {
        // Read the settings
        _nt = settings._trajectorySteps;
        _dt = settings._timeStep;
        _temperature = settings._temperature;
        _proposals = settings._proposals;
        _genMomKinetic = settings._genMomKinetic;
        _genMomPropose = settings._genMomPropose;
        _testBefore = settings._testBefore;
        _massMatrixType = settings._massMatrixType;
        _ergodic = settings._ergodic;
        _hmc = settings._hamiltonianMonteCarlo;
        _outputSamples = settings._outputSamples;
        _outputTrajectory = settings._outputTrajectory;
        _inputMatrix = settings._inputMatrix;
        _inputData = settings._inputData;
        _algorithmNew = settings._algorithmNew;
        _window = settings._window;
        _adaptTimestep = settings._adaptTimestep;

        // Initialise random number generator
        srand((unsigned int) time(nullptr));

        std::cout << std::endl << "Hamiltonian Monte Carlo Sampler" << std::endl << "Lars Gebraad, "
                "2017" << std::endl << "Use --help or -h to display the documentation." << std::endl << std::endl;

        // Load forward matrix
        auto startCPU = get_cpu_time();
        auto startWall = get_wall_time();
        std::cout << "Loading forward matrix ..." << std::endl;
        arma::mat forward_matrix;
        forward_matrix.load(settings._inputMatrix);
        _model = hmc::forward_model(forward_matrix);
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
        _data = hmc::data(_model, synthData, 0.5, true);
        arma::dcolvec means(forward_matrix.n_cols);
        arma::dcolvec std(means.size());
        // generate uniform prior
        for (int i = 0; i < means.size(); i++) {
            means[i] = settings._means;
            std[i] = settings._std_dev;
        }
        _prior = hmc::prior(means, std);

        // Start precomputation
        startCPU = std::clock();
        startWall = get_wall_time();

        // Set precomputed quantities
        _A = _prior._inv_cov_m + _model._g.t() * _data._inv_cov_d * _model._g;
        _bT = _prior._inv_cov_m * _prior._means + _model._g.t() * _data._inv_cov_d * _data._observedData;
        _c = arma::as_scalar(0.5 * (
                _prior._means.t() * _prior._inv_cov_m * _prior._means +
                _data._observedData.t() * _data._inv_cov_d * _data._observedData
        ));

        // Pre-compute mass matrix and other associated quantities
        switch (_massMatrixType) {
            case 0:
                _massMatrix = &_A;
                break;
            case 1:
                _optionalMassMatrixMemory = arma::diagmat(_A);
                _massMatrix = &_optionalMassMatrixMemory;
                break;
            case 2:
                _optionalMassMatrixMemory = arma::eye(size(_A));
                _massMatrix = &_optionalMassMatrixMemory;
                break;
            default:
                _massMatrix = &_A;
        }
        _inverseMassMatrixDiagonal = arma::diagmat(1.0 / (*_massMatrix).diag());

        // Prepare mass matrix decomposition and inverse
        _CholeskyLowerMassMatrix = arma::chol(*_massMatrix, "lower");
        _inverseMassMatrix = arma::inv(*_massMatrix);

        // Set starting proposal
        propose_momentum();
        _proposedModel = randn(_prior._means, arma::vec(arma::mat(_prior._covariance).diag()));

        // Set starting model
        _currentModel = _proposedModel;
        _currentMomentum = _proposedMomentum;

        // Do analysis of the product _A * _massMatrix to determine optimal time step
        if (_adaptTimestep) {
            double maxFrequency;
            arma::vec eigval;
            arma::mat eigvec;
            switch (_massMatrixType) {
                case 0:
                    _dt = (2.0 * PI / _nt);
                    break;
                case 1:
                case 2:

                    eig_sym(eigval, eigvec, _inverseMassMatrix * _A);
                    maxFrequency = arma::max(eigval);
                    eigval.clear();
                    eigvec.clear();
                    _dt = (2.0 * PI / _nt) * 1.0 / sqrt(maxFrequency);
                    break;

                default:
                    break;
            }
        }

        std::cout << "Set up time CPU: " << (std::clock() - startCPU) / (double) (CLOCKS_PER_SEC) << "s, "
                "wall: " << get_wall_time() - startWall << "s" << std::endl << std::endl;

        // Output settings
        std::cout << "Inversion of linear model using MCMC sampling." << std::endl;
        std::cout << "Selected method;      \033[1;34m" << (_hmc ? "Hamiltonian Monte Carlo" : "Metropolis-Hastings")
                  << "\033[0m with following options:" << std::endl;
        std::cout << "\t parameters:        \033[1;32m" << _currentModel.size() << "\033[0m" << std::endl;
        std::cout << "\t proposals:         \033[1;32m" << _proposals << "\033[0m" << std::endl;
        std::cout << "\t temperature:       \033[1;32m" << _temperature << "\033[0m" << std::endl;
        std::cout << "\t timestep:          \033[1;32m" << _dt << "\033[0m" << std::endl;
        std::cout << "\t number of steps:   \033[1;32m" << _nt << "\033[0m" << std::endl << std::endl;
        std::cout << "\t New algorithm:     \033[1;32m" << (_algorithmNew ? "true" : "false") << "\033[0m" << std::endl;
        std::cout << "\t Optimal timestep:  \033[1;32m" << (_adaptTimestep ? "true" : "false") << "\033[0m"
                  << std::endl;
        std::cout << "\t mass matrix type:  \033[1;32m" <<
                  (_massMatrixType == 0 ? "full optimal matrix" :
                   (_massMatrixType == 1 ? "diagonal optimal matrix" : "unit matrix"))
                  << "\033[0m" << std::endl << std::endl;
        std::cout << "\t input matrix:      \033[1;32m" << _inputMatrix << "\033[0m" << std::endl;
        std::cout << "\t input data:        \033[1;32m" << _inputData << "\033[0m" << std::endl;
        std::cout << "\t output samples:    \033[1;32m" << _outputSamples << "\033[0m" << std::endl;
        std::cout << "\t output trajectory: \033[1;32m" << _outputTrajectory << "\033[0m" << std::endl << std::endl;

        if (_testBefore) {
            std::cout << "\t - Exploiting conservation of energy by evaluating before propagation" << std::endl;
        }
        std::cout << "\t - Use generalised mass matrix with" << (_genMomPropose ? "" : "out")
                  << " off diagonal entries" << std::endl;
        if (_genMomKinetic) std::cout << "\t - Use generalised momentum for kinetic energy" << std::endl;
        std::cout << "\t - " << (_ergodic ? "E" : "Not e") << "nforcing ergodicity" << std::endl << std::endl;
    };

    void sampler::setStarting(arma::vec &model) {
        _currentModel = model;
        _proposedModel = model;
    }

    void sampler::propose_metropolis() {
        // This method might not be the actual Metropolis-Hastings method. It should draw from distributions with
        // means at the CURRENT model, not at the prior model.
        _proposedModel = randn(_prior._means, _prior._covariance.t());
    }

    void sampler::propose_momentum() {
        // Draw random prior momenta according to the distribution defined by the mass matrix.
        _proposedMomentum = (_genMomPropose && _massMatrixType == 0) ? randn_Cholesky(_CholeskyLowerMassMatrix) : randn(
                *_massMatrix);
    }

    double sampler::precomp_misfit() {
        return as_scalar(0.5 * _proposedModel.t() * (_A * _proposedModel) - _bT.t() * _proposedModel + _c);
    }

    arma::vec sampler::precomp_misfitGrad() {
        return (_A * _proposedModel - _bT);
    }

    double sampler::kineticEnergy() {
        // If we choose to only use kinetic energy from the diagonal we can use the inverse mass matrix diagonal.
        // In some choices, there's only non-zero numbers on the diagonal when explicitly choosing a specific
        // mass matrix, I hardcode these options to be more efficient.
        if (_genMomKinetic && _massMatrixType != 1 && _massMatrixType != 2) {
            return as_scalar(0.5 * _proposedMomentum.t() * _inverseMassMatrix * _proposedMomentum);
        } else {
            return as_scalar(0.5 * _proposedMomentum.t() * _inverseMassMatrixDiagonal * _proposedMomentum);
        }
    }

    double sampler::chi() {
        return precomp_misfit();
    }

    double sampler::energy() {
        return chi() + kineticEnergy();
    }

    void sampler::sample_new() {
        // Sample the distribution using the modified algorithm
        double x = _hmc ? energy() : chi();
        double x_new;
        int accepted = 1;
        int uturns = 0;

        // Open output file
        std::ofstream samplesfile;
        samplesfile.open(_outputSamples);
        samplesfile << _prior._means.size() << " " << _proposals << std::endl;

        write_sample(samplesfile, x);

        // Write progress to console
        std::cout << "[" << std::setw(3) << (int) (100.0 * double(0) / _proposals) << "%] "
                  << std::string(((unsigned long) ((_window.ws_col - 7) * 0 / _proposals)), *"=") <<
                  "\r" << std::flush;

        for (int it = 1; it < _proposals; it++) {

            // Write progress to console
            if (it % 100 == 0) {
                std::cout << "[" << std::setw(3) << (int) (100.0 * double(it) / _proposals) << "%] "
                          << std::string(((unsigned long) ((_window.ws_col - 7) * it / _proposals)), *"=") <<
                          "\r" << std::flush;
            }

            // Propose
            if (_hmc) {
                propose_momentum();
                if (!_testBefore) {
                    // Write only for the last trajectory  (or any else, modify as fit)
                    leap_frog(uturns, it == _proposals - 1);
                }
            } else {
                propose_metropolis();
            }

            // Evaluate new misfit
            x_new = (_hmc ? energy() : chi());
            double result_exponent = exp((x - x_new) / _temperature);
            if ((x_new < x) || (result_exponent > randf(0.0, 1.0))) {
                if (_testBefore) {
                    leap_frog(uturns, it == _proposals - 1);
                }
                accepted++;
                x = x_new;
                _currentModel = _proposedModel;
                write_sample(samplesfile, x);
            }
        }

        std::cout << "[" << 100 << "%] " << std::string((unsigned long) (_window.ws_col - 7), *"=") << "\r\n"
                  << std::flush;
        std::cout << "Number of accepted models: " << accepted << std::endl;
        std::cout << "Number of U-Turn terminations in propagation: " << uturns << std::endl;

        // Write result
        samplesfile << accepted << std::endl;
        samplesfile.close();
    }

    void sampler::sample_neal() {
        // Sample the distribution using the modified algorithm
        double x = _hmc ? energy() : chi();
        double x_new;
        int accepted = 1;
        int uturns = 0;

        // Open output file
        std::ofstream samplesfile;
        samplesfile.open(_outputSamples);
        samplesfile << _prior._means.size() << " " << _proposals << std::endl;

        write_sample(samplesfile, x);

        // Write progress to console
        std::cout << "[" << std::setw(3) << (int) (100.0 * double(0) / _proposals) << "%] "
                  << std::string(((unsigned long) ((_window.ws_col - 7) * 0 / _proposals)), *"=") <<
                  "\r" << std::flush;

        for (int it = 1; it < _proposals; it++) {

            // Write progress to console
            if (it % 100 == 0) {
                std::cout << "[" << std::setw(3) << (int) (100.0 * double(it) / _proposals) << "%] "
                          << std::string(((unsigned long) ((_window.ws_col - 7) * it / _proposals)), *"=") <<
                          "\r" << std::flush;
            }

            // Propose
            if (_hmc) {
                propose_momentum();
                x = energy(); // this differs from our new algorithm
                // Write only for the last trajectory  (or any else, modify as fit)
                leap_frog(uturns, it == _proposals - 1);
            } else {
                propose_metropolis();
            }

            // Evaluate new misfit
            x_new = (_hmc ? energy() : chi());

            double result_exponent = exp((x - x_new) / _temperature);

            if ((x_new < x) || (result_exponent > randf(0.0, 1.0))) {
                if (_testBefore) {
                    leap_frog(uturns, it == _proposals - 1);
                }
                accepted++;
                x = x_new;
                _currentModel = _proposedModel;
                write_sample(samplesfile, x);
            }
        }

        std::cout << "[" << 100 << "%] " << std::string((unsigned long) (_window.ws_col - 7), *"=") << "\r\n"
                  << std::flush;
        std::cout << "Number of accepted models: " << accepted << std::endl;
        std::cout << "Number of U-Turn terminations in propagation: " << uturns << std::endl;

        // Write result
        samplesfile << accepted << std::endl;
        samplesfile.close();
    }

    void sampler::leap_frog(int &uturns, bool writeTrajectory) {

        // start proposal at current momentum
        _proposedModel = _currentModel;
        // Acts as starting momentum
        _currentMomentum = _proposedMomentum;

        arma::vec misfitGrad;
        double angle1, angle2;

        std::ofstream trajectoryfile;
        if (writeTrajectory) {
            trajectoryfile.open(_outputTrajectory);
            trajectoryfile << _prior._means.size() << " " << _nt << std::endl;
        }

        unsigned long local_nt = _nt;
        double local_dt = _dt;

        if (_ergodic) {
            local_nt = static_cast<unsigned long>(_nt * randf(0.5, 1.5));
            local_dt = _dt * randf(0.5, 1.5);
        }

        for (int it = 0; it < local_nt; it++) {

            misfitGrad = precomp_misfitGrad();
            _proposedMomentum = _proposedMomentum - 0.5 * local_dt * misfitGrad;

            if (writeTrajectory) write_sample(trajectoryfile, chi());

            _proposedModel = _proposedModel + local_dt * (
                    (_genMomKinetic ? _inverseMassMatrix : _inverseMassMatrixDiagonal) * _proposedMomentum);

            misfitGrad = precomp_misfitGrad();
            _proposedMomentum = _proposedMomentum - 0.5 * local_dt * misfitGrad;

            /* Check no-U-turn criterion. */
            angle1 = arma::as_scalar(_proposedMomentum.t() * (_currentModel - _proposedModel));
            angle2 = arma::as_scalar(_currentMomentum.t() * (_proposedModel - _currentModel));

            if (angle1 > 0.0 && angle2 > 0.0) {
                uturns++;
                break;
            }
        }
        if (writeTrajectory) trajectoryfile.close();
    }

    void sampler::write_sample(std::ofstream &outfile, double misfit) {
        for (double j : _proposedModel) {
            outfile << j << "  ";
        }
        outfile << misfit;
        outfile << std::endl;

    }

    void sampler::sample() {
        // Start timers
        auto startCPU = std::clock();
        auto startWall = get_wall_time();

        if (_algorithmNew) {
            sample_new();
        } else {
            sample_neal();
        }

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