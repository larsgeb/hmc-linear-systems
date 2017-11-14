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
    sampler::sampler(prior &prior, data &data, forward_model &model, GenerateInversionSettings settings) :
            _data(data), _prior(prior), _model(model) {
        _nt = settings._trajectorySteps;
        _dt = settings._timeStep;
        _acceptanceFactor = settings._acceptanceFactor;
        _gravity = settings._gravity;
        _proposals = settings._proposals;
        _genMomKinetic = settings._genMomKinetic;
        _genMomPropose = settings._genMomPropose;
        _norMom = settings._norMom;
        _testBefore = settings._testBefore;
        _window = settings._window;
        _ergodic = settings._ergodic;
        _hmc = settings._hamiltonianMonteCarlo;
        _outfile = settings._outfile;

        /* Initialise random number generator. ----------------------------------------*/
        srand((unsigned int) time(nullptr));

        // Pre-compute mass matrix and other associated quantities
        _massMatrix = arma::conv_to<arma::mat>::from(
                _gravity * (_prior._inv_cov_m + ((_model._g.t() * _data._inv_cov_d) * _model._g))
        );
        _inverseMassMatrixDiagonal = arma::diagmat(1.0 / _massMatrix.diag());

        // Prepare mass matrix decomposition and inverse.
        _CholeskyLowerMassMatrix = arma::chol(_massMatrix, "lower");
        _inverseMassMatrix = arma::inv(_massMatrix);

        // Set starting proposal.
        _proposedMomentum = _genMomPropose ? randn_Cholesky(_CholeskyLowerMassMatrix) : randn(_massMatrix);
//        if (_norMom) { _proposedMomentum = normalise(_proposedMomentum); }
        _proposedModel = randn(_prior._means, arma::vec(arma::mat(_prior._covariance).diag()));

        // Set starting model.
        _currentModel = _proposedModel;
        _currentMomentum = _proposedMomentum;

        // Set precomputed quantities.
        _A = _prior._inv_cov_m + _model._g.t() * _data._inv_cov_d * _model._g;
        _bT = _prior._inv_cov_m * _prior._means + _model._g.t() * _data._inv_cov_d * _data._observedData;
        _c = arma::conv_to<double>::from(0.5 * (
                _prior._means.t() * _prior._inv_cov_m * _prior._means +
                _data._observedData.t() * _data._inv_cov_d * _data._observedData
        ));
    };

    void sampler::setStarting(arma::vec &model) {
        _currentModel = model;
        _proposedModel = model;
    }

    void sampler::propose_metropolis() {
        _proposedModel = randn(_prior._means, _prior._covariance.t());
    }

    void sampler::propose_hamilton(int &uturns) {
        /* Draw random prior momenta. */
        _proposedMomentum = _genMomPropose ? randn_Cholesky(_CholeskyLowerMassMatrix) : randn(_massMatrix);
//        if (_norMom) _proposedMomentum = sqrt(_currentMomentum * _currentMomentum) * normalise(_proposedMomentum);
    }

    double sampler::precomp_misfit() {
        return conv_to<double>::from(0.5 * _proposedModel.t() * (_A * _proposedModel) - _bT.t() * _proposedModel + _c);
    }

    arma::vec sampler::precomp_misfitGrad() {
        // Should actually be left multiply, but matrix is symmetric, so skipped that bit.
        return _A * _proposedModel - _bT;
    }

    double sampler::kineticEnergy() {
        return _genMomKinetic ?
               conv_to<double>::from(0.5 * _proposedMomentum.t() * _inverseMassMatrix * _proposedMomentum) :
               conv_to<double>::from(0.5 * _proposedMomentum.t() * _inverseMassMatrixDiagonal * _proposedMomentum);
    }

    double sampler::chi() {
        return precomp_misfit();
    }

    double sampler::energy() {
        return chi() + kineticEnergy();
    }

    void sampler::sample() {

        std::cout << "Inversion of linear model using MCMC sampling." << std::endl;
        std::cout << "Selected method; \033[1;34m" << (_hmc ? "hmc" : "Metropolis-Hastings")
                  << "\033[0m with following options:"
                  << std::endl;
        std::cout << "\t parameters:   \033[1;32m" << _currentModel.size() << "\033[0m" << std::endl;
        std::cout << "\t proposals:    \033[1;32m" << _proposals << "\033[0m" << std::endl;
        std::cout << "\t timestep:     \033[1;32m" << _dt << "\033[0m" << std::endl;
        std::cout << "\t gravity:      \033[1;32m" << _gravity << "\033[0m" << std::endl;

        if (_testBefore) {
            std::cout << "\t - Exploiting conservation of energy by evaluating before propagation" << std::endl;
        }
        std::cout << "\t - Use generalised mass matrix with" << (_genMomPropose ? "" : "out")
                  << " off diagonal entries" << std::endl;
        if (_genMomKinetic) std::cout << "\t - Use generalised momentum for kinetic energy" << std::endl;

        double x = _hmc ? energy() : chi();
        double x_new;
        int accepted = 1;
        int uturns = 0;

        std::ofstream samplesfile;
        samplesfile.open(_outfile);
        samplesfile << _prior._means.size() << " " << _proposals << std::endl;

        write_sample(samplesfile, x);

        std::cout << "[" << std::setw(3) << (int) (100.0 * double(0) / _proposals) << "%] "
                  << std::string(((unsigned long) ((_window.ws_col - 7) * 0 / _proposals)), *"=") <<
                  "\r" << std::flush;

        for (int it = 1; it < _proposals; it++) {

            if (it % 850 == 0) { // Display progress
                std::cout << "[" << std::setw(3) << (int) (100.0 * double(it) / _proposals) << "%] "
                          << std::string(((unsigned long) ((_window.ws_col - 7) * it / _proposals)), *"=") <<
                          "\r" << std::flush;
            }

            if (_hmc) {
                propose_hamilton(uturns);
                if (!_testBefore) {
                    leap_frog(uturns, it == _proposals - 1);
                }
            } else {
                propose_metropolis();
            }

            x_new = _acceptanceFactor * (_hmc ? energy() : chi());

            double result;
            result = x - x_new;
            double result_exponent;
            result_exponent = exp(result);

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
        // Write results
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
            trajectoryfile.open("OUTPUT/trajectory.txt");
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
            angle1 = arma::conv_to<double>::from(_proposedMomentum.t() * (_currentModel - _proposedModel));
            angle2 = arma::conv_to<double>::from(_currentMomentum.t() * (_proposedModel - _currentModel));

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

    arma::vec sampler::precomp_misfitGrad(arma::vec parameters) {
        // Should actually be left multiply, but matrix is symmetric, so skipped that bit.
        return _A * parameters - _bT;
    }
} // namespace hmc
