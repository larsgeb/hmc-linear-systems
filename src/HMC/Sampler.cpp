//
// Created by Lars Gebraad on 18-8-17.
//
#include <randomnumbers.hpp>
#include <cmath>
#include "ForwardModel.hpp"
#include "Data.hpp"
#include "Prior.hpp"
#include "Sampler.hpp"

namespace HMC {
    Sampler::Sampler(Prior &prior, Data &data, ForwardModel &model, GenerateInversionSettings settings) :
            _data(data),
            _prior(prior),
            _model(model) {
        _nt = settings._trajectorySteps;
        _dt = settings._timeStep;
        _gravity = settings._gravity;
        _proposals = settings._proposals;
        _genMomKinetic = settings._genMomKinetic;
        _genMomPropose = settings._genMomPropose;
        _norMom = settings._norMom;
        _testBefore = settings._testBefore;

        /* Initialise random number generator. ----------------------------------------*/
        srand((unsigned int) time(nullptr));

        // Pre-compute mass matrix and other associated quantities
        _massMatrix = _gravity * (_prior._inv_cov_m + ((_model._g.Transpose() * _data._inv_cov_d) * _model._g));
        _inverseMassMatrixDiagonal = AlgebraLib::VectorToDiagonal(_massMatrix.InvertMatrixElements(true).Trace());

        // Prepare mass matrix decomposition and inverse.
        _CholeskyLowerMassMatrix = _massMatrix.CholeskyDecompose();
        AlgebraLib::Matrix InverseCholeskyLowerMassMatrix = _CholeskyLowerMassMatrix.InvertLowerTriangular();
        _inverseMassMatrix = InverseCholeskyLowerMassMatrix.Transpose() * InverseCholeskyLowerMassMatrix;

        // Set starting proposal.
        _proposedMomentum = _genMomPropose ? randn_Cholesky(_CholeskyLowerMassMatrix) : randn(_massMatrix);
        _norMom ? _proposedMomentum = _proposedMomentum.Normalize() : AlgebraLib::Vector();
        _proposedModel = randn(_prior._means, _prior._covariance.Trace());

        // Set starting model.
        _currentModel = _proposedModel;
        _currentMomentum = _proposedMomentum;

        // Set precomputed quantities.
        _A = _prior._inv_cov_m + _model._g.Transpose() * _data._inv_cov_d * _model._g;
        _bT = _prior._inv_cov_m * _prior._means + _model._g.Transpose() * _data._inv_cov_d * _data._observedData;
        _c = 0.5 * (
                _prior._means.Transpose() * _prior._inv_cov_m * _prior._means +
                _data._observedData.Transpose() * _data._inv_cov_d * _data._observedData
        );
    };

    void Sampler::propose_metropolis() {
        _proposedModel = randn(_prior._means, _prior._covariance.Trace());
    }

    void Sampler::propose_hamilton(int &uturns) {
        /* Draw random prior momenta. */
        _proposedMomentum = _genMomPropose ? randn_Cholesky(_CholeskyLowerMassMatrix) : randn(_massMatrix);
        if (_norMom) _proposedMomentum = sqrt(_currentMomentum * _currentMomentum) * _proposedMomentum.Normalize();
    }

    double Sampler::precomp_misfit() {
        return 0.5 * _proposedModel * (_A * _proposedModel) - _bT * _proposedModel + _c;
    }

    AlgebraLib::Vector Sampler::precomp_misfitGrad() {
        // Should actually be left multiply, but matrix is symmetric, so skipped that bit.
        return _A * _proposedModel - _bT;
    }

    double Sampler::kineticEnergy() {
        return _genMomKinetic ?
               0.5 * _proposedMomentum.Transpose() * _inverseMassMatrix * _proposedMomentum :
               0.5 * _proposedMomentum.Transpose() * _inverseMassMatrixDiagonal * _proposedMomentum;
    }

    double Sampler::chi() {
        return precomp_misfit();
    }

    double Sampler::energy() {
        return chi() + kineticEnergy();
    }

    void Sampler::sample(bool hamilton) {
        double x = hamilton ? energy() : chi();
        double x_new;
        int accepted = 1;
        int uturns = 0;

        std::ofstream samplesfile;
        samplesfile.open("OUTPUT/samples.txt");
        samplesfile << _prior._means.size() << " " << _proposals << std::endl;

        write_sample(samplesfile, x);
        for (int it = 1; it < _proposals; it++) {

            if (hamilton) {
                propose_hamilton(uturns);
                if (!_testBefore) {
                    leap_frog(uturns, it == _proposals - 1);
                }
            } else {
                propose_metropolis();
            }

            x_new = (hamilton ? energy() : chi());

            double result;
            result = x - x_new;
            double result_exponent;
            result_exponent = exp(result);

//        if (true) {
            if ((x_new < x) || (result_exponent > randf(0.0, 1.0))) {
//            double Hamiltonian = energy();
//            std::cout<< Hamiltonian;
                if (_testBefore) {
                    leap_frog(uturns, it == _proposals - 1);
                }
//            Hamiltonian = energy();
//            std::cout<< Hamiltonian;
                accepted++;
                x = x_new;
                _currentModel = _proposedModel;
                write_sample(samplesfile, x);
            }
        }
        samplesfile << accepted << std::endl;
        samplesfile.close();

        std::cout << "Number of accepted models: " << accepted << std::endl;
        std::cout << "Number of U-Turn terminations in propagation: " << uturns;
    }

/* Leap-frog integration of Hamilton's equations. ---------------------------------*/
    void Sampler::leap_frog(int &uturns, bool writeTrajectory) {

        // start proposal at current momentum
        _proposedModel = _currentModel;
        // Acts as starting momentum
        _currentMomentum = _proposedMomentum;

        AlgebraLib::Vector misfitGrad;
        double angle1, angle2;

        std::ofstream trajectoryfile;
        if (writeTrajectory) {
            trajectoryfile.open("OUTPUT/trajectory.txt");
            trajectoryfile << _prior._means.size() << " " << _nt << std::endl;
        }

        for (int it = 0; it < _nt; it++) {
            misfitGrad = precomp_misfitGrad();
            _proposedMomentum = _proposedMomentum - 0.5 * _dt * misfitGrad;

            if (writeTrajectory) write_sample(trajectoryfile, chi());
            // Full step in position. Linear algebra does not allow for dividing by diagonal of matrix, hence the loop.
            _proposedModel = _proposedModel + _dt * (
                    (_genMomKinetic ? _inverseMassMatrix : _inverseMassMatrixDiagonal) * _proposedMomentum);
            // Second branch produces unnecessary overhead (lot of zeros).

            misfitGrad = precomp_misfitGrad();
            _proposedMomentum = _proposedMomentum - 0.5 * _dt * misfitGrad;

            /* Check no-U-turn criterion. */
            angle1 = _proposedMomentum * (_currentModel - _proposedModel);
            angle2 = _currentMomentum * (_proposedModel - _currentModel);

            if (angle1 > 0.0 && angle2 > 0.0) {
                uturns++;
                break;
            }
        }
        if (writeTrajectory) trajectoryfile.close();
    }

    void Sampler::write_sample(std::ofstream &outfile, double misfit) {
        for (double j : _proposedModel) {
            outfile << j << "  ";
        }
        outfile << misfit;
        outfile << std::endl;

    }

    AlgebraLib::Vector Sampler::precomp_misfitGrad(AlgebraLib::Vector parameters) {
        // Should actually be left multiply, but matrix is symmetric, so skipped that bit.
        return _A * parameters - _bT;
    }

}
