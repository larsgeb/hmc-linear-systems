//
// Created by Lars Gebraad on 7/11/17.
//

#include <ctime>
#include <utility>
#include <vector>
#include <cmath>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include <cstdlib>
#include "randomnumbers.hpp"
#include "linearalgebra.hpp"
#include <cstdio>
#include <iostream>

montecarlo::montecarlo(prior &in_prior, data &in_data, forwardModel in_model, int in_nt, double in_dt, int in_iterations,
                       bool useGeneralisedMomentumPropose, bool useGeneralisedMomentumKinetic, bool normalizeMomentum, bool
                       evaluateHamiltonianBeforeLeap) {
    _prior = in_prior;
    _data = in_data;
    _model = std::move(in_model);
    _nt = in_nt;
    _dt = in_dt;
    _iterations = in_iterations;
    _useGeneralisedMomentumKinetic = useGeneralisedMomentumKinetic;
    _useGeneralisedMomentumPropose = useGeneralisedMomentumPropose;
    _normalizeMomentum = normalizeMomentum;
    _evaluateHamiltonianBeforeLeap = evaluateHamiltonianBeforeLeap;

    /* Initialise random number generator. ----------------------------------------*/
    srand((unsigned int) time(nullptr));

    // Pre-compute mass matrix and other associated quantities
    _massMatrix = _prior._inverseCovarianceMatrix +
                  (TransposeMatrix(_model._designMatrix) * _data._inverseCD) * _model._designMatrix;
    _inverseMassMatrixDiagonal = VectorToDiagonal(MatrixTrace(InvertMatrixElements(_massMatrix)));

    // This is where the magic happens
    _CholeskyLowerMassMatrix = CholeskyDecompose(_massMatrix);
    std::vector<std::vector<double>> InverseCholeskyLowerMassMatrix = InvertLowerTriangular(_CholeskyLowerMassMatrix);
    _inverseMassMatrix = TransposeMatrix(InverseCholeskyLowerMassMatrix) * InverseCholeskyLowerMassMatrix;

    _proposedMomentum = _useGeneralisedMomentumPropose ?
                        randn_Cholesky(_CholeskyLowerMassMatrix) :
                        randn(_massMatrix);
    _proposedMomentum = _normalizeMomentum ? NormalizeVector(_proposedMomentum) : _proposedMomentum;
    _proposedModel = randn(_prior._mean, _prior._std);
    _currentModel = _proposedModel;
    _currentMomentum = _proposedMomentum;

    _A = _massMatrix;
    _bT = (_prior._inverseCovarianceMatrix * _prior._mean) +
          (TransposeMatrix(_model._designMatrix) * _data._inverseCD) * _data._observedData;
    _c = 0.5 * (_prior._mean * (_prior._inverseCovarianceMatrix * _prior._mean) +
                _data._observedData * (_data._inverseCD * _data._observedData));
};

//montecarlo::montecarlo(prior &in_prior, data &in_data, posterior &in_posterior, forwardModel &in_model, int in_nt,
//                       double in_dt, int in_iterations, bool useGeneralisedMomentumPropose,
//                       bool useGeneralisedMomentumKinetic, bool normalizeMomentum) :
//        montecarlo::montecarlo(in_prior, in_data, std::move(in_model), in_nt, in_dt, in_iterations,
//                               useGeneralisedMomentumPropose, useGeneralisedMomentumKinetic, normalizeMomentum) {
//    _posterior = in_posterior;
//}

montecarlo::~montecarlo() = default;

void montecarlo::propose_metropolis() {
    for (int i = 0; i < _prior._numberParameters; i++)
        _proposedModel[i] = randn(_prior._mean[i], _prior._std[i]);
}

void montecarlo::propose_hamilton(int &uturns) {
    /* Draw random prior momenta. */
    double normalizationFactor = sqrt(_proposedMomentum * _proposedMomentum);
    _proposedMomentum = _useGeneralisedMomentumPropose ?
                        randn_Cholesky(_CholeskyLowerMassMatrix) :
                        randn(_massMatrix);
    _proposedMomentum = _normalizeMomentum ? normalizationFactor * NormalizeVector(_proposedMomentum) : _proposedMomentum;
}

double montecarlo::precomp_misfit() {
    return 0.5 * _proposedModel * (_A * _proposedModel) - _bT * _proposedModel + _c;
}

double montecarlo::kineticEnergy() {
    return _useGeneralisedMomentumKinetic ?
           0.5 * (_proposedMomentum * (_inverseMassMatrix * _proposedMomentum)) :
           0.5 * _proposedMomentum * (_inverseMassMatrixDiagonal * _proposedMomentum);
}

std::vector<double> montecarlo::precomp_misfitGrad() {
    // Should actually be left multiply, but matrix is symmetric, so skipped that bit.
    return _A * _proposedModel - _bT;
}

double montecarlo::chi() {
//    return _posterior.misfit(_proposedModel, _prior, _data, _model);
    return precomp_misfit();
}

double montecarlo::energy() {
    double H = chi();
    H += kineticEnergy();
    return H;
}

void montecarlo::sample(bool hamilton) {
    double x = hamilton ? energy() : chi();
    double x_new;
    int accepted = 1;
    int uturns = 0;

    std::ofstream samplesfile;
    samplesfile.open("OUTPUT/samples.txt");
    samplesfile << _prior._numberParameters << " " << _iterations << std::endl;

    write_sample(samplesfile, x);
    for (int it = 1; it < _iterations; it++) {

        if (hamilton) {
            propose_hamilton(uturns);
            if (!_evaluateHamiltonianBeforeLeap) {
                leap_frog(uturns, it == _iterations - 1);
            }
        } else {
            propose_metropolis();
        }

        x_new = (hamilton ? energy() : chi());

        double result;
        result = x - x_new;
        double result_exponent;
        result_exponent = exp(result);

        if (true) {
//        if ((x_new < x) || (result_exponent > randf(0.0, 1.0))) {
//            double Hamiltonian = energy();
//            std::cout<< Hamiltonian;
            if (_evaluateHamiltonianBeforeLeap) {
                leap_frog(uturns, it == _iterations - 1);
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
void montecarlo::leap_frog(int &uturns, bool writeTrajectory) {

    // start proposal at current momentum
    _proposedModel = _currentModel;
    // Acts as starting momentum
    _currentMomentum = _proposedMomentum;

    std::vector<double> misfitGrad;
    double angle1, angle2;

    std::ofstream trajectoryfile;
    if (writeTrajectory) {
        trajectoryfile.open("OUTPUT/trajectory.txt");
        trajectoryfile << _prior._numberParameters << " " << _nt << std::endl;
    }

    for (int it = 0; it < _nt; it++) {
        misfitGrad = precomp_misfitGrad();
        _proposedMomentum = _proposedMomentum - 0.5 * _dt * misfitGrad;
        misfitGrad.clear();

        if (writeTrajectory) write_sample(trajectoryfile, chi());
        // Full step in position. Linear algebra does not allow for dividing by diagonal of matrix, hence the loop.
        _proposedModel = _proposedModel + _dt * (
                (_useGeneralisedMomentumKinetic ? _inverseMassMatrix : _inverseMassMatrixDiagonal) * _proposedMomentum);
        // Second branch produces unnecessary overhead (lot of zeros).

        misfitGrad = precomp_misfitGrad();
        _proposedMomentum = _proposedMomentum - 0.5 * _dt * misfitGrad;
        misfitGrad.clear();

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

void montecarlo::write_sample(std::ofstream &outfile, double misfit) {
    for (double j : _proposedModel) {
        outfile << j << "  ";
    }
    outfile << misfit;
    outfile << std::endl;

}
