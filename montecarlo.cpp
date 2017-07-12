//
// Created by Lars Gebraad on 7/11/17.
//

#include <time.h>
#include <vector>
#include <math.h>
#include "auxiliary.hpp"
#include "montecarlo.hpp"
#include <stdlib.h>
#include "randomnumbers.hpp"
#include "linearalgebra.hpp"

montecarlo::montecarlo(std::vector<double> startModel, prior &in_prior, data &in_data, posterior &in_posterior, int in_nt,
                       double
                       in_dt, int in_iterations) {
    _prior = in_prior;
    _data = in_data;
    _posterior = in_posterior;
    _misfitApproximation = taylorExpansion(startModel, 0.001, _prior, _data, _posterior);
    _nt = in_nt;
    _dt = in_dt;
    _iterations = in_iterations;

    /* Initialise random number generator. ----------------------------------------*/
    srand((unsigned int) time(0));

    _currentModel = startModel;
    _proposedModel = startModel;

    // Assigning a random moment to the momentum vectors
    for (int i = 0; i < _prior._numberParameters; i++) {
        // But this momentum assignment is as of yet only the diagonal.
        _proposedMomentum.push_back(randn(0.0, sqrt(_prior._massMatrix[i][i])));
        _proposedModel[i] = randn(_prior._mean[i], _prior._std[i]);
//        _currentModel[i] = _proposedModel[i]; Remnant of Andreas' code. Just to assign a random start model. This is
// specified in the current code, as it allows some control over the Taylor Expansion location.
    }
}

montecarlo::~montecarlo() {

}

/* Propose test model based on prior. ---------------------------------------------*/
void montecarlo::propose_metropolis() {
    for (int i = 0; i < _prior._numberParameters; i++)
        _proposedModel[i] = randn(_prior._mean[i], _prior._std[i]);
}

/* Proposal based on the solution of Hamilton's equation. -------------------------*/
void montecarlo::propose_hamilton() {
    /* Draw random prior momenta. */
    for (int i = 0; i < _prior._numberParameters; i++)
        _proposedMomentum[i] = (randn(0.0, sqrt(_prior._massMatrix[i][i]))); // only diagonal implemented!
    /* Integrate Hamilton's equations. */
    leap_frog();
}

double montecarlo::chi() {
    return _posterior.misfit(_proposedModel, _prior, _data);
}

/* Misfit for Hamiltonian Monte Carlo. --------------------------------------------*/
double montecarlo::energy() {
    double H = chi();
    for (int i = 0; i < _prior._numberParameters; i++) {
        H += 0.5 * _proposedMomentum[i] * _proposedMomentum[i] * _prior._massMatrix[i][i];
    }
    return H;
}

/* Leap-frog integration of Hamilton's equations. ---------------------------------*/
void montecarlo::leap_frog() {

    _proposedModel = _currentModel;
    std::vector<double> tempMisfitGrad;
    /* March forward. -------------------------------------------------------------*/
    for (int it = 0; it < _nt; it++) {

        tempMisfitGrad = _misfitApproximation.gradient(_proposedModel);
        /* First half step in momentum. */
        for (int i = 0; i < _prior._numberParameters; i++) {
            _proposedMomentum[i] = _proposedMomentum[i] - 0.5 * _dt * tempMisfitGrad[i];
        }
        tempMisfitGrad.clear();

        /* Full step in position. */
        for (int i = 0; i < _prior._numberParameters; i++) {
            _proposedModel[i] = _proposedModel[i] + _dt * _proposedMomentum[i] * _prior._massMatrix[i][i];
        }

        // Update misfit to new model position
        tempMisfitGrad = _misfitApproximation.gradient(_proposedModel);
        /* Second half step in momentum. */
        for (int i = 0; i < _prior._numberParameters; i++) {
            _proposedMomentum[i] = _proposedMomentum[i] - 0.5 * _dt * tempMisfitGrad[i];
        }
        tempMisfitGrad.clear();

        // TODO, no U-Turn criterion

    }
}