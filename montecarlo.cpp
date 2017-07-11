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
//    double Xi;
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
    /* Local variables and setup. -------------------------------------------------*/

//    double angle1, angle2;
//    std::vector<double> initialModel;
//    std::vector<double> out;
//    std::vector<double> p_half;
//    std::vector<double> p_init;
//    std::vector<double> p_start;


    /* Set initial values. --------------------------------------------------------*/
//    initialModel = _currentModel;
/*    for (int i = 0; i < _prior._numberParameters; i++) {
        p_init.push_back(_proposedMomentum[i]);
        p_start.push_back(_proposedMomentum[i]);
        p_half.push_back(0);
    }*/

    _proposedModel = _currentModel;
    std::vector<double> tempMisfit;
    /* March forward. -------------------------------------------------------------*/
    for (int it = 0; it < _nt; it++) {

        tempMisfit = _misfitApproximation.gradient(_proposedModel);
        /* First half step in momentum. */
        for (int i = 0; i < _prior._numberParameters; i++) {
            _proposedMomentum[i] = _proposedMomentum[i] - 0.5 * _dt * tempMisfit[i];
        }
        tempMisfit.clear();

        /* Full step in position. */
        for (int i = 0; i < _prior._numberParameters; i++) {
            _proposedModel[i] = _proposedModel[i] + _dt * _proposedMomentum[i] * (1 / _data._inverseCD[i][i]);
        }

        // Update misfit to new model position
        tempMisfit = _misfitApproximation.gradient(_proposedModel);
        /* Second half step in momentum. */
        for (int i = 0; i < _prior._numberParameters; i++) {
            _proposedMomentum[i] = _proposedMomentum[i] - 0.5 * _dt * tempMisfit[i];
        }
        tempMisfit.clear();

/*        *//* Check no-U-turn criterion. *//*
        angle1 = 0.0;
        angle2 = 0.0;
        for (int i = 0; i < m_Nq; i++) {
            angle1 += p_new[i] * (m_q_new.q[i] - m_q.q[i]);
            angle2 += p_start[i] * (m_q.q[i] - m_q_new.q[i]);
        }

        if (angle1 < 0.0 && angle2 < 0.0) {
            if (verbose) printf("steps: %d\n", it);
            break;
        }*/
    }
}