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
#include <stdio.h>
#include <iostream>

montecarlo::montecarlo(std::vector<double> startModel, prior &in_prior, data &in_data, posterior &in_posterior, int in_nt,
                       double
                       in_dt, int in_iterations) {
    _prior = in_prior;
    _data = in_data;
    _posterior = in_posterior;
    _misfitApproximation = taylorExpansion(startModel, 0.01, _prior, _data, _posterior);
    _nt = in_nt;
    _dt = in_dt;
    _iterations = in_iterations;

    /* Initialise random number generator. ----------------------------------------*/
    srand((unsigned int) time(0));

    _currentModel = startModel;
    _proposedModel = startModel;

    // These shouldn't be initialized yet, but better safe than sorry when using std::vector.push_back().
    _proposedMomentum.clear();
    _currentMomentum.clear();

    // Assigning a random moment to the momentum vectors
    for (int i = 0; i < _prior._numberParameters; i++) {
        // But this momentum assignment is as of yet only the diagonal.
        _proposedMomentum.push_back(randn(0.0, sqrt(_prior._massMatrix[i][i])));
        _proposedModel[i] = randn(_prior._mean[i], _prior._std[i]);
        _currentModel[i] = _proposedModel[i];
        _currentMomentum.push_back(_proposedMomentum[i]);
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
        _proposedMomentum[i] = (randn(0.0, 1 / sqrt(float(_prior._massMatrix[i][i])))); // only diagonal implemented!
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

void montecarlo::sample(bool hamilton) {
    double x = (hamilton ? energy() : chi());
    double x_new;
    int accepted = 0;

    FILE *pfile;
    pfile = fopen("OUTPUT/samples.txt", "w");
    write_sample(pfile, x, 0);
    for (int it = 1; it < _iterations; it++) {
        hamilton ? propose_hamilton() : propose_metropolis();
        x_new = (hamilton ? energy() : chi());

        if ((x_new < x) || (exp(x - x_new) > randf(0.0, 1.0))) {
            accepted++;
            x = x_new;
            _currentModel = _proposedModel;
            write_sample(pfile, x, it);
            if (accepted % 50 == 0) _misfitApproximation.updateExpansion(_currentModel);
        }
    }
    fprintf(pfile, "%i ", accepted);
    fprintf(pfile, "\n");
    fclose(pfile);
    std::cout << "Number of accepted models: " << accepted;
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
            _proposedModel[i] = _proposedModel[i] + _dt * _proposedMomentum[i];
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

void montecarlo::write_sample(FILE *pfile, double misfit, int iteration) {
    if (iteration == 0) fprintf(pfile, "%d %d\n", (int) _prior._numberParameters, _iterations);

    for (int i = 0; i < (int) _prior._numberParameters; i++) fprintf(pfile, "%lg ", _proposedModel[i]);
    fprintf(pfile, "%lg ", misfit);
    fprintf(pfile, "\n");
}