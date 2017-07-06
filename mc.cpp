#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mc.hpp"
#include "aux.hpp"

/********************************************************************************//**
 * Base class for Monte Carlo sampling.
************************************************************************************/

/* Constructor. -------------------------------------------------------------------*/
mc::mc(int iterations_in, int nt_in, double dt_in, int _Nq, bool verbose_in) {
    /* Set basic parameters. ------------------------------------------------------*/
    nt = nt_in;
    dt = dt_in;
    verbose = verbose_in;
    iterations = iterations_in;
    Nq = _Nq;
    /* Read starting model and prior data. */
    q.read_input("INPUT/parameters_starting_model.txt");
    q_new.read_input("INPUT/parameters_starting_model.txt");

    /* Read observed and compute prior data. */
    dat.read_data("DATA/synthetics.txt");

    q.tExpand(dat, 1, 1.000001);

    zeroth = q.t0U;
    first = q.t1U;

    /* Initialise random number generator. ----------------------------------------*/
    srand(time(0));

    /* Initialise models and momenta. ---------------------------------------------*/
    p = new double[Nq];
    p_new = new double[Nq];

    for (int i = 0; i < Nq; i++) {
        p_new[i] = randn(0.0, sqrt(this->q.iCM[i][i]));
        // Not sure about assigning random moment, now diagonal element of inverse covariance matrix
        q_new.q[i] = randn(q.mean_q[i], q.sigma_q[i]);
        q.q[i] = q_new.q[i];
    }
}

/* Destructor. --------------------------------------------------------------------*/
mc::~mc() {

}


/********************************************************************************//**
 * Proposals.
************************************************************************************/

/* Propose test model based on prior. ---------------------------------------------*/
void mc::propose_metropolis() {
    for (int i = 0; i < Nq; i++) q_new.q[i] = randn(q.mean_q[i], q.sigma_q[i]);
}

/* Proposal based on the solution of Hamilton's equation. -------------------------*/
void mc::propose_hamilton() {
    /* Draw random prior momenta. */
    for (int i = 0; i < Nq; i++) p_new[i] = randn(0.0, sqrt(q.iCM[i][i]));

    /* Integrate Hamilton's equations. */
    leap_frog(false);

}


/********************************************************************************//**
 * Misfits.
************************************************************************************/

/* Misfit for Metropolis Hastings. ------------------------------------------------*/
double mc::chi() {
    double Xi;
    Xi = dat.misfit(&q_new);
    return Xi;
}

/* Misfit for Hamiltonian Monte Carlo. --------------------------------------------*/
double mc::energy() {
    /* Local variables. */
    double H;
    /* Compute energy - model part. */
    H = chi();
    /* Compute energy - momentum part. */
    for (int i = 0; i < Nq; i++) {
        H += 0.5 * p_new[i] * p_new[i] * q_new.iCM[i][i];
    }

    return H;
}

/*=================================================================================*/
/* Right-hand sides of Hamilton equations. ----------------------------------------*/
/*=================================================================================*/

/********************************************************************************//**
 * Miscellaneous.
************************************************************************************/

/* Write a sample to an open file. ------------------------------------------------*/
void mc::write_sample(FILE *pfile, double misfit, int iteration) {
    if (iteration == 0) fprintf(pfile, "%d %d\n", Nq, iterations + 1);

    for (int i = 0; i < Nq; i++) fprintf(pfile, "%lg ", q.q[i]);
    fprintf(pfile, "%lg ", misfit);
    fprintf(pfile, "\n");
}

/* Leap-frog integration of Hamilton's equations. ---------------------------------*/
void mc::leap_frog(bool verbose) {
    /* Local variables and setup. -------------------------------------------------*/

    double *p_half, *p_init, *out;
    double angle1, angle2;
    parameters q_init;

    out = new double[Nq];
    p_half = new double[Nq];
    p_init = new double[Nq];

    FILE *pfile;
    if (verbose) pfile = fopen("OUTPUT/trajectory.txt", "w");

    /* Set initial values. --------------------------------------------------------*/

    q_init = q;
    q_new = q;
    for (int i = 0; i < Nq; i++) {
        p_init[i] = p_new[i];
    }

    /* March forward. -------------------------------------------------------------*/

    if (verbose) fprintf(pfile, "%d %d\n", 2 * Nq, nt);

    for (int it = 0; it < nt; it++) {
        /* Some output. */
        if (verbose) {
            for (int i = 0; i < Nq; i++) fprintf(pfile, "%lg ", q_init.q[i]);
            for (int i = 0; i < Nq; i++) fprintf(pfile, "%lg ", p_init[i]);
            fprintf(pfile, "\n");
        }

        /* First half step in momentum. */
        for (int i = 0; i < Nq; i++) {
            p_half[i] = p_init[i] - 0.5 * dt * q.t1U[i];
            // q.t1U is the gradient approximated up to first order evaluated at the starting model.
        }

        /* Full step in position. */
        for (int i = 0; i < Nq; i++) {
            q_new.q[i] = q_init.q[i] + dt * p_half[i] * q.iCM[i][i];
            // q.t1U is the gradient approximated up to first order evaluated at the starting model.
        }

        /* Second half step in momentum. */
        for (int i = 0; i < Nq; i++) {
            p_new[i] = p_half[i] - 0.5 * dt * q.t1U[i];
        }

        /* Update position and momentum. */
        for (int i = 0; i < Nq; i++) {
            p_init[i] = p_new[i];
            q_init.q[i] = q_new.q[i];
        }

        /* Check no-U-turn criterion. */
        angle1 = 0.0;
        angle2 = 0.0;
        for (int i = 0; i < Nq; i++) {
            angle1 += p_new[i] * (q_new.q[i] - q.q[i]);
            angle2 += p[i] * (q.q[i] - q_new.q[i]);
        }

        if (angle1 < 0.0 && angle2 < 0.0) {
            if (verbose) printf("steps: %d\n", it);
            break;
        }
    }

    /* Clean up. ------------------------------------------------------------------*/

    if (verbose) fclose(pfile);

    delete[] p_half;
    delete[] p_init;
    delete[] out;
}




