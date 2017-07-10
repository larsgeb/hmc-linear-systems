#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "montecarlo.hpp"
#include "aux.hpp"

/********************************************************************************//**
 * Base class for Monte Carlo sampling.
************************************************************************************/

/* Constructor. -------------------------------------------------------------------*/
montecarlo::montecarlo(int in_iterations, int in_nt, double in_dt, int in_Nq,
                       bool in_verbose) {
    /* Set basic parameters. ------------------------------------------------------*/
    m_nt = in_nt;
    m_dt = in_dt;
    m_verbose = in_verbose;
    m_iterations = in_iterations;
    m_Nq = in_Nq;

    /* Read starting model and prior data. */
    m_q.read_input("INPUT/parameters_starting_model.txt");
    m_q_new.read_input("INPUT/parameters_starting_model.txt");

    m_q_starting = m_q.q;

    /* Read observed data. */
    m_data.read_data("DATA/synthetics.txt");

    /* Do a Taylor expansion for fast misfit and fast leap-frog propagation */
    m_q.tExpand(m_data, 1, 1.000001);

//    zeroth = m_q.t0U;
//    first = m_q.t1U;

    /* Initialise random number generator. ----------------------------------------*/
    srand(time(0));

    /* Not sure about assigning random moment, now diagonal element of inverse
     * covariance matrix. Correlation could be implemented. The physical meaning
     * (within Hamiltonian Mechanics) of this is not immediately clear to me.
     * Some kind of connection between particles. */
    for (int i = 0; i < m_Nq; i++) {
        // But this momentum assignment is as of yet only the diagonal.
        p_new.push_back(randn(0.0, sqrt(this->m_q.iCM[i][i])));
        // Assigning new model vectors based on prior (For a random start model
        // within prior)
        m_q_new.q[i] = randn(m_q.mean_q[i], m_q.sigma_q[i]);
        m_q.q[i] = m_q_new.q[i];
    }
}

/* Destructor. --------------------------------------------------------------------*/
montecarlo::~montecarlo() {

}

/********************************************************************************//**
 * Proposals.
************************************************************************************/

/* Propose test model based on prior. ---------------------------------------------*/
void montecarlo::propose_metropolis() {
    for (int i = 0; i < m_Nq; i++)
        m_q_new.q[i] = randn(m_q.mean_q[i], m_q.sigma_q[i]);
}

/* Proposal based on the solution of Hamilton's equation. -------------------------*/
void montecarlo::propose_hamilton() {
    /* Draw random prior momenta. */
    for (int i = 0; i < m_Nq; i++)
        p_new[i] = randn(0.0, sqrt(m_q.iCM[i][i])); // Only with diagonal
    /* Integrate Hamilton's equations. */
    leap_frog(false);
}

/********************************************************************************//**
 * Misfits.
************************************************************************************/

/* Misfit for Metropolis Hastings. ------------------------------------------------*/
double montecarlo::chi() {
//    double Xi;
//    Xi = m_data.misfit(&m_q_new);
//    return Xi;
    return m_data.misfit(&m_q_new);
}

/* Misfit for Hamiltonian Monte Carlo. --------------------------------------------*/
double montecarlo::energy() {
    /* Local variables. */
    double H;
    /* Compute energy - model part. */
    H = chi();
    /* Compute energy - momentum part. */
    for (int i = 0; i < m_Nq; i++) {
        H += 0.5 * p_new[i] * p_new[i] * m_q_new.iCM[i][i];
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
void montecarlo::write_sample(FILE *pfile, double misfit, int iteration) {
    if (iteration == 0) fprintf(pfile, "%d %d\n", m_Nq, m_iterations + 1);

    for (int i = 0; i < m_Nq; i++) fprintf(pfile, "%lg ", m_q.q[i]);
    fprintf(pfile, "%lg ", misfit);
    fprintf(pfile, "\n");
}

/* Leap-frog integration of Hamilton's equations. ---------------------------------*/
void montecarlo::leap_frog(bool verbose) {
    /* Local variables and setup. -------------------------------------------------*/

    double angle1, angle2;
    parameters q_init;

    std::vector<double> out;
    std::vector<double> p_half;
    std::vector<double> p_init;
    std::vector<double> p_start;

    FILE *pfile;
    if (verbose) pfile = fopen("OUTPUT/trajectory.txt", "w");

    /* Set initial values. --------------------------------------------------------*/
    q_init = m_q;
    m_q_new = m_q;
    for (int i = 0; i < m_Nq; i++) {
        p_init.push_back(p_new[i]);
        p_start.push_back(p_new[i]);
        p_half.push_back(0);
    }

    /* March forward. -------------------------------------------------------------*/
    if (verbose) fprintf(pfile, "%d %d\n", 2 * m_Nq, m_nt);

    for (int it = 0; it < m_nt; it++) {
        /* Some output. */
        if (verbose) {
            for (int i = 0; i < m_Nq; i++) fprintf(pfile, "%lg ", q_init.q[i]);
            for (int i = 0; i < m_Nq; i++) fprintf(pfile, "%lg ", p_init[i]);
            fprintf(pfile, "\n");
        }

        /* First half step in momentum. */
        for (int i = 0; i < m_Nq; i++) {
            double secDer;
            secDer = 0;
            std::vector<double> qDiff;
            qDiff = VectorDifference(m_q_starting,m_q.q);
            for (int j = 0; j < m_Nq; j++){
                secDer += m_q.t2U[i][j]*2*(qDiff[j]);
            }
            qDiff.clear();
            p_half[i] = p_init[i] + 0.5 * m_dt * m_q.t1U[i];
            // m_q.t1U is the gradient approximated up to first order evaluated at the starting model.
        }

        /* Full step in position. */
        for (int i = 0; i < m_Nq; i++) {
            m_q_new.q[i] = q_init.q[i] + m_dt * p_half[i] * m_q.iCM[i][i];
            // m_q.t1U is the gradient approximated up to first order evaluated at the starting model.
        }

        /* Second half step in momentum. */
        for (int i = 0; i < m_Nq; i++) {
            p_new[i] = p_half[i] + 0.5 * m_dt * m_q.t1U[i];
        }

        /* Update position and momentum. */
        for (int i = 0; i < m_Nq; i++) {
            p_init[i] = p_new[i];
            q_init.q[i] = m_q_new.q[i];
        }

        /* Check no-U-turn criterion. */
        angle1 = 0.0;
        angle2 = 0.0;
        for (int i = 0; i < m_Nq; i++) {
            angle1 += p_new[i] * (m_q_new.q[i] - m_q.q[i]);
            angle2 += p_start[i] * (m_q.q[i] - m_q_new.q[i]);
        }

        if (angle1 < 0.0 && angle2 < 0.0) {
            if (verbose) printf("steps: %d\n", it);
            break;
        }
    }

    /* Clean up. ------------------------------------------------------------------*/

    if (verbose) fclose(pfile);

    p_half.clear();
    p_init.clear();
    out.clear();
}





