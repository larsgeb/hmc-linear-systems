#ifndef mc_hpp
#define mc_hpp

#include <stdio.h>
#include <vector>
#include "aux.hpp"

/********************************************************************************//**
 * Base class for Monte Carlo sampling.
************************************************************************************/

class montecarlo {
public:
    /* Regular member variables. -------------------------------------------------*/

    bool m_verbose;              /* Print and write output */

    int m_Nq;                    /**< Model space dimension. */
    int m_iterations;            /**< Number of iterations. */

    int m_nt;                    /**< Maximum number of time steps for numerical integration. */
    double m_dt;                 /**< Time increment for numerical integration. */

    data m_data;                 /* Observed data. */
    data syn;                    /* Synthetics. */

    std::vector<double> m_q_starting;     /**< Start model vector. */

    parameters m_q;              /**< Current model vector. */
    parameters m_q_new;          /**< Test model vector. */

    std::vector<double> p;     /**< Current momentum vector. */
    std::vector<double> p_new; /**< Test momentum vector. */

    // Taylor expansion objects (tensors that grow in dimension with every next order)
    double zeroth;
    std::vector<double> first;

    /* Constructor and destructor. ------------------------------------------------*/
    montecarlo(int in_iterations, int in_nt, double in_dt, int in_Nq,
               bool in_verbose);            /**< Constructor. */
    ~montecarlo();                                                          /**< Destructor. */

    /* Test model proposals. ------------------------------------------------------*/
    void
    propose_metropolis();      /**< Proposal based on prior, to be specified in this function. */
    void propose_hamilton();        /**< Proposal based on Hamiltonian mechanics. */

    /* Misfits. -------------------------------------------------------------------*/
    double chi();       /** Misfit functional for Metropolis Hastings. */
    double energy();    /** Total energy for Hamiltonian Monte Carlo. */

    /** Write a sample to an open file. */
    void write_sample(FILE *pfile,    /**< Pointer to open file. */
                      double misfit,  /**< Misfit value. */
                      int iteration   /**< Current iteration. */
    );

    /** Leap-frog integration of Hamilton's equations with initial positions and momenta. */
    void leap_frog(
            bool output               /**< File output or not. */
    );
};

#endif /* mc_hpp */