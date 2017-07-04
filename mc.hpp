
#ifndef mc_hpp
#define mc_hpp

#include <stdio.h>
#include "aux.hpp"

#endif /* mc_hpp */


/********************************************************************************//**
 * Base class for Monte Carlo sampling.
************************************************************************************/

class mc
{
public:
    
    /* Regular member variables. -------------------------------------------------*/
    
    bool verbose;                                   /**< Print and write output, including
                                                        + Hamiltonian trajectory written to OUTPUT/trajectory.txt,
                                                        + number of time steps along the Hamiltonian trajectory. 
                                                        + pre-computed scalar, vector and matrix. */
    
    int Nq;                                         /**< Model space dimension. */
    int iterations;                                 /**< Number of iterations. */
    
    int nt;                                         /**< Maximum number of time steps for numerical integration. */
    double dt;                                      /**< Time increment for numerical integration. */
    
    data dat;                                       /* Observed data. */
    data syn;                                       /* Synthetics. */
    
    parameters q;                                   /**< Current model vector. */
    parameters q_new;                               /**< Test model vector. */
    
    double *p;                                      /**< Current momentum vector. */
    double *p_new;                                  /**< Test momentum vector. */
    
    double *m;                                      /**< Mass matrix. */
    
    /* Constructor and destructor. ------------------------------------------------*/
    
    mc(int iterations, int nt, double dt, bool verbose);            /**< Constructor. */
    ~mc();                                                          /**< Destructor. */
    
    /* Precomputed matrices for fast misfit & derivative evaluation. --------------*/
    
    double **A;
    double *B;
    double C;
    
    /* Test model proposals. ------------------------------------------------------*/
    
    void propose_metropolis();      /**< Proposal based on prior, to be specified in this function. */
    void propose_hamilton();        /**< Proposal based on Hamiltonian mechanics. */
    
    /* Misfits. -------------------------------------------------------------------*/
    
    double chi();       /** Misfit functional for Metropolis Hastings. */
    double energy();    /** Total energy for Hamiltonian Monte Carlo. */
    
    /* Right-hand sides of Hamilton equations. ------------------------------------*/
    
    /** dH/dp. Needs to be adjusted as needed. */
    void dHdp(
              double *p,        /**< Momentum vector. */
              double *rhs       /**< Output right-hand side, i.e. dH/dp. Must be allocated. */
            );
    
    /** dH/dq. Needs to be adjusted as needed. */
    void dHdq(
              parameters &q,        /**< Model (position) vector. */
              double *rhs           /**< Output right-hand side, i.e. dH/dq. Must be allocated. */
            );
    
    /* Miscellaneous. -------------------------------------------------------------*/
    
    /** Write a sample to an open file. */
    void write_sample(FILE *pfile,                  /**< Pointer to open file. */
                      double misfit,                /**< Misfit value. */
                      int iteration                 /**< Current iteration. */
                      );
    
    /** Leap-frog integration of Hamilton's equations with initial positions and momenta. */
    void leap_frog(
                    bool output                                 /**< File output or not. */
                    );
};

