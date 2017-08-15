//
// Created by Lars Gebraad on 7/11/17.
//

#ifndef HMC_VSP_MONTECARLO_HPP
#define HMC_VSP_MONTECARLO_HPP

/*! \mainpage Hamiltonian Monte Carlo Documentation
 *
 * This set of mainly C++ files allow one to do probabilistic inversion of (Linear) models.
 * What this means is that it takes a function, or an algorithm, accompanied by (synthetic)
 * data and tries to invert for the input parameters.
 *
 * It does so in a probabilistic sense; there's is not one solution, only the solution. This
 * solution contains the probabilities for all possible combinations of parameters. This allows
 * for heavily underdetermined problems to still have a meaningful solution without the need
 * for stabilization.
 *
 * Typically (which might have been some decades ago, as methods have advanced), one would
 * sample our so-called posterior distrubtion by randomly proposing models (based on some prior
 * information) and calculate the misfit. This process is done many times, and each time a
 * model satisfies some criterion, it increases the likelihood of a solution in that
 * 'neighborhood' (since one can not infenitessimally explore the model space, we avarage over
 * cells). A typical algorithm for this is Metropolis-Hastings.
 *
 * However, if one has a very poor idea of where the parameters would be, the proposals will
 * mostly be rejeceted. This is where Hamiltonian Monte Carlo comes in. It not only uses the
 * prior information, but also the local curvature of the misfit functional (typically a scalar
 * function on many dimensions, \f$ \chi = f(x_1,x_2,x_3, \ldots , x_n) \f$) to direct the new
 * proposal towards a likeli model. This propagation takes no history of proposals, so it is
 * called a Markov Chain Monte Carlo (MCMC) method.
 *
 * For a more detailed explanation, I do strongly recommend taking a look at
 * <A href="https://larsgeb.github.io">larsgeb.github.io</A>
 * where a full technical document is present, with more mathematical rigour. The purpose of
 * this site is merely the documentation and clarification of code in a nice online format.
 *
 * I recommend taking a first look at the classes in montecarlo.hpp and auxiliary.hpp, which
 * define all the constructs for the theoretical parts. The rest is mostly tools for linear
 * algebra or random number generators. Interesting, but not essential.
 *
 */

#include <libio.h>
#include <fstream>

class montecarlo {
public:
    // Constructors and destructors
    montecarlo(prior &in_prior, data &in_data, forwardModel in_model, int in_nt, double in_dt, int in_iterations,
               bool useGeneralisedMomentumPropose, bool useGeneralisedMomentumKinetic, bool normalizeMomentum, bool
               evaluateHamiltonianBeforeLeap);

//    montecarlo(prior &in_prior, data &in_data, posterior &in_posterior, forwardModel &in_model, int in_nt, double in_dt,
//               int in_iterations, bool useGeneralisedMomentumPropose, bool useGeneralisedMomentumKinetic, bool normalizeMomentum);


    ~montecarlo();

    void sample(bool hamilton);

    std::vector<double> precomp_misfitGrad(std::vector<double> parameters);

private:
    // Fields
    prior _prior;
    data _data;
    forwardModel _model;
    posterior _posterior;
    int _nt; // Number of time steps for trajectory
    double _dt; // Time step for trajectory
    int _iterations; // Number of iterations for Monte Carlo sampling
    bool _useGeneralisedMomentumKinetic;
    bool _useGeneralisedMomentumPropose;
    bool _normalizeMomentum;
    bool _evaluateHamiltonianBeforeLeap;

    std::vector<double> _currentModel;
    std::vector<double> _proposedModel;
    std::vector<double> _currentMomentum;
    std::vector<double> _proposedMomentum;

    std::vector<std::vector<double>> _massMatrix;
    std::vector<std::vector<double>> _CholeskyLowerMassMatrix;
    std::vector<std::vector<double>> _inverseMassMatrix; // needed to write Hamilton's equations in vector form
    std::vector<std::vector<double>> _inverseMassMatrixDiagonal; // needed to write Hamilton's equations in vector form

    // Precomputed misfit function elements
    std::vector<std::vector<double>> _A;
    std::vector<double> _bT; // Because I haven't coded up the actual left multiplication of vector-matrices
    double _c;

    // Member functions
    void propose_metropolis();

    void propose_hamilton(int &uturns);

    void leap_frog(int &uturns, bool writeTrajectory);

    double chi();

    double energy();

    void write_sample(std::ofstream &outfile, double misfit);

    double precomp_misfit();

    std::vector<double> precomp_misfitGrad();

    double kineticEnergy();
};

#endif //HMC_VSP_MONTECARLO_HPP
