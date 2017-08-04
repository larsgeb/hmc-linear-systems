//
// Created by Lars Gebraad on 7/11/17.
//

#ifndef HMC_VSP_MONTECARLO_HPP
#define HMC_VSP_MONTECARLO_HPP

#include <libio.h>
#include <fstream>

class montecarlo {
public:
    // Constructors and destructors
    montecarlo(prior &in_prior, data &in_data, forwardModel in_model, int in_nt, double in_dt, int in_iterations,
               bool useGeneralisedMomentumPropose, bool useGeneralisedMomentumKinetic);

    montecarlo(prior &in_prior, data &in_data, posterior &in_posterior, forwardModel &in_model, int in_nt, double in_dt,
               int in_iterations, bool useGeneralisedMomentumPropose, bool useGeneralisedMomentumKinetic);


    ~montecarlo();

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

    void propose_hamilton(int &uturns, bool writeTrajectory);

    void leap_frog(int &uturns, bool writeTrajectory);

    double chi();

    double energy();

    void sample(bool hamilton);

    void write_sample(std::ofstream &outfile, double misfit);

    double precomp_misfit();

    std::vector<double> precomp_misfitGrad();

    double kineticEnergy();
};

#endif //HMC_VSP_MONTECARLO_HPP
