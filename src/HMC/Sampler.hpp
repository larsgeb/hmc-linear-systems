//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_SAMPLER_HPP
#define HMC_LINEAR_SYSTEM_SAMPLER_HPP

#include "Prior.hpp"
#include "Data.hpp"
#include "Posterior.hpp"

namespace HMC {
    struct GenerateInversionSettings {
        double _gravity = 1.0;
        unsigned long int _proposals = 10000;
        double _timeStep = 0.1;
        unsigned long int _trajectorySteps = 10;
        bool _genMomPropose = true; // Use generalized mass matrix to propose new momenta (true).
        bool _genMomKinetic = true; // Use generalized mass matrix to compute kinetic energy (true).
        bool _norMom = false; // Normalize momentum to previous value to keep constant energy level (true).
        bool _testBefore = true; // Decreases required computation time by order of magnitude, no other influence.
        bool _hamiltonianMonteCarlo = true; // Metropolis Hastings (false) or Hamiltonian Monte Carlo (true).

        GenerateInversionSettings &setSamples(unsigned long int samples) { _proposals =
                                                                                   samples; };
        GenerateInversionSettings &setTimeStep(double _timeStep) { _timeStep = _timeStep; };
        GenerateInversionSettings &
        setTrajectorySteps(unsigned long int _trajectorySteps) { _trajectorySteps = _trajectorySteps; };
        GenerateInversionSettings &setgenMomPropose(bool genMomPropose) { _genMomPropose = genMomPropose; }
        GenerateInversionSettings &setgenMomKinetic(bool genMomKinetic) { _genMomKinetic = genMomKinetic; }
        GenerateInversionSettings &setnorMom(bool norMom) { _norMom = norMom; }
        GenerateInversionSettings &setAcceptBeforeTraj(bool acceptBeforeTraj) { _testBefore = acceptBeforeTraj; }
        GenerateInversionSettings &setGravity(double gravity) { _gravity = gravity; }
        GenerateInversionSettings &sethamiltonianMonteCarlo
                (bool hamiltonianMonteCarlo) { _hamiltonianMonteCarlo = hamiltonianMonteCarlo; }

    };

    class Sampler {
    public:
        // Constructors and destructors
        Sampler(Prior &prior, Data &data, ForwardModel &model, GenerateInversionSettings settings);

        void sample(bool hamilton);

        AlgebraLib::Vector precomp_misfitGrad(AlgebraLib::Vector parameters);

    private:
        // Fields
        Prior _prior;
        Data _data;
        ForwardModel _model;
        Posterior _posterior;
        unsigned long _nt; // Number of time steps for trajectory
        double _dt; // Time step for trajectory
        double _gravity; // Global gravitational constant
        unsigned long _proposals; // Number of iterations for Monte Carlo sampling
        bool _genMomKinetic;
        bool _genMomPropose;
        bool _norMom;
        bool _testBefore;

        AlgebraLib::Vector _currentModel;
        AlgebraLib::Vector _proposedModel;
        AlgebraLib::Vector _currentMomentum;
        AlgebraLib::Vector _proposedMomentum;

        AlgebraLib::Matrix _massMatrix;
        AlgebraLib::Matrix _CholeskyLowerMassMatrix;
        AlgebraLib::Matrix _inverseMassMatrix; // needed to write Hamilton's equations in vector form
        AlgebraLib::Matrix _inverseMassMatrixDiagonal; // needed to write Hamilton's equations in vector form

        // Precomputed misfit function size
        AlgebraLib::Matrix _A;
        AlgebraLib::Vector _bT; // Because I haven't coded up the actual left multiplication of vector-matrices
        double _c;

        // Member functions
        void propose_metropolis();

        void propose_hamilton(int &uturns);

        void leap_frog(int &uturns, bool writeTrajectory);

        double chi();

        double energy();

        void write_sample(std::ofstream &outfile, double misfit);

        double precomp_misfit();

        AlgebraLib::Vector precomp_misfitGrad();

        double kineticEnergy();
    };
}

#endif //HMC_LINEAR_SYSTEM_SAMPLER_HPP
