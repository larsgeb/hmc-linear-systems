//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_SAMPLER_HPP
#define HMC_LINEAR_SYSTEM_SAMPLER_HPP

#include "prior.hpp"
#include "data.hpp"
#include "posterior.hpp"

#include <sys/ioctl.h>
#include <cstdio>
#include <unistd.h>

namespace hmc {
    struct GenerateInversionSettings {
        double _gravity = 1.0;
        unsigned long int _proposals = 10000;
        double _timeStep = 0.1;
        unsigned long int _trajectorySteps = 10;
        struct winsize window{};
        bool _genMomPropose = true; // Use generalized mass matrix to propose new momenta (true).
        bool _genMomKinetic = true; // Use generalized mass matrix to compute kinetic energy (true).
        bool _norMom = false; // Normalize momentum to previous value to keep constant energy level (true).
        bool _testBefore = true; // Decreases required computation time by order of magnitude, no other influence.
        bool _hamiltonianMonteCarlo = true; // Metropolis Hastings (false) or Hamiltonian Monte Carlo (true).

        GenerateInversionSettings() {
            ioctl(STDOUT_FILENO, TIOCGWINSZ, &window);
            if (window.ws_col < 5)
                window.ws_col = 20;
            if (window.ws_row < 5)
                window.ws_row = 20;
        }

        GenerateInversionSettings &setSamples(unsigned long int samples) { _proposals = samples; };

        GenerateInversionSettings &setTimeStep(double timeStep) { _timeStep = timeStep; };

        GenerateInversionSettings &
        setTrajectorySteps(unsigned long int trajectorySteps) { _trajectorySteps = trajectorySteps; };

        GenerateInversionSettings &setgenMomPropose(bool genMomPropose) { _genMomPropose = genMomPropose; }

        GenerateInversionSettings &setgenMomKinetic(bool genMomKinetic) { _genMomKinetic = genMomKinetic; }

        GenerateInversionSettings &setnorMom(bool norMom) { _norMom = norMom; }

        GenerateInversionSettings &setAcceptBeforeTraj(bool acceptBeforeTraj) { _testBefore = acceptBeforeTraj; }

        GenerateInversionSettings &setGravity(double gravity) { _gravity = gravity; }

        GenerateInversionSettings &sethamiltonianMonteCarlo
                (bool hamiltonianMonteCarlo) { _hamiltonianMonteCarlo = hamiltonianMonteCarlo; }

    };

    class sampler {
    public:
        // Constructors and destructors
        sampler(prior &prior, data &data, ForwardModel &model, GenerateInversionSettings settings);

        void sample();

        algebra_lib::vector precomp_misfitGrad(algebra_lib::vector parameters);

        void setStarting(algebra_lib::vector &model, algebra_lib::vector &momentum);

        algebra_lib::vector _currentModel;
        algebra_lib::vector _proposedModel;
        algebra_lib::vector _currentMomentum;
        algebra_lib::vector _proposedMomentum;

    private:
        // Fields
        prior _prior;
        data _data;
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
        bool _hmc;
        winsize _window;

        algebra_lib::matrix _massMatrix;
        algebra_lib::matrix _CholeskyLowerMassMatrix;
        algebra_lib::matrix _inverseMassMatrix; // needed to write Hamilton's equations in vector form
        algebra_lib::matrix _inverseMassMatrixDiagonal; // needed to write Hamilton's equations in vector form

        // Precomputed misfit function size
        algebra_lib::matrix _A;
        algebra_lib::vector _bT; // Because I haven't coded up the actual left multiplication of vector-matrices
        double _c;

        // Member functions
        void propose_metropolis();

        void propose_hamilton(int &uturns);

        void leap_frog(int &uturns, bool writeTrajectory);

        double chi();

        double energy();

        void write_sample(std::ofstream &outfile, double misfit);

        double precomp_misfit();

        algebra_lib::vector precomp_misfitGrad();

        double kineticEnergy();
    };
}

#endif //HMC_LINEAR_SYSTEM_SAMPLER_HPP
