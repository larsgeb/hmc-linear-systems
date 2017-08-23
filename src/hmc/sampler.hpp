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

using namespace algebra_lib;

namespace hmc {
    struct GenerateInversionSettings {
        double _gravity = 1.0;
        unsigned long int _proposals = 10000;
        double _timeStep = 0.1;
        unsigned long int _trajectorySteps = 10;
        struct winsize window{};
        char *_outfile = "samples.txt";
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

        GenerateInversionSettings &setOutfile(char *outputFile) { _outfile = outputFile; }

        GenerateInversionSettings &setGravity(double gravity) { _gravity = gravity; }

        GenerateInversionSettings &sethamiltonianMonteCarlo
                (bool hamiltonianMonteCarlo) { _hamiltonianMonteCarlo = hamiltonianMonteCarlo; }

    };

    class sampler {
    public:
        // Constructors and destructors
        sampler(prior &prior, data &data, forward_model &model, GenerateInversionSettings settings);

        void sample();

        vector precomp_misfitGrad(vector parameters);

        void setStarting(vector &model);

        vector _currentModel;
        vector _proposedModel;
        vector _currentMomentum;
        vector _proposedMomentum;

    private:
        // Fields
        prior _prior;
        data _data;
        forward_model _model;
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
        char *_outfile;
        winsize _window;

        matrix _massMatrix;
        matrix _CholeskyLowerMassMatrix;
        matrix _inverseMassMatrix; // needed to write Hamilton's equations in vector form
        matrix _inverseMassMatrixDiagonal; // needed to write Hamilton's equations in vector form

        // Precomputed misfit function size
        matrix _A;
        vector _bT; // Because I haven't coded up the actual left multiplication of vector-matrices
        double _c;

        // Member functions
        void propose_metropolis();

        void propose_hamilton(int &uturns);

        void leap_frog(int &uturns, bool writeTrajectory);

        double chi();

        double energy();

        void write_sample(std::ofstream &outfile, double misfit);

        double precomp_misfit();

        vector precomp_misfitGrad();

        double kineticEnergy();


    };
}

#endif //HMC_LINEAR_SYSTEM_SAMPLER_HPP
