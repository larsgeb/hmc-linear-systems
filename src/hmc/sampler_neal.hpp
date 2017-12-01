//
// Created by Lars Gebraad on 1-12-17.
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
    struct InversionSettings {
        const double PI = 3.14159265358979323846264338327;
        double _gravity = 1.0;
        double _temperature = 1.0;
        unsigned long int _proposals = 10000;
        double _timeStep = 0.1;
        double _acceptanceFactor = 1.0;
        unsigned long int _trajectorySteps = 10;
        struct winsize _window{};
        char *_outfile = const_cast<char *>("samples.txt");
        bool _genMomPropose = true; // Use generalized mass matrix to propose new momenta (true).
        bool _genMomKinetic = true; // Use generalized mass matrix to compute kinetic energy (true).
//        bool _norMom = false; // Normalize momentum to previous value to keep constant energy level (true).
        bool _testBefore = true; // Decreases required computation time by order of magnitude, no other influence.
        bool _ergodic = false;  // Randomizes trajectory length and step size
        bool _hamiltonianMonteCarlo = true; // Metropolis Hastings (false) or Hamiltonian Monte Carlo (true).

        InversionSettings() {
            ioctl(STDOUT_FILENO, TIOCGWINSZ, &_window);
            if (_window.ws_col < 5)
                _window.ws_col = 20;
            if (_window.ws_row < 5)
                _window.ws_row = 20;
        }

        InversionSettings &setSamples(unsigned long int samples) { _proposals = samples; };

        InversionSettings &setTimeStep(double timeStep) { _timeStep = timeStep; };

        InversionSettings &setTimeStepFromGrav_nSteps() {
            _timeStep = 1.0 / (_trajectorySteps * sqrt(1.0 / _gravity) / (2 * PI));
        };

        InversionSettings &
        setAcceptanceFactor(double acceptanceFactor) { _acceptanceFactor = acceptanceFactor; };

        InversionSettings &
        setTrajectorySteps(unsigned long int trajectorySteps) { _trajectorySteps = trajectorySteps; };

        InversionSettings &setGenMomPropose(bool genMomPropose) { _genMomPropose = genMomPropose; }

        InversionSettings &setGenMomKinetic(bool genMomKinetic) { _genMomKinetic = genMomKinetic; }

        InversionSettings &setErgodicity(bool ergodic) { _ergodic = ergodic; }

//        InversionSettings &setnorMom(bool norMom) { _norMom = norMom; }

        InversionSettings &setAcceptBeforeTraj(bool acceptBeforeTraj) { _testBefore = acceptBeforeTraj; }

        InversionSettings &setOutfile(char *outputFile) { _outfile = outputFile; }

        InversionSettings &setGravity(double gravity) { _gravity = gravity; }
        InversionSettings &setTemperature(double temperature) { _temperature = temperature; }

        InversionSettings &setHamiltonianMonteCarlo
                (bool hamiltonianMonteCarlo) { _hamiltonianMonteCarlo = hamiltonianMonteCarlo; }
    };

    class sampler {
    public:
        // Constructors and destructors
        sampler(prior &prior, data &data, forward_model &model, InversionSettings settings);

        void sample();

        arma::vec precomp_misfitGrad(arma::vec parameters);

        void setStarting(arma::vec &model);

        arma::vec _currentModel;
        arma::vec _proposedModel;
        arma::vec _currentMomentum;
        arma::vec _proposedMomentum;

    private:
        // Fields
        prior _prior;
        data _data;
        forward_model _model;
        Posterior _posterior;
        unsigned long _nt; // Number of time steps for trajectory
        double _dt; // Time step for trajectory
        double _gravity; // Global gravitational constant
        double _temperature; // Global gravitational constant
        unsigned long _proposals; // Number of iterations for Monte Carlo sampling
        bool _genMomKinetic;
        bool _genMomPropose;
//        bool _norMom;
        bool _testBefore;
        bool _hmc;
        char *_outfile;
        winsize _window;

        arma::mat _massMatrix;
        arma::mat _CholeskyLowerMassMatrix;
        arma::mat _inverseMassMatrix; // needed to write Hamilton's equations in vector form
        arma::mat _inverseMassMatrixDiagonal; // needed to write Hamilton's equations in vector form

        // Precomputed misfit function size
        arma::mat _A;
        arma::vec _bT; // Because I haven't coded up the actual left multiplication of vector-matrices
        double _c;

        // Member functions
        void propose_metropolis();

        void propose_hamilton(int &uturns);

        void leap_frog(int &uturns, bool writeTrajectory);

        double chi();

        double energy();

        void write_sample(std::ofstream &outfile, double misfit);

        double precomp_misfit();

        arma::vec precomp_misfitGrad();

        double kineticEnergy();


        bool _ergodic;
        double _acceptanceFactor;
    };
}

#endif //HMC_LINEAR_SYSTEM_SAMPLER_HPP
