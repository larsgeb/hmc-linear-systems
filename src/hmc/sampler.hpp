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

    struct InversionSettings {
        // Defaults
        const double PI = 3.14159265358979323846264338327;
        double _temperature = 1.0;
        unsigned long int _proposals = 1000;
        unsigned long int _trajectorySteps = 10;
        unsigned long int _massMatrixType = 0;
        struct winsize _window{};
        double _means = 1;
        double _std_dev = 1;
        char *_outputSamples = const_cast<char *>("OUTPUT/samples.txt");
        char *_outputTrajectory = const_cast<char *>("OUTPUT/trajectory.txt");
        char *_inputMatrix = const_cast<char *>("INPUT/matrix.txt");
        char *_inputData = const_cast<char *>("INPUT/data.txt");
        bool _algorithmNew = true;
        bool _genMomPropose = true; // Use generalized mass matrix to propose new momenta (true).
        bool _genMomKinetic = true; // Use generalized mass matrix to compute kinetic energy (true).
        bool _testBefore = true; // Decreases required computation time by order of magnitude, no other influence.
        bool _ergodic = true;  // Randomizes trajectory length and step size
        bool _hamiltonianMonteCarlo = true; // Metropolis Hastings (false) or Hamiltonian Monte Carlo (true).
        double _timeStep;

        InversionSettings &setTimeStepFromGrav_nSteps() {
            _timeStep = 2.0 * PI / _trajectorySteps;
        };


        // Parse command line options
        void parse_input(int argc, char *argv[]) {
            for (int i = 1; i < argc; i++) {

                if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                    std::cout << "Display Help!" << std::endl;
                    exit(EXIT_SUCCESS);
                }

                if (i + 1 != argc) {
                    if (strcmp(argv[i], "-im") == 0 || strcmp(argv[i], "--inputmatrix") == 0) {
                        _inputMatrix = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-id") == 0 || strcmp(argv[i], "--inputdata") == 0) {
                        _inputData = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-mtype") == 0 || strcmp(argv[i], "--massmatrixtype") == 0) {
                        parse_long_unsigned(argv, i, _massMatrixType);
                        i++;
                    } else if (strcmp(argv[i], "-os") == 0 || strcmp(argv[i], "--outputsamples") == 0) {
                        _outputSamples = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-ot") == 0 || strcmp(argv[i], "--outputtrajectory") == 0) {
                        _outputTrajectory = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--temperature") == 0) {
                        parse_double(argv, i, _temperature);
                        i++;
                    } else if (strcmp(argv[i], "-nt") == 0 || strcmp(argv[i], "--trajectorysteps") == 0) {
                        parse_long_unsigned(argv, i, _trajectorySteps);
                        i++;
                    } else if (strcmp(argv[i], "-dt") == 0 || strcmp(argv[i], "--timestep") == 0) {
                        parse_double(argv, i, _timeStep);
                        i++;
                    } else if (strcmp(argv[i], "-ns") == 0 || strcmp(argv[i], "--numberofsamples") == 0) {
                        parse_long_unsigned(argv, i, _proposals);
                        i++;
                    } else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--ergodic") == 0) {
                        parse_boolean(argv, i, _ergodic);
                        i++;
                    } else if (strcmp(argv[i], "-gmp") == 0 || strcmp(argv[i], "--generalizedmomentumpropose") == 0) {
                        parse_boolean(argv, i, _genMomPropose);
                        i++;
                    } else if (strcmp(argv[i], "-gmc") == 0 || strcmp(argv[i], "--generalizedmomentumcalculate") == 0) {
                        parse_boolean(argv, i, _genMomKinetic);
                        i++;
                    } else if (strcmp(argv[i], "-Hb") == 0 || strcmp(argv[i], "--hamiltonianbefore") == 0) {
                        parse_boolean(argv, i, _testBefore);
                        i++;
                    } else if (strcmp(argv[i], "-an") == 0 || strcmp(argv[i], "--algorithmnew") == 0) {
                        parse_boolean(argv, i, _algorithmNew);
                        i++;
                    } else if (strcmp(argv[i], "-means") == 0 || strcmp(argv[i], "--means") == 0) {
                        parse_double(argv, i, _means);
                        i++;
                    } else if (strcmp(argv[i], "-std") == 0 || strcmp(argv[i], "--standarddeviation") == 0) {
                        parse_double(argv, i, _std_dev);
                        i++;
                    }
                }
            }
        }

        // Check input for errors
        void parse_long_unsigned(char *argv[], int position, long unsigned int &to_set) {
            std::stringstream ss(argv[position + 1]);
            long unsigned x;
            if (!(ss >> x) || !ss.eof())
                std::cerr << "Invalid number " << argv[position + 1] << " for option " << argv[position] << std::endl;
            else
                to_set = x;
        }

        void parse_double(char *argv[], int position, double &to_set) {
            std::stringstream ss(argv[position + 1]);
            double x;
            if (!(ss >> x) || !ss.eof())
                std::cerr << "Invalid number " << argv[position + 1] << " for option " << argv[position] << std::endl;
            else
                to_set = x;
        }

        void parse_boolean(char *argv[], int position, bool &to_set) {
            std::stringstream ss(argv[position + 1]);
            bool x;
            if (!(ss >> x) || !ss.eof())
                std::cerr << "Invalid number " << argv[position + 1] << " for option " << argv[position] << std::endl;
            else
                to_set = x;
        }

        // Constructor
        InversionSettings(int argc, char *argv[]) {
            // Constructor for settings

            // Use standard timestep defined by gravity
            _timeStep = 0.05;
            setTimeStepFromGrav_nSteps();

            // Parse command line input
            parse_input(argc, argv);

            // Try to get the command window size for displaying results, otherwise default
            ioctl(STDOUT_FILENO, TIOCGWINSZ, &_window);
            if (_window.ws_col < 5)
                _window.ws_col = 20;
            if (_window.ws_row < 5)
                _window.ws_row = 20;
        }
    };

    class sampler {
    public:
        // Constructors and destructors
        sampler(prior &prior, data &data, forward_model &model, InversionSettings settings);

        // Sample the posterior and write samples out to file
        void sample();

        void sample_new();

        void sample_neal();

        // Set the starting model explicitly instead of prior-based
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
//        Posterior _posterior; // Necessary?
        unsigned long _nt; // Number of time steps for trajectory
        double _dt; // Time step for trajectory
        double _REMOVETHIS; // Global gravitational constant
        double _temperature; // Global gravitational constant
        unsigned long _proposals; // Number of iterations for Monte Carlo sampling
        unsigned long _massMatrixType; // Number of iterations for Monte Carlo sampling
        bool _genMomKinetic;
        bool _genMomPropose;
        bool _algorithmNew;
//        bool _norMom;
        bool _testBefore;
        bool _hmc;
        char *_outputSamples;
        char *_outputTrajectory;
        char *_inputMatrix;
        char *_inputData;
        winsize _window;
        bool _ergodic;

        arma::mat *_massMatrix;
        arma::mat _optionalMassMatrixMemory;
        arma::mat _CholeskyLowerMassMatrix;
        arma::mat _inverseMassMatrix; // needed to write Hamilton's equations in vector form
        arma::mat _inverseMassMatrixDiagonal; // needed to write Hamilton's equations in vector form

        // Precomputed misfit function size
        arma::mat _A;
        arma::vec _bT; // Because I haven't coded up the actual left multiplication of vector-matrices
        double _c;

        // Member functions
        void propose_metropolis();

        void propose_momentum();

        void leap_frog(int &uturns, bool writeTrajectory);

        double chi();

        double energy();

        void write_sample(std::ofstream &outfile, double misfit);

        double precomp_misfit();

        arma::vec precomp_misfitGrad();

        double kineticEnergy();
    };
}

#endif //HMC_LINEAR_SYSTEM_SAMPLER_HPP
