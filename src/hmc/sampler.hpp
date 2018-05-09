//
// Created by Lars Gebraad on 18-8-17.
//

#ifndef HMC_LINEAR_SYSTEM_SAMPLER_HPP
#define HMC_LINEAR_SYSTEM_SAMPLER_HPP

#include "prior.hpp"
#include "data.hpp"

#include <sys/ioctl.h>
#include <cstdio>
#include <unistd.h>

namespace hmc {

    struct InversionSettings {
        // Defaults
        const double PI = 3.14159265358979323846264338327;
        struct winsize _window{};

        // Output files
        char *_outputSamplesFile = const_cast<char *>("OUTPUT/samples.txt");
        char *_outputTrajectoryFile = const_cast<char *>("OUTPUT/trajectory.txt");

        // Forward model
        char *_inputMatrixFile = const_cast<char *>("INPUT/matrix.txt");

        // Prior stuff
        bool _fixedPrior = true;
        double _fixedPriorMean = 1;
        double _fixedPriorCovariance = 1;
        char *_priorCovarianceFile  = const_cast<char *>("");
        char *_priorMeanFile  = const_cast<char *>("");

        // Data stuff
        int _useFixedDataCovariance = 2; // 0 is not fixed, 1 is absolute, 2 is percentual
        double _fixedDataCovariance = 1; // being absolute or percentual
        char *_inputDataFile = const_cast<char *>("INPUT/data.txt");
        char *_dataCovarianceFile  = const_cast<char *>("");

        // Tuning parameters
        double _timeStep = 0.1;
        double _temperature = 1.0;
        unsigned long int _proposals = 1000;
        unsigned long int _trajectorySteps = 10;
        unsigned long int _massMatrixType = 0;

        // Other options
        bool _algorithmNew = true;
        bool _genMomPropose = true; // Use generalized mass matrix to propose new momenta (true).
        bool _genMomKinetic = true; // Use generalized mass matrix to compute kinetic energy (true).
        bool _testBefore = true; // Decreases required computation time by order of magnitude, no other influence.
        bool _ergodic = true;  // Randomizes trajectory length and step size
        bool _hamiltonianMonteCarlo = true; // Metropolis Hastings (false) or Hamiltonian Monte Carlo (true).
        bool _adaptTimestep = true; // adapt timestep for mass-matrix choice


        // Parse command line options
        void parse_input(int argc, char *argv[]) {

            int a = 1;
            if (argc < 2) {
                display_help();
                exit(EXIT_SUCCESS);
            }
            for (int i = 1; i < argc; i++) {
                if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                    display_help();
                    exit(EXIT_SUCCESS);
                }
                if (i + 1 != argc) {
                    if (strcmp(argv[i], "-im") == 0 || strcmp(argv[i], "--inputmatrix") == 0) {
                        _inputMatrixFile = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-id") == 0 || strcmp(argv[i], "--inputdata") == 0) {
                        _inputDataFile = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-icm") == 0 || strcmp(argv[i], "--inputparametercovariance") == 0) {
                        _dataCovarianceFile = (argv[i + 1]);
                        _fixedPriorCovariance = false;
                        i++;
                    } else if (strcmp(argv[i], "-icd") == 0 || strcmp(argv[i], "--inputdatacovariance") == 0) {
                        _priorCovarianceFile = (argv[i + 1]);
                        _useFixedDataCovariance = false;
                        i++;
                    } else if (strcmp(argv[i], "-mtype") == 0 || strcmp(argv[i], "--massmatrixtype") == 0) {
                        parse_long_unsigned(argv, i, _massMatrixType);
                        i++;
                    } else if (strcmp(argv[i], "-os") == 0 || strcmp(argv[i], "--outputsamples") == 0) {
                        _outputSamplesFile = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-ot") == 0 || strcmp(argv[i], "--outputtrajectory") == 0) {
                        _outputTrajectoryFile = (argv[i + 1]);
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
                    } else if (strcmp(argv[i], "-at") == 0 || strcmp(argv[i], "--adapttimestep") == 0) {
                        parse_boolean(argv, i, _adaptTimestep);
                        i++;
                    } else if (strcmp(argv[i], "-ns") == 0 || strcmp(argv[i], "--numberofsamples") == 0) {
                        parse_long_unsigned(argv, i, _proposals);
                        i++;
                    } else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--ergodic") == 0) {
                        parse_boolean(argv, i, _ergodic);
                        i++;
                    } else if (strcmp(argv[i], "-percd") == 0 || strcmp(argv[i], "--percentualdatacovariance") == 0) {
                        if(_useFixedDataCovariance != 1){
                            std::cout << "Option clash for fixed data covariance (percentual) together with input data "
                                    "covariance matrix, exiting...";
                            exit(EXIT_SUCCESS);
                        }else{
                            parse_double(argv, i, _fixedDataCovariance);
                        }
                        i++;
                    } else if (strcmp(argv[i], "-abscd") == 0 || strcmp(argv[i], "--absolutedatacovariance") == 0) {
                        if(_useFixedDataCovariance != 2){
                            std::cout << "Option clash for fixed data covariance (absolute) together with input data "
                                    "covariance matrix, exiting...";
                            exit(EXIT_SUCCESS);
                        }else{
                            parse_double(argv, i, _fixedDataCovariance);
                        }
                        i++;
                    } else if (strcmp(argv[i], "-gmp") == 0 || strcmp(argv[i], "--correlatedmomenta") == 0) {
                        parse_boolean(argv, i, _genMomPropose);
                        i++;
                    } else if (strcmp(argv[i], "-gmc") == 0 || strcmp(argv[i], "--generalkinetic") == 0) {
                        parse_boolean(argv, i, _genMomKinetic);
                        i++;
                    } else if (strcmp(argv[i], "-Hb") == 0 || strcmp(argv[i], "--hamiltonianbefore") == 0) {
                        parse_boolean(argv, i, _testBefore);
                        i++;
                    } else if (strcmp(argv[i], "-an") == 0 || strcmp(argv[i], "--algorithmnew") == 0) {
                        parse_boolean(argv, i, _algorithmNew);
                        i++;
                    } else if (strcmp(argv[i], "-means") == 0 || strcmp(argv[i], "--means") == 0) {
                        parse_double(argv, i, _fixedPriorMean);
                        i++;
                    } else if (strcmp(argv[i], "-std") == 0 || strcmp(argv[i], "--standarddeviation") == 0) {
                        parse_double(argv, i, _fixedPriorCovariance);
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
                std::cerr << "Invalid input " << argv[position + 1] << " for option " << argv[position] << std::endl;
            else
                to_set = x;
        }

        void parse_double(char *argv[], int position, double &to_set) {
            std::stringstream ss(argv[position + 1]);
            double x;
            if (!(ss >> x) || !ss.eof())
                std::cerr << "Invalid input " << argv[position + 1] << " for option " << argv[position] << std::endl;
            else
                to_set = x;
        }

        void parse_boolean(char *argv[], int position, bool &to_set) {
            std::stringstream ss(argv[position + 1]);
            bool x;
            if (!(ss >> x) || !ss.eof())
                std::cerr << "Invalid input " << argv[position + 1] << " for option " << argv[position] << std::endl;
            else
                to_set = x;
        }

        void display_help() {
            std::cout << std::endl << "Hamiltonian Monte Carlo Sampler" << std::endl
                      << "Lars Gebraad, 2017" << std::endl << "Displaying help ..." << std::endl << std::endl;

            std::cout << "\tFiles" << std::endl
                      << "\t\t \033[1;31m -im \033[0m (existing path to non-existing file, required)" << std::endl
                      << "\t\t input matrix file (G), every line of text file should be a matrix row, \r\n\t\t entries "
                              "separated"
                      << " by spaces" << std::endl
                      << "\t\t \033[1;31m -id \033[0m (existing path to non-existing file, required)" << std::endl
                      << "\t\t input data file (d), every datapoint should be a new line" << std::endl
                      << "\t\t \033[1;31m -os \033[0m (existing path to non-existing file, required)" << std::endl
                      << "\t\t output samples file" << std::endl
                      << "\t\t \033[1;31m -ot \033[0m (existing path to non-existing file, required)" << std::endl
                      << "\t\t output trajectory file" << std::endl << std::endl
                      << "\tPrior information" << std::endl
                      << "\t\t \033[1;31m -means \033[0m (double, required)" << std::endl
                      << "\t\t Prior means" << std::endl
                      << "\t\t \033[1;31m -std \033[0m (double, required)" << std::endl
                      << "\t\t Prior standard deviation (square root of variance)" << std::endl << std::endl
                      << "\tTuning parameters" << std::endl
                      << "\t\t \033[1;32m -dt \033[0m (double, default = adaptive)" << std::endl
                      << "\t\t size of time discretization steps" << std::endl
                      << "\t\t \033[1;32m -nt \033[0m (integer, default = 10)" << std::endl
                      << "\t\t number of time discretization steps" << std::endl
                      << "\t\t \033[1;32m -t \033[0m (double, default = 1)" << std::endl
                      << "\t\t temperature" << std::endl
                      << "\t\t \033[1;32m -mtype \033[0m (0, 1 or 2, default = 0)" << std::endl
                      << "\t\t mass matrix type: full ideal (0), diagonal ideal (1) or unit matrix (2)"
                      << std::endl
                      << std::endl
                      << "\tOther options" << std::endl
                      << "\t\t \033[1;32m -ns \033[0m (integer, default = 1000)" << std::endl
                      << "\t\t number of proposals" << std::endl
                      << "\t\t \033[1;32m -at \033[0m (boolean, default = 1) " << std::endl
                      << "\t\t adapt timestep to be stable, using eigen-decomposition of the term sqrt(Q^-1 A)"
                      << std::endl
                      << "\t\t \033[1;32m -e\033[0m (boolean, default = 1)" << std::endl
                      << "\t\t ensure ergodicity of the sampler by uniformly modifying nt and dt by \r\n\t\t 0.5-1.5 "
                              "(randomly) per sample" << std::endl
                      << "\t\t \033[1;32m -gmp\033[0m (boolean, default = 1)" << std::endl
                      << "\t\t use full mass matrix to propose new momenta (correlated samples)" << std::endl
                      << "\t\t \033[1;32m -gmc\033[0m (boolean, default = 1)" << std::endl
                      << "\t\t use full mass matrix to calculate kinetic energy instead of diagonal" << std::endl
                      << "\t\t \033[1;32m -Hb\033[0m (boolean, default = 1)" << std::endl
                      << "\t\t use conservation of energy to evaluate the Hamiltonian and acceptance \r\n\t\t     criterion "
                              "before "
                              "propagating, only for algorithm 1" << std::endl
                      << "\t\t \033[1;32m -an\033[0m (boolean, default = 1)" << std::endl
                      << "\t\t choose HMC algorithm; classic (0), new (1)" << std::endl << std::endl
                      << "\tFor examples, see test/" << std::endl << std::endl;
        }

        // Constructor
        InversionSettings(int argc, char *argv[]) {
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
        explicit sampler(InversionSettings settings);

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

        prior _prior;
        forward_model _model;
        data _data;
    private:
        // Fields

        unsigned long _nt; // Number of time steps for trajectory
        double _dt; // Time step for trajectory
        double _temperature; // Temperature for acceptance criterion
        unsigned long _proposals; // Number of iterations for Monte Carlo sampling
        unsigned long _massMatrixType; // Number of iterations for Monte Carlo sampling
        bool _genMomKinetic;
        bool _genMomPropose;
        bool _algorithmNew;
        bool _testBefore;
        bool _hmc;
        int _useFixedDataCovariance;
        char *_outputSamples;
        char *_dataCovarianceFile;
        char *_outputTrajectory;
        char *_inputMatrix;
        winsize _window;
        bool _ergodic;
        bool _adaptTimestep;
        double _fixedDataCovariance;

        // Data stuff
        char *_inputData;

        // Prior stuff
        bool _fixedPrior;
        double _fixedPriorCovariance;
        double _fixedPriorMean;
        char *_priorMeanFile;
        char *_priorCovarianceFile;

        arma::mat *_massMatrix;
        arma::mat _optionalMassMatrixMemory;
        arma::mat _CholeskyLowerMassMatrix;
        arma::mat _inverseMassMatrix; // needed to write Hamilton's equations in vector form
        arma::mat _inverseMassMatrixDiagonal; // needed to write Hamilton's equations in vector form TODO rewrite as sparse

        // Precomputed misfit function size
        arma::mat _A;
        arma::vec _bT;
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

double get_wall_time();

double get_cpu_time();

#endif //HMC_LINEAR_SYSTEM_SAMPLER_HPP
