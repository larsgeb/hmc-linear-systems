/*
 * Created by Lars Gebraad on 18-08-17.
 * Last modified 29-05-18.
 */

#ifndef HMC_LINEAR_SYSTEM_SAMPLER_HPP
#define HMC_LINEAR_SYSTEM_SAMPLER_HPP

#include <sys/ioctl.h>
#include <cstdio>
#include <unistd.h>
#include <armadillo>

using namespace arma;

namespace hmc {
    // settings structural object, mostly done to tidy up the linearSampler class
    struct InversionSettings {

        // Defaults
        const double PI = 3.14159265358979323846264338327;
        struct winsize _window{};

        // Output files
        char *_outputSamplesFile = const_cast<char *>("OUTPUT/samples.txt");
        char *_outputTrajectoryFile = const_cast<char *>("OUTPUT/trajectory.txt");

        // ABC-style
        char *A_file = const_cast<char *>("");
        char *B_file = const_cast<char *>("");
        char *C_file = const_cast<char *>("");

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
                    if (strcmp(argv[i], "-ia") == 0 || strcmp(argv[i], "--inputA") == 0) {
                        A_file = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-ib") == 0 || strcmp(argv[i], "--inputB") == 0) {
                        B_file = (argv[i + 1]);
                        i++;
                    } else if (strcmp(argv[i], "-ic") == 0 || strcmp(argv[i], "--inputC") == 0) {
                        C_file = (argv[i + 1]);
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

        // Parsing double input
        void parse_double(char *argv[], int position, double &to_set) {
            std::stringstream ss(argv[position + 1]);
            double x;
            if (!(ss >> x) || !ss.eof())
                std::cerr << "Invalid input " << argv[position + 1] << " for option " << argv[position] << std::endl;
            else
                to_set = x;
        }

        // Parsing boolean input
        void parse_boolean(char *argv[], int position, bool &to_set) {
            std::stringstream ss(argv[position + 1]);
            bool x;
            if (!(ss >> x) || !ss.eof())
                std::cerr << "Invalid input " << argv[position + 1] << " for option " << argv[position] << std::endl;
            else
                to_set = x;
        }

        // Display help for sampler
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
                      << "\t\t ensure ergodicity of the linearSampler by uniformly modifying nt and dt by \r\n\t\t 0.5-1.5 "
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
                      << "\tFor examples, see inversions/" << std::endl << std::endl;
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

    class linearSampler {
    public:
        /** \brief Constructor for a probabilistic sampler.
          * \param InversionSettings settings
          * \return linearSampler object
          * */
        explicit linearSampler(InversionSettings settings);

        /** \brief General wrapper for any sampler
          * \return void
          * */
        void sample();

        /** \brief Set the starting model explicitly instead of prior-based.
          * \param arma::vec startingModel
          * \return void
          * */
        void setStarting(arma::vec &model);

        /** \brief Method for sampling using the criterion as described in Neal's HMC introduction.
          * \return void
          * */
        void sample_neal();

    private:
        // States
        vec _currentModel; ///< State of markov chain describing coordinates of current point.
        vec _proposedModel; ///< State of markov chain describing coordinates of proposal.
        vec _proposedMomentum; ///< State of markov chain describing momentum of proposal.

        // Quadratic form
        mat A; ///< A in quadratic form.
        mat At; ///< A^t in quadratic form, speeds up proposals at the cost of memory.
        bool symmetricA = false; ///< boolean for symmetry of A in quadratic form, speeds up gradient calculations if true.
        colvec B; ///< B in quadratic form.
        double C; ///< C in quadratic form.

        // Mass matrices
        mat massMatrix; ///< Mass matrix for HMC.
        mat invMass; ///< Inverse mass matrix for calculation of kinetic energy in HMC.

        // Settings
        unsigned long nt; ///< Number of time steps for trajectory in HMC.
        double dt; ///< Time step for trajectory in HMC.
        double temperature; ///< Temperature for acceptance criterion in MCMC.
        unsigned long proposals; ///< Number of proposals for MCMC.
        unsigned long massMatrixType; ///< Number of iterations in HMC.
        winsize window; ///< Size of terminal for nice output.

        // Pointers to files
        char *A_file; ///< Pointer to character array of filename containing A in the quadratic form.
        char *B_file; ///< Pointer to character array of filename containing B in the quadratic form.
        char *C_file; ///< Pointer to character array of filename containing C in the quadratic form.
        char *_outputSamples; ///< Pointer to character array of filename to store samples in MCMC.
        char *_outputTrajectory; ///< Pointer to character array of filename to store trajectory samples from HMC.

        // Member methods

        //
        /** \brief Propose new momentum according to N(0,M), writes to \ref linearSampler::_proposedMomentum.
          * \return void
          * */
        void propose_momentum();

        // Integrate Hamilton's equations using a leapfrog scheme
        void leap_frog(bool writeTrajectory);

        // Evaluate misfit (chi)
        double chi();

        // Evaluate Hamiltonian (H)
        double energy();

        // Write sample to one line of opened filestream
        void write_sample(std::ofstream &outfile, double misfit);

        // Calculate misfit of quadratic form
        double misfit();

        // Calculate kinetic energy as 1/2 pt M^-1 p
        double kineticEnergy();

        arma::mat CholeskyLowerMassMatrix;
    };
}

double get_wall_time();

double get_cpu_time();

#endif //HMC_LINEAR_SYSTEM_SAMPLER_HPP
