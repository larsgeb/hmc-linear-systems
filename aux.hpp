//
// Created by Lars Gebraad on 7/3/17.
//

#ifndef HMC_VSP_AUX_HPP
#define HMC_VSP_AUX_HPP

/*=================================================================================*/
/* Inversion and model parameters. ------------------------------------------------*/
/*=================================================================================*/
class parameters;

class data;

class parameters {
    // This seems to be working satisfactorily.
public:
    int nLayers;
    int Nq;
    // Maybe implement a separate class for a priori information
    std::vector<double> q,
            sigma_q,
            mean_q;
    // m_q, sigma_q and mean_q contain layer speeds (0:nLayers, nLayers+1) and layer thicknesses (nLayers+1:nLayers*2 nLayers)
    // Usage of std::vector allows us to not define the model dimensions in the code, but rather in the parameters_synthetics_model.txt-file (dynamic allocation of memory)
    std::vector<std::vector<double> > iCM;


    std::vector<double> layerBase; // should be automatically computed from m_q within constructors/setters
    double t0U; // A scalar for the misfit at the starting model (for Taylor expansion)
    std::vector<double> t1U; // A vector containing all first derivatives of Xi evaluated at the starting model (for Taylor expansion).
    std::vector<std::vector<double> > t2U; // A matrix containing all second derivatives of Xi evaluated at the starting model (for Taylor expansion).

    parameters();         /**< Constructor. */
    ~parameters();        /**< Destructor. */
    parameters(const parameters &a);            // Copy constructor.
    parameters &operator=(const parameters &a); // Assignment operator.

    /** Member functions. */
    void
    read_input(const char *filename);  /**< Fill values by reading input file. */
    void tExpand(data d, int order, double ratioStep);

    void calculateLayerBase();

    void calculateInverseCM();
};

parameters operator+(const parameters &a, const parameters &b);

parameters operator-(const parameters &a, const parameters &b);

/*=================================================================================*/
/* VSP data class. ----------------------------------------------------------------*/
/*=================================================================================*/

class data {
public:

    /** data setup. */
    int recN;      /**< Number of receivers. */
    std::vector<double> recT; /**< Time observed. */
    std::vector<double> recZ; /**< Depth of receiver. */

    std::vector<std::vector<double> > iCD; // Inverse data covariance matrix

    /** Constructor and destructor. */
    data();         /**< Constructor. */
    ~data();

    /**< Destructor. */

    /** Member functions. */

    void make_synthetics(parameters &q);

    std::vector<double> forwardModel(parameters &q);

    void read_data(
            const char *filename);           /**< Fill data from previously computed and stored file. */
    void write(const char *filename);               /**< Write data to a file. */
    double misfit(parameters *starting_q);

    void calculateInverseCD(
            double MeasurementErrorStd); // Calculates data covariance matrix based upon a single standard deviation
    std::vector<double> misfitT1(parameters *starting_q, double ratioStep);

    std::vector<std::vector<double> >
    misfitT2(parameters *starting_q, double ratioStep);
};

// Some useful linear algebra functions
std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B);

std::vector<double>
MatrixVectorProduct(std::vector<std::vector<double> > M, std::vector<double> A);

double VectorVectorProduct(std::vector<double> A, std::vector<double> B);

/** Double-valued, uniformly distributed random numbers. */
double randf(
        double min,    /**< Minimum value. */
        double max     /**< Maximum value. */
);

/** Double-valued, normally distributed random numbers. */
void randn(
        double mean,        /**< Mean. */
        double stdv,        /**< Standard deviation. */
        double *x1,         /**< Pointer to first random number. */
        double *x2          /**< Pointer to second random number. */
);

double randn(
        double mean,       /**< Mean. */
        double stdv        /**< Standard deviation. */
);

#endif //HMC_VSP_AUX_HPP
