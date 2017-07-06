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
    std::vector<double> q,
            sigma_q,
            mean_q;
    // q, sigma_q and mean_q contain layer speeds (0:nLayers, nLayers+1) and layer thicknesses (nLayers+1:nLayers*2 nLayers)
    // Usage of std::vector allows us to not define the model dimensions in the code, but rather in the parameters.txt-file (dynamic allocation of memory)
    std::vector<std::vector<double>> iCM;


    std::vector<double> layerBase; // should be automatically computed from q within constructors/setters
    double t0U; // A scalar for the misfit at the starting model (for Taylor expansion)
    std::vector<double> t1U; // A vector containing all first derivatives of Xi evaluated at the starting model (for Taylor expansion).


    parameters();         /**< Constructor. */
    ~parameters();        /**< Destructor. */
    parameters(const parameters &a);            // Copy constructor.
    parameters &operator=(const parameters &a); // Assignment operator.

    /** Member functions. */
    void read_input(const char *filename);  /**< Fill values by reading input file. */
    void tExpand(data d, int order,double ratioStep);
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

    /** Constructor and destructor. */
    data();         /**< Constructor. */
    ~data();

    /**< Destructor. */

    /** Member functions. */

    void make_synthetics(parameters &q);
    std::vector<double> forwardModel(parameters &q);
    void read_data(const char *filename);           /**< Fill data from previously computed and stored file. */
    void write(const char *filename);               /**< Write data to a file. */
    double misfit(parameters *q);

    std::vector<double> misfitT1(parameters *q,double ratioStep);
};

// Some useful linear algebra functions
std::vector<double> VectorDifference(std::vector<double> A,std::vector<double> B);
std::vector<double> MatrixVectorProduct(std::vector<std::vector<double>> M, std::vector<double> A);
double VectorVectorProduct(std::vector<double> A, std::vector<double> B);

#endif //HMC_VSP_AUX_HPP
