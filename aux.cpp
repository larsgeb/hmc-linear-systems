//
// Created by Lars Gebraad on 7/3/17.
//

#include <time.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "aux.hpp"

const double PI = 3.14159265358979323846264338327;
/*=================================================================================*/
/* Inversion and model parameters class. ------------------------------------------*/
/*=================================================================================*/


/* Constructor. -----------------------------------------------------------------------------*/
parameters::parameters() {

}

/* Destructor. --------------------------------------------------------------------*/
parameters::~parameters() {

}

/* Copy constructor. --------------------------------------------------------------*/
parameters::parameters(const parameters &a) {
    for (int i = 0; i < a.nLayers; i++) {
        q[i] = a.q[i];
        sigma_q[i] = a.sigma_q[i];
        mean_q[i] = a.mean_q[i];
        layerBase[i] = a.layerBase[i];
    }
    nLayers = a.nLayers;
    Nq = a.Nq;
}

/* Assignment operator. ----------------------------------------------------------------------------------------------*/
parameters &parameters::operator=(const parameters &a) {
    /* Check for self-assignment. */
    if (this == &a) return *this;
    q.clear();
    mean_q.clear();
    sigma_q.clear();
    layerBase.clear();
    for (int i = 0; i < a.q.size(); i++) {
        q.push_back(a.q[i]);
        mean_q.push_back(a.mean_q[i]);
        sigma_q.push_back(a.sigma_q[i]);
        layerBase.push_back(a.layerBase[i]);
    }
//    calculateLayerBase(); // Not needed because assignment allows copying
    return *this;
}

/* Addition. ---------------------------------------------------------------------------------------------------------*/
parameters operator+(const parameters &a, const parameters &b) {
    // add values (in m_q), takes mean & sigma of a
    parameters c;

    for (int i = 0; i < a.q.size(); i++) {
        c.q.push_back(a.q[i] + b.q[i]);
        c.mean_q.push_back(a.mean_q[i]);
        c.sigma_q.push_back(a.sigma_q[i]);
    }
    c.calculateLayerBase();
    return c;
}

/* Subtraction. ------------------------------------------------------------------------------------------------------*/
parameters operator-(const parameters &a, const parameters &b) {
    parameters c;

    for (int i = 0; i < a.q.size(); i++) {
        c.q.push_back(a.q[i] - b.q[i]);
        c.mean_q.push_back(a.mean_q[i]);
        c.sigma_q.push_back(a.sigma_q[i]);
    }

    return c;
}

/* Read input from file. ---------------------------------------------------------------------------------------------*/
void parameters::read_input(const char *filename) {
    FILE *fid;
    char str[1000];

    fid = fopen(filename, "r");

    fgets(str, 1000, fid); // Move line
    fscanf(fid, "%d", &nLayers);
    fgets(str, 1000, fid);

    Nq = nLayers * 2 + 1;

    double temp;

    /*- Layer Speeds -------------------------------------------------------------------------------------------------*/
    fgets(str, 1000, fid);
    // Scan all priori mean layer speeds (number of layers + half space)
    for (int i = 0; i < nLayers + 1; i++) {
        fscanf(fid, "%lg", &temp);
        mean_q.push_back(temp);
    }
    fgets(str, 1000, fid);
    // Scan all priori sigma layer speeds (number of layers + half space)
    for (int i = 0; i < nLayers + 1; i++) {
        fscanf(fid, "%lg", &temp);
        sigma_q.push_back(temp);
    }
    fgets(str, 1000, fid);

    /*- Layer thicknesses --------------------------------------------------------------------------------------------*/
    fgets(str, 1000, fid);
    for (int i = nLayers + 1; i < nLayers * 2 + 1; i++) {
        fscanf(fid, "%lg", &temp);
        mean_q.push_back(temp);
    }
    fgets(str, 1000, fid);
    for (int i = nLayers + 1; i < nLayers * 2 + 1; i++) {
        fscanf(fid, "%lg", &temp);
        sigma_q.push_back(temp);
    }
    fgets(str, 1000, fid);

    /* Assign values equal to prior means. */
    for (int i = 0; i < mean_q.size(); i++) {
        q.push_back(mean_q[i]);
    }

    calculateLayerBase();
    calculateInverseCM();

    fclose(fid);
}

/* Calculate layer bases ---------------------------------------------------------------------------------------------*/
void parameters::calculateLayerBase() {
    for (int i = 0; i < nLayers; i++) {
        if (i == 0) {
            layerBase.push_back(q[nLayers + 1 + i]);
        } else {
            layerBase.push_back(layerBase[i - 1] + q[nLayers + 1 + i]);
        }
    }
}

/* Calculate taylor expansion of starting model ----------------------------------------------------------------------*/
void parameters::tExpand(data d, int order, double ratioStep) {
    if (order != 1) {
        return;
    }
    t0U = d.misfit(this);   // Calculate Xi
    t1U = d.misfitT1(this, ratioStep); // Calculate and evaluate first derivative
    t2U = d.misfitT2(this, ratioStep); // Calculate and evaluate first derivative
}

/* Calculate inverse covariance matrix -------------------------------------------------------------------------------*/
void parameters::calculateInverseCM() {
    // Note that only uncorrelated prior is as of yet implemented (unit matrices)
    std::vector<double> row_iCM((unsigned long) Nq, 0);
    iCM.insert(iCM.end(), (unsigned long) Nq, row_iCM);
    row_iCM.clear();

    for (int i = 0; i < Nq; i++) {
        iCM[i][i] = 1 / sigma_q[i];
    }
}

/*====================================================================================================================*/
/* VSP data class. ---------------------------------------------------------------------------------------------------*/
/*====================================================================================================================*/

data::data() {
    // Set up has to be always loaded from the same file, therefor it's hardcoded into the constructor.
    FILE *fid2;
    char str[1000];

    /* Read VSP setup. */
    fid2 = fopen("INPUT/setup.txt", "r");

    double receiverSpacing, initialLocation;

    fgets(str, 1000, fid2);
    fscanf(fid2, "%i %lf %lf", &recN, &receiverSpacing, &initialLocation);

    for (int i = 0; i < recN; i++) {
        recZ.push_back(initialLocation + receiverSpacing * i);
    }

    calculateInverseCD(500); // Hardcoded measurement error. Update as you see fit.
    fclose(fid2);
}

data::~data() {

}

void data::make_synthetics(parameters &q) {
    recT = forwardModel(q);
}

std::vector<double> data::forwardModel(parameters &q) {
    // This is the forward model code for a simple VSP where first arrivals are calculated from layer thicknesses and speeds.
    // Remember, positive Z is downward
    // m_q is filled in the following way: 0 - nLayers: speeds || nLayers+1 - nLayers*2:
    std::vector<double> travelTime;
    travelTime.clear();
    for (int i = 0; i < recN; i++) {
        // Looping through receivers
        travelTime.push_back(0.0);
        for (int j = 0; j < q.nLayers; j++) {
            // Looping through layers
            double thickness, base, top;
            thickness = q.q[j + 1 + q.nLayers];
            base = q.layerBase[j];
            top = j != 0 ? q.layerBase[j - 1] : 0;

            if (recZ[i] >= base) {
                travelTime[i] += thickness / q.q[j];
            } else if (recZ[i] > top) {
                travelTime[i] += (recZ[i] - top) / q.q[j];
            }
        }
        if (recZ[i] > q.layerBase[q.nLayers - 1]) {
            travelTime[i] += (recZ[i] - q.layerBase[q.nLayers - 1]) / q.q[q.nLayers];
        }
    }
    return travelTime;
}

void data::read_data(const char *filename) {
    // Be careful if data is already generated, will add more! (and will add a lot of debugging time..)
    double a;
    recT.clear(); // This actually mitigates the above described error
    std::ifstream infile(filename);
    for (int i = 0; i < recN; i++) {
        infile >> a;
        recT.push_back(a);
    }
}

void data::write(const char *filename) {
    FILE *fid;

    fid = fopen(filename, "w");

    for (int i = 0; i < recN; i++) {
        fprintf(fid, "%.15lg ", recT[i]);
        fprintf(fid, "\n");
    }
    fclose(fid);
}

double data::misfit(parameters *starting_q) {
    // This function only incorporates prior and forward model misfit. Added could be forward model errors.
    double Xi;

    // Maybe implement some extra functions for vector arithmetic
    std::vector<double> q_difference;
    q_difference = VectorDifference(starting_q->q, starting_q->mean_q);
    Xi = 0.5 * VectorVectorProduct(q_difference, MatrixVectorProduct(starting_q->iCM,q_difference));

    std::vector<double> d_difference;
    std::vector<double> d_forward;
    d_forward = forwardModel(*starting_q);
    d_difference = VectorDifference(d_forward, recT);
    Xi += 0.5 *
          VectorVectorProduct(d_difference, MatrixVectorProduct(iCD, d_difference));

    return Xi;
}

std::vector<double> data::misfitT1(parameters *starting_q, double ratioStep) {
    std::vector<double> firstDerivative;

    // Doing a Taylor expansion of the first order. This evaluates the derivatives at the starting model,
    // but note that it still would have to be multiplied by (m - m0) for the taylor function
    for (int i = 0; i < starting_q->Nq; i++) {
        double current_qi;
        current_qi = starting_q->q[i];
        parameters q_perturbed;
        q_perturbed.read_input("INPUT/parameters_synthetics_model.txt");
        q_perturbed.q[i] = current_qi * ratioStep;
        firstDerivative.push_back((this->misfit(&q_perturbed) - starting_q->t0U) /
                                  (current_qi * (1-ratioStep)));
    }
    return firstDerivative;
}

std::vector<std::vector<double> >
data::misfitT2(parameters *starting_q, double ratioStep) {
    std::vector<std::vector<double> > secondDerivative;
    std::vector<double> startingFirstDerivative;

    startingFirstDerivative = misfitT1(starting_q,ratioStep);
    for (int j = 0; j < starting_q->Nq; j++) {
        double current_qi;
        current_qi = starting_q->q[j];
        parameters q_perturbed;
        q_perturbed.read_input("INPUT/parameters_synthetics_model.txt");
        q_perturbed.q[j] = current_qi * ratioStep;
        secondDerivative.push_back(VectorDifference(misfitT1(&q_perturbed,ratioStep),
                            startingFirstDerivative));
        for (int i = 0; i<starting_q->Nq; i++){
            secondDerivative[j][i] *= 0.5 / (current_qi * (1-ratioStep));
        }

    }
    return secondDerivative;
}

void data::calculateInverseCD(double MeasurementErrorStd) {

    // A faster way of allocating a zero matrix.
    std::vector<double> row_iCD((unsigned long) recN, 0);
    iCD.insert(iCD.end(), (unsigned long) recN, row_iCD);
    row_iCD.clear();

    // Filling the diagonal
    for (int i = 0; i < recN; i++) {
        iCD[i][i] = 1 / MeasurementErrorStd;
    }
}

std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B) {

    std::vector<double> C;

    if (A.size() != B.size()) {
        // Get some Exception class to THROW.
        std::cout
                << "Vectors are not the same dimension! The code DIDN'T run successfully.";
        return C;
    }

    std::vector<double> q_difference;
    for (int i = 0; i < A.size(); i++) {
        // Prior misfit
        C.push_back(A[i] - B[i]);
    }
    return C;
}

std::vector<double>
MatrixVectorProduct(std::vector<std::vector<double> > M, std::vector<double> A) {
    std::vector<double> C;

    // Using std::vector<>.size() requires casting for clean compilation (seems unnecessary.. But oh well)
    // So watch out if you're working on 2^63 order problems.. ;) (Maximum <int>, not maximum <unsigned int>)
    int rowsM = static_cast<int>(M.size());
    int columnsM = static_cast<int>(M[0].size());
    int rowsA = static_cast<int>(A.size());

    if (columnsM != rowsA) {
        // Get some Exception class to THROW.
        std::cout
                << "Vector and matrix are not compatible in dimension! The code DIDN'T run successfully.";
        return C;
    }

    for (int i = 0; i < rowsM; i++) {
        C.push_back(0);
        for (int j = 0; j < columnsM; j++) {
            C[i] += M[i][j] * A[j];
        }
    }
    return C;
}

double VectorVectorProduct(std::vector<double> A, std::vector<double> B) {
    double C;

    if (A.size() != B.size()) {
        // Get some Exception class to THROW.
        std::cout
                << "Vectors are not compatible in dimension! The code DIDN'T run successfully.";
        return C;
    }
    C = 0;
    for (int i = 0; i < A.size(); i++) {
        C += (A[i] * B[i]);
    }
    return C;
}

/* Uniformly distributed, double-valued random numbers. ---------------------------*/
double randf(double min, double max) {
    return (max - min) * (double) rand() / RAND_MAX + min;
}

/* Normally distributed, double-valued random numbers. ----------------------------*/
/* This function implements the Box-Muller transform to obtain a pair of
 normally distributed random numbers with a given mean and standard deviation. */
void randn(double mean, double stdv, double *x1, double *x2) {
    double z1 = (double) rand() / RAND_MAX;
    double z2 = (double) rand() / RAND_MAX;

    *x1 = sqrt(-2.0 * log(z1)) * cos(2.0 * PI * z2);
    *x2 = sqrt(-2.0 * log(z1)) * sin(2.0 * PI * z2);

    *x1 = stdv * (*x1) + mean;
    *x2 = stdv * (*x2) + mean;
}

double randn(double mean, double stdv) {
    double x;

    double z1 = (double) rand() / RAND_MAX;
    double z2 = (double) rand() / RAND_MAX;

    x = sqrt(-2.0 * log(z1)) * cos(2.0 * PI * z2);
    x = stdv * x + mean;

    return x;
}
