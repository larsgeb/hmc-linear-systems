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

/*====================================================================================================================*/
/* Inversion and model parameters class. -----------------------------------------------------------------------------*/
/*====================================================================================================================*/


/* Constructor. ------------------------------------------------------------------------------------------------------*/
parameters::parameters() {

}

/* Destructor. -------------------------------------------------------------------------------------------------------*/
parameters::~parameters() {

}

/* Copy constructor. -------------------------------------------------------------------------------------------------*/
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
    // add values (in q), takes mean & sigma of a
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
}

/* Calculate inverse covariance matrix -------------------------------------------------------------------------------*/
void parameters::calculateInverseCM() {
    // Note that only uncorrelated prior is as of yet implemented (unit matrices)
    for (int i = 0; i < Nq; i++) {
        std::vector<double> row_iCM;
        // Generating the row for an identity matrix
        for (int j = 0; j < Nq; j++) {
            if (i == j) {
                row_iCM.push_back(1 / sigma_q[i]);
            } else {
                row_iCM.push_back(0);
            }
        }
        // Pushing the row into the matrix
        iCM.push_back(row_iCM);
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
    // q is filled in the following way: 0 - nLayers: speeds || nLayers+1 - nLayers*2: thicknesses
    std::vector<double> travelTime;

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
    double a;
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

double data::misfit(parameters *q) {
    // This function only incorporates prior and forward model misfit. Added could be forward model errors.
    double Xi;

    // Maybe implement some extra functions for vector arithmetic
    std::vector<double> q_difference;
    q_difference = VectorDifference(q->q, q->mean_q);
    Xi = 0.5 * VectorVectorProduct(q_difference, MatrixVectorProduct(q->iCM, q_difference));

    std::vector<double> d_difference;
    std::vector<double> d_forward;
    d_forward = forwardModel(*q);
    d_difference = VectorDifference(d_forward, recT);
    d_forward.clear();

    return Xi;
}

std::vector<double> data::misfitT1(parameters *q, double ratioStep) {
    std::vector<double> firstDerivative;

    // Doing a Taylor expansion of the first order. This evaluates the derivatives at the starting model,
    // but note that it still would have to be multiplied by (m - m0) for the taylor function
    for (int i = 0; i < q->Nq; i++) {
        double current_q;
        current_q = q->q[i];
        parameters q_perturbed;
        q_perturbed.read_input("INPUT/parameters.txt");
        q_perturbed.q[i] = current_q*ratioStep;
        firstDerivative.push_back((q_perturbed.t0U - q->t0U)/(current_q*ratioStep));
    }
    return firstDerivative;
}

std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B) {

    std::vector<double> C;

    if (A.size() != B.size()) {
        // Get some Exception class to THROW.
        std::cout << "Vectors are not the same dimension! The code DIDN'T run successfully.";
        return C;
    }

    std::vector<double> q_difference;
    for (int i = 0; i < A.size(); i++) {
        // Prior misfit
        C.push_back(A[i] - B[i]);
    }
    return C;
}

std::vector<double> MatrixVectorProduct(std::vector<std::vector<double>> M, std::vector<double> A) {
    std::vector<double> C;

    // Using std::vector<>.size() requires casting for clean compilation (seems unnecessary.. But oh well)
    // So watch out if you're working on 2^63 order problems.. ;) (Maximum <int>, not maximum <unsigned int>)
    int rowsM = static_cast<int>(M.size());
    int columnsM = static_cast<int>(M[0].size());
    int rowsA = static_cast<int>(A.size());

    if (columnsM != rowsA) {
        // Get some Exception class to THROW.
        std::cout << "Vector and matrix are not compatible in dimension! The code DIDN'T run successfully.";
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
        std::cout << "Vectors are not compatible in dimension! The code DIDN'T run successfully.";
        return C;
    }
    C = 0;
    for (int i = 0; i < A.size(); i++) {
        C += (A[i] * B[i]);
    }
    return C;
}

