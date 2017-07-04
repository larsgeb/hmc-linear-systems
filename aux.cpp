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

/*=================================================================================*/
/* Inversion and model parameters class. ------------------------------------------*/
/*=================================================================================*/

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
}

/* Assignment operator. -----------------------------------------------------------*/
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

/* Addition. ----------------------------------------------------------------------*/
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

/* Subtraction. -------------------------------------------------------------------*/
parameters operator-(const parameters &a, const parameters &b) {
    parameters c;

    for (int i = 0; i < a.q.size(); i++) {
        c.q.push_back(a.q[i] - b.q[i]);
        c.mean_q.push_back(a.mean_q[i]);
        c.sigma_q.push_back(a.sigma_q[i]);
    }

    return c;
}


/* Read input from file. ----------------------------------------------------------*/
void parameters::read_input(const char *filename) {
    FILE *fid;
    char str[1000];

    fid = fopen(filename, "r");

    fgets(str, 1000, fid); // Move line
    fscanf(fid, "%d", &nLayers);
    fgets(str, 1000, fid);

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

    fclose(fid);
}

void parameters::calculateLayerBase() {
    for (int i = 0; i < nLayers; i++) {
        if (i == 0) {
            layerBase.push_back(q[nLayers + 1 + i]);
        } else {
            layerBase.push_back(layerBase[i - 1] + q[nLayers + 1 + i]);
        }
    }
}

/*=================================================================================*/
/* VSP data class. ----------------------------------------------------------------*/
/*=================================================================================*/

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
                // checking whether receiver is completely below layers:
//                std::cout << "Receiver " << i + 1 << " at x=" << recZ[i] << " is below  layer " << j + 1
//                          << " with top " << top << " and base " << base << std::endl;
                travelTime[i] += thickness / q.q[j];
            } else if (recZ[i] > top) {
//                std::cout << "Receiver " << i + 1 << " at x=" << recZ[i] << " is within layer " << j + 1 << std::endl;
                travelTime[i] += (recZ[i] - top) / q.q[j];
            }
        }
        if (recZ[i] > q.layerBase[q.nLayers-1]) {
//            std::cout << "Receiver " << i + 1 << " at x=" << recZ[i] << " is in the lower halfspace" << std::endl;
            travelTime[i] += (recZ[i] - q.layerBase[q.nLayers-1]) / q.q[q.nLayers];
        }
    }
    return travelTime;
}

void data::read_data(const char *filename) {
    double a;
    std::ifstream infile(filename);
    for (int i = 0; i < recN; i++){
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






