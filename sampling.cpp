//
// Created by Lars Gebraad on 7/4/17.
//

#include <iostream>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "aux.hpp"
#include "mc.hpp"

int main(int argc, char *argv[]) {

    int nt, iterations;
    double dt;
    bool verbose;

    if (!strcmp(argv[1], "metropolis")) {
        nt = 0;
        dt = 0.0;
        iterations = atoi(argv[2]);
        if (argv[3] && !strcmp(argv[3], "verbose")) verbose = true;
    } else if (!strcmp(argv[1], "hamilton")) {
        iterations = atoi(argv[2]);
        nt = atoi(argv[3]);
        dt = atof(argv[4]);
        if (argv[5] && !strcmp(argv[5], "verbose")) verbose = true;
    } else {
        std::cout << "No valid sampling method is specified, terminating." << std::endl;
        return 0;
    }

    parameters parameters1;
    data data1;
    parameters1.read_input("INPUT/parameters.txt"); // Load subsurface parameters into q
    data1.read_data("DATA/synthetics.txt");  //

    mc m(iterations, nt, dt, verbose);
    clock_t start = clock();
    int accepted = 0;
    double x, x_new;

    /* Initial values. ----------------------------------------------------------------*/

    if (!strcmp(argv[1], "metropolis")) x = m.chi();
    else if (!strcmp(argv[1], "hamilton")) x = m.energy();

    m.write_sample(pfile, x, 0);


    return 0;
}