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
#include "montecarlo.hpp"

int main(int argc, char *argv[]) {

    int nt, iterations;
    double dt;
    bool verbose = false;

    if (!strcmp(argv[1], "metropolis")) {
        nt = 0;
        dt = 0.0;
        iterations = atoi(argv[2]);
        if (argv[3] && !strcmp(argv[3], "m_verbose")) verbose = true;
    } else if (!strcmp(argv[1], "hamilton")) {
        iterations = atoi(argv[2]);
        nt = atoi(argv[3]);
        dt = atof(argv[4]);
        if (argv[5] && !strcmp(argv[5], "m_verbose")) verbose = true;
    } else {
        std::cout << "No valid sampling method is specified, terminating."
                  << std::endl;
        return 0;
    }

    parameters start_model;
    data data1;
    start_model.read_input(
            "INPUT/parameters_starting_model.txt"); // Load subsurface parameters into m_q
    data1.read_data("DATA/synthetics.txt");

    start_model.tExpand(data1, 1,
                        1.000001); // Do a Taylor expansion of the misfit function to avoid MANY calculations

    montecarlo m(iterations, nt, dt, start_model.Nq, verbose);
    clock_t start = clock();
    int accepted = 0;
    double x, x_new;

    FILE *pfile;
    pfile = fopen("OUTPUT/samples.txt", "w");

    /* Initial values. ----------------------------------------------------------------*/
    if (!strcmp(argv[1], "metropolis")) x = m.chi();
    else if (!strcmp(argv[1], "hamilton")) x = m.energy();

    /* Random walk. -------------------------------------------------------------------*/
    for (int it = 0; it < m.m_iterations; it++) {
        /* Make a model proposition and compute misfit/energy. */
        if (!strcmp(argv[1], "metropolis")) {
            m.propose_metropolis();
            x_new = m.chi();
        } else if (!strcmp(argv[1], "hamilton")) {
            m.propose_hamilton();
            x_new = m.energy();
        }

        /* Check Metropolis rule. */
        if ((x_new < x) || (exp(x - x_new) > randf(0.0, 1.0))) {
            x = x_new;
            m.m_q = m.m_q_new;
            accepted++;
        }

    }

    printf("accepted: %d\n", accepted);
    printf("elapsed time: %f\n", (double) (clock() - start) / CLOCKS_PER_SEC);

    /* Clean up. ----------------------------------------------------------------------*/
    fclose(pfile);
    return 0;
}