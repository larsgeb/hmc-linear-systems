#include <iostream>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "aux.hpp"

int main() {

    parameters q;
    parameters q_start;
    data d; // Already loads setup.xt

    q.read_input("INPUT/parameters_synthetics_model.txt"); // Load subsurface parameters into q
    q_start.read_input("INPUT/parameters_starting_model.txt"); // Load subsurface parameters into q
    d.make_synthetics(q);
    d.write("DATA/synthetics.txt");
    d.read_data("DATA/synthetics.txt");
    q_start.tExpand(d, 1, 1.000001); //1.000001 seems sufficient for working precision. (~15 digits)
    return 0;
}