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
    data d; // Already loads setup.xt

    q.read_input("INPUT/parameters.txt"); // Load subsurface parameters into q
    d.make_synthetics(q);
    d.write("DATA/synthetics.txt");

    q.tExpand(d,1,1.000001); //1.000001 seems sufficient for working precision. (~15 digits)

    return 0;
}