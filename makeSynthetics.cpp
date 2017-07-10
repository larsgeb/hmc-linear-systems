//
// Created by Lars Gebraad on 7/10/17.
// Creates synthetic data based on setup.txt (for receivers) and parameters_synthetics_model.txt (for model).
//
#include <vector>
#include "auxiliary.hpp"

int main() {
    // Easy as can be.
    data synthetics; // Needed for location of receivers, hardcoded into class
    forwardVSP::generateSynthetics("INPUT/parameters_synthetics_model.txt",synthetics._depthReceivers);
    return 0;
}