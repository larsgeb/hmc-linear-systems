//
// Created by lars on 11/04/18.
//

#include <armadillo>
#include "../hmc/sampler.hpp"

using namespace arma;

int main(int argc, char *argv[]) {

    mat G = zeros(4, 4);

    G(0, 0) = 1;
    G(0, 3) = 1;
    G(1, 1) = 1;
    G(2, 2) = 1;
    G(3, 2) = 1;

    vec m = {1, 2, 3, 4};

    vec d_obs = G * m;

    d_obs.save("observed_data.bin");
    G.save("model_matrix.bin");

    return EXIT_SUCCESS;
}
