//
// Created by lars on 11/04/18.
//

#include <armadillo>
#include "../hmc/sampler.hpp"

using namespace arma;

int main(int argc, char *argv[]) {

    mat G = mat(10000, 10000, fill::eye);

    G(1, 1) = 5;
    G(0, 1) = 5;

    vec m = vec(10000, fill::ones);

    m(5) = 2;
    
    vec d_obs = G * m;

    d_obs.save("observed_data.bin");
    G.save("model_matrix.bin");
    G.save("model_matrix.txt", raw_ascii); // For plotting the matrix using matplotlib

    return EXIT_SUCCESS;
}
