//
// Created by Lars Gebraad on 18-8-17.
//

#include "forward_model.hpp"

using namespace arma;

hmc::forward_model::forward_model(const mat &_model_matrix)
        : _g(_model_matrix) {
}

hmc::forward_model::forward_model() {}

