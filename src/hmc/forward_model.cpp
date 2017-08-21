//
// Created by Lars Gebraad on 18-8-17.
//

#include "forward_model.hpp"

using namespace algebra_lib;

hmc::forward_model::forward_model(const matrix &_model_matrix)
        : _g(_model_matrix) {

}
