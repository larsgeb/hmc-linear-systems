//
// Created by Lars Gebraad on 18-8-17.
//

#include "forwardmodel.hpp"

hmc::ForwardModel::ForwardModel(const algebra_lib::matrix &_model_matrix)
        : _g(_model_matrix) {

}
