//
// Created by Lars Gebraad on 18-8-17.
//

#include "ForwardModel.hpp"

HMC::ForwardModel::ForwardModel(const AlgebraLib::Matrix &_model_matrix)
        : _g(_model_matrix) {

}
