//
// Created by Lars Gebraad on 7/19/17.
//

#include <fstream>
#include "tomographyForwardModel.hpp"
#include "linearalgebra.hpp"

/* -----------------------------------------------------------------------------------------------------------------------
 * Forward model class.
 * ----------------------------------------------------------------------------------------------------------------------- */
void forwardModel::constructDesignMatrix(int numberParameters) {

    // Loading model properties
    std::ifstream matrixPropertiesFile("straight-ray-solver/forward_matrix_properties.txt");
    matrixPropertiesFile >> _numberReceivers;
    matrixPropertiesFile >> _numberSources;
    matrixPropertiesFile >> _numberCellsX;
    matrixPropertiesFile >> _numberCellsY;
    matrixPropertiesFile.close();
    _numberParameters = _numberCellsX * _numberCellsY;
    _numberData = _numberSources * _numberReceivers;

    // Read file for observed data
    double b;
    std::ifstream infile("straight-ray-solver/forward_matrix.txt");
    for (int row = 0; row < _numberData; row++) {
        _designMatrix.push_back(std::vector<double> {});
        for (int column = 0; column < _numberParameters; column++) {
            infile >> b;
            _designMatrix[row].push_back(b);
        }
    }
    infile.close();
}

forwardModel::forwardModel(int numberParameters) {
    _numberParameters = numberParameters;
    constructDesignMatrix(numberParameters);
}

std::vector<double> forwardModel::calculateData(std::vector<double> parameters) {
    return MatrixVectorProduct(_designMatrix, parameters);
}

forwardModel::forwardModel() {}