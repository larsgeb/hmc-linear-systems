//
// Created by Lars Gebraad on 7/19/17.
//

#ifndef HMC_VSP_FORWARDMODEL_HPP
#define HMC_VSP_FORWARDMODEL_HPP

#include <vector>

class forwardModel;

class forwardModel {
public:
    // Constructors & destructors
    forwardModel(int numberParameters);

    forwardModel();

    // Fields
    int _numberParameters;
    int _numberData;
    int _numberSources;
    int _numberReceivers;
    int _numberCellsX;
    int _numberCellsY;
    std::vector <std::vector<double>> _designMatrix;

    // Methods
    void constructDesignMatrix(int numberParameters);

    std::vector<double> calculateData(std::vector<double> parameters);
};

#endif //HMC_VSP_FORWARDMODEL_HPP
