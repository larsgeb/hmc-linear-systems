//
// Created by Lars Gebraad on 7/10/17.
//

#ifndef HMC_VSP_LINEARALGEBRA_HPP
#define HMC_VSP_LINEARALGEBRA_HPP

std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorSum(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorScalarProduct(std::vector<double> A, double b);

std::vector<double>
MatrixVectorProduct(std::vector<std::vector<double>> M, std::vector<double> A);

double VectorVectorProduct(std::vector<double> A, std::vector<double> B);

#endif //HMC_VSP_LINEARALGEBRA_HPP
