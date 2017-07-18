//
// Created by Lars Gebraad on 7/10/17.
//

#ifndef HMC_VSP_LINEARALGEBRA_HPP
#define HMC_VSP_LINEARALGEBRA_HPP

std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorSum(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorScalarProduct(std::vector<double> A, double b);

std::vector<double> MatrixVectorProduct(std::vector<std::vector<double>> M, std::vector<double> A);

std::vector<std::vector<double>> MatrixMatrixProduct(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);

std::vector<double> GetMatrixRow(std::vector<std::vector<double>> M, int row);

std::vector<std::vector<double>> TransposeMatrix(std::vector<std::vector<double>> M);

std::vector<double> GetMatrixColumn(std::vector<std::vector<double>> M, int column);

double VectorVectorProduct(std::vector<double> A, std::vector<double> B);

#endif //HMC_VSP_LINEARALGEBRA_HPP
