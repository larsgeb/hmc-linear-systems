//
// Created by Lars Gebraad on 7/10/17.
//

#ifndef HMC_VSP_LINEARALGEBRA_HPP
#define HMC_VSP_LINEARALGEBRA_HPP

std::vector<double> MatrixVectorProduct(std::vector<std::vector<double>> M, std::vector<double> A);

std::vector<std::vector<double>> MatrixMatrixProduct(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);

std::vector<std::vector<double>> MatrixMatrixSum(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);

std::vector<double> GetMatrixRow(std::vector<std::vector<double>> M, int row);

std::vector<double> GetMatrixColumn(std::vector<std::vector<double>> M, int column);

std::vector<std::vector<double>> TransposeMatrix(std::vector<std::vector<double>> M);

std::vector<double> MatrixTrace(std::vector<std::vector<double>> M);

std::vector<std::vector<double>> InvertMatrix(std::vector<std::vector<double>> M);

std::vector<std::vector<double>> VectorToDiagonal(std::vector<double>);

double VectorVectorProduct(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorSum(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorScalarProduct(std::vector<double> A, double b);

/* ---------------------------------------------------------------------------------------------------------------------- *
 * Operators which wrap to the above functions. The trick in getting these to use less memory is std::move.
 * ---------------------------------------------------------------------------------------------------------------------- */
// Matrix multiplication
std::vector<std::vector<double>> operator*(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);
// In-product of two vectors
double operator*(std::vector<double> A, std::vector<double> B);
// Matrix-vector multiplication
std::vector<double> operator*(std::vector<std::vector<double>> M, std::vector<double> A);
// Two identical scalar-vector multiplcations
std::vector<double> operator*(double b, std::vector<double> A);
std::vector<double> operator*(std::vector<double> A, double b);
// Vector sum
std::vector<double> operator+(std::vector<double> A, std::vector<double> B);
// Vector difference
std::vector<double> operator-(std::vector<double> A, std::vector<double> B);
// Matrix sum
std::vector<std::vector<double>> operator+(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);

#endif //HMC_VSP_LINEARALGEBRA_HPP
