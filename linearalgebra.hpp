//
// Created by Lars Gebraad on 7/10/17.
//

#ifndef HMC_VSP_LINEARALGEBRA_HPP
#define HMC_VSP_LINEARALGEBRA_HPP

/**
 * @brief Function incorporating the standard matrix-vector product. The result is a column vector.
 * @param M A std::vector of a std::vector of doubles, containing the matrix row elements (inner dimension) in the matrix rows (outer dimension).
 * @param A std::vector<double> object containing the elements of the column vector.
 * @return std::vector<double> object containing the elements of the matrix-vector product.
 */
std::vector<double> MatrixVectorProduct(std::vector<std::vector<double>> M, std::vector<double> A);
/**
 * @brief Function incorporating the standard matrix-vector product, producing a new matrix. Matrix M should have as many columns as N has rows, otherwise an exception is thrown.
 * @param M A std::vector of a std::vector of doubles, containing the matrix row elements (inner dimension) in the matrix rows (outer dimension) of the left matrix. Dimensions are a x b.
 * @param N A std::vector of a std::vector of doubles, containing the matrix row elements (inner dimension) in the matrix rows (outer dimension) of the right matrix. Dimensions are b x c.
 * @return A std::vector of a std::vector of doubles, containing the matrix row elements (inner dimension) in the matrix rows (outer dimension) of the resulting matrix. Dimensions are a x c.
 */
std::vector<std::vector<double>> MatrixMatrixProduct(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);
/**
 * @brief Function incorporating the standard matrix-vector product, producing a new matrix. Matrix M should have as many columns as N has rows, otherwise an exception is thrown.
 * @param M
 * @param N
 * @return
 */
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

std::vector<std::vector<double>> CholeskyDecompose(std::vector<std::vector<double>> A);

/* ---------------------------------------------------------------------------------------------------------------------- *
 * Operators which wrap to the above functions. The trick in getting these to use less memory is std::move.
 * ---------------------------------------------------------------------------------------------------------------------- */
// Matrix multiplication
std::vector<std::vector<double>> operator*(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);
// In-product of two vectors
double operator*(std::vector<double> A, std::vector<double> B);
// Matrix-vector multiplication
std::vector<double> operator*(std::vector<std::vector<double>> M, std::vector<double> A);
// Two identical scalar-vector multiplications
std::vector<double> operator*(double b, std::vector<double> A);
std::vector<double> operator*(std::vector<double> A, double b);
// Vector sum
std::vector<double> operator+(std::vector<double> A, std::vector<double> B);
// Vector difference
std::vector<double> operator-(std::vector<double> A, std::vector<double> B);
// Matrix sum
std::vector<std::vector<double>> operator+(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);

#endif //HMC_VSP_LINEARALGEBRA_HPP
