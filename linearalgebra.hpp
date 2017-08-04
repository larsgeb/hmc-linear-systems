//
// Created by Lars Gebraad on 7/10/17.
//
/** @file
 *
 * @brief Linear algebra functions operating on standard library containers.
 *
 * A kind of library for all kinds of linear algebra functions. Incorporates basic operations such as multiplication and
 * dot products (with respective operator overloading) as well as some more specialized functions such as Cholesky
 * decomposition and lower triangular matrix inversion.
 *
 * This set of function operate on containers defined in the standard library. Most of the objects passed as input or
 * output within this set of functions are std::vector's, of either std::vector<doubles> or doubles itself. These
 * containers act as vectors and matrices. The benefit of these is that they are memory optimized, allow for (nearly)
 * perfect forwarding and other easy manipulation.
 *
 * Note that if there are nested containers, the convention is that the
 * first index (the outer container) represents rows, where the inner containers represent row elements. There is no easy
 * way to directly access a column, except when usen GetMatrixColumn().
 *
 * */
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
 * @param M std::vector of a std::vector of doubles containing the matrix, dimensions i x j.
 * @param N std::vector of doubles containing the product vector, length j.
 * @return std::vector of doubles containing the result vector, length i.
 */
std::vector<std::vector<double>> MatrixMatrixSum(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);
/**
 * @brief Function to get a row from a matrix.
 * @param M std::vector of a std::vector of doubles, which is the matrix from which to extract the row.
 * @param row int indicating the row. Numbering starts at 0.
 * @return std::vector of doubles containing the extracted matrix row.
 */
std::vector<double> GetMatrixRow(std::vector<std::vector<double>> M, int row);
/**
 * @brief Function to get a column from a matrix.
 * @param M std::vector of a std::vector of doubles, which is the matrix from which to extract the row.
 * @param row int indicating the column. Numbering starts at 0.
 * @return std::vector of doubles containing the extracted matrix column.
 */
std::vector<double> GetMatrixColumn(std::vector<std::vector<double>> M, int column);
/**
 * @brief Function to transpose a matrix M of size i x j into a matrix N of size j x i, where \f$M_{ij} = N_{ji}\f$
 * @param M Any std::vector of std::vector of double.
 * @return The transpose of M.
 */
std::vector<std::vector<double>> TransposeMatrix(std::vector<std::vector<double>> M);

std::vector<double> MatrixTrace(std::vector<std::vector<double>> M);

std::vector<std::vector<double>> InvertMatrixElements(std::vector<std::vector<double>> M);

std::vector<std::vector<double>> VectorToDiagonal(std::vector<double>);

double VectorVectorProduct(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorSum(std::vector<double> A, std::vector<double> B);

std::vector<double> VectorScalarProduct(std::vector<double> A, double b);

std::vector<double> SolveLowerTriangular(std::vector<std::vector<double>> L, std::vector<double> x);

std::vector<std::vector<double>> InvertLowerTriangular(std::vector<std::vector<double>> L);

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
