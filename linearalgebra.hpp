/** @file linearalgebra.hpp
 * @Author Lars Gebraad (larsgebraad@gmail.com)
 * @date   September, 2008
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
 * When I refer to matrices or vectors, the actual C++ objects will be these containers (of containers) of doubles.
 * The first elements have index zero.
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
 * @param M Any \f$ n \times m \f$ matrix \f$M\f$.
 * @param A Any \f$ m \f$ dimensional vector A.
 * @return A \f$ n \f$ dimensional vector B which is defined through: \f$ B_i = \sum_{k=1}^m M_{ik} A_k \f$.
 */
std::vector<double> MatrixVectorProduct(std::vector<std::vector<double>> M, std::vector<double> A);

/**
 * @brief Function incorporating the standard matrix-matrix product, producing a new matrix. Matrix M should have as many columns as N has rows, otherwise an exception is thrown.
 * @param M Any \f$ n \times m \f$ matrix \f$M\f$.
 * @param N Any \f$ m \times p \f$ matrix \f$N\f$.
 * @return A \f$ n \times p \f$ matrix P defined by \f$ P_{ij} = \sum_{k=1}^m M_{ik} N_{kj} \f$.
 */
std::vector<std::vector<double>> MatrixMatrixProduct(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);

/**
 * @brief A function to calculate the sum of the entries of two matrices.
 * @param M Any \f$ n \times m \f$ matrix \f$M\f$.
 * @param N Any \f$ n \times m \f$ matrix \f$N\f$.
 * @return A \f$ n \times m \f$ matrix \f$ S \f$ containing the sum of the entries of \f$M\f$ and \f$N\f$.
 */
std::vector<std::vector<double>> MatrixMatrixSum(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);

/**
 * @brief Function to get a row from a matrix.
 * @param M Any \f$ n \times m \f$ matrix \f$M\f$.
 * @param row Integer indicating the row. Numbering starts at 0.
 * @return An \f$ m \f$ dimensional vector containing the extracted matrix row.
 */
std::vector<double> GetMatrixRow(std::vector<std::vector<double>> M, int row);

/**
 * @brief Function to get a column from a matrix.
 * @param M Any \f$ n \times m \f$ matrix \f$M\f$.
 * @param row Integer indicating the column. Numbering starts at 0.
 * @return An \f$ n \f$ dimensional vector containing the extracted matrix column.
 */
std::vector<double> GetMatrixColumn(std::vector<std::vector<double>> M, int column);

/**
 * @brief Function to transpose a matrix M of size i x j into a matrix N of size j x i, where \f$M_{ij} = N_{ji}\f$.
 * @param M Any \f$ n \times m \f$ matrix \f$M\f$.
 * @return The transpose of \f$M\f$.
 */
std::vector<std::vector<double>> TransposeMatrix(std::vector<std::vector<double>> M);

/**
 * @brief Function to the trace of a square \f$n \times n\f$ matrix \f$M\f$.
 * @param M Any square matrix.
 * @return The trace of \f$ M \f$.
 */
std::vector<double> MatrixTrace(std::vector<std::vector<double>> M);

/**
 * @brief Function to take the inverse of the individual matrix elements
 * @param M Any \f$ n \times m \f$ matrix \f$M\f$.
 * @return The inverted elements \f$ 1/M_{ij}\f$ in a std::vector of std::vector of double.
 */
std::vector<std::vector<double>> InvertMatrixElements(std::vector<std::vector<double>> M);

/**
 * @brief Function which takes a std::vector of double to make a diagonal matrix of it, such that \f$ A_{i} = M_{ii} \f$.
 * @param A Any Any \f$ n \f$ dimensional vector.
 * @return Diagonal \f$ n \times n \f$ matrix \f$M\f$.
 */
std::vector<std::vector<double>> VectorToDiagonal(std::vector<double> A);

/**
 * @brief Dot product of vectors.
 * @param A Vector \f$ A \f$  of dimension \f$ n\f$.
 * @param B Vector \f$ B \f$  of dimension \f$ n\f$.
 * @return Dot product \f$ c \f$ of \f$ A \f$ and \f$ B \f$ where where \f$ c = A_i \cdot B_i \f$ (float).
 */
double VectorVectorProduct(std::vector<double> A, std::vector<double> B);

/**
 * @brief Vector difference between two vectors.
 * @param A Vector \f$ A \f$ of dimension \f$ n\f$.
 * @param B Vector \f$ B \f$ of dimension \f$ n\f$.
 * @return Vector \f$ C \f$ where \f$ C_i = A_i - B_i \f$.
 */
std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B);

/**
 * @brief Vector sum between two vectors.
 * @param A Vector \f$ A \f$ of dimension \f$ n\f$.
 * @param B Vector \f$ B \f$ of dimension \f$ n\f$.
 * @return Vector \f$ C \f$ where \f$ C_i = A_{i} + B_i \f$.
 */
std::vector<double> VectorSum(std::vector<double> A, std::vector<double> B);

/**
 * @brief Vector scalar prodcut of vector and scalar.
 * @param A Any vector \f$ A \f$.
 * @param B Scalar \f$ b \f$.
 * @return Vector \f$ C \f$ where \f$ C_{i} = A_{i} \cdot b \f$.
 */
std::vector<double> VectorScalarProduct(std::vector<double> A, double b);

/**
 * @brief Solve linear equation \f$ L y = x\f$ where \f$ L \f$ is a lower triangular \f$ n \times n\f$ matrix and \f$ x
 * \f$ and \f$ y \f$ are \f$ n\f$ dimensional vectors. Uses forward and backward substitution to iteratively solve the
 * system.
 * @param L Any lower triangular \f$ n \times n\f$ matrix.
 * @param x Any \f$ n\f$ dimensional vector.
 * @return Solution y of the system, an \f$ n\f$ dimensional vector.
 */
std::vector<double> SolveLowerTriangular(std::vector<std::vector<double>> L, std::vector<double> x);

/**
 * @brief Invert a lower triangular \f$ n \times n\f$ matrix by use of solving the system \f$ L L^{-1} = I\f$ per column
 * of \f$ I \f$ using SolveLowerTriangular(), the identity matrix.
 * @param L Any lower triangular \f$ n \times n\f$ matrix.
 * @return The inverse of \f$ L \f$, another lower triangular \f$ n \times n\f$ matrix.
 */
std::vector<std::vector<double>> InvertLowerTriangular(std::vector<std::vector<double>> L);

/**
 * @brief Cholesky-decomposition of a positive definite Hermitian \f$ n \times n\f$ matrix \f$ M \f$.
 * @param A A positive definite Hermitian \f$ n \times n\f$ matrix.
 * @return A lower triangular \f$ n \times n\f$ matrix for which holds \f$ L L^\dagger = M \f$, where \f$ L^\dagger \f$ is
 * the adjoint of \f$ L \f$ (which of course simplifies to transposition for real matrices).
 */
std::vector<std::vector<double>> CholeskyDecompose(std::vector<std::vector<double>> A);

/* ---------------------------------------------------------------------------------------------------------------------- *
 * Operators which wrap to the above functions. The trick in getting these to use less memory is std::move.
 * ---------------------------------------------------------------------------------------------------------------------- */
/**
 * @brief Operator form of MatrixMatrixProduct(), using std library forwarding.
 */
std::vector<std::vector<double>> operator*(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);
/**
 * @brief Operator form of VectorVectorProduct(), using std library forwarding.
 */
double operator*(std::vector<double> A, std::vector<double> B);
/**
 * @brief Operator form of MatrixVectorProduct(), using std library forwarding.
 */
std::vector<double> operator*(std::vector<std::vector<double>> M, std::vector<double> A);
/**
 * @brief Operator form of VectorScalarProduct(), using std library forwarding.
 */
std::vector<double> operator*(double b, std::vector<double> A);
/**
 * @brief Operator form of VectorScalarProduct(), using std library forwarding.
 */
std::vector<double> operator*(std::vector<double> A, double b);
/**
 * @brief Operator form of VectorSum(), using std library forwarding.
 */
std::vector<double> operator+(std::vector<double> A, std::vector<double> B);
/**
 * @brief Operator form of VectorDifference(), using std library forwarding.
 */
std::vector<double> operator-(std::vector<double> A, std::vector<double> B);
/**
 * @brief Operator form of MatrixMatrixSum(), using std library forwarding.
 */
std::vector<std::vector<double>> operator+(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N);

#endif //HMC_VSP_LINEARALGEBRA_HPP
