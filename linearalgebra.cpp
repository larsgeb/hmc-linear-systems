//
// Created by Lars Gebraad on 7/10/17.
// Simple Linear Algebra stuff.
//


#include <vector>
#include <iostream>
#include <cmath>


// Linear algebra functions
std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B) {

    std::vector<double> C;

    if (A.size() != B.size()) {
        // Get some Exception class to THROW.
        std::cout << "Vectors are not the same dimension! The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }

    std::vector<double> q_difference;
    for (int i = 0; i < A.size(); i++) {
        // Prior misfit
        C.push_back(A[i] - B[i]);
    }
    return C;
}

std::vector<double> VectorScalarProduct(std::vector<double> A, double b) {
    for (double &i : A) {
        i = i * b;
    }
    return A;
}

std::vector<double> VectorSum(std::vector<double> A, std::vector<double> B) {

    std::vector<double> C;

    if (A.size() != B.size()) {
        // Get some Exception class to THROW.
        std::cout << "Vectors are not the same dimension! The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }

    std::vector<double> q_difference;
    for (int i = 0; i < A.size(); i++) {
        // Prior misfit
        C.push_back(A[i] + B[i]);
    }
    return C;
}

std::vector<double> MatrixVectorProduct(std::vector<std::vector<double> > M, std::vector<double> A) {
    std::vector<double> C;

    // Using std::vector<>.size() requires casting for clean compilation (seems unnecessary.. But oh well)
    // So watch out if you're working on 2^63 order problems.. ;) (Maximum <int>, not maximum <unsigned int>)
    int rowsM = static_cast<int>(M.size());
    int columnsM = static_cast<int>(M[0].size());
    int rowsA = static_cast<int>(A.size());

    if (columnsM != rowsA) {
        // Get some Exception class to THROW.
        std::cout << "Vector and matrix are not compatible in dimension! The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }

    for (int i = 0; i < rowsM; i++) {
        C.push_back(0);
        for (int j = 0; j < columnsM; j++) {
            C[i] += M[i][j] * A[j];
        }
    }
    return C;
}

double VectorVectorProduct(std::vector<double> A, std::vector<double> B) {
    double C;

    if (A.size() != B.size()) {
        // Get some Exception class to THROW.
        std::cout << "Vectors are not compatible in dimension! The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }
    C = 0;
    for (int i = 0; i < A.size(); i++) {
        C += (A[i] * B[i]);
    }
    return C;
}

std::vector<double> GetMatrixColumn(std::vector<std::vector<double>> M, int column) {
    std::vector<double> A;
    unsigned long rows = M.size();

    if (column > M[0].size()) {
        // Get some Exception class to THROW.
        std::cout << "Matrix index out of range. The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }

    A.clear();
    for (unsigned long row = 0; row < rows; row++) {
        A.push_back(M[row][column]);
    }
    return A;
}

std::vector<double> GetMatrixRow(std::vector<std::vector<double>> M, int row) {
    std::vector<double> A;
    unsigned long columns = M[0].size();

    if (row > M.size()) {
        // Get some Exception class to THROW.
        std::cout << "Matrix index out of range. The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }

    A.clear();
    for (unsigned long column = 0; column < columns; column++) {
        A.push_back(M[row][column]);
    }
    return A;
}

std::vector<std::vector<double>> TransposeMatrix(std::vector<std::vector<double> > M) {
    unsigned long rows = M.size();
    unsigned long columns = M[0].size();
    std::vector<std::vector<double>> N;

    N.clear();
    std::vector<double> zeroRow(rows, 0);
    N.insert(N.end(), columns, zeroRow);

    for (int row = 0; row < rows; row++) {
        for (int column = 0; column < columns; column++) {
            N[column][row] = M[row][column];
        }
    }
    return N;
}

std::vector<double> MatrixTrace(std::vector<std::vector<double>> M) {
    if (M.size() != M[0].size()) {
        // Get some Exception class to THROW.
        std::cout << "Matrix is not square, and trace doesn't exist. The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }
    std::vector<double> trace;
    trace.clear();
    for (int i = 0; i < M.size(); i++) {
        trace.push_back(M[i][i]);
    }
    return trace;
}

std::vector<std::vector<double>> InvertMatrix(std::vector<std::vector<double>> M) {
    for (int i = 0; i < M.size(); i++) {
        for (int j = 0; j < M[i].size(); j++) {
            if (M[i][j] != 0) {
                M[i][j] = 1 / M[i][j];
            }
        }
    }
    return M;
}

std::vector<std::vector<double>> VectorToDiagonal(std::vector<double> A) {
    std::vector<std::vector<double>> M;

    M.clear();
    std::vector<double> zeroRow(A.size(), 0);
    M.insert(M.end(), A.size(), zeroRow);

    for (int i = 0; i < A.size(); i++) {
        M[i][i] = A[i];
    }

    return M;
}

std::vector<std::vector<double>>
MatrixMatrixProduct(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N) {
    std::vector<std::vector<double>> P;
    unsigned long columnsM = M[0].size();
    unsigned long rowsM = M.size();
    unsigned long columnsN = N[0].size();
    unsigned long rowsN = N.size();

    if (columnsM != rowsN) {
        // Get some Exception class to THROW.
        std::cout << "Matrices are not compatible in dimension. The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }

    P.clear();
    std::vector<double> zeroRow(columnsN, 0);
    P.insert(P.end(), rowsM, zeroRow);

    for (int rowM = 0; rowM < rowsM; rowM++) {
        for (int columnN = 0; columnN < columnsN; columnN++) {
            double sum = 0;
            for (int innerDim = 0; innerDim < rowsN; innerDim++) {
                sum += M[rowM][innerDim] * N[innerDim][columnN];
            }
            P[rowM][columnN] = sum;
        }
    }
    return P;
}

std::vector<std::vector<double>> MatrixMatrixSum(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N) {
    std::vector<std::vector<double>> S;

    unsigned long columnsM = M[0].size();
    unsigned long rowsM = M.size();
    unsigned long columnsN = N[0].size();
    unsigned long rowsN = N.size();

    if (columnsM != columnsN || rowsM != rowsN) {
        // Get some Exception class to THROW.
        std::cout << "Matrices are not compatible in dimension. The code DIDN'T run successfully." << std::endl;
        throw std::exception();
    }

    S.clear();
    std::vector<double> zeroRow(columnsN, 0);
    S.insert(S.end(), rowsN, zeroRow);

    for (int rowM = 0; rowM < rowsM; rowM++) {
        for (int columnM = 0; columnM < columnsM; columnM++) {
            S[rowM][columnM] = M[rowM][columnM] + N[rowM][columnM];
        }
    }

    return S;
}

std::vector<std::vector<double>> operator+(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N) {
    // Just a simple wrapper. Doesn't need much memory because of std::move
    return MatrixMatrixSum(std::move(M), std::move(N));
}

std::vector<std::vector<double>> operator*(std::vector<std::vector<double>> M, std::vector<std::vector<double>> N) {
    // Just a simple wrapper. Doesn't need much memory because of std::move
    return MatrixMatrixProduct(std::move(M), std::move(N));
}

double operator*(std::vector<double> A, std::vector<double> B) {
    return VectorVectorProduct(std::move(A), std::move(B));
}

std::vector<double> operator*(std::vector<std::vector<double>> M, std::vector<double> A) {
    // Just a simple wrapper. Doesn't need much memory because of std::move
    return MatrixVectorProduct(std::move(M), std::move(A));
}

std::vector<double> operator+(std::vector<double> A, std::vector<double> B) {
    return VectorSum(std::move(A), std::move(B));
}

std::vector<double> operator-(std::vector<double> A, std::vector<double> B) {
    return VectorDifference(std::move(A), std::move(B));
}

std::vector<double> operator*(double b, std::vector<double> A) {
    return VectorScalarProduct(std::move(A), b);
}

std::vector<double> operator*(std::vector<double> A, double b) {
    return VectorScalarProduct(std::move(A), b);
}

std::vector<std::vector<double>> CholeskyDecompose(std::vector<std::vector<double>> A) {
    std::vector<std::vector<double>> L;
    unsigned long rowsA = A.size();
    L.clear();
    std::vector<double> zeroRow(rowsA, 0);
    L.insert(L.end(), rowsA, zeroRow);

    L[0][0] = sqrt(A[0][0]);

    for (int row = 1; row < rowsA; ++row) {
        for (int column = 0; column <= row; ++column) {

            double sum = 0;
            for (int k = 0; k < column; k++) {
                sum += L[row][k]*L[column][k];
            }

            if (column == row) {
                L[row][column] = sqrt(A[row][column] - sum);
            } else {
                L[row][column] = (1/L[column][column]) * (A[row][column] - sum);
            }
        }
    }

    return L;
}
