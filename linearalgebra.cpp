//
// Created by Lars Gebraad on 7/10/17.
// Simple Linear Algebra stuff.
//


#include <vector>
#include <iostream>

// Linear algebra functions


std::vector<double> VectorDifference(std::vector<double> A, std::vector<double> B) {

    std::vector<double> C;

    if (A.size() != B.size()) {
        // Get some Exception class to THROW.
        std::cout
                << "Vectors are not the same dimension! The code DIDN'T run successfully.";
        return C;
    }

    std::vector<double> q_difference;
    for (int i = 0; i < A.size(); i++) {
        // Prior misfit
        C.push_back(A[i] - B[i]);
    }
    return C;
}

std::vector<double> VectorScalarProduct(std::vector<double> A, double b) {
    std::vector<double> B;
    for (int i = 0; i < A.size(); i++) {
        // Prior misfit
        B.push_back(A[i] * b);
    }
    return B;
}

std::vector<double> VectorSum(std::vector<double> A, std::vector<double> B) {

    std::vector<double> C;

    if (A.size() != B.size()) {
        // Get some Exception class to THROW.
        std::cout
                << "Vectors are not the same dimension! The code DIDN'T run successfully.";
        return C;
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
        std::cout
                << "Vector and matrix are not compatible in dimension! The code DIDN'T run successfully.";
        return C;
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
        std::cout
                << "Vectors are not compatible in dimension! The code DIDN'T run successfully.";
        return C;
    }
    C = 0;
    for (int i = 0; i < A.size(); i++) {
        C += (A[i] * B[i]);
    }
    return C;
}
