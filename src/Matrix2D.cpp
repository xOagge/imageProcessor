#include "Matrix2D.h"

Matrix2D::Matrix2D(std::vector<std::vector<double>> D) {
    // Verify if there aren't vectors with different sizes
    int size = D[0].size();
    for (int col = 0; col < D.size(); col++) { // Corrected loop condition
        if (D[col].size() != size) {
            std::cout << "It's not a matrix, there are vectors with different lengths\n";
            return;
        }
    }

    // Verify if proportions are odd
    if (D.size() % 2 == 0 || D[0].size() % 2 == 0) {
        std::cout << "The matrix sides need to be odd\n";
        return;
    }
    M = D;
}

void Matrix2D::Normalize() {
    // Calculate the total number of elements in the matrix
    double totalElements = 0;
    for (int i = 0; i < M.size(); i++) {
        for (int j = 0; j < M[0].size(); j++) {
            totalElements += M[i][j];
        }
    }

    std::cout << totalElements << "\n";

    // Normalize each element by dividing it by the total number of elements
    for (int i = 0; i < M.size(); i++) {
        for (int j = 0; j < M[0].size(); j++) {
            M[i][j] /= totalElements;
        }
    }
}

std::vector<std::vector<double>> Matrix2D::GetMatrix() {
    return M;
}