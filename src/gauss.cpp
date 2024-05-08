// C++ Libraries
#include "gauss.h"

// Gaussian function
double gaussian(double d, double s) {
    return (1.0 / (sqrt(2 * M_PI) * s)) * exp(-d * d / (2 * s * s));
}

// Create a Gaussian kernel matrix
std::vector<std::vector<double>> createGaussianMatrix(int size, double hardness) {
    //verify arguments
    if(size%2 == 0){
        std::cerr << "ERROR: The dimensions of the Gaussian matrix needs to be odd\n";
        exit(1);
    }

    // Needed variables
    std::vector<std::vector<double>> matrix(size, std::vector<double>(size));
    double sum = 0.0;
    int indx_center = size / 2;

    // Create matrix
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            //distances of index corrdinates to index center
            double x = i - indx_center;
            double y = j - indx_center;
            double d = sqrt(x*x+y*y);
            matrix[i][j] = gaussian(d, hardness);
            sum += matrix[i][j];
        }
    }

    // Normalize the matrix
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix[i][j] /= sum;
        }
    }

    return matrix;
}
