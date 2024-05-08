#ifndef __GAUSS__
#define __GAUSS__

// C++ lIBRARIES
#include <iostream>
#include <vector>
#include <cmath>

// Gaussian function
double gaussian(double x, double s);

// Create a Gaussian kernel matrix
std::vector<std::vector<double>> createGaussianMatrix(int size, double hardness);

#endif
