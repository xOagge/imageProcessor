#ifndef __MATRIX2D__
#define __MATRIX2D__

#include <iostream>
#include <vector>
#include <string>
#include <map>


class Matrix2D{
    public:
        Matrix2D() = default;
        Matrix2D(std::vector<std::vector<double>> D);
        ~Matrix2D() = default;

        void Normalize();
        std::vector<std::vector<double>> GetMatrix();

    private:
        std::vector<std::vector<double>> M;

};

#endif