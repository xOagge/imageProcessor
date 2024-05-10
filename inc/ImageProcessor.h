#ifndef __IMAGEPROCESSOR__
#define __IMAGEPROCESSOR__

//STL´s
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <random>
#include "Matrix2D.h"

#include "PlotsMaker.h"

struct Image;

class ImageProcessor{
    public:
        ImageProcessor();
        ~ImageProcessor() = default;

        //manipulate pgm files
        void CreateEmpty(std::string filename, int ncols, int nrows, int N);
        void ReadImage(std::string filename);
        void InvertColors(std::string stored, std::string filename);
        void ImageSum(std::string img1, std::string img2, std::string result);
        void ProduceImage(std::string stored);
        void addSegment(std::string stored, std::vector<int> start, std::vector<int> end, int color=0, int
size=2, std::string cluster_orientation="H"); //not working properly
        void SaltAndPepper(std::string stored, int p, std::string filename);
        void GaussianNoise(std::string stored, double stdv, std::string filename);
        void ApplyMatrix(std::string stored, std::string filename, Matrix2D matrix, bool Norm);
        void ApplyMedian(std::string stored, std::string filename, int dim);
        void Threshholding(std::string stored, std::string filename, double threshhold);
        void HoughTransform(std::string stored, std::string edged, std::string filename);
        void Mediana_Quadrado(std::string stored, std::string filename);

        //frequencies and root
        std::vector<int> GetColourFreq(std::string stored);
        std::vector<double> GetColourRelFreq(std::string stored);
        double MedColor(std::string stored);
        double Variance(std::string stored);
        void PlotAbsFreq(std::string stored, std::string filename);
        void PlotRelFreq(std::string stored, std::string filename);
        void Image2DHist(std::string stored, std::string filename);
        Matrix2D GetMatrix(std::string matrix);
        

    private:
        //auxiliar methods
        Image CopySettings(std::string stored); //makes new image with same settings as the stored
        bool IsStored(std::string stored);
        bool SameSettings(std::string img1, std::string img2);
        bool GenerateRandomBool(int p);
        int GetRandomExtreme(int N); //N to know what is the white color
        double Gaussian(double stdv);
        void PrintColorCoordinates(std::string stored, int color);
        std::vector<int> cart_to_matrix(int nrows, std::vector<int> c);
        int Median(std::vector<int> M);

        std::map<std::string, Image> ImageStorage; //filename->Image
        std::map<std::string, Matrix2D> MatrixStorage;
        PlotsMaker RProcessor;

};

struct Image {
    std::string cod; // codification (not using)
    int nrows, ncols; //numb of rows and columns
    int N; // max color value (white)
    std::vector<std::vector<int>> C; // matrix of pixel colors
};

#endif
