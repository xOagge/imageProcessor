#include <iostream>
#include "ImageProcessor.h"

using namespace std;

int main(){
    std:cout<< "\n\n\nPROGRAM STARTED:\n";

    ImageProcessor Proc;

    //read glassware, invert and sum original and inverted-------------------------------
    Proc.ReadImage("glassware_noisy.ascii.pgm");
    Proc.InvertColors("glassware_noisy.ascii.pgm", "Inverted_Glassware");
    Proc.ImageSum("glassware_noisy.ascii.pgm", "Inverted_Glassware", "White_Glassware");
    Proc.ProduceImage("Inverted_Glassware");
    Proc.ProduceImage("White_Glassware");

    //segments-------------------------------
    Proc.CreateEmpty("A", 256 , 256, 255);
    std::vector<int> P1 = {68,28};
    std::vector<int> P2 = {128,228};
    std::vector<int> P3 = {188,28};
    std::vector<int> P4 = {98,128};
    std::vector<int> P5 = {158,128};//
    Proc.addSegment("A", P1, P2);
    Proc.addSegment("A", P3, P2);
    Proc.addSegment("A", P1, P2);
    Proc.addSegment("A", P4, P5);
    Proc.ProduceImage("A");

    //Make A image with reduced noise
    Proc.ApplyMedian("glassware_noisy.ascii.pgm", "glassware_Median", 3);
    Proc.ProduceImage("glassware_Median");
    Proc.ApplyMatrix("glassware_Median", "glassware_median_sharp", Proc.GetMatrix("Sharp"), false);
    Proc.ApplyMatrix("glassware_Median", "glassware_median_gaussian", Proc.GetMatrix("Gaussian"), true);
    Proc.ApplyMatrix("glassware_median_sharp", "glassware_median_sharp_gaussian", Proc.GetMatrix("Gaussian"), true);
    Proc.ProduceImage("glassware_median_sharp");
    Proc.ProduceImage("glassware_median_gaussian");
    Proc.ProduceImage("glassware_median_sharp_gaussian");


    //SaltAndPepper-------------------------------
    Proc.SaltAndPepper("glassware_noisy.ascii.pgm", 30, "SaltAndPepper");
    Proc.ProduceImage("SaltAndPepper");

    //gaussian-------------------------------
    Proc.GaussianNoise("glassware_noisy.ascii.pgm", 20 , "GaussianNoise");
    Proc.ProduceImage("GaussianNoise");


    //frequencies and root-------------------------------
    Proc.PlotAbsFreq("glassware_noisy.ascii.pgm", "GlassWare_AbsFreq.pdf");
    Proc.PlotRelFreq("glassware_noisy.ascii.pgm", "GlassWare_RelFreq.pdf");
    Proc.Image2DHist("glassware_noisy.ascii.pgm", "GlassWare_Hist2D.pdf");
    Proc.Image2DHist("A", "A_2D.pdf");

    //Media
    std::cout << "ORIGINAL IMAGE VALUES" << std::endl;
    std::cout <<"Media: " <<  Proc.MedColor("glassware_noisy.ascii.pgm") << std::endl;
    //Variance
    std::cout <<"Variance: " <<  Proc.Variance("glassware_noisy.ascii.pgm") << std::endl;

    //tranformada Hough
    Proc.Mediana_Quadrado("glassware_noisy.ascii.pgm", "Mediana_Quadrado");
    Proc.ProduceImage("Mediana_Quadrado");
    Proc.Threshholding("Mediana_Quadrado", "threshholding", 10);
    Proc.ProduceImage("threshholding");
    //Hough
    Proc.HoughTransform("glassware_noisy.ascii.pgm", "threshholding", "Hough");
    
    return 0;
}
