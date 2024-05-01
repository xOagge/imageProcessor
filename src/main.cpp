#include <iostream>
#include "ImageProcessor.h"

int main(){
    ImageProcessor Proc;

    //read glassware, invert and sum original and inverted-------------------------------
    Proc.ReadImage("glassware_noisy.ascii.pgm");
    Proc.InvertColors("glassware_noisy.ascii.pgm", "Inverted_Glassware");
    Proc.ImageSum("glassware_noisy.ascii.pgm", "Inverted_Glassware", "White_Glassware");
    Proc.ProduceImage("Inverted_Glassware");
    Proc.ProduceImage("White_Glassware");

    //segments-------------------------------
    Proc.CreateEmpty("A", 300 , 400, 255);
    std::vector<int> P1 = {10,10};
    std::vector<int> P2 = {200,250};
    Proc.addSegment("A", P1, P2);
    Proc.ProduceImage("A");

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

    //Media
    std::cout << "ORIGINAL IMAGE VALUES" << std::endl;
    std::cout <<"Media: " <<  Proc.MedColor("glassware_noisy.ascii.pgm") << std::endl;
    //Variance
    std::cout <<"Variance: " <<  Proc.Variance("glassware_noisy.ascii.pgm") << std::endl;
    return 0;
}