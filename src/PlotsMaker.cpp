#include "PlotsMaker.h"
#include "PlotsMaker.h"
#include "TApplication.h"
#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>

void PlotsMaker::MakeHistogram(std::vector<int> info, std::string filename){
    std::string destination = "../res/" + filename;
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Histogram", 800, 600);
    // Create a histogram
    TH1F *hist = new TH1F("hist", filename.c_str(), info.size(), 0, info.size());
    // Fill the histogram with the info
    for (int i = 0; i < info.size(); ++i) {
        hist->SetBinContent(i+1, info[i]);
    }
    // Set histogram style
    hist->SetFillColor(38);
    // Draw the histogram
    hist->Draw();
    // Update the canvas to display the histogram
    c1->Update();
    c1->SaveAs(destination.c_str());
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    delete hist;
    delete c1;
}

void PlotsMaker::MakeHistogram(std::vector<double> info, std::string filename){
    std::string destination = "../res/" + filename;
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Histogram", 800, 600);
    // Create a histogram
    TH1F *hist = new TH1F("hist", filename.c_str(), info.size(), 0, info.size());
    // Fill the histogram with the info
    for (int i = 0; i < info.size(); ++i) {
        hist->SetBinContent(i+1, info[i]);
    }
    // Set histogram style
    hist->SetFillColor(38);
    // Draw the histogram
    hist->Draw();
    // Update the canvas to display the histogram
    c1->Update();
    c1->SaveAs(destination.c_str());
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    delete hist;
    delete c1;
}

void PlotsMaker::MakeHistogram2D(std::vector<std::vector<int>> info, std::string filename){
    std::string destination = "../res/" + filename;
    
    // Determine the dimensions of the histogram
    int nx = info.size(); // Number of bins along x-axis
    int ny = info[0].size(); // Number of bins along y-axis

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "2D Histogram", nx, ny);
    
    // Set grayscale color palette
    gStyle->SetPalette(kGreyScale);
    
    // Create a 2D histogram
    TH2F *hist2D = new TH2F("hist2D", filename.c_str(), nx, 0, nx, ny, 0, ny);
    hist2D->SetMinimum(-0.1);
    
    // Fill the histogram with the info
    for (int col = 0; col < nx; col++) {
        for (int row = 0; row < ny; row++) {
            hist2D->Fill(col+0.5, ny-row-0.5, info[col][row]);
        }
    }

    int total = nx*ny;
    int count = 0;
    for (int col = 0; col < nx; col++) {
        for (int row = 0; row < ny; row++) {
            if(info[col][row] == hist2D->GetBinContent(col+0.5, ny-row-0.5)){count += 1;}
        }
    }
    
    std::cout << "TOTAL EQUALITY: " << (double)count/(double)total << std::endl;


    // Draw the 2D histogram
    hist2D->Draw("COLZ");
    
    // Update the canvas to display the histogram
    c1->Update();
    
    // Save the canvas as an image file
    c1->SaveAs(destination.c_str());
    
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    
    // Clean up
    delete hist2D;
    delete c1;
}

void PlotsMaker::MakePointsPlot(std::vector<int>info , std::string filename){
    std::cout << "Empty" << std::endl;
}