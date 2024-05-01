#ifndef __PLOTSMAKER__
#define __PLOTSMAKER__

#include <iostream>
#include <vector>
#include <string>

class PlotsMaker{
    public:
        PlotsMaker() = default;
        ~PlotsMaker() = default;

        //make plots
        void MakeHistogram(std::vector<int> info, std::string filename);
        void MakeHistogram(std::vector<double> info, std::string filename);
        void MakeHistogram2D(std::vector<std::vector<int>> info, std::string filename);

        void MakePointsPlot(std::vector<int> info, std::string filename);

    private:

};

#endif