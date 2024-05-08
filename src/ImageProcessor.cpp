#include "ImageProcessor.h"
#include <algorithm>
#include "gauss.h"
using namespace std;

ImageProcessor::ImageProcessor(){
    //reducing noise
    std::vector<std::vector<double>> S = {
        {0.0, -1.0, 0.0},
        {-1.0, 5.0, -1.0},
        {0.0, -1.0, 0.0}
    };

    std::vector<std::vector<double>> B = {
        {1.0, 1.0, 1.0},
        {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0},
    };

    std::vector<std::vector<double>> M = {
        {1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0},
    };

    std::vector<std::vector<double>> G = createGaussianMatrix(3, 2);


    Matrix2D Sharp(S);
    Matrix2D Box(B);
    Matrix2D Media(M);
    Matrix2D Gaussian(G);

    MatrixStorage["Sharp"] = Sharp;
    MatrixStorage["Box"] = Box;
    MatrixStorage["Media"] = Media;
    MatrixStorage["Gaussian"] = Gaussian;
}

void ImageProcessor::CreateEmpty(std::string filename, int ncols, int nrows, int N){
    if(IsStored(filename)){
        std::cout << "A file with such name already exists" << std::endl;
        return;
    }
    Image image;
    image.nrows = nrows;
    image.ncols = ncols;
    image.N = N;
    image.C.resize(image.ncols, std::vector<int>(image.nrows, N));
    ImageStorage[filename] = image;
}

void ImageProcessor::ReadImage(std::string filename){
    Image image;
    std::ifstream inFile("../img/" + filename);

    //Make sure file is found
    if (!inFile) {
        std::cerr << "Error: Unable to open file for reading." << std::endl;
        return;
    }

    // Read PGM header
    inFile >> image.cod; // Read codification (P2)
    inFile >> image.ncols >> image.nrows; // Read number of columns and rows
    inFile >> image.N; // Read max color value (white)

    // Resize image matrix to match dimensions
    image.C.resize(image.ncols, std::vector<int>(image.nrows, 0));

    // Read pixel values
    for (int row = 0; row < image.nrows; row++) {
        for (int col = 0; col < image.ncols; col++) {
            inFile >> image.C[col][row];
        }
    }
    inFile.close();
    ImageStorage[filename] = image; //store image information
}

void ImageProcessor::InvertColors(std::string stored, std::string filename){
    //check if image to invert exists
    if(!IsStored(stored)){
        std::cout << "Key '" << stored << "' does not exist" << std::endl;
        return;
    }
    //get the image from storage and copy settings
    Image image = CopySettings(stored);

    //change it into a new file
    for (int col = 0; col < image.ncols; col++) {
        for (int row = 0; row < image.nrows; row++) {
            image.C[col][row] = image.N - ImageStorage[stored].C[col][row];
        }
    }
    ImageStorage[filename] = image;
}

void ImageProcessor::ImageSum(std::string img1, std::string img2, std::string result){
    //check if image to invert exists
    if(!IsStored(img1) || !IsStored(img2) || !SameSettings(img1,img2)){
        std::cout << IsStored(img1) << " || " << IsStored(img2) << " || " << SameSettings(img1,img2) << "\n";
        std::cout << "One if the Keys does not exist or images have different settings"<< std::endl;
        return;
    }
    Image image = CopySettings(img1); //doesnt matter if is img1 or img2, they have the same settings

    for (int col = 0; col < image.ncols; col++) {
        for (int row = 0; row < image.nrows; row++) {
            image.C[col][row] = ImageStorage[img1].C[col][row] + ImageStorage[img2].C[col][row];
        }
    }
    ImageStorage[result] = image;
}

void ImageProcessor::ProduceImage(std::string stored) {
    // Check if the image exists
    if (!IsStored(stored)) {return;}

    Image image = ImageStorage[stored];

    std::ofstream outFile("../res/" + stored + ".pgm", std::ios::binary);
    if (!outFile) {
        std::cerr << "Error: Unable to open file '" << stored << ".pgm' for writing." << std::endl;
        return;
    }

    // Write the PGM file header
    outFile << "P2" << std::endl;
    outFile << image.ncols << " " << image.nrows << std::endl;
    outFile << image.N << std::endl;

    // Write pixel data to the file
    for (int row = 0; row < image.nrows; row++) {
        for (int col = 0; col < image.ncols; col++) {
            outFile << image.C[col][row] << " ";
        }
        outFile << std::endl; // Newline after each row
    }

    // Close the file
    outFile.close();
}

void ImageProcessor::addSegment(std::string stored, std::vector<int> start, std::vector<int> end, int color, int
size, string cluster_orientation) {

    if (!IsStored(stored)) {return;}

    if (start.size() != 2 || end.size() != 2) {
        std::cout << "Error: Invalid start or end points." << std::endl;
        return;
    }

    //////////////// X IN THE MATRIX CORRESPONDS TO -Y CARTESIAN /////////////////////////////
    //////////////// Y IN THE MATRIX CORRESPONDES TO X CARTESIAN /////////////////////////////
    ////////////////  CORRESPONDES TO A ROTATION OF -90 DEGREES  /////////////////////////////
    // AN AUXILIARY FUNCTION WAS CREATED TO TAKE A POINT IN CARTESIAN AND DRAW IT IN MATRIX //

    int xi = start[0];
    int yi = start[1];
    int xf = end[0];
    int yf = end[1];
    if(xi>xf){
        xi = end[0];
        yi = end[1];
        xf = start[0];
        yf = start[1];
    }

    double slope=0;
    if(xf-xi!=0){slope= static_cast<double> (yf-yi)/(xf-xi);}

    double x = xi;
    double y = yi;
    if(abs(slope)<=1){
        while (x <= xf) {
            for(int z=x-size/2;z<=x+size/2;z++){
                std::vector<int> c = {z,y};
                std::vector<int> m = cart_to_matrix(ImageStorage[stored].nrows, c);
                ImageStorage[stored].C[m[1]][m[0]] = color;}
            x++;
            y=(x-xi)*slope+yi;
        }
    }
    if(abs(slope)>1){
        while (x <= xf) {
            for(int z=x-size/2;z<=x+size/2;z++){
                std::vector<int> c = {z,y};
                std::vector<int> m = cart_to_matrix(ImageStorage[stored].nrows, c);
                ImageStorage[stored].C[m[1]][m[0]] = color;}
            x+=1/abs(slope);
            y+=slope/abs(slope);
        }
    }
}

void ImageProcessor::SaltAndPepper(std::string stored, int p, std::string filename){
    if (!IsStored(stored)) {return;}
    Image image = CopySettings(stored);

    //iterate image to add SaltAndPepper
    for (int col = 0; col < image.ncols; col++) {
        for (int row = 0; row < image.nrows; row++) {
            if(!GenerateRandomBool(p)){image.C[col][row] = ImageStorage[stored].C[col][row];}
            else{image.C[col][row] = GetRandomExtreme(image.N);}
        }
    }

    //store image
    ImageStorage[filename] = image;
}

void ImageProcessor::GaussianNoise(std::string stored, double stdv, std::string filename){
    if (!IsStored(stored)) {return;}
    Image image = CopySettings(stored);

    //iterate image to add gaussian noise
    int color;
    for (int col = 0; col < image.ncols; col++) {
        for (int row = 0; row < image.nrows; row++) {
            int color = ImageStorage[stored].C[col][row] + (int)Gaussian(stdv);
            //std::cout << "Color: " << color << std::endl;
            if(color<0){image.C[col][row] = 0;}
            if(color>image.N){image.C[col][row] = image.N;}
            else{image.C[col][row] = color;}

        }
    }

    //store image
    ImageStorage[filename] = image;
}

void ImageProcessor::ApplyMatrix(std::string stored, std::string filename, Matrix2D matrix, bool Norm) {
    // Verifications
    if (!IsStored(stored)) return;
    if (matrix.GetMatrix().empty()) return;
    if(Norm){matrix.Normalize();}

    std::vector<std::vector<double>> M = matrix.GetMatrix();
    // Declare variables
    Image image = CopySettings(stored);
    int ni = M.size() / 2;
    int nj = M[0].size() / 2;
    
    // Iterate stored image to define new image pixels
    for (int col = 0; col < image.ncols; col++) {
        for (int row = 0; row < image.nrows; row++) {
            double sum = 0;

            // Iterate over neighborhood
            for (int i = -ni; i <= ni; i++) {
                for (int j = -nj; j <= nj; j++) {
                    int new_col = col + i;
                    int new_row = row + j;
                    
                    // Check bounds
                    if (new_col < 0 || new_col >= image.ncols || new_row < 0 || new_row >= image.nrows){continue;}
                    sum += (double)(ImageStorage[stored].C[new_col][new_row]) * M[i + ni][j + nj];
                }
            }
            image.C[col][row] = (int)(sum);
        }
    }
    
    ImageStorage[filename] = image;
}

void ImageProcessor::ApplyMedian(std::string stored, std::string filename, int dim) {
    // Verifications
    if (!IsStored(stored)) return;
    if(dim%2 == 0) return;

    // Declare variables
    Image image = CopySettings(stored);
    int ni = dim / 2;
    int nj = dim / 2;
    
    // Iterate stored image to define new image pixels
    for (int col = 0; col < image.ncols; col++) {
        for (int row = 0; row < image.nrows; row++) {
            std::vector<int> Values;

            // Iterate over neighborhood
            for (int i = -ni; i <= ni; i++) {
                for (int j = -nj; j <= nj; j++) {
                    int new_col = col + i;
                    int new_row = row + j;
                    // Check bounds
                    if (new_col < 0 || new_col >= image.ncols || new_row < 0 || new_row >= image.nrows){continue;}
                    Values.push_back(ImageStorage[stored].C[new_col][new_row]);
                }
            }
            int Color = Median(Values);
            image.C[col][row] = Color;
        }
    }
    
    ImageStorage[filename] = image;
}

//frequencies and root
std::vector<int> ImageProcessor::GetColourFreq(std::string stored){
    if (!IsStored(stored)) {return {};}
    //create vector to store frequencies
    std::vector<int> AbsFrequencies;
    AbsFrequencies.resize(ImageStorage[stored].N+1, 0);

    //iterate image
    for (int col = 0; col < ImageStorage[stored].ncols; col++){
        for (int row = 0; row < ImageStorage[stored].nrows; row++) {
            AbsFrequencies[ImageStorage[stored].C[col][row]] += 1;
        }
    }
    return AbsFrequencies;
}

std::vector<double> ImageProcessor::GetColourRelFreq(std::string stored){
    if (!IsStored(stored)) {return {};}
    //get absolute frequencies and total number of pixels
    std::vector<int> AbsFrequencies = GetColourFreq(stored);
    int Npixels = ImageStorage[stored].ncols * ImageStorage[stored].nrows;
    //to store relative frequencies
    std::vector<double> RelFrequencies;
    RelFrequencies.resize(ImageStorage[stored].N+1, 0);
    //iterate absolute frequencies to get relative frequencies
    for(int col = 0; col<ImageStorage[stored].N+1; col++){
        RelFrequencies[col] = (double)AbsFrequencies[col] / (double)Npixels;
    }

    return RelFrequencies;

}

double ImageProcessor::MedColor(std::string stored){
    if (!IsStored(stored)) {return {};}
    int Sum = 0;
    int Npixels = ImageStorage[stored].ncols * ImageStorage[stored].nrows;
    for (int col = 0; col < ImageStorage[stored].ncols; col++){
        for (int row = 0; row < ImageStorage[stored].nrows; row++) {
            Sum += ImageStorage[stored].C[col][row];
        }
    }
    return ((double)Sum / (double)Npixels);
}

double ImageProcessor::Variance(std::string stored){
    if (!IsStored(stored)) {return {};}
    double diff;
    double Media = MedColor(stored);
    int Npixels = ImageStorage[stored].ncols * ImageStorage[stored].nrows;
    double Sum = 0;
    for (int col = 0; col < ImageStorage[stored].ncols; col++){
        for (int row = 0; row < ImageStorage[stored].nrows; row++) {
            diff = (double)ImageStorage[stored].C[col][row] - Media;
            Sum += diff*diff;
        }
    }
    return (Sum/(double)(Npixels-1));
}

void ImageProcessor::PlotAbsFreq(std::string stored, std::string filename){
    if (!IsStored(stored)) {return;}
    // Get the absolute frequencies
    std::vector<int> AbsFrequencies = GetColourFreq(stored);
    //plot
    RProcessor.MakeHistogram(AbsFrequencies, filename);
}

void ImageProcessor::PlotRelFreq(std::string stored,std::string filename){
    if (!IsStored(stored)) {return;}
    // Get the relative frequencies
    std::vector<double> AbsFrequencies = GetColourRelFreq(stored);
    //plot
    RProcessor.MakeHistogram(AbsFrequencies, filename);
}

void ImageProcessor::Image2DHist(std::string stored, std::string filename){
    if (!IsStored(stored)) {return;}
    // Get the relative frequencies
    std::vector<std::vector<int>> Image = ImageStorage[stored].C;
    //plot
    RProcessor.MakeHistogram2D(Image, filename);
}

Matrix2D ImageProcessor::GetMatrix(std::string matrix){
    return MatrixStorage[matrix];
}

//auxiliar methods
Image ImageProcessor::CopySettings(std::string stored){
    Image image;
    image.nrows = ImageStorage[stored].nrows;
    image.ncols = ImageStorage[stored].ncols;
    image.N = ImageStorage[stored].N;
    image.C.resize(image.ncols, std::vector<int>(image.nrows, 0));
    return image;
}

bool ImageProcessor::IsStored(std::string stored){
    auto it = ImageStorage.find(stored);
    if (it == ImageStorage.end()) {
        // Key does not exist in the map
        //std::cout << "File '" << stored << "' does not exist." << std::endl;
        return false;
    }
    //std::cout << "File '" << stored << "' exist." << std::endl;
    return true;
}

bool ImageProcessor::SameSettings(std::string img1, std::string img2){
    //check if images exist
    if(!IsStored(img1) || !IsStored(img2)){
        std::cout << "One if the Keys does not exist"<< std::endl;
        return false;
    }
    //MELHORAR COM POINTERS EM VEZ DE FAZER COPIA 
    Image image1 = ImageStorage[img1];
    Image image2 = ImageStorage[img2];
    if(image1.nrows == image2.nrows && image1.ncols == image2.ncols && image1.N == image2.N){
        return true;
    }
    return false;
}

bool ImageProcessor::GenerateRandomBool(int p) {
    if (p <= 0 || p > 100) {
        std::cerr << "Error: Probability must be above 0 and equal or below 100." << std::endl;
        return false;
    }

    // Seed the random number generator
    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    // Generate a random integer between 0 and 99
    std::uniform_int_distribution<> dis(0, 100);
    int randNum = dis(gen);

    // Return true with probability p%
    return randNum <= p;
}

int ImageProcessor::GetRandomExtreme(int N) {
    // Seed the random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Generate a random boolean value
    std::uniform_int_distribution<> dis(0, 1);
    bool value = dis(gen);

    // Return 45 if the random boolean value is false (0), otherwise return 88
    return value ? 0 : N;
}

double ImageProcessor::Gaussian(double stdv) {
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, stdv);

    // Generate a random value from the normal distribution
    return distribution(gen);
}

void ImageProcessor::PrintColorCoordinates(std::string stored, int color){
    if (!IsStored(stored)) {return;}
    for (int col = 0; col < ImageStorage[stored].ncols; col++){
        for (int row = 0; row < ImageStorage[stored].nrows; row++) {
            if(ImageStorage[stored].C[col][row] == color){std::cout << "{ " << col << " , " << row << " } " << std::endl;}
        }
    }
}

std::vector<int> ImageProcessor::cart_to_matrix(int nrows, std::vector<int> c){
    // takes a point in cartesian and returns the corresponding point in the image matrix
    std::vector<int> m;
    m.push_back(nrows - c[1]);
    m.push_back(c[0]);
    return m;
}

int ImageProcessor::Median(std::vector<int> M){
    sort(M.begin(), M.end());
    int Mediana;
    if (M.size() % 2 == 0){
        Mediana = (M[M.size()/2]+M[((M.size()/2)-1)])/2;
      } else{
        Mediana = M[(M.size()-1)/2];
      }
    return Mediana;
}