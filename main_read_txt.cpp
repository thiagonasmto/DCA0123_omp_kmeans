#include "./lib/kMeansCalc.cpp"
//#include "./pascal-releases-master/include/pascalops.h"
//#include "pascalops.h"

/*
Here is an example of the primary way to use the kMeansCalc object;
with a text file.
Here is a main.cpp file to show how it works.
*/

int main(int argc, char **argv){
     if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <dimSize>\n";
        return 1;
    }

    int dimSize = std::atoi(argv[1]); // Receber o tamanho do problema como argumento
    int clusterCount = 20;
    int iterationCount = 1000;
    double threshold = 1;
    
    kMeansCalc<double> k {"./data/mnist_train.txt",dimSize};
    std::cout << "Avg fitness: " << k.doubleFindOptimalClusters(clusterCount,iterationCount,threshold,1) << "\n";

    k.voidWritePointsToFile("out.txt");
    //k.voidPrintSummary();
    //k.voidPrintClusterSummary();

    return 0;
}
