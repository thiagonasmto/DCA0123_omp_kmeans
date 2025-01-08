/*
kMeansCalculator.
This program takes in data and clusters it.
See all the possible arguments to pass into the kMeansObject.
For instance, can pass in the number of clusters to create, the number 
of iterations, the threshold to stop at, and the number of random seeds.

The program will output summaries related to points / centroids. It can
also output a text file.
*/
#ifndef POINT_READER_CPP
#define POINT_READER_CPP

#include "./PointReader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include <omp.h>

/*

PointReader class.

This is the implementation of the point reader class.
The point reader class is designed to read points for the kMeansCalc and return vectors
of the points and their stats.

*/


/*
**************************
CTOR, DTORS
**************************
*/
// default ctor
template <class T>
PointReader<T>::PointReader(){
    /*
    Inputs:
    void
    Outputs:
    void
    Utility:
    Create reader, set point dim to 0
    */
    this->intPointDimensions = 0;
}
// param ctor that accepts dim size
template <class T>
PointReader<T>::PointReader(const int &dimSize){
    this->intPointDimensions = dimSize;
}

// copy ctor
template <class T>
PointReader<T>::PointReader(const PointReader<T> &reader){
    /*
    Inputs:
    PointReader, PointReader being copied
    Outputs:
    void
    Utility:
    Create point reader in image of reader
    */
    this->myvectorPoints.voidSetEqualLhsRhs(reader.myvectorGetPoints());
    this->intPointDimensions = reader.intGetPointDimensions();
    this->myvectorStdDevs.voidSetEqualLhsRhs(reader.myvectorGetStdDevs());
    this->myvectorAvgs.voidSetEqualLhsRhs(reader.myvectorGetAvgs());
    this->myvectorMins.voidSetEqualLhsRhs(reader.myvectorGetMins());
    this->myvectorMaxes.voidSetEqualLhsRhs(reader.myvectorGetMaxes());
}

// param ctor that accepts 2d arr
template <class T>
PointReader<T>::PointReader(T **data, const int &rows, const int &cols) {
    myVector<Point<T>> points{rows};
    this->intPointDimensions = cols;

    // Vetores para cada propriedade
    myVector<double> stddevs{this->intPointDimensions, 0};
    myVector<double> avgs{this->intPointDimensions, 0};
    myVector<double> mins{this->intPointDimensions, __DBL_MAX__};
    myVector<double> maxes{this->intPointDimensions, __DBL_MIN__};

    // Paralelizar o loop principal - paralelizado, segunda versão
    #pragma omp parallel
    {
        // Variáveis locais para evitar conflitos entre threads
        myVector<double> local_avgs(cols, 0);
        myVector<double> local_mins(cols, __DBL_MAX__);
        myVector<double> local_maxes(cols, __DBL_MIN__);

        // - paralelizado, segunda versão
        #pragma omp for
        for (int i = 0; i < rows; ++i) {
            Point<T> newPoint{cols};
            for (int j = 0; j < cols; ++j) {
                double val = data[i][j];
                newPoint[j] = val;

                // Atualiza valores locais
                if (val >= local_maxes[j]) {
                    local_maxes[j] = val;
                }
                if (val <= local_mins[j]) {
                    local_mins[j] = val;
                }
                local_avgs[j] += val;
            }
            points[i] = newPoint;
        }

        // Redução para mins, maxes e avgs - paralelizado, segunda versão
        #pragma omp critical
        {
            for (int j = 0; j < cols; ++j) {
                if (local_mins[j] < mins[j]) mins[j] = local_mins[j];
                if (local_maxes[j] > maxes[j]) maxes[j] = local_maxes[j];
                avgs[j] += local_avgs[j];
            }
        }
    }

    // Divide somas para calcular médias
    int pointCount = points.intLength();
    for (int i = 0; i < this->intPointDimensions; ++i) {
        avgs[i] /= pointCount;
    }

    // Calcula desvios padrão (stddev) - paralelizado, segunda versão
    #pragma omp parallel for
    for (int i = 0; i < pointCount; ++i) {
        for (int j = 0; j < this->intPointDimensions; ++j) {
            double diff = this->myvectorPoints.tGetByReference(i).tGetValAtDimNByReference(j) - avgs[j];
            diff = (diff * diff);

            #pragma omp atomic
            stddevs[j] += diff;
        }
    }

    // Divide pelo número de pontos e calcula raiz quadrada
    for (int i = 0; i < this->intPointDimensions; ++i) {
        stddevs[i] = std::sqrt((stddevs[i] / pointCount));
    }

    // Agora atribuir os vetores aos atributos do objeto
    this->myvectorPoints.voidSetEqualLhsRhs(points);
    this->myvectorStdDevs.voidSetEqualLhsRhs(stddevs);
    this->myvectorAvgs.voidSetEqualLhsRhs(avgs);
    this->myvectorMins.voidSetEqualLhsRhs(mins);
    this->myvectorMaxes.voidSetEqualLhsRhs(maxes);
}

// param ctor
template <class T>
PointReader<T>::PointReader(const std::string &fileName,const int &intPointDimensions){
    /*
    Inputs:
    string, int, File name, number of dimensions
    Outputs:
    void
    Utility:
    Reads points based on a passed in file name and sets the dimensions to the passed in dimension amount
    */
    this->intPointDimensions = intPointDimensions;
    voidReadTxt(fileName);
}

// dtor
template <class T>
PointReader<T>::~PointReader(){
    /*
    Inputs:
    void
    Outputs:
    void    
    Utility:
    PointReader destructor
    */

}

/*
**************************
GETTERS, SETTERS
**************************
*/
// get vector of stddevs
template <class T>
myVector<double> PointReader<T>::myvectorGetStdDevs() const{
    /*
    Inputs: 
    void
    Outputs: 
    myVector, vector containing the standard deviations for each dimension of the points
    Utility: 
    Getter for get the standard deviation vector for all the points in the object
    */
    return this->myvectorStdDevs;
}
// get the vector of avgs
template <class T>
myVector<double> PointReader<T>::myvectorGetAvgs() const{
    /*
    Inputs:
    void
    Outputs: 
    myVector, vector containing the averages for each dimension of the points
    Utility: 
    Getter for get the averages vector for all the points in the object
    */
    return this->myvectorAvgs;
}
// get the vector of mins
template <class T>
myVector<double> PointReader<T>::myvectorGetMins() const{
    /*
    Inputs:
    void
    Outputs: 
    myVector, vector containing the minimums for each dimension of the points
    Utility: 
    Getter for get the minimums vector for all the points in the object
    */
    return this->myvectorMins;
}
// get the vector of maxes
template <class T>
myVector<double> PointReader<T>::myvectorGetMaxes() const{
    /*
    Inputs:
    void
    Outputs: 
    myVector, vector containing the maxes for each dimension of the points
    Utility: 
    Getter for get the maxes vector for all the points in the object
    */
    return this->myvectorMaxes;
}
// get the vector of points
template <class T>
myVector<Point<T>> PointReader<T>::myvectorGetPoints() const{
    /*
    Inputs:
    void
    Outputs:
    myVector, vector containing all of the points in the object
    Utility:
    returns a vector containing all of the points in the object
    */
    return this->myvectorPoints;
}
// get point count
template <class T>
int PointReader<T>::intGetPointCount() const{
    return this->myvectorPoints.intLength();
}
// get point dimension size 
template <class T>
int PointReader<T>::intGetPointDimensions() const{
    /*
    Inputs:
    void
    Outputs:
    int, number of dimensions 
    Utility:
    Returns the number of dimensions each point has in the object
    */
    return this->intPointDimensions;
}
// get the total lines of file
template <class T>
int PointReader<T>::intGetTotalLines(const std::string &fileName) {
    /*
    Inputs:
    string, file name
    Outputs:
    int, number of lines in the file
    Utility:
    Gets the total number of lines that are within the passed in file
    */
    // Abre o arquivo em modo binário para leitura rápida
    std::ifstream ifstreamFile(fileName, std::ios::in | std::ios::binary);
    if (!ifstreamFile) {
        std::cerr << "Error: file cannot be opened.\n";
        exit(1);
    }

    // Move para o final para determinar o tamanho total do arquivo
    ifstreamFile.seekg(0, std::ios::end);
    std::streampos fileSize = ifstreamFile.tellg();
    ifstreamFile.seekg(0, std::ios::beg);

    // Lê o arquivo em blocos paralelos
    int totalLines = 0;
    const int numThreads = omp_get_max_threads();  // Número de threads disponíveis
    const std::streamsize blockSize = fileSize / numThreads;

    #pragma omp parallel reduction(+:totalLines) // Paralelizado na segunda versão
    {
        const int threadID = omp_get_thread_num();
        std::streampos start = threadID * blockSize;
        std::streampos end = (threadID == numThreads - 1) ? fileSize : static_cast<std::streampos>((threadID + 1) * blockSize);

        // Buffer para leitura local
        std::vector<char> buffer(end - start);
        ifstreamFile.seekg(start, std::ios::beg);
        ifstreamFile.read(buffer.data(), buffer.size());

        // Conta quebras de linha dentro do buffer
        int localCount = 0;
        for (char c : buffer) {
            if (c == '\n') localCount++;
        }

        totalLines += localCount;
    }

    // Adiciona a linha final, se não terminar com '\n'
    ifstreamFile.clear();  // Limpa os flags de EOF
    ifstreamFile.seekg(-1, std::ios::end);
    char lastChar;
    ifstreamFile.get(lastChar);
    if (lastChar != '\n') totalLines++;

    return totalLines;
}

/*
**************************
void returning functions
**************************
*/
// read txt file
template <class T>
void PointReader<T>::voidReadTxt(const std::string &fileName) {
    /*
    Inputs:
    string, file name
    Outputs:
    void
    Utility:
    Reads the data from the passed-in file and populates the object with points being read from the file
    */

    // Calcula o total de linhas (número de pontos) e inicializa o vetor
    int totalLines = intGetTotalLines(fileName);
    this->myvectorPoints.voidSetSize(totalLines);

    // Vetores para armazenar resultados parciais
    myVector<double> stddevs(this->intPointDimensions, 0);
    myVector<double> avgs(this->intPointDimensions, 0);
    myVector<double> mins(this->intPointDimensions, __DBL_MAX__);
    myVector<double> maxes(this->intPointDimensions, __DBL_MIN__);

    // Abrimos o arquivo para leitura
    std::ifstream ifstreamFile(fileName, std::ios::in);
    if (!ifstreamFile) {
        std::cerr << "Error: file cannot be opened.\n";
        exit(1);
    }

    // Lê todas as linhas do arquivo e as armazena em memória
    std::vector<std::string> fileLines(totalLines);
    for (int i = 0; i < totalLines && std::getline(ifstreamFile, fileLines[i]); ++i) {}

    ifstreamFile.close();

    // Vetores locais para cada thread
    std::vector<myVector<double>> localAvgs(omp_get_max_threads(), myVector<double>(this->intPointDimensions, 0));
    std::vector<myVector<double>> localMins(omp_get_max_threads(), myVector<double>(this->intPointDimensions, __DBL_MAX__));
    std::vector<myVector<double>> localMaxes(omp_get_max_threads(), myVector<double>(this->intPointDimensions, __DBL_MIN__));

    // Processa os pontos em paralelo
    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        myVector<T> myvectorTCoords(this->intPointDimensions);
        Point<T> pointNew(myvectorTCoords);

        #pragma omp for
        for (int i = 0; i < totalLines; ++i) {
            std::stringstream ss(fileLines[i]);
            T tReadValueAtDimN;
            int intDimCounter = 0;

            while (ss >> tReadValueAtDimN && intDimCounter < this->intPointDimensions) {
                myvectorTCoords[intDimCounter] = tReadValueAtDimN;

                // Atualiza mínimos e máximos locais
                if (tReadValueAtDimN > localMaxes[threadID].tGetByReference(intDimCounter)) {
                    localMaxes[threadID][intDimCounter] = tReadValueAtDimN;
                }
                if (tReadValueAtDimN < localMins[threadID].tGetByReference(intDimCounter)) {
                    localMins[threadID][intDimCounter] = tReadValueAtDimN;
                }

                // Soma para calcular médias locais
                localAvgs[threadID][intDimCounter] += tReadValueAtDimN;
                intDimCounter++;
            }

            pointNew.voidSetCoordVector(myvectorTCoords);
            this->myvectorPoints[i] = pointNew;
        }
    }

    for (int i = 0; i < this->intPointDimensions; ++i) {
        for (int t = 0; t < omp_get_max_threads(); ++t) {
            avgs[i] += localAvgs[t][i];
            mins[i] = std::min(mins[i], localMins[t][i]);
            maxes[i] = std::max(maxes[i], localMaxes[t][i]);
        }
        avgs[i] /= totalLines;
    }

    // Calcula o desvio padrão em paralelo
    #pragma omp parallel for
    for (int i = 0; i < totalLines; ++i) {
        for (int j = 0; j < this->intPointDimensions; ++j) {
            double diff = this->myvectorPoints.tGetByReference(i).tGetValAtDimNByReference(j) - avgs[j];
            #pragma omp atomic
            stddevs[j] += diff * diff;
        }
    }

    for (int i = 0; i < this->intPointDimensions; ++i) {
        stddevs[i] = std::sqrt(stddevs[i] / totalLines);
    }

    // Atribui os valores calculados aos membros da classe
    this->myvectorStdDevs.voidSetEqualLhsRhs(stddevs);
    this->myvectorAvgs.voidSetEqualLhsRhs(avgs);
    this->myvectorMins.voidSetEqualLhsRhs(mins);
    this->myvectorMaxes.voidSetEqualLhsRhs(maxes);

}

// read csv file
template <class T>
void PointReader<T>::voidReadCSV(const std::string &strFileName){
    /*
    Inputs:
    string, file name
    Outputs: 
    void
    Utility: 
    reads a csv file to populate the point reader with points
    */
    std::ifstream fIn;    
    std::string line;
    this->myvectorPoints.voidSetSize(intGetTotalLines(strFileName));

    fIn.open(strFileName);
    if(!fIn){
        std::cout << "Error: file cannot be opened.\n";
        exit(1);
    }

    std::string strLine;
    int counter;
    int intLineCounter = 0;

    myVector<double> stddevs {this->intPointDimensions,0};
    myVector<double> avgs {this->intPointDimensions,0};
    myVector<double> mins {this->intPointDimensions,__DBL_MAX__};
    myVector<double> maxes {this->intPointDimensions,__DBL_MIN__};

    while(fIn.good()){
        while(getline(fIn,strLine)){
            myVector<T> pointTCoords {this->intPointDimensions};
            counter = 0;
            std::stringstream ss(strLine);
            std::string token;
            while(getline(ss,token,',')){
                double val = std::atof(token.c_str());
                // put the val in min, max, std, etc
                if(val >= maxes.tGetByReference(counter)){
                    maxes[counter] = val;
                }
                if(val <= mins.tGetByReference(counter)){
                    mins[counter] = val;
                }
                // add to avgs, later we divide by num points
                avgs[counter] += val;
                pointTCoords[counter] = val;
                counter ++;
            }
            Point<T> newPoint {pointTCoords};
            this->myvectorPoints[intLineCounter] = newPoint;
            intLineCounter ++ ;
        }
    }

    fIn.close();
    // divide total points by num points to get avgs
    for(int i=0;i<this->intPointDimensions;++i){
        avgs[i] /= this->myvectorPoints.intLength();
    }
    // s^2 = sum(xi - mean(x))^2 * (1/(n)) pop stddev
    for(int i=0;i<this->myvectorPoints.intLength();++i){
        for(int j=0;j<this->intPointDimensions;++j){
            double diff = (this->myvectorPoints.tGetByReference(i).tGetValAtDimNByReference(j) - avgs[j]);
            diff = (diff * diff);

            stddevs[j] += diff;
        }
    }
    int pointCount = this->myvectorPoints.intLength();
    // now divide by num points (pop std) and square root
    for(int i=0;i<this->intPointDimensions;++i){
        stddevs[i] = std::sqrt((stddevs.tGetByReference(i)/(pointCount)));
    }

    // now assign vectors
    this->myvectorStdDevs = stddevs;
    this->myvectorAvgs = avgs;
    this->myvectorMins = mins;
    this->myvectorMaxes = maxes;
}
// set the point dim size
template <class T>
void PointReader<T>::voidSetIntPointDimensions(const int &dimSize){
    /*
    Inputs:
    int, number of dimension to be set to
    outputs: 
    void
    Utility: 
    Sets the number of dimensions of the point reader
    */
    this->intPointDimensions = dimSize;
}
// set vector of points
template <class T>
void PointReader<T>::voidSetMyVectorPoints(const myVector<Point<T>>&points){
    /*
    Inputs:
    myVector, points vector to be set to 
    outputs: 
    void
    Utility: 
    Sets the number of points vector to the passed in points vector
    */
    this->myvectorPoints.voidSetEqualLhsRhs(points);
}
// set vector of avgs
template <class T>
void PointReader<T>::voidSetMyVectorAvgs(const myVector<double>&points){
    /*
    Inputs:
    myVector, points vector to have the averages set to 
    outputs: 
    void
    Utility: 
    Sets the averages of the current point reader to the averages of the passed in vector
    */
    this->myvectorAvgs.voidSetEqualLhsRhs(points);
}
// set min vector
template <class T>
void PointReader<T>::voidSetMyVectorMins(const myVector<double>&points){
    /*
    Inputs:
    myVector, points vector to have the minimums set to 
    outputs: 
    void
    Utility: 
    Sets the minimums of the current point reader to the minimums of the passed in vector
    */
    this->myvectorMins.voidSetEqualLhsRhs(points);
}
// set vector of maxes
template <class T>
void PointReader<T>::voidSetMyVectorMaxes(const myVector<double>&points){
    /*
    Inputs:
    myVector, points vector to have the maximums set to 
    outputs: 
    void
    Utility: 
    Sets the maxes of the current point reader to the maxes of the passed in vector
    */
    this->myvectorMaxes.voidSetEqualLhsRhs(points);
}
// set vector of stddevs
template <class T>
void PointReader<T>::voidSetMyVectorStdDevs(const myVector<double>&points){
    /*
    Inputs:
    myVector, points vector to have the standard deviations set to 
    outputs: 
    void
    Utility: 
    Sets the standard deviations of the current point reader to the standard deviations of the passed in vector
    */
    this->myvectorStdDevs.voidSetEqualLhsRhs(points);
}
// check that a point has dimensions
template <class T>
void PointReader<T>::voidCheckHasPointDimensions() const{
    /*
    Inputs: 
    void
    Outputs: 
    void
    Utility: 
    Checks to see if point dimensions have been initialized, exits if they haven't
    */
    if(this->intPointDimensions <= 0){
        std::cout << "Error: PointReader.cpp\n";
        std::cout << "Point dimensions uninitialized. Check arguments on ctor.\n";
        exit(1);
    }
}


/*
**************************
operator
**************************
*/
// stream extraction
template <class T>
std::istream &operator>>(std::istream &fileInput, PointReader<T> &aReader) {
    /*
    Inputs: 
    istream, PointReader, istream to have data pulled from, point reader to have data put into 
    Outputs: 
    istream
    Utility: 
    Stream extraction operator for populating a PointReader object
    */

    // Check if the object has point dimensions.
    aReader.voidCheckHasPointDimensions();

    // Read the file to determine the total number of lines
    std::vector<std::string> lines;
    std::string strLine;
    while (std::getline(fileInput, strLine)) {
        lines.push_back(strLine);
    }

    // Reset the input stream
    fileInput.clear();
    fileInput.seekg(0, std::ios::beg);

    int intLineCounter = lines.size();
    int intPointDimensions = aReader.intGetPointDimensions();

    // Initialize vectors for statistics
    myVector<double> stddevs(intPointDimensions, 0);
    myVector<double> avgs(intPointDimensions, 0);
    myVector<double> mins(intPointDimensions, __DBL_MAX__);
    myVector<double> maxes(intPointDimensions, __DBL_MIN__);

    myVector<Point<T>> points(intLineCounter);

    // Parallel processing of lines
    #pragma omp parallel
    {
        myVector<double> localAvgs(intPointDimensions, 0);
        myVector<double> localMins(intPointDimensions, __DBL_MAX__);
        myVector<double> localMaxes(intPointDimensions, __DBL_MIN__);
        myVector<double> localStddevs(intPointDimensions, 0);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < intLineCounter; ++i) {
            std::stringstream ss(lines[i]);
            myVector<T> myvectorTCoords(intPointDimensions);
            Point<T> pointNew(myvectorTCoords);

            T tValAtDimN = 0;
            int counter = 0;

            while (ss >> tValAtDimN && counter < intPointDimensions) {
                myvectorTCoords[counter] = tValAtDimN;

                // Update local statistics
                localAvgs[counter] += tValAtDimN;
                localMins[counter] = std::min(localMins[counter], tValAtDimN);
                localMaxes[counter] = std::max(localMaxes[counter], tValAtDimN);

                counter++;
            }

            pointNew.voidSetCoordVector(myvectorTCoords);
            points[i] = pointNew;
        }

        // Reduce results to global vectors
        #pragma omp critical
        {
            for (int j = 0; j < intPointDimensions; ++j) {
                avgs[j] += localAvgs[j];
                mins[j] = std::min(mins[j], localMins[j]);
                maxes[j] = std::max(maxes[j], localMaxes[j]);
            }
        }
    }

    // Calculate averages
    int pointCount = points.intLength();
    for (int i = 0; i < intPointDimensions; ++i) {
        avgs[i] /= pointCount;
    }

    // Compute standard deviations
    // #pragma omp parallel for schedule(dynamic) reduction(+ : stddevs[:intPointDimensions])
    for (int i = 0; i < pointCount; ++i) {
        for (int j = 0; j < intPointDimensions; ++j) {
            double diff = points[i].tGetValAtDimNByReference(j) - avgs[j];
            stddevs[j] += diff * diff;
        }
    }

    for (int i = 0; i < intPointDimensions; ++i) {
        stddevs[i] = std::sqrt(stddevs[i] / pointCount);
    }

    // Assign results to the reader object
    aReader.voidSetMyVectorAvgs(avgs);
    aReader.voidSetMyVectorMaxes(maxes);
    aReader.voidSetMyVectorMins(mins);
    aReader.voidSetMyVectorStdDevs(stddevs);
    aReader.voidSetMyVectorPoints(points);

    return fileInput;
}

// equality operator
template <class T>
bool PointReader<T>::operator == (const PointReader<T>&rhs) const{
    /*
    Inputs: 
    PointReader, point reader being tested for equality
    Outputs: 
    bool, true if point readers are equal, false if not
    Utility: 
    equality operator for testing if two point readers are equal
    */
    if(rhs.myvectorGetAvgs() != this->myvectorAvgs){
        return false;
    }
    if(rhs.myvectorGetMaxes() != this->myvectorMaxes){
        return false;
    }
    if(rhs.myvectorGetMins() != this->myvectorMins){
        return false;
    }
    if(rhs.myvectorGetAvgs() != this->myvectorAvgs){
        return false;
    }
    if(this->intPointDimensions != rhs.intGetPointDimensions()){
        return false;
    }
    if(this->myvectorPoints != rhs.myvectorGetPoints()){
        return false;
    }
    return true;
}
// inequality operator
template <class T>
bool PointReader<T>::operator != (const PointReader<T>&rhs)const{
    /*
    Inputs: 
    PointReader, point reader being tested for inequality
    Outputs: 
    bool, true if point readers are not equal, false if they are
    Utility: 
    inequality operator for testing if two point readers are not equal
    */
    return !(*this == rhs);
}
// assignment operator
template <class T>
const PointReader<T>& PointReader<T>::operator = (const PointReader<T> &rhs){
    /*
    Inputs: 
    PointReader, point reader to be set equal to
    Outputs: 
    PointReader, passed in point reader for cascadability
    Utility: 
    Assignment operator to set the current point reader to be equal to the passed in point reader
    */

    // set all attributes equal
    this->intPointDimensions = rhs.intGetPointDimensions();

    this->myvectorPoints.voidSetEqualLhsRhs(rhs.myvectorGetPoints());
    this->myvectorStdDevs.voidSetEqualLhsRhs(rhs.myvectorGetStdDevs());
    this->myvectorAvgs.voidSetEqualLhsRhs(rhs.myvectorGetAvgs());
    this->myvectorMins.voidSetEqualLhsRhs(rhs.myvectorGetMins());
    this->myvectorMaxes.voidSetEqualLhsRhs(rhs.myvectorGetMaxes());

    return *this;
}

#endif