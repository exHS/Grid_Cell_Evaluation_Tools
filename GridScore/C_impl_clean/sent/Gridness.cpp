#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <valarray>     // std::valarray
using namespace std;

// This class calculates the gridness score of a given spike data set.
class Gridness
{
public:
    long dataSize;
    float** spikingData;
    float** rateMap;
    float** autoCorrMap;
    float* gridnessScores;
    string filePath;

private:

    // get the size of the file
    long getFileSize(FILE *file)
    {
        long lCurPos, lEndPos;
        lCurPos = ftell(file);
        fseek(file, 0, 2);
        lEndPos = ftell(file);
        fseek(file, lCurPos, 0);
        return lEndPos;
    }

    // read data into float** spikingData
    void readData(string filePath)
    {
        float *fileBuf;			// Pointer to our buffered data
        FILE *file = NULL;		// File pointer
        this->filePath = filePath;

        // Open the file in binary mode using the "rb" format string
        // This also checks if the file exists and/or can be opened for reading correctly
        if ((file = fopen(filePath.c_str(), "rb")) == NULL)
            cout << "Could not open specified file" << endl;

        // Get the size of the file in bytes
        this->dataSize = getFileSize(file);


        // Allocate space in the buffer for the whole file
        fileBuf = new float[this->dataSize];

        // Read the file in to the buffer
        fread(fileBuf, this->dataSize, 1, file);

        // one float consists of 4 bytes
        this->dataSize /= 4;
        // 3 floats in one row
        this->dataSize /= 3;

        // init fitting array
        this->spikingData = new float*[this->dataSize];

        // bring data in correct format
        for (int i = 0; i < this->dataSize; i++)
        {
            this->spikingData[i] = new float[2];
            this->spikingData[i][0] = fileBuf[i*3 +1];
            this->spikingData[i][1] = fileBuf[i*3 +2];
        }

        // clean
        delete[]fileBuf;
        fclose(file);

    }


    // calculates the maxiumum of the given 2d array
    double maximum(float** arr,int sizeX, int sizeY)
    {

        double maximum = arr[0][0];
        long dataSize = this->dataSize;

        for (int i= 0; i < sizeX; i++)
        {
            for (int j= 0; j < sizeY; j++)
            {

                if (arr[i][j] > maximum)
                {
                    maximum  = arr[i][j];
                }

            }
        }

        return maximum;
    }

    // returns the minium of a given 2d array
    double minimum(float** arr,int sizeX, int sizeY)
    {

        double minimum = arr[0][0];
        long dataSize = this->dataSize;

        for (int i= 0; i < sizeX; i++)
        {
            for (int j= 0; j < sizeY; j++)
            {


                if (arr[i][j] < minimum)
                {
                    minimum  = arr[i][j];
                }

            }
        }

        return minimum;

    }


    // calculates the mean of given data; n is the size of data
    double mean(int n, float *data)
    {
        double m=0;
        int i=n/8;
        while (i>0)
        {
            m += data[0];
            m += data[1];
            m += data[2];
            m += data[3];
            m += data[4];
            m += data[5];
            m += data[6];
            m += data[7];
            data += 8;
            i--;
        }

        switch (n%8)
        {
        case 7:
            m+=data[6];
        case 6:
            m+=data[5];
        case 5:
            m+=data[4];
        case 4:
            m+=data[3];
        case 3:
            m+=data[2];
        case 2:
            m+=data[1];
        case 1:
            m+=data[0];
        }
        return m/n;
    }





    // returns a gaussian kernel with sigma standard deviation
    float** gaussianKernel(float sigma)
    {
        int kernelSize = 2*ceil(2.0*sigma) +1;

        float** kernel = new float*[kernelSize];

        float x, y;

        for (int i = 0; i < kernelSize; i++)
        {
            kernel[i] = new float[kernelSize];
            for (int j = 0; j < kernelSize; j++)
            {
                x = (float)(i - (int)(kernelSize / 2));
                y = (float)(j - (int)(kernelSize / 2));

                kernel[i][j] = exp((-1.0f) * (x * x + y * y) / (2 * sigma * sigma));
                kernel[i][j] = kernel[i][j] / (2 * M_PI * sigma * sigma);
            }
        }
        return kernel;
    }

    // apply a gaussian kernel to the given array
    void gaussianFilter(float** arr,int dimension,float sigma)
    {

        int kernelSize = 2*ceil(2.0*sigma) +1;

        // get kernel
        float** kernel = this->gaussianKernel(sigma);

        // initialize temporal array with values of arr
        float** tmp = new float*[dimension];
        for (int i = 0; i < dimension; i++)
        {
            tmp[i] = new float[dimension];
            for (int j = 0; j < dimension; j++)
            {
                tmp[i][j] = arr[i][j];
            }
        }

        // calculate normalization number of kernel
        float norm = 0.0;
        for (int k = 0 ; k <  kernelSize; k++)
        {
            for (int l = 0; l < kernelSize; l++)
            {
                norm += kernel[k][l];
            }
        }

        // do actual calculation
        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                float sum = 0.0;

                for (int k = -floor(kernelSize/2.0); k <  floor(kernelSize/2.0); k++)
                {
                    for (int l = -floor(kernelSize/2.0); l < floor(kernelSize/2.0); l++)
                    {

                        if( i+k >= 0 && i+k < dimension && j+l >= 0 && j+l < dimension )
                        {
                            //int test = floor(kernelSize/-2.0);
                            sum += kernel[ k + int(floor(kernelSize/2.0)) ][ l + int(floor(kernelSize/2.0)) ] * tmp[i+k][j+l];
                        }

                    }
                }

                arr[i][j] = sum / norm;

            }
        }
    }


    // calculates a rate map of the spikingData and stores it in rateMap variable
    void calcRateMap(int resolution,bool filter,float sigma)
    {
        float maxi = this->maximum(this->spikingData,this->dataSize,2);
        float mini = this->minimum(this->spikingData,this->dataSize,2);

        maxi = abs(ceil(maxi));
        mini = abs(floor(mini));

        if (maxi > mini)
            mini = maxi;

        // prepare the data for further computation
        for(int i = 0; i < this->dataSize; i++)
        {
            for(int j = 0; j < 2; j++)
            {
                this->spikingData[i][j] += mini;
                this->spikingData[i][j] /= (2.0*mini);
                // HACK !!
                this->spikingData[i][j] *= float(resolution-1);
                this->spikingData[i][j] = round(this->spikingData[i][j]);
            }
        }


        this->rateMap = new float*[resolution];
        // build map
        for(int i = 0; i < resolution; i++)
        {
            this->rateMap[i] = new float[resolution];

            for(int j = 0; j < resolution; j++)
            {
                this->rateMap[i][j] = 0;
            }
        }

        int j,k;
        for(int i= 0; i < this->dataSize; i++)
        {

            j = this->spikingData[i][0];
            k = this->spikingData[i][1];

            this->rateMap[k][j] += 1;


        }

        // use gaussian filter on data
        if (filter)
        {
            gaussianFilter(this->rateMap,resolution,sigma);
        }

        maxi = this->maximum(this->rateMap,resolution,resolution);


        // normalize data and mirror it
        for(int i = 0; i < resolution; i++)
        {
            for(int j = 0; j < resolution; j++)
            {
                this->rateMap[i][j] /= maxi;
            }
        }


    }

    // takes a float array and its dimension (quadratic) and builds a list from that data
    float* makeList(float** arr, int dim)
    {
        float* listFloat = new float[dim*dim];
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                listFloat[i*dim+j] = arr[i][j];
            }
        }
        return listFloat;
    }



    // returns the pearson correlation coefficient of the 2 given lists.
    // It is assumed that both lists are of length N
    double pearson(int N, float *data, float *data1)
    {
        double med = mean(N, data), med1 = mean(N, data1), rez = 0, rez1 = 0, rez2 = 0;

        for(int i=0; i<N; i++)
        {
            rez+=(data[i]-med)*(data1[i]-med1);
            rez1+=(data[i]-med)*(data[i]-med);
            rez2+=(data1[i]-med1)*(data1[i]-med1);
        }
        return rez/sqrt(rez1*rez2);
    }

    // calculates the convolution of to matrices
    float** conv2(float **a,float **b, int N){

        float** c = new float*[N+N-1];

        for(int i = 0; i < N+N-1; i++){
            c[i] = new float[N+N-1];
        }


        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++) {
                float w = b[i][j];
                if (w != 0.0){
                    for(int k = 0; k < N; k++){
                        for(int m = 0; m < N; m++){
                            c[i+k][j+m] += w*a[k][m];
                        }
                    }
                }
            }
        }

        return c;
    }

    // calculatess the pearsons correlation coefficient for each entry of the matrix
    float** xcorr2(float **arr,int N)
    {

        // rotate arr
        float** arrMirr = new float*[N];
        float** ones = new float*[N];
        float** arrQuad = new float*[N];

        for(int i = 0; i < N; i++)
        {
            arrMirr[i] = new float[N];
            ones[i] = new float[N];
            arrQuad[i] = new float[N];
            for(int j = 0; j < N; j++)
            {
                arrMirr[i][j] = arr[N-1-i][N-1-j];
                ones[i][j] = 1;
                arrQuad[i][j] = arr[i][j] * arr[i][j];
            }
        }

        float** c = conv2(arr,arrMirr,N);
        float** a = conv2(arrQuad,ones,N);

        float b =0;

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                b += arrQuad[i][j];
            }
        }

        for(int i = 0; i < N+N-1; i++)
        {
            for(int j = 0; j < N+N-1; j++)
            {
                c[i][j] = c[i][j] / sqrt(a[i][j]*b);
            }
        }


        return c;
    }





    // calculate autocorrelation map using pearson's correlation and the rateMap
    void calcAutocorrelationMap(int resolution)
    {

        int nRateMap = resolution*resolution;
        int nAutoCor = resolution+resolution -1;


        // initialize variables
        std::valarray<float> listVal(nRateMap) ;
        for (int i = 0; i < resolution; i++)
        {
            for (int j = 0; j < resolution; j++)
            {
                listVal[i*resolution+j] = this->rateMap[i][j];
            }
        }

        double mean = listVal.sum() / nRateMap;

        for (int i = 0; i < resolution; i++)
        {
            for (int j = 0; j < resolution; j++)
            {
                this->rateMap[i][j] -= mean;
            }
        }

        // initialize variables
        this->autoCorrMap = new float*[nAutoCor];
        float** tmp = new float*[nAutoCor];
        for (int i = 0; i < nAutoCor ; i++)
        {
            this->autoCorrMap[i] = new float[nAutoCor];
            tmp[i] = new float[nAutoCor];
            for (int j = 0; j < nAutoCor ; j++)
            {
                this->autoCorrMap[i][j] = 0.0;
                tmp[i][j] = 0.0;
            }
        }

        this->autoCorrMap = this->xcorr2(this->rateMap,resolution);

        // calc std
        float sd = 0.0;
        for( int i = 0; i < nRateMap; i++ )
        {
          sd += (listVal[i] - mean)*(listVal[i] - mean);
        }
        sd /= nRateMap;
        sd = sqrt( sd );

        // calc correct autoCorrMap
        for (int i = 0; i < nAutoCor ; i++)
        {
            for (int j = 0; j < nAutoCor ; j++)
            {
                this->autoCorrMap[i][j] = this->autoCorrMap[i][j] / (nRateMap*pow(sd,2));
            }
        }


    }


    // rotates the array. assuming the array is quadratic
    float** rotateArray(float** arr,int arrSize, int degree)
    {

        if (degree == 0)
            return arr;

        float** tmp = new float*[arrSize];
        for (int i = 0; i < arrSize; i++)
        {
            tmp[i] = new float[arrSize];
            for (int j = 0; j < arrSize; j++)
            {
                tmp[i][j] = 0;
            }
        }


        float angle = float(degree) * M_PI / 180.0;

        float centerX,centerY;
        float rotatedX,rotatedY;
        int resultX,resultY;

        float center = arrSize/2.0  ;

        float cosAngle = cos(angle);
        float sinAngle = sin(angle);

        for (int i = 0; i < arrSize; i++)
        {
            for (int j = 0; j < arrSize; j++)
            {

                // calc center point
                centerX = i - center;
                centerY = j - center;

                rotatedX = roundf( centerX * cosAngle - centerY * sinAngle );
                rotatedY = roundf( centerX * sinAngle + centerY * cosAngle );

                resultX = rotatedX + center;
                resultY = rotatedY + center;

                if ( resultX >= 0 && resultX < arrSize && resultY >=0 && resultY < arrSize)
                    tmp[i][j] = arr[resultX][resultY];
            }
        }
        return tmp;
    }


    // rotate the map in 6 degree steps and calculate the gridness score for each step
    void calcGridnessScores(int resolution,int inRad,int outRad)
    {

        // times 2 since the autocorrelation map is 2 times bigger
        int dim=resolution*2 -1;

        // maybe this parameter needs to be calculate as well
        int cntr = resolution;


        int** ringFilter = new int*[dim];
        float** ringMap = new float*[dim];
        for (int i = 0; i < dim; i++)
        {
            ringFilter[i] = new int[dim];
            ringMap[i] = new float[dim];
            for (int j =0; j < dim; j++)
            {
                ringFilter[i][j] = 0;
                ringMap[i][j] = 0;
            }
        }

        // build the filter in that the corraltion will be calculated
        float cntrI;
        float cntrJ;
        float dist;
        for (int i = 0; i < dim; i++)
        {
            for (int j =0; j < dim; j++)
            {
                cntrI = pow(cntr-i,2);
                cntrJ = pow(cntr-j,2);
                dist = cntrI+cntrJ;

                if ( pow(inRad,2) <= dist && dist <= pow(outRad,2))
                {
                    ringFilter[i][j] = 1;
                }

            }
        }

        // build ringmap
        for (int i = 0; i < dim; i++)
        {
            for (int j =0; j < dim; j++)
            {
                ringMap[i][j] = this->autoCorrMap[i][j] * ringFilter[i][j];
            }
        }



        // resize ringmap
        int startPos = cntr - outRad ;
        int endPos = cntr + outRad ;
        int dimRingMapResized = endPos-startPos;
        float** ringMapResized = new float*[dimRingMapResized];
        for (int i = 0; i < dimRingMapResized; i++)
        {
            ringMapResized[i] = new float[dimRingMapResized];
            for (int j =0; j < dimRingMapResized; j++)
            {
                ringMapResized[i][j] = 0;

                ringMapResized[i][j] = ringMap[i+startPos][j+startPos];
            }
        }


        this->gridnessScores = new float[30];


        float* ringMapList = this->makeList(ringMapResized,dimRingMapResized);

        // rotate ringMap in one degree steps
        for(int k=0; k< 30; k++)
        {

            float** rotatedArray = this->rotateArray(ringMapResized,dimRingMapResized,k*6);

            float* rotatedArrayList = this->makeList(rotatedArray,dimRingMapResized);

            this->gridnessScores[k] = pearson(dimRingMapResized*dimRingMapResized,ringMapList,rotatedArrayList);

            //cout << gridnessScores[k] << endl;
        }





    }

    // calculate the actual gridness score according to the paper of Sargolini
    float gridScore()
    {
        float gridScore = min(this->gridnessScores[10],this->gridnessScores[20]) - max(this->gridnessScores[5],max(this->gridnessScores[15],this->gridnessScores[25]));
        return gridScore;

    }

    public:






    float getGridScore(string filePath,int resolution,int innerBound, int outerBound,bool filter,float sigma){

        this->readData(filePath);

        this->calcRateMap(resolution,filter,sigma);

//        for (int i = 0; i < resolution; i++)
//        {
//            for (int j =0; j < resolution; j++)
//            {
//                cout << this->rateMap[i][j] << " ";
//            }
//            cout << "\n";
//        }



       this->calcAutocorrelationMap(resolution);


        this->calcGridnessScores(resolution,innerBound,outerBound);

        return this->gridScore();
    }


};
