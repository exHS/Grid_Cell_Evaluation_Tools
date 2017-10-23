#include <stdio.h>
#include <iostream>
#include "src/Gridness.cpp"
using namespace std;

int main()
{
    char *filePath = "resizedData/10704-07070407_T2C3.dat";


    int resolution = 60;
    bool filter = true;
    float sigma = 1.5;
    int innerBound = 18;
    int outerBound = 35;

    Gridness grid;


    cout << grid.getGridScore(filePath,resolution,innerBound,outerBound,filter,sigma) << endl;


    return 0;
}
