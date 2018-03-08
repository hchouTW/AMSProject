#ifndef SEGMENT_H
#define SEGMENT_H

#include <vector>
#include <TChain.h>
#include <iostream>
#include <TMath.h>
#include <TString.h>
#include "Cell.h"

class Segment{


public:
    Segment();
    Segment(Cell* start, double theta_threshold=0.9);

    int Add(Cell* _cell);


    int Finish();


    int GetNHits(){return hits.size();}
    int GetNCells(){return cells.size();}

    std::vector<Cell*> cells;
    std::vector<TRD2DHit*> hits;
    void Print();

    double GetLength(){
        return sqrt(pow(z2-z1,2)+pow(x2-x1,2));
    }

    void FitStraightLine(double &chi2, double &ndof);

    void Draw();

    float slope;
    float intercept;
    float x1,y1,z1;
    float x2,y2,z2;
    void Init();
};


#endif // SEGMENT_H
