#ifndef TRD2DHIT_H
#define TRD2DHIT_H



#include "SimpleHitState.h"

#include "TFile.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "math.h"
#include "TRandom3.h"
#include "TLine.h"
#include <TColor.h>
#include <TLatex.h>
#include <TText.h>
#include <vector>
#include <algorithm>
#include "root.h"
#include <iostream>

#include "Cell.h"
class Cell;

class TRD2DHit{

public:
    TRD2DHit();
    TRD2DHit(float _x, float _z,float _Amp,int _layer);

    TRD2DHit(SimpleHitState simple);
    TRD2DHit(TrdRawHitR* raw);

    ~TRD2DHit(){}

    bool operator()(TRD2DHit i,TRD2DHit j);
    bool operator()(TRD2DHit *i,TRD2DHit *j);


    float GetDistance(TRD2DHit* p2);

    float GetDistance_X(TRD2DHit* p2);

    float GetDistance_Z(TRD2DHit* p2);

    void AddConnectedUpper(Cell* cell);
    void AddConnectedLower(Cell* cell);

    bool MergeWith(TRD2DHit *hit);

    void ClearConnections();

    int type;  // 0 for XZ,  1 for YZ Plane

    double y;
    double z;
    double y_err;
    double z_err;
    double Amp;
    int layer;
    double weight;
    bool ismergedhit;
    float boundary[2];
    bool isused;

    std::vector<Cell*> v_connected_upper;
    std::vector<Cell*> v_connected_lower;


//    std::vector<TRD2DHit*> v_connectedhit_upper;
//    std::vector<TRD2DHit*> v_connectedhit_lower;


    double connect_upper_mindist;
    double connect_lower_mindist;

    virtual void Print();

    void Init();
    void AddConnectedHitUpper(TRD2DHit *hit);
    void AddConnectedHitLower(TRD2DHit *hit);
};




class TRD2DHitCollection{

public:
    TRD2DHitCollection(){
        Init();
    }
    ~TRD2DHitCollection(){
        ClearMemory();
    }

    void Add(SimpleHitState simplehit);



    void Add(TRD2DHit *hit);



    void Add_WidthMerging(SimpleHitState simplehit);

    void Add_WidthMerging(TRD2DHit *hit);

    int Size(int i);
    void Init();
    void ClearMemory();

    void Draw_XZ();
    void Draw_YZ();
    void Draw();

    int runid;
    int eventid;
    static const int NLayers=20;

    TRD2DHit* GetHit(int ilayer, int index);


    std::vector<TRD2DHit*> Hits[NLayers];


    TRD2DHit *MakeHitState(SimpleHitState simple);

    int NTotalHits;


    void Sort();
};




#endif // TRD2DHIT_H
