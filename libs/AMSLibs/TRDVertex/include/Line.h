#ifndef LINE_H
#define LINE_H

#include "TMinuit.h"
#include "TObject.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "iostream"
#include "math.h"
#include "TGraph.h"
#include "TLine.h"
#include "TRD2DHit.h"


class Line : public TObject{

public:

    int iside;
    double weight;
    Line();


    void Init(int side);

    Line(int side);

    void Draw();

    //====original measurement=========

    void AddHit(TRD2DHit* hit);


    std::vector<TRD2DHit*> v_hit;

    double x0;
    double z0;

    void SetLineAsPointAndSlopes(double _x0, double _z0, double slope);

    void SetCrossPoint(double _x0, double _z0);


    int Merge(Line* ll);

    static TMinuit *gMinuit;

    double LineFit();
    double LineFit_TMinuit();



    double LineFit_TMinuit_Robust();





    static void fcn_LineFit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

    double fun_LineFit(Double_t *par);

    //===================

    //x=a+bz;
    double a;
    double b;
    double a_err;
    double b_err;
    double chi2;
    int ndof;

    double _aa[2];
    double _nn[2];

    void SetMatrixRepresentation();

    TMatrixD aa;
    TMatrixD nn;
    TMatrixD nn_T;

    static TMatrixD I22;


    double GetDistanceError(double x, double z);
    double GetDistance(double x, double z);

    double GetDistance(TMatrixD p);



    double GetWeightedChi2(int ihit);

    double GetChi2(int ihit);

    double GetTotalWeightedChi2();

    double GetTotalChi2();

    double GetNhit(){return v_hit.size();}

    void Print();


};


bool IsSmallWeightLine(Line *a);


#endif // LINE_H
