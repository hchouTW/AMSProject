#ifndef VERTEX2D_H
#define VERTEX2D_H

#include "TObject.h"
#include "TMinuit.h"
#include "TObject.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "Line.h"
#include "math.h"
#include <vector>
#include <algorithm>
#include "stdlib.h"

class Vertex2D : public TObject{

public:
    Vertex2D();

    std::vector<Line*>v_lines;
    double weight;
    static TMatrixD I22;
    double chi2;
    double nhits;
    int iside;


    bool has3D;
    bool isremoved;

    double dir();

    void Add(Line *line);


    int GetNTracks(){return v_lines.size();}


    void Solve();

    TMatrixD p;

    void SolveOnce();

    void Draw();

    void SetP0andSlopes(double x, double z, double* slope);


    void SetP0(double x, double z);

    double GetTotalChi2();

    double GetTotalNHits();


    bool Contains(Line* line);


    bool ContainsAll(Vertex2D *vv);

    bool BelongsTo(Vertex2D* vv);


    int Merge(Vertex2D *vv);


    int Distance(Vertex2D *vv){
        return sqrt((p-vv->p).E2Norm());
    }

    void Print();



    double x,z;
    double x_err,z_err;

    double V; //probability in arbitrary unit


    double Z(){return z;}
    double XY(){return x;}

    static bool DEBUG;


    static TMinuit *gMinuit;



    static void fcn_CombinedFit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


    double fun_GetCombinedFitChi2(Double_t *par);

    void RemoveOutliers();

    void CombineFit();


    double fun_GetCombinedFitChi2_New(Double_t *par);
    static void fcn_CombinedFit_New(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);



    void CombineFit_New();



    void CombineFit_New_Adaptive();

};



class VertexFinder2D_ZVTOP : public TObject{


public:

    VertexFinder2D_ZVTOP(){
        v_track=new std::vector<Line*>;
        Niteration=0;
        iside=-1;
        x=0;
        z=0;
    }

    VertexFinder2D_ZVTOP( int side){
        v_track=new std::vector<Line*>;
        Niteration=0;
        iside=side;
        x=0;
        z=0;
    }

    ~VertexFinder2D_ZVTOP(){
        delete v_track;
    }

    void SetLineColection(std::vector<Line*> *v){

        *v_track=*v;
        return;
        //        v_track=v;
    }

    void SetVertexColection(std::vector<Vertex2D*> *v){
        v_v2d=v;
    }

    int iside;
    int Niteration;
    static bool DEBUG;

    void FindBestVertex();
    void FindVertex();

    double fun_ZVTOP_V(double* x,double* par);
    double fun_ZVTOP_f(Line* line, double *x);


    static void fcn_ZVTOP(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
    {
        VertexFinder2D_ZVTOP *v=(VertexFinder2D_ZVTOP*)gMinuit->GetObjectFit();
        f=v->fun_ZVTOP_V(par,NULL);
    }

    std::vector<Line*> *v_track;

    double x,z;
    double x_err,z_err;

    static TMinuit *gMinuit;

    std::vector<Vertex2D*> *v_v2d;


};


#endif // VERTEX2D_H
