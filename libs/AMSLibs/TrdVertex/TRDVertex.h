#ifndef TRDVERTEX_H
#define TRDVERTEX_H
#include <vector>
#include "SimpleHitState.h"
#include "root.h"
#include "TMultiGraph.h"
#include "TRDVertex.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TMinuit.h"
#include "Minuit2/ABObj.h"


class TRDVertexFit 
{
  public:
    TRDVertexFit();
    ~TRDVertexFit();
    vector<vector<double> > _cxy;
    vector<vector<double> > _cz;
    vector<vector<double> > _sxy;
    int _ntrack;
    int _ndof;
    double _chi2;
    double _vtxcc[2];
    vector<vector<double> > tinfo;   // a,b,chisq, ndof, normchi2
    void  Clear();
    int   AddTrack( vector<double> &trackx , vector<double> &trackz, vector<double> &errorx  );
    
    int FastIter();
    
    double VertexFit();
    double StepFit(vector<double> &trackx , vector<double> &trackz, vector<double> &errorx , vector<double> &param);
    double FastFit( int nlim = 3 );
    double GetChi2();
    double GetNormChi2();

};

class TRDVertex2DFit: public TObject{
  public:
    TRDVertexFit *fitx;
    TRDVertexFit *fity;
    TRDVertex2DFit(){ _ndof=0; _ntrack=0; _ntothit=0; _chi2 =0;_normchi2=0;};
    ~TRDVertex2DFit(){};
    void SetV( TRDVertexFit *_fit1, TRDVertexFit *_fit2  ) {  fitx = _fit1; fity=_fit2;};
    int _ndof;
    int _ntrack;
    int _ntothit;
    double _chi2;
    double _normchi2;
    AMSPoint _ccvtx;


    double func_Chi2Gain( Double_t *par);
    static void  fcn_TRDChi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

    double DoFit();
    static TMinuit *gMinuit;
    double GetChi2(){ return _chi2; };
    double GetNormChi2(){return _normchi2;};
    void GetVertex( AMSPoint &vcc )  {vcc = _ccvtx;};

};

#endif
