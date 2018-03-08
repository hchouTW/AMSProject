#ifndef TRDNEWREC_H
#define TRDNEWREC_H
#include "SimpleHitState.h"
#include "root.h"
#include "TMultiGraph.h"
#include "TRDVertex.h"
#include "TRDSTrack.h"
#include "TMinuit.h"
#include "Minuit2/ABObj.h"




class TRDNewRec
{
public:
    TRDNewRec();
    ~TRDNewRec();
    void Init(vector<SimpleHitState> *hits,int itree );
    void Init( AMSEventR *evt,int itree );
    void Clear();
    vector<SimpleHitState> *hitco; 
    int Process();
    AMSPoint _vcc[2];
    double   _vdd[2];
    int      _ntxtrack[2];
//    int AddVertex( AMSPoint cc);
    static double _maxcos;
    static double _maxchi2;
    int GetNTrack( int iside = 1 );
    vector<semiIndex> xcaTrack;
    vector<semiIndex> ycaTrack;
    vector<ccL> xccList;
    vector<ccL> yccList;
    int FillHits( vector<ccL> &ccList, int iside=0) ;
    
    void SaveVertex( TRDVertexFit *tvtx  , int iside , vector<int> &vtkidx);
    int BuildTrk(int iside = 0);
    int RefineTrk( semiIndex &trk, int iside = 1 );
    int BuildVertex( int method =0  );
    
    //new method for iteration way.
    int BuildIterVtx(   );
    int   ReGroup(  vector<double> &czg   );
    int   ReGroup(  vector<double> &cxg ,vector<double> &czg  );

   //old method 
    int BuildVtx(int iside = 0);
    double ReFitVtx();
    vector<int>  chntrkidx[2];

    
    AMSEventR *_evt;
    int _ntrack;
    int _nvertex;
private:
    static int  _trtrkh;
    int _itree;
};
#endif // TRDNEWREC_H
