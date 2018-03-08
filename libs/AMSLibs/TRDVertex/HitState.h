#ifndef _HITSTATE_H
#define _HITSTATE_H

//#include "TrRecHit.h"
//#include "ParticleState.h"
//#include "TrdKCluster.h"
#include "SimpleHitState.h"

#include "TFile.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "math.h"
#include "TRandom3.h"
#include "TLine.h"
#include <vector>
#include <algorithm>
#include <iostream>

class HitState : public TObject{
public:

    HitState();




    //    HitState(ParticleState *particle);


    //    HitState(TrdKHit* hit, AMSPoint* p0=0, AMSDir *dir=0);


    HitState(float xx,float yy,float zz);

    HitState(float xx[3]);

    //======For AMSRoot=======

    //    HitState(TrRecHitR *hit);

    //    HitState(TrRecHitR *hit,AMSEventR *ev);


    //    HitState(EcalHitR *hit);
    //    HitState(MCEventgR *hit);
    //    HitState(EcalShowerR *shower);



    ~HitState(){}

    bool operator() (HitState i,HitState j);



    //    void RegisterTrRecHit(TrRecHitR *hit);


    void Init();

    int FindLayer(float z);


    void Smear();


    bool IsRealTrTrackHit();
    bool IsTrTrackHit();
    bool IsTRDHit();


    SimpleHitState GetSimpleHitState();


    void Print(){


        std::cout<<"\t  Hit: " <<x<<", "<<y<<", "<<z<<", Amp: "<<AmpY<<",  res:"<<residule_x<<", "<<residule_y<<" err:"<<rx<<", "<<ry<<std::endl;


    }



public:
    //    static float Layer_z[10];
    double x,y,z;
    double tx,ty,p,q;
    double rx,ry,rtx,rty,rq;
    int layer_trk;
    int layer_trd;

    int hittype; // TrTrack:1 , ECALShower:2, TRDKHit:3
    int YOnly;
    int XOnly;
    int IsNull;


    float Amp;
    float AmpY;

    float Dist_Neibour;
    float N_Neibour;
    float Amp_Neibour;

    float pathlength;

    float residule_x;
    float residule_y;

    int hitno;





    HitState Merge(HitState hit2){

        if(AmpY != hit2.AmpY){
            std::cout<<"Error,  Not the same hit:"<<AmpY<<",   "<<hit2.AmpY<<std::endl;
            return *this;
        }


        // Keep XY Hit
        if(YOnly && !hit2.YOnly) return hit2;


        // Keep smallest residule_x
        if(fabs(residule_x)<fabs(hit2.residule_x))return *this;
        else return hit2;

    }


    ClassDef(HitState, 2);
};


class HitStateCollection{

public:

    HitStateCollection(){
        Init();
    }

    ~HitStateCollection(){
        Clear();
    }


    void Add(int i,HitState hit){
        Hits[i].push_back(hit);
        return;
    }

    int Size(int i){
        return Hits[i].size();
    }

    void Init(){
        for(int i=0;i<10;i++){
            Hits[i].clear();
            Hits[i].reserve(10);
        }
    }
    void Clear(){
        for(int i=0;i<10;i++)Hits[i].clear();
    }


    void CleanUp(int UseYInfoOnly);

    std::vector<HitState> Hits[10];

    std::vector<HitState> Hits_TRD;
    std::vector<HitState> Hits_ECAL;

    void Draw();

    int runid;
    int eventid;
};



class FittedState{

public:

    FittedState();
    FittedState(HitState hit, float Chi2y);


    FittedState(HitState hit, TMatrixD X,TMatrixD res, float Chi2y, TMatrixD R);


    bool IsNullHit(){return originalhit.IsNull;}
    bool IsTRDHit(){return originalhit.IsTRDHit();}
    bool IsTrTrackHit(){return originalhit.IsTrTrackHit();}
    bool IsRealTrTrackHit(){return originalhit.IsRealTrTrackHit();}


    ~FittedState(){}


    HitState originalhit;
    float Chi2;

    double ndof;

    double errx;
    double erry;

    double resx;
    double resy;
    int layer;
    double rig;
    double rigerr;

    float Chi2_x;
    float Chi2_y;
    bool isempty;

    float delta_theta;

    double z;
    double state[5];
    double cov[25];
    double err[5];
};

#endif
