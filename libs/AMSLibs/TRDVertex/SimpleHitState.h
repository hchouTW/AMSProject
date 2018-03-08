#ifndef SIMPLEHITSTATE_H
#define SIMPLEHITSTATE_H

#include "TObject.h"
#include "vector"

class SimpleHitState: public TObject{

public:

    SimpleHitState(){
        x=-999;
        y=-999;
        z=-999;
        layer=-1;
        hittype=-1; // TrTrack:1 , ECALShower:2, TRDKHit:3
        YOnly=-1;
        XOnly=-1;
        Amp=-1;
        AmpY=-1;
        mult=0;
	//        memset(xmul,0,10*sizeof(double));
//        memset(ymul,0,10*sizeof(double));
//        memset(zmul,0,10*sizeof(double));
    }
    SimpleHitState(float _x,float _y, float _z, int _layer, int _hittype, int _YOnly, int _XOnly,
                   float _Amp, float _AmpY, float _path=0):x(_x),y(_y),z(_z),layer(_layer),hittype(_hittype),YOnly(_YOnly),XOnly(_XOnly),Amp(_Amp),AmpY(_AmpY),pathlength(_path){
        mult=0;
	//        memset(xmul,0,10*sizeof(double));
//        memset(ymul,0,10*sizeof(double));
//        memset(zmul,0,10*sizeof(double));
    }
    ~SimpleHitState(){}
    bool operator() (SimpleHitState i,SimpleHitState j) { return (i.z>j.z);}

public:

    double x,y,z;
    int layer;
    int hittype; // TrTrack:1 , ECALShower:2, TRDKHit:3
    int YOnly;
    int XOnly;
    float Amp;
    float AmpY;
    float pathlength;

    int mult;
    std::vector<float> xmul;
//    double ymul[10];
//    double zmul[10];


    ClassDef(SimpleHitState, 4)
};

#endif // SIMPLEHITSTATE_H
