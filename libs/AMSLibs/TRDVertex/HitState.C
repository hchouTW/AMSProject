#include "HitState.h"
#include "iostream"

bool CompareY (HitState i,HitState j) {
    return (i.y>j.y);
}


HitState::HitState(){
    Init();
    x=-999;
    y=-999;
    z=-999;
    layer_trk=-999;
    layer_trd=-999;

    hittype=0;
    IsNull=0;
}

//HitState::HitState(ParticleState *particle){
//    //        cout<<"Hit from particle:"<<particle->x<<", "<<particle->y<<", "<<particle->z<<endl;
//    Init();
//    x=particle->x;
//    y=particle->y;
//    z=particle->z;
//    layer=FindLayer(z);
//    Smear();
//}

//HitState::HitState(TrdKHit *hit, AMSPoint *p0, AMSDir *dir){
//    Init();
//    rx=0.3;
//    ry=0.3;
//    x=hit->TRDHit_x;
//    y=hit->TRDHit_y;
//    z=hit->TRDHit_z;
//    if(hit->TRDHit_Direction==1)  {
//        x=0;
//        YOnly=1;
//        rx=10000;
//    }
//    if(hit->TRDHit_Direction==0)  {
//        y=0;
//        XOnly=1;
//        ry=10000;
//    }
//    layer_trd=hit->TRDHit_Layer*-1;
//    hittype=3;

//    //=====path length=========

//    if(p0 && dir){
//        pathlength=hit->Tube_Track_3DLength(p0,dir);
//    }


//    Amp=hit->TRDHit_Amp;
//    AmpY=hit->TRDHit_GainCorrection;

//    return;
//}

HitState::HitState(float xx, float yy, float zz){
    Init();
    x=xx;
    y=yy;
    z=zz;
    layer_trk=FindLayer(z);
    hittype=1;


    //debug
    rx=1e-6;
    ry=1e-6;



}

HitState::HitState(float xx[]){
    Init();
    x=xx[0];
    y=xx[1];
    z=xx[2];
    layer_trk=FindLayer(z);
    Smear();
}

//HitState::HitState(TrRecHitR *hit){ // Assuming XY Hit
//    RegisterTrRecHit(hit);
//}

//HitState::HitState(TrRecHitR *hit, AMSEventR *ev){ // Assuming XY Hit
//    RegisterTrRecHit(hit);

//    N_Neibour=0;
//    Amp_Neibour=0;
//    Dist_Neibour=999;

//    if(!ev)return;
//    int nhit=ev->NTrRecHit();

//    for(int i=0;i<nhit;i++){
//        TrRecHitR *myhit=ev->pTrRecHit(i);
//        if(!myhit)continue;
//        if(myhit->GetLayerJ()!=hit->GetLayer())continue;
//        if(myhit->GetYCluster()==hit->GetYCluster())continue;
//        float dist=fabs(myhit->GetCoord().y()-hit->GetCoord().y());
//        if( dist> 0.1)continue;
//        N_Neibour++;
//        Amp_Neibour+=myhit->Sum();
//        if(dist<Dist_Neibour)Dist_Neibour=dist;
//    }

//}

//HitState::HitState(EcalHitR *hit){ // Assuming XY Hit
//    Init();

//    x=hit->Coo[0];
//    y=hit->Coo[1];
//    z=hit->Coo[2];
//    if(hit->Proj==1)  {
//        x=0;
//        YOnly=1;
//        rx=1000;
//    }else {
//        y=0;
//        XOnly=1;
//        ry=1000;
//    }
//    layer=hit->Plane;
//    hittype=3;

//    AmpY=hit->EdCorr;
//    Amp=hit->Edep;
//}

//HitState::HitState(MCEventgR *hit){ // Assuming XY Hit

//    Init();

//    x = hit->Coo[0];
//    y=hit->Coo[1];
//    z=hit->Coo[2];


//    layer=hit->Particle;
//    hittype=5;

//    AmpY=hit->Charge;
//    Amp=hit->Momentum;


//}

//HitState::HitState(EcalShowerR *shower){ // Assuming XY Hit
//    Init();
//    x=shower->CofG[0];
//    y=shower->CofG[1];
//    z=shower->CofG[2];
//    layer=10;
//    rx=0.5;
//    ry=0.5;

//    tx=shower->Dir[0]/shower->Dir[2];
//    ty=shower->Dir[1]/shower->Dir[2];
//    p=shower->GetCorrectedEnergy();


//    rtx=0.03;
//    rty=0.03;

//    q=-1/p;
//    rq=3*q;

//    //        cout<<"Hit from EcalShower...."<<endl;
//    //        cout<<"Dir:"<<shower->Dir[0]<<","<<shower->Dir[1]<<","<<shower->Dir[2]<<endl;
//    //        cout<<"Coo:"<<x<<", "<<y<<", "<<z<<endl;
//    //        cout<<"Energy : "<<p<<endl;

//    hittype=2;
//}

bool HitState::operator()(HitState i, HitState j) { return (i.z>j.z);}

//void HitState::RegisterTrRecHit(TrRecHitR *hit){

//    Init();
//    x=hit->GetCoord()[0];
//    y=hit->GetCoord()[1];
//    z=hit->GetCoord()[2];
//    if(hit->OnlyY() || hit->GetResolvedMultiplicity()==-1)  {
//        x=0;
//        YOnly=1;
//        rx=1000;
//    }
//    if(hit->OnlyX())  {
//        y=0;
//        XOnly=1;
//        ry=1000;
//    }
//    layer=hit->GetLayerJ();
//    hittype=1;

//    AmpY=hit->Sum();
//    Amp=hit->GetTotSignal();

//    return;
//}

void HitState::Init(){
    //        rx=3*0.003;
    //        ry=3*0.001;

    rx =  25.e-4;  // Tuned with muon MC 25/Jan/2011 (SH)
    ry =  13.e-4;  // Tuned with muon MC 25/Jan/2011 (SH)

    tx=0;
    ty=0;
    q=0.01;

    rtx=1;
    rty=1;
    rq=10;

    YOnly=0;
    XOnly=0;

    hittype=0;

    IsNull=1;

    hitno=-1;

    return;
}

int HitState::FindLayer(float z){
    static float Layer_z[10]={
        200,
        158.91,
        53.06,
        29.2,
        27.5,
        1.7,
        -1.7,
        -27.5,
        -29.2,
        -135.48,
    };
    for(int i=0;i<10;i++){
        if( fabs(z-Layer_z[i]) < 0.5) return i;
    }

    if(z>160)return 0;

    return -1;
}

void HitState::Smear(){
    x+=gRandom->Gaus(0,rx);
    y+=gRandom->Gaus(0,ry);
    return;
}

bool HitState::IsTrTrackHit(){
    return hittype==1;
}

bool HitState::IsRealTrTrackHit(){
    return IsTrTrackHit()&&!IsNull;
}

bool HitState::IsTRDHit(){
    return hittype==3;
}

SimpleHitState HitState::GetSimpleHitState(){
    return SimpleHitState(x,y,z,layer_trk, hittype, YOnly, XOnly,Amp, AmpY,pathlength);
}


FittedState::FittedState(){
    isempty=1;
    Chi2_x=-999;
    Chi2_y=-999;
    Chi2=-999;
    resx=0;
    resy=0;
    errx=10;
    erry=10;
    layer=0;
    originalhit.IsNull=1;
    delta_theta=0;
}

FittedState::FittedState(HitState hit, float Chi2y){
    originalhit=hit;

    if(hit.IsTrTrackHit())layer=hit.layer_trk;
    if(hit.IsTRDHit())layer=10+hit.layer_trd;

    Chi2=Chi2y;
    delta_theta=0;
}

FittedState::FittedState(HitState hit, TMatrixD X, TMatrixD res, float myChi2, TMatrixD R){
    isempty=0;
    originalhit=hit;
    z=hit.z;
    if(hit.IsTrTrackHit())layer=hit.layer_trk;
    if(hit.IsTRDHit())layer=10+hit.layer_trd;

    Chi2=myChi2;

    double *_state=X.GetMatrixArray();
    double *rr=res.GetMatrixArray();

    resx=rr[0];
    resy=rr[1];

    double* RR=R.GetMatrixArray();

    Chi2_x=0;
    Chi2_y=0;
    if(!hit.YOnly) Chi2_x= pow(resx,2)/RR[0];
    if(!hit.XOnly) Chi2_y= pow(resy,2)/RR[6];

    delta_theta=0;



//    memcpy(err,rr,5*sizeof(double));
    memcpy(state,_state,5*sizeof(double));



//    err[0]=sqrt(RR[0]);
//    err[1]=sqrt(RR[6]);
//    err[2]=sqrt(RR[12]);
//    err[3]=sqrt(RR[18]);
//    err[4]=sqrt(RR[24]);


    //    std::cout<<"State: " <<state[0]<<", "<<state[1]<<", "<<state[2]<<", "<<state[3]<<", "<<state[4]<<std::endl;



}


void HitStateCollection::CleanUp(int UseYInfoOnly)
{

    std::vector<HitState> Hits_XY[10];
    for(int ilayer=1;ilayer<10;ilayer++){

        // Sort in y
        std::sort(Hits[ilayer].begin(), Hits[ilayer].end(), CompareY);

        for(int i=0;i<Size(ilayer);i++){
            HitState hit=Hits[ilayer].at(i);
            if(!Hits_XY[ilayer].size()){
                Hits_XY[ilayer].push_back(Hits[ilayer].at(0));
            }else{
                bool ismatched=0;
                for(int j=0;j<Hits_XY[ilayer].size();j++){
                    HitState xyhit=Hits_XY[ilayer].at(j);
                    if(hit.AmpY==xyhit.AmpY){
                        HitState newhit=xyhit.Merge(hit);
                        Hits_XY[ilayer][j]=newhit;
                        ismatched=1;
                    }
                }
                if(!ismatched) Hits_XY[ilayer].push_back(hit);
            }
        }
        Hits[ilayer]=Hits_XY[ilayer];
    }



    if(UseYInfoOnly){
        for(int ilayer=1;ilayer<=9;ilayer++){
            for(int i=0;i<Size(ilayer);i++){
                HitState *hit=&(Hits[ilayer].at(i));
                hit->YOnly=1;
            }
        }
    }


    return;
}

void HitStateCollection::Draw()
{


    //    c->Divide(2,1);

    //    c->cd(1);
    //    gPad->DrawFrame(-60,-160,60,160);

    //    for(int i=1;i<10;i++){
    //        for(int j=0;j<Hits[i].size();j++){
    //            HitState* hit=&(Hits[i].at(j));
    //            if(hit->YOnly)continue;
    //            TEllipse* point=new TEllipse(hit->x,hit->z,0.5,0.5);
    //            point->Draw();
    //        }
    //    }

    ////    if(Hits_ECAL.size()){
    ////        HitState* hit=&(Hits_ECAL.at(0));
    ////        TEllipse* point=new TEllipse(hit->x,hit->z,0.5,0.5);
    ////        double x1,y1,x2,y2;
    ////        y1=-145;
    ////        y2=-160;
    ////        x1=hit->x+(y1-hit->z)*hit->tx;
    ////        x2=hit->x+(y2-hit->z)*hit->tx;
    ////        TLine* line=new TLine(x1,y1,x2,y2);
    ////        point->Draw();
    ////        line->SetLineWidth(3);
    ////        line->SetLineColor(kRed);
    ////        line->Draw();


    ////        y1=-145;
    ////        y2=160;
    ////        x1=hit->x+(y1-hit->z)*hit->tx;
    ////        x2=hit->x+(y2-hit->z)*hit->tx;
    ////        line=new TLine(x1,y1,x2,y2);
    ////        point->Draw();
    ////        line->SetLineWidth(3);
    ////        line->SetLineStyle(2);
    ////        line->SetLineColor(kRed);
    ////        line->Draw();


    ////    }
    ////    if(Hits_TRD.size()){
    ////        HitState* hit=&(Hits_TRD.at(0));
    ////        TEllipse* point=new TEllipse(hit->x,hit->z,0.5,0.5);
    ////        double x1,y1,x2,y2;
    ////        y1=150;
    ////        y2=90;
    ////        x1=hit->x+(y1-hit->z)*hit->tx;
    ////        x2=hit->x+(y2-hit->z)*hit->tx;
    ////        TLine* line=new TLine(x1,y1,x2,y2);
    ////        point->Draw();
    ////        line->SetLineWidth(3);
    ////        line->SetLineColor(kGreen+2);
    ////        line->Draw();
    ////        y1=150;
    ////        y2=-150;
    ////        x1=hit->x+(y1-hit->z)*hit->tx;
    ////        x2=hit->x+(y2-hit->z)*hit->tx;
    ////        line=new TLine(x1,y1,x2,y2);
    ////        line->SetLineStyle(2);
    ////        line->SetLineColor(kGreen+2);
    ////        point->Draw();
    ////        line->Draw();
    ////    }


    //    if(Hits_TRD.size() && Hits_ECAL.size()){
    //        HitState* trd=&(Hits_TRD.at(0));
    //        HitState* ecal=&(Hits_ECAL.at(0));
    //        double x1,y1,x2,y2;
    //        y1=trd->z;
    //        y2=ecal->z;
    //        x1=trd->x;
    //        x2=ecal->x;
    //        TLine* line=new TLine(x1,y1,x2,y2);
    //        line->SetLineWidth(3);
    //        line->SetLineColor(kBlue);
    //        line->Draw();
    //    }


    //    c->cd(2);
    //    gPad->DrawFrame(-60,-160,60,160);

    for(int i=1;i<10;i++){
        for(int j=0;j<Hits[i].size();j++){
            HitState* hit=&(Hits[i].at(j));
            if(hit->XOnly)continue;
            TEllipse* point=new TEllipse(hit->y,hit->z,0.5,0.5);
            point->Draw();
        }
    }
    if(Hits_ECAL.size()){
        HitState* hit=&(Hits_ECAL.at(0));
        TEllipse* point=new TEllipse(hit->y,hit->z,0.5,0.5);

        //        double x1,y1,x2,y2;
        //        y1=-145;
        //        y2=-160;
        //        x1=hit->y+(y1-hit->z)*hit->ty;
        //        x2=hit->y+(y2-hit->z)*hit->ty;
        //        TLine* line=new TLine(x1,y1,x2,y2);
        point->Draw();
        //        line->SetLineWidth(3);
        //        line->SetLineColor(kRed);
        //        line->Draw();

        //        y1=-145;
        //        y2=160;
        //        x1=hit->y+(y1-hit->z)*hit->ty;
        //        x2=hit->y+(y2-hit->z)*hit->ty;
        //        line=new TLine(x1,y1,x2,y2);
        //        point->Draw();
        //        line->SetLineWidth(3);
        //        line->SetLineStyle(2);
        //        line->SetLineColor(kRed);
        //        line->Draw();
    }

    if(Hits_TRD.size()){
        HitState* hit=&(Hits_TRD.at(0));
        TEllipse* point=new TEllipse(hit->y,hit->z,0.5,0.5);

        //        double x1,y1,x2,y2;
        //        y1=150;
        //        y2=90;
        //        x1=hit->y+(y1-hit->z)*hit->ty;
        //        x2=hit->y+(y2-hit->z)*hit->ty;
        //        TLine* line=new TLine(x1,y1,x2,y2);
        //        line->SetLineWidth(3);
        //        line->SetLineColor(kGreen+2);
        point->Draw();
        //        line->Draw();

        //        y1=150;
        //        y2=-150;
        //        x1=hit->y+(y1-hit->z)*hit->ty;
        //        x2=hit->y+(y2-hit->z)*hit->ty;
        //        line=new TLine(x1,y1,x2,y2);
        //        line->SetLineStyle(2);
        //        line->SetLineColor(kGreen+2);
        //        point->Draw();
        //        line->Draw();


    }



    if(Hits_TRD.size() && Hits_ECAL.size()){
        HitState* trd=&(Hits_TRD.at(0));
        HitState* ecal=&(Hits_ECAL.at(0));
        double x1,y1,x2,y2;
        y1=trd->z;
        y2=ecal->z;
        x1=trd->y;
        x2=ecal->y;
        TLine* line=new TLine(x1,y1,x2,y2);
        line->SetLineWidth(3);
        line->SetLineColor(kBlue);
        line->Draw();
    }





    return;
}
