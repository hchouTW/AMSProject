#include "TRD2DHit.h"

using namespace std;

void TRD2DHit::Init()
{
    weight=1;
    y=0;
    z=0;
    y_err=0.3;
    z_err=0.3;
    Amp=0;
    layer=0;
    ismergedhit=0;
    boundary[0]=0;
    boundary[1]=0;
    connect_upper_mindist=999;
    connect_lower_mindist=999;
    type=-1;
    type=-1;

    v_connected_upper.clear();
    v_connected_lower.clear();

    v_connected_upper.resize(0);
    v_connected_lower.resize(0);

    //    v_connectedhit_upper.clear();
    //    v_connectedhit_lower.clear();
    //    v_connectedhit_upper.resize(0);
    //    v_connectedhit_lower.resize(0);


}

TRD2DHit::TRD2DHit(){
    Init();
}

TRD2DHit::TRD2DHit(float _x, float _z, float _Amp, int _layer){
    Init();
    y=_x;
    z=_z;
    Amp=_Amp;
    layer=_layer;
}

TRD2DHit::TRD2DHit(SimpleHitState simple){
    Init();
    if(simple.XOnly){
        y=simple.x;
        type=0;
    }else {
        y=simple.y;
        type=1;
    }
    z=simple.z;
    layer=fabs(simple.layer);
    Amp=simple.Amp;
    boundary[0]=y-0.3;
    boundary[1]=y+0.3;

    y_err=0.3;
    z_err=0.3;

}




TRD2DHit::TRD2DHit(TrdRawHitR *raw){
    Init();
    layer=raw->Layer;
    TRDHitRZD rzd(*raw);
    y=rzd.r;
    z=rzd.z;

    if(layer<=3 || layer>=16){
        type=0;
    }else {
        type=1;
    }
    Amp=raw->Amp;
    boundary[0]=y-0.3;
    boundary[1]=y+0.3;

    y_err=0.3;
    z_err=0.3;

}







bool TRD2DHit::operator()(TRD2DHit i, TRD2DHit j){
    return i.z < j.z;
}


bool TRD2DHit::operator()(TRD2DHit* i, TRD2DHit *j){
    return i->z < j->z;
}

float TRD2DHit::GetDistance(TRD2DHit *p2){
    return sqrt(pow(y-p2->y,2)+pow(z-p2->z,2));
}

float TRD2DHit::GetDistance_X(TRD2DHit *p2){
    return fabs(y-p2->y);
}

float TRD2DHit::GetDistance_Z(TRD2DHit *p2){
    return fabs(z-p2->z);
}

void TRD2DHit::AddConnectedUpper(Cell* cell)
{
    v_connected_upper.push_back(cell);
}

void TRD2DHit::AddConnectedLower(Cell* cell)
{
    v_connected_lower.push_back(cell);
}


void TRD2DHit::AddConnectedHitUpper(TRD2DHit* hit)
{
    //    v_connectedhit_upper.push_back(hit);
}

void TRD2DHit::AddConnectedHitLower(TRD2DHit* hit)
{
    //    v_connectedhit_lower.push_back(hit);
}




bool TRD2DHit::MergeWith(TRD2DHit *hit){
    //    return 0;
    if(layer!=hit->layer)return 0;
    if(z!=hit->z)return 0;
    double distance = hit->y-y;
    if(distance>0 && fabs(boundary[1]-hit->boundary[0])<0.2){
        boundary[1]=hit->boundary[1];
        y=(Amp*y+hit->Amp*hit->y)/(Amp+hit->Amp);
        y_err=fabs(y-boundary[1]);
        Amp+=hit->Amp;
        ismergedhit=1;
        return 1;
    }else if(distance<0 && fabs(boundary[0]-hit->boundary[1])<0.2){
        boundary[0]=hit->boundary[0];
        y=(Amp*y+hit->Amp*hit->y)/(Amp+hit->Amp);
        y_err=fabs(y-boundary[0]);
        Amp+=hit->Amp;
        ismergedhit=1;
        return 1;
    }


    return 0;
}

void TRD2DHit::ClearConnections()
{
    v_connected_lower.clear();
    v_connected_upper.clear();
    //    v_connectedhit_lower.clear();
    //    v_connectedhit_upper.clear();
}

void TRD2DHit::Print()
{
    std::cout<<"TRD 2D Hit: " <<y<<", "<<z<<", layer: " <<layer<<", type: " <<type<<std::endl;
}

void TRD2DHitCollection::Add(SimpleHitState simplehit){
    Add(new TRD2DHit(simplehit));
    return;
}

void TRD2DHitCollection::Add(TRD2DHit *hit){
    Hits[abs(hit->layer)].push_back(hit);
    NTotalHits++;
    return;
}

void TRD2DHitCollection::Add_WidthMerging(SimpleHitState simplehit){
    TRD2DHit* hit=new TRD2DHit(simplehit);
    Add_WidthMerging(hit);
    return;
}

void TRD2DHitCollection::Add_WidthMerging(TRD2DHit *hit){
    for(int i=0;i<Hits[abs(hit->layer)].size();i++){
        TRD2DHit* existinghit=(Hits[abs(hit->layer)].at(i));
        int result=existinghit->MergeWith(hit);
        if(result){
            delete hit;
            return;
        }
    }

    Hits[abs(hit->layer)].push_back(hit);
    NTotalHits++;
    return;
}

int TRD2DHitCollection::Size(int i){
    return Hits[i].size();
}

void TRD2DHitCollection::Init(){
    for(int i=0;i<NLayers;i++){
        Hits[i].clear();
    }
    NTotalHits=0;
}

void TRD2DHitCollection::ClearMemory(){
    for(int ilayer=0;ilayer<NLayers;ilayer++){
        for(int ihit=0;ihit<Hits[ilayer].size();ihit++)delete Hits[ilayer].at(ihit);
        Hits[ilayer].clear();
    }
    NTotalHits=0;
}


void TRD2DHitCollection::Draw(){
    for(int ilayer=0;ilayer<NLayers;ilayer++){
        if(!Hits[ilayer].size())continue;
        for(int ihit=0;ihit<Hits[ilayer].size();ihit++){
            TRD2DHit* hit=(Hits[ilayer].at(ihit));
            TEllipse* point=new TEllipse(hit->y,hit->z,0.3,0.3);
            int color=kRed;
            point->SetFillColor(color);
            point->Draw();
        }
    }
}



void TRD2DHitCollection::Draw_YZ(){

    for(int ilayer=0;ilayer<NLayers;ilayer++){
        if(!Hits[ilayer].size())continue;
        if(ilayer>3 && ilayer<16)continue;
        for(int ihit=0;ihit<Hits[ilayer].size();ihit++){
            TRD2DHit* hit=(Hits[ilayer].at(ihit));
            TEllipse* point=new TEllipse(hit->y,hit->z,hit->y_err,0.3);
            int color=TColor::GetColor(log10(hit->Amp)/log10(5000),float(0.),float(0.));
            point->SetFillColor(color);
            point->Draw();
        }
    }


}



void TRD2DHitCollection::Draw_XZ(){
    for(int ilayer=0;ilayer<NLayers;ilayer++){
        if(!Hits[ilayer].size())continue;
        if(ilayer<=3 || ilayer>=16)continue;
        for(int ihit=0;ihit<Hits[ilayer].size();ihit++){
            TRD2DHit* hit=(Hits[ilayer].at(ihit));
            TEllipse* point=new TEllipse(hit->y,hit->z,hit->y_err,0.3);
            int color=TColor::GetColor(log10(hit->Amp)/log10(5000),float(0.),float(0.));
            point->SetFillColor(color);
            point->Draw();
        }
    }
    return;
}

TRD2DHit *TRD2DHitCollection::GetHit(int ilayer, int index){
    if(ilayer<0 || ilayer>NLayers)return 0;
    if(index>=Hits[ilayer].size())return 0;
    return (Hits[ilayer].at(index));
}

void TRD2DHitCollection::Sort(){
    TRD2DHit dummy;
    for(int ilayer=0;ilayer<NLayers;ilayer++){
        std::sort(Hits[ilayer].begin(),Hits[ilayer].end(),dummy);
        std::reverse(Hits[ilayer].begin(),Hits[ilayer].end());
    }
}
