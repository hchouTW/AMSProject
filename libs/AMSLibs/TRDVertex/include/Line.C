#include "Line.h"
using namespace std;
TMatrixD Line::I22 = TMatrixD(2,2);
TMinuit *Line::gMinuit = NULL;

bool IsSmallWeightLine(Line *a){
    return a->weight<1e-5;
}

double Line::LineFit(){

    int nhit = v_hit.size();
    if(nhit < 3) return -1;
    std::vector<double> res;
    res.resize(nhit);
    double mtx[2][2] = { { 0, 0 }, { 0, 0 } }, minv[2][2];
    double vec[2] = { 0, 0 };
    for (int i = 0; i < nhit; i++) {
        TRD2DHit* hit=v_hit.at(i);
        double w = (hit->y_err > 0) ? 1/hit->y_err/hit->y_err : 0;
        mtx[0][0] += w;      mtx[0][1] += w*hit->z;
        mtx[1][0] += w*hit->z; mtx[1][1] += w*hit->z*hit->z;
        vec[0]    += w*hit->y; vec[1]    += w*hit->y*hit->z;
    }
    double det = mtx[0][0]*mtx[1][1]-mtx[0][1]*mtx[1][0];
    if (det == 0) return -1;
    minv[0][0] = mtx[1][1]/det; minv[0][1] = -mtx[0][1]/det;
    minv[1][1] = mtx[0][0]/det; minv[1][0] = -mtx[1][0]/det;

    a = minv[0][0]*vec[0]+minv[0][1]*vec[1];
    b = minv[1][0]*vec[0]+minv[1][1]*vec[1];

    chi2 = 0;
    double SSxx=0;
    double SSzz=0;
    double xmean=0;
    double zmean=0;
    for (int i = 0; i < nhit; i++) {
        TRD2DHit* hit=v_hit.at(i);
        xmean+=hit->y/nhit;
        zmean+=hit->z/nhit;
    }

    for (int i = 0; i < nhit; i++) {
        TRD2DHit* hit=v_hit.at(i);
        SSxx=pow(hit->y-xmean,2);
        SSzz=pow(hit->z-zmean,2);
    }


    for (int i = 0; i < nhit; i++) {
        TRD2DHit* hit=v_hit.at(i);
        res[i] = hit->y-(a+b*hit->z);
        SSxx+=pow(res[i],2);
        chi2 += (  hit->y_err > 0) ? res[i]*res[i]/hit->y_err/hit->y_err : 0;
    }
    ndof = nhit-2;
    double normchi2 = chi2/ndof;

    b_err=sqrt(normchi2)/sqrt(SSzz);
    a_err=sqrt(normchi2)*sqrt(1/(nhit)+pow(zmean,2)/SSzz);
    if(b_err<1e-10)b_err=1e-10;



    SetMatrixRepresentation();
    return normchi2;
}

double Line::LineFit_TMinuit_Robust(){

    double T=10;
    double chi2_c=3;



    if(!gMinuit)gMinuit=new TMinuit(2);
    gMinuit->SetFCN(fcn_LineFit);
    gMinuit->SetObjectFit(this);
    double arglist[1];
    arglist[0] = -1;
    gMinuit->SetPrintLevel(-1);


    while(T>0.5){
        T/=2;
        int ierr=0;
        gMinuit->mnparm(0,  "a",  a,   1,  -1000,1000,ierr);
        gMinuit->mnparm(1,  "b",  b,   1,  -10,10,ierr);

        arglist[0]=0;
        gMinuit->mnexcm("MIGRAD", arglist, 0,ierr);

        double bnd1, bnd2;
        int ivar;
        TString name;
        gMinuit->mnpout(0, name, a, a_err, bnd1, bnd2, ivar);
        gMinuit->mnpout(1, name, b, b_err, bnd1, bnd2, ivar);


        double wc=exp(-1*chi2_c/T);
        for(int i=0;i<GetNhit();i++){
            double w=exp(-1*GetChi2(i)/T);
            v_hit.at(i)->weight=w/(w+wc);
        }
    }


    vector<TRD2DHit*>::iterator iter;
    for(iter=v_hit.begin();iter!=v_hit.end();){
        if((*iter)->weight<0.1)v_hit.erase(iter);
        else iter++;
    }

    ndof = v_hit.size() - 2;
    chi2=GetTotalChi2();
    SetMatrixRepresentation();

    return chi2;
}

void Line::fcn_LineFit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    Line *v=(Line*)gMinuit->GetObjectFit();
    f=v->fun_LineFit(par);
}

double Line::fun_LineFit(Double_t *par){
    a=par[0];
    b=par[1];
    return GetTotalWeightedChi2();
}

void Line::SetMatrixRepresentation(){
    //    _aa[0]=0;
    //    _aa[1]=a;
    //    _nn[0]=1/sqrt(1+b*b);
    //    _nn[1]=b/sqrt(1+b*b);
    //    aa.SetMatrixArray(_aa);
    //    nn.SetMatrixArray(_nn);
    //    nn_T.Transpose(nn);


    _aa[0]=a;
    _aa[1]=0;
    double norm=sqrt(1+b*b);

    _nn[0]=b/norm;
    _nn[1]=1/norm;

    aa.SetMatrixArray(_aa);
    nn.SetMatrixArray(_nn);
    nn_T.Transpose(nn);

    return;
}





void Line::SetLineAsPointAndSlopes(double _x0, double _z0, double slope){
    x0=_x0;
    z0=_z0;
    b=slope;
    a=x0-b*z0;
    SetMatrixRepresentation();
    return;
}

void Line::SetCrossPoint(double _x0, double _z0){
    x0=_x0;
    z0=_z0;
    a=x0-b*z0;
    SetMatrixRepresentation();
    return;
}


void Line::Draw(){
    TGraph* g=new TGraph;
    g->SetMarkerStyle(24);
    g->SetMarkerSize(2);

    double zmax=0;
    double zmin=999;


    for(int i=0;i<v_hit.size();i++){
        TRD2DHit *hit=v_hit.at(i);
        g->SetPoint(i,hit->y,hit->z);
        if(hit->z < zmin)zmin=hit->z;
        if(hit->z > zmax)zmax=hit->z;
    }
    g->Draw("p");








    TLine *line = new TLine(a+b*zmax,zmax,a+b*zmin,zmin);
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
    line->Draw("same");
}

Line::Line(int side){
    const double a_I22[2][2]={
        {1,0},
        {0,1},
    };
    I22.SetMatrixArray((const double*)a_I22);

    aa.ResizeTo(2,1);
    nn.ResizeTo(2,1);
    nn_T.ResizeTo(1,2);

    Init(side);
}

void Line::Init(int side)
{
    weight=1;
    iside=side;
    a=0;
    b=0;
    a_err=1;
    b_err=1;
    chi2=0;
    ndof=0;
    SetMatrixRepresentation();
}

Line::Line(){

    Init(0);
}

void Line::AddHit(TRD2DHit *hit){
    v_hit.push_back(hit);
    return;
}

int Line::Merge(Line *ll){
    // for the moment , this assume no shared hits
    bool merge=0;

    double newchi2_1=0;
    for(int i=0;i<ll->GetNhit();i++) {
        TRD2DHit* hit=ll->v_hit.at(i);
        newchi2_1+=pow(GetDistance(hit->y,hit->z)/hit->y_err,2);
    }

    double newchi2_2=0;
    for(int i=0;i<GetNhit();i++) {
        TRD2DHit* hit=v_hit.at(i);
        newchi2_2+=pow(ll->GetDistance(hit->y,hit->z)/hit->y_err,2);
    }
    newchi2_1/=ll->GetNhit();
    newchi2_2/=GetNhit();
    //    std::cout<<"Simple Chi2 test: " <<newchi2_1<<",  " <<newchi2_2<<std::endl;
    if(newchi2_1>3 && newchi2_2>3) return 0;


    Line newline(iside);
    for(int i=0;i<GetNhit();i++)    newline.AddHit(v_hit.at(i));
    for(int i=0;i<ll->GetNhit();i++)    newline.AddHit(ll->v_hit.at(i));
    newline.LineFit_TMinuit();
    double newchi2=newline.GetTotalChi2();
    if(newchi2/GetNhit()<5)merge=1;
    //    std::cout<<"Full Chi2 test: " <<newchi2/GetNhit()<<std::endl;

    if(merge){
        merge=1;
        for(int i=0;i<ll->GetNhit();i++){
            v_hit.push_back(ll->v_hit[i]);
        }
        LineFit_TMinuit_Robust();
        return 1;
    }else return 0;
}

double Line::LineFit_TMinuit(){

    if(!gMinuit)gMinuit=new TMinuit(2);
    gMinuit->SetFCN(fcn_LineFit);
    gMinuit->SetObjectFit(this);
    double arglist[1];
    arglist[0] = -1;
    gMinuit->SetPrintLevel(-1);
    int ierr=0;
    gMinuit->mnparm(0,  "a",  a,   1,  -1000,1000,ierr);
    gMinuit->mnparm(1,  "b",  b,   1,  -10,10,ierr);

    arglist[0]=0;
    gMinuit->mnexcm("MIGRAD", arglist, 0,ierr);

    double bnd1, bnd2;
    int ivar;
    TString name;
    gMinuit->mnpout(0, name, a, a_err, bnd1, bnd2, ivar);
    gMinuit->mnpout(1, name, b, b_err, bnd1, bnd2, ivar);

    //        cout<<"Line Fit Result : " <<a<<"+/-"<<a_err<<",  b: " <<b<<"+/-"<<b_err<<", Chi2: " <<GetTotalChi2()<<",  nhit: " <<GetNhit()<<endl;

    SetMatrixRepresentation();
    return GetTotalChi2();
}

double Line::GetDistance(double x, double z){
    //x=a+bz;
    //a+bz+(-1)x=0;
    double distance=fabs(a+b*z-x)/sqrt(b*b+1);

    // TMatrixD ap(2,1);
    // double xx[2]={z,x};
    // ap.SetMatrixArray(xx);
    // double distance2= GetDistance(ap);

    // cout<<"distance 2: " <<distance2<<",  simple: " <<distance<<endl;

    return distance;

}

double Line::GetDistanceError(double x, double z){
    double sigma_r=sqrt(a_err*a_err/(1+b*b) + b_err*b_err*pow(-1*a*b+b*x+z,2)/pow(1+b*b,3));
    //    cout<<"Distance: " <<GetDistance(x,z)<<", err: " <<sigma_r<<endl;
    return sigma_r;
}

double Line::GetWeightedChi2(int ihit){
    if(ihit>GetNhit())return 0;
    return  v_hit.at(ihit)->weight*GetChi2(ihit);
}

double Line::GetChi2(int ihit){
    if(ihit>GetNhit())return 0;
    double dist=v_hit.at(ihit)->y-(a+b*v_hit.at(ihit)->z);
    return  dist*dist/(v_hit.at(ihit)->y_err*v_hit.at(ihit)->y_err);
}

double Line::GetTotalWeightedChi2(){
    chi2=0;
    for(int i=0;i<v_hit.size();i++){
        chi2+=GetWeightedChi2(i);
    }
    return chi2;
}

double Line::GetTotalChi2(){
    chi2=0;
    for(int i=0;i<v_hit.size();i++){
        chi2+=GetChi2(i);
    }
    return chi2;
}

void Line::Print(){
    std::cout<<"Line: " <<a<<", "<<b<<", "<<a_err<<",  "<<b_err<<"  nhit: " <<GetNhit()<<",  Chi2: " <<GetTotalChi2()<<",  wweight: " <<weight<<std::endl;
    for(int i=0;i<GetNhit();i++){
        v_hit.at(i)->Print();
        //        std::cout<<"Hit: " <<i<<": chi2: "<<GetChi2(i)<<std::endl;
    }
    //        std::cout<<"Matrix: " <<std::endl;
    //        std::cout<<" aa: " ; aa.Print();
    //        std::cout<<" nn: " ; nn.Print();
}

double Line::GetDistance(TMatrixD p){
    TMatrixD ap(2,1);
    ap=aa-p;
    TMatrixD ap_T(1,2);
    ap_T.Transpose(ap);
    TMatrixD D=ap_T*(I22-nn*nn_T)*ap;
    double dist2=D.GetMatrixArray()[0];
    if(dist2<1e-6)dist2=1e-6;
    return sqrt(dist2);
}


