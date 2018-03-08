#include "Vertex3D.h"
bool Vertex3D::DEBUG = 0;
TMinuit *Vertex3D::gMinuit = NULL;

using namespace std;


void Vertex3D::CombineFit_New(){


    if(gMinuit){
        delete gMinuit;
        gMinuit=0;
    }
    gMinuit=new TMinuit(3);

    gMinuit->SetFCN(fcn_CombinedFit_New);
    gMinuit->SetObjectFit(this);
    double arglist[1];
    arglist[0] = -1;
    gMinuit->SetPrintLevel(-1);
    int ierr=0;
    gMinuit->mnparm(0,  "x0",  x,   1,  -380,380,ierr);
    gMinuit->mnparm(1,  "y0",  y,   1,  -380,380,ierr);
    gMinuit->mnparm(2,  "z0",  z,   1,  -380,380,ierr);

    arglist[0]=0;

    gMinuit->mnexcm("MIGRAD", arglist, 0,ierr);


    double bnd1, bnd2;
    int ivar;
    TString name;
    gMinuit->mnpout(0, name, x, x_err, bnd1, bnd2, ivar);
    gMinuit->mnpout(1, name, y, y_err, bnd1, bnd2, ivar);
    gMinuit->mnpout(2, name, z, z_err, bnd1, bnd2, ivar);

    return;
}

void Vertex3D::fcn_CombinedFit_New(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Vertex3D *v=(Vertex3D*)gMinuit->GetObjectFit();
    f=v->fun_GetCombinedFitChi2_New(par);
}

double Vertex3D::fun_GetCombinedFitChi2_New(Double_t *par)
{
    double sum_x=0;
    double sum_y=0;
    //=====XZ
    for(int i=0;i<v_x->GetNTracks();i++){
        double dist=v_x->v_lines.at(i)->GetDistance(par[0],par[2]);
        double sigma=v_x->v_lines.at(i)->GetDistanceError(par[0],par[2]);
        sum_x+=pow(dist/sigma,2);
    }
    //    cout<<"Chi2 X:  "<<sum_x<<endl;


    //=====YZ
    for(int i=0;i<v_y->GetNTracks();i++){
        double dist=v_y->v_lines.at(i)->GetDistance(par[1],par[2]);
        double sigma=v_y->v_lines.at(i)->GetDistanceError(par[1],par[2]);
        sum_y+=pow(dist/sigma,2);
    }
    //    cout<<"Chi2 Y:  "<<sum_y<<endl;
    return sum_x+sum_y;
}

void Vertex3D::RemoveSmallWeight(vector<Line*> &v)
{
    vector<Line*>::iterator it;
    for(it=v.begin();it!=v.end();){
        if(IsSmallWeightLine(*it)){
            delete *it;
            it=v.erase(it);
        }else{
            it++;
        }
    }
}

int Vertex3D::CombineFit_3D_Adaptive()
{

    if(!HasYZ && !HasXZ)return -1;


    if(!v_x->v_lines.size()){
        Is2D=1;
        HasXZ=0;
    }
    if(!v_y->v_lines.size()){
        Is2D=1;
        HasYZ=0;
    }


    double T=1024;
    double chi2_c=5;

    if(!gMinuit)gMinuit=new TMinuit(3);
    gMinuit->SetFCN(fcn_3D_New);
    gMinuit->SetObjectFit(this);
    double arglist[1];
    arglist[0] = -1;
    gMinuit->SetPrintLevel(-1);

    x=0;
    y=0;
    z=0;

    if(HasXZ)v_x->SolveOnce();
    if(HasYZ)v_y->SolveOnce();

    if(HasXZ) x=v_x->x;
    if(HasYZ) y=v_y->x;

    if(HasXZ)z+=v_x->z;
    if(HasYZ)z+=v_y->z;
    z/=HasXZ+HasYZ;

    double scale=2;
    if(Is2D)scale=10;

    if(DEBUG)cout<<"Initial Position : "<<x<<", "<<y<<", "<<z<<endl;
    while(T>0.5){


        gMinuit->Release(0);
        gMinuit->Release(1);

        if(DEBUG){
            cout<<"Temperature: " <<T<<endl;
            cout<<"Position : "<<x<<", "<<y<<", "<<z<<endl;
            cout<<"Has XZ  : " <<HasXZ<<",   YZ: " <<HasYZ<<endl;
        }

        T/=scale;
        int ierr=0;
        gMinuit->mnparm(0,  "x0",  x,   1,  -380,380,ierr);
        gMinuit->mnparm(1,  "y0",  y,   1,  -380,380,ierr);
        gMinuit->mnparm(2,  "z0",  z,   1,  -380,380,ierr);

        if(!HasXZ)gMinuit->FixParameter(0);
        if(!HasYZ)gMinuit->FixParameter(1);

        arglist[0]=0;
        gMinuit->mnexcm("MIGRAD", arglist, 0,ierr);

        double bnd1, bnd2;
        int ivar;
        TString name;
        gMinuit->mnpout(0, name, x, x_err, bnd1, bnd2, ivar);
        gMinuit->mnpout(1, name, y, y_err, bnd1, bnd2, ivar);
        gMinuit->mnpout(2, name, z, z_err, bnd1, bnd2, ivar);

        if(DEBUG){
            cout<<"Fitted Position : "<<x<<", "<<y<<", "<<z<<endl;
            cout<<"Fitted Position Err: "<<x_err<<", "<<y_err<<", "<<z_err<<endl;
        }

        double wc=exp(-1*chi2_c/T);
        if(HasXZ){
            if(DEBUG) cout<<"========XZ Projection==================="<<endl;
            for(int i=0;i<v_x->GetNTracks();i++){
                double chi2=GetChi2_PointToLine_XZ(i);
                double w=exp(-1*chi2/T);
                v_x->v_lines.at(i)->weight=w/(w+wc);
                if(DEBUG) {

                    cout<<"Track : "<<i<<",  Chi2: " <<chi2 <<"   ";
                    v_x->v_lines.at(i)->Print();
                }
            }

            v_x->v_lines.erase( std::remove_if(v_x->v_lines.begin(), v_x->v_lines.end(), IsSmallWeightLine), v_x->v_lines.end() );

            if(v_x->GetNTracks()==1)v_x->v_lines.at(0)->weight=0;
            v_x->v_lines.erase( std::remove_if(v_x->v_lines.begin(), v_x->v_lines.end(), IsSmallWeightLine), v_x->v_lines.end() );

            if(!v_x->v_lines.size()){
                Is2D=1;
                HasXZ=0;
            }
        }


        if(HasYZ){
            if(DEBUG)cout<<"========YZ Projection==================="<<endl;
            for(int i=0;i<v_y->GetNTracks();i++){
                double chi2=GetChi2_PointToLine_YZ(i);
                double w=exp(-1*chi2/T);
                v_y->v_lines.at(i)->weight=w/(w+wc);
                if(DEBUG) {
                    cout<<"Track : "<<i<<",  Chi2: " <<chi2 <<"   ";
                    v_y->v_lines.at(i)->Print();
                }
            }
            v_y->v_lines.erase( std::remove_if(v_y->v_lines.begin(), v_y->v_lines.end(), IsSmallWeightLine), v_y->v_lines.end() );
            if(v_y->GetNTracks()==1)v_y->v_lines.at(0)->weight=0;
            v_y->v_lines.erase( std::remove_if(v_y->v_lines.begin(), v_y->v_lines.end(), IsSmallWeightLine), v_y->v_lines.end() );
            if(!v_y->v_lines.size()){
                Is2D=1;
                HasYZ=0;
            }

        }

        if(!v_x->GetNTracks() && v_y->GetNTracks()<=2) return -1;
        if(!v_y->GetNTracks() && v_x->GetNTracks()<=2) return -1;

    }

//    cout<<"NTrack : "<<v_x->GetNTracks()<<", "<<v_y->GetNTracks()<<",  "<<Is2D<<",  Position: " <<x<<", "<<y<<", "<<z<<endl;



    return 1;
}
