#ifndef VERTEX3D_H
#define VERTEX3D_H

#include "TMinuit.h"
#include "TObject.h"
#include "Vertex2D.h"
#include "iostream"

class Vertex3D : public TObject{



public:

    Vertex3D(){
        Is2D=0;
        HasXZ=0;
        HasYZ=0;

        v_x = new Vertex2D;
        v_y = new Vertex2D;


    }
    ~Vertex3D(){
      delete v_x;
      delete v_y;
    }

    Vertex2D *v_x;
    Vertex2D *v_y;

    double x;
    double y;
    double z_x;
    double z_y;
    double z;
    double x_err;
    double y_err;
    double z_err;

    double Is2D;
    double HasXZ;
    double HasYZ;

    void SetX(Vertex2D* v2d){
        v_x=v2d;
        x=v_x->XY();
        z_x=v_x->Z();
        if(!v2d->GetNTracks()){
            Is2D=1;
            HasXZ=0;
        }else HasXZ=1;
    }

    void SetY(Vertex2D* v2d){
        v_y=v2d;
        y=v_y->XY();
        z_y=v_y->Z();
        if(!v2d->GetNTracks()){
            Is2D=1;
            HasYZ=1;
        }else HasYZ=1;
    }

    void AddTrack_XZ(Line* line){
        HasXZ=1;
        v_x->Add(line);
    }

    void AddTrack_YZ(Line* line){
        HasYZ=1;
        v_y->Add(line);
    }




    int HasSharedTracks(Vertex3D* v2){

        int nshared=0;

        for(int i=0;i<v_x->GetNTracks();i++){
            for(int j=0;j<v2->v_x->GetNTracks();j++){
                if(v_x->v_lines.at(i)==v2->v_x->v_lines.at(j))nshared++;
            }
        }

        for(int i=0;i<v_y->GetNTracks();i++){
            for(int j=0;j<v2->v_y->GetNTracks();j++){
                if(v_y->v_lines.at(i)==v2->v_y->v_lines.at(j))nshared++;
            }
        }

        return nshared;


    }
    double GetNTracksX(){
        return v_x->GetNTracks();
    }

    double GetNTracksY(){
        return v_y->GetNTracks();
    }


    double GetNTracks(){
        return v_x->GetNTracks()+v_y->GetNTracks();
    }

    void Print(){
        std::cout<<"====== Vertex ============"<<std::endl;
        std::cout<<"Position : "<<x<<", "<<y<<", "<<z<<std::endl;
        std::cout<<"Position Err: "<<x_err<<", "<<y_err<<", "<<z_err<<std::endl;
        std::cout<<"NTracks: " <<v_x->GetNTracks()+v_y->GetNTracks()<<std::endl;
        std::cout<<"Chi2 X: " <<GetChi2X()<<std::endl;
        for(int i=0;i<v_x->GetNTracks();i++){
            v_x->v_lines[i]->Print();
        }
        std::cout<<"Chi2 Y: " <<GetChi2Y()<<std::endl;
        for(int i=0;i<v_y->GetNTracks();i++){
            v_y->v_lines[i]->Print();
        }



    }

    static bool DEBUG;


    void RemoveOutliers(){
        std::vector<Line*>::iterator it;
        std::vector<Line*>::iterator it_max;
        while(v_x->v_lines.size()>=3 || v_y->v_lines.size()>=3){
            double chi2max=0;

            int removed=0;

            if(v_x->v_lines.size()>=3){
                for(it=v_x->v_lines.begin();it!=v_x->v_lines.end();it++){
                    double c2=(*it)->GetTotalChi2()/(*it)->GetNhit();
                    if(c2 > chi2max){
                        chi2max=c2;
                        it_max=it;
                    }
                }
                if(chi2max > 5){
                    v_x->v_lines.erase(it_max);
                    removed++;
                }
            }

            chi2max=0;
            if(v_y->v_lines.size()>=3){
                for(it=v_y->v_lines.begin();it!=v_y->v_lines.end();it++){
                    double c2=(*it)->GetTotalChi2()/(*it)->GetNhit();
                    if(c2 > chi2max){
                        chi2max=c2;
                        it_max=it;
                    }
                }
                if(chi2max > 5){
                    v_y->v_lines.erase(it_max);
                    removed++;
                }
            }

            if(removed){
                CombineFit();
            }else break;

        }
    }




    void CombineFit(){

        int ntrack_x=v_x->GetNTracks();
        int ntrack_y=v_y->GetNTracks();

        z=0.5*(z_x+z_y);

        if(DEBUG)std::cout<<"Set P0 XZ.........."<<std::endl;
        v_x->SetP0(x,z);
        if(DEBUG)std::cout<<"Set P0 YZ.........."<<std::endl;
        v_y->SetP0(y,z);

        int chi2_x=v_x->GetTotalChi2();
        int chi2_y=v_y->GetTotalChi2();

        int nhit_x=v_x->GetTotalNHits();
        int nhit_y=v_y->GetTotalNHits();

        //        std::cout<<"Before Fit: XZ: " <<chi2_x<<"/"<<nhit_x<<",   YZ: " <<chi2_y<<"/"<<nhit_y<<std::endl;

        //        if(ntrack_x + ntrack_y>=4){
        if(gMinuit){
            delete gMinuit;
            gMinuit=0;
        }
        gMinuit=new TMinuit(3+ntrack_x+ntrack_y);

        gMinuit->SetFCN(fcn_CombinedFit);
        gMinuit->SetObjectFit(this);
        double arglist[1];
        arglist[0] = -1;
        gMinuit->SetPrintLevel(-1);
        int ierr=0;
        gMinuit->mnparm(0,  "x0",  x,   1,  -380,380,ierr);
        gMinuit->mnparm(1,  "y0",  y,   1,  -380,380,ierr);
        gMinuit->mnparm(2,  "z0",  z,   1,  -380,380,ierr);

        for(int i=0;i<v_x->GetNTracks();i++){
            gMinuit->mnparm(3+i,  Form("slope_x_%i",i),  v_x->v_lines[i]->b,   0.1,  -10,10,ierr);  // This range need to set properly
        }
        for(int i=0;i<v_y->GetNTracks();i++){
            gMinuit->mnparm(3+v_x->GetNTracks()+i,  Form("slope_y_%i",i),  v_y->v_lines[i]->b,   0.1,  -10,10,ierr);  // This range need to set properly
        }

        arglist[0]=0;

        gMinuit->mnexcm("MIGRAD", arglist, 0,ierr);


        double bnd1, bnd2;
        int ivar;
        TString name;
        gMinuit->mnpout(0, name, x, x_err, bnd1, bnd2, ivar);
        gMinuit->mnpout(1, name, y, y_err, bnd1, bnd2, ivar);
        gMinuit->mnpout(2, name, z, z_err, bnd1, bnd2, ivar);

        double slope_x[20];
        double slope_err_x[20];
        for(int i=0;i<v_x->GetNTracks();i++){
            gMinuit->mnpout(i+3, name, slope_x[i], slope_err_x[i], bnd1, bnd2, ivar);
        }

        double slope_y[20];
        double slope_err_y[20];
        for(int i=0;i<v_y->GetNTracks();i++){
            gMinuit->mnpout(3+v_x->GetNTracks()+i, name, slope_y[i], slope_err_y[i], bnd1, bnd2, ivar);
        }

        v_x->SetP0andSlopes(x,z,slope_x);
        v_y->SetP0andSlopes(y,z,slope_y);


    }




    void CombineFit_New();








    static TMinuit *gMinuit;


    double GetChi2X(){
        //        return v_x->GetTotalChi2();
        double chi2=0;
        for(int i=0;i<v_x->GetNTracks();i++){
            chi2+=pow(v_x->v_lines.at(i)->GetDistance(x,z)/v_x->v_lines.at(i)->GetDistanceError(x,z),2);
        }
        return chi2;
    }
    double GetChi2Y(){
        //        return v_y->GetTotalChi2();
        double chi2=0;
        for(int i=0;i<v_y->GetNTracks();i++){
            chi2+=pow(v_y->v_lines.at(i)->GetDistance(y,z)/v_y->v_lines.at(i)->GetDistanceError(y,z),2);
        }
        return chi2;
    }


    double GetNhitX(){
        return v_x->GetTotalNHits();
    }
    double GetNhitY(){
        return v_y->GetTotalNHits();
    }


    static void fcn_CombinedFit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
    {
        Vertex3D *v=(Vertex3D*)gMinuit->GetObjectFit();
        f=v->fun_GetCombinedFitChi2(par);
    }


    double fun_GetCombinedFitChi2(Double_t *par)
    {
        v_x->SetP0andSlopes(par[0],par[2],par+3);
        v_y->SetP0andSlopes(par[1],par[2],par+3+v_x->GetNTracks());
        double chi2_x=v_x->GetTotalChi2();
        double chi2_y=v_y->GetTotalChi2();
        return chi2_x+chi2_y;
    }




    static void fcn_CombinedFit_New(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

    double fun_GetCombinedFitChi2_New(Double_t *par);





    static void fcn_3D_New(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
    {
        Vertex3D *v=(Vertex3D*)gMinuit->GetObjectFit();
        f=v->fun_GetCombinedFitChi2_3D(par);
    }

    double GetChi2_PointToLine_XZ(int i){
        double dist=v_x->v_lines.at(i)->GetDistance(x,z);
        double sigma=v_x->v_lines.at(i)->GetDistanceError(x,z);
        return pow(dist/sigma,2);
    }

    double GetChi2_PointToLine_YZ(int i){
        double dist=v_y->v_lines.at(i)->GetDistance(y,z);
        double sigma=v_y->v_lines.at(i)->GetDistanceError(y,z);
        return pow(dist/sigma,2);
    }


    double GetWeightedChi2_PointToLine_XZ(int i){
        double dist=v_x->v_lines.at(i)->GetDistance(x,z);
        double sigma=v_x->v_lines.at(i)->GetDistanceError(x,z);
        return v_x->v_lines.at(i)->weight*pow(dist/sigma,2);
    }

    double GetWeightedChi2_PointToLine_YZ(int i){
        double dist=v_y->v_lines.at(i)->GetDistance(y,z);
        double sigma=v_y->v_lines.at(i)->GetDistanceError(y,z);
        return v_y->v_lines.at(i)->weight*pow(dist/sigma,2);
    }


    double fun_GetCombinedFitChi2_3D(Double_t *par)
    {
        x=par[0];
        y=par[1];
        z=par[2];

        //        std::cout<<"Fitting par: " <<x<<", "<<y<<", "<<z<<std::endl;

        double sum_x=0;
        for(int i=0;i<v_x->GetNTracks();i++){
            //            std::cout<<"Chi2*weight xz: " <<i<<",  "<<GetChi2_PointToLine_XZ(i)<<",  weight: " <<v_x->v_lines[i]->weight<<std::endl;
            sum_x+=GetWeightedChi2_PointToLine_XZ(i);
        }
        double sum_y=0;
        for(int i=0;i<v_y->GetNTracks();i++){
            //            std::cout<<"Chi2*weight yz: " <<i<<",   "<<GetChi2_PointToLine_YZ(i)<<",  weight: " <<v_y->v_lines[i]->weight<<std::endl;
            sum_y+=GetWeightedChi2_PointToLine_YZ(i);
        }


        //        std::cout<<"Total Chi2: " <<sum_x+sum_y<<std::endl;

        return sum_x+sum_y;
    }


    int CombineFit_3D_Adaptive();








    void RemoveSmallWeight(std::vector<Line*> &v);
};






#endif // VERTEX3D_H
