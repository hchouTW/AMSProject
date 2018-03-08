#include "Vertex2D.h"

TMinuit *VertexFinder2D_ZVTOP::gMinuit = NULL;
bool VertexFinder2D_ZVTOP::DEBUG = 0;

bool Vertex2D::DEBUG = 0;
TMinuit *Vertex2D::gMinuit = NULL;
TMatrixD Vertex2D::I22 = TMatrixD(2,2);

using namespace std;





bool comp2DVertexByProbability(Vertex2D *a, Vertex2D *b){
    return a->V>b->V;
}

bool Vertex2D::Contains(Line *line){
    for(int i=0;i<GetNTracks();i++){
        if(v_lines[i]==line) return 1;;
    }
    return 0;
}

int Vertex2D::Merge(Vertex2D *vv){
    static int cut_maxdist=15;

    if(V >1 && vv->V<=1){
        // remove 2-track vertex which have shared track with current vertex
        int shared=0;
        for(int i=0;i<GetNTracks();i++){
            for(int j=0;j<vv->GetNTracks();j++){
                if(v_lines[i] == vv->v_lines[j])shared=1;
            }

        }

        if(shared)return 1;
    }




    if(ContainsAll(vv)){
        return 1;
    }else if(BelongsTo(vv)){
        for(int i=0;i<vv->GetNTracks();i++){
            Add(vv->v_lines.at(i));
        }
        SolveOnce();
        //CombineFit();
        return 1;
    }else{
        if(Distance(vv)>cut_maxdist)return 0;

        double totalchi2=GetTotalChi2();
        double totalhit=GetTotalNHits();

        std::vector<Line*> v_lines_original=v_lines;
        for(int i=0;i<vv->GetNTracks();i++)Add(vv->v_lines.at(i));

        SolveOnce();
        //CombineFit();

        double totalchi2_new=GetTotalChi2();
        double totalhit_new=GetTotalNHits();

        // cout<<"Chi2 Change from  : "<<totalchi2<<"  to   "<<totalchi2_new<<endl;
        //  cout<<"Total Nhit chagne from  : "<<totalhit <<" to " <<totalhit_new<<endl;


        if((totalchi2_new-totalchi2)/(totalhit_new-totalhit) > 3){
            v_lines = v_lines_original;
            // SolveOnce();
            //CombineFit();
            return 0;
        }


        return 1;


        //            std::vector<Line*>v_lines_original=v_lines;

        //            for(int i=0;i<vv->GetNTracks();i++)Add(vv->v_lines.at(i));
        //            SolveOnce();
        //            CombineFit();
        //            if(GetTotalChi2()/GetTotalNHits()<5)return 1;
        //            else {
        //                v_lines=v_lines_original;
        //                SolveOnce();
        //                return 0;
        //            }


    }
}


void Vertex2D::Print(){
    std::cout<<"Vertex 2D: " <<XY()<<", "<<Z()<<",  ntracks: " <<GetNTracks()<<",   Chi2: " <<GetTotalChi2()<<"/"<<GetTotalNHits()<<",  Probability: " <<V<<std::endl;
    for(int i=0;i<v_lines.size();i++){
        std::cout<<"\t";v_lines[i]->Print();
        std::cout<<"\t Distance: " <<v_lines[i]->GetDistance(XY(),Z())<<std::endl;
    }
}

bool Vertex2D::ContainsAll(Vertex2D *vv){
    for(int i=0;i<vv->GetNTracks();i++){
        bool matched=0;
        for(int j=0;j<GetNTracks();j++){
            if(vv->v_lines[i]==v_lines[j])matched=1;
        }
        if(!matched)return 0;
    }
    return 1;
}

bool Vertex2D::BelongsTo(Vertex2D *vv){
    for(int i=0;i<GetNTracks();i++){
        bool matched=0;
        for(int j=0;j<vv->GetNTracks();j++){
            if(v_lines[i]==vv->v_lines[j])matched=1;
        }
        if(!matched)return 0;
    }
    return 1;
}

double Vertex2D::GetTotalNHits(){
    nhits=0;
    for(int i=0;i<GetNTracks();i++){
        nhits+=v_lines[i]->GetNhit();
    }
    return nhits;
}

double Vertex2D::GetTotalChi2(){
    chi2=0;
    for(int i=0;i<GetNTracks();i++){
        chi2+=v_lines[i]->GetTotalChi2();
    }
    return chi2;
}

void Vertex2D::SetP0(double x, double z){
    if(DEBUG)std::cout<<"line size: " <<GetNTracks()<<std::endl;
    for(int i=0;i<GetNTracks();i++){
        if(DEBUG)std::cout<<"i : " <<i<<std::endl;
        v_lines[i]->SetCrossPoint(x,z);
    }
    return;
}

void Vertex2D::SetP0andSlopes(double x, double z, double *slope){
    for(int i=0;i<GetNTracks();i++){
        v_lines[i]->SetLineAsPointAndSlopes(x,z,slope[i]);
    }
    return;
}

void Vertex2D::Draw(){
    cout<<"vertex 2D,  draw   ntracks : " <<GetNTracks()<<endl;
    for(int i=0;i<GetNTracks();i++){
        v_lines[i]->Print();
        v_lines[i]->Draw();
    }
}

void Vertex2D::SolveOnce(){



    TMatrixD R(2,2);
    TMatrixD R_inv(2,2);
    TMatrixD q(2,1);

    for(int i=0;i<GetNTracks();i++){
        double c=1./pow(v_lines[i]->b_err,2);  // use error of slope as weight
        TMatrixD nn=v_lines[i]->nn;
        TMatrixD nn_T(1,2);
        nn_T.Transpose(nn);
        R+=c*(I22-nn*nn_T);
        q+=c*(I22-nn*nn_T)*(v_lines[i]->aa);
    }
    R_inv=R;
    R_inv.Invert();
    p=R_inv*q;

//    z=p.GetMatrixArray()[0];
//    x=p.GetMatrixArray()[1];

    x=p.GetMatrixArray()[0];
    z=p.GetMatrixArray()[1];



    return;
}

void Vertex2D::Solve(){
    int maxiter=5;
    int niter=0;
    while(niter<maxiter){
        SolveOnce();
        std::vector<Line*>::iterator it=v_lines.begin();
        std::vector<Line*>::iterator it_maxdist;
        double dist_max=0;
        for(;it!=v_lines.end();it++){
            double distance=(*it)->GetDistance(p);
            if(distance>dist_max){
                dist_max=distance;
                it_maxdist=it;
            }
        }
        if(dist_max>2){
            it=v_lines.erase(it_maxdist);
            niter++;
        }else{
            break;
        }
    }
    //        p.Print();
    return;
}

void Vertex2D::Add(Line *line){
    for(int i=0;i<GetNTracks();i++){
        if(v_lines[i]==line){
            return;
        }
    }
    v_lines.push_back(line);
}

Vertex2D::Vertex2D(){
    v_lines.reserve(100);
    const double a_I22[2][2]={
        {1,0},
        {0,1},
    };
    I22.SetMatrixArray((const double*)a_I22);
    p.ResizeTo(2,1);
    double pp[2]={0,0};
    p.SetMatrixArray(pp);
    has3D=0;
    isremoved=0;
    weight=1;
}

double Vertex2D::dir(){
    double direction=0;
    for(int i=0;i<GetNTracks();i++){
        Line* line=v_lines.at(i);
        for(int j=0;j<v_lines[i]->GetNhit();j++){
            direction+=line->v_hit.at(i)->z-Z();
        }
    }
    return direction/fabs(direction);
}

void Vertex2D::fcn_CombinedFit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Vertex2D *v=(Vertex2D*)gMinuit->GetObjectFit();
    f=v->fun_GetCombinedFitChi2(par);
}

double Vertex2D::fun_GetCombinedFitChi2(Double_t *par)
{
    SetP0andSlopes(par[0],par[1],par+2);
    return GetTotalChi2();
}

void Vertex2D::RemoveOutliers(){
    std::vector<Line*>::iterator it;
    std::vector<Line*>::iterator it_max;
    while(v_lines.size()>3){
        double chi2max=0;
        for(it=v_lines.begin();it!=v_lines.end();it++){
            double c2=(*it)->GetTotalChi2()/(*it)->GetNhit();
            if(c2 > chi2max){
                chi2max=c2;
                it_max=it;
            }
        }
        if(chi2max > 5){
            v_lines.erase(it_max);
            CombineFit();
            continue;
        }else break;
    }
}

void Vertex2D::CombineFit(){
    int ntrack=GetNTracks();
    z=Z();
    x=XY();

    if(gMinuit){
        delete gMinuit;
        gMinuit=0;
    }
    gMinuit=new TMinuit(2+ntrack);

    gMinuit->SetFCN(fcn_CombinedFit);
    gMinuit->SetObjectFit(this);
    double arglist[1];
    arglist[0] = -1;
    gMinuit->SetPrintLevel(-1);
    int ierr=0;
    gMinuit->mnparm(0,  "x0",  x,   1,  -380,380,ierr);
    gMinuit->mnparm(1,  "z0",  z,   1,  -380,380,ierr);

    for(int i=0;i<GetNTracks();i++){
        gMinuit->mnparm(2+i,  Form("slope_x_%i",i),  v_lines[i]->b,   0.1,  -10,10,ierr);  // This range need to set properly
    }

    arglist[0]=0;

    gMinuit->mnexcm("MIGRAD", arglist, 0,ierr);


    double bnd1, bnd2;
    int ivar;
    TString name;
    gMinuit->mnpout(0, name, x, x_err, bnd1, bnd2, ivar);
    gMinuit->mnpout(1, name, z, z_err, bnd1, bnd2, ivar);

    double slope[20];
    double slope_err_x[20];
    for(int i=0;i<GetNTracks();i++){
        gMinuit->mnpout(i+2, name, slope[i], slope_err_x[i], bnd1, bnd2, ivar);
    }

    SetP0andSlopes(x,z,slope);


    double pp[2]={z,x};
    p.SetMatrixArray(pp);


    //        std::cout<<"Fit Result: " <<x<<", "<<y<<", "<<z<<std::endl;
    //        std::cout<<"NTracks: " <<v_x.GetNTracks()+v_y.GetNTracks()<<std::endl;
    //        std::cout<<"Chi2 X: " <<v_x.GetTotalChi2()<<"/"<<v_x.GetTotalNHits()<<std::endl;
    //        std::cout<<"Chi2 Y: " <<v_y.GetTotalChi2()<<"/"<<v_y.GetTotalNHits()<<std::endl;

    delete gMinuit;
    gMinuit=0;


}

double Vertex2D::fun_GetCombinedFitChi2_New(Double_t *par)
{
    double sum_x=0;
    for(int i=0;i<GetNTracks();i++){
        double dist=v_lines.at(i)->GetDistance(par[0],par[1]);
        double sigma=v_lines.at(i)->GetDistanceError(par[0],par[1]);
        sum_x+=v_lines.at(i)->weight*pow(dist/sigma,2);
    }
    return sum_x;
}

void Vertex2D::fcn_CombinedFit_New(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Vertex2D *v=(Vertex2D*)gMinuit->GetObjectFit();
    f=v->fun_GetCombinedFitChi2_New(par);
}

void Vertex2D::CombineFit_New(){

    if(gMinuit){
        delete gMinuit;
        gMinuit=0;
    }
    gMinuit=new TMinuit(2);

    gMinuit->SetFCN(fcn_CombinedFit_New);
    gMinuit->SetObjectFit(this);
    double arglist[1];
    arglist[0] = -1;
    gMinuit->SetPrintLevel(-1);
    int ierr=0;
    gMinuit->mnparm(0,  "x0",  x,   1,  -380,380,ierr);
    gMinuit->mnparm(2,  "z0",  z,   1,  -380,380,ierr);

    arglist[0]=0;

    gMinuit->mnexcm("MIGRAD", arglist, 0,ierr);


    double bnd1, bnd2;
    int ivar;
    TString name;
    gMinuit->mnpout(0, name, x, x_err, bnd1, bnd2, ivar);
    gMinuit->mnpout(2, name, z, z_err, bnd1, bnd2, ivar);

    return;
}

void Vertex2D::CombineFit_New_Adaptive(){



    double T=10;
    double chi2_c=3;

    if(!gMinuit)gMinuit=new TMinuit(2);
    gMinuit->SetFCN(fcn_CombinedFit_New);
    gMinuit->SetObjectFit(this);
    double arglist[1];
    arglist[0] = -1;
    gMinuit->SetPrintLevel(-1);

    SolveOnce();

//    cout<<"Initial Position : "<<x<<", "<<z<<endl;
    while(T>0.1){
//        cout<<"Temperature: " <<T<<endl;
//        cout<<"Position : "<<x<<", "<<z<<endl;
        T/=2;
        int ierr=0;
        gMinuit->mnparm(0,  "x0",  x,   1,  -380,380,ierr);
        gMinuit->mnparm(1,  "z0",  z,   1,  -380,380,ierr);

        arglist[0]=0;
        gMinuit->mnexcm("MIGRAD", arglist, 0,ierr);

        double bnd1, bnd2;
        int ivar;
        TString name;
        gMinuit->mnpout(0, name, x, x_err, bnd1, bnd2, ivar);
        gMinuit->mnpout(1, name, z, z_err, bnd1, bnd2, ivar);


        double wc=exp(-1*chi2_c/T);
        for(int i=0;i<GetNTracks();i++){

            double dist=v_lines.at(i)->GetDistance(x,z);
            double sigma=v_lines.at(i)->GetDistanceError(x,z);
            double chi2=pow(dist/sigma,2);

            double w=exp(-1*chi2/T);
            v_lines.at(i)->weight=w/(w+wc);


            v_lines.at(i)->Print();
        }
    }

    v_lines.erase( std::remove_if(v_lines.begin(), v_lines.end(), IsSmallWeightLine), v_lines.end() );
    //    vector<Line*>::iterator iter;
    //    for(iter=v_lines.begin();iter!=v_lines.end();){
    //        if((*iter)->weight<0.1)v_lines.erase(iter);
    //        else iter++;
    //    }
    return;

}





void VertexFinder2D_ZVTOP::FindBestVertex(){

    while(Niteration<5){
        if(DEBUG)cout<<"Start iteration: " <<Niteration<<",    Ntracks: " <<v_track->size()<<endl;
        Niteration++;
        double Vmax=1;
        Vertex2D *candidate;
        for(int i=0;i<v_track->size();i++){
            for(int j=i+1;j<v_track->size();j++){
                Vertex2D *v = new Vertex2D;
                v->iside=iside;
                v->Add((v_track->at(i)));
                v->Add((v_track->at(j)));
                v->SolveOnce();
                double xx[2]={v->XY(),v->Z()};
                //  v->Print();
                v->V=fun_ZVTOP_V(xx,0);
                if(fun_ZVTOP_V(xx,0)>Vmax){
                    Vmax=fun_ZVTOP_V(xx,0);
                    candidate=v;
                    x=v->XY();
                    z=v->Z();
                }else {
                    delete v;
                }
            }
        }
        if(Vmax>1){
            vector<Line*>::iterator iter_line;
            for(iter_line=v_track->begin();iter_line!=v_track->end();iter_line++){
                if((*iter_line)->GetDistance(candidate->XY(),candidate->Z()) < 10){
                    candidate->Add(*iter_line);
                }
            }

            //            candidate->SolveOnce();
            //            candidate->CombineFit();
            //            candidate->RemoveOutliers();

            for(iter_line=v_track->begin();iter_line!=v_track->end();){
                if(candidate->V>1 && candidate->Contains(*iter_line)){
                    v_track->erase(iter_line);
                } else iter_line++;
            }

            if(DEBUG){cout<<"Candidate: ";candidate->Print();}
            v_v2d->push_back(candidate);
            if(v_track->size()<2)break;
        }else {
            break;
        }
    }

    //=====Two-Track vertex
    for(int i=0;i<v_track->size();i++){
        for(int j=i+1;j<v_track->size();j++){
            Vertex2D *v=new Vertex2D;
            v->iside=iside;
            v->Add((v_track->at(i)));
            v->Add((v_track->at(j)));
            v->SolveOnce();
            double xx[2]={v->XY(),v->Z()};
            //  v->Print();
            v->V=fun_ZVTOP_V(xx,0);
            v_v2d->push_back(v);
        }
    }



    sort(v_v2d->begin(),v_v2d->end(),comp2DVertexByProbability);
    //   cout<<" v_2d candidate Sorted:  "<<endl;
    //  for(int i=0;i<v_v2d->size();i++)v_v2d->at(i).Print();



    int maxiter_v2d=10;
    int niter_v2d=0;
    int oldsize_v2d=v_v2d->size();
    vector<Vertex2D*>::iterator it_v2d;
    vector<Vertex2D*>::iterator it2_v2d;
    while(niter_v2d<maxiter_v2d){
        for(it_v2d=v_v2d->begin();it_v2d!=v_v2d->end();it_v2d++){
            if((*it_v2d)->V<=1)continue;
            for(it2_v2d=it_v2d+1;it2_v2d!=v_v2d->end();){
                if((*it_v2d)->Merge(*it2_v2d)){
                    double xx[2]={(*it_v2d)->XY(),(*it_v2d)->Z()};
                    (*it_v2d)->V=fun_ZVTOP_V(xx,0);
                    delete (*it2_v2d);
                    v_v2d->erase(it2_v2d);
                }
                else it2_v2d++;
            }
        }
        if(v_v2d->size()==oldsize_v2d){
            break;
        }
        oldsize_v2d=v_v2d->size();
        niter_v2d++;
    }

    //    vector<Vertex2D> v_v2d_candidate;
    //    double Vmin=999;
    //    Vertex2D candidate;
    //    for(int i=0;i<v_track->size();i++){
    //        for(int j=i+1;j<v_track->size();j++){
    //            Vertex2D v;
    //            v->iside=iside;
    //            v->Add((v_track->at(i)));
    //            v->Add((v_track->at(j)));
    //            v->SolveOnce();
    //            double xx[2]={v->XY(),v->Z()};
    //            if(fun_ZVTOP_V(xx,0)<Vmin){
    //                Vmin=fun_ZVTOP_V(xx,0);
    //                candidate=v;
    //                x=v->XY();
    //                z=v->Z();
    //            }
    //            v_v2d_candidate->push_back(v);
    //        }
    //    }

    //    vector<Vertex2D>::iterator it_v2d;
    //    for(it_v2d=v_v2d_candidate->begin();it_v2d!=v_v2d_candidate->end();){
    //        if(candidate->Merge(&(*it_v2d))){
    //            v_v2d_candidate->erase(it_v2d);
    //        }
    //        else it_v2d++;
    //    }

    //    candidate->Print();
    //    v_v2d->push_back(candidate);
    if(DEBUG){
        cout<<"FindBestVertex :"  <<Niteration<<"  finished,  N vertex: " <<v_v2d->size()<<endl;
        for(int i=0;i<v_v2d->size();i++)v_v2d->at(i)->Print();
    }
    return;

}


void VertexFinder2D_ZVTOP::FindVertex(){
    FindBestVertex();
}

double VertexFinder2D_ZVTOP::fun_ZVTOP_V(double* x, double* par)
{
    double sum=0;
    double sum2=0;
    double f;
    for(int i=0;i<v_track->size();i++){
        f=fun_ZVTOP_f(v_track->at(i),x);
        //        cout<<"F: " <<f<<endl;
        sum+=f;
        sum2+=f*f;
    }
    //        cout<<"par :  "<<x[0]<<",  "<<x[1]<<",  sum:"<<sum<<", sum2: " <<sum2<<",  V: "<<(sum-sum2/sum)<<endl;
    //    if(!sum)return 0;
    return sum-sum2/sum;
}

double VertexFinder2D_ZVTOP::fun_ZVTOP_f(Line* line, double *x)
{
    // use sigma=1 for test purpose, full treatment needs error propagation to around (x,z).
    static double sigma=1;
    double dist=line->GetDistance(x[0],x[1]);
    return exp(-0.5*pow(dist/sigma,2));
}


