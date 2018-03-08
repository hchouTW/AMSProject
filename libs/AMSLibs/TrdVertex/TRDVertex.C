#ifndef TRDVERTEX_C
#define TRDVERTEX_C
#include "TRDVertex.h"

TRDVertexFit::TRDVertexFit(){Clear();}
TRDVertexFit::~TRDVertexFit(){
  Clear();
}
void TRDVertexFit::Clear(){
  for(int i=0;i<_cxy.size();i++)   _cxy[i].clear(); 
  _cxy.clear();
  for(int i=0;i<_cz.size();i++)   _cz[i].clear(); 
  _cz.clear();
  for(int i=0;i<_sxy.size();i++)   _sxy[i].clear(); 
  _sxy.clear();
  for(int i=0;i<tinfo.size();i++)   tinfo[i].clear();
  tinfo.clear();
  _ntrack = 0;
  _ndof = 0;
  _chi2 = 0;
  _vtxcc[0] =0; _vtxcc[1] = 0;
}

int TRDVertexFit::AddTrack( vector<double> &trackx , vector<double> &trackz, vector<double> &errorx){
  if( trackx.size() != trackz.size()   ) return -1;
  if( trackx.size() != errorx.size()   ) return -1;
  _cxy.push_back(trackx);
  _cz.push_back(trackz );
  _sxy.push_back(errorx);
  _ntrack = _sxy.size();
  return 1;
}

int TRDVertexFit::FastIter(){
  if( _ntrack !=2   )  {Clear();  return -1;}
//  cout<<" ntrack is "<<_ntrack<<endl;  
  tinfo.clear();
  // get all results of all tracks  
  for(int itk =0;itk<_ntrack; itk++ )  {
    vector<double> params;//a, b, chisq, ndof, normchi2, itk 
    if(   StepFit(  _cxy[itk], _cz[itk],_sxy[itk], params ) <0) continue  ;
    //    if( params[2] > 100  ) continue;
    if( params[3] <= 0  ) continue;
    params.push_back( (double)itk );  
    tinfo.push_back(params);
  }

  double z0 = ( -tinfo[0][0] + tinfo[1][0] ) /   ( tinfo[0][1] - tinfo[1][1] )   ;
  double x0 = tinfo[0][0] + tinfo[0][1] * z0  ;
   _vtxcc[0] = x0;
   _vtxcc[1] = z0;
  return 1;  
}

double TRDVertexFit::FastFit(int nlim  ){
  
  if( _ntrack < nlim ) return -1;
  // get all results of all tracks  
  if( nlim == 2 ) {
    tinfo.clear();
    // get all results of all tracks  
    for(int itk =0;itk<_ntrack; itk++ )  {
      vector<double> params;//a, b, chisq, ndof, normchi2, itk 
      if(   StepFit( _cxy[itk], _cz[itk],_sxy[itk], params ) <0) continue  ;
      //    if( params[2] > 100  ) continue;
      if( params[3] <= 0  ) continue;
      params.push_back( (double)itk );  
      tinfo.push_back(params);
    }
  }
  
  tinfo.clear();
  for(int itk =0;itk<_ntrack; itk++ )  {
    vector<double> params;//a, b, chisq, ndof, normchi2, itk 
    if(   StepFit(  _cxy[itk], _cz[itk],_sxy[itk], params ) <0) continue  ;
    //    if( params[2] > 100  ) continue;
    if( params[3] <= 0  ) continue;
    params.push_back( (double)itk );  
    tinfo.push_back(params);
  }

  // get the rough estimation of vertex.  
  vector<double> a;
  vector<double> b;
  vector<double> s;
  for(int i=0;i<tinfo.size();i++  ) {
    a.push_back(tinfo[i][0] );
    b.push_back(-tinfo[i][1] );
    s.push_back(1);
  }
  vector<double> params;
  
  if(  StepFit( a,b,s,params ) < 0 ) return -1; // vertex fit failure;
  double vxy = params[0];
  double vz  = params[1];
  // recalculate the Chi2 with fix on vertex.
  _vtxcc[0] = vxy;
  _vtxcc[1] = vz;
  

  double reChi2 = 0;
  
  double atchi2  = 0;
  int    atnhit = 0;
  for(int icandi=0;icandi<tinfo.size();icandi++)  {
    int idx = (int) tinfo[icandi][5];
    //    if(idx < 0 || idx >=_cxy[idx].size() ) continue;
    double f1=0;
    double f2=0;
    for(int ih=0;ih<_cxy[idx].size();ih++ ) {
      f1 +=( ( vxy - _cxy[idx][ih]) *(vz - _cz[idx][ih]) ) ; 
      f2 += ( (vz - _cz[idx][ih])*(vz - _cz[idx][ih]) );
    }
    double tb = 1.0 * f1 / f2;
    double ta = vxy - tb  * vz;
    // now need to calculate residual of current track
    
    for(int ih=0;ih<_cxy[idx].size();ih++ ) {
      double res = _cxy[idx][ih] - ( ta + tb * _cz[idx][ih] )  ;
      atchi2 += res*res/_sxy[idx][ih]/_sxy[idx][ih]; 
      atnhit ++;
    } 
  }
  _chi2 = atchi2;
  _ndof = (atnhit - tinfo.size() - 2 );
  return  _chi2/_ndof;
}


double TRDVertexFit::StepFit(   vector<double> &trackx , vector<double> &trackz, vector<double> &errorx , vector<double> &param ){
    
  if( trackx.size() != trackz.size()   ) return -1;
  if( trackx.size() != errorx.size()   ) return -1;

  int nhit = trackx.size();
  if(nhit < 2) return -1;
  vector<double> res;
  res.resize(nhit);
  double mtx[2][2] = { { 0, 0 }, { 0, 0 } }, minv[2][2];
  double vec[2] = { 0, 0 };
  for (int i = 0; i < nhit; i++) {
    double w = (errorx[i] > 0) ? 1/errorx[i]/errorx[i] : 0;  
    mtx[0][0] += w;      mtx[0][1] += w*trackz[i];
    mtx[1][0] += w*trackz[i]; mtx[1][1] += w*trackz[i]*trackz[i];
    vec[0]    += w*trackx[i]; vec[1]    += w*trackx[i]*trackz[i];
  }
   double det = mtx[0][0]*mtx[1][1]-mtx[0][1]*mtx[1][0];
   if (det == 0) return -1;
   minv[0][0] = mtx[1][1]/det; minv[0][1] = -mtx[0][1]/det;
   minv[1][1] = mtx[0][0]/det; minv[1][0] = -mtx[1][0]/det;
   
   double a = minv[0][0]*vec[0]+minv[0][1]*vec[1];
   double b = minv[1][0]*vec[0]+minv[1][1]*vec[1];
//  for the params
//  x = a + b * z; 
   param.clear();
   param.push_back( a ) ;
   param.push_back( b ) ;
    
   double chisq = 0;
   for (int i = 0; i < nhit; i++) {
     res[i] = trackx[i]-(a+b*trackz[i]);
     chisq += (  errorx[i] > 0) ? res[i]*res[i]/errorx[i]/errorx[i] : 0;
   }
   param.push_back( chisq );
   double ndof = nhit-2;
   param.push_back(ndof);
   double normchi2 = 1.0* chisq/ndof;
   param.push_back( normchi2 );
   return normchi2;
}

double TRDVertexFit::GetChi2(){
  return _chi2;
}
double TRDVertexFit::GetNormChi2(){
  return _chi2/_ndof;
}


Double_t TRDVertex2DFit::func_Chi2Gain(Double_t *par){

  int ntrackx = fitx->tinfo.size();
  int ntracky = fity->tinfo.size();
  
  double atchi2  = 0;
  int    atnhit = 0;
  
  for(int iside=0;iside<2;iside++){
    TRDVertexFit *tfit = (iside ==0) ? fitx:fity ;
    for(int icandi=0;icandi<tfit->tinfo.size();icandi++)  {
      int idx = (int) (tfit->tinfo[icandi][5]);
      double z0  = par[2];
      double xy0 = (iside==0) ? par[0]:par[1]  ;
      double tb  = (iside==0) ? par[3+icandi] : par[3+ntrackx+icandi]   ;
      double ta =  xy0 - tb * z0 ;
//     cout<<__LINE__<<" , "<<idx<<" , "<<z0<<" , "<<xy0<<" , "<<ta<<" , "<<tb<< " , " << endl;
      for(int ih=0;ih<tfit->_cxy[idx].size();ih++ ) {
        double res = tfit->_cxy[idx][ih] - ( ta + tb * tfit->_cz[idx][ih] )  ;
//        cout<<"res "<<res<<endl;
        atchi2 += res*res/tfit->_sxy[idx][ih]/tfit->_sxy[idx][ih]; 
        atnhit ++;
      } 
    } 
  }
  return atchi2;
}

//void TRDVertex2DFit::fcn_TRDChi2( int &npar, double *gin, double &f, double *par, int iflag  ){
void TRDVertex2DFit::fcn_TRDChi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//  double chisq = 0;
  TRDVertex2DFit *tvtx = (TRDVertex2DFit *) gMinuit->GetObjectFit(); 
  f = tvtx->func_Chi2Gain( par  );
//  cout<<"step checking "<<endl; 
//  for(int i=0;i<npar;i++) cout<<Form("%3.3f ", par[i] );
//  cout<<f<<endl;
}

TMinuit *TRDVertex2DFit::gMinuit = 0;

double TRDVertex2DFit::DoFit(){
  if(fitx->FastFit(2) <0   ) return -1;
  if(fity->FastFit(2) <0   ) return -1;
  if(fabs( fitx->_vtxcc[1] - fity->_vtxcc[1] ) > 30. ) return -1;
  double x0 =   fitx->_vtxcc[0];
  double y0 =   fity->_vtxcc[0];
  double z0 = ( fitx->_vtxcc[1] + fity->_vtxcc[1] ) / 2. ;
  
  int ntrackx = fitx->tinfo.size();
  int ntracky = fity->tinfo.size();
  int npar = ntrackx + ntracky + 3;
  double *vstart = new double[npar];
  vstart[0] = x0;
  vstart[1] = y0;
  vstart[2] = z0;
  for(int i=0; i< ntrackx;i++ )    vstart[i+3] = fitx->tinfo[i][1];
  for(int i=0; i< ntracky;i++ )    vstart[i+3+ntrackx] = fity->tinfo[i][1];
  
//  cout<<"printing the as"<<endl;
//  for(int i=0; i< ntrackx;i++ )    {
//    for(int j=0;j<6;j++)
//      cout<<fitx->tinfo[i][j]<<" , ";
//    cout<<endl;
//  }
//  for(int i=0; i< ntracky;i++ )    {
//    for(int j=0;j<6;j++)
//      cout<<fity->tinfo[i][j]<<" , ";
//    cout<<endl;
//  }
     
  double *vstep = new double[npar];
  vstep[0] = 0.1;
  vstep[1] = 0.1;
  vstep[2] = 0.9;
  for(int i=0; i< ntrackx;i++ )    vstep[i+3] = 0.1;
  for(int i=0; i< ntracky;i++ )  vstep[i+3+ntrackx] = 0.1;

  if( !gMinuit  )  
    gMinuit = new TMinuit(50);
  gMinuit->mncler();
  gMinuit->SetFCN(fcn_TRDChi2);
  gMinuit->SetObjectFit(this); 
  gMinuit->SetPrintLevel(-1);
  double arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  for(int i=0;i<3;i++)
    gMinuit->mnparm(i, Form("c_%i",i+1), vstart[i], vstep[i], 0, 0,ierflg);
  
  for(int i=3;i<3+fitx->_ntrack;i++)
    gMinuit->mnparm(i, Form("trkx_%i",i+1), vstart[i], vstep[i], -10,10,ierflg);

  for(int i=3+fitx->_ntrack;i<3+fitx->_ntrack+fity->_ntrack;i++)
    gMinuit->mnparm(i, Form("trky_%i",i+1), vstart[i], vstep[i], -10,10,ierflg);


  delete vstep;
  delete vstart;

//  cout<<"checking pars "<<endl;
//  for(int i=0;i<npar;i++)
//    cout<<Form( "%3.3f %3.3f",vstart[i],vstep[i]   )<<endl;
  for(int i=npar;i<50;i++  )
    gMinuit->FixParameter(i); 
//fixparameter
  gMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);
  double cc[3], err[3];
  TString s_name[3]={"c_1","c_2","c_3"};
  for(int i=0;i<3;i++) {
    double bnd1, bnd2;
    int ivar;
    gMinuit->mnpout(i, s_name[i], cc[i], err[i], bnd1, bnd2, ivar);
  }
//  for(int i=0;i<3;i++) {
//    cout<<   cc[i]  <<"  ,  "  <<  err[i]<<" :  "<<endl  ;
//  }
//mnpout
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
//  cout<<Form("%3.3f, %3.3f, %3.3f, %i, %i, %i",  amin,edm,errdef,nvpar,nparx,icstat )<<endl;
//  delete gMinuit;
  _chi2 = amin;
  _ndof = fitx->_ndof + fity->_ndof -1;
  _normchi2 = 1.0*_chi2/_ndof;
  _ccvtx.setp(cc[0],cc[1],cc[2] );
  _ntrack = fitx->tinfo.size() +  fity->tinfo.size();
  return  _normchi2; 
}



#endif 
