#ifndef TRDNEWREC_C
#define TRDNEWREC_C
#include "TRDNewRec.h"
#include "TRDVertex.h"
#include "TRDVertex.C"



int TRDNewRec::_trtrkh = 0;
double TRDNewRec::_maxcos  = -0.99;
double TRDNewRec::_maxchi2 = 20;

TRDNewRec::TRDNewRec(){
  hitco = 0;
  _evt = 0;
  Clear();
}
TRDNewRec::~TRDNewRec(){
  Clear();
}

int TRDNewRec::GetNTrack(int iside){
  if( iside == 1 ) return xcaTrack.size();
  if( iside == 2 ) return ycaTrack.size();
  return 0;
}

//int TRDNewRec::AddVertex( AMSPoint cc, int iside = 1 ){
//  if( iside !=1 && iside != 2  ) return -1;
//  _vcc[iside-1].push_back( cc );
//  return _vcc[iside-1].size();
//}

void TRDNewRec::Init( vector<SimpleHitState> *hits, int itree ){
  hitco = hits;
  _itree = itree;
}

void TRDNewRec::Init( AMSEventR *evt, int itree ){
  _evt   = evt;
  _itree = itree;
}



void TRDNewRec::Clear(){
  if( hitco &&  hitco->size() > 0)  hitco->clear();
  _evt = 0;
  hitco =0;
  _ntrack = 0;
  _nvertex = 0;
  xcaTrack.clear();
  ycaTrack.clear();
  _vcc[0].setp(0,0,0);
  _vcc[1].setp(0,0,0);
  _vdd[0] = 0;
  _vdd[1] = 0;
}

int TRDNewRec::Process(){
  TRDNewRec::_maxcos = -0.99;
  BuildTrk(1);  //x track
  TRDNewRec::_maxcos = -0.994;
  BuildTrk(2); // y track
  BuildVertex( 0 );
  return 1;  
}

int TRDNewRec::BuildVertex( int imethod  ){
  if( imethod == 1  ) {
    if( xcaTrack.size() >= 2 )  TRDNewRec::BuildVtx(1); 
    if( ycaTrack.size() >= 2 )  TRDNewRec::BuildVtx(2); 
    ReFitVtx(); 
  }
  else {
    BuildIterVtx(); 
  }

}


struct trkcomb{
  int    idx[2];
  double cxy;
  double cz;
  void  setp( int id1, int id2, float tcxy, float tcz ) {
    idx[0] = id1;
    idx[1] = id2;
    cxy = tcxy;
    cz  = tcz ;
  }
};

struct  cgp{
  vector<double> cent;
  vector<vector<double> > cz;
  vector<vector<int> >    iz;
  void AddHit( double tc, int tidx ) {
    if(  cent.size() >0  ){
      bool isfd = false;
      for( int i=0;i<cent.size();i++  ) {
        if( fabs(  tc - cent[i] ) < 2 ) {
          cz[i].push_back( tc  );
          iz[i].push_back( tidx);
          double ctot=0;
          for( int ic=0;ic<cz[i].size();ic++  ) 
            ctot += cz[i][ic];
          ctot  /= cz[i].size();
          cent[i] = ctot;
          isfd = true;
          break;
        }
      }
      if( !isfd  ) {
        vector<double> tmpv;
        tmpv.push_back( tc );
        cz.push_back(tmpv);
        
        vector<int>   tmidx;
        tmidx.push_back(tidx);
        iz.push_back( tmidx );
        
        cent.push_back(tc);
      }
    }
    else {
     vector<double> tmpv;
     tmpv.push_back( tc );
     cz.push_back(tmpv);
     vector<int>   tmidx;
     tmidx.push_back(tidx);
     iz.push_back( tmidx );
     cent.push_back( tc );
    }
  }
  int NGroup() {
    return cent.size();
  }
  int GetGroup( int ig, vector<double> &cg, vector<int> &tidxg ) {
    if( ig>= NGroup() || ig<0  ) return -1;
    cg = cz[ig];
    tidxg = iz[ig];
    return cg.size();
  }
};


int TRDNewRec::BuildIterVtx(){
  
  if(  xcaTrack.size() <2 || ycaTrack.size()<2   ) return -1;
  TRDVertexFit *tvtx = new TRDVertexFit();
  vector<trkcomb>  icomb[2];
//  cout<<" track size is "<<xcaTrack.size()<<" , "<<ycaTrack.size()<<endl;
  for( int iside=1;iside<=2;iside++  ) {
    vector<semiIndex> *caTrack = 0;
    if( iside == 1 )  caTrack = &xcaTrack;
    if( iside == 2 )  caTrack = &ycaTrack;
    for(int ic1=0;ic1 < caTrack->size()-1;ic1++ ) {
      for(int ic2 = ic1+1;ic2 < caTrack->size();ic2++ ) {
          if(  ic1 ==  ic2  ) continue;
//quick way to get the combination          
          tvtx->Clear();
          for( int ic=0;ic<2;ic++  ){
            int idk = (ic==0) ? ic1:ic2;
            semiIndex tcandi = caTrack->at(idk);
            vector<double> errorx;
            vector<double> trackx;
            vector<double> trackz;
            for(int i=0;i<tcandi.nLength;i++ ) {
              errorx.push_back( 0.6 );
              trackx.push_back( (double) tcandi.cxy[i] );
              trackz.push_back( (double) tcandi.cz[i]  );
            }
            tvtx->AddTrack( trackx, trackz, errorx);
          }
          if(tvtx->FastIter() <=0  ) continue  ;
          trkcomb tcomb;
          tcomb.setp( ic1,ic2,  tvtx-> _vtxcc[0], tvtx-> _vtxcc[1]   );
          icomb[iside-1].push_back( tcomb  );
      }
    } 
  }
  vector<double>  czg;
  vector<int>     izg;
  for(int iside=1;iside<=2;iside++) {
    for(int itk=0;itk<icomb[iside-1].size();itk++){
      czg.push_back( icomb[iside-1][itk].cz);
      izg.push_back( (iside-1)*100 + itk  );
    }
  }
   
  vector<int> reID;
  
  vector<int> igroup;
  cgp  hitcont;
  for( int i=0;i< czg.size();i++  )
    hitcont.AddHit( czg[i] ,i  );

//  cout<<"geting groups "<<hitcont.NGroup()<<endl;
  for( int ig=0;ig<hitcont.NGroup();ig++  ) {
    vector<double> cg;
    vector<int>    tigz;
    int ires = hitcont.GetGroup( ig, cg , tigz);
    if( ires <= 1  ) continue;
    int ngood[3] = {0,0,0};
    for( int i=0;i<tigz.size();i++  ){
      ngood[0] ++;
      if( izg[tigz[i]] <100   ) ngood[1]++;
      else    ngood[2] ++;
    }
//  now we can start the fitting procedure.
    if( ngood[1]*ngood[2] <= 0 )  continue;
    for( int igood=1;igood<=2;igood++ ) {
      if( ngood[igood] <=1 ) {
        if(  ngood[igood] <=0  ) continue;
        for( int i=0;i<tigz.size();i++  ){ 
          if( igood == 1 &&  izg[tigz[i]] < 100 ) 
            reID.push_back( izg[tigz[i]]  );
          if( igood == 2 &&  izg[tigz[i]] >= 100 ) 
            reID.push_back( izg[tigz[i]]  );
        }
        continue;
      }
      vector<double>  cyg;
      vector<int>     iyg;
      // else, we should consider the xy direction compitable
      for( int i=0;i<tigz.size();i++  ){
          int iside = 0;
          int tidx = izg[tigz[i]] ;
          if( tidx<100 ) iside = 1;
          else iside =2;
          if( iside != igood  ) continue;
          int itk = (tidx % 100 ) ;
          cyg.push_back( icomb[iside-1][itk].cxy);
          iyg.push_back( (iside-1)*100 + itk  );
      }
      
     
      vector<int> igroupy;
      cgp  hitconty;
      for( int i=0;i< cyg.size();i++  )
        hitconty.AddHit( cyg[i] ,i  );
     
      int maxN = 0, imaxN = -1;
      for( int igy=0;igy<hitconty.NGroup();igy++ ) {
        vector<double> cgy;
        vector<int>    tigy;
        int iresy = hitconty.GetGroup( igy, cgy , tigy);
        if(  iresy >maxN  ) {
          maxN = iresy;
          imaxN = igy;
        }
      }// igy  NGroup
     
      if( imaxN >=0  ){
        vector<double> cgy;
        vector<int>    tigy;
        int iresy = hitconty.GetGroup( imaxN, cgy , tigy);
        for(  int i=0;i< tigy.size();i++  )  {
          reID.push_back(  iyg[tigy[i]] );
        }
      }
    } // iside, igood
  }

//  for(int i=0;i<reID.size();i++) 
//    cout<<reID[i]<<" , ";
//   cout<<endl;
 //reID stores all the inform of icomb.index;
   chntrkidx[0].clear(); 
   chntrkidx[1].clear(); 
  for(int ire=0;ire<reID.size();ire++  ) {
    int tidx = reID[ire];
    int index = tidx %100;
//    cout<<"index is "<< index<<endl;
    if(tidx <100 )  {
      bool isfd[2] = {false,false};
      for( int i=0;i<chntrkidx[0].size();i++   ){
        if(  icomb[0][index].idx[0] == chntrkidx[0][i]  )    isfd[0] = true; 
        if(  icomb[0][index].idx[1] == chntrkidx[0][i]  )    isfd[1] = true; 
      }
      if(!isfd[0] )  chntrkidx[0].push_back( icomb[0][index].idx[0] ) ;   
      if(!isfd[1] )  chntrkidx[0].push_back( icomb[0][index].idx[1] ) ;   
    }
    if( tidx >= 100  ) {
      bool isfd[2] = {false,false};
      for( int i=0;i<chntrkidx[1].size();i++   ){
        if(  icomb[1][index].idx[0] == chntrkidx[1][i]  )    isfd[0] = true; 
        if(  icomb[1][index].idx[1] == chntrkidx[1][i]  )    isfd[1] = true; 
      }
      if(!isfd[0] )  chntrkidx[1].push_back( icomb[1][index].idx[0] ) ;   
      if(!isfd[1] )  chntrkidx[1].push_back( icomb[1][index].idx[1] ) ;   
    } 
  }
  
//  cout<<"checking the final results "<<endl;
//  for(int i=0;i<2;i++)  {
//    for(int j=0;j<  chntrkidx[i].size();j++ ) {
//      cout<< chntrkidx[i][j]<<" , ";
//    }
//    cout<<endl;
//  }
  ReFitVtx();  
//  cout<<
    
  
  return 0; 
}  


int  TRDNewRec::ReGroup( vector<double> &czg ){
//  vector<int> igroup;
//  cgp  hitcont;
//  for( int i=0;i< czg.size();i++  )
//    hitcont.AddHit( czg[i] ,i );
//  cout<<"geting groups "<<hitcont.NGroup()<<endl;
//
//  for( int ig=0;ig<hitcont.NGroup();ig++  ) {
//    vector<double> cg;
//    vector<int>    tigz;
//    int ires = hitcont.GetGroup( ig, cg, tigz );
//    if( ires <= 1  ) continue;
//     for(int i=0;i<cg.size();i++) cout<<cg[i]<<" , ";
//     cout<<endl;
//  }

  return 0;
}

int  TRDNewRec::ReGroup( vector<double> &cxg, vector<double> &czg ) {
  
  
  return 0;
}












double TRDNewRec::ReFitVtx(){
  
//  if( _vcc[0].z() ==0 || _vcc[1].z() == 0   ) return -1;
//  if(  fabs(  _vcc[0].z() -  _vcc[1].z() ) >35  ) return -1;
//  cout<<"##### for x side tracks "<< chntrkidx[0].size()<<endl;
//  cout<<"##### for y side tracks "<< chntrkidx[1].size()<<endl;
  TRDVertexFit *tvtx = new TRDVertexFit();
  TRDVertexFit *tvty = new TRDVertexFit();
  for(int ic=0;ic < xcaTrack.size();ic++ ) {
    bool isexist = false;
    for(int i=0;i<chntrkidx[0].size();i++  )
      if( ic == chntrkidx[0][i] ){
        isexist = true;
        break;
      }
    if(!isexist) continue;
    semiIndex tcandi = xcaTrack[ic];
    vector<double> errorx;
    vector<double> trackx;
    vector<double> trackz;
    for(int i=0;i<tcandi.nLength;i++ ) {
      errorx.push_back( 0.6 );
      trackx.push_back( (double) tcandi.cxy[i] );
      trackz.push_back( (double) tcandi.cz[i]  );
    }
    tvtx->AddTrack( trackx, trackz, errorx);
  }
  for(int ic=0;ic < ycaTrack.size();ic++ ) {
    bool isexist = false;
    for(int i=0;i<chntrkidx[1].size();i++  )
      if( ic == i ){
        isexist = true;
        break;
      }
    if(!isexist) continue;
    semiIndex tcandi = ycaTrack[ic];
    vector<double> errorx;
    vector<double> trackx;
    vector<double> trackz;
    for(int i=0;i<tcandi.nLength;i++ ) {
      errorx.push_back( 0.6 );
      trackx.push_back( (double) tcandi.cxy[i] );
      trackz.push_back( (double) tcandi.cz[i]  );
    }
    tvty->AddTrack( trackx, trackz, errorx);
  }
   
  tvtx->FastFit(2);
//  cout<<__LINE__<<"   fitting results "<<tvtx->_vtxcc[0]<<" , "<<tvtx->_vtxcc[1]<<" , "<<tvtx->_chi2<<" , "<<tvtx->_ndof<<endl;
tvty->FastFit();
//  cout<<__LINE__<<"   fitting results "<<tvty->_vtxcc[0]<<" , "<<tvty->_vtxcc[1]<<" , "<<tvty->_chi2<<" , "<<tvty->_ndof<<endl;

  TRDVertex2DFit *fit2d = new TRDVertex2DFit();
  fit2d->SetV( tvtx, tvty );
  if(fit2d->DoFit() >0)  {
    fit2d->GetVertex( _vcc[0] ); 
    fit2d->GetVertex( _vcc[1] ); 
    _vdd[0] = fit2d->GetNormChi2();
    _vdd[1] = fit2d->GetNormChi2();
  }
//  cout<<_vcc[0].x()<<" , "<<_vcc[0].y()<< " , " <<_vcc[0].z()<<endl; 
//  cout<<_vcc[1].x()<<" , "<<_vcc[1].y()<< " , " <<_vcc[1].z()<<endl; 

  delete tvtx; 
  delete tvty;
  return 0;
}


int TRDNewRec::FillHits(   vector<ccL> &ccList, int iside   ){
  if( _evt  ) {
//    for(int ihit=0;ihit<_evt->   ) 
    
    return -1;  
  }
  else if(hitco) {
    
    for(int ihit=0;ihit<hitco->size(); ihit++) {
      SimpleHitState &ssht = hitco->at(ihit);
      //    if(  ssht.layer > -16  )  continue;
      double cxy = -9999;
      if( iside == 1 && ssht.x != 0  ) cxy = ssht.x; 
      else if( iside == 2 && ssht.y != 0  ) cxy = ssht.y; 
      else continue;

      ccL tmpccL;
      bool isfd = false;
      for(int ic=0;ic<ccList.size();ic++){
        if( ccList[ic]. addAj(ihit, ssht.layer,  cxy , ssht.z , ssht.Amp, ssht.AmpY) >0       ){
          isfd = true;
          break;
        }
      }
      if( !isfd  ){
        tmpccL.addCC(ihit, ssht.layer,  cxy , ssht.z , ssht.Amp, ssht.AmpY );
        ccList.push_back( tmpccL);
      }
    }

    
  }
  else return -1;
  return -1;
}



int TRDNewRec::BuildTrk( int iside  ){
  //  if(iside==1)   cout<<"building xz trdtrack"<<endl;
  //  if(iside==2)   cout<<"building yz trdtrack"<<endl;

  vector<ccL> ccList;	
  FillHits( ccList, iside );
//  cout<<"iside is "<<iside<<endl;
//  for(int i=0;i<ccList.size();i++)  {
//    ccL tc = ccList[i];
//    cout<< tc.nLayer<<" , "<<tc.cxy<<" , "<<tc.cz<<endl;
//  }
//  for(int ihit=0;ihit<hitco->size(); ihit++) {
//    SimpleHitState &ssht = hitco->at(ihit);
//    //    if(  ssht.layer > -16  )  continue;
//    float cxy = -9999;
//    if( iside == 1 && ssht.x != 0  ) cxy = ssht.x; 
//    else if( iside == 2 && ssht.y != 0  ) cxy = ssht.y; 
//    else continue;
//
//    ccL tmpccL;
//    bool isfd = false;
//    for(int ic=0;ic<ccList.size();ic++){
//      if( ccList[ic]. addAj(ihit, ssht.layer,  cxy , ssht.z , ssht.Amp, ssht.AmpY) >0       ){
//        isfd = true;
//        break;
//      }
//    }
//    if( !isfd  ){
//      tmpccL.addCC(ihit, ssht.layer,  cxy , ssht.z , ssht.Amp, ssht.AmpY );
//      ccList.push_back( tmpccL);
//    }
//  }
  
  
  //  if( ccList.size() > 40 ) {
  ////    iend = lmax-1;
  //    double weightR = 1. + 300.*190./pow( 380./(1.0*ccList.size() - 40 ), 3 );
  //    double kiklimt = 0.01 /weightR ;
  //    TRDNewRec::_maxcos = kiklimt - 1;
  //    cout<<"size "<<ccList.size()<<"   "<<TRDNewRec::_maxcos<<endl;
  //  }
  //  cout<<"total hits "<<ccList.size()<<endl;
  //  for(int ihit=0;ihit<ccList.size();ihit++) {
  //    cout<<Form("layer %i, %3.1f, %3.1f, %3.1f, %s, %i",ccList[ihit].nLayer, ccList[ihit].cxy, ccList[ihit].cz, ccList[ihit].amp,(  ccList[ihit].isWide?"wide core":"narrow core"  ) , ccList[ihit].nadj  )<<endl;
  //    if( ccList[ihit].nadj > 1 )
  //    {
  //      for(int i=0;i<ccList[ihit].nadj;i++){
  //        cout<<"   #### "<<ccList[ihit].ajcxy[i] <<"   , "<<ccList[ihit].ajamp[i]   <<endl;
  //        
  //      }
  //    }
  //  }
  int ncl = ccList.size();
  ccL tmp;
  for(int i=0;i<ncl-1;i++){
    for(int j=0;j<ncl-i-1;j++){
      if(ccList[j].cz < ccList[j+1].cz){
        tmp = ccList[j];
        ccList[j] = ccList[j+1];	
        ccList[j+1] = tmp;
      }
    }
  }
  int nNeuron=0; 
  vector<vector<int>  > sNeu;
  sNeu.resize(ncl);
  for(int i=0;i<ncl;i++)
    sNeu[i].resize(ncl);
  for(int i=0;i<ncl;i++) {
    for(int j=0;j<ncl;j++){
      sNeu[i][j]=-9999;
      if(ccList[i].nLayer < ccList[j].nLayer)      {
        sNeu[i][j]=1;
        nNeuron++;
      }
    } 
  }

  ///////////////////////////////////////////////////////////
  /////////  start the CATS iteration process   ////////////
  int niter = (iside == 1) ?11: 7;

  vcos headcos;
  for(int isel=0;isel<niter;isel++){
    for(int i=0;i<ncl-2;i++) {
      for(int j =i+1;j<ncl-1;j++) {
        if(ccList[i].nLayer >= ccList[j].nLayer) continue;
        //        if( ccList[i].isWide ) continue;
        int flgS=1;
        int lmax=1;
        for(int k=j+1;k<ncl;k++)    {
          if(ccList[j].nLayer >= ccList[k].nLayer) continue;
          if(sNeu[i][j]!=sNeu[j][k]) continue;
          if(headcos.getcosyz(ccList[i],ccList[j],ccList[k])> TRDNewRec::_maxcos ) continue;
          //          cout<<headcos.getcosyz(ccList[i],ccList[j],ccList[k])<<endl;
          lmax=sNeu[j][k];
          flgS=0;
        }
        if(flgS==0){
          sNeu[i][j]=lmax+1;
        } }
    }
  }// end of iteration

  //  return 0;

  vector<vector<int> > tmin;
  tmin.resize( nNeuron );
  for(int i=0;i<nNeuron;i++ )
    tmin[i].resize( 3  );
  int ntmin=0;
  int nS[11] = {0,0,0,0,0,0, 0,0,0,0,0};
  for (int i=0;i<nNeuron;i++) {
    for(int j=i+1;j<3;j++) {
      tmin[i][j] = -9999;
    }}
  for (int i=0;i<ncl-1;i++) {
    for(int j=i+1;j<ncl;j++) {
      if(sNeu[i][j]==-9999) continue;
      //      cout<<sNeu[i][j]<<endl;
      tmin[ntmin][0] = sNeu[i][j];
      tmin[ntmin][1] = i;
      tmin[ntmin][2] = j;
      nS[niter-sNeu[i][j]]++;
      ntmin++;
    }
  }
  // 
  int tmp_tmin[3];
  for(int i=0;i<ntmin-1;i++)
  {
    for(int j=0;j<ntmin-i-1;j++)
    {
      if(tmin[j][0]  <  tmin[j+1][0]){
        tmp_tmin[0] = tmin[j][0];
        tmp_tmin[1] = tmin[j][1];
        tmp_tmin[2] = tmin[j][2];

        tmin[j][0] =tmin[j+1][0];
        tmin[j][1] =tmin[j+1][1];
        tmin[j][2] =tmin[j+1][2];

        tmin[j+1][0] = tmp_tmin[0];	
        tmin[j+1][1] = tmp_tmin[1];	
        tmin[j+1][2] = tmp_tmin[2];	
      }
    }
  }
  if(tmin[0][0]<2)
  {
    ccList.clear();
    return 0;
  }
  //  for(int i=0;i<ntmin;i++)
  //    cout<<i<<"th "<<tmin[i][0]  <<" , "<<tmin[i][1]<<" , "<<tmin[i][2]<<endl;
  //  cout<<endl;
  vector<int> stmin;
  stmin.resize( tmin.size());
  for(int i=0;i<stmin.size();i++) stmin[i] = 0;
  //  return 0;

  int  lmax =  tmin[0][0];
  CombIter tcomb( tmin, ccList , _maxcos );
  //  tcomb.selfiter( lmax, lmax );
  int iend = 3;
  //  if( ccList.size() >=70  ) iend = lmax-1;
  for(int ik = lmax; ik >= iend; ik--){
    tcomb.selfiter( ik, ik );
  }

  for(int itk=0;itk<tcomb.candi.size();itk++) {
    tcomb.candi[itk].isBad = false;
  }  
  tcomb.CheckOverlap();

  for(int itk=0;itk<tcomb.candi.size();itk++) {
    semiIndex tcandi = tcomb.candi[itk]; 
    if(tcandi.isBad ) continue;

  //  semiIndex &trk = tcandi;
  //  for( int il=0;il<trk.nLength;il++) {
  //    for(int in=0;in<trk.nLength; in++)  {
  //      int tid  = trk.trkIndex[in];
  //      cout<<trk.trkIndex[in]<<" , "<<trk.levelIndex[in]<<" , "<<trk.cxy[in]<<" , "<<trk.cz[in]<<" , "<< ccList[tid].nLayer<<" , "<< ccList[tid].cxy<<" , "<<ccList[tid].cz<<     endl;   
  //    }
  //  }
//    tcandi.Stat(); 
    for(int i=0;i<tcandi.nLength; i++) {
      int tid = tcandi.trkIndex[i];
      ccList[tid].isused = true;
    }
    if( iside == 1   )  xcaTrack.push_back( tcandi   ); 
    if( iside == 2   )  ycaTrack.push_back( tcandi   );
  }
  if(iside ==1 ) xccList = ccList;
  if(iside ==2 ) yccList = ccList;
  
  if( iside == 1  )  {
    for(int i=0;i<xcaTrack.size();i++)
    RefineTrk( xcaTrack[i],iside ) ;
  } 
  if (iside == 2  )  {
    for(int i=0;i<ycaTrack.size();i++)
    RefineTrk( ycaTrack[i],iside ) ;
  }
  return 0; 

}
int TRDNewRec::RefineTrk( semiIndex &trk , int iside  ){
  vector<ccL> &ccList = (iside == 1) ? xccList: yccList;
  // step 1.
  // solve the lost hit problem;
  
  
  for( int ilk=0;ilk>-20;ilk--   ) {
    if(iside ==1)   if( ilk>-4 && ilk<-15  ) continue;
    else if( iside ==2  )  if( ilk>0 || (ilk >-3&& ilk >-16) || ilk <-19 ) continue;
    else break;
    bool isexist = false;
    for(int ic=0;ic<trk.nLength; ic++){
      if( ilk == trk.levelIndex[ic] ) {isexist =true; break;}
    }
    if( isexist  ) continue;
    int imin = -1; double qmin = 200;

    for( int ic=0;ic<ccList.size();ic++  ) {
      if(  ccList[ic].isused  )  continue;
      if(  ccList[ic].nLayer != ilk )   continue;
      double qua = trk.AFit( ccList[ic].cxy, ccList[ic].cz  );
      if( qua <= 0 ) continue;
      if( qua > 10  ) continue;
      if( qua <qmin  ) {
        qmin = qua;
        imin = ic;
      }
    }
    if( imin>= 0 && qmin <10 ) {
//      cout<<"adding "<< ccList[imin].nLayer  <<" , "<< ccList[imin].cz<<endl;
      trk.AddHit(  imin, ccList[imin].nLayer,  ccList[imin].cxy, ccList[imin].cz,1 );
    }
  }
  trk.Quality();

// step 2.
  // solve the merging problem;
  for( int il=0;il<trk.nLength;il++) {
 
    int tid = trk.trkIndex[il];
    ccL tcandi = ccList[tid];
//    cout<<" nadj is "<<tcandi.nadj<<endl;
    if( tcandi.nadj ==0  ) continue;
    vector<float> ql;
    ql.resize( tcandi.nadj );
    for( int ic=0;ic<tcandi.nadj; ic++ ) {
      // refit the current track, to find the best position
      ql[ic] = trk.LFit( il, tcandi.ajcxy[ic] );   
    }
    float mql= 9999;
    int   mic = -1;
    for(int i=0;i<ql.size();i++) {
      if(ql[i]<=0) continue;
      if( ql[i] <  mql  ){
        mql = ql[i];
        mic = i;
      }
    }
    if( mic >= 0 && mql !=0 && mql <9999 ) {
//      trk.SetHit( tcandi.ajcxy[mic], il);
//      cout<<Form("int %i th, orignal %3.3f, new %3.3f ", il+1, trk.cxy[il]  , tcandi.ajcxy[mic]  )<<endl;
    }
  }


}
int TRDNewRec::BuildVtx( int iside  ){


  vector<vector<double> > aterrorx;
  vector<vector<double> > attrackx;
  vector<vector<double> > attrackz;

  vector<semiIndex> caTrack ;
  if(iside ==1 ) {
    caTrack = xcaTrack;
  }
  else if(iside ==2 ) {
    caTrack = ycaTrack;
  }
  else return 0;

  if( caTrack.size() >2 ) {
    int iend = caTrack.size();
    TRDVertexFit *tvtx = new TRDVertexFit();
    vector<vector<double> > aterrorx;
    vector<vector<double> > attrackx;
    vector<vector<double> > attrackz;
    vector<int>  trkindex;
    for(int ic=0;ic < caTrack.size();ic++ ) {
      semiIndex tcandi = caTrack[ic];
      vector<double> errorx;
      vector<double> trackx;
      vector<double> trackz;
      for(int i=0;i<tcandi.nLength;i++ ) {
        errorx.push_back( 0.6 );
        trackx.push_back( (double) tcandi.cxy[i] );
        trackz.push_back( (double) tcandi.cz[i]  );
      }
      aterrorx.push_back(errorx);
      attrackx.push_back(trackx);
      attrackz.push_back(trackz);
    }
    tvtx->Clear();
    trkindex.clear();
    double tchi2 = -1;

    for( int ic=0; ic<aterrorx.size();ic++ ) {
      tvtx->AddTrack( attrackx[ic], attrackz[ic], aterrorx[ic] );
      trkindex.push_back( ic );
    }
    if(  caTrack.size() == 2 ) {
      if( tvtx->FastFit(2) > 0 ) {
        tchi2 = tvtx->_chi2;
        if( tchi2 >0 && tchi2 <TRDNewRec::_maxchi2  ){
          // we got the correct combination
          SaveVertex(  tvtx, iside, trkindex  ) ;             
        } 
        return tchi2;
      }
      return 0;
    }  // ntrack ==2 situation
    if(  tvtx->FastFit() > 0  )   { 
      tchi2 = tvtx->_chi2;    
    }
    
    if( tchi2 >0 && tchi2 <TRDNewRec::_maxchi2  ){
      // we got the correct combination
      SaveVertex(  tvtx, iside, trkindex  ) ;             
    } 
    else {
      // we start to find the combinatin procedure
      vector<int> godidx;
      vector<int> badidx;
      while (true){
        if(  caTrack.size() - badidx.size() < 3 ) break;
        vector<double> qvalue;
        for( int ik=0;ik< aterrorx.size();ik++  ) qvalue.push_back(-1);
        for(int ik=0;ik< aterrorx.size();ik++ ){
          bool iskp0 = false;
          for(int j0=0;j0<badidx.size();j0++)  {
            if(  ik == badidx[j0]  ) {iskp0 = true; break;}
          }
          if(iskp0) continue;
//          cout<<__LINE__<<"   fitting results "<<tvtx->_vtxcc[0]<<" , "<<tvtx->_vtxcc[1]<<" , "<<tvtx->_chi2<<" , "<<tvtx->_ndof<<" , "<<tvtx->_ntrack<<endl;
          tvtx->Clear();
          trkindex.clear();
          for(int ic=0; ic<aterrorx.size();ic++) {
            if(ik == ic ) continue;
            bool iskp = false;
            for(int j=0;j<badidx.size();j++)  {
              if(  ic== badidx[j]  ) {iskp = true; break;}
            }
            if(iskp) continue;
            tvtx->AddTrack( attrackx[ic], attrackz[ic], aterrorx[ic] );
            
          }//adding ic
          if(  tvtx->FastFit() > 0  )   { 
            tchi2     = tvtx->_chi2;    
            qvalue[ik] = tvtx->_chi2;
          }
          else {
            qvalue[ik] = -1;
          }
//          cout<<__LINE__<<"   fitting results "<<tvtx->_vtxcc[0]<<" , "<<tvtx->_vtxcc[1]<<" , "<<tvtx->_chi2<<" , "<<tvtx->_ndof<<" , "<<tvtx->_ntrack<<endl;
        }//ik
        double minqv = 9e6;
        int   iminqv = -1;
        for(int ick=0;ick<qvalue.size(); ick++ ) {
          if( qvalue[ick] <= 0 ) continue;
          if( minqv > qvalue[ick]  )  { minqv = qvalue[ick]; iminqv = ick;}
        }
        if( minqv < TRDNewRec::_maxchi2 ) {
          // should break here, with track index,
          trkindex.clear();
          tvtx->Clear(); 
          for(int ic=0;ic< aterrorx.size();ic++)  {
            if( ic == iminqv   ) continue;
            bool iskp = false;
            for(int j=0;j<badidx.size();j++)  {
              if(  ic== badidx[j]  ) {iskp = true; break;}
            }
            if(iskp) continue;
            tvtx->AddTrack( attrackx[ic], attrackz[ic], aterrorx[ic] ); 
            trkindex.push_back( ic  );
          }
          tvtx->FastFit() ;             
          SaveVertex(  tvtx, iside, trkindex ) ;             
          break;
        }
        else if( minqv >= TRDNewRec::_maxchi2 ) {
          badidx.push_back( iminqv ); 
          continue;  
        }
      }
    }
    delete tvtx;
  }   // itrack >= 3



  return 0;
}




void  TRDNewRec::SaveVertex( TRDVertexFit *tvtx , int iside, vector<int> &vtkidx )
{
  if(iside == 1){
    _vcc[0].setp( tvtx->_vtxcc[0] ,0, tvtx->_vtxcc[1] ); 
    if(tvtx->_ndof !=0 ) _vdd[0] = 1.0*tvtx->_chi2/tvtx->_ndof ;
    _ntxtrack[0] = tvtx->_ntrack;
    chntrkidx[0] = vtkidx;  
  }
  if(iside == 2 ) {
    _vcc[1].setp( 0,  tvtx->_vtxcc[0] , tvtx->_vtxcc[1] );  
    if( tvtx->_ndof !=0 )  _vdd[1] = 1.0 * tvtx->_chi2/ tvtx->_ndof;
    _ntxtrack[1] = tvtx->_ntrack;
    chntrkidx[1] = vtkidx;
  }
}





#endif 
