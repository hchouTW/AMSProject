#ifndef TRDSTrack_H
#define TRDSTrack_H
#include <vector>


class fastFitQu{
  public:
    double cx[12];
    double cy[12];
    double cz[12];
    int nhit;
    fastFitQu(){
      nhit = 0;
    }
    void addYZ(double y, double z){
      if(nhit>=12) return;
      cy[nhit] = y;	
      cz[nhit] = z;
      nhit ++;	
    }
    void sort(){
      for(int i=0;i<nhit-1;i++){
        for(int j=0; j<nhit-i-1;j++){
          if(cz[j] > cz[j+1]){
            double tmp;
            tmp = cz[j];	cz[j] = cz[j+1];     	cz[j+1] = tmp;
            tmp = cy[j];	cy[j] = cy[j+1];     	cy[j+1] = tmp;
            tmp = cx[j];	cx[j] = cx[j+1];     	cx[j+1] = tmp;
          }}}
    }
    double guessY( double tz ) {
        return  Intpol1(cz[0],cz[nhit-1],cy[0],cy[nhit-1],tz ); 
    }
    double Intpol1(double x1, double x2, 
        double y1, double y2, double x) {
      return (x-x2)/(x1-x2)*y1
        +(x-x1)/(x2-x1)*y2;
    }
    void bkstack() {  if(nhit > 0) nhit--;  }	
    double getCsq() {
      sort();
      double sum =0;
      double sums = 0;
      for(int i=1;i<nhit-1;i++)
      {
        double tvec = getvec( cy[i-1],cz[i-1],cy[i],cz[i],cy[i+1],cz[i+1]);
        sum += tvec;
        sums += tvec*tvec;
      }
      return sqrt( sums/(nhit-2) - (sum/(nhit-2))*(sum/(nhit-2))); 
    }
    void clear(){nhit = 0;}
    double getMean(){
      sort();
      double sum =0;
      double sums = 0;
      for(int i=1;i<nhit-1;i++) {
        double tvec = getvec( cy[i-1],cz[i-1],cy[i],cz[i],cy[i+1],cz[i+1]);
        sum += tvec;
        sums += tvec*tvec;
      }
      return sum/(nhit-2);	
    }

    double getvec(double ay,double az,double by, double bz,double cy, double cz )		
    {
      double v1y = (ay-by)/sqrt((az-bz)*(az-bz) + (ay-by)*(ay-by) );
      double v1z = (az-bz)/sqrt((az-bz)*(az-bz) + (ay-by)*(ay-by) );
      double v2y = (cy-by)/sqrt((cz-bz)*(cz-bz) + (cy-by)*(cy-by) );
      double v2z = (cz-bz)/sqrt((cz-bz)*(cz-bz) + (cy-by)*(cy-by) );
      return	v1z*v2y - v1y*v2z;
    }
};


class semiIndex{
  public:
    int nLength;
    bool isBad ;
    bool isXTrack ;
    bool isVertex;
    int trkIndex[12];
    int levelIndex[12];
    float cxy[12];
    float cz[12];
    float quality;
    float param[2];
    semiIndex(){
      isVertex = false;
      nLength = 0;isBad = false;
      quality = -1;
      for(int i=0;i<12;i++) trkIndex[i]=-9999;
      for(int i=0;i<12;i++) levelIndex[i]=-9999;
      for(int i=0;i<12;i++) cxy[i]=-9999;
      for(int i=0;i<12;i++) cz[i]=-9999;
      isXTrack = false;
    }

    void AddHit(int iccL, int tlevel, float tcxy, float tcz, int sort=0){
      if(nLength>=12) return;
      trkIndex[nLength] = iccL; 
      levelIndex[nLength] = tlevel; 
      cxy[nLength] = tcxy;
      cz[nLength]  = tcz;
      nLength++;
      if(sort >0) SortHit();
    }
    void SortHit() {// sort hit by z
      for(int i=0;i<nLength-1; i++) {
        for(int j=0;j<nLength-i-1;j++) {
          int it; float ct;
          if( levelIndex[j]  > levelIndex[j+1] ) {
            it = levelIndex[j]; levelIndex[j] = levelIndex[j+1]; levelIndex[j+1]=it;
            it = trkIndex[j]; trkIndex[j] = trkIndex[j+1]; trkIndex[j+1]=it;
            ct = cxy[j]; cxy[j] = cxy[j+1]; cxy[j+1] = ct;
            ct = cz[j]; cz[j] = cz[j+1]; cz[j+1] = ct;
          }
        }
      }
    }
    
    
    void SetHit(float tcxy,  int iLength){
      cxy[iLength] = tcxy;
    }
    
    void PopBack(){
       if( nLength >0  )  {
         nLength--;
         trkIndex[nLength] = -9999;
         levelIndex[nLength] = -9999;
         cxy[nLength] = -9999;
         cz[nLength] = -9999;
       }
       return;    
    }
    void Clear()
    {
      for(int i=0;i<12;i++){ trkIndex[i] = -9999;}
      for(int i=0;i<12;i++){ levelIndex[i] = -9999;}
      for(int i=0;i<12;i++){ cxy[i] = -9999;}
      for(int i=0;i<12;i++){ cz [i] = -9999;}
      nLength = 0;
    }
    double Quality() {
//      if(quality >=0 ) return quality;
      if( ! check() ) return false;
      TrFit trfitk;
      fastFitQu trfit1; 
      for(int i=0;i<nLength;i++){
        trfitk.Add(0.5, cxy[i],cz[i] ,-1,0.6,0.6 );
        trfit1.addYZ( cxy[i],cz[i]  );
      }
      float ly = trfitk.LinearFit(2);
      float lm  = fabs( trfit1.getMean())  ;
      float lc  = fabs( trfit1.getCsq() )  ;
//      cout<<"for current candidate "<<Form(" %4.2f, %4.2f, %4.2f", ly, lm, lc ) <<endl;
      if( ly <= 0  ) {
        quality = 0;
        isBad = true;
        return quality;
      }
      else{
        quality = 20./(ly+20.)  + 100.* nLength ;
        return quality; 
      }
    }
    bool check(  float  qscut =9000. ) {
      fastFitQu trfit1; 
      TrFit trfitk;
      for(int i=0;i<nLength;i++){
        trfitk.Add(0.5, cxy[i],cz[i] ,-1,0.6,0.6 );
        trfit1.addYZ( cxy[i],cz[i]  );
      }
      float ly = trfitk.LinearFit(2);
      float lm  = fabs( trfit1.getMean())  ;
      float lc  = fabs( trfit1.getCsq() )  ;
//      cout<<"for current candidate "<<Form(" %4.2f, %4.2f, %4.2f", ly, lm, lc ) <<endl;
      for(int i=0;i<nLength;i++) {
        double cc = trfitk.GetYr(i) ;
        if(fabs(cc) >2.0 ) return false;
      }
      param[0] = trfitk.GetParam(0);
      param[1] = trfitk.GetParam(1);
      ly = 1.0*ly/nLength;
      if( ly<=0  ||  ly > qscut ) {
        isBad = true;
        return false;
      }
      return true;
    }
    float LFit( int iLength, float  tcxy ) {
      TrFit trfitk;
      for(int i=0;i<nLength;i++){
        if( i == iLength  )   trfitk.Add(0.5, tcxy,cz[i] ,-1,0.6,0.6 );
        else   trfitk.Add(0.5, cxy[i],cz[i] ,-1,0.6,0.6 );
      }
      float ly = trfitk.LinearFit(2);
      ly = 1.0*ly/nLength;
      return ly; 
    }
    float AFit( float  tcxy , float tcz ) {
      TrFit trfitk;
      if(  tcz > cz[0]  )   trfitk.Add(0.5,tcxy,tcz ,-1,0.6,0.6 );
      for(int i=0;i<nLength;i++){
        if( tcz < cz[i] && tcz > cz[i+1] && i<nLength-1)  trfitk.Add(0.5,tcxy,tcz ,-1,0.6,0.6 );
        trfitk.Add(0.5, cxy[i],cz[i] ,-1,0.6,0.6 );
      }
      if(  tcz > cz[nLength-1]  )   trfitk.Add(0.5,tcxy,tcz ,-1,0.6,0.6 );
      
      float ly = trfitk.LinearFit(2);
      ly = 1.0*ly/nLength;
      if( ly <=0  ) return -1;
      float res = 0;
      for(int i=0;i<nLength+1;i++){
        float cc = trfitk.GetYr(i);
        if( fabs(cc) == 0 ) return -1;
        if(  fabs(cc)  > res   ) res = fabs(cc);
      }
      if( res > 1 )  return -1;
      return ly; 
    }
    
    void Stat() {
      fastFitQu trfit1; 
      TrFit trfitk;
      for(int i=0;i<nLength;i++){
        trfitk.Add(0.5, cxy[i],cz[i] ,-1,0.6,0.6 );
        trfit1.addYZ( cxy[i],cz[i]  );
      }
      float ly = trfitk.LinearFit(2);
      float lm  = fabs( trfit1.getMean())  ;
      float lc  = fabs( trfit1.getCsq() )  ;
      cout<<"for current candidate "<<Form(" %4.2f, %4.2f, %4.2f", ly, lm, lc ) <<endl;
      for(int i=0;i<nLength;i++) {
        double cc = trfitk.GetYr(i) ;
        cout<<fabs( cc )<<" , ";
      }
      cout<<endl;
    }




};





class ccL{//coordinate list
  public:
    double cz,   cxy;
    int nLayer,  idx;
    double amp, ampy;
    bool isused;
    bool isWide;
    double cst,  ced;
    int nadj;
    int ajamp[50];
    int ajampy[50];
    int ajidx[50];
    double ajcxy[50];
     

    ccL(){  cxy = -9999; cz = -9999; nLayer = -9999; idx = -1; amp = -9999; ampy = -9999;isused = false; cst = -9999; ced = -9999;  nadj=0; isWide = false; }
    void  addCC( int tidx, int tnLayer, double tcxy, double tcz, double tamp = -9999, double tampy = -9999 ) {
        idx = tidx; nLayer = tnLayer;
        cxy = tcxy; cz = tcz;
        amp = tamp; ampy = tampy;
        cst = cxy;  ced = cxy;
        
        ajidx [nadj] = tidx;
        ajamp [nadj] = tamp;
        ajampy[nadj] = tampy;
        ajcxy [nadj] = tcxy;
        nadj++;

    }
    int addAj(  int tidx, int tnLayer, double tcxy, double tcz, double tamp = -9999, double tampy=-9999 ){

      float kl = 0.8 *2.0;
//       if( tnLayer != nLayer || tcxy < cst -0.8 || tcxy > ced + 0.8  ) return -1;
       if( tnLayer != nLayer || tcxy < cst - kl || tcxy > ced + kl  ) return -1;
//      if( tnLayer != nLayer || tcxy < cst - 3 || tcxy > ced + 3  ) return -1;
      else  {
        if(  tcxy < cst   ) cst = tcxy;
        if(  tcxy > ced   ) ced = tcxy;
        ajidx [nadj] = tidx;
        ajamp [nadj] = tamp;
        ajampy[nadj] = tampy;
        ajcxy [nadj] = tcxy;
        nadj ++; 
    
        if( fabs(cst-ced) > 1.5 ) isWide = true;
        // process the amp and ampy;
        float totamp = 0;
        for(int i=0;i<nadj;i++) totamp += ajamp[i];
        float totcxy = 0;
        for( int i=0;i<nadj;i++ ) {
          totcxy +=  ajcxy[i] * ajamp[i]/totamp;
        }
        cxy = totcxy;
        return 1;
      }
      return 0;
    }
};

class vcos{
  public:
    vcos(){}
    double getcosyz(double *a,double *b,double *c)
    {
      return ((a[1]-b[1])*(c[1]-b[1])+(a[2]-b[2])*(c[2]-b[2]))/
        (sqrt((a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))*
         sqrt((c[1]-b[1])*(c[1]-b[1])+(c[2]-b[2])*(c[2]-b[2])));
    }
    double getcosyz(ccL  a, ccL b, ccL c)
    {
      return
        sqrt((a.cxy-b.cxy)*(a.cxy-b.cxy)+(a.cz-b.cz)*(a.cz-b.cz)) == 0 || 
        sqrt((c.cxy-b.cxy)*(c.cxy-b.cxy)+(c.cz-b.cz)*(c.cz-b.cz)) ==0  ? -1:
        ((a.cxy-b.cxy)*(c.cxy-b.cxy)+(a.cz-b.cz)*(c.cz-b.cz))/
        sqrt((a.cxy-b.cxy)*(a.cxy-b.cxy)+(a.cz-b.cz)*(a.cz-b.cz))/
        sqrt((c.cxy-b.cxy)*(c.cxy-b.cxy)+(c.cz-b.cz)*(c.cz-b.cz));
    }
};


class CombIter{
  public:
    CombIter( vector<vector<int> > &tcomb,  vector<ccL> &tccList, float tmaxcos ){ 
      cpos[0] = 0; cpos[1] =0;
      comb = tcomb; 
      tbecon = 0;
      ccList = tccList;
      for(int i=0;i<3;i++) 
        ibecon[i] = -1;
      trdebug = false;
      maxcos = tmaxcos;
    }
    ~CombIter(){ 
      for(int i=0;i<comb.size();i++) comb[i].clear(); 
      comb.clear(); 
    }
    int nCandi;
    vector<semiIndex> candi;
    semiIndex tsemi;
    int tbecon;
    vector<vector<int> > comb;
    vector<ccL> ccList;
    bool trdebug;
    float maxcos; 
    int ibecon[3];
    int lmax;
    vcos headcos;
    int cpos[2];
    int  CheckOverlap(   int ilevel = -1   ) {
      // before checking ,we should order the quality first.
      // tracks can only be check overlap from top to bottom.
      if(  candi.size() < 2  ) return -1;
      vector<pair<double,int> >  canidx;
      for(int i=0;i<candi.size();i++)
      vector<pair<double,int> >tindex;
//      cout<<"candi size is "<<candi.size()<<endl;
      for(int i=0;i<candi.size();i++)
        canidx.push_back( make_pair<double,int>( candi[i].quality,i) );
      for(int i=0;i<canidx.size()-1;i++){
        for(int j=0;j<canidx.size()-i-1;j++){
          if( canidx[j].first < canidx[j+1].first ) {
            pair<double,int> ttp = canidx[j];
            canidx[j]  = canidx[j+1];
            canidx[j+1]= ttp;
          }
        }
      }
//      for(int i=0;i<canidx.size();i++)
//         cout<<"  %%%%%%    "<<  canidx[i].first<<"  ,  "<< canidx[i].second<<endl;
//      cout<<endl; 
      for(int i=0;i<canidx.size()-1; i++) {
        int ii = canidx[i].second;
        if(candi[ii].isBad  ) continue;
        for(int j=i+1;j<canidx.size();j++) {
          if( j<=i) continue;
          int jj = canidx[j].second;
          if(candi[jj].isBad  ) continue;
          if(CheckOverlap(  candi[ii], candi[jj] )  <0 ) {
             MaskGhost( candi[ii], candi[jj] );
          }
        }
      }
      return 0;
    }
    int CheckOverlap(  int idx0, int idx1 ) {
      for(int itk=0;itk< candi.size(); itk++ )  {
        if(candi[itk].isBad ) continue;
        for(int i=0;i< candi[i].nLength; i++ ) {
          if(  idx0 == candi[itk].trkIndex[i] ||  idx1 == candi[itk].trkIndex[i]  )
            return -1;
        }
      }
      return 1;
    }
    int CheckOverlap(  semiIndex &gsemi1, semiIndex &gsemi2 ) {
      for(int i=0;i< gsemi1.nLength; i++ ) {
        int ilevel = gsemi1.levelIndex[i];
        for(int j=0;j<gsemi2.nLength;j++ )  {
          int jlevel = gsemi2.levelIndex[j];
          if(ilevel != jlevel) continue;
          if( gsemi1.trkIndex[i] == gsemi2.trkIndex[j]  ) return -1;
        }
      }
      return 1;
    }

    int MaskGhost(  semiIndex &gsemi1, semiIndex &gsemi2 ) {
      int ires = 0; 
      if( gsemi1.isBad || gsemi2.isBad )   return ires;
      if( gsemi1.Quality() >= gsemi2.Quality() )   {
        gsemi2.isBad = true;
        ires =1;
      }
      else {
        ires =2;
        gsemi1.isBad = true;
      }
      return ires;
    }

    int selfiter( int ilevel, int tlmax ) {
      lmax = tlmax;
//      if( lmax <= 5 ) return 0;
      if( ilevel != lmax  ){
        if(ibecon[1]>=0 )  ibecon[0] = ibecon[1];
        if(ibecon[2]>=0 )  ibecon[1] = ibecon[2];
      }
      if( ilevel == 1 ) {
        int nres = 0;
        for(int i=0;i< comb.size(); i++) {
          if( comb[i][0] != ilevel   ) continue; 
          if( comb[i][1] != ibecon[1]   ) continue;
          if( ibecon[0] <0 ||  ibecon[1] <0 ) break;
          if( headcos.getcosyz( ccList[ibecon[0]],  ccList[ibecon[1]] , ccList[comb[i][2]] ) >maxcos  )  continue;
          ibecon[2] = comb[i][2];
          nres ++;
          // here we need to add something to current candidate.
          tsemi.AddHit( ibecon[2] , ccList[ibecon[2]].nLayer, ccList[ibecon[2]].cxy, ccList[ibecon[2]].cz);
          if(tsemi.Quality() > 0)  candi.push_back( tsemi );    
          tsemi.PopBack();
        }
        return nres;
      }
      else if( ilevel > 1 )   {
        int nres = 0 ;
        for(int i=0;i< comb.size(); i++) {
          if( comb[i][0] != ilevel   ) continue; 
          if(ilevel == lmax  ) {
            for(int ii=0;ii<3;ii++) 
              ibecon[ii] = -1;
            ibecon[0] = comb[i][1];
            ibecon[1] = comb[i][2];
            tsemi.Clear();
            tsemi.AddHit( ibecon[0], ccList[ibecon[0]].nLayer,ccList[ibecon[0]].cxy, ccList[ibecon[0]].cz);
            tsemi.AddHit( ibecon[1], ccList[ibecon[1]].nLayer ,ccList[ibecon[1]].cxy, ccList[ibecon[1]].cz);
            int nn =  selfiter(ilevel-1,lmax);
            if( nn <= 0  ) continue;
            continue;
          }
          int tres = 0;
          
          if(   comb[i][1] != ibecon[1]   ) continue;
          if( ibecon[0] <0 ||  ibecon[1] <0 ) break;
          if( headcos.getcosyz( ccList[ibecon[0]],  ccList[ibecon[1]] , ccList[comb[i][2]] ) > maxcos  )  continue;
          ibecon[2] = comb[i][2];
          nres ++;
          tsemi.AddHit( ibecon[2], ccList[ibecon[2]].nLayer, ccList[ibecon[2]].cxy, ccList[ibecon[2]].cz );
          int nlower = selfiter(ilevel-1,lmax);
          tsemi.PopBack();
          if( nlower <=0  ) continue;
        }// hit combination 
        if(0)
        if(ilevel == lmax ) {  // we need to remove all the overlap when iteration is over
          if( candi.size() >1  )
          for(int ic=0;ic<candi.size()-1;ic++) {  
            if( candi[ic].nLength <  lmax + 1  ) { candi[ic].isBad = true;  continue;}
            if( candi[ic].nLength != lmax + 1  ) continue;
            for(int jc=ic+1;jc<candi.size();jc++) {
              if( jc<=ic  ) continue;
              if( candi[jc].nLength != lmax + 1  ) continue;
              if( candi[jc].nLength <  lmax + 1  ) { candi[jc].isBad = true;  continue;}
              if(  CheckOverlap(  candi[ic], candi[jc] )  < 0  )  // should be improved, sort the order first
                MaskGhost( candi[ic], candi[jc] );
            }
          }
        }
        return nres;
      }
      return 0;
    } //self iter

    
};










#endif
