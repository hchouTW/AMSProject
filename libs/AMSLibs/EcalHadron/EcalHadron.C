#ifndef __EcalHadron_C__
#define __EcalHadron_C__

#include "EcalHadron.h"

static double EdepThreshold = 1./1000.; // GeV / 1 to 10 MeV per cell is OK

double EcalHadron::ECALPAR[3][8]; //[18]
double EcalHadron::EcalEdepCell[18][72];
double EcalHadron::EcalEdepPla[18];
double EcalHadron::EcalS13Pla[18];
double EcalHadron::EcalS35Pla[18];

int    EcalHadron::EcalApex;
double EcalHadron::EcalEdepRaw;
double EcalHadron::EcalEdep;
double EcalHadron::EcalRigidity;


void EcalHadron::InitParameters(){ 
	// APEX-INDEPENDENT & LOGPARABOLIC-LIKE
	double CA[3]={1.05,1.00,0.95};
	double CB[3]={0.00,0.17,0.00};
	double PAR0_A1[8]={ 1.605973e+00, 1.423663e+00, 1.246288e+00, 1.110446e+00, 1.102703e+00, 1.026874e+00, 9.107874e-01, 9.131400e-01 }; 
	double PAR1_A1[8]={ 9.997233e-01, 9.281403e-01, 8.214123e-01, 9.916973e-01, 9.077081e-01, 9.885551e-01, 1.005648e+00, 9.812549e-01 }; 
	double PAR2_A1[8]={ 1.068166e-01, 1.283885e-01, 1.489230e-01, 1.344360e-01, 1.762905e-01, 2.244266e-01, 2.377332e-01, 1.540307e-01 };   
	for(int dd=0;dd<18;dd++){
		for(int zz=0;zz<8;zz++){
			ECALPAR[0][zz]= PAR0_A1[zz]*CA[0]+CB[0];
			ECALPAR[1][zz]= PAR1_A1[zz]*CA[1]+CB[1];
			ECALPAR[2][zz]= PAR2_A1[zz]*CA[2]+CB[2];
		}
	}
}


void EcalHadron::InitEvent(){

	for(int pp=0;pp<18;pp++){
		for(int cc=0;cc<72;cc++){
			EcalEdepCell[pp][cc] = 0.;
		}

		EcalEdepPla[pp]  = 0.;
		EcalS13Pla[pp]   = 0.;
		EcalS35Pla[pp]   = 0.;
	}

	EcalApex = 0;     // apex index
	EcalEdepRaw = 0.; // energyD GeV
	EcalEdep = 0.;    // erec below apex
	EcalRigidity= 0.; // estimated rig
}



void EcalHadron::SetShowerProperties(EcalShowerR* shower){
	EcalEdepRaw = shower->EnergyD/1000.; // GeV

	// --- shower decomposition --- 
	int N2D= shower->NEcal2DCluster(); // shower dec: 2D clusters
	for(int ecl2=0;ecl2<N2D;ecl2++){ 
		Ecal2DClusterR* e2d= (Ecal2DClusterR*)shower->pEcal2DCluster(ecl2); 
		if(!e2d) continue; 

		int N1D= e2d->NEcalCluster(); // shower dec: 1D clusters    
		for(int ecl1=0;ecl1<N1D;ecl1++){ 
			EcalClusterR* e1d= (EcalClusterR*)e2d->pEcalCluster(ecl1);

			int NHITS= e1d->NEcalHit(); // shower dec: raw hits
			for(int hi=0;hi<NHITS;hi++){
				EcalHitR* hit= (EcalHitR*)e1d->pEcalHit(hi);
				if(!hit) continue; 
				int Plane= hit->Plane;
				int Cell = hit->Cell;
				double Edep = hit->Edep;
				EcalEdepCell[Plane][Cell] = Edep/1000.; // GeV
			}  
		}
	}


	// --- recomposition: energy per plane ---
	for(int pp=0;pp<18;pp++){
		EcalEdepPla[pp] = 0.;
		for(int cc=0;cc<72;cc++){    
			if(EcalEdepCell[pp][cc] > EdepThreshold ){ 
				EcalEdepPla[pp]  += EcalEdepCell[pp][cc];
			}
		}
	}


	// --- find max-edep cell per plane ---
	int EcalCellSeed[18];
	for(int pp=0;pp<18;pp++){
		EcalCellSeed[pp] = -1; 
		double MaxEdep = 0.;
		if(EcalEdepPla[pp]==0) continue; 
		for(int cc=0;cc<72;cc++){    
			if(EcalEdepCell[pp][cc] > EdepThreshold ){ 

				if(EcalEdepCell[pp][cc] > MaxEdep){
					MaxEdep= EcalEdepCell[pp][cc];
					EcalCellSeed[pp]= cc;
				}
			}
		}
	}


	// --- ratio S13 and S15 ---
	for(int pp=0;pp<18;pp++){
		if(EcalEdepPla[pp]==0.) continue; 
		int Seed= EcalCellSeed[pp]; 
		if(Seed<0) continue;
		double S1 = EcalEdepCell[pp][ Seed ];
		double S3 = S1;
		if(Seed> 0) S3+= EcalEdepCell[pp][Seed - 1];
		if(Seed<71) S3+= EcalEdepCell[pp][Seed + 1];

		double S5= S3;
		if(Seed> 1) S5+= EcalEdepCell[pp][Seed - 2];
		if(Seed<70) S5+= EcalEdepCell[pp][Seed + 2];

		if(S3>0) EcalS13Pla[pp] = S1/S3;
		if(S5>0) EcalS35Pla[pp] = S3/S5;
	}

}


void EcalHadron::SetApex(){

	// dynamic threshold [from ~5MeV @1GV to ~100MeV @500GV]
	double EdepThresholdApex = EcalEdepRaw/1000.; 

	EcalApex = 19; // MIP default
	for(int pp=0;pp<18;pp++){   
		bool IsApex= (EcalS13Pla[pp]<0.930 && EcalS35Pla[pp]<0.980 && EcalEdepPla[pp]>EdepThresholdApex); 
		if(!IsApex) continue; 

		if(pp>0 && pp<17 && EcalEdepPla[pp]>EcalEdepPla[pp-1] && EcalEdepPla[pp]<0.75*EcalEdepPla[pp+1]){
			EcalApex = pp+1;
			break;
		}

		if(pp==0 && EcalEdepPla[pp]<0.75*EcalEdepPla[pp+1]){
			EcalApex = pp+1;
			break;
		}

		if(pp==17 && EcalEdepPla[pp]>EcalEdepPla[pp-1]){
			EcalApex = pp+1;
			break;
		}
	}


	// edep below apex
	EcalEdep=0.;
	if(EcalApex<1 || EcalApex>18){ EcalEdep = 0.;}
	else{
		for(int pp=EcalApex-1;pp<18;pp++){
			EcalEdep += EcalEdepPla[pp];
		}
	}

}



void EcalHadron::SetRigidity(int Z){
	// apex & edep/pla must already be set
	EcalRigidity=0.;
	if(EcalApex<1 || EcalApex>18) return;
	if(EcalEdep<=0.) return;

	int indZ = TMath::Min( TMath::Abs(Z),8 ) -1;
	//int indA = EcalApex -1; 
	double ShowerDepth= 18. - EcalApex + 1.;
	double EdepPlaNT= EcalEdep/ShowerDepth;

	double P[3];
	P[0]=ECALPAR[0][indZ]; //[indA];
	P[1]=ECALPAR[1][indZ]; //[indA];
	P[2]=ECALPAR[2][indZ]; //[indA];
	double LogEdep=log10(EdepPlaNT);
	double LogRig= P[0] + P[1]*LogEdep + P[2]*LogEdep*LogEdep;
	EcalRigidity = TMath::Power(10., LogRig);
}


void EcalHadron::Build(EcalShowerR* shower, int Z){
	InitEvent();
	SetShowerProperties(shower);
	SetApex();
	SetRigidity(Z);
}

#endif // __EcalHadron_C__
