#ifndef __MgntRecon_H__
#define __MgntRecon_H__

#include <stdarg.h>

#include "MgntPhySt.h"
#include "MgntHitSt.h"
#include "MgntMag.h"
#include "MgntMat.h"
#include "MgntProp.h"

//---- Track ----//
class Track {
	public :
		enum MethodOption {
			ALCARAZ = 0, CHOUTKO = 1, CHIKANIAN = 2,
			ANALYTIC = 3, SIMPLE = 4, PHYSICAL = 5
		};
		
		enum PhyOption {
			NONE = 0, MSCAT = 1, ENGLS = 2, MSCAT_ENGLS = 3
		};

		enum DirOption {
			UPWARD = 0, DOWNWARD = 1
		};

	public :
		Track(PhySt::ParticleList part = PhySt::Proton, DirOption dirOpt = Track::DOWNWARD) 
			{ MgntMag::Load(); MgntMat::Load(); setPart(part, dirOpt); }
		~Track() {}

		void setPart(PhySt::ParticleList part = PhySt::Proton, DirOption dirOpt = Track::DOWNWARD) 
			{ fPartSt.setPart(part); fDirOpt = dirOpt; }
		
		void addHits(std::vector<HitSt>& hits) 
			{ fHits = hits; HitSt::SortHits(fHits, fDirOpt); }
		
		Bool_t fit(MethodOption method = Track::SIMPLE, PhyOption phyOpt = Track::MSCAT_ENGLS);
		
		void print() {
			COUT("******************************** Track ********************************\n");
			fPartSt.print();
			for (Int_t ih = 0; ih < fHits.size(); ++ih) fHits.at(ih).print();
			COUT("***********************************************************************\n");
		}

		PhySt PartStatus() { return fPartSt; }

	//protected :
	public :
		Bool_t officalFit(MethodOption method, PhyOption phyOpt);
		Bool_t analyticFit();
		Bool_t simpleFit();
		Bool_t physicalFit();

	//protected :
	public :
		MethodOption fMethod;
		PhyOption    fPhyOpt;
		DirOption    fDirOpt;
		PhySt        fPartSt; // physical status at (Z = outside of top or bottom hit)
		
		std::vector<HitSt> fHits; // all hits
		std::vector<MatPropParam> fMatPars; // physical parameters between two hit

		Int_t    fNdf;
		Double_t fChisq;
		Double_t fNormChisq;

		Double_t fChisqHX;
		Double_t fChisqHY;
		Double_t fChisqDV;
		Double_t fChisqDW;
		Double_t fChisqIN;
		Double_t fChisqBR;

		// Physical algorithm (Negative Log-Likelihood)
		SVecI<2> fNPMHit; // (HitX, HitY)
		SVecD<2> fNLLHit; // (HitX, HitY) 
		SVecI<4> fNPMInt; // (MscatV, MscatW, EnglsI, EnglsB)
		SVecD<4> fNLLInt; // (MscatV, MscatW, EnglsI, EnglsB)
		
		Int_t           fNDF;     // (fNPMHit + fNPMInt - 5)
		Double_t        fNLL;     // (fNLLHit + fNLLInt)
		Double_t        fNormNLL; // (NLL / NDF)

		Int_t           fNDFX;     // (fNPMHit + fNPMInt - 2)
		Double_t        fNLLX;     // (fNLLHit + fNLLInt)
		Double_t        fNormNLLX; // (NLLX / NDFX)

		Int_t           fNDFY;     // (fNPMHit + fNPMInt - 3)
		Double_t        fNLLY;     // (fNLLHit + fNLLInt)
		Double_t        fNormNLLY; // (NLLY / NDFY)

		void NLLZero() {
			fNPMHit = SVecI<2>();
			fNLLHit = SVecD<2>();
			fNPMInt = SVecI<4>();
			fNLLInt = SVecD<4>();
			fNDF     = 0;
			fNLL     = ZERO;
			fNormNLL = ZERO;
			
			fNDFX     = 0;
			fNLLX     = ZERO;
			fNormNLLX = ZERO;
			fNDFY     = 0;
			fNLLY     = ZERO;
			fNormNLLY = ZERO;

			fNdf = 0;
			fChisq = ZERO;
			fNormChisq = ZERO;
		}

		//std::vector<PhySt>        fPartSt;
		//std::vector<MatPropParam> fIntSt; // physical parameters between two hit

	//protected :
	public :
		static const Double_t ZERO = 0.;
		static const Double_t ONE =  1.;
		static const Double_t NEG = -1.;

		// Num of Hit requirement
		static const Int_t    NUM_MIN_HITX = 3;
		static const Int_t    NUM_MIN_HITY = 4;
		
		// Momentum
		static const Double_t EPSILON_MIN_MOM = 1.0e-4;
		static const Double_t EPSILON_RAT_MOM = 0.15;

		// Coordinate
		static const Double_t EPSILON_MAX_COO = 500.;

		// Minimization
		static const Int_t    NUM_MAX_ITER_FACT_ANALYTIC = 10;
		static const Int_t    NUM_MAX_ITER_FACT_SIMPLE   = 10;
		static const Int_t    NUM_MAX_ITER_FACT_PHYSICAL = 10;
		static const Double_t CONVG_EPSILON   = 1.0e-3;
		static const Double_t CONVG_TOLERANCE = 1.0e-1;

		// Levenberg-Marquardt algorithm
		static const Int_t    LM_IGNORE_ITER  = 12; // best 5
		static const Int_t    LM_NUM_MAX_ITER = 15;
		static const Double_t LM_LAMBDA_0     = 1.00;
		static const Double_t LM_LAMBDA_UP    = 1.00;
		static const Double_t LM_LAMBDA_DN    = 0.25;
};





















/*
//---- Vertex ----//
class Vertex {
	public :
		Vertex() { init(); }
		~Vertex() {}

		void init();
		void analyticsFit();
		void sampleFit();
		void physicalFit();

	protected :
		void minimizer();

	protected :
		std::vector<Track &> fOrgTracks;
		std::vector<Track>   fRecTracks;
};


//---- MgntRecon ----//
class MgntRecon {
  public :
    MgntRecon() {}
    ~MgntRecon() {}

	protected :
		std::vector<Track>  fTracks;
		std::vector<Vertex> fVertexs;
};
*/
#endif // __MgntRecon_H__
