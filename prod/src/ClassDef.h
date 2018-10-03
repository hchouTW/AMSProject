/**********************************************************
 * Author        : Hsin-Yi Chou
 * Email         : hchou@cern.ch
 * Last modified : 2015-07-22 16:29
 * Filename      : ClassDef.h
 * Description   :
 * *******************************************************/
#ifndef __ClassDef_H__
#define __ClassDef_H__

#include <iostream>
#include <vector>
#include <TObject.h>
#include <TString.h>

// RunTagInfo
class RunTagInfo : public TObject {
	public :
		RunTagInfo() { init(); }
		~RunTagInfo() {}

		void init() {
			run = 0;
			eventFT = 0;
			eventLT = 0;
			numOfSelEvent = 0;
			numOfTrgEvent = 0;
			file.clear();
			
			dateUTC = 0;
			timeUTC = 0;
		}

	public :
		UInt_t run;
		UInt_t eventFT;            // first trigger event
		UInt_t eventLT;            // last trigger event
		UInt_t numOfSelEvent;      // number of select  events
		UInt_t numOfTrgEvent;      // number of trigger events
		std::vector<TString> file; // file list

		UInt_t  dateUTC;           // UTC date
		UInt_t  timeUTC;           // UTC time

	ClassDef(RunTagInfo, 6)
};


// SegPARTMCInfo
class SegPARTMCInfo : public TObject {
	public :
		SegPARTMCInfo() { init(); }
		~SegPARTMCInfo() {}

		void init() {
			lay = -1;
			bta = -1;
			mom = -1;
			std::fill_n(coo, 3, 0);
			std::fill_n(dir, 3, 0);
		}

	public :
        Short_t lay;
		Float_t bta;
		Float_t mom;
		Float_t coo[3];
		Float_t dir[3];
	
    ClassDef(SegPARTMCInfo, 3)
};

struct SegPARTMCInfo_sort {
	bool operator() (const SegPARTMCInfo & hit1, const SegPARTMCInfo & hit2) {
		if      (hit1.lay < hit2.lay) return true;
		else if (hit1.lay > hit2.lay) return false;
		return false;
	}
};


// HitTRKMCInfo
class HitTRKMCInfo : public TObject {
	public :
		HitTRKMCInfo() { init(); }
		~HitTRKMCInfo() {}

		void init() {
			layJ = 0;
			edep = -1;
			mom  = -1;
			std::fill_n(coo, 3, 0);
            std::fill_n(smr, 2, 0);
			std::fill_n(loc, 2, -1);
		}

	public :
		Short_t layJ;    // layerJ
		Float_t edep;    // edep
		Float_t mom;     // momentum
		Float_t coo[3];
        Float_t smr[2];
        Float_t loc[2];

	ClassDef(HitTRKMCInfo, 7)
};

struct HitTRKMCInfo_sort {
	bool operator() (const HitTRKMCInfo & hit1, const HitTRKMCInfo & hit2) {
		if      (hit1.layJ < hit2.layJ) return true;
		else if (hit1.layJ > hit2.layJ) return false;
		else {
			if      (hit1.mom < hit2.mom) return true;
			else if (hit1.mom > hit2.mom) return false;
		}
		return false;
	}
};


// HitTOFMCInfo
class HitTOFMCInfo : public TObject {
	public :
		HitTOFMCInfo() { init(); }
		~HitTOFMCInfo() {}

		void init() {
			lay = 0;
			bta = 0.;
			std::fill_n(coo, 3, 0);
		}

	public :
		Short_t lay;    // layer
		Float_t bta;    // beta
		Float_t coo[3]; // coordinate
	
    ClassDef(HitTOFMCInfo, 1)
};

struct HitTOFMCInfo_sort {
	bool operator() (const HitTOFMCInfo & hit1, const HitTOFMCInfo & hit2) {
		if      (hit1.lay < hit2.lay) return true;
		else if (hit1.lay > hit2.lay) return false;
		return false;
	}
};


// HitTRDMCInfo
class HitTRDMCInfo : public TObject {
	public :
		HitTRDMCInfo() { init(); }
		~HitTRDMCInfo() {}

		void init() {
			lay = 0;
            mom = -1;
			std::fill_n(coo, 3, 0);
		}

	public :
		Short_t lay;    // layer
        Float_t mom;
		Float_t coo[3];

	ClassDef(HitTRDMCInfo, 6)
};

struct HitTRDMCInfo_sort {
	bool operator() (const HitTRDMCInfo & hit1, const HitTRDMCInfo & hit2) {
		if      (hit1.lay < hit2.lay) return true;
		else if (hit1.lay > hit2.lay) return false;
		return false;
	}
};


// PartMCInfo
class PartMCInfo : public TObject {
	public :
		PartMCInfo() { init(); }
		~PartMCInfo() {}

		void init() {
			chrg = 0;
			mass = 0;
            bta  = 0;
			mom  = 0;
			ke   = 0;
			std::fill_n(coo, 3, 0);
			std::fill_n(dir, 3, 0);
            
            segsTk.clear();
            segsTf.clear();
            segsTd.clear();
            segsEc.clear();
            segsRh.clear();
			
            hitsTk.clear();
            hitsTf.clear();
		}

	public :
		Float_t chrg;
		Float_t mass;
        Float_t bta;
		Float_t mom;
		Float_t ke; // kinetic energy
		Float_t coo[3];
		Float_t dir[3];

        std::vector<SegPARTMCInfo> segsTk;
        std::vector<SegPARTMCInfo> segsTf;
        std::vector<SegPARTMCInfo> segsTd;
        std::vector<SegPARTMCInfo> segsEc;
        std::vector<SegPARTMCInfo> segsRh;
		
        std::vector<HitTRKMCInfo>  hitsTk;
        std::vector<HitTOFMCInfo>  hitsTf;

	ClassDef(PartMCInfo, 9)
};

struct PartMCInfo_sort {
	bool operator() (const PartMCInfo & par1, const PartMCInfo & par2) {
		if      (par1.coo[2] > par2.coo[2]) return true;
		else if (par1.coo[2] < par2.coo[2]) return false;
		else {
			if      (par1.mom > par2.mom) return true;
			else if (par1.mom < par2.mom) return false;
		}
		return false;
	}
};


// VertexMCInfo
class VertexMCInfo : public TObject {
	public :
		VertexMCInfo() { init(); }
		~VertexMCInfo() {}

		void init() {
			status = false;
			std::fill_n(coo, 3, 0.);
            numOfPart = 0;
            partId = -1;
            partKe = -1;
		}

	public :
		Bool_t  status;
		Float_t coo[3];

        Short_t numOfPart;
        Short_t partId;
        Float_t partKe;
	
	ClassDef(VertexMCInfo, 7)
};


// HitTRKInfo
class HitTRKInfo : public TObject {
	public :
		HitTRKInfo() { init(); }
		~HitTRKInfo() {}

		void init() {
			std::fill_n(clsId,  2, -1);
			layJ = -1;
			tkid =  0;
			sens = -1;
			mult = -1;
			std::fill_n(side, 2, false);
			std::fill_n(coo, 3, 0);
			std::fill_n(loc, 2, -1);
			std::fill_n(chrg, 2, -1);

            std::fill_n(nsr, 2, 0);
            std::fill_n(sig[0], 2*5, 0);
		}

	public :
		Short_t clsId[2]; // cluster id (x, y)
		Short_t layJ;     // layerJ
		Short_t tkid;     // tkID
		Short_t sens;     // sensor
		Short_t mult;     // multiplicity
		Bool_t  side[2];  // side, x, y
		Float_t coo[3];   // coordinate
		Float_t loc[2];   // (elc) cofg loc
		Float_t chrg[2];  // (elc) chrg

        Short_t nsr[2];    // num of strip
        Float_t sig[2][5]; // index 2 -> seed

	ClassDef(HitTRKInfo, 10)
};

struct HitTRKInfo_sort {
	bool operator() (const HitTRKInfo & hit1, const HitTRKInfo & hit2) {
		if      (hit1.layJ < hit2.layJ) return true;
		else if (hit1.layJ > hit2.layJ) return false;
		else {
			if      (hit1.tkid < hit2.tkid) return true;
			else if (hit1.tkid > hit2.tkid) return false;
		}
		return false;
	}
};


// HitTRDInfo
class HitTRDInfo : public TObject {
	public :
		HitTRDInfo() { init(); }
		~HitTRDInfo() {}

		void init() {
			lay  = -1;
			side = 0;
			amp  = 0;
			len  = 0;
			std::fill_n(coo, 3, 0);

            dEdx  = -1;
            mcMom = -1;
		}

	public :
		Short_t lay;    // layer
		Short_t side;   // side, 1 x, 2 y, 3 xy
		Float_t amp;    // (elc) amp
		Float_t len;    // (elc) len
		Float_t coo[3];

        Float_t dEdx;   // dEdx = 0.01 * (amp/len)
        Float_t mcMom;  // (MC Info) mom

	ClassDef(HitTRDInfo, 7)
};

struct HitTRDInfo_sort {
	bool operator() (const HitTRDInfo & hit1, const HitTRDInfo & hit2) {
		if      (hit1.lay < hit2.lay) return true;
		else if (hit1.lay > hit2.lay) return false;
        else {
		    if      (hit1.amp < hit2.amp) return true;
		    else if (hit1.amp > hit2.amp) return false;
        }
		return false;
	}
};


// HitRICHInfo
class HitRICHInfo : public TObject {
	public :
		HitRICHInfo() { init(); }
		~HitRICHInfo() {}

		void init() {
            status = 0;
            cross = false;
			channel = -1;
            npe = -1;
            
			std::fill_n(coo, 3, 0);
            dist = -1.;
            dtha = -1.;
            dphi = -1.;
            rtha = -1.;
            rphi = -1.;
		}

	public :
        Short_t status;
        Bool_t  cross;   // cross
		Short_t channel; // channel = 16 * PMT + pixel
        Float_t npe;     // ADC counts above the pedestal/gain of the channel.
		Float_t coo[3];  // x, y, z
        Float_t dist;    // distance to particle
        Float_t dtha;    // direct theta
        Float_t dphi;    // direct phi (pmt)
        Float_t rtha;    // reflect theta
        Float_t rphi;    // reflect phi (pmt)

	ClassDef(HitRICHInfo, 1)
};

struct HitRICHInfo_sort {
	bool operator() (const HitRICHInfo & hit1, const HitRICHInfo & hit2) {
		if      (hit1.channel < hit2.channel) return true;
		else if (hit1.channel > hit2.channel) return false;
		return false;
	}
};


// HitECALInfo
class HitECALInfo : public TObject {
	public :
		HitECALInfo() { init(); }
		~HitECALInfo() {}

		void init() {
			id   = 0;
			side = 0;
			edep = -1;
			std::fill_n(coo, 3, 0);
		}

	public :
		Short_t id;     // plane * 100 + cell
		Bool_t  side;   // side 0 x, 1 y
		Float_t edep;   // (phy) edep
		Float_t coo[3]; // x, y, z

	ClassDef(HitECALInfo, 2)
};

struct HitECALInfo_sort {
	bool operator() (const HitECALInfo & hit1, const HitECALInfo & hit2) {
		if      (hit1.id < hit2.id) return true;
		else if (hit1.id > hit2.id) return false;
		return false;
	}
};


// ClsACCInfo
class ClsACCInfo : public TObject {
	public :
		ClsACCInfo() { init(); }
		~ClsACCInfo() {}

		void init() {
			sector = 0;
			time = 0;
			rawQ = -1;
			std::fill_n(coo, 3, 0);
		}

	public :
		Short_t sector; // sector (1~8)
		Float_t time;   // time
		Float_t rawQ;   // raw Charge
		Float_t coo[3]; // x, y, z

	ClassDef(ClsACCInfo, 1)
};

struct ClsACCInfo_sort {
	bool operator() (const ClsACCInfo & cls1, const ClsACCInfo & cls2) {
		return (cls1.time < cls2.time);
	}
};


// CKTrackInfo
class CKTrackInfo : public TObject {
	public :
		CKTrackInfo() { init(); }
		~CKTrackInfo() {}

		void init() {
            status = false;
            std::fill_n(ndof, 2, 0);
            std::fill_n(nchi, 2, 0);
            
            rig = 0;
            bta = 0;

            statusTop = false;
            std::fill_n(stateTop, 6, 0);
            
            statusBtm = false;
            std::fill_n(stateBtm, 6, 0);
            
            std::fill_n(statusLJ, 9, false);
            std::fill_n(stateLJ[0], 9*6, 0);

            cpuTime = 0;
        }

    public :
        Bool_t  status;
        Short_t ndof[2];
		Float_t nchi[2];
        
        Float_t rig;
        Float_t bta;
        
        Bool_t  statusTop;
        Float_t stateTop[6]; // (cx cy cz ux uy uz)
        
        Bool_t  statusBtm;
        Float_t stateBtm[6]; // (cx cy cz ux uy uz)

        Bool_t  statusLJ[9];
        Float_t stateLJ[9][6]; // track state at layerJ (1 2 3 4 5 6 7 8 9)

        Float_t cpuTime; // [ms]
	
    ClassDef(CKTrackInfo, 1)
};


// KFTrackInfo
class KFTrackInfo : public TObject {
	public :
		KFTrackInfo() { init(); }
		~KFTrackInfo() {}

		void init() {
            status = false;
            std::fill_n(ndof, 2, 0);
            std::fill_n(nchi, 2, 0);
            
            std::fill_n(rig, 4, 0);
            std::fill_n(bta, 4, 0);

            statusTop = false;
            std::fill_n(stateTop, 8, 0);
            
            statusBtm = false;
            std::fill_n(stateBtm, 8, 0);
            
            std::fill_n(statusLJ, 9, false);
            std::fill_n(stateLJ[0], 9*8, 0);

            cpuTime = 0;
        }

    public :
        Bool_t  status;
        Short_t ndof[2];
		Float_t nchi[2];
       
        Float_t rig[4]; // z = 195, 0, -70, -136
        Float_t bta[4]; // z = 195, 0, -70, -136

        Bool_t  statusTop;
        Float_t stateTop[8]; // (cx cy cz ux uy uz rig bta)
        
        Bool_t  statusBtm;
        Float_t stateBtm[8]; // (cx cy cz ux uy uz rig bta)

        Bool_t  statusLJ[9];
        Float_t stateLJ[9][8]; // track state at layerJ (1 2 3 4 5 6 7 8 9)
	
        Float_t cpuTime; // [ms]

    ClassDef(KFTrackInfo, 1)
};


// HCTrackInfo
class HCTrackInfo : public TObject {
	public :
		HCTrackInfo() { init(); }
		~HCTrackInfo() {}

		void init() {
            status = false;
            going = 0;
            chrg = 0;
            mass = 0;
            
            std::fill_n(ndof, 2, 0);
            std::fill_n(nchi, 2, 0);
            std::fill_n(quality, 2, 0);
            
            std::fill_n(state, 8, 0);
            
            std::fill_n(rig, 4, 0);
            std::fill_n(bta, 4, 0);

            statusTop = false;
            std::fill_n(stateTop, 8, 0);
            
            statusCen = false;
            std::fill_n(stateCen, 8, 0);
            
            statusRh = false;
            std::fill_n(stateRh, 8, 0);
            
            statusBtm = false;
            std::fill_n(stateBtm, 8, 0);
            
            std::fill_n(statusLJ, 9, false);
            std::fill_n(stateLJ[0], 9*8, 0);
           
            cpuTime = 0;
        }
	
    public :
        Bool_t  status;
        Short_t going; // 0, no  1, up-going  -1, down-going
        Short_t chrg;
        Float_t mass;
        
        Short_t ndof[2];
        Float_t nchi[2];
        Float_t quality[2];
        
        Float_t state[8]; // (cx cy cz ux uy uz rig bta)
        
        Float_t rig[4]; // z = 195, 0, -70, -136
        Float_t bta[4]; // z = 195, 0, -70, -136
        
        Bool_t  statusTop; // track at top of detector (z = 195.)
        Float_t stateTop[8];
        
        Bool_t  statusCen; // track at bottom of detector (z = 0.)
        Float_t stateCen[8];
        
        Bool_t  statusRh; // track at bottom of detector (z = -70.)
        Float_t stateRh[8];
        
        Bool_t  statusBtm; // track at bottom of detector (z = -136.)
        Float_t stateBtm[8];

        Bool_t  statusLJ[9];
        Float_t stateLJ[9][8]; // track state at layerJ (1 2 3 4 5 6 7 8 9)

        Float_t cpuTime; // [ms]

        ClassDef(HCTrackInfo, 6)
};


// ShowerInfo
class ShowerInfo : public TObject {
	public :
		ShowerInfo() { init(); }
		~ShowerInfo() {}
		
		void init() {
			energyD = -1;
			energyE = -1;
			PisaBDT = -2;
			Q = -1;

			hadronApex = -1;
			hadronEnergy = -1;
		}

	public :
		Float_t energyD;  // energy deposit [GeV]
		Float_t energyE;  // Pisa
		Float_t PisaBDT;  // Pisa
		Float_t Q;        // ecal charge

		Short_t hadronApex; // reject high-apex events (say, Apex > 12)
		Float_t hadronEnergy;

		ClassDef(ShowerInfo, 4)
};


// LIST
class LIST : public TObject {
	public :
		LIST() { init(); }
		~LIST() {}

		void init() {
			run    = 0;
			event  = 0;
			entry  = 0;
			weight = 1;
		}

	public :
		UInt_t   run;
		UInt_t   event;
		UInt_t   entry;
		Float_t  weight;

	ClassDef(LIST, 2)
};


// G4MC
class G4MC : public TObject {
	public :
		G4MC() { init(); }
		~G4MC() {}

		void init() {
			primPart.init();
			primVtx.init();
		}

	public :
		PartMCInfo   primPart;
		VertexMCInfo primVtx;

	ClassDef(G4MC, 7)
};


// RTI
class RTI : public TObject {
	public :
		RTI() { init(); }
		~RTI() {}

		void init() {
			flagRun = true;
			isGoodSecond = true;
			zenith = -1;
			std::fill_n(cfStormer, 4, -1);
			std::fill_n(cfIGRF, 4, -1);
			radiusGTOD = -1;
			thetaGTOD = -99;
			phiGTOD = -99;
		    latMAG = -99;
			longMAG = -99;
			latGAT = -99;
			longGAT = -99;
			std::fill_n(rptISS, 3, 0);
			std::fill_n(velISS, 3, 0);
			std::fill_n(yprISS, 3, 0);
			isInSAA = false;
			uTime = 0;
			dateUTC = 0;
			timeUTC = 0;
			liveTime = -1;
			std::fill_n(trackerAlign[0], 4, 0);
			trackerTemp = 0;
			isInShadow  = -1;
		}

	public :
		Bool_t  flagRun;             // true, if good
		Bool_t  isGoodSecond;        // true, if good
		Float_t zenith;              // ams zenith (z-axis) angle [degrees]
		Float_t cfStormer[4];        // max Stormer cutoff for 25, 30, 35, 40 degrees [GeV]
		Float_t cfIGRF[4];           // max IGRF cutoff for 25, 30, 35, 40 degrees [GeV]
		Float_t radiusGTOD;          // distance from earth to ams [cm]
		Float_t thetaGTOD;           // earth coordinate [rad]
		Float_t phiGTOD;             // earth coordinate [rad]
		Float_t latMAG;              // geomagnetic latitude [rad]
		Float_t longMAG;             // geomagnetic longitude [rad]
		Float_t latGAT;              // ams pointing galatic latitude [rad]
		Float_t longGAT;             // ams pointing galatic longitude [rad]
		Float_t rptISS[3];           // ISS coordinates (R, Phi, Theta) (GTOD)
		Float_t velISS[3];           // ISS velocity (Vel rad/sec, VelPhi rad, VelTheta rad)
		Float_t yprISS[3];           // ISS attitude (Yaw, Pitch, Roll)
		Bool_t  isInSAA;             // true, if ams in south atlantic anomaly
		UInt_t  uTime;               // unix time
		UInt_t  dateUTC;             // UTC date
		UInt_t  timeUTC;             // UTC time
		Float_t liveTime;            // fraction of "non-busy" time
		Float_t trackerAlign[2][2];  // L1 x,y L9 x,y
		Float_t trackerTemp;         // tracker temperature (Sensor A)
		
		Short_t isInShadow;          // particle pass through the ISS Solar Array
		                             // return
									 //  -1, no particle information
									 //   0, not in shadow
									 //   1, in shadow

	ClassDef(RTI, 8)
};


// TRG
class TRG : public TObject {
	public :
		TRG() { init(); }
		~TRG() {}

		void init() {
			bit = 0;
			physicalPatt = -1;
			logicPatt = -1;
		}

	public :
		Short_t bit;          // level 1 trigger : bit := extTrg * 1 + unBiasTrgTOF * 2 + unBiasTrgECAL * 4 + physTrg * 8;
		Int_t   physicalPatt; // level 1 trigger : physical
		Int_t   logicPatt;    // level 1 trigger : logic

	ClassDef(TRG, 2)
};


// TOF
class TOF : public TObject {
	public :
		TOF() { init(); }
		~TOF() {}

		void init() {
			numOfCls = 0;
			numOfClsH = 0;
			numOfBeta = 0;
			numOfBetaH = 0;

			numOfInTimeCls = -1;

			statusBeta = false;
			beta = 0;

			statusBetaH = false;
			betaHBit = 0;
			betaHPatt = 0;
			betaH = 0;
			normChisqT = -1;
			normChisqC = -1;
			std::fill_n(coo[0], 4*3, 0);
			std::fill_n(T, 4, -1);
			std::fill_n(Q, 4, -1);
			Qall = -1;
            Zall = -1;

			std::fill_n(numOfExtCls, 4, 0);
		}

	public :
		Short_t numOfCls;
		Short_t numOfClsH;
		Short_t numOfBeta;
		Short_t numOfBetaH;

		Short_t numOfInTimeCls;

		Bool_t  statusBeta;
		Float_t beta;

		Bool_t  statusBetaH;
		Short_t betaHBit;
		Short_t betaHPatt;
		Float_t betaH;
		Float_t normChisqT;
		Float_t normChisqC;
        Float_t coo[4][3];
		Float_t T[4];
		Float_t Q[4];
		Float_t Qall;
        Short_t Zall;

		// extern clusters
		Short_t numOfExtCls[4];

	ClassDef(TOF, 10)
};


// ACC
class ACC : public TObject {
	public :
		ACC() { init(); }
		~ACC() {}

		void init() {
			numOfCls = 0;
			clusters.clear();
		}

	public :
		Short_t numOfCls;

		std::vector<ClsACCInfo> clusters;

	ClassDef(ACC, 2)
};


// TRK
class TRK : public TObject {
	public :
		TRK() { init(); }
		~TRK() {}

		void init() {
            numOfTrack = 0;
			
            status = false;
            bitPatt = 0;
			bitPattXY = 0;
			
            QIn = -1;
			QL2 = -1;
			QL1 = -1;
			QL9 = -1;

            hits.clear();

            ckTr = std::vector<CKTrackInfo>(4);
            kfTr = std::vector<KFTrackInfo>(4);
            hcTr = std::vector<HCTrackInfo>(4);
            
            hcTrTF = std::vector<HCTrackInfo>(4);
		}

	public :
        Short_t numOfTrack;

		// bitPatt   := hasInner   * 1 + hasL2   * 2 + hasL1   * 4 + hasL9   * 8
		// bitPattXY := hasInnerXY * 1 + hasL2XY * 2 + hasL1XY * 4 + hasL9XY * 8
		// Inner    (XY) := ((bitPattXY&  1)==  1)
		// InnerL1  (XY) := ((bitPattXY&  5)==  5)
		// InnerL9  (XY) := ((bitPattXY&  9)==  9)
		// FullSpan (XY) := ((bitPattXY& 13)== 13)
        Bool_t  status;
		Short_t bitPatt;
		Short_t bitPattXY;

		// Track Charge
		Float_t QIn;
		Float_t QL1;
		Float_t QL2;
		Float_t QL9;
		
        // Track Hits
        std::vector<HitTRKInfo> hits;

        // Choutko [Inn InnL1 InnL9 FS]
        std::vector<CKTrackInfo> ckTr;

        // Choutko [Inn InnL1 InnL9 FS]
        std::vector<KFTrackInfo> kfTr;
        
        // HYChou [Inn InnL1 InnL9 FS]
        std::vector<HCTrackInfo> hcTr;
        
        std::vector<HCTrackInfo> hcTrTF;

	ClassDef(TRK, 9)
};


// TRD
class TRD : public TObject {
	public :
		TRD() { init(); }
		~TRD() {}

		void init() {
            numOfCls = 0;
            numOfSegment = 0;
			numOfTrack = 0;
			numOfHTrack = 0;
			
            trackStatus = false;
			std::fill_n(trackState, 6, 0);

			std::fill_n(statusKCls, 2, false);
			std::fill_n(LLRep, 2, -1);
			std::fill_n(LLReh, 2, -1);
			std::fill_n(LLRph, 2, -1);
			std::fill_n(LLRnhit, 2, -1);
			std::fill_n(Q, 2, -1);
            
            hits[0].clear();
            hits[1].clear();

            //std::fill_n(vtxNum, 3, 0);
            //vtxNTrk = 0;
            //vtxNHit = 0;
            //vtxChi2 = 0;
            //std::fill_n(vtxCoo, 3, 0);
		}

	public :
        Short_t numOfCls;
        Short_t numOfSegment;
		Short_t numOfTrack;
		Short_t numOfHTrack;
		
        // TrdTrack or TrdHTrack (first TrdH, second Trd)
		Bool_t  trackStatus;
		Float_t trackState[6]; // coo, dir

		// (TrdHTrack or TrdTrack) and TrTrack
		Bool_t  statusKCls[2]; // true, rebuild success (Trd, Trk)
		Float_t LLRep[2];
		Float_t LLReh[2];
		Float_t LLRph[2];
		Short_t LLRnhit[2];
		Float_t Q[2];
        
        std::vector<HitTRDInfo> hits[2];

        // TRDVertex
        //Short_t vtxNum[3]; // (3d, 2d_y, 2d_x)
        //Short_t vtxNTrk;
        //Short_t vtxNHit;
        //Float_t vtxChi2;
        //Float_t vtxCoo[3];

	ClassDef(TRD, 8)
};


// RICH
class RICH : public TObject {
	public :
		RICH() { init(); }
		~RICH() {}

		void init() {
            numOfHit = 0;
            
            status = false;
			isGood = false;
            kind = -1;
            tile = -1;
            refz =  0;
			beta = -1;
			Q = -1;
            numOfExpPE = -1;
            eftOfColPE = -1;
			
            vetoKind = -1;
			vetoTile = -1;
			vetoRfrIndex = -1;
			vetoDistToBorder = -1;
			vetoIsInFiducialVolume = false;

            //std::fill_n(vetoNumOfExpPE, 3, 0);
            //
            //vetoNumOfCrsHit = 0;
            //std::fill_n(vetoNumOfCkvHit[0], 3*3, 0);
            //
            //hits.clear();
		}

	public :
        Short_t numOfHit;
        
        // Official RichRingR
		Bool_t  status;
		Bool_t  isGood;
        Short_t kind;
        Short_t tile;
        Float_t refz;
		Float_t beta;
		Float_t Q;
        Float_t numOfExpPE;
        Float_t eftOfColPE;

		// Rich Veto
		Short_t vetoKind;          // -1, None, 0, Aerogel 1, NaF
		Short_t vetoTile;          // tile id
		Float_t vetoRfrIndex;      // refractive index
		Float_t vetoDistToBorder;  // dist To Border Edge
		Bool_t  vetoIsInFiducialVolume;
       
        // [0] electron [1] proton [2] deuterium
		//Float_t vetoNumOfExpPE[3];     // number of photoelectrons expected for a given track, beta and charge.
        //Short_t vetoNumOfCrsHit;       // Particle Cross PMT
        //Short_t vetoNumOfCkvHit[3][3]; // Photo
                                       // [0] In  Ring
                                       // [1] Out Ring
                                       // [2] Non Physics

        // Rich Hits
        //std::vector<HitRICHInfo> hits;
	
    ClassDef(RICH, 9)
};


// ECAL
class ECAL : public TObject {
	public :
		ECAL() { init(); }
		~ECAL() {}

		void init() {
			numOfShower = 0;
            shower.init();
		}

	public :
		Short_t    numOfShower;
		ShowerInfo shower;

	ClassDef(ECAL, 5)
};

#endif // __ClassDef_H__
