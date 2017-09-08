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
			//file.clear();
			
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


// HitTRKMCInfo
class HitTRKMCInfo : public TObject {
	public :
		HitTRKMCInfo() { init(); }
		~HitTRKMCInfo() {}

		void init() {
			layJ    = 0;
			tkid    = 0;
			edep    = -1;
			mom     = -1;
			std::fill_n(coo, 3, 0);
			std::fill_n(dir, 3, 0);
		}

	public :
		Short_t layJ;    // layerJ
		Short_t tkid;    // tkID
		Float_t edep;    // edep
		Float_t mom;     // momentum
		Float_t coo[3];
		Float_t dir[3];

	ClassDef(HitTRKMCInfo, 3)
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


// HitTRDMCInfo
class HitTRDMCInfo : public TObject {
	public :
		HitTRDMCInfo() { init(); }
		~HitTRDMCInfo() {}

		void init() {
			lay     = 0;
			sub     = 0;
			edep    = -1;
			std::fill_n(coo, 3, 0);
		}

	public :
		Short_t lay;     // layer
		Short_t sub;     // Ladder * 100 + Tube
		Float_t edep;    // (elc) edep
		Float_t coo[3];

	ClassDef(HitTRDMCInfo, 3)
};

struct HitTRDMCInfo_sort {
	bool operator() (const HitTRDMCInfo & hit1, const HitTRDMCInfo & hit2) {
		if      (hit1.lay < hit2.lay) return true;
		else if (hit1.lay > hit2.lay) return false;
		else {
			if      (hit1.sub < hit2.sub) return true;
			else if (hit1.sub > hit2.sub) return false;
			else {
				if      (hit1.edep < hit2.edep) return true;
				else if (hit1.edep > hit2.edep) return false;
			}
		}
		return false;
	}
};


// PartMCInfo
class PartMCInfo : public TObject {
	public :
		PartMCInfo() { init(); }
		~PartMCInfo() {}

		void init() {
			partID   = 0;
			chrg     = 0;
			mass     = 0;
			mom      = 0;
			kEng     = 0;
			std::fill_n(coo, 3, 0);
			std::fill_n(dir, 3, 0);
			hits.clear();
		}

	public :
		Short_t partID;
		Float_t chrg;
		Float_t mass;
		Float_t mom;
		Float_t kEng; // kinetic energy
		Float_t coo[3];
		Float_t dir[3];

		std::vector<HitTRKMCInfo> hits;

	ClassDef(PartMCInfo, 5)
};

struct PartMCInfo_sort {
	bool operator() (const PartMCInfo & par1, const PartMCInfo & par2) {
		if      (par1.coo[2] > par2.coo[2]) return true;
		else if (par1.coo[2] < par2.coo[2]) return false;
		else {
			if      (par1.partID < par2.partID) return true;
			else if (par1.partID > par2.partID) return false;
			else {
				if      (par1.mom > par2.mom) return true;
				else if (par1.mom < par2.mom) return false;
			}
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
			std::fill_n(coo, 3, 0.0);
		}

	public :
		Bool_t             status;
		Float_t            coo[3];
	
	ClassDef(VertexMCInfo, 3)
};


// HitTRKInfo
class HitTRKInfo : public TObject {
	public :
		HitTRKInfo() { init(); }
		~HitTRKInfo() {}

		void init() {
			std::fill_n(clsId,  2, -1);
			layJ =  0;
			tkid =  0;
			sens = -1;
			mult = -1;
			side =  0;
			std::fill_n(coo,  3, 0);
			std::fill_n(chrg, 2, -1);

			cofgX = 0;
			seedAddrX = -1;
			seedIndxX = -1;
			stripSigX.clear();
			stripSgmX.clear();

			cofgY = 0;
			seedAddrY = -1;
			seedIndxY = -1;
			stripSigY.clear();
			stripSgmY.clear();
		}

	public :
		Short_t clsId[2]; // cluster id (x, y)
		Short_t layJ;     // layerJ
		Short_t tkid;     // tkID
		Short_t sens;     // sensor
		Short_t mult;     // multiplicity
		Short_t side;     // side, 1 x, 2 y, 3 xy
		Float_t coo[3];   // coordinate
		Float_t chrg[2];  // charge
		
		Float_t              cofgX;     // (elc) cofg	
		Float_t              seedAddrX; // (elc) seed strip address
		Float_t              seedIndxX; // (elc) seed strip index
		std::vector<Float_t> stripSigX; // (elc) strip signal (value)
		std::vector<Float_t> stripSgmX; // (elc) strip signal (sigma)
		
		Float_t              cofgY;     // (elc) cofg	
		Float_t              seedAddrY; // (elc) seed strip address
		Float_t              seedIndxY; // (elc) seed strip index
		std::vector<Float_t> stripSigY; // (elc) strip signal (value)
		std::vector<Float_t> stripSgmY; // (elc) strip signal (sigma)

	ClassDef(HitTRKInfo, 6)
};

struct HitTRKInfo_sort {
	bool operator() (const HitTRKInfo & hit1, const HitTRKInfo & hit2) {
		if      (hit1.layJ < hit2.layJ) return true;
		else if (hit1.layJ > hit2.layJ) return false;
		else {
			if      (hit1.tkid < hit2.tkid) return true;
			else if (hit1.tkid > hit2.tkid) return false;
			else {
				if      (hit1.side < hit2.side) return true;
				else if (hit1.side > hit2.side) return false;
			}
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
			lay  = 0;
			sub  = 0;
			side = 0;
			amp  = -1;
			std::fill_n(coo, 3, 0);
		}

	public :
		Short_t lay;    // layer
		Short_t sub;    // ladder * 100 + tube
		Short_t side;   // side, 1 x, 2 y, 3 xy
		Float_t amp;    // (elc) amp
		Float_t coo[3];

	ClassDef(HitTRDInfo, 2)
};

struct HitTRDInfo_sort {
	bool operator() (const HitTRDInfo & hit1, const HitTRDInfo & hit2) {
		if      (hit1.lay < hit2.lay) return true;
		else if (hit1.lay > hit2.lay) return false;
		else {
			if      (hit1.sub < hit2.sub) return true;
			else if (hit1.sub > hit2.sub) return false;
			else {
				if      (hit1.amp < hit2.amp) return true;
				else if (hit1.amp > hit2.amp) return false;
			}
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


// TrackInfo
class TrackInfo : public TObject {
	public :
		TrackInfo() { init(); }
		~TrackInfo() {}

		void init() {
			bitPattJ = 0;
			bitPattXYJ = 0;
			bitPatt = 0;
			QIn = -1;
			QL2 = -1;
			QL1 = -1;
			QL9 = -1;
			
			std::fill_n(status[0], 2 * 4, false);
			std::fill_n(rigidity[0], 2 * 4, 0);
			std::fill_n(chisq[0][0], 2 * 4 * 2, -1);
			std::fill_n(state[0][0], 2 * 4 * 6, 0);
			
			std::fill_n(stateLJ[0][0][0], 2 * 4 * 9 * 6, 0);
			//std::fill_n(localID[0][0][0], 2 * 4 * 9 * 3, 0);
			//std::fill_n(localLJ[0][0][0], 2 * 4 * 9 * 2, -1);
			
			hits.clear();
		}

	public :
		UShort_t bitPattJ;
		UShort_t bitPattXYJ;
		
		// bitPatt := hasInner *   1 + hasInnerXY *   2 +
		//            hasL2    *   4 + hasL2XY    *   8 +
		//            hasL1    *  16 + hasL1XY    *  32 +
		//            hasL9    *  64 + hasL9XY    * 128
		// Inner    (XY) := ((bitPatt&   2)==  2)
		// InnerL1  (XY) := ((bitPatt&  34)== 34)
		// InnerL9  (XY) := ((bitPatt& 130)==130)
		// FullSpan (XY) := ((bitPatt& 162)==162)
		UShort_t bitPatt;

		// Track Charge
		Float_t QIn;
		Float_t QL1;
		Float_t QL2;
		Float_t QL9;

		// Algorithm     (CHOUTKO, CHIKANIANF)
		// Track Pattern (Inn, InnL1, InnL9, FS)
		Bool_t  status[2][4];
		Float_t rigidity[2][4];
		Float_t chisq[2][4][2];
		Float_t state[2][4][6];
		
		Float_t stateLJ[2][4][9][6]; // (x y z dirx diry dirz)
		//Short_t localID[2][4][9][3]; // (tkid sens mult)
		//Float_t localLJ[2][4][9][2]; // (xloc yloc)
	
		// Track Hits
		std::vector<HitTRKInfo> hits;

	ClassDef(TrackInfo, 6)
};


// VertexInfo
class VertexInfo : public TObject {
	public :
		VertexInfo() { init(); }
		~VertexInfo() {}

		void init() {
			std::fill_n(vtxZqu, 5, 0.0);
			std::fill_n(vtxR, 8, 0.0);
			vtxR[7] = -1.0;
			trackID.clear();
		}

	public :
		Float_t vtxZqu[5]; // vertex (coo, dx, dy) based on Zqu
		Float_t vtxR[8]; // vertex (coo, dir, mom, chi/ndf) based on VertexR
		std::vector<Int_t> trackID;
	
	ClassDef(VertexInfo, 2)
};


// ShowerInfo
class ShowerInfo : public TObject {
	public :
		ShowerInfo() { init(); }
		~ShowerInfo() {}
		
		void init() {
			energyD = -1;
			energyE = -1;
			energyP = -1;
			PisaBDT = -2;
			Q = -1;
			std::fill_n(showerAxis, 6, 0);

			hadronApex = -1;
			hadronEnergy = -1;
		}

	public :
		Float_t energyD;  // energy deposit [GeV]
		Float_t energyE;  // Pisa
		Float_t energyP;  // Choutko - base on energyC
		Float_t PisaBDT;  // Pisa
		Float_t Q;        // ecal charge
		Float_t showerAxis[6];

		Short_t hadronApex; // reject high-apex events (say, Apex > 12)
		Float_t hadronEnergy;

		ClassDef(ShowerInfo, 2)
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
			beamID = -1;

			primPart.init();
			secParts.clear();
			primVtx.init();
		}

	public :
		Short_t beamID; // only for MC Beam Test (400GeV proton)
		PartMCInfo                 primPart;
		std::vector<PartMCInfo>    secParts;
		VertexMCInfo               primVtx;

	ClassDef(G4MC, 6)
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
			std::fill_n(cutoffStormer, 4, -1);
			std::fill_n(cutoffIGRF, 4, -1);
			radiusGTOD = -1;
			thetaGTOD = -99;
			phiGTOD = -99;
			thetaMAG = -99;
			phiMAG = -99;
			latGXY = -999;
			longGXY = -999;
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
			//isFromSpace = -1;
			//std::fill_n(backtrace[0], 2*3, -1);
		}

	public :
		Bool_t  flagRun;             // true, if good
		Bool_t  isGoodSecond;        // true, if good
		Float_t zenith;              // ams zenith (z-axis) angle [degrees]
		Float_t cutoffStormer[4];    // max Stormer cutoff for 25, 30, 35, 40 degrees [GeV]
		Float_t cutoffIGRF[4];       // max IGRF cutoff for 25, 30, 35, 40 degrees [GeV]
		Float_t radiusGTOD;          // distance from earth to ams [cm]
		Float_t thetaGTOD;           // earth coordinate [rad]
		Float_t phiGTOD;             // earth coordinate [rad]
		Float_t thetaMAG;            // geomagnetic latitude [rad]
		Float_t phiMAG;              // geomagnetic longitude [rad]
		Float_t latGXY;              // ams pointing galatic latitude [rad]
		Float_t longGXY;             // ams pointing galatic longitude [rad]
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
		//Short_t isFromSpace;     // particle coming from space
		                             // return
									 //  -1, no particle information
									 //   0, particle isnot coming from space
									 //   1, particle is coming from space (weak)
									 //   2, particle is coming from space (middle)
									 //   3, particle is coming from space (strong)
		//Short_t backtrace[2][3];   // charge { 1, -1 }  (stable fact 1.00, 1.15, 1.30)
		                             // return
		                             //   0, unercutoff (i.e. atmospheric origin), 
							       //   1, over cutoff (i.e. coming from space), 
							       //   2, trapped, 
							       //  -1, error

	ClassDef(RTI, 7)
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
			numOfCluster = 0;
			numOfClusterH = 0;
			numOfBeta = 0;
			numOfBetaH = 0;

			numOfInTimeCluster = -1;

			statusBeta = false;
			beta = 0;

			statusBetaH = false;
			betaHBit = 0;
			betaHPatt = 0;
			betaHGoodTime = 0;
			betaH = 0;
			normChisqT = -1;
			normChisqC = -1;
			std::fill_n(T, 4,  0);
			std::fill_n(Q, 4, -1);
			Qall = -1;

			std::fill_n(extClsN, 2, 0);
			std::fill_n(extClsL[0], 2*2, -1);
			std::fill_n(extClsQ[0], 2*2, -1);
			std::fill_n(extClsT[0], 2*2,  0);

			//statusBetaHs = false;
			//betaHBits = 0;
			//betaHPatts = 0;
			//betaHGoodTimes = 0;
			//betaHs = 0;
			//normChisqTs = -1;
			//normChisqCs = -1;
			//std::fill_n(Ts, 4,  0);
			//std::fill_n(Qs, 4, -1);
			//Qalls = -1;
			//std::fill_n(betaHStates, 6, 0);
		}

	public :
		Short_t numOfCluster;
		Short_t numOfClusterH;
		Short_t numOfBeta;
		Short_t numOfBetaH;

		Short_t numOfInTimeCluster;

		Bool_t  statusBeta;
		Float_t beta;

		Bool_t  statusBetaH;
		Short_t betaHBit;
		Short_t betaHPatt;
		Short_t betaHGoodTime;
		Float_t betaH;
		Float_t normChisqT;
		Float_t normChisqC;
		Float_t T[4];
		Float_t Q[4];
		Float_t Qall;

		// extern clusters (near betaH hit by time)
		Short_t extClsN[2];
		Short_t extClsL[2][2];
		Float_t	extClsQ[2][2];
		Float_t extClsT[2][2];

		// TRK Track independent
		//Bool_t  statusBetaHs;
		//Short_t betaHBits;
		//Short_t betaHPatts;
		//Short_t betaHGoodTimes;
		//Float_t betaHs;
		//Float_t normChisqTs;
		//Float_t normChisqCs;
		//Float_t Ts[4];
		//Float_t Qs[4];
		//Float_t Qalls;
		//Float_t betaHStates[6];

	ClassDef(TOF, 4)
};


// ACC
class ACC : public TObject {
	public :
		ACC() { init(); }
		~ACC() {}

		void init() {
			numOfCluster = 0;
			clusters.clear();
		}

	public :
		Short_t numOfCluster;

		std::vector<ClsACCInfo> clusters;

	ClassDef(ACC, 2)
};


// TRK
class TRK : public TObject {
	public :
		TRK() { init(); }
		~TRK() {}

		void init() {
			beamID = -1;
			beamDist = -1;

			tracks.clear();
			vertices.clear();

			otherHits.clear();
		}

	public :
		Short_t beamID;
		Float_t beamDist;

		std::vector<TrackInfo> tracks;
		std::vector<VertexInfo> vertices;

		std::vector<HitTRKInfo> otherHits;

	ClassDef(TRK, 4)
};


// TRD
class TRD : public TObject {
	public :
		TRD() { init(); }
		~TRD() {}

		void init() {
			numOfTrack = 0;
			numOfHTrack = 0;
			std::fill_n(numOfHSeg, 2, 0);
			
			vtxSide = 0;
			std::fill_n(vtxCoo, 3, 0);
			std::fill_n(vtxTrCoo, 3, 0);

			std::fill_n(statusKCls, 2, false);
			std::fill_n(Q, 2, -1);
			std::fill_n(LLR[0], 6, -1);
			std::fill_n(LLR_nhit, 2, 0);

			trackStatus = false;
			std::fill_n(trackState, 6, 0);
		}

	public :
		Short_t numOfTrack;
		Short_t numOfHTrack;
		Short_t numOfHSeg[2];

		// (TrTrack, TrdSegment) Vertex
		Short_t vtxSide;
		Float_t vtxCoo[3];
		Float_t vtxTrCoo[3];

		// (TrdHTrack or TrdTrack) and TrTrack
		Bool_t  statusKCls[2]; // true, rebuild success (Trd, Trk)
		Float_t Q[2];
		Float_t LLR[2][3]; // eP eH PH
		Short_t LLR_nhit[2];

		// TrdTrack or TrdHTrack (first TrdH, second Trd)
		Bool_t  trackStatus;
		Float_t trackState[6]; // coo, dir

	ClassDef(TRD, 5)
};


// RICH
class RICH : public TObject {
	public :
		RICH() { init(); }
		~RICH() {}

		void init() {
			kindOfRad  = -1;
			tileOfRad  = -1;
			isGoodTile = false;
			rfrIndex   = -1;
			distToBorder = -1;
			isInFiducialVolume = false;
			std::fill_n(emission, 6, 0);
			std::fill_n(receiving, 6, 0);

			std::fill_n(numOfExpPE, 5, -1);
			std::fill_n(theta, 5, -1);

            hits.clear();
		
			std::fill_n(numOfCrossHit, 2, 0);
			std::fill_n(numOfRingHit[0], 5*3, 0);

			status = false;
			isGoodRecon = false;
			beta = -1;
			Q = -1;
		}

	public :
		// Rich Veto
		Short_t kindOfRad;     // -1, None, 0, Aerogel 1, NaF
		Short_t tileOfRad;     // tile id
		Bool_t isGoodTile;
		Float_t rfrIndex;      // refractive index
		Float_t distToBorder;  // dist To Border Edge
		Bool_t isInFiducialVolume;
		Float_t emission[6];
        Float_t receiving[6];
		
        // [0] electron [1] pion [2] kaon [3] proton [4] deuterium
		Float_t numOfExpPE[5]; // number of photoelectrons expected for a given track, beta and charge.
        Float_t theta[5];      // theta for a given track, beta and charge.

		// Rich Hits
        std::vector<HitRICHInfo> hits;

		Short_t numOfCrossHit[2];   // CrossHit[selected, others]
		Short_t numOfRingHit[5][3]; // RingHit[particle][selected, inside, outside]
                                    // [0] electron
                                    // [1] pion
                                    // [2] kaon
                                    // [3] proton
                                    // [4] deuterium

		// Official RichRingR
		Bool_t  status;
		Bool_t  isGoodRecon;
		Float_t beta;
		Float_t Q;

	ClassDef(RICH, 6)
};


// ECAL
class ECAL : public TObject {
	public :
		ECAL() { init(); }
		~ECAL() {}

		void init() {
			numOfShower = 0;

			showers.clear();

			//rawHits.clear();
		}

	public :
		Short_t numOfShower;
		
		std::vector<ShowerInfo> showers;

		//std::vector<HitECALInfo> rawHits;

	ClassDef(ECAL, 4)
};

#endif // __ClassDef_H__
