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
#include <string>
#include "TObject.h"

// RunTagInfo
class RunTagInfo : public TObject {
	public :
		RunTagInfo() { init(); }
		~RunTagInfo() {}

		void init() {
			runID = 0;
			eventFT = 0;
			eventLT = 0;
			numOfSelEvent = 0;
			numOfTrgEvent = 0;
			file.clear();
			
			UTCyr = 0;
			UTCyd = 0;
			UTCmo = 0;
			UTCmd = 0;
			UTChr = 0;
			UTCmn = 0;
			UTCsc = 0;
		}

	public :
		UInt_t runID;
		UInt_t eventFT; // first trigger event
		UInt_t eventLT; // last trigger event
		UInt_t numOfSelEvent; // number of select  events
		UInt_t numOfTrgEvent; // number of trigger events
		std::vector<std::string> file; // file list

		UInt_t UTCyr; // UTC year
		UInt_t UTCyd; // UTC year day
		UInt_t UTCmo; // UTC month
		UInt_t UTCmd; // UTC month day
		UInt_t UTChr; // UTC hour
		UInt_t UTCmn; // UTC min
		UInt_t UTCsc; // UTC second

	ClassDef(RunTagInfo, 4)
};


// HitTRKMCInfo
class HitTRKMCInfo : public TObject {
	public :
		HitTRKMCInfo() { init(); }
		~HitTRKMCInfo() {}

		void init() {
			layJ    = 0;
			tkid    = 0;
			edep    = 0;
			mom     = 0;
			std::fill_n(coo, 3, 0);
			std::fill_n(dir, 3, 0);
		}

	public :
		Short_t layJ;    // layerJ
		Short_t tkid;    // tkID
		Float_t edep;    // (elc) edep
		Float_t mom;     // (phy) momentum
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
			edep    = 0;
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
			trackID  = 0;
			parentID = 0;
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
		Short_t trackID;
		Short_t parentID;
		Short_t partID;
		Float_t chrg;
		Float_t mass;
		Float_t mom;
		Float_t kEng; // kinetic energy
		Float_t coo[3];
		Float_t dir[3];

		std::vector<HitTRKMCInfo> hits;
		//std::vector<HitTRDMCInfo> hits;

	ClassDef(PartMCInfo, 3)
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
			type = 0;
			std::fill_n(coo, 3, 0.0);
			partID.clear();
		}

	public :
		Int_t              type;
		Float_t            coo[3];
		std::vector<Int_t> partID;
	
	ClassDef(VertexMCInfo, 2)
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
			std::fill_n(cofg, 2, -1);
			std::fill_n(chrg, 2, -1);
			std::fill_n(coo,  3, 0);
		}

	public :
		Short_t clsId[2]; // cluster id (x, y)
		Short_t layJ;     // layerJ
		Short_t tkid;     // tkID
		Short_t sens;     // sensor
		Short_t mult;     // multiplicity
		Short_t side;     // side, 1 x, 2 y, 3 xy
		Float_t cofg[2];  // (elc) strip index
		Float_t chrg[2];  // (elc) charge
		Float_t coo[3];

	ClassDef(HitTRKInfo, 4)
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
			amp  = 0;
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


// HitECALInfo
class HitECALInfo : public TObject {
	public :
		HitECALInfo() { init(); }
		~HitECALInfo() {}

		void init() {
			id   = 0;
			side = 0;
			edep = 0;
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


// TrackInfo
class TrackInfo : public TObject {
	public :
		TrackInfo() { init(); }
		~TrackInfo() {}

		void init() {
			bitPattJ = 0;
			bitPattXYJ = 0;

			bitPatt = 0;
			Qinner = -1;
			
			std::fill_n(status[0], 3 * 6, false);
			std::fill_n(rigidity[0], 3 * 6, 0);
			std::fill_n(chisq[0][0], 3 * 6 * 2, -1);
			std::fill_n(state[0][0], 3 * 6 * 6, 0);
			std::fill_n(stateL1[0][0], 3 * 6 * 6, 0);
			std::fill_n(stateL9[0][0], 3 * 6 * 6, 0);
			
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
		Float_t  Qinner;

		// Algorithm     (CHOUTKO, ALCARAZ, CHIKANIANF)
		// Track Pattern (InnU, InnL, Inn, InnL1, InnL9, FS)
		Bool_t  status[3][6];
		Float_t rigidity[3][6];
		Float_t chisq[3][6][2];
		Float_t state[3][6][6];
		Float_t stateL1[3][6][6];
		Float_t stateL9[3][6][6];
		
		std::vector<HitTRKInfo> hits;

	ClassDef(TrackInfo, 4)
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
			runID   = 0;
			eventID = 0;
			entryID = 0;
			weight  = 1;
		}

	public :
		UInt_t   runID;
		UInt_t   eventID;
		UInt_t   entryID;
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

	ClassDef(G4MC, 5)
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
			xTime = 0;
			UTCyr = 0;
			UTCyd = 0;
			UTCmo = 0;
			UTCmd = 0;
			UTChr = 0;
			UTCmn = 0;
			UTCsc = 0;
			liveTime = -1;
			std::fill_n(trackerAlign[0], 4, 0);
			isInShadow  = -1;
			isFromSpace = -1;
			std::fill_n(backtrace[0], 2*3, -1);
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
		Float_t xTime;               // unix time + frac
		UInt_t  UTCyr;               // UTC year
		UInt_t  UTCyd;               // UTC year day
		UInt_t  UTCmo;               // UTC month
		UInt_t  UTCmd;               // UTC month day
		UInt_t  UTChr;               // UTC hour
		UInt_t  UTCmn;               // UTC min
		UInt_t  UTCsc;               // UTC second
		Float_t liveTime;            // fraction of "non-busy" time
		Float_t trackerAlign[2][2];  // L1 x,y L9 x,y
		
		Short_t isInShadow;          // particle pass through the ISS Solar Array
		                             // return
																 //  -1, no particle information
																 //   0, false
																 //   1, true
		Short_t isFromSpace;         // particle coming from space
		                             // return
																 //  -1, no particle information
																 //   0, particle isnot coming from space
																 //   1, particle is coming from space (weak)
																 //   2, particle is coming from space (middle)
																 //   3, particle is coming from space (strong)
		Short_t backtrace[2][3];     // charge { 1, -1 }  (stable fact 1.00, 1.15, 1.30)
		                             // return
		                             //   0, unercutoff (i.e. atmospheric origin), 
													       //   1, over cutoff (i.e. coming from space), 
													       //   2, trapped, 
													       //  -1, error

	ClassDef(RTI, 5)
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
			numOfRawCluster = 0;
			numOfRawSide = 0;
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
			std::fill_n(Q, 4, 0);
			Qall = -1;

			statusBetaHs = false;
			betaHBits = 0;
			betaHPatts = 0;
			betaHGoodTimes = 0;
			betaHs = 0;
			normChisqTs = -1;
			normChisqCs = -1;
			std::fill_n(Qs, 4, 0);
			Qalls = -1;
			std::fill_n(betaHStates, 6, 0);
		}

	public :
		Short_t numOfRawCluster;
		Short_t numOfRawSide;
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
		Float_t Q[4];
		Float_t Qall;

		// TRK Track independent
		Bool_t  statusBetaHs;
		Short_t betaHBits;
		Short_t betaHPatts;
		Short_t betaHGoodTimes;
		Float_t betaHs;
		Float_t normChisqTs;
		Float_t normChisqCs;
		Float_t Qs[4];
		Float_t Qalls;
		Float_t betaHStates[6];

	ClassDef(TOF, 2)
};


// ACC
class ACC : public TObject {
	public :
		ACC() { init(); }
		~ACC() {}

		void init() {
			numOfCluster = 0;
		}

	public :
		Short_t numOfCluster;

	ClassDef(ACC, 1)
};


// TRK
class TRK : public TObject {
	public :
		TRK() { init(); }
		~TRK() {}

		void init() {
			beamID = -1;
			beamDist = -1;

			numOfTrack = 0;

			tracks.clear();
			vertices.clear();

			for (UInt_t it = 0; it < 4; ++it)
				maxQHit[it].init();
		}

	public :
		Short_t beamID;
		Float_t beamDist;

		Short_t  numOfTrack;
		std::vector<TrackInfo> tracks;
		std::vector<VertexInfo> vertices;

		HitTRKInfo maxQHit[4]; // L1, L2, L3-8, L9

	ClassDef(TRK, 3)
};


// TRD
class TRD : public TObject {
	public :
		TRD() { init(); }
		~TRD() {}

		void init() {
			numOfRawHit = 0;
			numOfCluster = 0;
			numOfSegment = 0;
			numOfTrack = 0;
			numOfHSegment = 0;
			numOfHTrack = 0;

			numOfVertexWithTrTrack = -1;

			std::fill_n(statusKCls, 2, false);
			std::fill_n(Q, 2, -1);
			std::fill_n(LLR[0], 6, -1);
			std::fill_n(LLR_nhit, 2, 0);

			trackStatus = false;
			std::fill_n(trackState, 6, 0);
		}

	public :
		Short_t numOfRawHit;
		Short_t numOfCluster;
		Short_t numOfSegment;
		Short_t numOfTrack;
		Short_t numOfHSegment;
		Short_t numOfHTrack;

		Short_t numOfVertexWithTrTrack;

		// (TrdHTrack or TrdTrack) and TrTrack
		Bool_t  statusKCls[2]; // true, rebuild success (Trd, Trk)
		Float_t Q[2];
		Float_t LLR[2][3]; // eP eH PH
		Short_t LLR_nhit[2];

		// TrdTrack or TrdHTrack (first TrdH, second Trd)
		Bool_t  trackStatus;
		Float_t trackState[6]; // coo, dir

	ClassDef(TRD, 3)
};


// RICH
class RICH : public TObject {
	public :
		RICH() { init(); }
		~RICH() {}

		void init() {
			numOfRing = 0;
			numOfHit = 0;

			kindOfRad = -1;
			tileOfRad = -1;
			rfrIndex  = -1;
			std::fill_n(emission, 6, 0);
			distToBorder = -1;
			std::fill_n(numOfExpPE, 5, -1);
			isGoodTile = false;
			isInFiducialVolume = false;
			
			status = false;
			isGoodRecon = false;
			beta = -1;
			Q = -1;
		}

	public :
		// official RichRingR
		Short_t numOfRing;
		Short_t numOfHit;

		Short_t kindOfRad;     // -1, None, 0, Aerogel 1, NaF
		Short_t tileOfRad;     // tile id
		Float_t rfrIndex;      // refractive index
		Float_t emission[6];   //  
		Float_t distToBorder;  // dist To Border Edge
		Float_t numOfExpPE[5]; // number of photoelectrons expected for a given track, beta and charge.
		                       // [0] electron
													 // [1] pion
													 // [2] kaon
													 // [3] proton
													 // [4] deuterium
		Bool_t isGoodTile;
		Bool_t isInFiducialVolume;

		Bool_t  status;
		Bool_t  isGoodRecon;
		Float_t beta;
		Float_t Q;

	ClassDef(RICH, 3)
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

	ClassDef(ECAL, 3)
};

#endif // __ClassDef_H__
