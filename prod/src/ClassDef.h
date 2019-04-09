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
			std::fill_n(chrg, 3, -1);
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
		Float_t chrg[3];  // (elc) chrg (x y xy)

	ClassDef(HitTRKInfo, 11)
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
            cr = 0;
            cz = 0;

            dEdX  = -1;
            mcMom = -1;
		}

	public :
		Short_t lay;    // layer
		Short_t side;   // side, 1 x, 2 y
		Float_t amp;    // (elc) amp
		Float_t len;    // (elc) len
		Float_t cr;
		Float_t cz;

        Float_t dEdX;   // dEdX = 0.01 * (amp/len)
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


class HitTRDSimpleInfo : public TObject {
	public :
		HitTRDSimpleInfo() { init(); }
		~HitTRDSimpleInfo() {}

		void init() {
            cz  = 0;
			amp = 0;
			len = 0;

            dEdX  = -1;
            mcMom = -1;
        }

    public :
		Float_t cz;
		Float_t amp;    // (elc) amp
		Float_t len;    // (elc) len
        
        Float_t dEdX;   // dEdX = 0.01 * (amp/len)
        Float_t mcMom;  // (MC Info) mom

	ClassDef(HitTRDSimpleInfo, 1)
};


// HitRICHInfo
class HitRICHInfo : public TObject {
	public :
		HitRICHInfo() { init(); }
		~HitRICHInfo() {}

		void init() {
            channel = 0;
            type = 0;
            bta[0] = 0.;
            bta[1] = 0.;
            npe = 0.;
            cx = 0;
            cy = 0;
            dist = -1;
		}

	public :
        Short_t channel;
        Short_t type; // 0, direct   1, reflected
        Float_t bta[2];
        Float_t npe; // ADC counts above the pedestal/gain of the channel.
		Float_t cx;
        Float_t cy;
        Float_t dist;

	ClassDef(HitRICHInfo, 3)
};

struct HitRICHInfo_sort {
	bool operator() (const HitRICHInfo & hit1, const HitRICHInfo & hit2) {
		if      (hit1.channel < hit2.channel) return true;
		else if (hit1.channel > hit2.channel) return false;
		return false;
	}
};


class ChHitInfo : public TObject {
    public :
        ChHitInfo() { init(); }
        ~ChHitInfo() {}
        
        void init() {
            chann = -1;
            pmtid = -1;

            cls  = -1;
            mode = -1;
            beta = -1;
            npe  = -1;
            cx   =  0;
            cy   =  0;
        }

    public :
        Short_t chann;
        Short_t pmtid;

        Short_t cls;

        Short_t mode;
        Float_t beta;
        
        Float_t npe;

        Float_t cx;
        Float_t cy;
    
    ClassDef(ChHitInfo, 1)
};


class ChStoneInfo : public TObject {
    public :
        ChStoneInfo() { init(); }
        ~ChStoneInfo() {}
        
        void init() {
            status = false;
            nhit   = 0;
            npmt   = 0;
            cx     = 0;
            cy     = 0;
            dist   = 0;
            npe    = 0;
            cnt    = 0;
            nchi   = 0;
            chic   = 0;
        }

    public :
        Bool_t  status;
        Short_t nhit;
        Short_t npmt;
        
        Float_t cx;
        Float_t cy;
        Float_t dist;
        Float_t npe;

        Float_t cnt;
        Float_t nchi;
        Float_t chic;

    ClassDef(ChStoneInfo, 1)
};


class ChCloudInfo : public TObject {
    public :
        ChCloudInfo() { init(); }
        ~ChCloudInfo() {}
		
        void init() {
            status   = false;
            nhit     = 0;
            npmt     = 0;

            nhit_dir = 0;
            nhit_rfl = 0;

            border   = 0;
            trace    = 0;
            accuracy = 1;
            uniform  = 1;
            crrch    = 0;
            
            beta     = 0;
            cbta     = 0;
            npe      = 0;
            expnpe   = 0;

            cnt      = 0;
            nchi     = 0;

            clcrad   = 0;
            misjudge = 0;
        }

    public :
        Bool_t  status;
        Short_t nhit;
        Short_t npmt;
        Short_t nhit_dir;
        Short_t nhit_rfl;

        Float_t border;
        Float_t trace;
        Float_t accuracy;
        Float_t uniform;
        Float_t crrch;
        
        Float_t beta;
        Float_t cbta;
        Float_t npe;
        Float_t expnpe;

        Float_t cnt;
        Float_t nchi;

        Float_t clcrad;
        Float_t misjudge;

    ClassDef(ChCloudInfo, 1)
};


class ChTumorInfo : public TObject {
    public :
        ChTumorInfo() { init(); }
        ~ChTumorInfo() {}
		
        void init() {
            status  = false;
            nhit    = 0;
            npmt    = 0;
            beta    = 0;
            cbta    = 0;
            npe     = 0;
        }

    public :
        Bool_t  status;
        Short_t nhit;
        Short_t npmt;
        Float_t beta;
        Float_t cbta;
        Float_t npe;

    ClassDef(ChTumorInfo, 1)
};


class ChGhostInfo : public TObject {
    public :
        ChGhostInfo() { init(); }
        ~ChGhostInfo() {}
		
        void init() {
            status  = false;
            nhit    = 0;
            npmt    = 0;
            beta    = 0;
            cbta    = 0;
            npe     = 0;
        }

    public :
        Bool_t  status;
        Short_t nhit;
        Short_t npmt;
        Float_t beta;
        Float_t cbta;
        Float_t npe;

    ClassDef(ChGhostInfo, 1)
};


class ChFitInfo : public TObject {
    public :
        ChFitInfo() { init(); }
        ~ChFitInfo() {}

        void init() {
            status = false;
            kind = 0;
            tile = 0;

            index = 0;
            dist = 0;

            is_good_geom = false;
            is_bad_tile = false;

            is_in_pmt_plane = false;

            std::fill_n(radp, 0., 3);
            std::fill_n(radd, 0., 3);
            std::fill_n(pmtp, 0., 3);

            nstn = 0;
            ncld = 0;
            ntmr = 0;
            ngst = 0;
            
            stone.init();
            cloud.init();
            tumor.init();
            ghost.init();
        
            hits.clear();

            nhit_ttl = 0;
            nhit_stn = 0;
            nhit_cld = 0;
            nhit_tmr = 0;
            nhit_gst = 0;
            nhit_emy = 0;
            nhit_oth = 0;
            nhit_oth_inn = 0;
            nhit_oth_out = 0;
        
            npe_ttl = 0;
            npe_stn = 0;
            npe_cld = 0;
            npe_tmr = 0;
            npe_gst = 0;
            npe_emy = 0;
            npe_oth = 0;
            npe_oth_inn = 0;
            npe_oth_out = 0;

            bta_crr = 1;

            cpuTime = 0;
        }

    public :
        Bool_t  status;
        Short_t kind;
        Short_t tile;

        Float_t index;
        Float_t dist;

        Bool_t  is_good_geom;
        Bool_t  is_bad_tile;

        Bool_t  is_in_pmt_plane;

        Float_t radp[3];
        Float_t radd[3];
        Float_t pmtp[3];

        Short_t nstn;
        Short_t ncld;
        Short_t ntmr;
        Short_t ngst;

        ChStoneInfo stone;
        ChCloudInfo cloud;
        ChTumorInfo tumor;
        ChGhostInfo ghost;

        std::vector<ChHitInfo> hits;

        Short_t nhit_ttl;
        Short_t nhit_stn;
        Short_t nhit_cld;
        Short_t nhit_tmr;
        Short_t nhit_gst;
        Short_t nhit_emy;
        Short_t nhit_oth;
        Short_t nhit_oth_inn;
        Short_t nhit_oth_out;

        Float_t npe_ttl;
        Float_t npe_stn;
        Float_t npe_cld;
        Float_t npe_tmr;
        Float_t npe_gst;
        Float_t npe_emy;
        Float_t npe_oth;
        Float_t npe_oth_inn;
        Float_t npe_oth_out;

        Float_t bta_crr;

        Float_t cpuTime;

    ClassDef(ChFitInfo, 1)
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
            
            statusCen = false;
            std::fill_n(stateCen, 6, 0);

            statusBtm = false;
            std::fill_n(stateBtm, 6, 0);
            
            statusTd = false;
            std::fill_n(stateTd, 6, 0);
            
            statusRh = false;
            std::fill_n(stateRh, 6, 0);
            
            //std::fill_n(statusLJ, 9, false);
            //std::fill_n(stateLJ[0], 9*6, 0);

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
        
        Bool_t  statusCen;
        Float_t stateCen[6]; // (cx cy cz ux uy uz)

        Bool_t  statusBtm;
        Float_t stateBtm[6]; // (cx cy cz ux uy uz)

        Bool_t  statusTd;
        Float_t stateTd[6]; // (cx cy cz ux uy uz)
        
        Bool_t  statusRh;
        Float_t stateRh[6]; // (cx cy cz ux uy uz)
        
        //Bool_t  statusLJ[9];
        //Float_t stateLJ[9][6]; // track state at layerJ (1 2 3 4 5 6 7 8 9)

        Float_t cpuTime; // [ms]
	
    ClassDef(CKTrackInfo, 2)
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
            
            std::fill_n(rig, 5, 0);
            std::fill_n(bta, 5, 0);

            statusTop = false;
            std::fill_n(stateTop, 8, 0);
            
            statusCen = false;
            std::fill_n(stateCen, 8, 0);
            
            statusBtm = false;
            std::fill_n(stateBtm, 8, 0);
            
            statusTd = false;
            std::fill_n(stateTd, 8, 0);
            
            statusRh = false;
            std::fill_n(stateRh, 8, 0);

            //std::fill_n(statusLJ, 9, false);
            //std::fill_n(stateLJ[0], 9*8, 0);

            cpuTime = 0;
        }

    public :
        Bool_t  status;
        Short_t ndof[2];
		Float_t nchi[2];
       
        Float_t rig[5]; // z = 195, 0, -136, 115, -75
        Float_t bta[5]; // z = 195, 0, -136, 115, -75

        Bool_t  statusTop;
        Float_t stateTop[8]; // (cx cy cz ux uy uz rig bta)
        
        Bool_t  statusCen;
        Float_t stateCen[8]; // (cx cy cz ux uy uz rig bta)
        
        Bool_t  statusBtm;
        Float_t stateBtm[8]; // (cx cy cz ux uy uz rig bta)
        
        Bool_t  statusTd;
        Float_t stateTd[8]; // (cx cy cz ux uy uz rig bta)
        
        Bool_t  statusRh;
        Float_t stateRh[8]; // (cx cy cz ux uy uz rig bta)

        //Bool_t  statusLJ[9];
        //Float_t stateLJ[9][8]; // track state at layerJ (1 2 3 4 5 6 7 8 9)
	
        Float_t cpuTime; // [ms]

    ClassDef(KFTrackInfo, 2)
};


// ShowerInfo
class ShowerInfo : public TObject {
	public :
		ShowerInfo() { init(); }
		~ShowerInfo() {}
		
		void init() {
            status = false;

			energyD = -1;
			energyE = -1;
			PisaBDT = -2;
			Q = -1;

			hadronApex = -1;
			hadronEnergy = -1;
		}

	public :
        Bool_t status;

		Float_t energyD;  // energy deposit [GeV]
		Float_t energyE;  // Pisa
		Float_t PisaBDT;  // Pisa
		Float_t Q;        // ecal charge

		Short_t hadronApex; // reject high-apex events (say, Apex > 12)
		Float_t hadronEnergy;

		ClassDef(ShowerInfo, 5)
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
            maxCfStormer = -1;
            maxCfIGRF = -1;
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
		Float_t cfStormer[4];        // Stormer cutoff for 25, 30, 35, 40 degrees [GV]
		Float_t cfIGRF[4];           // IGRF cutoff for 25, 30, 35, 40 degrees [GV]
        Float_t maxCfStormer;        // Stormer max cutoff [GV]
        Float_t maxCfIGRF;           // IGRF max cutoff [GV]
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
            noiseExtCls = 0;
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
        Short_t noiseExtCls; // noiseALL*1 + noiseUTOF*2 + noiseLTOF*4

	ClassDef(TOF, 11)
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
            QInMin = -1;

            hits.clear();

            cktr = std::vector<CKTrackInfo>(4);
            //kftr = std::vector<KFTrackInfo>(4);
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
		Float_t QInMin;
		
        // Track Hits
        std::vector<HitTRKInfo> hits;

        // Choutko [Inn InnL1 InnL9 FS]
        std::vector<CKTrackInfo> cktr;

        // Kalman Fitter [Inn InnL1 InnL9 FS]
        //std::vector<KFTrackInfo> kftr;

	ClassDef(TRK, 11)
};


// TRD
class TRD : public TObject {
	public :
		TRD() { init(); }
		~TRD() {}

		void init() {
            numOfCls = 0;
            numOfSegment = 0;
            numOfHSegment = 0;
			numOfTrack = 0;
			numOfHTrack = 0;
			
            trackStatus = false;
			std::fill_n(trackState, 6, 0);

			std::fill_n(LLRstatus, 2, false);
			std::fill_n(LLRnh, 2, -1);
			std::fill_n(LLRep, 2, -1);
			std::fill_n(LLReh, 2, -1);
			std::fill_n(LLRph, 2, -1);
			std::fill_n(Q, 2, -1);
            
			std::fill_n(ITstatus, 2, false);
			std::fill_n(ITnh, 2, 0);
			std::fill_n(ITcz, 2, 0);
			std::fill_n(ITamp, 2, 0);
			std::fill_n(ITlen, 2, 0);
			std::fill_n(ITdEdX, 2, 0);
			std::fill_n(ITMcMom, 2, 0);

            IThits[0].clear();
            IThits[1].clear();

            VTXstatus = false;
            VTXncls = 0;
            VTXnseg = 0;
            VTXcx = 0;
            VTXcy = 0;
            VTXcz = 0;
            
            OTHERncls = 0;
            OTHERnseg = 0;
		}

	public :
        Short_t numOfCls;
        Short_t numOfSegment;
		Short_t numOfTrack;
        Short_t numOfHSegment;
		Short_t numOfHTrack;
		
        // TrdTrack or TrdHTrack (first TrdH, second Trd)
		Bool_t  trackStatus;
		Float_t trackState[6]; // coo, dir

		// (TrdHTrack or TrdTrack) and TrTrack
		Bool_t  LLRstatus[2]; // true, rebuild success (Trd, Trk)
		Short_t LLRnh[2];
		Float_t LLRep[2];
		Float_t LLReh[2];
		Float_t LLRph[2];
		Float_t Q[2];

        // TRDRec
        Bool_t  ITstatus[2];
        Short_t ITnh[2];
        Float_t ITcz[2];
        Float_t ITamp[2];
        Float_t ITlen[2];
        Float_t ITdEdX[2];
        Float_t ITMcMom[2];

        std::vector<HitTRDSimpleInfo> IThits[2];

        // vertex seg
        Bool_t  VTXstatus;
        Short_t VTXncls;
        Short_t VTXnseg;
        Float_t VTXcx;
        Float_t VTXcy;
        Float_t VTXcz;
        
        // other seg
        Short_t OTHERncls;
        Short_t OTHERnseg;

	ClassDef(TRD, 13)
};


// RICH
class RICH : public TObject {
	public :
		RICH() { init(); }
		~RICH() {}

		void init() {
            numOfHit = 0;
            numOfCls = 0;
            numOfRing = 0;
            
            status = false;
            kind =  0;
            tile = -1;
            refz =  0;
            pmtz =  0;
			beta = -1;
			Q = -1;

			isGood = false;
            dist = -1;
            nhit = -1;
            npmt = -1;
            prob = -1.0;
            cstcb = -1.0;
            cstcq = -1.0;
            numOfClsPE = -1;
            numOfExpPE = -1;
            numOfColPE = -1;
            eftOfColPE = -1;
            
            betaCrr = 1;

            clsbta = -1;
            //uhits.clear();
            //ohits.clear();

            chfit.init();

            //vetoKind = -1;
			//vetoTile = -1;
			//vetoRfrIndex = -1;
			//vetoDistToBorder = -1;
			//vetoIsInFiducialVolume = false;

            //std::fill_n(vetoNumOfExpPE, 3, 0);
            //
            //vetoNumOfCrsHit = 0;
            //std::fill_n(vetoNumOfCkvHit[0], 3*3, 0);
            //
            //hits.clear();
		}

	public :
        Short_t numOfHit;
        Short_t numOfCls;
        Short_t numOfRing;

        // Official RichRingR
		Bool_t  status;
        Short_t kind;
        Short_t tile;
        Float_t refz;
        Float_t pmtz;
        Float_t dist;
		Float_t beta;
		Float_t Q;
        
		Bool_t  isGood;
        Short_t nhit;
        Short_t npmt;
        Float_t prob;
        Float_t cstcb; // consistency beta
        Float_t cstcq; // consistency charge
        Float_t numOfClsPE;
        Float_t numOfExpPE;
        Float_t numOfColPE;
        Float_t eftOfColPE;
    
        Float_t betaCrr;

        // Hits
        Float_t clsbta;
        //std::vector<HitRICHInfo> uhits;
        //std::vector<HitRICHInfo> ohits;

        ChFitInfo chfit;

		// Rich Veto
		//Short_t vetoKind;          // -1, None, 0, Aerogel 1, NaF
		//Short_t vetoTile;          // tile id
		//Float_t vetoRfrIndex;      // refractive index
		//Float_t vetoDistToBorder;  // dist To Border Edge
		//Bool_t  vetoIsInFiducialVolume;
       
        // [0] electron [1] proton [2] deuterium
		//Float_t vetoNumOfExpPE[3];     // number of photoelectrons expected for a given track, beta and charge.
        //Short_t vetoNumOfCrsHit;       // Particle Cross PMT
        //Short_t vetoNumOfCkvHit[3][3]; // Photo
                                       // [0] In  Ring
                                       // [1] Out Ring
                                       // [2] Non Physics

        // Rich Hits
        //std::vector<HitRICHInfo> hits;
	
    ClassDef(RICH, 11)
};


// ECAL
class ECAL : public TObject {
	public :
		ECAL() { init(); }
		~ECAL() {}

		void init() {
			numOfShower = 0;
            shower.init();

            numOfOther = 0;
            other.init();
		}

	public :
		Short_t    numOfShower;
		ShowerInfo shower;

        Short_t    numOfOther;
        ShowerInfo other;

	ClassDef(ECAL, 6)
};


// HCTrInfo
class HCTrInfo : public TObject {
	public :
		HCTrInfo() { init(); }
		~HCTrInfo() {}

		void init() {
            status = false;
            going = 0;
            chrg = 0;
            mass = 0;
            
            std::fill_n(ndof, 3, 0);
            std::fill_n(nchi, 3, 0);
            std::fill_n(qlt, 3, 0);
            
            std::fill_n(state, 8, 0);
            
            std::fill_n(rig, 5, 0);
            std::fill_n(bta, 5, 0);

            statusTop = false;
            std::fill_n(stateTop, 8, 0);
            
            statusCen = false;
            std::fill_n(stateCen, 8, 0);
            
            statusBtm = false;
            std::fill_n(stateBtm, 8, 0);
            
            statusTd = false;
            std::fill_n(stateTd, 8, 0);
            
            statusRh = false;
            std::fill_n(stateRh, 8, 0);
            
            //std::fill_n(statusLJ, 9, false);
            //std::fill_n(stateLJ[0], 9*8, 0);
           
            cpuTime = 0;
        }
	
    public :
        Bool_t  status;
        Short_t going; // 0, no  1, up-going  -1, down-going
        Short_t chrg;
        Float_t mass;
        
        Short_t ndof[3]; // x, y-beta, all
        Float_t nchi[3];
        Float_t qlt[3];
        
        Float_t state[8]; // (cx cy cz ux uy uz rig bta)
        
        Float_t rig[5]; // z = 195, 0, -136, 115, -75
        Float_t bta[5]; // z = 195, 0, -136, 115, -75
        
        Bool_t  statusTop; // track at (z = 195.)
        Float_t stateTop[8];
        
        Bool_t  statusCen; // track at (z = 0.)
        Float_t stateCen[8];
        
        Bool_t  statusBtm; // track at (z = -136.)
        Float_t stateBtm[8];
        
        Bool_t  statusTd; // track at (z = -75.)
        Float_t stateTd[8];
        
        Bool_t  statusRh; // track at (z = 115.)
        Float_t stateRh[8];

        //Bool_t  statusLJ[9];
        //Float_t stateLJ[9][8]; // track state at layerJ (1 2 3 4 5 6 7 8 9)

        Float_t cpuTime; // [ms]

        ClassDef(HCTrInfo, 6)
};


// HCBtaInfo
class HCBtaInfo : public TObject {
	public :
		HCBtaInfo() { init(); }
		~HCBtaInfo() {}

		void init() {
            status = false;
            going = 0;
            chrg = 0;
            mass = 0;

            ndof = 0;
            nchi = 0;
            qlt = 0;
            
            rerr = 0;
            
            std::fill_n(state, 8, 0);
            
            std::fill_n(rig, 5, 0);
            std::fill_n(bta, 5, 0);
            
            cpuTime = 0;
        }
    
    public :
        Bool_t  status;
        Short_t going; // 0, no  1, up-going  -1, down-going
        Short_t chrg;
        Float_t mass;
        
        Short_t ndof;
        Float_t nchi;
        Float_t qlt;
        
        Float_t rerr;

        Float_t state[8]; // (cx cy cz ux uy uz rig bta)
        
        Float_t rig[5]; // z = 195, 0, -136, 115, -75
        Float_t bta[5]; // z = 195, 0, -136, 115, -75
        
        Float_t cpuTime; // [ms]

        ClassDef(HCBtaInfo, 2)
};


// HCMuInfo
class HCMuInfo : public TObject {
	public :
		HCMuInfo() { init(); }
		~HCMuInfo() {}

		void init() {
            status = false;
            going = 0;
            chrg = 0;
            mass = 0;
            
            std::fill_n(muNdof, 3, 0);
            std::fill_n(muNchi, 3, 0);
            std::fill_n(muQlt, 3, 0);
            
            std::fill_n(ndof, 3, 0);
            std::fill_n(nchi, 3, 0);
            std::fill_n(qlt,  3, 0);
            
            std::fill_n(state, 8, 0);
            
            std::fill_n(rig, 5, 0);
            std::fill_n(bta, 5, 0);

            statusTop = false;
            std::fill_n(stateTop, 8, 0);
            
            statusCen = false;
            std::fill_n(stateCen, 8, 0);
            
            statusBtm = false;
            std::fill_n(stateBtm, 8, 0);
            
            statusTd = false;
            std::fill_n(stateTd, 8, 0);
            
            statusRh = false;
            std::fill_n(stateRh, 8, 0);
           
            cpuTime = 0;
        }
	
    public :
        Bool_t  status;
        Short_t going; // 0, no  1, up-going  -1, down-going
        Short_t chrg;
        Float_t mass;
        
        Short_t muNdof[3]; // x, y, beta
        Float_t muNchi[3];
        Float_t muQlt[3];
        
        Short_t ndof[3]; // x, y-beta, all
        Float_t nchi[3];
        Float_t qlt[3];
        
        Float_t state[8]; // (cx cy cz ux uy uz rig bta)
        
        Float_t rig[5]; // z = 195, 0, -136, 115, -75
        Float_t bta[5]; // z = 195, 0, -136, 115, -75
        
        Bool_t  statusTop; // track at (z = 195.)
        Float_t stateTop[8];
        
        Bool_t  statusCen; // track at (z = 0.)
        Float_t stateCen[8];
        
        Bool_t  statusBtm; // track at (z = -136.)
        Float_t stateBtm[8];
        
        Bool_t  statusTd; // track at (z = -75.)
        Float_t stateTd[8];
        
        Bool_t  statusRh; // track at (z = 115.)
        Float_t stateRh[8];

        Float_t cpuTime; // [ms]

        ClassDef(HCMuInfo, 2)
};


// HYC
class HYC : public TObject {
	public :
		HYC() { init(); }
		~HYC() {}

		void init() {
            trPrm    = std::vector<HCTrInfo>(4);
            trPrmAll = std::vector<HCTrInfo>(4);
            btaPrm.init();
            massPrm = 0.0;

            trSecIn.init();
            trSecAllIn.init();
            btaSec.init();
            massSec = 0.0;

            mutr.init();

            //trM1 = std::vector<HCTrInfo>(4);
            //trM2 = std::vector<HCTrInfo>(4);
            //
            //trM1All = std::vector<HCTrInfo>(4);
            //trM2All = std::vector<HCTrInfo>(4);
         
            //btaM1.init();
            //btaM2.init();

            //mutr.init();

            //massM1 = 0.0;
            //massM2 = 0.0;
		}
	
    public :
        // Two main particles (P/D or He4/He3)
        
        // first (P or He4)
        std::vector<HCTrInfo> trPrm;
        std::vector<HCTrInfo> trPrmAll;
        HCBtaInfo             btaPrm;
        Float_t               massPrm;

        // second (D or He3)
        HCTrInfo  trSecIn;
        HCTrInfo  trSecAllIn;
        HCBtaInfo btaSec;
        Float_t   massSec;
        
        // Mass Fit [Inn]
        HCMuInfo mutr;

        // Two main particles (P/D or He4/He3)
        // M1: P or He4
        // M2: D or He3

        // Track Fit [Inn InnL1 InnL9 FS]
        //std::vector<HCTrInfo> trM1;
        //std::vector<HCTrInfo> trM2;

        //std::vector<HCTrInfo> trM1All;
        //std::vector<HCTrInfo> trM2All;
        
        // Beta Fit [Inn]
        //HCBtaInfo btaM1;
        //HCBtaInfo btaM2;

        // Mass Fit [Inn]
        //HCMuInfo mutr;
        
        // Mass from Track&Beta Fit [Inn]
        //Float_t massM1;
        //Float_t massM2;

	ClassDef(HYC, 5)
};


#endif // __ClassDef_H__
