/**********************************************************
 * Author        : Hsin-Yi Chou
 * Email         : hchou@cern.ch
 * Last modified : 2015-07-22 16:29
 * Filename      : YiProdNtuple.h
 * Description   :
 * *******************************************************/
#ifndef __YiProdNtuple_H__
#define __YiProdNtuple_H__

// C++ library
#include <iostream>
#include <fstream>
#include <cstring>
#include <cctype>
#include <locale>
#include <vector>
#include <utility>
#include <cmath>
#include <functional>
#include <algorithm>
#include <cstdio>

// ROOT library
#include <TString.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TF1.h>
#include <TF2.h>

// AMS library
#include <root.h>
#include <amschain.h>
#include <root_setup.h>
#include <Tofrec02_ihep.h>
#include <TrdKCluster.h>
#include <TrExtAlignDB.h>
#include <TrTrack.h>
#include <tkdcards.h>
#include <root_RVSP.h>
#include <TofTrack.h>
#include <TrFit.h>
#include <bcorr.h>
#include <EcalChi2CY.h>
#include <TkDBc.h>
#include <TkSens.h>
#include <TrCharge.h>
#include <richradidOff.h>
#include <richtrrecOff.h>
#include <GeoMagField.h>
#include <GeoMagTrace.h>
#include <TrReconQ.h>
#include <MagField.h>
#include <TrMass.h>

// User defination library
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>

#define __HAS_AMS_OFFICE_LIBS__
#include <TRACKSys.h>

#include <EcalHadron.h>
#include <EcalHadron.C>
#include <TRDVertex.h>

#include "ClassDef.h"

// User defination macro
#define Debug true


//---- RecEvent ----//
class RecEvent {
	public :
		RecEvent() { init(); }
		~RecEvent() {}

		void init();
		bool rebuild(AMSEventR * event);
		
		inline float time() { return fStopwatch.time(); }

	public :
		int iBeta;
		int iBetaH;
		int iTrTrack;
		int iEcalShower;
		int iTrdTrack;
		int iTrdHTrack;
		int iRichRing;

	public :
		double trackerZJ[9];

    // Track Refit Option
	// Algorithm     (CHOUTKO, KALMAN, HCHOU)
	// Track Pattern (Inn, InnL1, InnL9, FS)
    public :
        inline void initTrRft() { std::fill_n(_trRft[0], 3*4, 22); }
        inline int trRft(int i, int j) { 
            if (_trRft[i][j]==22) { _trRft[i][j] = 21; return 22; }
            else return _trRft[i][j];
        }

    protected :
        int _trRft[3][4];

	protected :
		MGClock::HrsStopwatch fStopwatch;
};

static RecEvent recEv;


//---- EventBase ----//
class EventBase {
	public :
		enum MODE {
			ISS, BT, MC
		};
		static enum MODE eventMode;
		static void setEventMode(EventBase::MODE mode) {EventBase::eventMode = mode;}
		static bool checkEventMode(EventBase::MODE mode) {return (EventBase::eventMode == mode);}

		enum VERSION {
			B620, B950
		};
		static enum VERSION eventVersion;
		static void setEventVersion(EventBase::VERSION version) {EventBase::eventVersion = version;}
		static bool checkEventVersion(EventBase::VERSION version) {return (EventBase::eventVersion == version);}

	public :
		EventBase();
		~EventBase();

		virtual void initEvent() = 0;
		virtual void setEventTree(TTree * evTree = nullptr) = 0;
		virtual void setEnvironment() = 0;
		virtual bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr) = 0;
		virtual bool selectEvent(AMSEventR * event = nullptr) = 0;
		inline void fill();

		inline float time() { return fStopwatch.time(); }

	protected :
		inline void setTree(const std::string& name, const std::string& title);
		inline void deleteTree();

	protected :
		bool isTreeSelf;
		TTree * tree;

	protected :
		MGClock::HrsStopwatch fStopwatch;
};


//---- EventList ----//
class EventList : virtual public EventBase {
	public :
		EventList();
		~EventList();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
		LIST fList;
		G4MC fG4mc;

		static Float_t Weight;
};
Float_t EventList::Weight = 1;


//---- EventRti ----//
class EventRti : virtual public EventBase {
	public :
		EventRti();
		~EventRti();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
		RTI fRti;

	protected :
		static UInt_t CurrUTime; // current utime
};

UInt_t EventRti::CurrUTime = 0;


//---- EventTrg ----//
class EventTrg : virtual public EventBase {
	public :
		EventTrg();
		~EventTrg();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
		TRG fTrg;
};


//---- EventTof ----//
class EventTof : virtual public EventBase {
	public :
		EventTof();
		~EventTof();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
		TOF fTof;
};


//---- EventAcc ----//
class EventAcc : virtual public EventBase {
	public :
		EventAcc();
		~EventAcc();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
	ACC fAcc;
};


//---- EventTrk ----//
class EventTrk : virtual public EventBase {
	public :
		EventTrk();
		~EventTrk();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
		TRK fTrk;
};


//---- EventTrd ----//
class EventTrd : virtual public EventBase {
	public :
		EventTrd();
		~EventTrd();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
		TRD fTrd;
};


//---- EventRich ----//
class EventRich : virtual public EventBase {
	public :
		EventRich();
		~EventRich();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
		RICH fRich;
};

//---- EventEcal ----//
class EventEcal : virtual public EventBase {
	public :
		EventEcal();
		~EventEcal();

		void initEvent();
		void setEventTree(TTree * evTree = nullptr);
		void setEnvironment();
		bool processEvent(AMSEventR * event = nullptr, AMSChain * chain = nullptr);
		bool selectEvent(AMSEventR * event = nullptr);

	public :
		ECAL fEcal;
};


//---- DataSelection ----//
class DataSelection {
	public :
		enum SWITCH {
			ON, OFF
		};
		enum OPTION {
			LIST, RTI, TRG, TOF, ACC, TRK, TRD, RICH, ECAL, NUMBER
		};
		static enum SWITCH option[DataSelection::NUMBER];
		static void setOption(DataSelection::OPTION opt, DataSelection::SWITCH sw) {option[opt] = sw;}
		static bool checkOption(DataSelection::OPTION opt) {return (DataSelection::option[opt] == DataSelection::ON);}

	public :
		DataSelection();
		~DataSelection();

		void setMultiTree(bool isMTree = false);
		void setEventTree();
		void setEnvironment();
		void fill();

		int processEvent(AMSEventR * event, AMSChain * chain);
		int preselectEvent(AMSEventR * event, const std::string& officialDir = "");
		int selectEvent(AMSEventR * event);
		int analysisEvent(AMSEventR * event);

	protected :
		bool isMultiTree;
		TTree * evTree;
		EventList list;
		EventRti rti;
		EventTrg trg;
		EventTof tof;
		EventAcc acc;
		EventTrk trk;
		EventTrd trd;
		EventRich rich;
		EventEcal ecal;

	public :
		static Float_t gScaleFact;
		//static TF1     gScaleFunc1D; // (chrg)
		static TF1     gScaleFunc1D; // (rig)
		static TF2     gScaleFunc2D; // (rig, chrg)
};

Float_t DataSelection::gScaleFact = 0.02;
//TF1     DataSelection::gScaleFunc1D("gScaleFunc1D", "0.5*((1.0+[0])+(1.0-[0])*TMath::Erf(x*x-6.0))", 0, 10);
TF1     DataSelection::gScaleFunc1D("gScaleFunc1D", "0.5*((1.0+[0])+(1.0-[0])*TMath::Erf(0.75*(TMath::Log(TMath::Abs(x))-4.5)))*(x>0)+(x<=0)", -2000, 2000);
TF2     DataSelection::gScaleFunc2D("gScaleFunc2D", "0.5*((1.0+0.5*((1.0+[0])+(1.0-[0])*TMath::Erf(y*y-6.0)))+(1.0-0.5*((1.0+[0])+(1.0-[0])*TMath::Erf(y*y-6.0)))*TMath::Erf(0.75*(TMath::Log(TMath::Abs(x))-4.5)))*(x>0)+(x<=0)", -2000, 2000, 0, 10);


//---- RunTagOperator ----//
class RunTagOperator {
	public :
		RunTagOperator();
		~RunTagOperator();
		
		void init();
		bool processEvent(AMSEventR * event, AMSChain * chain);
		void save(TFile * file);

	protected :
		std::map<UInt_t, RunTagInfo> fRunTag;
};


//---- YiNtuple ----//
// We have some Bug in AMSChain->AddFromFile().
// Limit of Files -> Maybe we can set 30 files.
class YiNtuple {
	public :
		enum MODE {
			NORM, COPY
		};
		static enum MODE selectionMode;
		static void setSelectionMode(YiNtuple::MODE mode) {YiNtuple::selectionMode = mode;}
		static bool checkSelectionMode(YiNtuple::MODE mode) {return (YiNtuple::selectionMode == mode);}

	public :
		YiNtuple();
		~YiNtuple();

		inline void init();
		inline void setOutputFile(const std::string& file_name = "YiNtuple.root", const std::string& path = ".", bool isMultiTree = false);
		void readDataFrom(const std::string& file_list = "fileList.txt", Long64_t group_id = 0, Long64_t group_size = -1);
        void saveInputFileList(TFile * file);
		void loopEventChain();

	protected :
		std::pair<Long64_t, Long64_t> fGroup;
		std::vector<std::string>      fFileList;
		std::string                   fFileDir;

		std::string      fFileName;
		AMSChain       * fChain;
		DataSelection  * fData;
		RunTagOperator * fRunTagOp; 

	protected :
		MGClock::HrsStopwatch fStopwatch;
};


#endif // __YiProdNtuple_H__
