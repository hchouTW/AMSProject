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

// Root library
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TSystem.h"

// AMS library
#include "root.h"
#include "amschain.h"
#include "root_setup.h"
#include "Tofrec02_ihep.h"
#include "TrdKCluster.h"
#include "TrExtAlignDB.h"
#include "TrTrack.h"
#include "tkdcards.h"
#include "root_RVSP.h"
#include "TofTrack.h"
#include "TrFit.h"
#include "bcorr.h"
#include "EcalChi2CY.h"
#include "TkDBc.h"
#include "TkSens.h"
#include "TrCharge.h"
#include "richradidOff.h"
#include "richtrrecOff.h"
#include "GeoMagField.h"
#include "GeoMagTrace.h"
#include "TrReconQ.h"

// User defination library
#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/ROOT/selflib.h"

#include "/afs/cern.ch/user/h/hchou/private/YiService/src/EcalHadron/EcalHadron.h"
#include "/afs/cern.ch/user/h/hchou/private/YiService/src/EcalHadron/EcalHadron.cpp"

#include "/afs/cern.ch/user/h/hchou/private/YiService/src/MgntFit/Include.h"

#include "ClassDef.h"
#include "ClassDef.C"

// User defination macro
#define Debug true


//---- RecEvent ----//
class RecEvent {
	public :
		RecEvent() { init(); }
		~RecEvent() {}

		void init();
		bool rebuild(AMSEventR * event);
		
		inline float time() { return timer.time(); }

	public :
		int iBeta;
		int iBetaH;
		int iTrTrack;
		int iEcalShower;
		int iTrdTrack;
		int iTrdHTrack;
		int iRichRing;
	
	protected :
		MgntClock::HrsTimer timer;
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
		virtual void setEventTree(TTree * evTree = 0) = 0;
		virtual void setEnvironment() = 0;
		virtual bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0) = 0;
		virtual bool selectEvent(AMSEventR * event = 0) = 0;
		inline void fill();

		inline float time() { return timer.time(); }

	protected :
		inline void setTree(const char * name, const char * title);
		inline void deleteTree();

	protected :
		bool isTreeSelf;
		TTree * tree;

	protected :
		MgntClock::HrsTimer timer;
};


//---- EventList ----//
class EventList : virtual public EventBase {
	public :
		EventList();
		~EventList();

		void initEvent();
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

	public :
		LIST fList;
		G4MC fG4mc;

		static Float_t gWeight;
};
Float_t EventList::gWeight = 1;


//---- EventRti ----//
class EventRti : virtual public EventBase {
	public :
		EventRti();
		~EventRti();

		void initEvent();
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

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
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

	public :
		TRG fTrg;
};


//---- EventTof ----//
class EventTof : virtual public EventBase {
	public :
		EventTof();
		~EventTof();

		void initEvent();
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

	public :
		TOF fTof;
};


//---- EventAcc ----//
class EventAcc : virtual public EventBase {
	public :
		EventAcc();
		~EventAcc();

		void initEvent();
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

	public :
	ACC fAcc;
};


//---- EventTrk ----//
class EventTrk : virtual public EventBase {
	public :
		EventTrk();
		~EventTrk();

		void initEvent();
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

	public :
		TRK fTrk;
};


//---- EventTrd ----//
class EventTrd : virtual public EventBase {
	public :
		EventTrd();
		~EventTrd();

		void initEvent();
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

	public :
		TRD fTrd;
};


//---- EventRich ----//
class EventRich : virtual public EventBase {
	public :
		EventRich();
		~EventRich();

		void initEvent();
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

	public :
		RICH fRich;
};

//---- EventEcal ----//
class EventEcal : virtual public EventBase {
	public :
		EventEcal();
		~EventEcal();

		void initEvent();
		void setEventTree(TTree * evTree = 0);
		void setEnvironment();
		bool processEvent(AMSEventR * event = 0, AMSChain * chain = 0);
		bool selectEvent(AMSEventR * event = 0);

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
		int preselectEvent(AMSEventR * event);
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
		static TRandom3 gRandom;
		static Float_t  gScaleFact;
		static TF1      gScaleFunc;
};

TRandom3 DataSelection::gRandom(0);
Float_t  DataSelection::gScaleFact = 0.01;
TF1      DataSelection::gScaleFunc("gScaleFunc", "0.5*((1.0+[0])+(1.0-[0])*TMath::Erf(0.75*(TMath::Log(TMath::Abs(x))-4.5)))*(x>0)+(x<=0)", -2000, 2000);


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
		void readDataFrom(const std::string& file_list = "fileList.txt", Long64_t group_th = 0, Long64_t group_size = -1);
    void saveInputFileList(TFile * file);
		void loopEventChain();

	protected :
		std::pair<Long64_t, Long64_t> fGroup;
		std::vector<std::string>      fFileList;

		std::string      fFileName;
		AMSChain       * fChain;
		DataSelection  * fData;
		RunTagOperator * fRunTagOp; 

	protected :
		MgntClock::HrsTimer fTimer;
		
	public :
		static std::string FileDir;
};

std::string YiNtuple::FileDir = "";


#endif // __YiProdNtuple_H__
