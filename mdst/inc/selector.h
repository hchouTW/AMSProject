#ifndef __Selector_H__
#define __Selector_H__

// STL c++
#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <cstdarg>
#include <vector>
#include <fstream>
#include <chrono>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

// AMS
#include <root.h>
#include <amschain.h>
#include <root_setup.h>
#include <Tofrec02_ihep.h>
#include <TrdKHit.h>
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

#include "EcalHadron/EcalHadron.h"
#include "EcalHadron/EcalHadron.C"

#include "AmsRich/AmsRich.h"
#include "AmsRich/AmsRich.C"

#include "AmsRich/CherenkovMeas.h"
#include "AmsRich/CherenkovMeas.C"

#include "AmsRich/ChSearch.h"
#include "AmsRich/ChSearch.C"

#include "TrSys.h"

// User defination library
#include "ClassDef.h"

#include <gflags/gflags.h>
#include <glog/logging.h>


template<typename... Args>
inline std::string Fmt(const std::string& fmt, Args... args) {
	std::vector<char> buf(1+std::snprintf(nullptr, 0, fmt.c_str(), args...));
	std::snprintf(buf.data(), buf.size(), fmt.c_str(), args...);
	std::string&& str = std::string(buf.begin(), buf.end());
	str.erase(std::remove_if(str.begin(), str.end(), ([](const char& ch)->bool{return (ch==char('\0'));}) ), str.end());
	return str;
}



class Stopwatch {
    public :
        using Clock    = std::chrono::high_resolution_clock;
        using Time     = std::chrono::high_resolution_clock::time_point;
        using Duration = std::chrono::high_resolution_clock::duration; 
        
        using FloatSeconds  = std::chrono::duration<double>;
        using Nanoseconds   = std::chrono::nanoseconds;
        using Microseconds  = std::chrono::microseconds;
        using Milliseconds  = std::chrono::milliseconds;
        using Seconds       = std::chrono::seconds;
        using Minutes       = std::chrono::minutes;
        using Hours         = std::chrono::hours; 

    public :
        Stopwatch() { init(); }
        ~Stopwatch() {}
        
        inline void init() { times_.first = Clock::now(); times_.second = Clock::now(); }

        inline Time& start() { times_.first  = Clock::now(); return times_.first;  } 
        inline Time& stop()  { times_.second = Clock::now(); return times_.second; }

        inline Duration duration() const { return (times_.second - times_.first); }
        inline double   time()     const { return std::chrono::duration<double>( (times_.second - times_.first) ).count(); }

        inline std::string str() const {
            Duration&& durt = duration();
            double time = std::chrono::duration<double>(durt).count();
            Hours        hours   = std::chrono::duration_cast<Hours>(durt);   durt -= hours;
            Minutes      minutes = std::chrono::duration_cast<Minutes>(durt); durt -= minutes;
            FloatSeconds seconds = std::chrono::duration_cast<FloatSeconds>(durt);

            std::string outstr;
            std::time_t utime1 = Clock::to_time_t(times_.first);
            std::time_t utime2 = Clock::to_time_t(times_.second);
            outstr += Fmt("========================  Stopwatch  ==========================\n");
            outstr += Fmt("==  START TIME : Unix( %ld )    %s", utime1, std::asctime(std::gmtime(&utime1)));
            outstr += Fmt("==  STOP  TIME : Unix( %ld )    %s", utime2, std::asctime(std::gmtime(&utime2)));
            outstr += Fmt("==  Duration   : %-3d HR %-2d MIN %6.3f   SEC (%17.6f)\n", hours.count(), minutes.count(), seconds.count(), time);
            outstr += Fmt("===============================================================\n");
            return outstr;
        }

    private :
        std::pair<Time, Time> times_;
};


inline std::vector<std::string> ReadListFile(const std::string& path) {
    if (std::system((Fmt("test -f \"%s\"", path.c_str())).c_str()) != 0) return std::vector<std::string>();
	
    std::fstream fstr;
	fstr.open(path, std::ios::in);
    if (!fstr.is_open()) return std::vector<std::string>();
	fstr.seekg(0);
    
    std::vector<std::string> list;
	for (std::string line; std::getline(fstr, line); ) {
        list.push_back(line);
    }
    fstr.close();

	return list;
}


class Selector {
    public :
		enum class Type {
			ISS, BT, MC
		};
		static Type kType;
		static void SetType(Type type) { kType = type; }
		static bool CheckType(Type type) {return (kType == type); }
		
    public :
        static double GetTkZJ(int layJ) { return (CheckType(Type::ISS) ? TkDBc::Head->GetZlayerAJ(layJ) :TkDBc::Head->GetZlayerJ(layJ)); }

    public :
        static TRandom3 RndmGenerator;

    private :
        AMSChain*    amsch;
        AMSEventR*   event;
        ParticleR*   pParticle;
        BetaHR*      pBetaH;
        TrTrackR*    pTrTrack;
        TrdTrackR*   pTrdTrack;
        TrdHTrackR*  pTrdHTrack;
        EcalShowerR* pEcalShower;
        RichRingR*   pRichRing;

        double dist_tk_td;
        double dist_tk_tdh;
        double dist_tk_ecal;
        
        TFile* file;
        TTree* tree;
        
        TTree* tree_runenv;
        RUNENV data_runenv;
        
        TTree* tree_decenv;
        DECENV data_decenv;
        
        int    data_zin;
        double data_rin;
        double data_mass;
        double data_beta;
        double data_msqr;

        LIST data_list;
        G4MC data_g4mc;
        RTI  data_rti;
        TRG  data_trg;
        ACC  data_acc;
        TOF  data_tof;
        TRK  data_trk;
        TRD  data_trd;
        ECAL data_ecal;
        RICH data_rich;
        HYC  data_hyc;

    public :
        Selector(AMSChain* ams) { init(); amsch = ams; }
        ~Selector() { init(); }

        inline void set_output(const std::string& outpath);
        inline void write();
        inline void close();

        void set_environment();
        void process_events();

    protected :
        inline void init();

        inline void process_init();
       
        bool process_runenv();
        bool process_decenv(UInt_t beg_event = 0, UInt_t end_event = 0);
        bool process_prefix();
        bool process_data();
        bool process_prd();
        bool process_sel();

        bool process_list();
        bool process_g4mc();
        bool process_rti();
        bool process_trg();
        bool process_acc();
        bool process_tof();
        bool process_trk();
        bool process_trd();
        bool process_ecal();
        bool process_rich();
        bool process_hyc();
};

Selector::Type Selector::kType = Selector::Type::ISS;
TRandom3 Selector::RndmGenerator = TRandom3(0);


void Selector::init() {
    amsch = nullptr;

    file = nullptr;
    tree = nullptr;
    
    tree_runenv = nullptr;
    tree_decenv = nullptr;

    process_init();
}

void Selector::process_init() {
    event       = nullptr;
    pParticle   = nullptr;
    pBetaH      = nullptr;
    pTrTrack    = nullptr;
    pTrdTrack   = nullptr;
    pTrdHTrack  = nullptr;
    pEcalShower = nullptr;
    pRichRing   = nullptr;

    dist_tk_td   = 0.0;
    dist_tk_tdh  = 0.0;
    dist_tk_ecal = 0.0;

    data_zin  = 1.0;
    data_rin  = 0.0;
    data_mass = TrFit::Mproton;
    data_beta = 1.0;
    data_msqr = 0.0;

    data_runenv.init();
    data_decenv.init();

    data_list.init();
    data_g4mc.init();
    data_rti .init();
    data_trg .init();
    data_acc .init();
    data_tof .init();
    data_trk .init();
    data_trd .init();
    data_ecal.init();
    data_rich.init();
    data_hyc .init();
}
        

void Selector::set_output(const std::string& outpath) {
    file = new TFile(outpath.c_str(), "RECREATE");

    if (file->IsZombie()) {
        std::cerr << Fmt("Error opening file: %s\n", outpath.c_str());
        LOG(FATAL) << Fmt("Error opening file: %s", outpath.c_str());
        file = nullptr;
        return;
    }
    
    tree_runenv = new TTree("runenv", "runenv");
    tree_runenv->Branch("runenv", &data_runenv);
    
    tree_decenv = new TTree("decenv", "decenv");
    tree_decenv->Branch("decenv", &data_decenv);
    
    tree = new TTree("mdst", "mdst");
    tree->Branch("list", &data_list);
    //if (CheckType(Type::MC )) tree->Branch("g4mc", &data_g4mc);
    //if (CheckType(Type::ISS)) tree->Branch("rti" , &data_rti );
    //tree->Branch("trg" , &data_trg );
    //tree->Branch("acc" , &data_acc );
    //tree->Branch("tof" , &data_tof );
    tree->Branch("trk" , &data_trk );
    //tree->Branch("trd" , &data_trd );
    //tree->Branch("ecal", &data_ecal);
    //tree->Branch("rich", &data_rich);
    //tree->Branch("hyc" , &data_hyc );
    
    std::string statement;
    statement += Fmt("\n****  Output Info ****\n");
    statement += Fmt("File: %s\n", outpath.c_str());
    statement += Fmt("Tree: %s\n\n", tree->GetName());
    std::cout << statement;
    LOG(INFO) << statement;
}


void Selector::write() {
    if (tree == nullptr || file == nullptr) return;
    file->cd();
    file->Write();
}

void Selector::close() {
    if (tree == nullptr || file == nullptr) return;
    file->Close();
    file = nullptr;
    tree = nullptr;
    tree_runenv = nullptr;
    tree_decenv = nullptr;
}

#endif // __Selector_H__
