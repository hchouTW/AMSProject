#ifndef __Analyzer_H__
#define __Analyzer_H__

// STL c++
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdarg>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>

// ROOT
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

// User defination library
#include "ClassDef.h"

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <CPPLibs.h>
#include <ROOTLibs.h>
#include <TrSys.h>

using namespace MGROOT;

template<typename... Args>
inline std::string Format(const std::string& fmt, Args... args) {
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
            outstr += Format("========================  Stopwatch  ==========================\n");
            outstr += Format("==  START TIME : Unix( %ld )    %s", utime1, std::asctime(std::gmtime(&utime1)));
            outstr += Format("==  STOP  TIME : Unix( %ld )    %s", utime2, std::asctime(std::gmtime(&utime2)));
            outstr += Format("==  Duration   : %-3d HR %-2d MIN %6.3f   SEC (%17.6f)\n", hours.count(), minutes.count(), seconds.count(), time);
            outstr += Format("===============================================================\n");
            return outstr;
        }

    private :
        std::pair<Time, Time> times_;
};


inline std::vector<std::string> ReadListFile(const std::string& path) {
    if (std::system((Format("test -f \"%s\"", path.c_str())).c_str()) != 0) return std::vector<std::string>();
	
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


class Analyzer {
    public :
		enum class Type {
			ISS, BT, MC
		};
		static Type kType;
		static void SetType(Type type) { kType = type; }
		static bool CheckType(Type type) {return (kType == type); }
		
    public :
        static TRandom3 RndmGenerator;

    private :
        std::unordered_map<UInt_t, UInt_t> runev;
        UInt_t map_run;
        UInt_t map_ev;
        
        TChain* mdst;
        LIST*   list;
        G4MC*   g4mc;
        RTI*    rti;
        TRG*    trg;
        ACC*    acc;
        TOF*    tof;
        TRK*    trk;
        TRD*    trd;
        ECAL*   ecal;
        RICH*   rich;
        HYC*    hyc;

        TFile* file;
        TTree* tree;
        TTree* clone_tree;

    public :
        Analyzer(TChain* chain) { init(); set_branch(chain); }
        ~Analyzer() { init(); }

        inline void set_output(const std::string& outpath);
        inline void write();
        inline void close();

        void set_environment();
        void process_events();

    protected :
        inline void init();
        inline void set_branch(TChain* chain);

        bool build_tree();
        bool build_hist();

        bool process_presel();
        bool process_data();
        
        bool process_data_prof();
        bool process_data_l();
        bool process_data_m();
        bool process_data_i();
        bool process_data_h();
};

Analyzer::Type Analyzer::kType = Analyzer::Type::ISS;
TRandom3 Analyzer::RndmGenerator = TRandom3(0);

void Analyzer::init() {
    runev.clear();
    map_run = 0;
    map_ev = 0;

    mdst = nullptr;
    list = nullptr;
    g4mc = nullptr;
    rti  = nullptr;
    trg  = nullptr;
    acc  = nullptr;
    tof  = nullptr;
    trk  = nullptr;
    trd  = nullptr;
    ecal = nullptr;
    rich = nullptr;
    hyc  = nullptr;

    file = nullptr;
    tree = nullptr;
    clone_tree = nullptr;
}

inline void Analyzer::set_branch(TChain* chain) {
    if (chain == nullptr) return;
    mdst = chain;

    list = new LIST;
    if (CheckType(Type::MC))  g4mc = new G4MC;
    if (CheckType(Type::ISS)) rti  = new RTI;
    trg  = new TRG;
    acc  = new ACC;
    tof  = new TOF;
    trk  = new TRK;
    trd  = new TRD;
    ecal = new ECAL;
    rich = new RICH;
    hyc  = new HYC;

    mdst->SetBranchAddress("list", &list);
    if (CheckType(Type::MC))  mdst->SetBranchAddress("g4mc", &g4mc);
    if (CheckType(Type::ISS)) mdst->SetBranchAddress("rti",  &rti);
    mdst->SetBranchAddress("trg",  &trg);
    mdst->SetBranchAddress("acc",  &acc);
    mdst->SetBranchAddress("tof",  &tof);
    mdst->SetBranchAddress("trk",  &trk);
    mdst->SetBranchAddress("trd",  &trd);
    mdst->SetBranchAddress("ecal", &ecal);
    mdst->SetBranchAddress("rich", &rich);
    mdst->SetBranchAddress("hyc",  &hyc);
}


void Analyzer::set_output(const std::string& outpath) {
    file = new TFile(outpath.c_str(), "RECREATE");

    if (file->IsZombie()) {
        std::cerr << Format("Error opening file: %s\n", outpath.c_str());
        LOG(FATAL) << Format("Error opening file: %s", outpath.c_str());
        file = nullptr;
        return;
    }

    tree = new TTree("ana", "data");
    tree->Branch("run", &map_run);
    tree->Branch("event", &map_ev);

    clone_tree = mdst->CloneTree(0);

    std::string statement;
    statement += Format("\n****  Output Info ****\n");
    statement += Format("File : %s\n", outpath.c_str());
    statement += Format("Tree : %s\n\n", tree->GetName());
    statement += Format("Clone: %s\n\n", clone_tree->GetName());
    std::cout << statement;
    LOG(INFO) << statement;
}


void Analyzer::write() {
    if (tree != nullptr && runev.size() != 0) {
        for (auto&& elem : runev) {
            map_run = elem.first;
            map_ev  = elem.second;
            tree->Fill();
        }
    }

    if (tree == nullptr || file == nullptr) return;
    file->cd();
    file->Write();
}

void Analyzer::close() {
    if (tree == nullptr || file == nullptr) return;
    file->Close();
    file = nullptr;
    tree = nullptr;
}

#endif // __Analyzer_H__
