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

#include <TMVA/Tools.h>
#include <TMVA/Reader.h>

// User defination library
#include <gflags/gflags.h>
#include <glog/logging.h>

#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "Variables.h"

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
        static UInt_t FlxUTime;

    private :
        TChain* runlist;
        VarsLTF varsLTF;
        VarsLRH varsLRH;
        VarsIIN varsIIN;
        VarsIEX varsIEX;
        VarsHEX varsHEX;
        VarsHFS varsHFS;

        TFile* file;
        TTree* tree;

        TTree*  aptree;
        Float_t apcc;
        
    public :
        Analyzer(TChain* list, TChain* chainLTF, TChain* chainLRH, TChain* chainIIN, TChain* chainIEX, TChain* chainHEX, TChain* chainHFS) { init(); set_branch(list, chainLTF, chainLRH, chainIIN, chainIEX, chainHEX, chainHFS); }
        ~Analyzer() { init(); }

        inline void set_output(const std::string& outpath);
        inline void write();
        inline void close();

        void set_environment();
        void process_events();

    protected :
        inline void init();
        inline void set_branch(TChain* list, TChain* chainLTF, TChain* chainLRH, TChain* chainIIN, TChain* chainIEX, TChain* chainHEX, TChain* chainHFS);

        bool build_hist();
        bool process_data(TChain* mdst, bool (Analyzer::*process)());

        bool process_data_ltf();
        bool process_data_lrh();
        bool process_data_iin();
        bool process_data_iex();
        bool process_data_hex();
        bool process_data_hfs();
};

Analyzer::Type Analyzer::kType = Analyzer::Type::ISS;
UInt_t Analyzer::FlxUTime = 0;

void Analyzer::init() {
    runlist = nullptr;

    varsLTF.init();
    varsLRH.init();
    varsIIN.init();
    varsIEX.init();
    varsHEX.init();
    varsHFS.init();

    file = nullptr;
    tree = nullptr;

    aptree = nullptr;
    apcc = 0.0;
}

inline void Analyzer::set_branch(TChain* list, TChain* chainLTF, TChain* chainLRH, TChain* chainIIN, TChain* chainIEX, TChain* chainHEX, TChain* chainHFS) {
    if (list    != nullptr) runlist = list;
    if (runlist != nullptr) tree = runlist->CopyTree("");
    
    if (chainHFS != nullptr) aptree = chainHFS->CloneTree(0);
    if (aptree) aptree->Branch("cc", &apcc);

    if (chainLTF != nullptr) varsLTF.set_tree(chainLTF);
    if (chainLRH != nullptr) varsLRH.set_tree(chainLRH);
    if (chainIIN != nullptr) varsIIN.set_tree(chainIIN);
    if (chainIEX != nullptr) varsIEX.set_tree(chainIEX);
    if (chainHEX != nullptr) varsHEX.set_tree(chainHEX);
    if (chainHFS != nullptr) varsHFS.set_tree(chainHFS);
    
    if (chainHEX != nullptr) varsHEX.set_tmva_reader("BDTG", "/afs/cern.ch/user/h/hchou/AMSProject/apflux/mva/MVAex/weights/TMVAClassification_BDTG.weights.xml");
    if (chainHFS != nullptr) varsHFS.set_tmva_reader("BDTG", "/afs/cern.ch/user/h/hchou/AMSProject/apflux/mva/MVAfs/weights/TMVAClassification_BDTG.weights.xml");
    //if (chainHEX != nullptr) varsHEX.set_tmva_reader("BDTG", "/afs/cern.ch/user/h/hchou/AMSProject/apflux/mva/test_MVAex/weights/TMVAClassification_BDTG.weights.xml");
    //if (chainHFS != nullptr) varsHFS.set_tmva_reader("BDTG", "/afs/cern.ch/user/h/hchou/AMSProject/apflux/mva/test_MVAfs/weights/TMVAClassification_BDTG.weights.xml");
}


void Analyzer::set_output(const std::string& outpath) {
    file = new TFile(outpath.c_str(), "RECREATE");

    if (file->IsZombie()) {
        std::cerr << Format("Error opening file: %s\n", outpath.c_str());
        LOG(FATAL) << Format("Error opening file: %s", outpath.c_str());
        file = nullptr;
        return;
    }

    std::string statement;
    statement += Format("\n****  Output Info ****\n");
    statement += Format("File : %s\n", outpath.c_str());
    std::cout << statement;
    LOG(INFO) << statement;
}


void Analyzer::write() {
    if (file == nullptr || tree == nullptr) return;
    file->cd();
    tree->Write();
    if (aptree) aptree->Write();
    file->Write();
}

void Analyzer::close() {
    if (file == nullptr || tree == nullptr) return;
    file->Close();
    file = nullptr;
    tree = nullptr;
    aptree = nullptr;
}

#endif // __Analyzer_H__
