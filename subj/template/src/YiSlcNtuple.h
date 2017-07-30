#ifndef __YiSlcNtuple_H__
#define __YiSlcNtuple_H__

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
#include <map>

// Root library
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile2D.h>
#include <TSystem.h>

// User defination library
#include "CPPLibs/CPPLibs.h"
#include "ROOTLibs/ROOTLibs.h"
#include "TRACKLibs/TRACKLibs.h"

// User defination macro
#define Debug true

//---- YiNtuple ----//
class YiNtuple {
	public :
		enum MODE {
			ISS, BT, MC
		};
		static enum MODE eventMode;
		static void SetEventMode(YiNtuple::MODE mode) { YiNtuple::eventMode = mode; }
		static bool CheckEventMode(YiNtuple::MODE mode) { return (YiNtuple::eventMode == mode); }

	public :
		YiNtuple();
		~YiNtuple();

		void init();
		void setOutputFile(const std::string& file_name = "YiNtuple.root", const std::string& path = ".");
		void readDataFrom(const std::string& file_list = "fileList.txt", Long64_t group_th = 0, Long64_t group_size = -1);
		void loopEventChain();

		virtual void setBranchAddress() = 0;
		virtual void analyzeEvent() = 0;

	protected :
		std::pair<long, long> group;
		std::vector<std::string> fileList;

		TChain * fRunChain;
		TChain * fDataChain;
		TFile * fFile;

	protected :
		MGClock::HrsStopwatch fStopwatch;
};


#endif // __YiSlcNtuple_H__
