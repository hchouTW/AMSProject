#ifndef __YiSlcNtuple_TCC__
#define __YiSlcNtuple_TCC__
#include "YiSlcNtuple.h"

//---- YiNtuple ----//
YiNtuple::MODE YiNtuple::eventMode = YiNtuple::ISS;

YiNtuple::YiNtuple() {
#if Debug == true
    std::cerr << "Debug : Now, YiNtuple::YiNtuple()\n";
#endif

    fRunChain = 0;
    fDataChain = 0;
    fFile = 0;
    init();

    fStopwatch.start();
}

YiNtuple::~YiNtuple() {
#if Debug == true
    std::cerr << "Debug : Now, YiNtuple::~YiNtuple()\n";
#endif
    init();

    fStopwatch.stop();
    fStopwatch.print();
}

inline void YiNtuple::init() {
#if Debug == true
    std::cerr << "Debug : Now, YiNtuple::init()\n";
#endif

    group.first = 0;
    group.second = -1;
    fileList.clear();
    if (fRunChain != 0) delete fRunChain;
    if (fDataChain != 0) delete fDataChain;
    if (fFile != 0) delete fFile;
    fRunChain = 0;
    fDataChain = 0;
    fFile = 0;
}

inline void YiNtuple::setOutputFile(const std::string& file_name, const std::string& path) {
#if Debug == true
    std::cerr << "Debug : Now, YiNtuple::setOutputFile()\n";
#endif

    std::string name = std::string(path) + "/" + std::string(file_name);
    fFile = new TFile(name.c_str(), "RECREATE");
    fFile->cd();
}

void YiNtuple::readDataFrom(const std::string& file_list, Long64_t group_th, Long64_t group_size) {
#if Debug == true
    std::cerr << "Debug : Now, YiNtuple::readDataFrom()\n";
#endif

    MGConfig::ReadList(file_list, group_th, group_size);
    COUT("\n**--------------------------------------------**\n");
    COUT("\n**    Read Data Form Source File List Info    **\n");
    COUT("\n**--------------------------------------------**\n");

    // start check sourceFileList.txt
    std::vector<std::string>&& flist = MGIO::ReadFileContent(file_list);
    if (flist.size() == 0)
        MGSys::ShowErrorAndExit(LOC_ADDR(), MGSys::Message("ROOT file list cannot be opend! Exiting ..."));
    // end check sourceFileList.txt

    // start load data with group
    if (group_th == 0 && group_size == -1) group_size = flist.size();
    if (group_size <= 0 || group_size > flist.size() || group_th < 0 || group_th >= flist.size())
        MGSys::ShowErrorAndExit(LOC_ADDR(), "Group format has error(1)! Exiting ...");

    Long64_t begin = group_th * group_size;
    Long64_t end   = (group_th + 1) * group_size;
    if (begin >= 0 && begin < flist.size() && end > flist.size()) {
        end = flist.size();
    }
    else if (begin < 0 || begin >= flist.size() || end < 1 || end > flist.size())
        MGSys::ShowErrorAndExit(LOC_ADDR(), "Group format has error(2)! Exiting ...");

    group = std::make_pair(group_th, group_size);
    for (int i = begin; i < end; i++) {
        fileList.push_back(flist.at(i).c_str());
    }

    COUT("\n---- Loading Root Files ----\n");
    COUT("Group : %ld th   [%ld files/group],    Total of Load Files : %ld \n", group.first, group.second, fileList.size());
    for (Long64_t  it = 0; it < fileList.size(); it++) {
        COUT("    Number : %ld,   %s\n", it, fileList.at(it).c_str());
    }
    // end load data with group

    // start read source file list
    fFile->cd();
    fRunChain = new TChain("runTag", "");
    fDataChain = new TChain("data", "");
	for (Long64_t it = 0; it < fileList.size(); it++) {
		fRunChain->Add(fileList.at(it).c_str());
		fDataChain->Add(fileList.at(it).c_str());
	}
	COUT("Totally : %ld data events.\n", fDataChain->GetEntries());
	// end read source file list

	COUT("\n**-------------------------------------------**\n");
	COUT("\n**    Read Data Form Source File List End    **\n");
	COUT("\n**-------------------------------------------**\n");
}

void YiNtuple::loopEventChain() {
	COUT("\n**-----------------------------**\n");
	COUT("\n**    Loop Event Chain Info    **\n");
	COUT("\n**-----------------------------**\n");

	// check event type
	if (fDataChain->GetEntries() <= 0)
        MGSys::ShowErrorAndExit(LOC_ADDR(), "Don't have event! Exiting ...");
	fFile->cd();

	setBranchAddress();
	analyzeEvent();

	fFile->cd();
	fFile->Write();
	fFile->Close();

	COUT("\n**----------------------------**\n");
	COUT("\n**    Loop Event Chain End    **\n");
	COUT("\n**----------------------------**\n");
}

#endif // __YiSlcNtuple_TCC__
