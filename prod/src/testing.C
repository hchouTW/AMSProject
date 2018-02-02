// ROOT library
#include <TString.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

// AMS library
#include <root.h>
#include <amschain.h>

// User defination library
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>

int main(int argc, const char ** argv) {
	COUT("\n\n");
	COUT("Usage : testing eventMode fileList begId endId\n");
	COUT("    Parameters : \n");
	COUT("    eventMode [ISS BT MC]\n");
	COUT("    fileList\n");
	COUT("    begId\n");
	COUT("    endId\n");
	COUT("\n\n");

	if (argc != 5 && argc != 6)
		MGSys::ShowErrorAndExit(LOC_ADDR(), "Number of argument is not conform! Exiting ...");

	std::string eventMode = argv[1];
	std::string fileList = argv[2];
	Long64_t begId = atol(argv[3]);
	Long64_t endId = atol(argv[4]);

    // Load Filelist
    bool stagedonly = true;
	unsigned int timeout = 10;
	AMSChain fChain("AMSRoot");
    int fileStatus = fChain.AddFromFile(fileList.c_str(), begId, endId, stagedonly, timeout);
    if (fileStatus == -1) 
        MGSys::ShowErrorAndExit(LOC_ADDR(), "ROOT file list cannot be opend! Exiting ..."); 
	COUT("FileStatus : %d\n", fileStatus);
    COUT("Totally : %ld data events.\n", fChain.GetEntries());
	
    // Event-Loop
    for (Long64_t ientry = 0; ientry < fChain.GetEntries(); ++ientry){
		AMSEventR* event = fChain.GetEvent(ientry);
        if (event == nullptr) continue;
        if (event->NTrTrack() != 1) continue;
    }

    return 1;
}
