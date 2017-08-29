#ifndef __YiAnaNtuple_C__
#define __YiAnaNtuple_C__

#include "YiSlcNtuple.h"
#include "YiSlcNtuple.tcc"

//#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jul27/src/ClassDef.h"
//#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jul27/src/ClassDef.C"
#include <vdev/src/ClassDef.h>
#include <vdev/src/ClassDef.C>

using namespace MGROOT;

// User defination macro
#ifdef Debug
    #undef Debug
#endif
#define Debug false

/*-------*/
/*  DST  */
/*-------*/
class DST {
    public :
        DST() : fDST(0) {}
        ~DST() {}

        void initDST(TFile * file, TChain * runChain = 0, TChain * dataChain = 0) {
            if (file == 0 || !gSaveDST) return;
            resetDST();
            file->cd();

            COUT("\n<<  init DST Info  >>\n");

            if (runChain != nullptr) {
                UInt_t runID = 0;
                UInt_t trgEV = 0;
                UInt_t selEV = 0;
                TTree * DSTRun = new TTree("runTag", "DST Run Info");
                DSTRun->Branch("runID", &runID);
                DSTRun->Branch("trgEV", &trgEV);
                DSTRun->Branch("selEV", &selEV);
                
                RunTagInfo * runTag = new RunTagInfo;
                runChain->SetBranchAddress("runTag", &runTag);
                for (Long64_t it = 0; it < runChain->GetEntries(); ++it) {
                    runChain->GetEntry(it);
                    runID = runTag->run;
                    trgEV = runTag->numOfTrgEvent;
                    selEV = runTag->numOfSelEvent;
                    DSTRun->Fill();
                }
                delete runTag;
            }

            if (dataChain != 0 && gSaveDSTClone) fDST = dataChain->CloneTree(0);
            if (fDST == 0      && gSaveDSTTree)  {
                fDST = new TTree("data", "MinDST data");
                //BranchMinDST(fDST, fMDst);
            }
            file->cd();
            
            Axis AXrso("Rigidity Resolution", 400, -1.5, 1.5);

            // rigidity binning
            AXnr = Axis("Rigidity [GV]", BinList( 
                    {   1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
                        3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
                        8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
                        19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
                        41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
                        93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } ) );
            AXir = Axis("1/Rigidity [1/GV]", AXnr, 1, true);
        
            // dx = 0.5 * (sqrt(x*x + 4 * dA / TWO_PI) - x)
            const Double_t AreaFT = 180.;
            std::vector<Double_t> radius; radius.push_back(0.);
            while (radius.at(radius.size()-1) < 52.) {
                Double_t r0 = radius.at(radius.size()-1);
                Double_t dr = 0.5 * (std::sqrt(r0*r0 + 2.*TMath::InvPi()*AreaFT) - r0);
                radius.push_back( (r0+dr) );
            }
            Axis AXCOMradius("Radius", radius);
            
            const Double_t CosFT = 0.03;
            std::vector<Double_t> aglcos; aglcos.push_back(1.);
            std::vector<Double_t> agldeg; agldeg.push_back(0.);
            while (agldeg.at(agldeg.size()-1) < 40.) {
                Double_t cos0 = aglcos.at(aglcos.size()-1);
                Double_t dcos = 0.5 * (std::sqrt(cos0*cos0 + 2.*TMath::InvPi()*CosFT) - cos0);
                Double_t cosv = (cos0 - dcos);
                Double_t degv = std::acos(cosv) * TMath::RadToDeg();
                aglcos.push_back(cosv);
                agldeg.push_back(degv);
            }
            Axis AXCOMangle("Angle", agldeg);
            
            COUT("\n<<  init DST End  >>\n");
            file->cd();
            return;
        }

        void resetDST() {
#if Debug == true
    std::cerr << "DST::resetDST()\n";
#endif
            //InitMinDST(fMDst);
        }

        void fillDST() {
#if Debug == true
    std::cerr << "DST::fillDST()\n";
#endif
            if (!gSaveDST) return;
            if (fDST == 0 || !gSaveDSTTree) return;
            fDST->Fill();
        }

        void finishDST() {
            Hist::Write();
        }

    public :
        TTree *  fDST;
        //MinDST   fMDst;
    
    public :
        static Bool_t gSaveDST;
        static Bool_t gSaveDSTTree;
        static Bool_t gSaveDSTClone;
        static UInt_t gUTime[2]; // (pre, cur)

    public :
        static Axis AXnr;
        static Axis AXir;
        static constexpr Float_t CfStableFT = 1.0; // 1.2
};

Bool_t DST::gSaveDST = true;
Bool_t DST::gSaveDSTTree = true;
Bool_t DST::gSaveDSTClone = false;
UInt_t DST::gUTime[2] = {0, 0};

Axis DST::AXnr;
Axis DST::AXir;


/*---------*/
/*  YiAna  */
/*---------*/
class YiAna : public YiNtuple, public DST {
    public :
        YiAna();
        ~YiAna();

        void setBranchAddress();
        void freeBranchAddress();
        void analyzeEvent();
        bool analyzeRunInfo();
        bool analyzeCore();

    public :
        LIST * fList;
        G4MC * fG4mc;
        RTI  * fRti;
        TRG  * fTrg;
        TOF  * fTof;
        ACC  * fAcc;
        TRK  * fTrk;
        TRD  * fTrd;
        RICH * fRich;
        ECAL * fEcal;
};

YiAna::YiAna() {
    DST::fDST = 0;
    YiNtuple::fRunChain = 0;
    YiNtuple::fDataChain = 0;
    YiNtuple::init();
    YiNtuple::fStopwatch.start();
    fList = 0;
    fG4mc = 0;
    fRti  = 0;
    fTrg  = 0;
    fTof  = 0;
    fAcc  = 0;
    fTrk  = 0;
    fTrd  = 0;
    fRich = 0;
    fEcal = 0;
}

YiAna::~YiAna() {
    freeBranchAddress();
}


void YiAna::setBranchAddress() {
    COUT("\n<<  Set Branch Address Info  >>\n");

    /*******************/
    /**  USER DEFINE  **/
    /*******************/
    freeBranchAddress();
    fList = new LIST;
    fG4mc = (YiNtuple::CheckEventMode(YiNtuple::MC)) ? (new G4MC) : 0;
    fRti  = (YiNtuple::CheckEventMode(YiNtuple::ISS)) ? (new RTI) : 0;
    fTrg  = new TRG;
    fTof  = new TOF;  
    fAcc  = new ACC;  
    fTrk  = new TRK;  
    fTrd  = new TRD;  
    fRich = new RICH;
    fEcal = new ECAL;

    fDataChain->SetBranchAddress("list", &fList);
    if (YiNtuple::CheckEventMode(YiNtuple::MC))
        fDataChain->SetBranchAddress("g4mc", &fG4mc);
    if (YiNtuple::CheckEventMode(YiNtuple::ISS))
        fDataChain->SetBranchAddress("rti",  &fRti);
    fDataChain->SetBranchAddress("trg",  &fTrg);
    fDataChain->SetBranchAddress("tof",  &fTof);
    fDataChain->SetBranchAddress("acc",  &fAcc);
    fDataChain->SetBranchAddress("trk",  &fTrk);
    fDataChain->SetBranchAddress("trd",  &fTrd);
    fDataChain->SetBranchAddress("rich", &fRich);
    fDataChain->SetBranchAddress("ecal", &fEcal);
    /*******************/
    /**               **/
    /*******************/

    COUT("\n<<  Set Branch Address End  >>\n");
}

void YiAna::freeBranchAddress() {
    if (fList != 0) { delete fList; fList = 0; }
    if (fG4mc != 0) { delete fG4mc; fG4mc = 0; }
    if (fRti  != 0) { delete fRti ; fRti  = 0; }
    if (fTrg  != 0) { delete fTrg ; fTrg  = 0; }
    if (fTof  != 0) { delete fTof ; fTof  = 0; }
    if (fAcc  != 0) { delete fAcc ; fAcc  = 0; }
    if (fTrk  != 0) { delete fTrk ; fTrk  = 0; }
    if (fTrd  != 0) { delete fTrd ; fTrd  = 0; }
    if (fRich != 0) { delete fRich; fRich = 0; }
    if (fEcal != 0) { delete fEcal; fEcal = 0; }
}

void YiAna::analyzeEvent() {
    COUT("\n<<  Analyze Event Info  >>\n");

    /*******************/
    /**  USER DEFINE  **/
    /*******************/
    initDST(fFile, fRunChain, fDataChain);
    fFile->cd();
    /*******************/
    /**               **/
    /*******************/

    Long64_t nentries = fDataChain->GetEntries();
    Long64_t npassed = 0;
    Long64_t nprocessed = 0;
    Long64_t printRate = nentries / 200;
    if (printRate < 250000) printRate = 250000;
    if (printRate > 2500000) printRate = 2500000;

    for (Long64_t ientry = 0; ientry < nentries; ientry++) {
        if (nprocessed%printRate == 0) {
            const UInt_t MemSize = 1024;
            ProcInfo_t procinfo;
            gSystem->GetProcInfo(&procinfo);
            Long64_t memRes = procinfo.fMemResident / MemSize;
            Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
            fStopwatch.stop();

            COUT("Info :: %lf %\n", 100. * float(nprocessed)/float(nentries));
            COUT("        Processed       : %ld / %ld\n", nprocessed, nentries);
            COUT("        Passed          : %ld / %ld\n", npassed, nprocessed);
            COUT("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
            COUT("        Real Time       : %9.2f (second)\n", fStopwatch.time());
            COUT("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fStopwatch.time());
            COUT("        Cpu    System   : %4.1f %\n", procinfo.fCpuSys);
            COUT("               User     : %4.1f %\n", procinfo.fCpuUser);
            COUT("        Memory Resident : %2ld GB %4ld MB\n", memRes / MemSize, memRes % MemSize);
            COUT("               Virtual  : %2ld GB %4ld MB\n", memVrl / MemSize, memVrl % MemSize);
        }
        nprocessed++;

        fDataChain->GetEntry(ientry);

        /*******************/
        /**  USER DEFINE  **/
        /*******************/
        resetDST();
        //if (!analyzeRunInfo()) continue;
        if (!analyzeCore()) continue;
        fillDST();
        /*******************/
        /**               **/
        /*******************/

        npassed++;
    }

    /*******************/
    /**  USER DEFINE  **/
    /*******************/
    finishDST();
    /*******************/
    /**               **/
    /*******************/

    if (nprocessed == nentries) {
        const UInt_t MemSize = 1024;
        ProcInfo_t procinfo;
        gSystem->GetProcInfo(&procinfo);
        Long64_t memRes = procinfo.fMemResident / MemSize;
        Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
        fStopwatch.stop();
        
        COUT("Info :: %lf %\n", 100. * float(nprocessed)/float(nentries));
        COUT("        Processed       : %ld / %ld\n", nprocessed, nentries);
        COUT("        Passed          : %ld / %ld\n", npassed, nprocessed);
        COUT("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
        COUT("        Real Time       : %9.2f (second)\n", fStopwatch.time());
        COUT("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fStopwatch.time());
        COUT("        Cpu    System   : %4.1f %\n", procinfo.fCpuSys);
        COUT("               User     : %4.1f %\n", procinfo.fCpuUser);
        COUT("        Memory Resident : %2ld GB %4ld MB\n", memRes / MemSize, memRes % MemSize);
        COUT("               Virtual  : %2ld GB %4ld MB\n", memVrl / MemSize, memVrl % MemSize);
        COUT("Info :: Root Files Processed Successfully Finished.\n");
    }
    else {
        COUT("Info :: Root Files Processed Seems Failed.\n");
        COUT("        Processed %ld in %ld\n", nprocessed, nentries);
    }

    COUT("\n<<  Analyze Event End  >>\n");
}


/*******************/
/**  USER DEFINE  **/
/*******************/
bool YiAna::analyzeRunInfo() {
#if Debug == true
    std::cerr << "YiAna::analyzeRunInfo()\n";
#endif

    if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
        // preselection (Mini-RTI Requirement)
        if (!fRti->flagRun) return false;
        if (!fRti->isGoodSecond) return false;
        if (fRti->isInSAA) return false;
        if (fRti->liveTime < 0.5) return false;

        // Depend on Events
        if (fRti->isInShadow) return false;
    }

    return true;
}
/*******************/
/**               **/
/*******************/


/*******************/
/**  USER DEFINE  **/
/*******************/
bool YiAna::analyzeCore() {
#if Debug == true
    std::cerr << "YiAna::analyzeCore()\n";
#endif
    
    return true;
}
/*******************/
/**               **/
/*******************/


/*-----------------*/
/*  Main Function  */
/*-----------------*/
int main(int argc, const char ** argv) {
    COUT("\n**--------------------------**\n");
    COUT("\n**    YiAnaNtuple START     **\n");
    COUT("\n**--------------------------**\n");
    MGROOT::Style::LoadDefaultEnvironment();

    COUT("\n\n");
    COUT("Usage : YiAnaNtuple event_mode file_list group_th group_size (path)\n");
    COUT("    Parameters : \n");
    COUT("    event_mode [ISS BT MC]\n");
    COUT("    file_list\n");
    COUT("    group_th\n");
    COUT("    group_size\n");
    COUT("    (path)\n");
    COUT("\n\n");

    if (argc != 5 && argc != 6)
        MGSys::ShowErrorAndExit(LOC_ADDR(), MGSys::Message("Number of argument is not conform! Exiting ..."));

    std::string event_mode = argv[1];
    std::string file_list = argv[2];
    Long64_t group_th = atol(argv[3]);
    Long64_t group_size = atol(argv[4]);

    std::use_facet<std::ctype<char> >(std::locale()).toupper(&event_mode[0], &event_mode[0] + event_mode.size());
    if (event_mode == "ISS") YiNtuple::SetEventMode(YiNtuple::ISS);
    else if (event_mode == "BT") YiNtuple::SetEventMode(YiNtuple::BT);
    else if (event_mode == "MC") YiNtuple::SetEventMode(YiNtuple::MC);
    else MGSys::ShowErrorAndExit(LOC_ADDR(), MGSys::Message("Can't find event mode (ISS, BT, MC)! Exiting ..."));

    std::string outputFile = STR_FMT("YiAnalytics_%s.%07ld.root", event_mode.c_str(), group_th);
    std::string path = ".";
    if (argc == 6) path = argv[5];
    
    YiAna ntuple;
    ntuple.setOutputFile(outputFile.c_str(), path.c_str());
    ntuple.readDataFrom(file_list.c_str(), group_th, group_size);
    ntuple.loopEventChain();

    //YiAna * ntuple = new YiAna();
    //ntuple->setOutputFile(outputFile.c_str(), path.c_str());
    //ntuple->readDataFrom(file_list.c_str(), group_th, group_size);
    //ntuple->loopEventChain();
    //if (ntuple != nullptr) delete ntuple;
    //ntuple = nullptr;

    COUT("\n**------------------------**\n");
    COUT("\n**    YiAnaNtuple END     **\n");
    COUT("\n**------------------------**\n");
    return 0;
}
#endif // __YiAnaNtuple_C__
