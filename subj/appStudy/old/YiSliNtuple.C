#ifndef __YiAnaNtuple_C__
#define __YiAnaNtuple_C__

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.h"
#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.tcc"

#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jun05/src/ClassDef.h"
#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jun05/src/ClassDef.C"

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/src/minDST.h"

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
                BranchMinDST(fDST, fMDst);
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
        
            //TString sbnPH = "/afs/cern.ch/user/h/hchou/public/DATABASE/physics/binning/phe_bin2.root";
            //TFile * fPHBin = TFile::Open(sbnPH);
            //TH1D  * hPHBn0 = (TH1D*)fPHBin->Get("hist2");
            //Axis AXPHnr("Rigidity [GV]", hPHBn0, Axis::kX);
            //Axis AXPHir = MgntROOT::Axis::Invert("Inverse Rigidity [1/GV]", AXPHnr);
            //Axis AXnr = AXPHnr;
            //Axis AXir = AXPHir;
            //fPHBin->Close();
            //file->cd();
            
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
            InitMinDST(fMDst);
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
        MinDST   fMDst;
    
    public :
        static Bool_t gSaveDST;
        static Bool_t gSaveDSTTree;
        static Bool_t gSaveDSTClone;
        static UInt_t gUTime[2]; // (pre, cur)

    public :
        static Axis AXnr;
        static Axis AXir;
        static const Float_t CfStableFT = 1.0; // 1.2

    public :
        static Float_t GetStdMDR(Int_t ipatt);
        static Float_t GetStdSqrMDR(Int_t ipatt);
        static Float_t GetMDR(Int_t ipatt);
        static Float_t GetSqrMDR(Int_t ipatt);
        static Float_t GetRigSgm(Int_t ipatt, Float_t rig, Float_t mass);
};

Bool_t DST::gSaveDST = true;
Bool_t DST::gSaveDSTTree = true;
Bool_t DST::gSaveDSTClone = false;
UInt_t DST::gUTime[2] = {0, 0};

Axis DST::AXnr;
Axis DST::AXir;

Float_t DST::GetStdMDR(Int_t ipatt) {
    Float_t MDR = 1800.; // (GeV)
    if      (ipatt == 2) MDR =  240.;
    if      (ipatt == 3) MDR =  540.;
    if      (ipatt == 4) MDR =  750.;
    if      (ipatt == 5) MDR = 1800.;
    return MDR;
}

Float_t DST::GetStdSqrMDR(Int_t ipatt) {
    Float_t SqrMDR = 3.24; // (TV)
    if      (ipatt == 2) SqrMDR = 5.760000e-02;
    else if (ipatt == 3) SqrMDR = 2.916000e-01;
    else if (ipatt == 4) SqrMDR = 5.625000e-01;
    else if (ipatt == 5) SqrMDR = 3.240000e+00;
    return SqrMDR;
}

Float_t DST::GetMDR(Int_t ipatt) {
    Float_t MDR = 2.0; // (TV)
    if      (ipatt == 0) MDR = 1.398205e-01;
    else if (ipatt == 1) MDR = 1.619846e-01;
    else if (ipatt == 2) MDR = 2.885081e-01;
    else if (ipatt == 3) MDR = 8.571684e-01;
    else if (ipatt == 4) MDR = 1.079867e+00;
    else if (ipatt == 5) MDR = 2.000000e+00;
    return MDR;
}

Float_t DST::GetSqrMDR(Int_t ipatt) {
    Float_t SqrMDR = 4.0; // (TV)
    if      (ipatt == 0) SqrMDR = 1.398205e-01 * 1.398205e-01;
    else if (ipatt == 1) SqrMDR = 1.619846e-01 * 1.619846e-01;
    else if (ipatt == 2) SqrMDR = 2.885081e-01 * 2.885081e-01;
    else if (ipatt == 3) SqrMDR = 8.571684e-01 * 8.571684e-01;
    else if (ipatt == 4) SqrMDR = 1.079867e+00 * 1.079867e+00;
    else if (ipatt == 5) SqrMDR = 2.000000e+00 * 2.000000e+00;
    return SqrMDR;
}

Float_t DST::GetRigSgm(Int_t ipatt, Float_t nrig, Float_t mass) {
    Float_t trSgm = 1.0;
    if      (ipatt == 0) {
        Float_t trParIU[3] = { 1.44615e-02, 1.58178e-03, 5.11515e-05 };
        Float_t trSgmIU = std::sqrt(trParIU[0] * (1.0 + mass*mass/nrig/nrig) + trParIU[1] * nrig + trParIU[2] * nrig * nrig) / nrig;
        trSgm = trSgmIU;
    }
    else if (ipatt == 1) {
        Float_t trParIL[3] = { 9.60391e-03, 7.54520e-04, 3.81112e-05 };
        Float_t trSgmIL = std::sqrt(trParIL[0] * (1.0 + mass*mass/nrig/nrig) + trParIL[1] * nrig + trParIL[2] * nrig * nrig) / nrig;
        trSgm = trSgmIL;
    }
    else if (ipatt == 2) {
        Float_t trParIn[3] = { 8.39123e-03, 3.36288e-04, 1.20139e-05 };
        Float_t trSgmIn = std::sqrt(trParIn[0] * (1.0 + mass*mass/nrig/nrig) + trParIn[1] * nrig + trParIn[2] * nrig * nrig) / nrig;
        trSgm = trSgmIn;
    }
    else if (ipatt == 3) {
        Float_t trParL1[3] = { 7.64252e-03, 5.34561e-04, 1.36103e-06 };
        Float_t trSgmL1 = std::sqrt(trParL1[0] * (1.0 + mass*mass/nrig/nrig) + trParL1[1] * nrig + trParL1[2] * nrig * nrig) / nrig;
        trSgm = trSgmL1;
    }
    else if (ipatt == 4) {
        Float_t trParL9[3] = { 8.23647e-03, 2.97335e-04, 8.57550e-07 };
        Float_t trSgmL9 = std::sqrt(trParL9[0] * (1.0 + mass*mass/nrig/nrig) + trParL9[1] * nrig + trParL9[2] * nrig * nrig) / nrig;
        trSgm = trSgmL9;
    }
    else if (ipatt == 5) {
        Float_t trParFs[3] = { 9.29537e-03, 1.14609e-04, 2.50000e-07 };
        Float_t trSgmFs = std::sqrt(trParFs[0] * (1.0 + mass*mass/nrig/nrig) + trParFs[1] * nrig + trParFs[2] * nrig * nrig) / nrig;
        trSgm = trSgmFs;
    }
    return trSgm;
}


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
    
    Long64_t runID   = fList->run;
    Long64_t eventID = fList->event;
    Long64_t entryID = fList->entry;
    Double_t weight  = fList->weight;

    // Monte Carlo
    Float_t MCTRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
                     0.0 : (fG4mc->primPart.mom / fG4mc->primPart.chrg);
    Float_t MCNRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
                     0.0 : std::fabs(MCTRig);
    Float_t MCIRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ?
                     0.0 : (1. / MCTRig);
    Int_t   MCSign = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ?
                     0 : MGNumc::Compare(MCTRig);
    
    VertexMCInfo * mcVtx = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
                           nullptr : &fG4mc->primVtx;

    Float_t MCFrac = -1;
    if (mcVtx != nullptr && fG4mc->secParts.size() > 0) {
        PartMCInfo * cand = nullptr;
        for (auto&& part : fG4mc->secParts) {
            if (cand == nullptr) cand = &part;
            if (MGNumc::Compare(cand->mom, part.mom) >= 0) continue;
            cand = &part;
        }
        MCFrac = (cand->mom / fG4mc->primPart.mom);
    }


#if Debug == true
    std::cerr << "YiAna::analyzeCore()   [][] COM Selection [][]\n";
#endif
    
    // preselection (Mini-RTI Requirement)
    if (!YiAna::analyzeRunInfo()) return false;
    
    // preselection (Mini-Track Requirement)
    TrackInfo * track = (fTrk->tracks.size() == 1) ? (&fTrk->tracks.at(0)) : nullptr;
    if (track == nullptr) return false;
    if ((track->bitPatt&1)!=1 || !track->status[0][0]) return false; // Inner XY
    
    // preselection (Min-Track Charge Requirement) 
    Bool_t trQL2Cut = ((MGNumc::Compare(track->QL2) > 0) && (track->QL2 < 0.75));
    Bool_t trQL1Cut = ((MGNumc::Compare(track->QL1) > 0) && (track->QL1 < 0.75));
    Bool_t trQL9Cut = ((MGNumc::Compare(track->QL9) > 0) && (track->QL9 < 0.75));
    if (trQL2Cut) return false;
    if (trQL1Cut) return false;
    if (trQL9Cut) return false;

    // preselection (Min-Track Pattern Requirement)
    Bool_t isTrIn = ((track->bitPatt&  5)==  5 && 
                     track->status[0][0]);                         // Inner Y + L2             (  5 = 1 + 4)
    Bool_t isTrL1 = ((track->bitPatt& 37)== 37 && 
                     track->status[0][0] && track->status[0][1]);  // Inner Y + L2 + L1 XY     ( 37 = 1 + 4 + 32)
    Bool_t isTrL9 = ((track->bitPatt&133)==133 && 
                     track->status[0][0] && track->status[0][2]);  // Inner Y + L2 + L9 XY     (133 = 1 + 4 + 128)
    Bool_t isTrFs = ((track->bitPatt&161)==161 && 
                     track->status[0][3] && 
                     track->status[0][1] && track->status[0][2]);  // Inner Y + L1 XY + L9 XY  (161 = 1 + 32 + 128)
    Bool_t isTrPt = (isTrIn || isTrL1 || isTrL9 || isTrFs);
    if (!isTrPt) return false;

    // Track
    const Float_t massPr = 0.938272297; 
    const Float_t massPi = 0.139570180; 
    Float_t trigIn   = track->rigidity[0][0];
    Float_t nrigIn   = std::fabs(track->rigidity[0][0]);
    Float_t trXY[2]  = { track->state[0][0][0], track->state[0][0][1] };
    Float_t trRad    = std::sqrt(trXY[0]*trXY[0]+trXY[1]*trXY[1]);
    Float_t trAgl    = std::acos( std::fabs(track->stateLJ[0][0][0][5]) ) * TMath::RadToDeg();
    Float_t trMassPr = (massPr/nrigIn)*(massPr/nrigIn);
    Float_t trMassPi = (massPi/nrigIn)*(massPi/nrigIn);
    
    Short_t trSign[4] = {0, 0, 0, 0};
    Float_t trNRig[4] = {0, 0, 0, 0};
    if (isTrIn) { trSign[0] = MGNumc::Compare(track->rigidity[0][0]); trNRig[0] = std::fabs(track->rigidity[0][0]); }
    if (isTrL1) { trSign[1] = MGNumc::Compare(track->rigidity[0][1]); trNRig[1] = std::fabs(track->rigidity[0][1]); }
    if (isTrL9) { trSign[2] = MGNumc::Compare(track->rigidity[0][2]); trNRig[2] = std::fabs(track->rigidity[0][2]); }
    if (isTrFs) { trSign[3] = MGNumc::Compare(track->rigidity[0][3]); trNRig[3] = std::fabs(track->rigidity[0][3]); }
    
    Short_t     trPt = 0;
    std::string trNm = "";
    if      (isTrFs) { trPt = 3; trNm = "Fs"; }
    else if (isTrL9) { trPt = 2; trNm = "L9"; }
    else if (isTrL1) { trPt = 1; trNm = "L1"; }
    else if (isTrIn) { trPt = 0; trNm = "In"; }
    
    Float_t trig = track->rigidity[0][trPt];
    Int_t   sign = MGNumc::Compare(trig);
    Float_t nrig = std::fabs(trig);
    Float_t irig = 1. / trig;
    
    // preselection (IGRF Cutoff Requirement) 
    Int_t   patAgl = -1;
    if      (trAgl < 25.0) patAgl = 0; 
    else if (trAgl < 30.0) patAgl = 1; 
    else if (trAgl < 35.0) patAgl = 2; 
    else if (trAgl < 40.0) patAgl = 3;
    else return false;
    
    // preselection (Mini-FriducialVolume Requirement)
    const Float_t trFiducialVolume = 47.0;
    if (trRad > trFiducialVolume) return false;
    
    // preselection (Mini-Acc Requirement)
    if (fAcc->clusters.size() != 0) return false;
    
    // preselection (Mini-Trigger Requirement)
    if ((fTrg->bit&8) != 8) return false;

    // preselection (Min-Track Charge Requirement) 
    if (track->QIn < 0.8 || track->QIn > 1.4) return false;
    
    // preselection (Mini-Tracker External Hit Requirement) 
    Short_t trNHL1 = 0;
    Short_t trNHL9 = 0;
    for (auto&& hit : fTrk->otherHits) {
        if (hit.side != 3) continue;
        if (hit.layJ != 1 && hit.layJ != 9) continue;
        Float_t chrg = 0.5 * (hit.chrg[0] + hit.chrg[1]);
        if (chrg < 0.75) continue;
        if      (hit.layJ == 1) trNHL1++;
        else if (hit.layJ == 9) trNHL9++;
    }
    //if (trNHL1 >= 2) return false;
    //if (trNHL9 >= 1) return false;
    //if (trNHL9 >= 2) return false; // testcode

    // preselection (Mini-TOF Requirement)
    if (fTof->numOfBetaH != 1) return false;
    if (!fTof->statusBetaH) return false;
    if (fTof->betaHPatt != 15) return false;
    if (fTof->betaHGoodTime != 15) return false;
    if (fTof->betaH < 0.3 || fTof->betaH > 1.3) return false;

    Float_t tofQ    = fTof->Qall;
    Float_t tofQu   = 0.50 * (fTof->Q[0] + fTof->Q[1]);
    Float_t tofQl   = 0.50 * (fTof->Q[2] + fTof->Q[3]);
    
    // preselection (Min-TOF Charge Requirement) 
    if (tofQ  < 0.8) return false;
    if (tofQu > 1.5) return false;
    if (tofQl > 2.0) return false;
    
    // preselection (Min-TOF Diff-Charge Requirement) 
    Float_t tofDQu = std::fabs(fTof->Q[1] - fTof->Q[0]);
    Float_t tofDQl = std::fabs(fTof->Q[3] - fTof->Q[2]);
    if (tofDQu > 0.85) return false;
    if (tofDQl > 0.85) return false;

    // preselection (Min-TOF Quality Requirement) 
    if (fTof->normChisqT > 10) return false;
    if (fTof->normChisqC > 10) return false;
    if (fTof->numOfInTimeCluster > 4) return false;
   
    // preselection (Mini-TOF External Hit Requirement)
    if ((fTof->extClsN[0] + fTof->extClsN[1]) > 1) return false;

    Float_t tofMSgm = 0.08;
    Float_t tofM    = ((1.0/fTof->betaH/fTof->betaH - 1.0) - trMassPr) / tofMSgm;

    // preselection (Mini-TRD Track && Vertex Requirement)
    Bool_t isTrdVtx = (fTrd->vtxSide == 15 && fTrd->vtxCoo[2] > 50 && fTrd->vtxCoo[2] < 200);
    //Bool_t isCleanTrdTrack = ((fTrd->numOfTrack == 1 || fTrd->numOfHTrack == 1) && !isTrdVtx);
    Bool_t isCleanTrdTrack = ((fTrd->numOfTrack == 1 || fTrd->numOfHTrack == 1));
    if (!isCleanTrdTrack) return false;
    
    // preselection (Mini-TRD Estimator Requirement)
    Int_t trdEstPt = -1; // trd 0, trk 1
    if      (fTrd->statusKCls[0] && fTrd->LLR_nhit[0] >= 8) trdEstPt = 0;
    else if (fTrd->statusKCls[1] && fTrd->LLR_nhit[1] >= 8) trdEstPt = 1;
    Bool_t  hasTrdEst = (trdEstPt != -1);
    Bool_t  istrdHe   = (hasTrdEst) ? (fTrd->LLR[trdEstPt][2] > 0.3) : false;
    Float_t trdEst    = (hasTrdEst) ? (fTrd->LLR[trdEstPt][0]) : -1.0;
    if (!hasTrdEst) return false;
    if (istrdHe) return false;
    
    // preselection (Mini-Shower Requirement)
    ShowerInfo * shower = (fEcal->showers.size() >= 1) ? (&fEcal->showers.at(0)) : nullptr;
    Bool_t  hasEcal  = (shower != nullptr);
    Float_t ecalbdt  = (!hasEcal) ? -2.0  : shower->PisaBDT;
    
    // preselection (Mini-Rich Requirement)
    Bool_t  hasRich    = (fRich->status && fRich->kindOfRad != -1);
    Float_t richMSgm   = (fRich->kindOfRad == 0) ? 0.0025 : 0.0065;
    Float_t richMCrrPr = richMSgm * (1.0 + 100.0*(std::sqrt(1.0+trMassPr) - 1.0));
    Float_t richMCrrPi = richMSgm * (1.0 + 100.0*(std::sqrt(1.0+trMassPi) - 1.0));
    Float_t richMPr = (!hasRich) ? 0.0 : 
                      ((1.0/fRich->beta/fRich->beta - 1.0) - trMassPr) / richMCrrPr;
    Float_t richMPi = (!hasRich) ? 0.0 : 
                      ((1.0/fRich->beta/fRich->beta - 1.0) - trMassPi) / richMCrrPi;

    // preselection (Mini-Rich Hits Requirement)
    if ((fRich->numOfCrossHit[0] + fRich->numOfCrossHit[1]) >= 15) return false; // Particle Hits >= 15
    //if (fRich->numOfRingHit[0][2] >=  8) return false; // Outside Ring Hits >= 8
    //if (fRich->numOfRingHit[3][1] >= 12) return false; // Non-Selected Ring Hits >= 12
    
    
    // Charge Confusion ---- Chisq
    Float_t lchixIn = (track->status[0][0]) ? std::log(track->chisq[0][0][0]) : 0.0;
    Float_t lchixL1 = (track->status[0][1]) ? std::log(track->chisq[0][1][0]) : 0.0;
    Float_t lchixL9 = (track->status[0][2]) ? std::log(track->chisq[0][2][0]) : 0.0;
    Float_t lchixFs = (track->status[0][3]) ? std::log(track->chisq[0][3][0]) : 0.0;
    
    Float_t lchiyIn = (track->status[0][0]) ? std::log(track->chisq[0][0][1]) : 0.0;
    Float_t lchiyL1 = (track->status[0][1]) ? std::log(track->chisq[0][1][1]) : 0.0;
    Float_t lchiyL9 = (track->status[0][2]) ? std::log(track->chisq[0][2][1]) : 0.0;
    Float_t lchiyFs = (track->status[0][3]) ? std::log(track->chisq[0][3][1]) : 0.0;
    
    // Charge Confusion ---- AsymRig
    const Float_t lasymLMT = 1.0e-8;
    const Float_t lasymSGM = 1.8;
    Float_t trSgmIn = GetRigSgm(2, nrig, massPr);
    Float_t trSgmL1 = GetRigSgm(3, nrig, massPr);
    Float_t trSgmL9 = GetRigSgm(4, nrig, massPr);
    Float_t trSgmFs = GetRigSgm(5, nrig, massPr);

    Float_t crrPar1I[5] = { 4.93251e-01, 1.90417e+00, 9.18502e+01, 4.86389e+00, 3.29747e-01 };
    Float_t crrSgm1I    = crrPar1I[0]*std::erf(crrPar1I[1]*(std::log(nrig+crrPar1I[2])-crrPar1I[3]))+crrPar1I[4];
    Float_t trSgm1I     = std::sqrt(trSgmIn * trSgmIn + trSgmL1 * trSgmL1) * crrSgm1I;
    Float_t asym1I      = (1.0/track->rigidity[0][0] - 1.0/track->rigidity[0][1]) / trSgm1I; 
    Float_t lasym1I     = (track->status[0][0] && track->status[0][1]) ? std::log(asym1I * asym1I + lasymLMT) / lasymSGM : 0.0;
    
    Float_t crrPar9I[5] = { 4.17294e-01, 1.15603e+00, 1.62194e+01, 3.81780e+00, 3.92404e-01 };
    Float_t crrSgm9I    = crrPar9I[0]*std::erf(crrPar9I[1]*(std::log(nrig+crrPar9I[2])-crrPar9I[3]))+crrPar9I[4];
    Float_t trSgm9I     = std::sqrt(trSgmIn * trSgmIn + trSgmL9 * trSgmL9) * crrSgm9I;
    Float_t asym9I      = (1.0/track->rigidity[0][0] - 1.0/track->rigidity[0][2]) / trSgm9I;
    Float_t lasym9I     = (track->status[0][0] && track->status[0][2]) ? std::log(asym9I * asym9I + lasymLMT) / lasymSGM : 0.0;

    Float_t crrPar91[5] = { 5.25352e-01, 9.88059e-01, 1.21200e+01, 3.72027e+00, 4.99424e-01 };
    Float_t crrSgm91    = crrPar91[0]*std::erf(crrPar91[1]*(std::log(nrig+crrPar91[2])-crrPar91[3]))+crrPar91[4];
    Float_t trSgm91     = std::sqrt(trSgmL1 * trSgmL1 + trSgmL9 * trSgmL9) * crrSgm91;
    Float_t asym91      = (1.0/track->rigidity[0][1] - 1.0/track->rigidity[0][2]) / trSgm91;
    Float_t lasym91     = (track->status[0][1] && track->status[0][2]) ? std::log(asym91 * asym91 + lasymLMT) / lasymSGM : 0.0;

    // MinDST
    fMDst.run    = runID;
    fMDst.event  = eventID;
    fMDst.weight = weight;

    if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
        fMDst.uTime = fRti->uTime;
        fMDst.cfRig = fRti->cutoffIGRF[patAgl];
    }

    if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
        fMDst.mcSign = MCSign;
        fMDst.mcNRig = MCNRig;
        fMDst.mcInt  = mcVtx->status;
        std::copy(mcVtx->coo, mcVtx->coo+3, fMDst.mcVtx);
        fMDst.mcFrac = MCFrac;
    }
    
    fMDst.trPt = trPt;
    fMDst.trNHL1  = trNHL1;
    fMDst.trNHL9  = trNHL9;
    std::copy(trSign, trSign+4, fMDst.sign);
    std::copy(trNRig, trNRig+4, fMDst.nrig);
    
    fMDst.lchix[0] = lchixIn;
    fMDst.lchix[1] = lchixL1;
    fMDst.lchix[2] = lchixL9;
    fMDst.lchix[3] = lchixFs;
    
    fMDst.lchiy[0] = lchiyIn;
    fMDst.lchiy[1] = lchiyL1;
    fMDst.lchiy[2] = lchiyL9;
    fMDst.lchiy[3] = lchiyFs;
    
    fMDst.lasym1I = lasym1I;
    fMDst.lasym9I = lasym9I;
    fMDst.lasym91 = lasym91;

    fMDst.tofM = tofM;
    for (Int_t it = 0; it < 2; ++it) {
        fMDst.tofCN[it] = fTof->extClsN[it];
        for (Int_t jt = 0; jt < 2; ++jt) {
            fMDst.tofCL[it][jt] = fTof->extClsL[it][jt];
            if (fTof->extClsL[it][jt] == 0) continue;
            fMDst.tofCQ[it][jt] = fTof->extClsQ[it][jt];
            fMDst.tofCT[it][jt] = (fTof->extClsQ[it][jt] < 0.9) ? -10 : fTof->extClsT[it][jt];
        }
    }


    fMDst.trdEst = trdEst;
    fMDst.trdVtx = isTrdVtx;
    fMDst.trdSeg[0] = fTrd->numOfHSeg[0];  
    fMDst.trdSeg[1] = fTrd->numOfHSeg[1];

    if (hasEcal) {
        fMDst.hasEcal = hasEcal;
        fMDst.ecalEst = ecalbdt;
    }

    fMDst.richRad  = fRich->kindOfRad;
    fMDst.richPhEl = fRich->numOfExpPE[0];
    fMDst.richPhPr = fRich->numOfExpPE[3];
    fMDst.richCh[0] = fRich->numOfCrossHit[0];
    fMDst.richCh[1] = fRich->numOfCrossHit[1];
    fMDst.richRhEl[0] = fRich->numOfRingHit[0][0];
    fMDst.richRhEl[1] = fRich->numOfRingHit[0][1];
    fMDst.richRhEl[2] = fRich->numOfRingHit[0][2];
    fMDst.richRhPr[0] = fRich->numOfRingHit[3][0];
    fMDst.richRhPr[1] = fRich->numOfRingHit[3][1];
    fMDst.richRhPr[2] = fRich->numOfRingHit[3][2];
    if (hasRich) {
        fMDst.hasRich  = hasRich;
        fMDst.richMPr  = richMPr;
        fMDst.richMPi  = richMPi;
    }


/*
    VertexMCInfo & vtx = fG4mc->primVtx;
    if (vtx.status && std::max(std::fabs(vtx.coo[0]), std::fabs(vtx.coo[1])) < 40) {
        std::cout << Form("MGRIG %8.2f VTX %8.2f %8.2f %8.2f\n", MCNRig, vtx.coo[0], vtx.coo[1], vtx.coo[2]);
        std::cout << Form("RIG %8.2f %8.2f %8.2f %8.2f\n", trSign[0] * trNRig[0], trSign[1] * trNRig[1], trSign[2] * trNRig[2], trSign[3] * trNRig[3]);
    //    std::cout << Form("NTRK %ld\n", fG4mc->secParts.size());
    //    for (Int_t it = 0; it < fG4mc->secParts.size(); ++it) {
    //        PartMCInfo & part = fG4mc->secParts.at(it);
    //        if (MGNumc::EqualToZero(part.chrg)) continue;
    //        std::cout << Form("IT%02d ID %03d MOM %8.2f CHRG %6.2f\n", it, part.partID, part.mom, part.chrg);
    //    }
    }
*/
        
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

    YiAna * ntuple = new YiAna();
    ntuple->setOutputFile(outputFile.c_str(), path.c_str());
    ntuple->readDataFrom(file_list.c_str(), group_th, group_size);
    ntuple->loopEventChain();
    if (ntuple != 0) delete ntuple;
    ntuple = 0;

    COUT("\n**------------------------**\n");
    COUT("\n**    YiAnaNtuple END     **\n");
    COUT("\n**------------------------**\n");
    return 0;
}
#endif // __YiAnaNtuple_C__
