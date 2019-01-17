#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/19Jan12/src/ClassDef.h"

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

    FLAGS_logtostderr = true;
    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(google::GLOG_FATAL);

    TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");
    
    TrackSys::MagMgnt::Load();
    TrackSys::MatMgnt::Load();

    MGConfig::JobOpt opt(argc, argv);

    TChain * dst = new TChain("data");
    for (auto&& file : opt.flist()) dst->Add(file.c_str());

    LIST * fList = new LIST;
    G4MC * fG4mc = (opt.mode() == MGConfig::JobOpt::MODE::MC ) ? new G4MC : nullptr;
    RTI  * fRti  = (opt.mode() == MGConfig::JobOpt::MODE::ISS) ? new RTI  : nullptr;
    TRG  * fTrg  = new TRG ;
    TOF  * fTof  = new TOF ;
    ACC  * fAcc  = new ACC ;
    TRK  * fTrk  = new TRK ;
    TRD  * fTrd  = new TRD ;
    RICH * fRich = new RICH;
    ECAL * fEcal = new ECAL;
    HYC  * fHyc  = new HYC ;

    dst->SetBranchAddress("list", &fList);
    if (opt.mode() == MGConfig::JobOpt::MODE::MC)
        dst->SetBranchAddress("g4mc", &fG4mc);
    if (opt.mode() == MGConfig::JobOpt::MODE::ISS)
        dst->SetBranchAddress("rti",  &fRti);
    dst->SetBranchAddress("trg",  &fTrg);
    dst->SetBranchAddress("tof",  &fTof);
    dst->SetBranchAddress("acc",  &fAcc);
    dst->SetBranchAddress("trk",  &fTrk);
    dst->SetBranchAddress("trd",  &fTrd);
    dst->SetBranchAddress("rich", &fRich);
    dst->SetBranchAddress("ecal", &fEcal);
    dst->SetBranchAddress("hyc",  &fHyc);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    TFile * ofle = new TFile(Form("%s/fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");

    PartInfo info(PartType::Proton);
    //PartInfo info(PartType::Helium4);
    
    PartInfo::SetDefault(info.type());
    PhyArg::SetOpt(true, true);
   
    std::vector<Double_t> vmom;
    if (info.type() == PartType::Proton) {
        vmom = std::vector<Double_t>( {
               0.50,   0.80,
               1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
               3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
               8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
               19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
               41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
               93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00,
               800.00 } );
    }
    if (info.type() == PartType::Helium4) {
        vmom = std::vector<Double_t>( {
               0.50,   0.80,
               1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
               3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
               8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
               19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
               41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
               93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00,
               800.00 } );
    }
    
    std::vector<Double_t> vrig;
    for (auto&& val : vmom) vrig.push_back(std::fabs(val / info.chrg()));

    Axis AXmom("Momentum [GeV]", vmom);
    Axis AXrig("Rigidity [GV]", vrig);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    
    Axis AXqlt("Quality", 200, -2.0, 10.0);
    Hist* hPRqltx = Hist::New("hPRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltx = Hist::New("hNRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqlty = Hist::New("hPRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqlty = Hist::New("hNRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqltxAll = Hist::New("hPRqltxAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltxAll = Hist::New("hNRqltxAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqltyAll = Hist::New("hPRqltyAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltyAll = Hist::New("hNRqltyAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    
    Hist* hPRqltxC1 = Hist::New("hPRqltxC1", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltxC1 = Hist::New("hNRqltxC1", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqltyC1 = Hist::New("hPRqltyC1", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltyC1 = Hist::New("hNRqltyC1", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqltxC1All = Hist::New("hPRqltxC1All", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltxC1All = Hist::New("hNRqltxC1All", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqltyC1All = Hist::New("hPRqltyC1All", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltyC1All = Hist::New("hNRqltyC1All", HistAxis(AXrig, AXqlt, "Events/Bin"));
    
    Axis AXtrd("TRD Estimator", 200, 0.0, 1.8);
    Hist* hPRTRDllr = Hist::New("hPRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hNRTRDllr = Hist::New("hNRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hPRTRDllrAll = Hist::New("hPRTRDllrAll", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hNRTRDllrAll = Hist::New("hNRTRDllrAll", HistAxis(AXrig, AXtrd, "Events/Bin"));
   
    Hist* hIRflux = Hist::New("hIRflux", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRflux = Hist::New("hPRflux", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRflux = Hist::New("hNRflux", HistAxis(AXrig, "Events/Bin"));
    Hist* hIRfluxAll = Hist::New("hIRfluxAll", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRfluxAll = Hist::New("hPRfluxAll", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRfluxAll = Hist::New("hNRfluxAll", HistAxis(AXrig, "Events/Bin"));

    Hist* hIRfluxCF10 = Hist::New("hIRfluxCF10", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRfluxCF10 = Hist::New("hPRfluxCF10", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRfluxCF10 = Hist::New("hNRfluxCF10", HistAxis(AXrig, "Events/Bin"));
    Hist* hIRfluxCF10All = Hist::New("hIRfluxCF10All", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRfluxCF10All = Hist::New("hPRfluxCF10All", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRfluxCF10All = Hist::New("hNRfluxCF10All", HistAxis(AXrig, "Events/Bin"));
    
    Hist* hIRfluxCF10C1 = Hist::New("hIRfluxCF10C1", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRfluxCF10C1 = Hist::New("hPRfluxCF10C1", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRfluxCF10C1 = Hist::New("hNRfluxCF10C1", HistAxis(AXrig, "Events/Bin"));
    Hist* hIRfluxCF10C1All = Hist::New("hIRfluxCF10C1All", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRfluxCF10C1All = Hist::New("hPRfluxCF10C1All", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRfluxCF10C1All = Hist::New("hNRfluxCF10C1All", HistAxis(AXrig, "Events/Bin"));
    
    Hist* hIRfluxCF10C2 = Hist::New("hIRfluxCF10C2", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRfluxCF10C2 = Hist::New("hPRfluxCF10C2", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRfluxCF10C2 = Hist::New("hNRfluxCF10C2", HistAxis(AXrig, "Events/Bin"));
    Hist* hIRfluxCF10C2All = Hist::New("hIRfluxCF10C2All", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRfluxCF10C2All = Hist::New("hPRfluxCF10C2All", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRfluxCF10C2All = Hist::New("hNRfluxCF10C2All", HistAxis(AXrig, "Events/Bin"));

    Hist* hPRmuqxAll = Hist::New("hPRmuqxAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRmuqxAll = Hist::New("hNRmuqxAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRmuqyAll = Hist::New("hPRmuqyAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRmuqyAll = Hist::New("hNRmuqyAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRmuqbAll = Hist::New("hPRmuqbAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRmuqbAll = Hist::New("hNRmuqbAll", HistAxis(AXrig, AXqlt, "Events/Bin"));
    
    Axis AXmass("Mass", 2000, -5.0, 5.0);
    Hist* hmass = Hist::New("hmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hmassC1 = Hist::New("hmassC1", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hmassC2 = Hist::New("hmassC2", HistAxis(AXrig, AXmass, "Events/Bin"));
    
    Hist* hmassAll = Hist::New("hmassAll", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hmassC1All = Hist::New("hmassC1All", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hmassC2All = Hist::New("hmassC2All", HistAxis(AXrig, AXmass, "Events/Bin"));
    
    Hist* hLmass = Hist::New("hLmass", HistAxis(AXmass, "Events/Bin"));
    Hist* hLmassC1 = Hist::New("hLmassC1", HistAxis(AXmass, "Events/Bin"));
    Hist* hLmassC2 = Hist::New("hLmassC2", HistAxis(AXmass, "Events/Bin"));
    
    Hist* hLmassAll = Hist::New("hLmassAll", HistAxis(AXmass, "Events/Bin"));
    Hist* hLmassC1All = Hist::New("hLmassC1All", HistAxis(AXmass, "Events/Bin"));
    Hist* hLmassC2All = Hist::New("hLmassC2All", HistAxis(AXmass, "Events/Bin"));
    
    Hist* hPRM1llx = Hist::New("hPRM1llx", HistAxis(AXqlt, "Events/Bin"));
    Hist* hNRM1llx = Hist::New("hNRM1llx", HistAxis(AXqlt, "Events/Bin"));
    Hist* hPRM1lly = Hist::New("hPRM1lly", HistAxis(AXqlt, "Events/Bin"));
    Hist* hNRM1lly = Hist::New("hNRM1lly", HistAxis(AXqlt, "Events/Bin"));
    
    Hist* hPRM2llx = Hist::New("hPRM2llx", HistAxis(AXqlt, "Events/Bin"));
    Hist* hNRM2llx = Hist::New("hNRM2llx", HistAxis(AXqlt, "Events/Bin"));
    Hist* hPRM2lly = Hist::New("hPRM2lly", HistAxis(AXqlt, "Events/Bin"));
    Hist* hNRM2lly = Hist::New("hNRM2lly", HistAxis(AXqlt, "Events/Bin"));
    
    Hist* hM1llx = Hist::New("hM1llx", HistAxis(AXmass, AXqlt, "Events/Bin"));
    Hist* hM1lly = Hist::New("hM1lly", HistAxis(AXmass, AXqlt, "Events/Bin"));
    
    Hist* hM2llx = Hist::New("hM2llx", HistAxis(AXmass, AXqlt, "Events/Bin"));
    Hist* hM2lly = Hist::New("hM2lly", HistAxis(AXmass, AXqlt, "Events/Bin"));
   

    // Choutko
    Hist* hCKPRTRDllr = Hist::New("hCKPRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hCKNRTRDllr = Hist::New("hCKNRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    
    Hist* hCKPRqltx = Hist::New("hCKPRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hCKNRqltx = Hist::New("hCKNRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hCKPRqlty = Hist::New("hCKPRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hCKNRqlty = Hist::New("hCKNRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    
    Hist* hCKIRflux = Hist::New("hCKIRflux", HistAxis(AXirig, "Events/Bin"));
    Hist* hCKPRflux = Hist::New("hCKPRflux", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKNRflux = Hist::New("hCKNRflux", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKIRfluxCF10 = Hist::New("hCKIRfluxCF10", HistAxis(AXirig, "Events/Bin"));
    Hist* hCKPRfluxCF10 = Hist::New("hCKPRfluxCF10", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKNRfluxCF10 = Hist::New("hCKNRfluxCF10", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKIRfluxCF10C1 = Hist::New("hCKIRfluxCF10C1", HistAxis(AXirig, "Events/Bin"));
    Hist* hCKPRfluxCF10C1 = Hist::New("hCKPRfluxCF10C1", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKNRfluxCF10C1 = Hist::New("hCKNRfluxCF10C1", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKIRfluxCF10C2 = Hist::New("hCKIRfluxCF10C2", HistAxis(AXirig, "Events/Bin"));
    Hist* hCKPRfluxCF10C2 = Hist::New("hCKPRfluxCF10C2", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKNRfluxCF10C2 = Hist::New("hCKNRfluxCF10C2", HistAxis(AXrig, "Events/Bin"));
    
    Hist* hCKmass = Hist::New("hCKmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hCKmassC1 = Hist::New("hCKmassC1", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hCKmassC2 = Hist::New("hCKmassC2", HistAxis(AXrig, AXmass, "Events/Bin"));
    
    Hist* hCKLmass = Hist::New("hCKLmass", HistAxis(AXmass, "Events/Bin"));
    Hist* hCKLmassC1 = Hist::New("hCKLmassC1", HistAxis(AXmass, "Events/Bin"));
    Hist* hCKLmassC2 = Hist::New("hCKLmassC2", HistAxis(AXmass, "Events/Bin"));
    
    Axis AXq("Q", 800, 0.0, 3.0);
    Hist* hQL1   = Hist::New("hQL1",   HistAxis(AXmass, AXq, "Events/Bin"));
    Hist* hPRQL1 = Hist::New("hPRQL1", HistAxis(AXq, "Events/Bin"));
    Hist* hNRQL1 = Hist::New("hNRQL1", HistAxis(AXq, "Events/Bin"));
    
    Hist* hQL2   = Hist::New("hQL2",   HistAxis(AXmass, AXq, "Events/Bin"));
    Hist* hPRQL2 = Hist::New("hPRQL2", HistAxis(AXq, "Events/Bin"));
    Hist* hNRQL2 = Hist::New("hNRQL2", HistAxis(AXq, "Events/Bin"));
    
    Hist* hQL9   = Hist::New("hQL9",   HistAxis(AXmass, AXq, "Events/Bin"));
    Hist* hPRQL9 = Hist::New("hPRQL9", HistAxis(AXq, "Events/Bin"));
    Hist* hNRQL9 = Hist::New("hNRQL9", HistAxis(AXq, "Events/Bin"));

    
    Long64_t printRate = dst->GetEntries();
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        dst->GetEntry(entry);
        
        // Geometry (TRK)
        if (fTrk->numOfTrack != 1) continue;
     
        // Geometry (TOF)
        if (fTof->numOfBetaH != 1) continue;
        if (!fTof->statusBetaH) continue;
        
        // Geometry (TRD)
        if (!fTrd->statusKCls[0]) continue;
        if (fTrd->LLRnhit[0] < 8) continue;
        
        // Geometry (ACC)
        if (fAcc->clusters.size() != 0) continue;
        
        // Down-going
        if (fTof->betaH < 0.) continue;

        // Charge
        if (std::abs(info.chrg()) == 1) {
            if (fTof->Qall < 0.8 || fTof->Qall > 1.3) continue;
            if (fTrk->QIn < 0.8 || fTrk->QIn > 1.3) continue;
        }
        if (std::abs(info.chrg()) == 2) {
            if (fTof->Qall < 1.7 || fTof->Qall > 2.4) continue;
            if (fTrk->QIn < 1.7 || fTrk->QIn > 2.4) continue;
        }

        // TOF
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        
        // TRD
        //if (std::abs(info.chrg()) == 1 && fTrd->LLRph[0] > 0.3) continue;
        //if (fTrd->LLRep[0] < 0.7) continue;
        
        // Track In
        CKTrackInfo& cktrIn = fTrk->cktr.at(0);
        HCTrInfo&    hctrIn = fHyc->trM1.at(0);
        HCTrInfo&    hctrInAll = fHyc->trM1All.at(0);
        HCTrInfo&    hctrInAll2 = fHyc->trM2All.at(0);
        
        // Track L1
        CKTrackInfo& cktrL1 = fTrk->cktr.at(1);
        HCTrInfo&    hctrL1 = fHyc->trM1.at(1);
        HCTrInfo&    hctrL1All = fHyc->trM1All.at(1);
        HCTrInfo&    hctrL1All2 = fHyc->trM2All.at(1);
        
        // Track L9
        CKTrackInfo& cktrL9 = fTrk->cktr.at(2);
        HCTrInfo&    hctrL9 = fHyc->trM1.at(2);
        HCTrInfo&    hctrL9All = fHyc->trM1All.at(2);
        HCTrInfo&    hctrL9All2 = fHyc->trM2All.at(2);
        
        // Track Fs
        CKTrackInfo& cktrFs = fTrk->cktr.at(3);
        HCTrInfo&    hctrFs = fHyc->trM1.at(3);
        HCTrInfo&    hctrFsAll = fHyc->trM1All.at(3);
        HCTrInfo&    hctrFsAll2 = fHyc->trM2All.at(3);
   
        // Mass
        HCMuInfo& mutr = fHyc->mutr;
        Short_t musign = (mutr.status ? (mutr.rig[0] > 0 ? 1 : -1) : 0);
        Bool_t  mugood = (mutr.status && mutr.muQlt[0] < 2.0 && mutr.muQlt[1] < 2.0 && mutr.muQlt[2] < 2.0 && mutr.qlt[0] < 2.0 && mutr.qlt[1] < 2.0);

        // IGRF RTI
        double maxCF = 0;
        if (opt.mode() == MGConfig::JobOpt::MODE::ISS && hctrIn.status) {
            double maxStormer = (*std::max_element(fRti->cfStormer, fRti->cfStormer+4));
            double maxIGRF    = (*std::max_element(fRti->cfIGRF,    fRti->cfIGRF+4));
            maxCF = std::max(maxStormer, maxIGRF);
            maxCF = maxIGRF;
        }
        
        // CK Track
        Short_t      ckpatt = -1;
        CKTrackInfo* cktr = nullptr; 
        if      (cktrFs.status) { ckpatt = 3; cktr = &cktrFs; } 
        else if (cktrL9.status) { ckpatt = 2; cktr = &cktrL9; }
        else if (cktrL1.status) { ckpatt = 1; cktr = &cktrL1; }
        else if (cktrIn.status) { ckpatt = 0; cktr = &cktrIn; }
        if (cktr == nullptr || !cktr->status) continue;
        
        Short_t  cksign = TrackSys::Numc::Compare(cktr->rig);
        Double_t ckarig = std::fabs(cktr->rig);
        Double_t ckirig = 1.0 / cktr->rig;
        Double_t ckqltx = std::log(cktr->nchi[0]);
        Double_t ckqlty = std::log(cktr->nchi[1]);
        Double_t ckmass = (fTof->betaH < 1.0) ? std::fabs(info.chrg() * cktrIn.rig * std::sqrt(1.0/fTof->betaH/fTof->betaH-1.0)): -1.0;
        
        Bool_t ckllrc = fTrd->LLRep[0] > (0.2 + 0.6 * cktr->bta * cktr->bta);
       
        Bool_t ckislow = (fTof->betaH > 0.5 && fTof->betaH < 0.8);

        if (cksign > 0) hCKPRqltx->fillH2D(ckarig, ckqltx, fList->weight);
        if (cksign < 0) hCKNRqltx->fillH2D(ckarig, ckqltx, fList->weight);
        if (cksign > 0) hCKPRqlty->fillH2D(ckarig, ckqlty, fList->weight);
        if (cksign < 0) hCKNRqlty->fillH2D(ckarig, ckqlty, fList->weight);
        
        if (cksign > 0) hCKPRTRDllr->fillH2D(ckarig, fTrd->LLRep[0], fList->weight);
        if (cksign < 0) hCKNRTRDllr->fillH2D(ckarig, fTrd->LLRep[0], fList->weight);
        
        hCKIRflux->fillH1D(ckirig, fList->weight);
        if (cksign > 0) hCKPRflux->fillH1D(ckarig, fList->weight);
        if (cksign < 0) hCKNRflux->fillH1D(ckarig, fList->weight);
        
        if (ckarig > 1.0 * maxCF) hCKIRfluxCF10->fillH1D(ckirig, fList->weight);
        if (cksign > 0 && ckarig > 1.0 * maxCF) hCKPRfluxCF10->fillH1D(ckarig, fList->weight);
        if (cksign < 0 && ckarig > 1.0 * maxCF) hCKNRfluxCF10->fillH1D(ckarig, fList->weight);
        
        if (ckllrc && ckarig > 1.0 * maxCF) hCKIRfluxCF10C1->fillH1D(ckirig, fList->weight);
        if (ckllrc && cksign > 0 && ckarig > 1.0 * maxCF) hCKPRfluxCF10C1->fillH1D(ckarig, fList->weight);
        if (ckllrc && cksign < 0 && ckarig > 1.0 * maxCF) hCKNRfluxCF10C1->fillH1D(ckarig, fList->weight);
        
        if (ckllrc && ckqltx < 2.0 && ckqlty < 2.0 && ckarig > 1.0 * maxCF) hCKIRfluxCF10C2->fillH1D(ckirig, fList->weight);
        if (ckllrc && ckqltx < 2.0 && ckqlty < 2.0 && cksign > 0 && ckarig > 1.0 * maxCF) hCKPRfluxCF10C2->fillH1D(ckarig, fList->weight);
        if (ckllrc && ckqltx < 2.0 && ckqlty < 2.0 && cksign < 0 && ckarig > 1.0 * maxCF) hCKNRfluxCF10C2->fillH1D(ckarig, fList->weight);
        
        if (ckmass > 0.0) hCKmass->fillH2D(ckarig, cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        if (ckmass > 0.0 && ckqltx < 2.0 && ckqlty < 2.0) hCKmassC1->fillH2D(ckarig, cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        if (ckmass > 0.0 && ckllrc && ckqltx < 2.0 && ckqlty < 2.0) hCKmassC2->fillH2D(ckarig, cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        
        if (ckislow && ckmass > 0.0) hCKLmass->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        if (ckislow && ckmass > 0.0 && ckqltx < 2.0 && ckqlty < 2.0) hCKLmassC1->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        if (ckislow && ckmass > 0.0 && ckllrc && ckqltx < 2.0 && ckqlty < 2.0) hCKLmassC2->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        
        // HC Track
        Short_t   patt = -1;
        HCTrInfo* hctr = nullptr; 
        if      (hctrFs.status) { patt = 3; hctr = &hctrFs; } 
        else if (hctrL9.status) { patt = 2; hctr = &hctrL9; }
        else if (hctrL1.status) { patt = 1; hctr = &hctrL1; }
        else if (hctrIn.status) { patt = 0; hctr = &hctrIn; }
        if (hctr == nullptr || !hctr->status) continue;
        
        Short_t  sign = TrackSys::Numc::Compare(hctr->rig[0]);
        Double_t arig = std::fabs(hctr->rig[0]);
        Double_t irig = (TrackSys::Numc::ONE<> / hctr->rig[0]);
        Double_t qltx = hctr->qlt[0];
        Double_t qlty = hctr->qlt[1];
        Double_t qltb = fHyc->btaM1.qlt;

        Bool_t llrc = fTrd->LLRep[0] > (0.2 + 0.6 * hctr->bta[0] * hctr->bta[0]);

        Bool_t islow = (fHyc->btaM1.status && fHyc->btaM1.bta[0] > 0.5 && fHyc->btaM1.bta[0] < 0.8);

        if (sign > 0) hPRqltx->fillH2D(arig, qltx, fList->weight);
        if (sign < 0) hNRqltx->fillH2D(arig, qltx, fList->weight);
        if (sign > 0) hPRqlty->fillH2D(arig, qlty, fList->weight);
        if (sign < 0) hNRqlty->fillH2D(arig, qlty, fList->weight);
        
        if (llrc && sign > 0) hPRqltxC1->fillH2D(arig, qltx, fList->weight);
        if (llrc && sign < 0) hNRqltxC1->fillH2D(arig, qltx, fList->weight);
        if (llrc && sign > 0) hPRqltyC1->fillH2D(arig, qlty, fList->weight);
        if (llrc && sign < 0) hNRqltyC1->fillH2D(arig, qlty, fList->weight);
        
        if (sign > 0) hPRTRDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
        if (sign < 0) hNRTRDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
        
        hIRflux->fillH1D(irig, fList->weight);
        if (sign > 0) hPRflux->fillH1D(arig, fList->weight);
        if (sign < 0) hNRflux->fillH1D(arig, fList->weight);
        
        if (arig > 1.0 * maxCF) hIRfluxCF10->fillH1D(irig, fList->weight);
        if (sign > 0 && arig > 1.0 * maxCF) hPRfluxCF10->fillH1D(arig, fList->weight);
        if (sign < 0 && arig > 1.0 * maxCF) hNRfluxCF10->fillH1D(arig, fList->weight);
        
        if (llrc && arig > 1.0 * maxCF) hIRfluxCF10C1->fillH1D(irig, fList->weight);
        if (llrc && sign > 0 && arig > 1.0 * maxCF) hPRfluxCF10C1->fillH1D(arig, fList->weight);
        if (llrc && sign < 0 && arig > 1.0 * maxCF) hNRfluxCF10C1->fillH1D(arig, fList->weight);
        
        if (llrc && qltx < 2.0 && qlty < 2.0 && arig > 1.0 * maxCF) hIRfluxCF10C2->fillH1D(irig, fList->weight);
        if (llrc && qltx < 2.0 && qlty < 2.0 && sign > 0 && arig > 1.0 * maxCF) hPRfluxCF10C2->fillH1D(arig, fList->weight);
        if (llrc && qltx < 2.0 && qlty < 2.0 && sign < 0 && arig > 1.0 * maxCF) hNRfluxCF10C2->fillH1D(arig, fList->weight);
        
        if (fHyc->massM1 > 0) hmass->fillH2D(arig, sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
        if (fHyc->massM1 > 0 && qltx < 2.0 && qlty < 2.0 && qltb < 2.0) hmassC1->fillH2D(arig, sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
        if (fHyc->massM1 > 0 && llrc && qltx < 2.0 && qlty < 2.0 && qltb < 2.0) hmassC2->fillH2D(arig, sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
        
        if (islow && fHyc->massM1 > 0) hLmass->fillH1D(sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
        if (islow && fHyc->massM1 > 0 && qltx < 2.0 && qlty < 2.0 && qltb < 2.0) hLmassC1->fillH1D(sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
        if (islow && fHyc->massM1 > 0 && llrc && qltx < 2.0 && qlty < 2.0 && qltb < 2.0) hLmassC2->fillH1D(sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
       

        Short_t   pattAll = -1;
        HCTrInfo* hctrAll = nullptr; 
        if      (hctrFsAll.status) { pattAll = 3; hctrAll = &hctrFsAll; } 
        else if (hctrL9All.status) { pattAll = 2; hctrAll = &hctrL9All; }
        else if (hctrL1All.status) { pattAll = 1; hctrAll = &hctrL1All; }
        else if (hctrInAll.status) { pattAll = 0; hctrAll = &hctrInAll; }
        if (hctrAll == nullptr || !hctrAll->status) continue;

        Short_t  signAll = TrackSys::Numc::Compare(hctrAll->rig[0]);
        Double_t arigAll = std::fabs(hctrAll->rig[0]);
        Double_t irigAll = (TrackSys::Numc::ONE<> / hctrAll->rig[0]);
        Double_t qltxAll = hctrAll->qlt[0];
        Double_t qltyAll = hctrAll->qlt[1];
        
        Bool_t llrcAll = fTrd->LLRep[0] > (0.2 + 0.6 * mutr.bta[0] * mutr.bta[0]);

        Bool_t islowAll = (mutr.status && mutr.bta[0] > 0.5 && mutr.bta[0] < 0.8 && (mutr.rig[0]*hctr->rig[0] > 0));

        if (signAll > 0) hPRqltxAll->fillH2D(arigAll, qltxAll, fList->weight);
        if (signAll < 0) hNRqltxAll->fillH2D(arigAll, qltxAll, fList->weight);
        if (signAll > 0) hPRqltyAll->fillH2D(arigAll, qltyAll, fList->weight);
        if (signAll < 0) hNRqltyAll->fillH2D(arigAll, qltyAll, fList->weight);
        
        if (llrcAll && signAll > 0) hPRqltxC1All->fillH2D(arigAll, qltxAll, fList->weight);
        if (llrcAll && signAll < 0) hNRqltxC1All->fillH2D(arigAll, qltxAll, fList->weight);
        if (llrcAll && signAll > 0) hPRqltyC1All->fillH2D(arigAll, qltyAll, fList->weight);
        if (llrcAll && signAll < 0) hNRqltyC1All->fillH2D(arigAll, qltyAll, fList->weight);
        
        if (signAll > 0) hPRTRDllrAll->fillH2D(arigAll, fTrd->LLRep[0], fList->weight);
        if (signAll < 0) hNRTRDllrAll->fillH2D(arigAll, fTrd->LLRep[0], fList->weight);

        hIRfluxAll->fillH1D(irigAll, fList->weight);
        if (signAll > 0) hPRfluxAll->fillH1D(arigAll, fList->weight);
        if (signAll < 0) hNRfluxAll->fillH1D(arigAll, fList->weight);

        if (arigAll > 1.0 * maxCF) hIRfluxCF10All->fillH1D(irigAll, fList->weight);
        if (signAll > 0 && arigAll > 1.0 * maxCF) hPRfluxCF10All->fillH1D(arigAll, fList->weight);
        if (signAll < 0 && arigAll > 1.0 * maxCF) hNRfluxCF10All->fillH1D(arigAll, fList->weight);
        
        if (llrcAll && arigAll > 1.0 * maxCF) hIRfluxCF10C1All->fillH1D(irigAll, fList->weight);
        if (llrcAll && signAll > 0 && arigAll > 1.0 * maxCF) hPRfluxCF10C1All->fillH1D(arigAll, fList->weight);
        if (llrcAll && signAll < 0 && arigAll > 1.0 * maxCF) hNRfluxCF10C1All->fillH1D(arigAll, fList->weight);
        
        if (llrcAll && qltxAll < 2.0 && qltyAll < 2.0 && arigAll > 1.0 * maxCF) hIRfluxCF10C2All->fillH1D(irigAll, fList->weight);
        if (llrcAll && qltxAll < 2.0 && qltyAll < 2.0 && signAll > 0 && arigAll > 1.0 * maxCF) hPRfluxCF10C2All->fillH1D(arigAll, fList->weight);
        if (llrcAll && qltxAll < 2.0 && qltyAll < 2.0 && signAll < 0 && arigAll > 1.0 * maxCF) hNRfluxCF10C2All->fillH1D(arigAll, fList->weight);
        
        if (mutr.status && signAll > 0) hPRmuqxAll->fillH2D(arigAll, mutr.muQlt[0], fList->weight);
        if (mutr.status && signAll < 0) hNRmuqxAll->fillH2D(arigAll, mutr.muQlt[0], fList->weight);
        if (mutr.status && signAll > 0) hPRmuqyAll->fillH2D(arigAll, mutr.muQlt[1], fList->weight);
        if (mutr.status && signAll < 0) hNRmuqyAll->fillH2D(arigAll, mutr.muQlt[1], fList->weight);
        if (mutr.status && signAll > 0) hPRmuqbAll->fillH2D(arigAll, mutr.muQlt[2], fList->weight);
        if (mutr.status && signAll < 0) hNRmuqbAll->fillH2D(arigAll, mutr.muQlt[2], fList->weight);

        if (mutr.status) hmassAll->fillH2D(arigAll, signAll * mutr.mass, fList->weight * (signAll>0?0.001:1.0));
        if (mutr.status && mugood) hmassC1All->fillH2D(arigAll, signAll * mutr.mass, fList->weight * (signAll>0?0.001:1.0));
        if (mutr.status && mugood && llrcAll) hmassC2All->fillH2D(arigAll, signAll * mutr.mass, fList->weight * (signAll>0?0.001:1.0));
        
        if (islowAll && mutr.status) hLmassAll->fillH1D(signAll * mutr.mass, fList->weight * (signAll>0?0.001:1.0));
        if (islowAll && mutr.status && mugood) hLmassC1All->fillH1D(signAll * mutr.mass, fList->weight * (signAll>0?0.001:1.0));
        if (islowAll && mutr.status && mugood && llrcAll) hLmassC2All->fillH1D(signAll * mutr.mass, fList->weight * (signAll>0?0.001:1.0));
        
        Short_t   pattAll2 = -1;
        HCTrInfo* hctrAll2 = nullptr; 
        if      (hctrFsAll2.status) { pattAll2 = 3; hctrAll2 = &hctrFsAll2; } 
        else if (hctrL9All2.status) { pattAll2 = 2; hctrAll2 = &hctrL9All2; }
        else if (hctrL1All2.status) { pattAll2 = 1; hctrAll2 = &hctrL1All2; }
        else if (hctrInAll2.status) { pattAll2 = 0; hctrAll2 = &hctrInAll2; }
        if (hctrAll2 == nullptr || !hctrAll2->status) continue;
        
        if (islowAll && mutr.status && mugood && mutr.rig[0] > 0) hPRM1llx->fillH1D(hctrAll->qlt[0], fList->weight);
        if (islowAll && mutr.status && mugood && mutr.rig[0] < 0) hNRM1llx->fillH1D(hctrAll->qlt[0], fList->weight);
        if (islowAll && mutr.status && mugood && mutr.rig[0] > 0) hPRM1lly->fillH1D(hctrAll->qlt[1], fList->weight);
        if (islowAll && mutr.status && mugood && mutr.rig[0] < 0) hNRM1lly->fillH1D(hctrAll->qlt[1], fList->weight);
        
        if (islowAll && mutr.status && mugood && mutr.rig[0] > 0) hPRM2llx->fillH1D(hctrAll2->qlt[0], fList->weight);
        if (islowAll && mutr.status && mugood && mutr.rig[0] < 0) hNRM2llx->fillH1D(hctrAll2->qlt[0], fList->weight);
        if (islowAll && mutr.status && mugood && mutr.rig[0] > 0) hPRM2lly->fillH1D(hctrAll2->qlt[1], fList->weight);
        if (islowAll && mutr.status && mugood && mutr.rig[0] < 0) hNRM2lly->fillH1D(hctrAll2->qlt[1], fList->weight);
        
        if (islowAll && mutr.status && mugood) hM1llx->fillH2D(musign * mutr.mass, hctrAll->qlt[0], fList->weight * (musign>0?0.001:1.0));
        if (islowAll && mutr.status && mugood) hM1lly->fillH2D(musign * mutr.mass, hctrAll->qlt[1], fList->weight * (musign>0?0.001:1.0));
        
        if (islowAll && mutr.status && mugood) hM2llx->fillH2D(musign * mutr.mass, hctrAll2->qlt[0], fList->weight * (musign>0?0.001:1.0));
        if (islowAll && mutr.status && mugood) hM2lly->fillH2D(musign * mutr.mass, hctrAll2->qlt[1], fList->weight * (musign>0?0.001:1.0));
        
        if (islowAll && mutr.status && mugood && cktrL1.status && fTrk->QL1 > 0) hQL1->fillH2D(musign * mutr.mass, fTrk->QL1, fList->weight * (musign>0?0.001:1.0));
        if (islowAll && mutr.status && mugood && cktrL1.status && fTrk->QL1 > 0 && mutr.rig[0] > 0) hPRQL1->fillH1D(fTrk->QL1, fList->weight);
        if (islowAll && mutr.status && mugood && cktrL1.status && fTrk->QL1 > 0 && mutr.rig[0] < 0) hNRQL1->fillH1D(fTrk->QL1, fList->weight);
        
        if (islowAll && mutr.status && mugood && fTrk->QL2 > 0) hQL2->fillH2D(musign * mutr.mass, fTrk->QL2, fList->weight * (musign>0?0.001:1.0));
        if (islowAll && mutr.status && mugood && fTrk->QL2 > 0 && mutr.rig[0] > 0) hPRQL2->fillH1D(fTrk->QL2, fList->weight);
        if (islowAll && mutr.status && mugood && fTrk->QL2 > 0 && mutr.rig[0] < 0) hNRQL2->fillH1D(fTrk->QL2, fList->weight);
        
        if (islowAll && mutr.status && mugood && cktrL9.status && fTrk->QL9 > 0) hQL9->fillH2D(musign * mutr.mass, fTrk->QL9, fList->weight * (musign>0?0.001:1.0));
        if (islowAll && mutr.status && mugood && cktrL9.status && fTrk->QL9 > 0 && mutr.rig[0] > 0) hPRQL9->fillH1D(fTrk->QL9, fList->weight);
        if (islowAll && mutr.status && mugood && cktrL9.status && fTrk->QL9 > 0 && mutr.rig[0] < 0) hNRQL9->fillH1D(fTrk->QL9, fList->weight);
    }
    ofle->Write();
    ofle->Close();

    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    if (fList) { delete fList; fList = nullptr; }
    if (fG4mc) { delete fG4mc; fG4mc = nullptr; }
    if (fRti ) { delete fRti ; fRti  = nullptr; }
    if (fTrg ) { delete fTrg ; fTrg  = nullptr; }
    if (fTof ) { delete fTof ; fTof  = nullptr; }
    if (fAcc ) { delete fAcc ; fAcc  = nullptr; }
    if (fTrk ) { delete fTrk ; fTrk  = nullptr; }
    if (fTrd ) { delete fTrd ; fTrd  = nullptr; }
    if (fRich) { delete fRich; fRich = nullptr; }
    if (fEcal) { delete fEcal; fEcal = nullptr; }
    if (fHyc ) { delete fHyc;  fHyc  = nullptr; }

    google::ShutdownGoogleLogging();
    return 0;
}
