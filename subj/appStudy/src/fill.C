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
               0.20,   0.50,   0.80,
               1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
               3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
               8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
               19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
               41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
               93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } );
    }
    if (info.type() == PartType::Helium4) {
        vmom = std::vector<Double_t>( {
               0.20,   0.50,   0.80,
               1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
               3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
               8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
               19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
               41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
               93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } );
    }
    
    std::vector<Double_t> vrig;
    for (auto&& val : vmom) vrig.push_back(std::fabs(val / info.chrg()));

    Axis AXmom("Momentum [GeV]", vmom);
    Axis AXrig("Rigidity [GV]", vrig);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    
    Axis AXtrd("TRD Estimator", 200, 0.0, 1.8);
    Hist* hPRTRDllr = Hist::New("hPRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hNRTRDllr = Hist::New("hNRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    
    Axis AXqlt("Quality", 400, -2.0, 10.0);
    Hist* hPRqltx = Hist::New("hPRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltx = Hist::New("hNRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqlty = Hist::New("hPRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqlty = Hist::New("hNRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    
    Axis AXmass("Mass", 400, 0.03, 6.0);
    Hist* hCKPRmass = Hist::New("hCKPRmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hCKNRmass = Hist::New("hCKNRmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    
    Hist* hPRmass = Hist::New("hPRmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hNRmass = Hist::New("hNRmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    
    Hist* hPRmass2 = Hist::New("hPRmass2", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hNRmass2 = Hist::New("hNRmass2", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hPRmass3 = Hist::New("hPRmass3", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hNRmass3 = Hist::New("hNRmass3", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hPRmass4 = Hist::New("hPRmass4", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hNRmass4 = Hist::New("hNRmass4", HistAxis(AXrig, AXmass, "Events/Bin"));
    
    Hist* hIRflux = Hist::New("hIRflux", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRflux = Hist::New("hPRflux", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRflux = Hist::New("hNRflux", HistAxis(AXrig, "Events/Bin"));
    
    Hist* hIRflux2 = Hist::New("hIRflux2", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRflux2 = Hist::New("hPRflux2", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRflux2 = Hist::New("hNRflux2", HistAxis(AXrig, "Events/Bin"));
    
    Hist* hIRflux3 = Hist::New("hIRflux3", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRflux3 = Hist::New("hPRflux3", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRflux3 = Hist::New("hNRflux3", HistAxis(AXrig, "Events/Bin"));
    
    Hist* hIRflux4 = Hist::New("hIRflux4", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRflux4 = Hist::New("hPRflux4", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRflux4 = Hist::New("hNRflux4", HistAxis(AXrig, "Events/Bin"));
    
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
        //if (fTrd->numOfTrack != 1 && fTrd->numOfHTrack != 1) continue;
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
        if (std::abs(info.chrg()) == 1 && fTrd->LLRph[0] > 0.3) continue;
        //if (fTrd->LLRep[0] < 0.7) continue;
        
        // Track In
        CKTrackInfo& cktrIn = fTrk->cktr.at(0);
        HCTrInfo&    hctrIn = fHyc->trM1.at(0);
        HCTrInfo&    hctrInAll = fHyc->trM1All.at(0);
        
        // Track L1
        CKTrackInfo& cktrL1 = fTrk->cktr.at(1);
        HCTrInfo&    hctrL1 = fHyc->trM1.at(1);
        HCTrInfo&    hctrL1All = fHyc->trM1All.at(1);
        
        // Track L9
        CKTrackInfo& cktrL9 = fTrk->cktr.at(2);
        HCTrInfo&    hctrL9 = fHyc->trM1.at(2);
        HCTrInfo&    hctrL9All = fHyc->trM1All.at(2);
        
        // Track Fs
        CKTrackInfo& cktrFs = fTrk->cktr.at(3);
        HCTrInfo&    hctrFs = fHyc->trM1.at(3);
        HCTrInfo&    hctrFsAll = fHyc->trM1All.at(3);
   
        // Mass
        HCMuInfo& mutr = fHyc->mutr;
        bool goodmu = (mutr.status && mutr.bta[0] < 1.0 && mutr.muQlt[0] < 2.5 && mutr.muQlt[1] < 2.5 && mutr.muQlt[2] < 2.5);

        // IGRF RTI
        double maxCF = 0;
        if (opt.mode() == MGConfig::JobOpt::MODE::ISS && hctrIn.status) {
            double maxStormer = (*std::max_element(fRti->cfStormer, fRti->cfStormer+4));
            double maxIGRF    = (*std::max_element(fRti->cfIGRF,    fRti->cfIGRF+4));
            maxCF = std::max(maxStormer, maxIGRF);
            maxCF = maxIGRF;

            //CERR("%14.8f STORMER %14.8f IGRF %14.8f\n", maxCF, maxStormer, maxIGRF);
            //if (std::fabs(hctrIn.rig[0]) < 1.2 * maxCF) continue;
        }
        
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
        
        HCTrInfo* hctrAll = nullptr; 
        if      (hctrFsAll.status) { patt = 3; hctrAll = &hctrFsAll; } 
        else if (hctrL9All.status) { patt = 2; hctrAll = &hctrL9All; }
        else if (hctrL1All.status) { patt = 1; hctrAll = &hctrL1All; }
        else if (hctrInAll.status) { patt = 0; hctrAll = &hctrInAll; }
        if (hctrAll == nullptr || !hctrAll->status) continue;

        Short_t  signAll = TrackSys::Numc::Compare(hctrAll->rig[0]);
        Double_t arigAll = std::fabs(hctrAll->rig[0]);
        Double_t irigAll = (TrackSys::Numc::ONE<> / hctrAll->rig[0]);
        Double_t qltxAll = hctrAll->qlt[0];
        Double_t qltyAll = hctrAll->qlt[1];
        
        Double_t ckmass = (fTof->betaH < 1.0) ? std::fabs(info.chrg() * cktrIn.rig * std::sqrt(1.0/fTof->betaH/fTof->betaH-1.0)): -1.0;
        
        if (cktrIn.rig > 0 && ckmass > 0 && cktrIn.nchi[0] < 8 && cktrIn.nchi[1] < 8) hCKPRmass->fillH2D(std::fabs(cktrIn.rig), ckmass, fList->weight);
        if (cktrIn.rig < 0 && ckmass > 0 && cktrIn.nchi[0] < 8 && cktrIn.nchi[1] < 8) hCKNRmass->fillH2D(std::fabs(cktrIn.rig), ckmass, fList->weight);

        if (sign > 0) hPRqltx->fillH2D(arig, qltx, fList->weight);
        if (sign < 0) hNRqltx->fillH2D(arig, qltx, fList->weight);
        if (sign > 0) hPRqlty->fillH2D(arig, qlty, fList->weight);
        if (sign < 0) hNRqlty->fillH2D(arig, qlty, fList->weight);
        
        if (sign > 0) hPRTRDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
        if (sign < 0) hNRTRDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);

        if (goodmu) hIRflux->fillH1D(irig, fList->weight);
        if (sign > 0 && goodmu) hPRflux->fillH1D(arig, fList->weight);
        if (sign < 0 && goodmu) hNRflux->fillH1D(arig, fList->weight);
        
        if (sign > 0 && goodmu) hPRmass->fillH2D(arig, mutr.mass, fList->weight);
        if (sign < 0 && goodmu) hNRmass->fillH2D(arig, mutr.mass, fList->weight);
        
        if (goodmu && qltx < 2.0 && qlty < 2.0) hIRflux2->fillH1D(irig, fList->weight);
        if (sign > 0 && goodmu && qltx < 2.0 && qlty < 2.0) hPRflux2->fillH1D(arig, fList->weight);
        if (sign < 0 && goodmu && qltx < 2.0 && qlty < 2.0) hNRflux2->fillH1D(arig, fList->weight);
        
        if (sign > 0 && goodmu && qltx < 2.0 && qlty < 2.0) hPRmass2->fillH2D(arig, mutr.mass, fList->weight);
        if (sign < 0 && goodmu && qltx < 2.0 && qlty < 2.0) hNRmass2->fillH2D(arig, mutr.mass, fList->weight);
       
        if (goodmu) hIRflux3->fillH1D(irigAll, fList->weight);
        if (signAll > 0 && goodmu) hPRflux3->fillH1D(arigAll, fList->weight);
        if (signAll < 0 && goodmu) hNRflux3->fillH1D(arigAll, fList->weight);
        
        if (signAll > 0 && goodmu) hPRmass3->fillH2D(arigAll, mutr.mass, fList->weight);
        if (signAll < 0 && goodmu) hNRmass3->fillH2D(arigAll, mutr.mass, fList->weight);
        
        if (goodmu && qltxAll < 2.0 && qltyAll < 2.0) hIRflux4->fillH1D(irigAll, fList->weight);
        if (signAll > 0 && goodmu && qltxAll < 2.0 && qltyAll < 2.0) hPRflux4->fillH1D(arigAll, fList->weight);
        if (signAll < 0 && goodmu && qltxAll < 2.0 && qltyAll < 2.0) hNRflux4->fillH1D(arigAll, fList->weight);
        
        if (signAll > 0 && goodmu && qltxAll < 2.0 && qltyAll < 2.0) hPRmass4->fillH2D(arigAll, mutr.mass, fList->weight);
        if (signAll < 0 && goodmu && qltxAll < 2.0 && qltyAll < 2.0) hNRmass4->fillH2D(arigAll, mutr.mass, fList->weight);
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
