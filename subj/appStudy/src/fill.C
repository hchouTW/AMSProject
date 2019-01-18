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
    
    Axis AXtrd("TRD Estimator", 200, 0.0, 1.8);
    Hist* hCKPRTRDllr = Hist::New("hCKPRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hCKNRTRDllr = Hist::New("hCKNRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hPRTRDllr = Hist::New("hPRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hNRTRDllr = Hist::New("hNRTRDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    
    Axis AXqlt("Quality", 200, -2.0, 10.0);
    Hist* hCKPRqltx = Hist::New("hCKPRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hCKNRqltx = Hist::New("hCKNRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hCKPRqlty = Hist::New("hCKPRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hCKNRqlty = Hist::New("hCKNRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqltx = Hist::New("hPRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqltx = Hist::New("hNRqltx", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hPRqlty = Hist::New("hPRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hNRqlty = Hist::New("hNRqlty", HistAxis(AXrig, AXqlt, "Events/Bin"));
    
    Hist* hCKIRflux = Hist::New("hCKIRflux", HistAxis(AXirig, "Events/Bin"));
    Hist* hCKPRflux = Hist::New("hCKPRflux", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKNRflux = Hist::New("hCKNRflux", HistAxis(AXrig, "Events/Bin"));
    Hist* hIRflux = Hist::New("hIRflux", HistAxis(AXirig, "Events/Bin"));
    Hist* hPRflux = Hist::New("hPRflux", HistAxis(AXrig, "Events/Bin"));
    Hist* hNRflux = Hist::New("hNRflux", HistAxis(AXrig, "Events/Bin"));

    Axis AXmass("Mass", 500, -5.0, 5.0);
    Hist* hCKmass   = Hist::New("hCKmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hCKmassC  = Hist::New("hCKmassC", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hCKLmass  = Hist::New("hCKLmass", HistAxis(AXmass, "Events/Bin"));
    Hist* hCKLmassC = Hist::New("hCKLmassC", HistAxis(AXmass, "Events/Bin"));
   
    Hist* hmass   = Hist::New("hmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hmassC  = Hist::New("hmassC", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hLmass  = Hist::New("hLmass", HistAxis(AXmass, "Events/Bin"));
    Hist* hLmassC = Hist::New("hLmassC", HistAxis(AXmass, "Events/Bin"));
    
    Hist* hHCmass   = Hist::New("hHCmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hHCmassC  = Hist::New("hHCmassC", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hHCLmass  = Hist::New("hHCLmass", HistAxis(AXmass, "Events/Bin"));
    Hist* hHCLmassC = Hist::New("hHCLmassC", HistAxis(AXmass, "Events/Bin"));

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
            
            if (fTrk->QL2 > 0 && fTrk->QL2 < 0.7) continue;
            if (fTrk->QL1 > 0 && fTrk->QL1 < 0.7) continue;
            if (fTrk->QL9 > 0 && fTrk->QL9 < 0.7) continue;
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
        HCTrInfo&    hctrM1In = fHyc->trM1All.at(0);
        HCTrInfo&    hctrM2In = fHyc->trM2All.at(0);
        
        // Track L1
        CKTrackInfo& cktrL1 = fTrk->cktr.at(1);
        HCTrInfo&    hctrL1 = fHyc->trM1.at(1);
        HCTrInfo&    hctrM1L1 = fHyc->trM1All.at(1);
        HCTrInfo&    hctrM2L1 = fHyc->trM2All.at(1);
        
        // Track L9
        CKTrackInfo& cktrL9 = fTrk->cktr.at(2);
        HCTrInfo&    hctrL9 = fHyc->trM1.at(2);
        HCTrInfo&    hctrM1L9 = fHyc->trM1All.at(2);
        HCTrInfo&    hctrM2L9 = fHyc->trM2All.at(2);
        
        // Track Fs
        CKTrackInfo& cktrFs = fTrk->cktr.at(3);
        HCTrInfo&    hctrFs = fHyc->trM1.at(3);
        HCTrInfo&    hctrM1Fs = fHyc->trM1All.at(3);
        HCTrInfo&    hctrM2Fs = fHyc->trM2All.at(3);
   
        // Mass
        HCMuInfo& mutr = fHyc->mutr;
        Short_t musign = (mutr.status ? (mutr.rig[0] > 0 ? 1 : -1) : 0);
        Float_t muarig = (mutr.status ? std::fabs(mutr.rig[0]) : 0);
        Float_t muirig = (mutr.status ? 1.0/mutr.rig[0] : 0);
        Bool_t  mugood = (mutr.status && mutr.muQlt[0] < 2.0 && mutr.muQlt[1] < 2.0 && mutr.muQlt[2] < 2.0 && mutr.qlt[0] < 2.0 && mutr.qlt[1] < 2.0);
        Bool_t  mulow  = (mutr.status && mutr.bta[0] > 0.5 && mutr.bta[0] < 0.8);
        Bool_t  mullrc = fTrd->LLRep[0] > (0.2 + 0.6 * mutr.bta[0] * mutr.bta[0]);

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
        Bool_t   ckllrc = fTrd->LLRep[0] > (0.2 + 0.6 * cktr->bta * cktr->bta);
        Bool_t   cklow  = (fTof->betaH > 0.5 && fTof->betaH < 0.8);
        Bool_t   ckcf   = (ckarig > 1.0 * maxCF);

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
        Bool_t   llrc = fTrd->LLRep[0] > (0.2 + 0.6 * hctr->bta[0] * hctr->bta[0]);
        Bool_t   low  = (fHyc->btaM1.status && fHyc->btaM1.bta[0] > 0.5 && fHyc->btaM1.bta[0] < 0.8);
        Bool_t   cf   = (arig > 1.0 * maxCF);

        Short_t   pattM1 = -1;
        HCTrInfo* hctrM1 = nullptr; 
        if      (hctrM1Fs.status) { pattM1 = 3; hctrM1 = &hctrM1Fs; } 
        else if (hctrM1L9.status) { pattM1 = 2; hctrM1 = &hctrM1L9; }
        else if (hctrM1L1.status) { pattM1 = 1; hctrM1 = &hctrM1L1; }
        else if (hctrM1In.status) { pattM1 = 0; hctrM1 = &hctrM1In; }
        if (hctrM1 == nullptr || !hctrM1->status) continue;

        Short_t  signM1 = TrackSys::Numc::Compare(hctrM1->rig[0]);
        Double_t arigM1 = std::fabs(hctrM1->rig[0]);
        Double_t irigM1 = (TrackSys::Numc::ONE<> / hctrM1->rig[0]);
        Double_t qltxM1 = hctrM1->qlt[0];
        Double_t qltyM1 = hctrM1->qlt[1];
        Bool_t   cfM1   = (arigM1 > 1.0 * maxCF);

        Short_t   pattM2 = -1;
        HCTrInfo* hctrM2 = nullptr; 
        if      (hctrM2Fs.status) { pattM2 = 3; hctrM2 = &hctrM2Fs; } 
        else if (hctrM2L9.status) { pattM2 = 2; hctrM2 = &hctrM2L9; }
        else if (hctrM2L1.status) { pattM2 = 1; hctrM2 = &hctrM2L1; }
        else if (hctrM2In.status) { pattM2 = 0; hctrM2 = &hctrM2In; }
        if (hctrM2 == nullptr || !hctrM2->status) continue;

        Short_t  signM2 = TrackSys::Numc::Compare(hctrM2->rig[0]);
        Double_t arigM2 = std::fabs(hctrM2->rig[0]);
        Double_t irigM2 = (TrackSys::Numc::ONE<> / hctrM2->rig[0]);
        Double_t qltxM2 = hctrM2->qlt[0];
        Double_t qltyM2 = hctrM2->qlt[1];
        Bool_t   cfM2   = (arigM2 > 1.0 * maxCF);

        // Choutko
        if (ckcf && cksign > 0) hCKPRTRDllr->fillH2D(ckarig, fTrd->LLRep[0], fList->weight);
        if (ckcf && cksign < 0) hCKNRTRDllr->fillH2D(ckarig, fTrd->LLRep[0], fList->weight);

        if (ckllrc && ckcf && cksign > 0) hCKPRqltx->fillH2D(ckarig, ckqltx, fList->weight);
        if (ckllrc && ckcf && cksign < 0) hCKNRqltx->fillH2D(ckarig, ckqltx, fList->weight);
        if (ckllrc && ckcf && cksign > 0) hCKPRqlty->fillH2D(ckarig, ckqlty, fList->weight);
        if (ckllrc && ckcf && cksign < 0) hCKNRqlty->fillH2D(ckarig, ckqlty, fList->weight);

        if (ckllrc && ckcf && ckqltx < 2.0 && ckqlty < 2.0) hCKIRflux->fillH1D(ckirig, fList->weight);
        if (ckllrc && ckcf && ckqltx < 2.0 && ckqlty < 2.0 && cksign > 0) hCKPRflux->fillH1D(ckarig, fList->weight);
        if (ckllrc && ckcf && ckqltx < 2.0 && ckqlty < 2.0 && cksign < 0) hCKNRflux->fillH1D(ckarig, fList->weight);

        if (ckmass > 0 && ckllrc) hCKmass->fillH2D(ckarig, cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        if (ckmass > 0 && ckllrc && ckqltx < 2.0 && ckqlty < 2.0) hCKmassC->fillH2D(ckarig, cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        
        if (cklow && ckmass > 0 && ckllrc) hCKLmass->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
        if (cklow && ckmass > 0 && ckllrc && ckqltx < 2.0 && ckqlty < 2.0) hCKLmassC->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.001:1.0));
       
        // HYChou
        if (cf && sign > 0) hPRTRDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
        if (cf && sign < 0) hNRTRDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
        
        if (llrc && cf && sign > 0) hPRqltx->fillH2D(arig, qltx, fList->weight);
        if (llrc && cf && sign < 0) hNRqltx->fillH2D(arig, qltx, fList->weight);
        if (llrc && cf && sign > 0) hPRqlty->fillH2D(arig, qlty, fList->weight);
        if (llrc && cf && sign < 0) hNRqlty->fillH2D(arig, qlty, fList->weight);
        
        if (llrc && cf && qltx < 2.0 && qlty < 2.0) hIRflux->fillH1D(irig, fList->weight);
        if (llrc && cf && qltx < 2.0 && qlty < 2.0 && sign > 0) hPRflux->fillH1D(arig, fList->weight);
        if (llrc && cf && qltx < 2.0 && qlty < 2.0 && sign < 0) hNRflux->fillH1D(arig, fList->weight);
        
        if (fHyc->massM1 > 0 && llrc) hmass->fillH2D(arig, sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
        if (fHyc->massM1 > 0 && llrc && qltx < 2.0 && qlty < 2.0 && qltb < 2.0) hmassC->fillH2D(arig, sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
        
        if (low && fHyc->massM1 > 0 && llrc) hLmass->fillH1D(sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
        if (low && fHyc->massM1 > 0 && llrc && qltx < 2.0 && qlty < 2.0 && qltb < 2.0) hLmassC->fillH1D(sign * fHyc->massM1, fList->weight * (sign>0?0.001:1.0));
       
        if (mutr.status && mullrc) hHCmass->fillH2D(muarig, musign * mutr.mass, fList->weight * (musign>0?0.001:1.0));
        if (mutr.status && mullrc && mugood) hHCmassC->fillH2D(muarig, musign * mutr.mass, fList->weight * (musign>0?0.001:1.0));

        if (low && mutr.status && mullrc) hHCLmass->fillH1D( musign * mutr.mass, fList->weight * (musign>0?0.001:1.0));
        if (low && mutr.status && mullrc && mugood) hHCLmassC->fillH1D(musign * mutr.mass, fList->weight * (musign>0?0.001:1.0));
        
        if (mulow && mutr.status && mugood) hM1llx->fillH2D(musign * mutr.mass, qltxM1, fList->weight * (musign>0?0.001:1.0));
        if (mulow && mutr.status && mugood) hM1lly->fillH2D(musign * mutr.mass, qltyM1, fList->weight * (musign>0?0.001:1.0));
        
        if (mulow && mutr.status && mugood) hM2llx->fillH2D(musign * mutr.mass, qltxM2, fList->weight * (musign>0?0.001:1.0));
        if (mulow && mutr.status && mugood) hM2lly->fillH2D(musign * mutr.mass, qltyM2, fList->weight * (musign>0?0.001:1.0));
        
        if (mulow && mutr.status && mugood && mutr.rig[0] > 0) hPRM1llx->fillH1D(qltxM1, fList->weight);
        if (mulow && mutr.status && mugood && mutr.rig[0] < 0) hNRM1llx->fillH1D(qltxM1, fList->weight);
        if (mulow && mutr.status && mugood && mutr.rig[0] > 0) hPRM1lly->fillH1D(qltyM1, fList->weight);
        if (mulow && mutr.status && mugood && mutr.rig[0] < 0) hNRM1lly->fillH1D(qltyM1, fList->weight);
        
        if (mulow && mutr.status && mugood && mutr.rig[0] > 0) hPRM2llx->fillH1D(qltxM2, fList->weight);
        if (mulow && mutr.status && mugood && mutr.rig[0] < 0) hNRM2llx->fillH1D(qltxM2, fList->weight);
        if (mulow && mutr.status && mugood && mutr.rig[0] > 0) hPRM2lly->fillH1D(qltyM2, fList->weight);
        if (mulow && mutr.status && mugood && mutr.rig[0] < 0) hNRM2lly->fillH1D(qltyM2, fList->weight);
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
