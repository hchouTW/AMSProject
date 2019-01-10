#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/19Jan08/src/ClassDef.h"

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
    
    Int_t    nmom     = 20;
    Double_t mombd[2] = { 1., 1000. };
    if (info.type() == PartType::Proton)   { mombd[0] = 0.55; mombd[1] = 30.0; }
    if (info.type() == PartType::Helium4)  { mombd[0] = 2.20; mombd[1] = 60.0; }
    Axis AXmom("Momentum [GeV]", nmom, mombd[0], mombd[1], AxisScale::kLog);
    Axis AXrig("Rigidity [GV]", AXmom.nbin(), mombd[0]/std::fabs(info.chrg()), mombd[1]/std::fabs(info.chrg()), AxisScale::kLog);
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
        if (fTof->betaHPatt != 15) continue;
        
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
        
        if (fTof->numOfInTimeCls > 4) continue;
        if ((fTof->numOfExtCls[0]+fTof->numOfExtCls[1]) > 0 || 
            (fTof->numOfExtCls[2]+fTof->numOfExtCls[3]) > 1) continue; 
    
        // TRD
        if (std::abs(info.chrg()) == 1 && fTrd->LLRph[0] > 0.3) continue;
        //if (fTrd->LLRep[0] < 0.7) continue;
        
        // Track In
        CKTrackInfo& cktrIn = fTrk->cktr.at(0);
        HCTrInfo&    hctrIn = fHyc->trM1.at(0);
        
        // Track L1
        CKTrackInfo& cktrL1 = fTrk->cktr.at(1);
        HCTrInfo&    hctrL1 = fHyc->trM1.at(1);
        
        // Track L9
        CKTrackInfo& cktrL9 = fTrk->cktr.at(2);
        HCTrInfo&    hctrL9 = fHyc->trM1.at(2);
        
        // Track Fs
        CKTrackInfo& cktrFs = fTrk->cktr.at(3);
        HCTrInfo&    hctrFs = fHyc->trM1.at(3);
   
        // Mass
        HCMuInfo& mutr = fHyc->mutr;
        bool goodmu = (mutr.status && mutr.bta[0] < 0.99 && mutr.muQuality[0] < 2.5 && mutr.muQuality[1] < 2.5 && mutr.muQuality[2] < 2.5);

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
        Double_t qltx = hctr->quality[0];
        Double_t qlty = hctr->quality[1];
       
        Double_t ckmass = (fTof->betaH < 1.0) ? std::fabs(info.chrg() * cktrIn.rig * std::sqrt(1.0/fTof->betaH/fTof->betaH-1.0)): -1.0;

        if (sign > 0) hPRqltx->fillH2D(arig, qltx, fList->weight);
        if (sign < 0) hNRqltx->fillH2D(arig, qltx, fList->weight);
        if (sign > 0) hPRqlty->fillH2D(arig, qlty, fList->weight);
        if (sign < 0) hNRqlty->fillH2D(arig, qlty, fList->weight);

        hIRflux->fillH1D(irig, fList->weight);
        if (sign > 0) hPRflux->fillH1D(arig, fList->weight);
        if (sign < 0) hNRflux->fillH1D(arig, fList->weight);
        
        if (sign > 0) hPRTRDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
        if (sign < 0) hNRTRDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
        
        if (sign > 0 && goodmu) hPRmass->fillH2D(arig, mutr.mass, fList->weight);
        if (sign < 0 && goodmu) hNRmass->fillH2D(arig, mutr.mass, fList->weight);
        
        if (cktrIn.rig > 0 && ckmass > 0) hPRmass2->fillH2D(std::fabs(cktrIn.rig), ckmass, fList->weight);
        if (cktrIn.rig < 0 && ckmass > 0) hNRmass2->fillH2D(std::fabs(cktrIn.rig), ckmass, fList->weight);
        
        if (sign > 0 && fHyc->massM1 > 0 && fHyc->btaM1.bta[0] < 0.99) hPRmass3->fillH2D(arig, fHyc->massM1, fList->weight);
        if (sign < 0 && fHyc->massM1 > 0 && fHyc->btaM1.bta[0] < 0.99) hNRmass3->fillH2D(arig, fHyc->massM1, fList->weight);
        if (sign > 0 && fHyc->massM2 > 0 && fHyc->btaM2.bta[0] < 0.99) hPRmass4->fillH2D(arig, fHyc->massM2, fList->weight);
        if (sign < 0 && fHyc->massM2 > 0 && fHyc->btaM2.bta[0] < 0.99) hNRmass4->fillH2D(arig, fHyc->massM2, fList->weight);
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
