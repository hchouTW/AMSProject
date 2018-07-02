#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

//#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/18May15/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18May19/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18May27/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18Jun10/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18Jun18/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18Jun23/src/ClassDef.h"
#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/18Jun28/src/ClassDef.h"

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(google::GLOG_FATAL);

    //TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");
    
    TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    
    //TrackSys::Sys::ShowMsg( TrackSys::Sys::GetEnv("TRACKSys_MagBox") );
    //TrackSys::Sys::ShowMsg( TrackSys::Sys::GetEnv("TRACKSys_MatBox") );

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
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    TFile * ofle = new TFile(Form("%s/fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXrig("Rigidity [GV]", 100, 0.55, 1000., AxisScale::kLog);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    
    // Fit Eff
    Hist* hMCnum   = Hist::New("hMCnum",   HistAxis(AXrig, "Events/Bin"));
    Hist* hCKnum   = Hist::New("hCKnum",   HistAxis(AXrig, "Events/Bin"));
    Hist* hKFnum   = Hist::New("hKFnum",   HistAxis(AXrig, "Events/Bin"));
    Hist* hHCnum   = Hist::New("hHCnum",   HistAxis(AXrig, "Events/Bin"));
    Hist* hHCnumMU = Hist::New("hHCnumMU", HistAxis(AXrig, "Events/Bin"));

    Hist* hMCtme   = Hist::New("hMCtme",   HistAxis(AXrig, ""));
    Hist* hCKtme   = Hist::New("hCKtme",   HistAxis(AXrig, "Mean Time"));
    Hist* hKFtme   = Hist::New("hKFtme",   HistAxis(AXrig, "Mean Time"));
    Hist* hHCtme   = Hist::New("hHCtme",   HistAxis(AXrig, "Mean Time"));
    Hist* hHCtmeMU = Hist::New("hHCtmeMU", HistAxis(AXrig, "Mean Time"));
    
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 1600, -1.0, 1.0);
    Hist* hCKRrso   = Hist::New("hCKRrso",   HistAxis(AXrig, AXRrso));
    Hist* hKFRrso   = Hist::New("hKFRrso",   HistAxis(AXrig, AXRrso));
    Hist* hHCRrso   = Hist::New("hHCRrso",   HistAxis(AXrig, AXRrso));
    Hist* hHCRrsoMU = Hist::New("hHCRrsoMU", HistAxis(AXrig, AXRrso));
    
    Axis AXBrso("(Bm - Bt) [1]", 1600, -0.1, 0.1);
    Hist* hCKBrso   = Hist::New("hCKBrso",   HistAxis(AXrig, AXBrso));
    Hist* hKFBrso   = Hist::New("hKFBrso",   HistAxis(AXrig, AXBrso));
    Hist* hHCBrso   = Hist::New("hHCBrso",   HistAxis(AXrig, AXBrso));
    Hist* hHCBrsoMU = Hist::New("hHCBrsoMU", HistAxis(AXrig, AXBrso));
    
    Hist* hTFBrso   = Hist::New("hTFBrso",   HistAxis(AXrig, AXBrso));
    Hist* hAGLBrso  = Hist::New("hAGLBrso",  HistAxis(AXrig, AXBrso));
    Hist* hNAFBrso  = Hist::New("hNAFBrso",  HistAxis(AXrig, AXBrso));
    
    Axis AXMrso("Mass [GeV]", 400, 0.03, 5.0);
    Hist* hCKMrso   = Hist::New("hCKMrso",   HistAxis(AXrig, AXMrso));
    Hist* hKFMrso   = Hist::New("hKFMrso",   HistAxis(AXrig, AXMrso));
    Hist* hHCMrso   = Hist::New("hHCMrso",   HistAxis(AXrig, AXMrso));
    Hist* hHCMrsoMU = Hist::New("hHCMrsoMU", HistAxis(AXrig, AXMrso));
       
    Hist* hHCMrsoMUCO = Hist::New("hHCMrsoMUCO", HistAxis(AXrig));
    Hist* hHCMrsoMUEO = Hist::New("hHCMrsoMUEO", HistAxis(AXrig));
    Hist* hHCMrsoMUSO = Hist::New("hHCMrsoMUSO", HistAxis(AXrig, AXMrso));
   
    Hist* hHCMrsoMUC[3] = { nullptr };
    Hist* hHCMrsoMUE[3] = { nullptr };
    Hist* hHCMrsoMUS[3] = { nullptr };
    std::vector<Double_t> cuts({ 0.5, 0.7, 1.0 });
    for (Int_t is = 0; is < 3; ++is) {
        hHCMrsoMUC[is] = Hist::New(Form("hHCMrsoMUC%d", is), HistAxis(AXrig));
        hHCMrsoMUE[is] = Hist::New(Form("hHCMrsoMUE%d", is), HistAxis(AXrig));
        hHCMrsoMUS[is] = Hist::New(Form("hHCMrsoMUS%d", is), HistAxis(AXrig, AXMrso));
    }
    
    Axis AXchi("Log-Chi-square [1]", 400, -3.0, 8.0);
    Hist* hCKchix   = Hist::New("hCKchix",   HistAxis(AXrig, AXchi));
    Hist* hKFchix   = Hist::New("hKFchix",   HistAxis(AXrig, AXchi));
    Hist* hHCchix   = Hist::New("hHCchix",   HistAxis(AXrig, AXchi));
    Hist* hHCchixMU = Hist::New("hHCchixMU", HistAxis(AXrig, AXchi));
    
    Hist* hCKchiy   = Hist::New("hCKchiy",   HistAxis(AXrig, AXchi));
    Hist* hKFchiy   = Hist::New("hKFchiy",   HistAxis(AXrig, AXchi));
    Hist* hHCchiy   = Hist::New("hHCchiy",   HistAxis(AXrig, AXchi));
    Hist* hHCchiyMU = Hist::New("hHCchiyMU", HistAxis(AXrig, AXchi));
    
    //Axis AXu("Cos", 400, -0.05, 0.05);
    //Hist* hCKux = Hist::New("hCKux", HistAxis(AXrig, AXu));
    //Hist* hKFux = Hist::New("hKFux", HistAxis(AXrig, AXu));
    //Hist* hHCux = Hist::New("hHCux", HistAxis(AXrig, AXu));
    //Hist* hCKuy = Hist::New("hCKuy", HistAxis(AXrig, AXu));
    //Hist* hKFuy = Hist::New("hKFuy", HistAxis(AXrig, AXu));
    //Hist* hHCuy = Hist::New("hHCuy", HistAxis(AXrig, AXu));

    MGClock::HrsStopwatch hrssw; hrssw.start();
    Long64_t printRate = static_cast<Long64_t>(0.04 * dst->GetEntries());
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) {
            hrssw.stop();
            COUT("Entry %lld/%lld Time %14.8f\n", entry, dst->GetEntries(), hrssw.time());
        }
        dst->GetEntry(entry);

        CKTrackInfo& ckTr = fTrk->ckTr.at(0);
        KFTrackInfo& kfTr = fTrk->kfTr.at(0);
        //HCTrackInfo& hcTr = fTrk->hcPrTr.at(0); // Tracker
        //HCTrackInfo& hcMu = fTrk->hcMuInTr.at(0); // Tracker + TOF
        HCTrackInfo& hcTr = fTrk->hcPrInTr.at(0); // Tracker + TOF
        HCTrackInfo& hcMu = fTrk->hcMuInTr.at(0); // Tracker + TOF
        //HCTrackInfo& hcTr = fTrk->hcPrInTr.at(1); // Tracker + TOF + RICH
        //HCTrackInfo& hcMu = fTrk->hcMuInTr.at(1); // Tracker + TOF + RICH
    
        // Reweight (MC)
        Double_t wgt = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? 1.0 : std::pow(fG4mc->primPart.mom/AXrig.min(), -1.7));

        // Geometry (TOF)
        if (fTof->numOfBetaH != 1) continue;
        if (!fTof->statusBetaH) continue;
        if (fTof->betaHPatt != 15) continue;
       
        // Geometry (TRD)
        if (fTrd->numOfTrack != 1 && fTrd->numOfHTrack != 1) continue;
        if (!fTrd->statusKCls[0]) continue;
        if (fTrd->LLRnhit[0] < 8) continue;
        
        // Geometry (ACC)
        if (fAcc->clusters.size() != 0) continue;
        
        // Down-going
        if (fTof->betaH < 0.) continue;

        // Charge
        if (fTof->Qall < 0.8 || fTof->Qall > 1.3) continue;
        if (fTrk->QIn < 0.8 || fTrk->QIn > 1.3) continue;

        // TOF
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        
        if (fTof->numOfInTimeCls > 4) continue;
        if ((fTof->numOfExtCls[0]+fTof->numOfExtCls[1]) > 0 || 
            (fTof->numOfExtCls[2]+fTof->numOfExtCls[3]) > 1) continue; 

        // TRD
        if (fTrd->LLRep[0] < 0.7) continue;

        // RICH
        //if (!fRich->status || !fRich->isGood) continue;
        //if (fRich->kind != 0) continue; // AGL

        Bool_t status = (ckTr.status && kfTr.status && hcTr.status && hcMu.status);
        if (!status) continue;
        if (hcTr.cpuTime > 1000.) continue; // rmove material loading event
        
        Short_t ckSign   = (ckTr.rig > 0) ? 1 : -1;
        Short_t kfSign   = (kfTr.rig[0] > 0) ? 1 : -1;
        Short_t hcSign   = (hcTr.stateTop[6] > 0) ? 1 : -1;
        Short_t hcSignMU = (hcMu.stateTop[6] > 0) ? 1 : -1;
        
        Double_t ckRig   = ckTr.rig;
        Double_t kfRig   = kfTr.rig[0];
        Double_t hcRig   = hcTr.stateTop[6];
        Double_t hcRigMU = hcMu.stateTop[6];
        
        Double_t ckIRig   = Numc::ONE<> / ckRig;
        Double_t kfIRig   = Numc::ONE<> / kfRig;
        Double_t hcIRig   = Numc::ONE<> / hcRig;
        Double_t hcIRigMU = Numc::ONE<> / hcRigMU;
        
        Double_t ckBta   = ckTr.bta;
        Double_t kfBta   = kfTr.bta[0];
        Double_t hcBta   = hcTr.stateTop[7];
        Double_t hcBtaMU = hcMu.stateTop[7];
        
        Double_t ck_chix   = std::log(ckTr.nchi[0]); 
        Double_t kf_chix   = std::log(kfTr.nchi[0]); 
        Double_t hc_chix   = hcTr.quality[0]; 
        Double_t hc_chixMU = hcMu.quality[0]; 
        
        Double_t ck_chiy   = std::log(ckTr.nchi[1]); 
        Double_t kf_chiy   = std::log(kfTr.nchi[1]); 
        Double_t hc_chiy   = hcTr.quality[1]; 
        Double_t hc_chiyMU = hcMu.quality[1]; 
        
        Double_t ckMass   = ((fTof->betaH >= 1.0) ? -1.0 : std::fabs(ckTr.rig           * std::sqrt(1.0/fTof->betaH/fTof->betaH - 1.0)));
        Double_t kfMass   = ((fTof->betaH >= 1.0) ? -1.0 : std::fabs(kfTr.rig[1]        * std::sqrt(1.0/fTof->betaH/fTof->betaH - 1.0)));
        Double_t hcMass   = ((fTof->betaH >= 1.0) ? -1.0 : std::fabs(hcTr.stateLJ[3][6] * std::sqrt(1.0/fTof->betaH/fTof->betaH - 1.0)));
        Double_t hcMassMU = (hcMu.mass); 
        
        //Double_t ckMass   = ((fRich->beta >= 1.0) ? -1.0 : std::fabs(ckTr.rig           * std::sqrt(1.0/fRich->beta/fRich->beta - 1.0)));
        //Double_t kfMass   = ((fRich->beta >= 1.0) ? -1.0 : std::fabs(kfTr.rig[2]        * std::sqrt(1.0/fRich->beta/fRich->beta - 1.0)));
        //Double_t hcMass   = ((fRich->beta >= 1.0) ? -1.0 : std::fabs(hcTr.stateLJ[3][6] * std::sqrt(1.0/fRich->beta/fRich->beta - 1.0)));
        //Double_t hcMassMU = (hcMu.mass); 

        Double_t mom  = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? hcTr.stateTop[6] : fG4mc->primPart.mom);
        Double_t imom = Numc::ONE<> / mom;
        Double_t cen  = std::sqrt(AXrig.center(AXrig.find(mom), AxisScale::kLog));
        
        Double_t bta = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? std::fabs(hcTr.stateTop[7]) : fG4mc->primPart.bta);

        hMCnum->fillH1D(mom, wgt);
        if (ckSign   > 0) hCKnum  ->fillH1D(ckRig,   wgt);
        if (kfSign   > 0) hKFnum  ->fillH1D(kfRig,   wgt);
        if (hcSign   > 0) hHCnum  ->fillH1D(hcRig,   wgt);
        if (hcSignMU > 0) hHCnumMU->fillH1D(hcRigMU, wgt);
        
        hMCtme->fillH1D(mom);
        hCKtme  ->fillH1D(mom, ckTr.cpuTime);
        hKFtme  ->fillH1D(mom, kfTr.cpuTime);
        hHCtme  ->fillH1D(mom, hcTr.cpuTime);
        hHCtmeMU->fillH1D(mom, hcMu.cpuTime);
        
        hCKRrso  ->fillH2D(mom, cen * (ckIRig   - imom));
        hKFRrso  ->fillH2D(mom, cen * (kfIRig   - imom));
        hHCRrso  ->fillH2D(mom, cen * (hcIRig   - imom));
        hHCRrsoMU->fillH2D(mom, cen * (hcIRigMU - imom));
        
        hCKBrso  ->fillH2D(mom, (ckBta   - bta));
        hKFBrso  ->fillH2D(mom, (kfBta   - bta));
        hHCBrso  ->fillH2D(mom, (hcBta   - bta));
        hHCBrsoMU->fillH2D(mom, (hcBtaMU - bta));
        
        hTFBrso->fillH2D(mom, (fTof->betaH - bta));
        if (fRich->kind == 0 && fRich->status) hAGLBrso->fillH2D(mom, (fRich->beta - bta));
        if (fRich->kind == 1 && fRich->status) hNAFBrso->fillH2D(mom, (fRich->beta - bta));
        
        hCKMrso  ->fillH2D(mom, ckMass);
        hKFMrso  ->fillH2D(mom, kfMass);
        hHCMrso  ->fillH2D(mom, hcMass);
        hHCMrsoMU->fillH2D(mom, hcMassMU);
  
        if (hcMu.mass > 1.5) hHCMrsoMUCO->fillH1D(mom);
        hHCMrsoMUEO->fillH1D(mom);
        hHCMrsoMUSO->fillH2D(mom, hcMassMU);
        for (Int_t is = 0; is < 3; ++is) {
            if (hcMu.quality[0] > 2.0) continue;
            if (hcMu.quality[1] > cuts.at(is)) continue;
            if (hcMu.mass > 1.5) hHCMrsoMUC[is]->fillH1D(mom);
            hHCMrsoMUE[is]->fillH1D(mom);
            hHCMrsoMUS[is]->fillH2D(mom, hcMassMU);
        }
        
        hCKchix  ->fillH2D(mom, ck_chix);
        hKFchix  ->fillH2D(mom, kf_chix);
        hHCchix  ->fillH2D(mom, hc_chix);
        hHCchixMU->fillH2D(mom, hc_chixMU);
        
        hCKchiy  ->fillH2D(mom, ck_chiy);
        hKFchiy  ->fillH2D(mom, kf_chiy);
        hHCchiy  ->fillH2D(mom, hc_chiy);
        hHCchiyMU->fillH2D(mom, hc_chiyMU);

        //hCKux->fillH2D(mom, cen * (ckTr.stateTop[3] - fG4mc->primPart.dir[0]));
        //hKFux->fillH2D(mom, cen * (kfTr.stateTop[3] - fG4mc->primPart.dir[0]));
        //hHCux->fillH2D(mom, cen * (hcTr.stateTop[3] - fG4mc->primPart.dir[0]));
        //
        //hCKuy->fillH2D(mom, cen * (ckTr.stateTop[4] - fG4mc->primPart.dir[1]));
        //hKFuy->fillH2D(mom, cen * (kfTr.stateTop[4] - fG4mc->primPart.dir[1]));
        //hHCuy->fillH2D(mom, cen * (hcTr.stateTop[4] - fG4mc->primPart.dir[1]));
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

    google::ShutdownGoogleLogging();
    return 0;
}
