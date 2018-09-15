#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

//#include "/ams_home/hchou/AMSCore/prod/18Jul04/src/ClassDef.h"
#include "/ams_home/hchou/AMSCore/prod/18Sep15/src/ClassDef.h"
//#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/18Jul04/src/ClassDef.h"

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(google::GLOG_FATAL);

    TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");
    
    //TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    
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
    //dst->SetBranchAddress("ecal", &fEcal);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    TFile * ofle = new TFile(Form("%s/fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXrig("Rigidity [GV]", 80, 0.55, 2000., AxisScale::kLog);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    
    // Fit Eff
    Hist* hMCnum = Hist::New("hMCnum", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKnum = Hist::New("hCKnum", HistAxis(AXrig, "Events/Bin"));
    Hist* hKFnum = Hist::New("hKFnum", HistAxis(AXrig, "Events/Bin"));
    Hist* hHCnum = Hist::New("hHCnum", HistAxis(AXrig, "Events/Bin"));

    Hist* hMCtme = Hist::New("hMCtme", HistAxis(AXrig, ""));
    Hist* hCKtme = Hist::New("hCKtme", HistAxis(AXrig, "Mean Time"));
    Hist* hKFtme = Hist::New("hKFtme", HistAxis(AXrig, "Mean Time"));
    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXrig, "Mean Time"));
    
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 1600, -1.0, 1.0);
    Hist* hCKRrso = Hist::New("hCKRrso", HistAxis(AXrig, AXRrso));
    Hist* hKFRrso = Hist::New("hKFRrso", HistAxis(AXrig, AXRrso));
    Hist* hHCRrso = Hist::New("hHCRrso", HistAxis(AXrig, AXRrso));
    
    Axis AXqlt("Quality [1]", 400, -2.5, 8.0);
    Hist* hCKqltx = Hist::New("hCKqltx", HistAxis(AXrig, AXqlt));
    Hist* hKFqltx = Hist::New("hKFqltx", HistAxis(AXrig, AXqlt));
    Hist* hHCqltx = Hist::New("hHCqltx", HistAxis(AXrig, AXqlt));
    
    Hist* hCKqlty = Hist::New("hCKqlty", HistAxis(AXrig, AXqlt));
    Hist* hKFqlty = Hist::New("hKFqlty", HistAxis(AXrig, AXqlt));
    Hist* hHCqlty = Hist::New("hHCqlty", HistAxis(AXrig, AXqlt));
    
    Axis AXc("Cos", 1000, -0.5, 0.5);
    Hist* hCKcx = Hist::New("hCKcx", HistAxis(AXrig, AXc));
    Hist* hKFcx = Hist::New("hKFcx", HistAxis(AXrig, AXc));
    Hist* hHCcx = Hist::New("hHCcx", HistAxis(AXrig, AXc));
    Hist* hCKcy = Hist::New("hCKcy", HistAxis(AXrig, AXc));
    Hist* hKFcy = Hist::New("hKFcy", HistAxis(AXrig, AXc));
    Hist* hHCcy = Hist::New("hHCcy", HistAxis(AXrig, AXc));
    
    Axis AXu("Cos", 400, -0.05, 0.05);
    Hist* hCKux = Hist::New("hCKux", HistAxis(AXrig, AXu));
    Hist* hKFux = Hist::New("hKFux", HistAxis(AXrig, AXu));
    Hist* hHCux = Hist::New("hHCux", HistAxis(AXrig, AXu));
    Hist* hCKuy = Hist::New("hCKuy", HistAxis(AXrig, AXu));
    Hist* hKFuy = Hist::New("hKFuy", HistAxis(AXrig, AXu));
    Hist* hHCuy = Hist::New("hHCuy", HistAxis(AXrig, AXu));
    
    MGClock::HrsStopwatch hrssw; hrssw.start();
    Long64_t printRate = static_cast<Long64_t>(0.1 * dst->GetEntries());
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) {
            hrssw.stop();
            COUT("Entry %lld/%lld Time %14.8f\n", entry, dst->GetEntries(), hrssw.time());
        }
        dst->GetEntry(entry);

        Int_t patt = 0;

        CKTrackInfo& ckTr = fTrk->ckTr.at(patt);
        KFTrackInfo& kfTr = fTrk->kfTr.at(patt);
        HCTrackInfo& hcTr = fTrk->hcTr.at(patt); // Tracker
    
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

        Bool_t status = (ckTr.status && kfTr.status && hcTr.status);
        if (!status) continue;
        
        Double_t bta  = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? std::fabs(hcTr.stateTop[7]) : fG4mc->primPart.bta);
        Double_t mom  = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? hcTr.stateTop[6] : fG4mc->primPart.mom);
        Double_t imom = Numc::ONE<> / mom;
        Double_t scl  = std::sqrt(AXrig.center(AXrig.find(mom), AxisScale::kLog));
        
        Short_t ckSign = (ckTr.rig > 0) ? 1 : -1;
        Short_t kfSign = (kfTr.rig[0] > 0) ? 1 : -1;
        Short_t hcSign = (hcTr.stateTop[6] > 0) ? 1 : -1;
        
        Double_t ckRig = ckTr.rig;
        Double_t kfRig = kfTr.rig[0];
        Double_t hcRig = hcTr.stateTop[6];
        
        Double_t ckIRig = Numc::ONE<> / ckRig;
        Double_t kfIRig = Numc::ONE<> / kfRig;
        Double_t hcIRig = Numc::ONE<> / hcRig;
        
        Double_t ck_qltx = std::log(ckTr.nchi[0]); 
        Double_t kf_qltx = std::log(kfTr.nchi[0]); 
        Double_t hc_qltx = hcTr.quality[0]; 
        
        Double_t ck_qlty = std::log(ckTr.nchi[1]); 
        Double_t kf_qlty = std::log(kfTr.nchi[1]); 
        Double_t hc_qlty = hcTr.quality[1]; 
       
        hMCnum->fillH1D(mom, wgt);
        if (ckSign > 0) hCKnum->fillH1D(ckRig, wgt);
        if (kfSign > 0) hKFnum->fillH1D(kfRig, wgt);
        if (hcSign > 0) hHCnum->fillH1D(hcRig, wgt);
        
        hMCtme->fillH1D(mom);
        hCKtme->fillH1D(mom, ckTr.cpuTime);
        hKFtme->fillH1D(mom, kfTr.cpuTime);
        hHCtme->fillH1D(mom, hcTr.cpuTime);
        
        hCKRrso->fillH2D(mom, scl * (ckIRig - imom));
        hKFRrso->fillH2D(mom, scl * (kfIRig - imom));
        hHCRrso->fillH2D(mom, scl * (hcIRig - imom));
        
        hCKqltx->fillH2D(mom, ck_qltx);
        hKFqltx->fillH2D(mom, kf_qltx);
        hHCqltx->fillH2D(mom, hc_qltx);
        
        hCKqlty->fillH2D(mom, ck_qlty);
        hKFqlty->fillH2D(mom, kf_qlty);
        hHCqlty->fillH2D(mom, hc_qlty);
        
        hCKcx->fillH2D(mom, scl * (ckTr.stateTop[0] - fG4mc->primPart.coo[0]));
        hKFcx->fillH2D(mom, scl * (kfTr.stateTop[0] - fG4mc->primPart.coo[0]));
        hHCcx->fillH2D(mom, scl * (hcTr.stateTop[0] - fG4mc->primPart.coo[0]));
        
        hCKcy->fillH2D(mom, scl * (ckTr.stateTop[1] - fG4mc->primPart.coo[1]));
        hKFcy->fillH2D(mom, scl * (kfTr.stateTop[1] - fG4mc->primPart.coo[1]));
        hHCcy->fillH2D(mom, scl * (hcTr.stateTop[1] - fG4mc->primPart.coo[1]));
        
        hCKux->fillH2D(mom, scl * (ckTr.stateTop[3] - fG4mc->primPart.dir[0]));
        hKFux->fillH2D(mom, scl * (kfTr.stateTop[3] - fG4mc->primPart.dir[0]));
        hHCux->fillH2D(mom, scl * (hcTr.stateTop[3] - fG4mc->primPart.dir[0]));
        
        hCKuy->fillH2D(mom, scl * (ckTr.stateTop[4] - fG4mc->primPart.dir[1]));
        hKFuy->fillH2D(mom, scl * (kfTr.stateTop[4] - fG4mc->primPart.dir[1]));
        hHCuy->fillH2D(mom, scl * (hcTr.stateTop[4] - fG4mc->primPart.dir[1]));
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
