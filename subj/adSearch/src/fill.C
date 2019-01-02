#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/18Dec23/src/ClassDef.h"
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
    HYC  * fHyc  = new HYC;

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
    
    Int_t    nmom     = 100;
    Double_t mombd[2] = { 1., 1000. };
    if (info.type() == PartType::Proton)   { mombd[0] = 0.55; mombd[1] = 50.0; }
    if (info.type() == PartType::Helium4)  { mombd[0] = 2.20; mombd[1] = 100.0; }
    Axis AXmom("Momentum [GeV]", nmom, mombd[0], mombd[1], AxisScale::kLog);
    Axis AXrig("Rigidity [GV]", AXmom.nbin(), mombd[0]/std::fabs(info.chrg()), mombd[1]/std::fabs(info.chrg()), AxisScale::kLog);
    
    // Num (InnTr)
    Hist* hMCnum = Hist::New("hMCnum", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKnum = Hist::New("hCKnum", HistAxis(AXrig, "Events/Bin"));
    Hist* hHCnum = Hist::New("hHCnum", HistAxis(AXrig, "Events/Bin"));
   
    // RigDef
    Axis AXRDrso("(1/R_L9 - 1/R_L1) [1/GV]", 1600, -1.0, 1.0);
    Hist* hCKRDrso = Hist::New("hCKRDrso", HistAxis(AXrig, AXRDrso));
    Hist* hHCRDrso = Hist::New("hHCRDrso", HistAxis(AXrig, AXRDrso));
    
    // Time
    Axis AXtme("Time [ms]", 1600, 0., 1000.);
    Hist* hJFevt = Hist::New("hJFevt", HistAxis(AXrig));
    Hist* hHCevt = Hist::New("hHCevt", HistAxis(AXrig));
    Hist* hJFtme = Hist::New("hJFtme", HistAxis(AXrig, "Tme"));
    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXrig, "Tme"));
    
    // Fit M Res
    Axis AXM("Mass", 1000, 0.03, 8.0);
    Hist* hJFM = Hist::New("hJFM", HistAxis(AXrig, AXM));
    Hist* hHCM = Hist::New("hHCM", HistAxis(AXrig, AXM));
    
    Hist* hHCM2 = Hist::New("hHCM2", HistAxis(AXrig, AXM));
    
    //Axis AXMqlt("Quality [1]", 800, -2.0, 4.0);
    //Hist* hHCMqltx = Hist::New("hHCMqltx", HistAxis(AXmom, AXMqlt));
    //Hist* hHCMqlty = Hist::New("hHCMqlty", HistAxis(AXmom, AXMqlt));
    //Hist* hHCMqltb = Hist::New("hHCMqltb", HistAxis(AXmom, AXMqlt));
    
    MGClock::HrsStopwatch hrssw; hrssw.start();
    Long64_t printRate = static_cast<Long64_t>(0.1 * dst->GetEntries());
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) {
            hrssw.stop();
            COUT("Entry %lld/%lld Time %14.8f\n", entry, dst->GetEntries(), hrssw.time());
        }
        dst->GetEntry(entry);
    
        // Reweight (MC)
        Double_t wgt = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? 1.0 : std::pow(fG4mc->primPart.mom/AXmom.min(), -1.7));
        
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
        if (fTrd->LLRep[0] < 0.7) continue;
        
        // Track In
        CKTrackInfo& ckTrIn = fTrk->ckTr.at(0);
        HCTrInfo&    hcTrIn = fHyc->trM1.at(0);
        
        // Track L1
        CKTrackInfo& ckTrL1 = fTrk->ckTr.at(1);
        HCTrInfo&    hcTrL1 = fHyc->trM1.at(1);
        
        // Track L9
        CKTrackInfo& ckTrL9 = fTrk->ckTr.at(2);
        HCTrInfo&    hcTrL9 = fHyc->trM1.at(2);
        
        // Track Fs
        CKTrackInfo& ckTrFs = fTrk->ckTr.at(3);
        HCTrInfo&    hcTrFs = fHyc->trM1.at(3);
        
        if (opt.mode() == MGConfig::JobOpt::MODE::MC) hMCnum->fillH1D(fG4mc->primPart.mom/info.chrg(), wgt);
        if (ckTrIn.status) hCKnum->fillH1D(ckTrIn.rig,    wgt);
        if (hcTrIn.status) hHCnum->fillH1D(hcTrIn.rig[0], wgt);
        
        if (ckTrL1.status & ckTrL9.status && ckTrFs.status) {
            double rigx = std::fabs(ckTrFs.rig);
            double ckrd = (1.0/ckTrL9.rig - 1.0/ckTrL1.rig);
            double sclx = std::sqrt(AXrig.center(AXrig.find(rigx), AxisScale::kLog));
            hCKRDrso->fillH2D(ckTrFs.rig, sclx * ckrd);
        }
        
        if (hcTrL1.status & hcTrL9.status && hcTrFs.status) {
            double rigx = std::fabs(hcTrFs.rig[0]);
            double hcrd = (1.0/hcTrL9.rig[0] - 1.0/hcTrL1.rig[0]);
            double sclx = std::sqrt(AXrig.center(AXrig.find(rigx), AxisScale::kLog));
            hHCRDrso->fillH2D(hcTrFs.rig[0], sclx * hcrd);
        }
    
        bool jfStatus = hcTrIn.status && (fTof->JFbtaT > 0 && fTof->JFbtaT < 1.0);
        bool hcStatus = hcTrIn.status && (fHyc->btaM1T.status && fHyc->btaM1T.bta[0] < 1.0);

        double jfM = std::sqrt((hcTrIn.rig[0] * hcTrIn.rig[0]) * (1.0 / fTof->JFbtaT / fTof->JFbtaT - 1.0));
        double hcM = std::sqrt((hcTrIn.rig[0] * hcTrIn.rig[0]) * (1.0 / fHyc->btaM1T.bta[0] / fHyc->btaM1T.bta[0] - 1.0));

        if (jfStatus) hJFevt->fillH1D(hcTrIn.rig[0]);
        if (hcStatus) hHCevt->fillH1D(hcTrIn.rig[0]);
        if (jfStatus) hJFtme->fillH1D(hcTrIn.rig[0], fTof->JFT_cpuTime);
        if (hcStatus) hHCtme->fillH1D(hcTrIn.rig[0], fHyc->btaM1T.cpuTime);
        
        if (jfStatus) hJFM->fillH2D(hcTrIn.rig[0], jfM);
        if (hcStatus) hHCM->fillH2D(hcTrIn.rig[0], hcM);
        
        if (fHyc->mutrT.status) hHCM2->fillH2D(fHyc->mutrT.rig[0], fHyc->mutrT.mass);
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
