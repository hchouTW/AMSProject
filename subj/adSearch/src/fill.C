#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

//#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/18May15/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18May19/src/ClassDef.h"
#include "/ams_home/hchou/AMSCore/prod/18May23/src/ClassDef.h"

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
    dst->SetBranchAddress("ecal", &fEcal);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    TFile * ofle = new TFile(Form("%s/fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    //Axis AXrig("Rigidity [GV]", 100, 0.5, 4000., AxisScale::kLog);
    Axis AXrig("Rigidity [GV]", 100, 0.5, 5., AxisScale::kLog);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);

    // Fit Eff
    Hist* hMCnum = Hist::New("hMCnum", HistAxis(AXrig, "Events/Bin"));
    Hist* hCKnum = Hist::New("hCKnum", HistAxis(AXrig, "Events/Bin"));
    Hist* hKFnum = Hist::New("hKFnum", HistAxis(AXrig, "Events/Bin"));
    Hist* hHCnum = Hist::New("hHCnum", HistAxis(AXrig, "Events/Bin"));
    
    Axis AXMass("Mass", 100, 0., 4.);
    Hist* hMpos = Hist::New("hMpos", HistAxis(AXrig, AXMass));
    Hist* hMneg = Hist::New("hMneg", HistAxis(AXrig, AXMass));
    
    Axis AXQlt("Quality", 100, -2., 8.);
    Hist* hQpos = Hist::New("hQpos", HistAxis(AXrig, AXQlt));
    Hist* hQneg = Hist::New("hQneg", HistAxis(AXrig, AXQlt));
    
    Hist* hMQpos = Hist::New("hMQpos", HistAxis(AXMass, AXQlt));
    Hist* hMQneg = Hist::New("hMQneg", HistAxis(AXMass, AXQlt));
    
    Hist* hM2pos = Hist::New("hM2pos", HistAxis(AXrig, AXMass));
    Hist* hM2neg = Hist::New("hM2neg", HistAxis(AXrig, AXMass));
    
    Hist* hQ2pos = Hist::New("hQ2pos", HistAxis(AXrig, AXQlt));
    Hist* hQ2neg = Hist::New("hQ2neg", HistAxis(AXrig, AXQlt));
    
    Hist* hMQ2pos = Hist::New("hMQ2pos", HistAxis(AXMass, AXQlt));
    Hist* hMQ2neg = Hist::New("hMQ2neg", HistAxis(AXMass, AXQlt));
    
    Axis AXMres("res", 100, -6., 6.);
    Hist* hMres = Hist::New("hMres", HistAxis(AXrig, AXMres));
    
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
        HCTrackInfo& hcTr = fTrk->hcTr.at(0);
        HCTrackInfo& hcMu = fTrk->hcMu;
    
        // Reweight (MC)
        Double_t wgt = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? 1.0 : std::pow(fG4mc->primPart.mom/AXrig.min(), -1.7));

        // Geometry (TOF)
        if (fTof->numOfBetaH != 1) continue;
        if (!fTof->statusBetaH) continue;
        if (fTof->betaHPatt != 15) continue;
       
        // Geometry (TRD)
        if (fTrd->numOfTrack != 1 && fTrd->numOfHTrack != 1) continue;
        if (!fTrd->statusKCls[0]) continue;
        if (fTrd->LLRnhit[0] < 10) continue;
        
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

        if (!ckTr.status && !kfTr.status && !hcTr.status && !hcMu.status) continue;
        
        Short_t ckSign = (ckTr.rig > 0) ? 1 : -1;
        Short_t kfSign = (kfTr.rig[0] > 0) ? 1 : -1;
        Short_t hcSign = (hcMu.stateTop[6] > 0) ? 1 : -1;

        Double_t mom = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? std::fabs(hcMu.stateTop[6]) : fG4mc->primPart.mom);
       
        Double_t TOFibta = 1.0/fTof->betaH;
        Double_t TOFmass = ((TOFibta > 1.0) ? std::sqrt(TOFibta*TOFibta-1.0)*std::fabs(kfTr.rig[1]) : 0.0);

        hMCnum->fillH1D(mom, wgt);
        if (ckSign > 0) hCKnum->fillH1D(ckTr.rig, wgt);
        if (kfSign > 0) hKFnum->fillH1D(kfTr.rig[0], wgt);
        if (hcTr.stateTop[6] > 0) hHCnum->fillH1D(hcTr.stateTop[6], wgt);

        if (hcSign > 0) hMpos->fillH2D(mom, hcMu.mass, wgt);
        if (hcSign < 0) hMneg->fillH2D(mom, hcMu.mass, wgt);
        
        if (hcSign > 0) hQpos->fillH2D(mom, hcMu.quality[1], wgt);
        if (hcSign < 0) hQneg->fillH2D(mom, hcMu.quality[1], wgt);
        
        if (hcSign > 0 && mom < 3.0) hMQpos->fillH2D(hcMu.mass, hcMu.quality[1], wgt);
        if (hcSign < 0 && mom < 3.0) hMQneg->fillH2D(hcMu.mass, hcMu.quality[1], wgt);
        
        if (kfSign > 0 && (TOFibta > 1.0)) hM2pos->fillH2D(mom, TOFmass, wgt);
        if (kfSign < 0 && (TOFibta > 1.0)) hM2neg->fillH2D(mom, TOFmass, wgt);
        
        if (kfSign > 0 && (TOFibta > 1.0)) hQ2pos->fillH2D(mom, std::log(kfTr.nchi[1]), wgt);
        if (kfSign < 0 && (TOFibta > 1.0)) hQ2neg->fillH2D(mom, std::log(kfTr.nchi[1]), wgt);
        
        if (kfSign > 0 && mom < 3.0 && (TOFibta > 1.0)) hMQ2pos->fillH2D(TOFmass, std::log(kfTr.nchi[1]), wgt);
        if (kfSign < 0 && mom < 3.0 && (TOFibta > 1.0)) hMQ2neg->fillH2D(TOFmass, std::log(kfTr.nchi[1]), wgt);
        
        if (hcSign > 0) hMres->fillH2D(mom, ((hcMu.mass - hcTr.mass) / hcTr.error[6]), wgt);
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
