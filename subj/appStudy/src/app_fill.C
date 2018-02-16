//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/ams_home/hchou/AMSCore/prod/18Feb13/src/ClassDef.h"
#include "/ams_home/hchou/AMSCore/prod/18Feb13/src/ClassDef.C"

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

    //google::InitGoogleLogging(argv[0]);
    //google::ParseCommandLineFlags(&argc, &argv, true);
    //google::SetStderrLogging(google::GLOG_FATAL);

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
    TFile * ofle = new TFile(Form("%s/app_fill%05ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXrig("Rigidity [GV]", 100, 0.5, 4000., AxisScale::kLog);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    
    
    Axis AXtrd("TRD Estimator", 200, 0.0, 1.6);
    Hist * hTRDllr_PR = Hist::New("hTRDllr_PR", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist * hTRDllr_NR = Hist::New("hTRDllr_NR", HistAxis(AXrig, AXtrd, "Events/Bin"));
    
    Axis AXchi("Log Chi-square X", 400, -4.0, 8.0);
    Hist * hTRKchix_PR = Hist::New("hTRKchix_PR", HistAxis(AXrig, AXchi, "Events/Bin"));
    Hist * hTRKchix_NR = Hist::New("hTRKchix_NR", HistAxis(AXrig, AXchi, "Events/Bin"));
    Hist * hTRKchiy_PR = Hist::New("hTRKchiy_PR", HistAxis(AXrig, AXchi, "Events/Bin"));
    Hist * hTRKchiy_NR = Hist::New("hTRKchiy_NR", HistAxis(AXrig, AXchi, "Events/Bin"));
    
    Hist * hCKIRflux = Hist::New("hCKIRflux", HistAxis(AXirig, "Events/Bin"));
    Hist * hCNIRflux = Hist::New("hCNIRflux", HistAxis(AXirig, "Events/Bin"));
    Hist * hKFIRflux = Hist::New("hKFIRflux", HistAxis(AXirig, "Events/Bin"));
    Hist * hHCIRflux = Hist::New("hHCIRflux", HistAxis(AXirig, "Events/Bin"));
    
    Hist * hCKPRflux = Hist::New("hCKPRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hCNPRflux = Hist::New("hCNPRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hKFPRflux = Hist::New("hKFPRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hHCPRflux = Hist::New("hHCPRflux", HistAxis(AXrig, "Events/Bin"));
    
    Hist * hCKNRflux = Hist::New("hCKNRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hCNNRflux = Hist::New("hCNNRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hKFNRflux = Hist::New("hKFNRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hHCNRflux = Hist::New("hHCNRflux", HistAxis(AXrig, "Events/Bin"));
    
    Long64_t printRate = dst->GetEntries();
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        dst->GetEntry(entry);
     
        TrackInfo& track = fTrk->track;
        
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
        if (track.QIn < 0.8 || track.QIn > 1.3) continue;

        // TOF
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        
        if (fTof->numOfInTimeCls > 4) continue;
        if ((fTof->extClsN[0]+fTof->extClsN[1]) > 0 || 
            (fTof->extClsN[2]+fTof->extClsN[3]) > 1) continue;
    
        // TRD
        if (fTrd->LLRph[0] > 0.3) continue;
    
        // IGRF RTI
        if (opt.mode() == MGConfig::JobOpt::MODE::ISS && track.status[0][0]) {
            if (std::fabs(track.rig[0][0]) < 1.2 * fRti->cfIGRF[0]) continue;
        }
       
        Int_t patt = 0;

        Bool_t ck_succ = track.status[0][patt];
        Bool_t cn_succ = track.status[1][patt];
        Bool_t kf_succ = track.status[2][patt];
        Bool_t hc_succ = track.status[3][patt];
        
        Double_t ck_irig = (ck_succ ? MGMath::ONE/track.rig[0][patt] : 0.);
        Double_t cn_irig = (cn_succ ? MGMath::ONE/track.rig[1][patt] : 0.);
        Double_t kf_irig = (kf_succ ? MGMath::ONE/track.rig[2][patt] : 0.);
        Double_t hc_irig = (hc_succ ? MGMath::ONE/track.rig[3][patt] : 0.);
        
        Short_t ck_sign = (ck_succ ? (ck_irig>0?1:-1) : 0);
        Short_t cn_sign = (cn_succ ? (cn_irig>0?1:-1) : 0);
        Short_t kf_sign = (kf_succ ? (kf_irig>0?1:-1) : 0);
        Short_t hc_sign = (hc_succ ? (hc_irig>0?1:-1) : 0);
        
        Double_t ck_rig = (ck_succ ? std::fabs(1.0/ck_irig) : 0.);
        Double_t cn_rig = (cn_succ ? std::fabs(1.0/cn_irig) : 0.);
        Double_t kf_rig = (kf_succ ? std::fabs(1.0/kf_irig) : 0.);
        Double_t hc_rig = (hc_succ ? std::fabs(1.0/hc_irig) : 0.);
        
        Double_t ck_chix = (ck_succ ? std::log(track.chisq[0][patt][0]) : 0.); 
        Double_t cn_chix = (cn_succ ? std::log(track.chisq[1][patt][0]) : 0.); 
        Double_t kf_chix = (kf_succ ? std::log(track.chisq[2][patt][0]) : 0.); 
        Double_t hc_chix = (hc_succ ? std::log(track.chisq[3][patt][0]) : 0.); 
        
        Double_t ck_chiy = (ck_succ ? std::log(track.chisq[0][patt][1]) : 0.); 
        Double_t cn_chiy = (cn_succ ? std::log(track.chisq[1][patt][1]) : 0.); 
        Double_t kf_chiy = (kf_succ ? std::log(track.chisq[2][patt][1]) : 0.); 
        Double_t hc_chiy = (hc_succ ? std::log(track.chisq[3][patt][1]) : 0.); 
        
        if (ck_succ && ck_sign > 0) hTRDllr_PR->fillH2D(ck_rig, fTrd->LLRep[0], fList->weight);
        if (ck_succ && ck_sign < 0) hTRDllr_NR->fillH2D(ck_rig, fTrd->LLRep[0], fList->weight);
        if (ck_succ && ck_sign > 0) hTRKchix_PR->fillH2D(ck_rig, ck_chix, fList->weight);
        if (ck_succ && ck_sign < 0) hTRKchix_NR->fillH2D(ck_rig, ck_chix, fList->weight);
        if (ck_succ && ck_sign > 0) hTRKchiy_PR->fillH2D(ck_rig, ck_chiy, fList->weight);
        if (ck_succ && ck_sign < 0) hTRKchiy_NR->fillH2D(ck_rig, ck_chiy, fList->weight);

        if (ck_succ) hCKIRflux->fillH1D(ck_irig, fList->weight);
        if (cn_succ) hCNIRflux->fillH1D(cn_irig, fList->weight);
        if (kf_succ) hKFIRflux->fillH1D(kf_irig, fList->weight);
        if (hc_succ) hHCIRflux->fillH1D(hc_irig, fList->weight);
        
        if (ck_succ && ck_sign > 0) hCKPRflux->fillH1D(ck_rig, fList->weight);
        if (cn_succ && cn_sign > 0) hCNPRflux->fillH1D(cn_rig, fList->weight);
        if (kf_succ && kf_sign > 0) hKFPRflux->fillH1D(kf_rig, fList->weight);
        if (hc_succ && hc_sign > 0) hHCPRflux->fillH1D(hc_rig, fList->weight);
        
        if (ck_succ && ck_sign < 0) hCKNRflux->fillH1D(ck_rig, fList->weight);
        if (cn_succ && cn_sign < 0) hCNNRflux->fillH1D(cn_rig, fList->weight);
        if (kf_succ && kf_sign < 0) hKFNRflux->fillH1D(kf_rig, fList->weight);
        if (hc_succ && hc_sign < 0) hHCNRflux->fillH1D(hc_rig, fList->weight);
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

    return 0;
}
