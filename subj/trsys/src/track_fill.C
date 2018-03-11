#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/18Feb27/src/ClassDef.h"

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

    google::InitGoogleLogging(argv[0]);

    //TrackSys::Sys::PutEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::PutEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSProject/subj/trsys/AMS02Mag.bin");
    TrackSys::Sys::PutEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");

    TrackSys::MagGeoBoxAms::Output();

    TrackSys::MagMgnt::Load();
    for (Int_t it = 0; it < 50; ++it) {
        MagFld&& mag = MagMgnt::Get(TrackSys::SVecD<3>(0., 0., it*8-200));
        COUT("%14.8f %14.8f %14.8f\n", mag.x(), mag.y(), mag.z());
    }

    TrackSys::Sys::ShowMsg( TrackSys::Sys::GetEnv("TRACKSys_MagBox") );
    TrackSys::Sys::ShowMsg( TrackSys::Sys::GetEnv("TRACKSys_MatBox") );

    //TrackSys::MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree("/ams_home/hchou/AMSData/new_material", "/ams_home/hchou/AMSData/material/g4mscan.root");
/*
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
    //dst->SetBranchAddress("rich", &fRich);
    //dst->SetBranchAddress("ecal", &fEcal);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    PartType type = PartType::Proton;
    //PartType type = PartType::Electron;
    PhyArg::SetOpt(true, true);
    Bool_t optL1 = false;
    Bool_t optL9 = false;
    
    TFile * ofle = new TFile(Form("%s/track_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 100, 0.35, 4000., AxisScale::kLog);
    
    Axis AXrig("Rigidity [GV]", 100, 0.35, 4000., AxisScale::kLog);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);

    // Time
    Axis AXtme("Time [ms]", 800, 0., 4.);
    Hist* hCKtme = Hist::New("hCKtme", HistAxis(AXmom, AXtme));
    Hist* hKFtme = Hist::New("hKFtme", HistAxis(AXmom, AXtme));
    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXmom, AXtme));
    
    // Fit Eff
    Hist* hCKnum = Hist::New("hCKnum", HistAxis(AXmom, "Events/Bin"));
    Hist* hKFnum = Hist::New("hKFnum", HistAxis(AXmom, "Events/Bin"));
    Hist* hHCnum = Hist::New("hHCnum", HistAxis(AXmom, "Events/Bin"));

    // Fit Res
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 2000, -10.0, 10.0);
    Hist* hCKRrso = Hist::New("hCKRrso", HistAxis(AXmom, AXRrso));
    Hist* hKFRrso = Hist::New("hKFRrso", HistAxis(AXmom, AXRrso));
    Hist* hHCRrso = Hist::New("hHCRrso", HistAxis(AXmom, AXRrso));
    
    Axis AXRchi("Log-Chi-square [1]", 800, -8.0, 8.0);
    Hist* hCKRchix = Hist::New("hCKRchix", HistAxis(AXmom, AXRchi));
    Hist* hKFRchix = Hist::New("hKFRchix", HistAxis(AXmom, AXRchi));
    Hist* hHCRchix = Hist::New("hHCRchix", HistAxis(AXmom, AXRchi));
    
    Hist* hCKRchiy = Hist::New("hCKRchiy", HistAxis(AXmom, AXRchi));
    Hist* hKFRchiy = Hist::New("hKFRchiy", HistAxis(AXmom, AXRchi));
    Hist* hHCRchiy = Hist::New("hHCRchiy", HistAxis(AXmom, AXRchi));
    
    Axis AXRres("Residual [#mum]", 4000, -800.0, 800.0);
    std::vector<Hist*> hCKresx(9, nullptr);
    std::vector<Hist*> hKFresx(9, nullptr);
    std::vector<Hist*> hHCresx(9, nullptr);
    for (UInt_t it = 0; it < hCKresx.size(); ++it) hCKresx[it] = Hist::New(STR_FMT("hCKresxL%d", it+1), HistAxis(AXmom, AXRres));
    for (UInt_t it = 0; it < hKFresx.size(); ++it) hKFresx[it] = Hist::New(STR_FMT("hKFresxL%d", it+1), HistAxis(AXmom, AXRres));
    for (UInt_t it = 0; it < hHCresx.size(); ++it) hHCresx[it] = Hist::New(STR_FMT("hHCresxL%d", it+1), HistAxis(AXmom, AXRres));
    
    std::vector<Hist*> hCKresy(9, nullptr);
    std::vector<Hist*> hKFresy(9, nullptr);
    std::vector<Hist*> hHCresy(9, nullptr);
    for (UInt_t it = 0; it < hCKresy.size(); ++it) hCKresy[it] = Hist::New(STR_FMT("hCKresyL%d", it+1), HistAxis(AXmom, AXRres));
    for (UInt_t it = 0; it < hKFresy.size(); ++it) hKFresy[it] = Hist::New(STR_FMT("hKFresyL%d", it+1), HistAxis(AXmom, AXRres));
    for (UInt_t it = 0; it < hHCresy.size(); ++it) hHCresy[it] = Hist::New(STR_FMT("hHCresyL%d", it+1), HistAxis(AXmom, AXRres));
    
    Axis AXRcos("Residual [10^{-4}]", 3000, -30, 30);
    std::vector<Hist*> hCKcosx(9, nullptr);
    std::vector<Hist*> hKFcosx(9, nullptr);
    std::vector<Hist*> hHCcosx(9, nullptr);
    for (UInt_t it = 0; it < hCKcosx.size(); ++it) hCKcosx[it] = Hist::New(STR_FMT("hCKcosxL%d", it+1), HistAxis(AXmom, AXRcos));
    for (UInt_t it = 0; it < hKFcosx.size(); ++it) hKFcosx[it] = Hist::New(STR_FMT("hKFcosxL%d", it+1), HistAxis(AXmom, AXRcos));
    for (UInt_t it = 0; it < hHCcosx.size(); ++it) hHCcosx[it] = Hist::New(STR_FMT("hHCcosxL%d", it+1), HistAxis(AXmom, AXRcos));
    
    std::vector<Hist*> hCKcosy(9, nullptr);
    std::vector<Hist*> hKFcosy(9, nullptr);
    std::vector<Hist*> hHCcosy(9, nullptr);
    for (UInt_t it = 0; it < hCKcosy.size(); ++it) hCKcosy[it] = Hist::New(STR_FMT("hCKcosyL%d", it+1), HistAxis(AXmom, AXRcos));
    for (UInt_t it = 0; it < hKFcosy.size(); ++it) hKFcosy[it] = Hist::New(STR_FMT("hKFcosyL%d", it+1), HistAxis(AXmom, AXRcos));
    for (UInt_t it = 0; it < hHCcosy.size(); ++it) hHCcosy[it] = Hist::New(STR_FMT("hHCcosyL%d", it+1), HistAxis(AXmom, AXRcos));

    Long64_t printRate = static_cast<Long64_t>(0.02 * dst->GetEntries());
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

        Bool_t hasMCL1 = false;
        Bool_t hasMCL9 = false;
        for (auto&& mchit : fG4mc->primPart.hits) {
            if (mchit.layJ == 1) hasMCL1 = true;
            if (mchit.layJ == 9) hasMCL9 = true;
        }
        
        Bool_t hasL1 = false;
        Bool_t hasL9 = false;
        TrFitPar fitPar(type);
        for (auto&& hit : track.hits) {
            HitSt mhit(hit.side%2==1, hit.side/2==1, hit.layJ);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
            mhit.set_nsr(hit.nsr[0], hit.nsr[1]);
            mhit.set_adc(hit.sig[0], hit.sig[1]);
          
            if (hit.layJ >= 2 && hit.layJ <= 8) fitPar.addHit(mhit);
            else {
                if (optL1 && hit.layJ == 1) { hasL1 = true; fitPar.addHit(mhit); }
                if (optL9 && hit.layJ == 9) { hasL9 = true; fitPar.addHit(mhit); }
            }
        }
        Short_t cutNHit = 4 + optL1 + optL9;
        if (fitPar.numOfHit() <= cutNHit) continue;

        if (optL1 && !(hasL1 && hasMCL1)) continue;
        if (optL9 && !(hasL9 && hasMCL9)) continue;
        Short_t patt = (optL1 + optL9 * 2);

        SegPARTMCInfo* mcs[9] = { nullptr };
        HitTRKMCInfo*  mch[9] = { nullptr };
        HitTRKInfo*    msh[9] = { nullptr };
        for (auto&& seg : fG4mc->primPart.segs) mcs[seg.lay -1] = &seg;
        for (auto&& hit : fG4mc->primPart.hits) mch[hit.layJ-1] = &hit;
        for (auto&& hit :           track.hits) msh[hit.layJ-1] = &hit;

        Bool_t hasLay[9] = { false };
        for (Int_t it = 0; it < 9; ++it) hasLay[it] = (mcs[it] && mch[it] && msh[it]);

        SegPARTMCInfo* topmc = nullptr;
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec == 0 && seg.lay >= 2-optL1) { topmc = &seg; break; } }
        if (topmc == nullptr) continue;

        Double_t mc_mom  = topmc->mom;
        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t bincen  = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
        
        //-------------------------------------//
        MGClock::HrsStopwatch sw; sw.start();
        PhyTrFit tr(fitPar);
        sw.stop();
        Bool_t hc_succ = tr.status();
        Double_t hc_irig = tr.part().irig();
        Double_t hc_tme  = sw.time()*1.0e3;
        Double_t hc_coo[9][3]; std::fill_n(hc_coo[0], 9*3, 0.);
        Double_t hc_dir[9][2]; std::fill_n(hc_dir[0], 9*2, 0.);
        for (Int_t it = 0; it < 9; ++it) {
            const PhySt* stt = tr.stts(it+1);
            if (stt == nullptr) continue;
            hc_coo[it][0] = stt->cx();
            hc_coo[it][1] = stt->cy();
            hc_coo[it][2] = stt->cz();
            hc_dir[it][0] = stt->ux();
            hc_dir[it][1] = stt->uy();
        }
        //-------------------------------------//
        
        Bool_t ck_succ = track.status[0][patt];
        Bool_t kf_succ = track.status[2][patt];
        //Bool_t hc_succ = track.status[3][patt];

        if (ck_succ) hCKnum->fillH1D(mc_mom);
        if (kf_succ) hKFnum->fillH1D(mc_mom);
        if (hc_succ) hHCnum->fillH1D(mc_mom);
        
        if (ck_succ) hCKtme->fillH2D(mc_mom, track.cpuTime[0][patt]);
        if (kf_succ) hKFtme->fillH2D(mc_mom, track.cpuTime[2][patt]*0.01);
        //if (hc_succ) hHCtme->fillH2D(mc_mom, track.cpuTime[3][patt]*0.1);
        if (hc_succ) hHCtme->fillH2D(mc_mom, hc_tme*0.1);

        Double_t ck_irig = (ck_succ ? MGMath::ONE/track.rig[0][patt] : 0.);
        Double_t kf_irig = (kf_succ ? MGMath::ONE/track.rig[2][patt] : 0.);
        //Double_t hc_irig = (hc_succ ? MGMath::ONE/track.rig[3][patt] : 0.);

        Double_t ck_rig = (ck_succ ? 1.0/ck_irig : 0.);
        Double_t kf_rig = (kf_succ ? 1.0/kf_irig : 0.);
        Double_t hc_rig = (hc_succ ? 1.0/hc_irig : 0.);
        
        Double_t ck_chix = (ck_succ ? std::log(track.chisq[0][patt][0]) : 0.); 
        Double_t kf_chix = (kf_succ ? std::log(track.chisq[2][patt][0]) : 0.); 
        //Double_t hc_chix = (hc_succ ? std::log(track.chisq[3][patt][0]) : 0.); 
        Double_t hc_chix = (hc_succ ? std::log(tr.nchix())              : 0.); 
        
        Double_t ck_chiy = (ck_succ ? std::log(track.chisq[0][patt][1]) : 0.); 
        Double_t kf_chiy = (kf_succ ? std::log(track.chisq[2][patt][1]) : 0.); 
        //Double_t hc_chiy = (hc_succ ? std::log(track.chisq[3][patt][1]) : 0.); 
        Double_t hc_chiy = (hc_succ ? std::log(tr.nchiy())              : 0.); 
        
        if (ck_succ) hCKRrso->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (kf_succ) hKFRrso->fillH2D(mc_mom, bincen * (kf_irig - mc_irig));
        if (hc_succ) hHCRrso->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_succ) hCKRchix->fillH2D(mc_mom, ck_chix);
        if (kf_succ) hKFRchix->fillH2D(mc_mom, kf_chix);
        if (hc_succ) hHCRchix->fillH2D(mc_mom, hc_chix);
        
        if (ck_succ) hCKRchiy->fillH2D(mc_mom, ck_chiy);
        if (kf_succ) hKFRchiy->fillH2D(mc_mom, kf_chiy);
        if (hc_succ) hHCRchiy->fillH2D(mc_mom, hc_chiy);

        for (Int_t it = 0; it < 9; ++it) {
            if (!hasLay[it]) continue;
            if (msh[it]->side!=3) continue;
            constexpr Double_t CM2UM = 1.0e4;
            constexpr Double_t COS   = 1.0e4;
            
            Double_t ck_resx = (ck_succ ? CM2UM * (track.stateLJ[0][patt][it][0]-mch[it]->coo[0]) : 0.);
            Double_t kf_resx = (kf_succ ? CM2UM * (track.stateLJ[2][patt][it][0]-mch[it]->coo[0]) : 0.);
            Double_t hc_resx = (hc_succ ? CM2UM * (                hc_coo[it][0]-mch[it]->coo[0]) : 0.);

            if (ck_succ) hCKresx.at(it)->fillH2D(mc_mom, ck_resx);
            if (kf_succ) hKFresx.at(it)->fillH2D(mc_mom, kf_resx);
            if (hc_succ) hHCresx.at(it)->fillH2D(mc_mom, hc_resx);
            
            Double_t ck_resy = (ck_succ ? CM2UM * (track.stateLJ[0][patt][it][1]-mch[it]->coo[1]) : 0.);
            Double_t kf_resy = (kf_succ ? CM2UM * (track.stateLJ[2][patt][it][1]-mch[it]->coo[1]) : 0.);
            Double_t hc_resy = (hc_succ ? CM2UM * (                hc_coo[it][1]-mch[it]->coo[1]) : 0.);
            
            if (ck_succ) hCKresy.at(it)->fillH2D(mc_mom, ck_resy);
            if (kf_succ) hKFresy.at(it)->fillH2D(mc_mom, kf_resy);
            if (hc_succ) hHCresy.at(it)->fillH2D(mc_mom, hc_resy);
            
            Double_t ck_cosx = (ck_succ ? COS * (track.stateLJ[0][patt][it][3]-mcs[it]->dir[0]) : 0.);
            Double_t kf_cosx = (kf_succ ? COS * (track.stateLJ[2][patt][it][3]-mcs[it]->dir[0]) : 0.);
            Double_t hc_cosx = (hc_succ ? COS * (                hc_dir[it][0]-mcs[it]->dir[0]) : 0.);
            
            if (ck_succ) hCKcosx.at(it)->fillH2D(mc_mom, ck_cosx);
            if (kf_succ) hKFcosx.at(it)->fillH2D(mc_mom, kf_cosx);
            if (hc_succ) hHCcosx.at(it)->fillH2D(mc_mom, hc_cosx);
            
            Double_t ck_cosy = (ck_succ ? COS * (track.stateLJ[0][patt][it][4]-mcs[it]->dir[1]) : 0.);
            Double_t kf_cosy = (kf_succ ? COS * (track.stateLJ[2][patt][it][4]-mcs[it]->dir[1]) : 0.);
            Double_t hc_cosy = (hc_succ ? COS * (                hc_dir[it][1]-mcs[it]->dir[1]) : 0.);
            
            if (ck_succ) hCKcosy.at(it)->fillH2D(mc_mom, ck_cosy);
            if (kf_succ) hKFcosy.at(it)->fillH2D(mc_mom, kf_cosy);
            if (hc_succ) hHCcosy.at(it)->fillH2D(mc_mom, hc_cosy);
        }
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
    */
    return 0;
}
