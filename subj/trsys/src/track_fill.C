#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

//#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/18Mar23/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18Mar23/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18May19/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18May27/src/ClassDef.h"
#include "/ams_home/hchou/AMSCore/prod/18Jun18/src/ClassDef.h"
//#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/18Jun18/src/ClassDef.h"

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
    
    //TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/afs/cern.ch/work/h/hchou/public/DATABASE/DB/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/afs/cern.ch/work/h/hchou/public/DATABASE/DB/material");

    //TrackSys::Sys::ShowMsg( TrackSys::Sys::GetEnv("TRACKSys_MagBox") );
    //TrackSys::Sys::ShowMsg( TrackSys::Sys::GetEnv("TRACKSys_MatBox") );

    //PhyArg::SetOpt(true, true);
    //PhySt st;
    //st.set_state_with_uxy(0, 0, 55, 0, 0, -1);
    //st.set_mom(3);
    //PropMgnt::PropToZ(-55, st);
    //st.print();

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
    //PartInfo::SetDefault(PartType::Electron);
    PartInfo::SetDefault(PartType::Proton);
    PhyArg::SetOpt(true, true);
    //PhyArg::SetOpt(false, false);
    Bool_t optL1 = false;
    Bool_t optL9 = false;
    
    TFile * ofle = new TFile(Form("%s/track_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 100, 0.5, 1000., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 100, 3.0, 1000., AxisScale::kLog); // RICH AGL
    
    Axis AXrig("Rigidity [GV]", 100, 0.5, 1000., AxisScale::kLog);
    //Axis AXrig("Rigidity [GV]", 100, 3.0, 1000., AxisScale::kLog); // RICH AGL
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    
    Double_t mass = PartInfo(PartType::Proton).mass();
    Axis AXeta("1/GammaBeta [1]", AXmom.nbin(), mass/AXmom.max(), mass/AXmom.min(), AxisScale::kLog);

    Double_t lbta = 1.0/std::sqrt(1.0+AXeta.max()*AXeta.max());
    Double_t ubta = 1.0/std::sqrt(1.0+AXeta.min()*AXeta.min());
    Axis AXbta("Beta [1]", AXmom.nbin(), lbta, ubta, AxisScale::kLog);

    // Time
    Axis AXtme("Time [ms]", 1600, 0., 1000.);
    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXmom, AXtme));
    
    // Fit Eff
    Hist* hCKnum = Hist::New("hCKnum", HistAxis(AXmom, "Events/Bin"));
    Hist* hKFnum = Hist::New("hKFnum", HistAxis(AXmom, "Events/Bin"));
    Hist* hHCnum = Hist::New("hHCnum", HistAxis(AXmom, "Events/Bin"));

    // Fit R Res
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 6000, -4.0, 4.0);
    Hist* hCKRrso = Hist::New("hCKRrso", HistAxis(AXmom, AXRrso));
    Hist* hKFRrso = Hist::New("hKFRrso", HistAxis(AXmom, AXRrso));
    Hist* hHCRrso = Hist::New("hHCRrso", HistAxis(AXmom, AXRrso));
    Hist* hHCRrso2 = Hist::New("hHCRrso2", HistAxis(AXmom, AXRrso));
    
    // Fit B Res
    Axis AXBrso("(Bm/Bt - 1) [1]", 1000, -0.08, 0.08);
    Hist* hCKBrso = Hist::New("hCKBrso", HistAxis(AXbta, AXBrso));
    Hist* hKFBrso = Hist::New("hKFBrso", HistAxis(AXbta, AXBrso));
    Hist* hHCBrso = Hist::New("hHCBrso", HistAxis(AXbta, AXBrso));
    Hist* hHCBrso2 = Hist::New("hHCBrso2", HistAxis(AXbta, AXBrso));
   
    // Fit M 
    Axis AXMrso("Mass [GeV]", 400, 0, 4);
    Hist* hHCMrso = Hist::New("hHCMrso", HistAxis(AXbta, AXMrso));
    
    Axis AXRchi("Log-Chi-square [1]", 800, -3.0, 8.0);
    Hist* hCKRchix = Hist::New("hCKRchix", HistAxis(AXmom, AXRchi));
    Hist* hKFRchix = Hist::New("hKFRchix", HistAxis(AXmom, AXRchi));
    Hist* hHCRchix = Hist::New("hHCRchix", HistAxis(AXmom, AXRchi));
    
    Hist* hCKRchiy = Hist::New("hCKRchiy", HistAxis(AXmom, AXRchi));
    Hist* hKFRchiy = Hist::New("hKFRchiy", HistAxis(AXmom, AXRchi));
    Hist* hHCRchiy = Hist::New("hHCRchiy", HistAxis(AXmom, AXRchi));
    
    //Axis AXRres("Residual [#mum]", 4000, -800.0, 800.0);
    //std::vector<Hist*> hCKresx(9, nullptr);
    //std::vector<Hist*> hKFresx(9, nullptr);
    //std::vector<Hist*> hHCresx(9, nullptr);
    //for (UInt_t it = 0; it < hCKresx.size(); ++it) hCKresx[it] = Hist::New(STR("hCKresxL%d", it+1), HistAxis(AXmom, AXRres));
    //for (UInt_t it = 0; it < hKFresx.size(); ++it) hKFresx[it] = Hist::New(STR("hKFresxL%d", it+1), HistAxis(AXmom, AXRres));
    //for (UInt_t it = 0; it < hHCresx.size(); ++it) hHCresx[it] = Hist::New(STR("hHCresxL%d", it+1), HistAxis(AXmom, AXRres));
    //
    //std::vector<Hist*> hCKresy(9, nullptr);
    //std::vector<Hist*> hKFresy(9, nullptr);
    //std::vector<Hist*> hHCresy(9, nullptr);
    //for (UInt_t it = 0; it < hCKresy.size(); ++it) hCKresy[it] = Hist::New(STR("hCKresyL%d", it+1), HistAxis(AXmom, AXRres));
    //for (UInt_t it = 0; it < hKFresy.size(); ++it) hKFresy[it] = Hist::New(STR("hKFresyL%d", it+1), HistAxis(AXmom, AXRres));
    //for (UInt_t it = 0; it < hHCresy.size(); ++it) hHCresy[it] = Hist::New(STR("hHCresyL%d", it+1), HistAxis(AXmom, AXRres));
    //
    //Axis AXRcos("Residual [10^{-4}]", 3000, -30, 30);
    //std::vector<Hist*> hCKcosx(9, nullptr);
    //std::vector<Hist*> hKFcosx(9, nullptr);
    //std::vector<Hist*> hHCcosx(9, nullptr);
    //for (UInt_t it = 0; it < hCKcosx.size(); ++it) hCKcosx[it] = Hist::New(STR("hCKcosxL%d", it+1), HistAxis(AXmom, AXRcos));
    //for (UInt_t it = 0; it < hKFcosx.size(); ++it) hKFcosx[it] = Hist::New(STR("hKFcosxL%d", it+1), HistAxis(AXmom, AXRcos));
    //for (UInt_t it = 0; it < hHCcosx.size(); ++it) hHCcosx[it] = Hist::New(STR("hHCcosxL%d", it+1), HistAxis(AXmom, AXRcos));
    //
    //std::vector<Hist*> hCKcosy(9, nullptr);
    //std::vector<Hist*> hKFcosy(9, nullptr);
    //std::vector<Hist*> hHCcosy(9, nullptr);
    //for (UInt_t it = 0; it < hCKcosy.size(); ++it) hCKcosy[it] = Hist::New(STR("hCKcosyL%d", it+1), HistAxis(AXmom, AXRcos));
    //for (UInt_t it = 0; it < hKFcosy.size(); ++it) hKFcosy[it] = Hist::New(STR("hKFcosyL%d", it+1), HistAxis(AXmom, AXRcos));
    //for (UInt_t it = 0; it < hHCcosy.size(); ++it) hHCcosy[it] = Hist::New(STR("hHCcosyL%d", it+1), HistAxis(AXmom, AXRcos));

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
        HCTrackInfo& hcTr = fTrk->hcPrInTr.at(0);
        
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

        Bool_t hasMCL1 = false;
        Bool_t hasMCL9 = false;
        for (auto&& mchit : fG4mc->primPart.hits) {
            if (mchit.layJ == 1) hasMCL1 = true;
            if (mchit.layJ == 9) hasMCL9 = true;
        }
       
        Int_t topLay = 1000;
        Bool_t hasL1 = false;
        Bool_t hasL9 = false;
        TrFitPar fitPar(PartType::Proton);
        //TrFitPar fitPar(PartType::Electron);
        for (auto&& hit : fTrk->hits) {
            HitStTRK mhit(hit.side[0], hit.side[1], hit.layJ);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
            mhit.set_nsr(hit.nsr[0], hit.nsr[1]);
            mhit.set_q(hit.adc[0], hit.adc[1]);
         
            if (hit.layJ >= 2 && hit.layJ <= 8) { fitPar.add_hit(mhit); topLay = std::min(topLay, hit.layJ-1); }
            else {
                if (optL1 && hit.layJ == 1) { hasL1 = true; fitPar.add_hit(mhit); topLay = 0; }
                if (optL9 && hit.layJ == 9) { hasL9 = true; fitPar.add_hit(mhit); }
            }
        }
        Short_t cutNHit = 4 + optL1 + optL9;

        Int_t cntCX = 0;
        Int_t cntCY = 0;
        for (auto&& hit : fitPar.hitsTRK()) { cntCX += hit.scx(); cntCY += hit.scy(); }
        if (cntCX < 4) continue;
        if (cntCY < 5) continue;

        for (Int_t il = 0; il < 4; ++il) {
            HitStTOF mhit(il);
            mhit.set_coo(fTof->coo[il][0], fTof->coo[il][1], fTof->coo[il][2]);
            mhit.set_q(fTof->Q[il]);
            mhit.set_t(fTof->T[il]*HitStTOF::TRANS_NS_TO_CM);
            fitPar.add_hit(mhit);
        }

        //if (!fRich->status) continue;
        //if (fRich->kind != 0) continue;
        //HitStRICH richHit( (fRich->kind == 0 ? HitStRICH::Radiator::AGL : HitStRICH::Radiator::NAF) );
        //richHit.set_coo(Numc::ZERO<>, Numc::ZERO<>, fRich->refz);
        //richHit.set_ib(Numc::ONE<> / fRich->beta);
        //fitPar.add_hit(richHit);

        //if (fTrd->hits[0].size() >= 1) {
        //    std::vector<std::pair<Double_t, std::pair<Int_t, Double_t>>> sigs;
        //    for (auto&& hit : fTrd->hits[0]) {
        //        if (hit.len < 0.3) continue;
        //        sigs.push_back(std::make_pair(static_cast<Double_t>(hit.dEdx), std::make_pair(hit.lay, static_cast<Double_t>(hit.coo[2]))));
        //    }
        //    std::sort(sigs.begin(), sigs.end());
        //    if (sigs.size() <= 5) continue;
        //    for (UInt_t it = 2; it < sigs.size(); ++it) {
        //        HitStTRD trdHit(sigs.at(it).second.first);
        //        trdHit.set_coo(0, 0, sigs.at(it).second.second);
        //        trdHit.set_el(sigs.at(it).first);
        //        //fitPar.add_hit(trdHit);
        //    }
        //}
        //else continue;

        if (!fitPar.check()) continue;

        if (optL1 && !(hasL1 && hasMCL1)) continue;
        if (optL9 && !(hasL9 && hasMCL9)) continue;
        Short_t patt = (optL1 + optL9 * 2);

        SegPARTMCInfo* mcs[9] = { nullptr };
        HitTRKMCInfo*  mch[9] = { nullptr };
        HitTRKInfo*    msh[9] = { nullptr };
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec == 0) mcs[seg.lay] = &seg; }
        for (auto&& hit : fG4mc->primPart.hits) mch[hit.layJ-1] = &hit;
        for (auto&& hit :           fTrk->hits) msh[hit.layJ-1] = &hit;

        Bool_t hasLay[9] = { false };
        for (Int_t it = 0; it < 9; ++it) hasLay[it] = (mcs[it] && mch[it] && msh[it]);

        SegPARTMCInfo* topmc = nullptr;
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec == 0 && seg.lay == topLay) { topmc = &seg; break; } }
        if (topmc == nullptr) continue;

        Double_t tmom = 0;
        Double_t tbta = fTof->mcBeta[0];
        if (tbta > 0) {
            tmom = mass / std::sqrt(1.0 / tbta / tbta - 1.0);
            if (!Numc::Valid(tmom)) tmom = topmc->mom;
        }
        else tmom = topmc->mom;
        Double_t mc_mom  = tmom;

        mc_mom = fG4mc->primPart.mom;
        //Double_t mc_mom  = topmc->mom; // Layer 2
        Double_t mc_eta  = mass/mc_mom;
        Double_t mc_bta  = 1.0/std::sqrt(1.0+mc_eta*mc_eta);
        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t bincen  = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
    
        //if (mc_mom < 1.0 || mc_mom > 10.0) continue; // testcode
        //if (mc_mom > 0.8) continue; // testcode
        //if (mc_mom < 300.0) continue; // testcode
        //-------------------------------------//
        MGClock::HrsStopwatch sw; sw.start();
        //PhyTrFit tr(fitPar, PhyTrFit::MuOpt::kFixed);
        PhyTrFit tr(fitPar, PhyTrFit::MuOpt::kFree);
        sw.stop();
        Bool_t hc_succ = tr.status();
        Double_t hc_irig = tr.part().irig();
        Double_t hc_tme  = sw.time()*1.0e3;
        //Double_t hc_coo[9][3]; std::fill_n(hc_coo[0], 9*3, 0.);
        //Double_t hc_dir[9][2]; std::fill_n(hc_dir[0], 9*2, 0.);
        //Double_t hc_lay_irig[9]; std::fill_n(hc_lay_irig, 9, 0.);
        //for (Int_t it = 0; hc_succ && it < 9; ++it) {
        //    PhySt&& stt = tr.interpolate_to_z(ckTr.stateLJ[it][2]);
        //    hc_coo[it][0] = stt.cx();
        //    hc_coo[it][1] = stt.cy();
        //    hc_coo[it][2] = stt.cz();
        //    hc_dir[it][0] = stt.ux();
        //    hc_dir[it][1] = stt.uy();
        //    hc_lay_irig[it] = stt.irig();
        //    //CERR("Lay%d Z %6.2f RIG %14.8f\n", it, hc_coo[it][2], 1.0/hc_lay_irig[it]);
        //}
        //hc_irig = hc_lay_irig[topLay];
       
        //CERR("FINAL FIT (MC MOM %14.8f) == RIG %14.8f MASS %14.8f QLT %14.8f\n", mc_mom, 1.0/hc_irig, tr.part().info().mass(), tr.quality(1));
        PhySt&& sttTop = tr.interpolate_to_z(195.0);
        if (Numc::EqualToZero(sttTop.mom())) continue;
        hc_irig = sttTop.irig();
        //CERR("FINAL FIT (MC MOM %14.8f) == RIG %14.8f MASS %14.8f QLT %14.8f TIME %14.8f  (Z %6.1f)\n", mc_mom, 1.0/hc_irig, tr.part().info().mass(), tr.quality(1), sw.time(), tr.part().cz());
        //CERR("FINAL FIT (MC MOM %14.8f) == RIG %14.8f MASS %14.8f QLT %14.8f\n", mc_mom, 1.0/hc_irig, tr.part().info().mass(), tr.quality(1));
        //-------------------------------------//
        
        Bool_t ck_succ = ckTr.status;
        Bool_t kf_succ = kfTr.status;
        //Bool_t hc_succ = hcTr.status;

        if (ck_succ) hCKnum->fillH1D(mc_mom);
        if (kf_succ) hKFnum->fillH1D(mc_mom);
        if (hc_succ) hHCnum->fillH1D(mc_mom);
        
        //if (hc_succ) hHCtme->fillH2D(mc_mom, hcTr.cpuTime);
        if (hc_succ) hHCtme->fillH2D(mc_mom, hc_tme);

        Double_t ck_irig = (ck_succ ? MGMath::ONE/ckTr.rig : 0.);
        //Double_t kf_irig = (kf_succ ? MGMath::ONE/kfTr.stateTop[6] : 0.);
        Double_t kf_irig = (kf_succ ? MGMath::ONE/kfTr.rig[0] : 0.);
        //Double_t hc_irig = (hc_succ ? MGMath::ONE/hcTr.state[6] : 0.);
        
        Double_t ck_bta = (ck_succ ? 1.0/std::sqrt(1.0+mass*ck_irig*mass*ck_irig) : 0.);
        Double_t kf_bta = (kf_succ ? 1.0/std::sqrt(1.0+mass*kf_irig*mass*kf_irig) : 0.);
        Double_t hc_bta = (hc_succ ? 1.0/std::sqrt(1.0+mass*hc_irig*mass*hc_irig) : 0.);
      
        if (hc_succ) hHCMrso->fillH2D(mc_bta, tr.part().mass());
        //if (hc_succ) hHCMrso->fillH2D(mc_bta, 1.0/hcTr.mass);

        Double_t ck_chix = (ck_succ ? std::log(ckTr.nchi[0]) : 0.); 
        Double_t kf_chix = (kf_succ ? std::log(kfTr.nchi[0]) : 0.); 
        //Double_t hc_chix = (hc_succ ? std::log(hcTr.nchi[0]) : 0.); 
        Double_t hc_chix = (hc_succ ? tr.quality(0) : 0.); 
        
        Double_t ck_chiy = (ck_succ ? std::log(ckTr.nchi[1]) : 0.); 
        Double_t kf_chiy = (kf_succ ? std::log(kfTr.nchi[1]) : 0.); 
        //Double_t hc_chiy = (hc_succ ? std::log(hcTr.nchi[1]) : 0.); 
        Double_t hc_chiy = (hc_succ ? tr.quality(1) : 0.); 
        
        if (ck_succ) hCKRrso->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (kf_succ) hKFRrso->fillH2D(mc_mom, bincen * (kf_irig - mc_irig));
        if (hc_succ) hHCRrso->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (hc_succ && tr.quality(1) < 1.0) hHCRrso2->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_succ) hCKBrso->fillH2D(mc_bta, (ck_bta/mc_bta - 1.0));
        if (kf_succ) hKFBrso->fillH2D(mc_bta, (kf_bta/mc_bta - 1.0));
        if (hc_succ) hHCBrso->fillH2D(mc_bta, (hc_bta/mc_bta - 1.0));
        
        if (hc_succ && tr.quality(1) < 1.0) hHCBrso2->fillH2D(mc_bta, (hc_bta/mc_bta - 1.0));
        
        if (ck_succ) hCKRchix->fillH2D(mc_mom, ck_chix);
        if (kf_succ) hKFRchix->fillH2D(mc_mom, kf_chix);
        if (hc_succ) hHCRchix->fillH2D(mc_mom, hc_chix);
        
        if (ck_succ) hCKRchiy->fillH2D(mc_mom, ck_chiy);
        if (kf_succ) hKFRchiy->fillH2D(mc_mom, kf_chiy);
        if (hc_succ) hHCRchiy->fillH2D(mc_mom, hc_chiy);
        
        //for (Int_t it = 0; it < 9; ++it) {
        //    if (!hasLay[it]) continue;
        //    if (!(msh[it]->side[0] && msh[it]->adc[0]>0)) continue;
        //    if (!(msh[it]->side[1] && msh[it]->adc[1]>0)) continue;
        //    constexpr Double_t CM2UM = 1.0e4;
        //    constexpr Double_t COS   = 1.0e4;
        //    
        //    Double_t ck_resx = (ck_succ ? CM2UM * (track.stateLJ[0][patt][it][0]-mch[it]->coo[0]) : 0.);
        //    Double_t kf_resx = (kf_succ ? CM2UM * (track.stateLJ[1][patt][it][0]-mch[it]->coo[0]) : 0.);
        //    Double_t hc_resx = (hc_succ ? CM2UM * (                hc_coo[it][0]-mch[it]->coo[0]) : 0.);

        //    if (ck_succ) hCKresx.at(it)->fillH2D(mc_mom, ck_resx);
        //    if (kf_succ) hKFresx.at(it)->fillH2D(mc_mom, kf_resx);
        //    if (hc_succ) hHCresx.at(it)->fillH2D(mc_mom, hc_resx);
        //    
        //    Double_t ck_resy = (ck_succ ? CM2UM * (track.stateLJ[0][patt][it][1]-mch[it]->coo[1]) : 0.);
        //    Double_t kf_resy = (kf_succ ? CM2UM * (track.stateLJ[1][patt][it][1]-mch[it]->coo[1]) : 0.);
        //    Double_t hc_resy = (hc_succ ? CM2UM * (                hc_coo[it][1]-mch[it]->coo[1]) : 0.);
        //    
        //    if (ck_succ) hCKresy.at(it)->fillH2D(mc_mom, ck_resy);
        //    if (kf_succ) hKFresy.at(it)->fillH2D(mc_mom, kf_resy);
        //    if (hc_succ) hHCresy.at(it)->fillH2D(mc_mom, hc_resy);
        //    
        //    Double_t ck_cosx = (ck_succ ? COS * (track.stateLJ[0][patt][it][3]-mcs[it]->dir[0]) : 0.);
        //    Double_t kf_cosx = (kf_succ ? COS * (track.stateLJ[1][patt][it][3]-mcs[it]->dir[0]) : 0.);
        //    Double_t hc_cosx = (hc_succ ? COS * (                hc_dir[it][0]-mcs[it]->dir[0]) : 0.);
        //    
        //    if (ck_succ) hCKcosx.at(it)->fillH2D(mc_mom, ck_cosx);
        //    if (kf_succ) hKFcosx.at(it)->fillH2D(mc_mom, kf_cosx);
        //    if (hc_succ) hHCcosx.at(it)->fillH2D(mc_mom, hc_cosx);
        //    
        //    Double_t ck_cosy = (ck_succ ? COS * (track.stateLJ[0][patt][it][4]-mcs[it]->dir[1]) : 0.);
        //    Double_t kf_cosy = (kf_succ ? COS * (track.stateLJ[1][patt][it][4]-mcs[it]->dir[1]) : 0.);
        //    Double_t hc_cosy = (hc_succ ? COS * (                hc_dir[it][1]-mcs[it]->dir[1]) : 0.);
        //    
        //    if (ck_succ) hCKcosy.at(it)->fillH2D(mc_mom, ck_cosy);
        //    if (kf_succ) hKFcosy.at(it)->fillH2D(mc_mom, kf_cosy);
        //    if (hc_succ) hHCcosy.at(it)->fillH2D(mc_mom, hc_cosy);
        //}
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
