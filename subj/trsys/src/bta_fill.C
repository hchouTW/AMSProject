#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/18Oct17/src/ClassDef.h"

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
    TFile * ofle = new TFile(Form("%s/bta_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    PartInfo info(PartType::Proton);
    //PartInfo info(PartType::Helium4);
    
    PartInfo::SetDefault(info.type());
    PhyArg::SetOpt(true, true);
    Bool_t optL1 = false;
    Bool_t optL9 = false;
    
    Int_t    nmom     = 80;
    Double_t mombd[2] = { 1., 1000. };
    if (info.type() == PartType::Proton)   { mombd[0] = 0.55; mombd[1] = 20.0; }
    if (info.type() == PartType::Helium4)  { mombd[0] = 2.20; mombd[1] = 40.0; }
    Axis AXmom("Momentum [GeV]", nmom, mombd[0], mombd[1], AxisScale::kLog);
    
    Axis AXrig("Rigidity [GV]", AXmom.nbin(), mombd[0]/std::fabs(info.chrg()), mombd[1]/std::fabs(info.chrg()), AxisScale::kLog);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    
    Axis AXigb("1/GammaBeta [1]", AXmom.nbin(), info.mass()/AXmom.max(), info.mass()/AXmom.min(), AxisScale::kLog);
    Double_t lbta = std::sqrt(1.0+AXigb.min()*AXigb.min());
    Double_t ubta = std::sqrt(1.0+AXigb.max()*AXigb.max());
    Axis AXib("1/Beta [1]", AXigb.nbin(), lbta, ubta, AxisScale::kLog);
    Axis AXbta("Beta [1]", AXigb.nbin(), 1.0/ubta, 1.0/lbta, AxisScale::kLog);

    // Time
    Axis AXtme("Time [ms]", 1600, 0., 1000.);
    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXmom, AXtme));
    
    // Fit R Res
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 2000, -1.8, 1.8);
    Hist* hCKRrso = Hist::New("hCKRrso", HistAxis(AXmom, AXRrso));
    Hist* hKFRrso = Hist::New("hKFRrso", HistAxis(AXmom, AXRrso));
    Hist* hHCRrso = Hist::New("hHCRrso", HistAxis(AXmom, AXRrso));
    
    // Fit B Res
    Axis AXTFBrso("(Bm/Bt - 1) [1]", 1000, -0.20, 0.20);
    Axis AXRHBrso("(Bm/Bt - 1) [1]", 1000, -0.05, 0.05);
    Hist* hTFBrso   = Hist::New("hTFBrso", HistAxis(AXmom, AXTFBrso));
    Hist* hRHBrso   = Hist::New("hRHBrso", HistAxis(AXmom, AXRHBrso));
    Hist* hHCTFBrso = Hist::New("hHCTFBrso", HistAxis(AXmom, AXTFBrso));
    Hist* hHCRHBrso = Hist::New("hHCRHBrso", HistAxis(AXmom, AXRHBrso));
    
    Axis AXBchi("Log-Chi-square [1]", 800, -3.0, 8.0);
    Hist* hHCBchi = Hist::New("hHCBchi", HistAxis(AXmom, AXBchi));
    
    // Fit M Res
    Axis AXM("Mass", 1000, 0., 10.);
    Hist* hCKM = Hist::New("hCKM", HistAxis(AXmom, AXM));
    Hist* hKFM = Hist::New("hKFM", HistAxis(AXmom, AXM));
    Hist* hHCM = Hist::New("hHCM", HistAxis(AXmom, AXM));
    
    Axis AXMchi("Log-Chi-square [1]", 800, -3.0, 5.0);
    Hist* hHCMchix = Hist::New("hHCMchix", HistAxis(AXmom, AXMchi));
    Hist* hHCMchiy = Hist::New("hHCMchiy", HistAxis(AXmom, AXMchi));
    Hist* hHCMchib = Hist::New("hHCMchib", HistAxis(AXmom, AXMchi));
    Hist* hHCMchi0 = Hist::New("hHCMchi0", HistAxis(AXmom, AXMchi));
    Hist* hHCMchi1 = Hist::New("hHCMchi1", HistAxis(AXmom, AXMchi));
    
    MGClock::HrsStopwatch hrssw; hrssw.start();
    Long64_t passEntry = 0;
    Long64_t printRate = static_cast<Long64_t>(0.1 * dst->GetEntries());
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) { // testcode
        if (entry%printRate==0) {
            hrssw.stop();
            COUT("Entry (%lld) %lld/%lld Time %14.8f\n", passEntry, entry, dst->GetEntries(), hrssw.time());
        }
        dst->GetEntry(entry);
     
        // No Interaction (testcode)
        //if (opt.mode() == MGConfig::JobOpt::MODE::MC)
        //    if (fG4mc->primVtx.status && fG4mc->primVtx.coo[2] > -120) continue;
        
        //if (entry > 100) break;
        if (fG4mc->primPart.mom > mombd[1]+5.0) continue; // testcode
        //if (fG4mc->primPart.mom > 0.7) continue; // testcode
        
        Int_t trPatt = optL1 + optL9 * 2;
        CKTrackInfo& ckTr = fTrk->ckTr.at(trPatt);
        KFTrackInfo& kfTr = fTrk->kfTr.at(trPatt);
        HCTrackInfo& hcTr = fTrk->hcTr.at(trPatt);
        
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

        Bool_t hasMCL1 = false;
        Bool_t hasMCL9 = false;
        for (auto&& mchit : fG4mc->primPart.hitsTk) {
            if (mchit.layJ == 1) hasMCL1 = true;
            if (mchit.layJ == 9) hasMCL9 = true;
        }
       
        Bool_t hasL1 = false;
        Bool_t hasL9 = false;
        TrFitPar fitPar(info.type());
        TrFitPar btaPar(info.type());
        TrFitPar muPar(info.type());
        for (auto&& hit : fTrk->hits) {
            Bool_t isInnTr = (hit.layJ >= 2 && hit.layJ <= 8);
            HitStTRK mhit(hit.side[0], hit.side[1], hit.layJ, isInnTr);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
         
            if (isInnTr) { fitPar.add_hit(mhit); muPar.add_hit(mhit); }
            else {
                if (optL1 && hit.layJ == 1) { hasL1 = true; fitPar.add_hit(mhit); muPar.add_hit(mhit); }
                if (optL9 && hit.layJ == 9) { hasL9 = true; fitPar.add_hit(mhit); muPar.add_hit(mhit); }
            }
           
            //if (hit.chrg[0] > 0 && hit.chrg[1] > 0 && hit.chrg[2] > 0) {
            //    HitStTRK mhit_q(false, false, hit.layJ, isInnTr);
            //    mhit_q.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
            //    mhit_q.set_q(hit.chrg[2], hit.chrg[0], hit.chrg[1], info.chrg());
            //    if (isInnTr) { btaPar.add_hit(mhit_q); }
            //    else {
            //        if (optL1 && hit.layJ == 1) { btaPar.add_hit(mhit_q); }
            //        if (optL9 && hit.layJ == 9) { btaPar.add_hit(mhit_q); }
            //    }
            //}
        }

        for (Int_t il = 0; il < 4; ++il) {
            HitStTOF mhit(il);
            mhit.set_coo(fTof->coo[il][0], fTof->coo[il][1], fTof->coo[il][2]);
            //mhit.set_q(fTof->Q[il], info.chrg());
            mhit.set_t(fTof->T[il]*HitStTOF::TRANS_NS_TO_CM);
            btaPar.add_hit(mhit);
            muPar.add_hit(mhit);
        }

        if (!fitPar.check()) continue;
        if (optL1 && !(hasL1 && hasMCL1)) continue;
        if (optL9 && !(hasL9 && hasMCL9)) continue;
        Short_t patt = (optL1 + optL9 * 2);

        SegPARTMCInfo* mcs[9] = { nullptr };
        HitTRKMCInfo*  mch[9] = { nullptr };
        HitTRKInfo*    msh[9] = { nullptr };
        for (auto&& seg : fG4mc->primPart.segsTk) { mcs[seg.lay] = &seg; }
        for (auto&& hit : fG4mc->primPart.hitsTk) mch[hit.layJ-1] = &hit;
        for (auto&& hit :             fTrk->hits) msh[hit.layJ-1] = &hit;

        Bool_t hasLay[9] = { false };
        for (Int_t it = 0; it < 9; ++it) hasLay[it] = (mcs[it] && mch[it] && msh[it]);

        Double_t mc_mom  = fG4mc->primPart.mom;
        Double_t mc_igb  = fG4mc->primPart.mass / mc_mom;
        Double_t mc_ibta = std::hypot(1.0, mc_igb);
        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t bincen  = std::sqrt(AXmom.center(AXmom.find(mc_mom), AxisScale::kLog));
        
        //-------------------------------------//
        PhyTrFit tr(fitPar);
        if (!tr.status()) COUT("HC TR FAILURE.\n");
        if (!tr.status()) continue;

        continue; // testcode

        PhyBtaFit trbta(btaPar, tr.part());
        if (!trbta.status()) COUT("HC BTA FAILURE.\n");
        if (!trbta.status()) continue;

        MGClock::HrsStopwatch sw; sw.start();
        PhyMuFit mufit(muPar);
        sw.stop();

        if (!mufit.status()) COUT("HC MU FAILURE.\n");
        if (!mufit.status()) continue;
        
        PhyTrFit trfit = mufit.fit();
        if (!trfit.status()) continue;
       
        //CERR("TIME %14.8f MASS %14.8f QLT %14.8f %14.8f %14.8f\n", sw.time(), mufit.part().mass(), mufit.quality(0), mufit.quality(1), mufit.quality(2));

        Bool_t   hc_succ = tr.status();
        Double_t hc_irig = tr.part().irig();
        Double_t hc_ibta = tr.part().ibta();
        Double_t hc_mass = tr.part().mass();
        Double_t hc_tme  = sw.time()*1.0e3;
  
        PhySt&& sttTop = tr.interpolate_to_z(195.0);
        PhySt&& sttBta = trbta.interpolate_to_z(195.0);
        hc_succ = (hc_succ ? !Numc::EqualToZero(sttTop.mom()) : false);
        if (hc_succ) hc_irig = sttTop.irig();
        if (hc_succ && trbta.status()) hc_ibta = sttBta.ibta();
        hc_mass = (sttTop.mom() * sttBta.igb()) * (sttTop.mom() * sttBta.igb());
        hc_mass = mufit.part().mass();
        //-------------------------------------//
        Bool_t ck_succ = ckTr.status;
        Bool_t kf_succ = kfTr.status;
        
        if (hc_succ) hHCtme->fillH2D(mc_mom, hc_tme);

        Double_t ck_irig = (ck_succ ? MGMath::ONE/ckTr.rig : 0.);
        Double_t kf_irig = (kf_succ ? MGMath::ONE/kfTr.rig[0] : 0.);
      
        if (ck_succ) hCKRrso->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (kf_succ) hKFRrso->fillH2D(mc_mom, bincen * (kf_irig - mc_irig));
        if (hc_succ) hHCRrso->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        Double_t wgt = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? 1.0 : std::pow(fG4mc->primPart.mom/AXrig.min(), -1.7));
        
        Bool_t TFstatus = (fTof->statusBetaH && fTof->betaH < 1.0);
        Bool_t RHstatus = (fRich->status && fRich->kind == 0 && fRich->beta < 1.0);
        Bool_t HCstatus = (hc_succ && trbta.status() && !mufit.is_like_el());
        
        Double_t TFibta = 1.0 / fTof->betaH;
        Double_t RHibta = 1.0 / fRich->beta;
        if (TFstatus) hTFBrso->fillH2D(mc_mom, (mc_ibta/TFibta - 1.0));
        if (RHstatus) hRHBrso->fillH2D(mc_mom, (mc_ibta/RHibta - 1.0));
        if (HCstatus) hHCTFBrso->fillH2D(mc_mom, (mc_ibta/hc_ibta - 1.0));
        if (HCstatus) hHCRHBrso->fillH2D(mc_mom, (mc_ibta/hc_ibta - 1.0));
        
        if (HCstatus) hHCBchi->fillH2D(mc_mom, trbta.quality());
        
        Double_t ck_mass = ((ck_succ && TFstatus) ? (TFibta*TFibta-1.0) * (ckTr.rig * ckTr.rig) : 0.0);
        Double_t kf_mass = ((kf_succ && TFstatus) ? (TFibta*TFibta-1.0) * (kfTr.rig[1] * kfTr.rig[1]) : 0.0);

        if (ck_succ && TFstatus) hCKM->fillH2D(mc_mom, ck_mass);
        if (kf_succ && TFstatus) hKFM->fillH2D(mc_mom, kf_mass);
        if (hc_succ && HCstatus) hHCM->fillH2D(mc_mom, hc_mass);
        
        if (HCstatus) hHCMchix->fillH2D(mc_mom, mufit.quality(0));
        if (HCstatus) hHCMchiy->fillH2D(mc_mom, mufit.quality(1));
        if (HCstatus) hHCMchib->fillH2D(mc_mom, mufit.quality(2));
        
        if (HCstatus) hHCMchi0->fillH2D(mc_mom, trfit.quality(0));
        if (HCstatus) hHCMchi1->fillH2D(mc_mom, trfit.quality(1));

        passEntry++;
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
