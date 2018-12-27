#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/18Dec23/src/ClassDef.h"

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
    TFile * ofle = new TFile(Form("%s/mu_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
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
    
    // Time
    Axis AXtme("Time [ms]", 1600, 0., 1000.);
    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXmom, AXtme));
    
    // Fit R Res
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 2000, -1.8, 1.8);
    Hist* hCKRrso = Hist::New("hCKRrso", HistAxis(AXmom, AXRrso));
    Hist* hHCRrso = Hist::New("hHCRrso", HistAxis(AXmom, AXRrso));
    
    Axis AXRqlt("Quality [1]", 800, -2.0, 4.0);
    Hist* hHCRqltx = Hist::New("hHCRqltx", HistAxis(AXmom, AXRqlt));
    Hist* hHCRqlty = Hist::New("hHCRqlty", HistAxis(AXmom, AXRqlt));
    
    Hist* hHCRrsoPr  = Hist::New("hHCRrsoPr", HistAxis(AXmom, AXRrso));
    Hist* hHCRqltxPr = Hist::New("hHCRqltxPr", HistAxis(AXmom, AXRqlt));
    Hist* hHCRqltyPr = Hist::New("hHCRqltyPr", HistAxis(AXmom, AXRqlt));
    
    Hist* hHCRrsoD  = Hist::New("hHCRrsoD", HistAxis(AXmom, AXRrso));
    Hist* hHCRqltxD = Hist::New("hHCRqltxD", HistAxis(AXmom, AXRqlt));
    Hist* hHCRqltyD = Hist::New("hHCRqltyD", HistAxis(AXmom, AXRqlt));
    
    Axis AXRscat("Scat Quality [1]", 600, 0.0, 6.0);
    Hist* hHCRscatx = Hist::New("hHCRscatx", HistAxis(AXRscat, AXRscat));
    Hist* hHCRscaty = Hist::New("hHCRscaty", HistAxis(AXRscat, AXRscat));
    Hist* hHCRsx = Hist::New("hHCRsx", HistAxis(AXmom, AXRscat));
    Hist* hHCRsy = Hist::New("hHCRsy", HistAxis(AXmom, AXRscat));
    Hist* hHCRst = Hist::New("hHCRst", HistAxis(AXmom, AXRscat));
    Hist* hHCRsr = Hist::New("hHCRsr", HistAxis(AXmom, AXRscat));
    
    Hist* hHCRsxMax = Hist::New("hHCRsxMax", HistAxis(AXmom, AXRscat));
    Hist* hHCRsyMax = Hist::New("hHCRsyMax", HistAxis(AXmom, AXRscat));
    Hist* hHCRstMax = Hist::New("hHCRstMax", HistAxis(AXmom, AXRscat));
    Hist* hHCRsrMax = Hist::New("hHCRsrMax", HistAxis(AXmom, AXRscat));
    
    // Fit B Res
    Axis AXTFBrso("(Bm/Bt - 1) [1]", 1000, -0.20, 0.20);
    Axis AXRHBrso("(Bm/Bt - 1) [1]", 1000, -0.05, 0.05);
    Hist* hTFBrso   = Hist::New("hTFBrso", HistAxis(AXmom, AXTFBrso));
    Hist* hRHBrso   = Hist::New("hRHBrso", HistAxis(AXmom, AXRHBrso));
    Hist* hHCTFBrso = Hist::New("hHCTFBrso", HistAxis(AXmom, AXTFBrso));
    Hist* hHCRHBrso = Hist::New("hHCRHBrso", HistAxis(AXmom, AXRHBrso));
    
    Hist* hHCTFBrsoPr = Hist::New("hHCTFBrsoPr", HistAxis(AXmom, AXTFBrso));
    Hist* hHCRHBrsoPr = Hist::New("hHCRHBrsoPr", HistAxis(AXmom, AXRHBrso));
    
    Axis AXBqlt("Quality [1]", 800, -2.0, 4.0);
    Hist* hHCBqltbPr = Hist::New("hHCBqltbPr", HistAxis(AXmom, AXBqlt));
    
    // Fit M Res
    Axis AXM("Mass", 800, 0.03, 8.0);
    Hist* hCKM = Hist::New("hCKM", HistAxis(AXmom, AXM));
    Hist* hHCM = Hist::New("hHCM", HistAxis(AXmom, AXM));
    
    Hist* hCKMcut = Hist::New("hCKMcut", HistAxis(AXmom, AXM));
    Hist* hHCMcut = Hist::New("hHCMcut", HistAxis(AXmom, AXM));
    
    Axis AXMqlt("Quality [1]", 800, -2.0, 4.0);
    Hist* hHCMqltx = Hist::New("hHCMqltx", HistAxis(AXmom, AXMqlt));
    Hist* hHCMqlty = Hist::New("hHCMqlty", HistAxis(AXmom, AXMqlt));
    Hist* hHCMqltb = Hist::New("hHCMqltb", HistAxis(AXmom, AXMqlt));
    
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
        TrFitPar fitParD(PartType::Deuterium);
        for (auto&& hit : fTrk->hits) {
            Bool_t isInnTr = (hit.layJ >= 2 && hit.layJ <= 8);
            HitStTRK mhit(hit.side[0], hit.side[1], hit.layJ, isInnTr);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
            mhit.set_nsr(hit.nsr[0], hit.nsr[1]);
            //mhit.set_q(hit.chrg[2], hit.chrg[0], hit.chrg[1], info.chrg());
         
            if (isInnTr) { fitPar.add_hit(mhit); fitParD.add_hit(mhit); }
            else {
                if (optL1 && hit.layJ == 1) { hasL1 = true; fitPar.add_hit(mhit); fitParD.add_hit(mhit); }
                if (optL9 && hit.layJ == 9) { hasL9 = true; fitPar.add_hit(mhit); fitParD.add_hit(mhit); }
            }
        }

        for (Int_t il = 0; il < 4; ++il) {
            HitStTOF mhit(il);
            mhit.set_coo(fTof->coo[il][0], fTof->coo[il][1], fTof->coo[il][2]);
            //mhit.set_q(fTof->Q[il], info.chrg());
            mhit.set_t(fTof->T[il]*HitStTOF::TRANS_NS_TO_CM);
            fitPar.add_hit(mhit);
            fitParD.add_hit(mhit);
        }

        if (!fitPar.check()) continue;
        if (!fitParD.check()) continue;
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
        Bool_t   ck_succ = ckTr.status;
        Double_t ck_irig = (ck_succ ? MGMath::ONE/ckTr.rig : 0.0);
        Double_t ck_mom  = (ck_succ ? std::fabs(ckTr.rig * info.chrg()) : 0.0);
        if (!ck_succ) { COUT("CK TR FAILURE.\n"); continue; }

        MGClock::HrsStopwatch sw; sw.start();
        PhyMuFit mufit(fitPar);
        if (!mufit.status()) { COUT("HC MU FAILURE.\n"); continue; }
        PhyTrFit trfit = mufit.fit();
        if (!trfit.status()) { COUT("HC TR FAILURE.\n"); continue; }
        sw.stop();
        
        PhyTrFit prtr(fitPar);
        if (!prtr.status()) { COUT("HC PR TR FAILURE.\n"); continue; }
        
        PhyBtaFit prbta(fitPar, prtr.part());
        if (!prbta.status()) { COUT("HC PR BTA FAILURE.\n"); continue; }
        
        PhyTrFit dtr(fitParD);
        
        //CERR("TIME %14.8f MASS %14.8f QLT %14.8f %14.8f %14.8f\n", sw.time(), mufit.part().mass(), mufit.quality(0), mufit.quality(1), mufit.quality(2));

        Bool_t   hc_succ = trfit.status();
        Double_t hc_irig = trfit.part().irig();
        Double_t hc_ibta = trfit.part().ibta();
        Double_t hc_mass = trfit.part().mass();
        Double_t hc_tme  = sw.time()*1.0e3;
 
        PhySt&& sttTop = trfit.interpolate_to_z(195.0);
        hc_succ = (hc_succ ? !Numc::EqualToZero(sttTop.mom()) : false);
        if (hc_succ) hc_irig = sttTop.irig();
        if (hc_succ) hc_ibta = sttTop.ibta();
        
        if (!hc_succ) { COUT("HC PROP FAILURE.\n"); continue; }
        //-------------------------------------//
        Double_t wgt = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? 1.0 : std::pow(fG4mc->primPart.mom/AXrig.min(), -1.7));
        
        hHCtme->fillH2D(mc_mom, hc_tme);
      
        hCKRrso->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        hHCRrso->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));

        hHCRqltx->fillH2D(mc_mom, trfit.quality(0));
        hHCRqlty->fillH2D(mc_mom, trfit.quality(1));
        
        PhySt&& sttTopPr = prtr.interpolate_to_z(195.0);
        if (prtr.status() && !Numc::EqualToZero(sttTopPr.mom())) {
            hHCRrsoPr->fillH2D(mc_mom, bincen * (sttTopPr.irig() - mc_irig));
            hHCRqltxPr->fillH2D(mc_mom, prtr.quality(0));
            hHCRqltyPr->fillH2D(mc_mom, prtr.quality(1));
        }
        
        Double_t maxx = 0, maxy = 0, maxt = 0, maxr = 0;
        for (auto&& loc : prtr.lscat()) {
            if (loc.scx()) hHCRscatx->fillH2D(loc.chic(0), loc.chis(0));
            if (loc.scy()) hHCRscaty->fillH2D(loc.chic(1), loc.chis(1));
            if (loc.scx()) hHCRsx->fillH2D(mc_mom, loc.chic(0));
            if (loc.scy()) hHCRsy->fillH2D(mc_mom, loc.chic(1));
            if (loc.scx()) hHCRst->fillH2D(mc_mom, loc.chis(0));
            if (loc.scy()) hHCRsr->fillH2D(mc_mom, loc.chis(1));
            maxx = std::max(maxx, loc.chic(0));
            maxy = std::max(maxy, loc.chic(1));
            maxt = std::max(maxt, loc.chis(0));
            maxr = std::max(maxr, loc.chis(1));
        }
        hHCRsxMax->fillH2D(mc_mom, maxx);
        hHCRsyMax->fillH2D(mc_mom, maxy);
        hHCRstMax->fillH2D(mc_mom, maxt);
        hHCRsrMax->fillH2D(mc_mom, maxr);
        
        PhySt&& sttTopD = dtr.interpolate_to_z(195.0);
        if (dtr.status() && !Numc::EqualToZero(sttTopD.mom())) {
            hHCRrsoD->fillH2D(mc_mom, bincen * (sttTopD.irig() - mc_irig));
            hHCRqltxD->fillH2D(mc_mom, dtr.quality(0));
            hHCRqltyD->fillH2D(mc_mom, dtr.quality(1));
        }

        Bool_t TFstatus = (fTof->betaH < 1.0 && ckTr.nchi[0] < 10.0 && ckTr.nchi[1] < 10.0);
        Bool_t RHstatus = (fRich->status && fRich->kind == 0 && fRich->beta < 1.0);
        Bool_t HCstatus = (!mufit.is_like_el());
        
        Double_t TFibta = 1.0 / fTof->betaH;
        Double_t RHibta = 1.0 / fRich->beta;
        if (TFstatus) hTFBrso->fillH2D(mc_mom, (mc_ibta/TFibta - 1.0));
        if (RHstatus) hRHBrso->fillH2D(mc_mom, (mc_ibta/RHibta - 1.0));
        if (HCstatus) hHCTFBrso->fillH2D(mc_mom, (mc_ibta/hc_ibta - 1.0));
        if (HCstatus) hHCRHBrso->fillH2D(mc_mom, (mc_ibta/hc_ibta - 1.0));
        
        PhySt&& sttTopBta = prbta.interpolate_to_z(195.0);
        if (prbta.status() && !Numc::EqualToZero(sttTopBta.mom()) && (sttTopBta.ibta() > prbta.rerr()+1.0)) {
            hHCTFBrsoPr->fillH2D(mc_mom, (mc_ibta/sttTopBta.ibta() - 1.0));
            hHCRHBrsoPr->fillH2D(mc_mom, (mc_ibta/sttTopBta.ibta() - 1.0));
            hHCBqltbPr->fillH2D(mc_mom, prbta.quality());
        }

        Double_t ck_mass = (TFstatus ? ck_mom * std::sqrt(TFibta * TFibta - 1.0) : 0.0);

        if (TFstatus) hCKM->fillH2D(mc_mom, ck_mass);
        if (HCstatus) hHCM->fillH2D(mc_mom, hc_mass);
        
        if (HCstatus) hHCMqltx->fillH2D(mc_mom, mufit.quality(0));
        if (HCstatus) hHCMqlty->fillH2D(mc_mom, mufit.quality(1));
        if (HCstatus) hHCMqltb->fillH2D(mc_mom, mufit.quality(2));
        
        if (TFstatus && HCstatus && mufit.quality(0) < 2.0 && mufit.quality(1) < 2.0 && mufit.quality(2) < 2.0) hCKMcut->fillH2D(mc_mom, ck_mass);
        if (HCstatus && mufit.quality(0) < 2.0 && mufit.quality(1) < 2.0 && mufit.quality(2) < 2.0)             hHCMcut->fillH2D(mc_mom, hc_mass);

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
