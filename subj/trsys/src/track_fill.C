#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

//#include "/ams_home/hchou/AMSCore/prod/18Sep21/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18Oct17/src/ClassDef.h"
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
    TFile * ofle = new TFile(Form("%s/track_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    //PartInfo info(PartType::Electron);
    PartInfo info(PartType::Proton);
    //PartInfo info(PartType::Helium4);
    
    PartInfo::SetDefault(info.type());
    PhyArg::SetOpt(true, true);
    Bool_t optL1 = true;
    Bool_t optL9 = true;
    
    Int_t    nmom     = 80;
    Double_t mombd[2] = { 1., 1000. };
    if (info.type() == PartType::Electron) { mombd[0] = 0.30; mombd[1] = 450.0; }
    if (info.type() == PartType::Proton)   { mombd[0] = 0.55; mombd[1] = 3800.0; }
    if (info.type() == PartType::Helium4)  { mombd[0] = 2.20; mombd[1] = 3800.0; }
    Axis AXmom("Momentum [GeV]", nmom, mombd[0], mombd[1], AxisScale::kLog);
    
    Axis AXrig("Rigidity [GV]", AXmom.nbin(), mombd[0]/std::fabs(info.chrg()), mombd[1]/std::fabs(info.chrg()), AxisScale::kLog);
    
    // Time
    Axis AXtme("Time [ms]", 1600, 0., 1000.);
    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXmom, AXtme));
    
    // Fit Eff
    Hist* hCKnum = Hist::New("hCKnum", HistAxis(AXmom, "Events/Bin"));
    Hist* hHCnum = Hist::New("hHCnum", HistAxis(AXmom, "Events/Bin"));

    // Fit R Res
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 2000, -1.3, 1.3);
    Hist* hCKRrso = Hist::New("hCKRrso", HistAxis(AXmom, AXRrso));
    Hist* hHCRrso = Hist::New("hHCRrso", HistAxis(AXmom, AXRrso));
    
    //Axis AXRDrso("(1/R_L9 - 1/R_L1) [1/GV]", 2000, -1.8, 1.8);
    //Hist* hCKRDrso = Hist::New("hCKRDrso", HistAxis(AXmom, AXRDrso));
    //Hist* hHCRDrso = Hist::New("hHCRDrso", HistAxis(AXmom, AXRDrso));
    //
    //Hist* hCKRD2rso = Hist::New("hCKRD2rso", HistAxis(AXmom, AXRDrso));
    //Hist* hHCRD2rso = Hist::New("hHCRD2rso", HistAxis(AXmom, AXRDrso));
    
    Axis AXRqlt("Quality [1]", 800, -2.0, 4.0);
    Hist* hCKRqltx = Hist::New("hCKRqltx", HistAxis(AXmom, AXRqlt));
    Hist* hHCRqltx = Hist::New("hHCRqltx", HistAxis(AXmom, AXRqlt));
    
    Hist* hCKRqlty = Hist::New("hCKRqlty", HistAxis(AXmom, AXRqlt));
    Hist* hHCRqlty = Hist::New("hHCRqlty", HistAxis(AXmom, AXRqlt));
   
    MGClock::HrsStopwatch hrssw;
    Long64_t passEntry = 0;
    Long64_t printRate = static_cast<Long64_t>(0.1 * dst->GetEntries());
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) { // testcode
        if (entry%printRate==0) {
            hrssw.stop();
            COUT("Entry %lld/%lld (%lld) Time %14.8f\n", entry, dst->GetEntries(), passEntry, hrssw.time());
        }
        dst->GetEntry(entry);
     
        // No Interaction (testcode)
        //if (opt.mode() == MGConfig::JobOpt::MODE::MC)
        //    if (fG4mc->primVtx.status && fG4mc->primVtx.coo[2] > -120) continue;
        
        //if (entry > 100) break;
        
        Int_t trPatt = optL1 + optL9 * 2;
        CKTrackInfo& ckTr = fTrk->ckTr.at(trPatt);
        KFTrackInfo& kfTr = fTrk->kfTr.at(trPatt);
        //HCTrackInfo& hcTr = fTrk->hcTr.at(trPatt);
        
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
       
        Int_t topLay = 1000;
        Bool_t hasL1 = false;
        Bool_t hasL9 = false;
        TrFitPar fitPar(info.type());
        TrFitPar fitParL1(info.type());
        TrFitPar fitParL9(info.type());
        for (auto&& hit : fTrk->hits) {
            Bool_t isInnTr = (hit.layJ >= 2 && hit.layJ <= 8);
            HitStTRK mhit(hit.side[0], hit.side[1], hit.layJ);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
            //mhit.set_q(hit.chrg[2], hit.chrg[0], hit.chrg[1]);
         
            if (isInnTr) { fitPar.add_hit(mhit); fitParL1.add_hit(mhit); fitParL9.add_hit(mhit); topLay = std::min(topLay, hit.layJ-1); }
            else {
                if (optL1 && hit.layJ == 1) { hasL1 = true; fitPar.add_hit(mhit); fitParL1.add_hit(mhit); topLay = 0; }
                if (optL9 && hit.layJ == 9) { hasL9 = true; fitPar.add_hit(mhit); fitParL9.add_hit(mhit); }
            }
        }

        for (Int_t il = 0; il < 4; ++il) {
            HitStTOF mhit(il);
            mhit.set_coo(fTof->coo[il][0], fTof->coo[il][1], fTof->coo[il][2]);
            mhit.set_q(fTof->Q[il]);
            mhit.set_t(fTof->T[il]*HitStTOF::TRANS_NS_TO_CM);
            //fitPar.add_hit(mhit);
        }

        //if (!fRich->status) continue;
        //if (fRich->kind != 0) continue;
        //HitStRICH richHit( (fRich->kind == 0 ? HitStRICH::Radiator::AGL : HitStRICH::Radiator::NAF) );
        //richHit.set_coo(Numc::ZERO<>, Numc::ZERO<>, fRich->refz);
        //richHit.set_ib(Numc::ONE<> / fRich->beta);
        //btaPar.add_hit(richHit);

        // TRD
        if (fTrd->statusKCls[0] && fTrd->recHits.size() >= 5) {
            HitStTRD mhit;
            mhit.set_coo(0, 0, fTrd->recCz);
            mhit.set_el(fTrd->recMen);
            fitPar.add_hit(mhit);
        }
        else continue;

        if (!fitPar.check()) continue;
        //if (!fitParL1.check()) continue;
        //if (!fitParL9.check()) continue;

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

        SegPARTMCInfo* topmc = nullptr;
        for (auto&& seg : fG4mc->primPart.segsTk) { if (seg.lay == topLay) { topmc = &seg; break; } }
        if (topmc == nullptr) continue;

        Double_t mc_mom  = fG4mc->primPart.mom;
        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t bincen  = std::sqrt(AXmom.center(AXmom.find(mc_mom), AxisScale::kLog));
        
        //-------------------------------------//
        MGClock::HrsStopwatch sw; sw.start();
        PhyTrFit tr(fitPar);
        sw.stop();
        
        //if (tr.status()) COUT("MC INFO   MOM %14.8f   IGB %14.8f    RIG %14.8f\n", fG4mc->primPart.mom, info.mass()/fG4mc->primPart.mom, tr.part().rig()); // testcode

        Bool_t   hc_succ = tr.status();
        Double_t hc_irig = tr.part().irig();
        Double_t hc_tme  = sw.time()*1.0e3;
  
        if (!hc_succ) COUT("HC PHY FAILURE.\n");
        //if (!trbta.status()) COUT("HC BTA FAILURE.\n");

        PhySt&& sttTop = tr.interpolate_to_z(195.0);
        hc_succ = (hc_succ ? !Numc::EqualToZero(sttTop.mom()) : false);
        if (hc_succ) hc_irig = sttTop.irig();
        
        //PhyTrFit trL1(fitParL1);
        //PhyTrFit trL9(fitParL9);
        //if (!trL1.status() || !trL9.status() || !tr.status()) continue;
        //PhySt&& sttTopL1 = trL1.interpolate_to_z(195.0);
        //PhySt&& sttTopL9 = trL9.interpolate_to_z(195.0);
        //-------------------------------------//
        
        Bool_t ck_succ = ckTr.status;
        if (ck_succ) hCKnum->fillH1D(mc_mom);
        if (hc_succ) hHCnum->fillH1D(mc_mom);
        
        if (hc_succ) hHCtme->fillH2D(mc_mom, hc_tme);

        Double_t ck_irig = (ck_succ ? MGMath::ONE/ckTr.rig : 0.);
      
        Double_t ck_qltx = (ck_succ ? std::log(ckTr.nchi[0]) : 0.); 
        Double_t hc_qltx = (hc_succ ? tr.quality(0) : 0.); 
        
        Double_t ck_qlty = (ck_succ ? std::log(ckTr.nchi[1]) : 0.); 
        Double_t hc_qlty = (hc_succ ? tr.quality(1) : 0.); 
        
        if (ck_succ) hCKRrso->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (hc_succ) hHCRrso->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_succ) hCKRqltx->fillH2D(mc_mom, ck_qltx);
        if (hc_succ) hHCRqltx->fillH2D(mc_mom, hc_qltx);
        
        if (ck_succ) hCKRqlty->fillH2D(mc_mom, ck_qlty);
        if (hc_succ) hHCRqlty->fillH2D(mc_mom, hc_qlty);
        
        //// Track L1
        //CKTrackInfo& ckTrL1 = fTrk->ckTr.at(1);
        //HCTrackInfo& hcTrL1 = fTrk->hcTr.at(1);
        //
        //// Track L9
        //CKTrackInfo& ckTrL9 = fTrk->ckTr.at(2);
        //HCTrackInfo& hcTrL9 = fTrk->hcTr.at(2);
        //
        //Double_t wgt = ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? 1.0 : std::pow(fG4mc->primPart.mom/AXrig.min(), -1.7));
        //
        //if (ckTrL1.status & ckTrL9.status && ckTr.status) {
        //    double ckrd = (1.0/ckTrL9.rig - 1.0/ckTrL1.rig);
        //    hCKRDrso->fillH2D(mc_mom, bincen * ckrd);
        //    
        //    Double_t ckcen = std::sqrt(AXmom.center(AXmom.find(std::fabs(2*ckTr.rig)), AxisScale::kLog));
        //    hCKRD2rso->fillH2D(2*ckTr.rig, ckcen * ckrd, wgt);
        //}
        //double hcrd = (1.0/sttTopL9.rig() - 1.0/sttTopL1.rig());
        //hHCRDrso->fillH2D(mc_mom, bincen * hcrd);
        //
        //Double_t hccen = std::sqrt(AXmom.center(AXmom.find(std::fabs(2*sttTop.rig())), AxisScale::kLog));
        //hHCRD2rso->fillH2D(2*sttTop.rig(), hccen * hcrd, wgt);
        //
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
