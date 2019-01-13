#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/19Jan09/src/ClassDef.h"

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();
    
    TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");
    
    //TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/material");

    MGConfig::JobOpt opt(argc, argv);
    PhyArg::SetOpt(true, true);

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
    HYC  * fHyc  = new HYC ;

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
    TFile * ofle = new TFile(Form("%s/hit_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    //PartInfo info(PartType::Proton);
    PartInfo info(PartType::Helium4);
    
    Double_t mombd[2] = { 1., 1000. };
    if (info.type() == PartType::Proton)  { mombd[0] = 0.55; mombd[1] = 3800.0; }
    if (info.type() == PartType::Helium4) { mombd[0] = 2.20; mombd[1] = 45000.0; }
    if (info.type() == PartType::Carbon12) { mombd[0] = 5.00; mombd[1] = 10000.0; }
    Axis AXmom("Momentum [GeV]", 150, mombd[0], mombd[1], AxisScale::kLog);
    
    Axis AXigb("1/GammaBeta [1]", AXmom.nbin(), info.mass()/AXmom.max(), info.mass()/AXmom.min(), AxisScale::kLog);
    Double_t lbta = std::sqrt(1.0+AXigb.min()*AXigb.min());
    Double_t ubta = std::sqrt(1.0+AXigb.max()*AXigb.max());
    Axis AXib("1/Beta [1]", AXigb.nbin(), lbta, ubta, AxisScale::kLinear);

    // Cut
    Axis AXcut("Cut", 11, 0., 11.);
    Hist* hCut = Hist::New("hCut", HistAxis(AXmom, AXcut));
    Hist* hEvt = Hist::New("hEvt", HistAxis(AXmom, AXcut));
  
    // Coo
    Axis AXres("res [#mum]", 1500, -300., 300.);
    
    Hist* hMrxInn   = Hist::New("hMrxInn",   HistAxis(AXigb, AXres));
    Hist* hMrxInnNN = Hist::New("hMrxInnNN", HistAxis(AXres, "Events/Bin"));
  
    Hist* hMryInn   = Hist::New("hMryInn",   HistAxis(AXigb, AXres));
    Hist* hMryInnNN = Hist::New("hMryInnNN", HistAxis(AXres, "Events/Bin"));

    Axis AXTKq("TKq", 1200, 0.6 * info.chrg() * info.chrg(), 10.0 * info.chrg() * info.chrg());
    Hist* hTKqxy = Hist::New("hTKqxy", HistAxis(AXigb, AXTKq));
    
    Axis AXTFq("TFq", 1200, 0.5 * info.chrg() * info.chrg(), 8.0 * info.chrg() * info.chrg());
    Hist* hTFq = Hist::New("hTFq", HistAxis(AXigb, AXTFq));
    
    Axis AXTDn("TDn", 30, 0., 30.);
    Axis AXTDc("TDc", 60, 85., 145.);
    Axis AXTDq("TDq", 400, 0.3 * info.chrg() * info.chrg(), 20.0 * info.chrg() * info.chrg());
    Hist* hTDn  = Hist::New("hTDn", HistAxis(AXigb, AXTDn));
    Hist* hTDc  = Hist::New("hTDc", HistAxis(AXigb, AXTDc));
    Hist* hTDq  = Hist::New("hTDq", HistAxis(AXigb, AXTDq));
    
    Axis AXTFtme("TFtme", 800, -25, 25);
    Hist* hTFtme = Hist::New("hTFtme", HistAxis(AXigb, AXTFtme));
    
    Axis AXAGLib("AGLib", 800, -0.005, 0.005);
    Hist* hAGLib = Hist::New("hAGLib", HistAxis(AXigb, AXAGLib));
    
    Axis AXNAFib("NAFib", 800, -0.015, 0.015);
    Hist* hNAFib = Hist::New("hNAFib", HistAxis(AXigb, AXNAFib));
    
    Axis AXAGLibElem("AGLib", 800*5, -0.005*10, 0.005*10);
    Axis AXNAFibElem("NAFib", 800*5, -0.015*10, 0.015*10);
    Hist* hAGLibElem = Hist::New("hAGLibElem", HistAxis(AXigb, AXAGLibElem));
    Hist* hNAFibElem = Hist::New("hNAFibElem", HistAxis(AXigb, AXNAFibElem));
    
    Axis AXnpe("npe", 100, 0.0, 20.0);
    Hist* hAGLibnpe = Hist::New("hAGLibnpe", HistAxis(AXAGLibElem, AXnpe));
    Hist* hNAFibnpe = Hist::New("hNAFibnpe", HistAxis(AXNAFibElem, AXnpe));

    Long64_t printRate = static_cast<Long64_t>(0.05*dst->GetEntries());
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        dst->GetEntry(entry);
       
        for (Int_t ic = 0; ic < AXcut.nbin(); ++ic) hEvt->fillH2D(fG4mc->primPart.mom, ic);
        hCut->fillH2D(fG4mc->primPart.mom, 0);
       
        // Geometry (TRK)
        if (fTrk->numOfTrack != 1) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 1);
        
        // Geometry (TOF)
        if (fTof->numOfBetaH != 1) continue;
        if (!fTof->statusBetaH) continue;
        if (fTof->betaHPatt != 15) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 2);
        
        // Geometry (TRD)
        //if (fTrd->numOfTrack != 1 && fTrd->numOfHTrack != 1) continue;
        if (!fTrd->statusKCls[0]) continue;
        if (fTrd->LLRnhit[0] < 8) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 3);
        
        // Geometry (ACC)
        if (fAcc->clusters.size() != 0) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 4);
        
        // Down-going
        if (fTof->betaH < 0.) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 5);

        // Charge
        if (std::abs(info.chrg()) == 1) {
            if (fTof->Qall < 0.8 || fTof->Qall > 1.3) continue;
            if (fTrk->QIn < 0.8 || fTrk->QIn > 1.3) continue;
        }
        else {
            if (fTof->Qall < 1.7 || fTof->Qall > 2.4) continue;
            if (fTrk->QIn < 1.7 || fTrk->QIn > 2.4) continue;
        }
        
        hCut->fillH2D(fG4mc->primPart.mom, 6);

        // TOF
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 7);
        
        if (fTof->numOfInTimeCls > 4) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 8);
        
        if ((fTof->numOfExtCls[0]+fTof->numOfExtCls[1]) > 0 || 
            (fTof->numOfExtCls[2]+fTof->numOfExtCls[3]) > 1) continue; 
        hCut->fillH2D(fG4mc->primPart.mom, 9);

        // No Interaction
        if (fG4mc->primVtx.status && fG4mc->primVtx.coo[2] > -120) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 10);

        // REC hit
        HitTRKInfo * rcTk[9]; std::fill_n(rcTk, 9, nullptr);
        for (auto&& hit : fTrk->hits) { rcTk[hit.layJ-1] = &hit; }

        // MC hit
        HitTRKMCInfo* mhTk[9]; std::fill_n(mhTk, 9, nullptr);
        for (auto&& hit : fG4mc->primPart.hitsTk) { mhTk[hit.layJ-1] = &hit; }
        
        HitTOFMCInfo* mhTf[4]; std::fill_n(mhTf, 4, nullptr);
        for (auto&& hit : fG4mc->primPart.hitsTf) { mhTf[hit.lay] = &hit; }
        
        // MC part
        SegPARTMCInfo* mpTk[9]; std::fill_n(mpTk, 9, nullptr);
        for (auto&& seg : fG4mc->primPart.segsTk) { mpTk[seg.lay] = &seg; }
        
        SegPARTMCInfo* mpTf[4]; std::fill_n(mpTf, 4, nullptr);
        for (auto&& seg : fG4mc->primPart.segsTf) { mpTf[seg.lay] = &seg; }
        
        SegPARTMCInfo* mpRh = nullptr;
        for (auto&& seg : fG4mc->primPart.segsRh) { mpRh = &seg; }
        
        for (Int_t it = 2; it <= 7; ++it) {
            if (!rcTk[it] || !mhTk[it] || !mpTk[it]) continue;
            Double_t igb = (fG4mc->primPart.mass / mpTk[it]->mom);
            if (igb < AXigb.min() || igb > AXigb.max()) continue;
            igb = AXigb.center(AXigb.find(igb), AxisScale::kLog);
            Double_t res[2] = { rcTk[it]->coo[0] - mhTk[it]->coo[0], rcTk[it]->coo[1] - mhTk[it]->coo[1] };
            if (!(rcTk[it]->side[0] && rcTk[it]->side[1])) continue;
            if (rcTk[it]->chrg[0] < 0.775) continue;
            if (rcTk[it]->chrg[1] < 0.775) continue;
            if (rcTk[it]->chrg[2] < 0.0) continue;

            constexpr Double_t CM2UM = 1.0e4;
            hMrxInn->fillH2D(igb, CM2UM * res[0]);
            hMryInn->fillH2D(igb, CM2UM * res[1]);
            if (mpTk[it]->mom/fG4mc->primPart.chrg > 40.0) {
               hMrxInnNN->fillH1D(CM2UM * res[0]);
               hMryInnNN->fillH1D(CM2UM * res[1]);
            }
            hTKqxy->fillH2D(igb, rcTk[it]->chrg[2] * rcTk[it]->chrg[2]);
        }
        
        for (Int_t it = 0; it < 4; ++it) {
            if (!mpTf[it] || fTof->Q[it]<=0) continue;
            Double_t igb = info.mass() / mpTf[it]->mom;
            if (igb < AXigb.min() || igb > AXigb.max()) continue;
            hTFq->fillH2D(igb, fTof->Q[it] * fTof->Q[it]);
        }
       
        if (mhTf[0] && mhTf[1] && mhTf[2] && mhTf[3]) {
            for (Int_t sl = 0; sl < 2; ++sl) {
                if (mhTf[sl*2+0]->bta >= 1.0 || mhTf[sl*2+1]->bta >= 1.0) continue;
                TrackSys::SVecD<3> vlen( (mhTf[sl*2+0]->coo[0] - mhTf[sl*2+1]->coo[0]), (mhTf[sl*2+0]->coo[1] - mhTf[sl*2+1]->coo[1]), (mhTf[sl*2+0]->coo[2] - mhTf[sl*2+1]->coo[2]) );
                Double_t rat = (fTof->coo[sl*2+0][2] - fTof->coo[sl*2+1][2]) / (mhTf[sl*2+0]->coo[2] - mhTf[sl*2+1]->coo[2]);
                Double_t len = rat * TrackSys::LA::Mag(vlen);
                
                Double_t ibta = 0.5 * (1.0/mhTf[sl*2+0]->bta + 1.0/mhTf[sl*2+1]->bta);
                Double_t igb = std::sqrt(ibta * ibta - 1.0);

                Double_t tme = len * ibta;
                Double_t mes = 2.99792458e+01 * (fTof->T[sl*2+1] - fTof->T[sl*2+0]);
                Double_t dlt = (mes - tme) * Numc::INV_SQRT_TWO;
                hTFtme->fillH2D(igb, dlt);
            }
        }

        if (mpRh && fRich->status && fRich->isGood) {
            PhySt st(info.type());
            st.set_state_with_cos(mpRh->coo[0], mpRh->coo[1], mpRh->coo[2], mpRh->dir[0], mpRh->dir[1], mpRh->dir[2]);
            st.set_mom(mpRh->mom);
            TrackSys::PropMgnt::PropToZ(fRich->refz, st);
            Double_t dlt = (1.0/fRich->beta - 1.0/st.bta());
            Double_t igb = st.igb();
            if (fRich->kind == 0) hAGLib->fillH2D(igb, dlt);
            if (fRich->kind == 1) hNAFib->fillH2D(igb, dlt);

            for (auto&& hit : fRich->uhits) {
                Double_t dlt = (1.0/hit.bta - 1.0/st.bta());
                if (fRich->kind == 0) hAGLibElem->fillH2D(igb, dlt);
                if (fRich->kind == 1) hNAFibElem->fillH2D(igb, dlt);
                if (fRich->kind == 0 && igb < 0.1) hAGLibnpe->fillH2D(dlt, hit.npe);
                if (fRich->kind == 1 && igb < 0.1) hNAFibnpe->fillH2D(dlt, hit.npe);
            }
        }
        
        // TRD
        if (fTrd->ITstatus[0]) {
            hTDn->fill(info.mass()/fTrd->ITMcMom[0], fTrd->ITnh[0]);
            hTDc->fill(info.mass()/fTrd->ITMcMom[0], fTrd->ITcz[0]);
            hTDq->fill(info.mass()/fTrd->ITMcMom[0], fTrd->ITdEdX[0]);
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
    if (fHyc)  { delete fHyc ; fHyc  = nullptr; }

    return 0;
}
