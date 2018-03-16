#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/18Mar12/src/ClassDef.h"

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();
    
    TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");
    
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
    TFile * ofle = new TFile(Form("%s/hit_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 150, 0.35, 4000., AxisScale::kLog);
    
    Double_t mass = 0.938272297;
    Axis AXeta("1/GammaBeta [1]", AXmom.nbin(), mass/AXmom.max(), mass/AXmom.min(), AxisScale::kLog);

    Double_t lbta = 1.0/std::sqrt(1.0+AXeta.max()*AXeta.max());
    Double_t ubta = 1.0/std::sqrt(1.0+AXeta.min()*AXeta.min());
    Axis AXbta("Beta [1]", AXmom.nbin(), lbta, ubta, AxisScale::kLog);

    // Cut
    Axis AXcut("Cut", 9, 0., 9.);
    Hist* hCut = Hist::New("hCut", HistAxis(AXmom, AXcut));
    Hist* hEvt = Hist::New("hEvt", HistAxis(AXmom, AXcut));
    
    // Coo
    Axis AXres("res [#mum]", 800, -200., 200.);
    Hist* hMrx = Hist::New("hMrx", HistAxis(AXmom, AXres));
    Hist* hMrxNN = Hist::New("hMrxNN", HistAxis(AXres, "Events/Bin"));
    Hist* hMrxN1 = Hist::New("hMrxN1", HistAxis(AXres, "Events/Bin"));
    Hist* hMrxN2 = Hist::New("hMrxN2", HistAxis(AXres, "Events/Bin"));
    Hist* hMrxN3 = Hist::New("hMrxN3", HistAxis(AXres, "Events/Bin"));
    
    Hist* hMry = Hist::New("hMry", HistAxis(AXmom, AXres));
    Hist* hMryNN = Hist::New("hMryNN", HistAxis(AXres, "Events/Bin"));
    Hist* hMryN1 = Hist::New("hMryN1", HistAxis(AXres, "Events/Bin"));
    Hist* hMryN2 = Hist::New("hMryN2", HistAxis(AXres, "Events/Bin"));
    Hist* hMryN3 = Hist::New("hMryN3", HistAxis(AXres, "Events/Bin"));
    Hist* hMryN4 = Hist::New("hMryN4", HistAxis(AXres, "Events/Bin"));
    
    Axis AXedep("Edep", 800, 0., 1.0);
    Hist* hMedep = Hist::New("hMedep", HistAxis(AXeta, AXedep));
    
    Axis AXadc("ADC", 800, 0., 400.);
    Hist* hMadcx = Hist::New("hMadcx", HistAxis(AXeta, AXadc));
    Hist* hMadcy = Hist::New("hMadcy", HistAxis(AXeta, AXadc));
    
    Axis AXTFres("TFres [cm]", 800, -10., 10.);
    Hist* hMTFrx[4] = { nullptr };
    Hist* hMTFry[4] = { nullptr };
    for (int it = 0; it < 4; ++it) {
        hMTFrx[it] = Hist::New(Form("hMTF%drx", it), HistAxis(AXmom, AXTFres));
        hMTFry[it] = Hist::New(Form("hMTF%dry", it), HistAxis(AXmom, AXTFres));
    }
    
    Axis AXTFadc("TFadc", 800, 0., 3.);
    Hist* hMTFadc = Hist::New("hMTFadc", HistAxis(AXeta, AXTFadc));

    Long64_t printRate = static_cast<Long64_t>(0.02*dst->GetEntries());
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        dst->GetEntry(entry);
        
        TrackInfo& track = fTrk->track;

        for (Int_t ic = 0; ic < AXcut.nbin(); ++ic) hEvt->fillH2D(fG4mc->primPart.mom, ic);
        hCut->fillH2D(fG4mc->primPart.mom, 0);
       
        // Geometry (TOF)
        if (fTof->numOfBetaH != 1) continue;
        if (!fTof->statusBetaH) continue;
        if (fTof->betaHPatt != 15) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 1);
        
        // Geometry (TRD)
        if (fTrd->numOfTrack != 1 && fTrd->numOfHTrack != 1) continue;
        if (!fTrd->statusKCls[0]) continue;
        if (fTrd->LLRnhit[0] < 10) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 2);
        
        // Geometry (ACC)
        if (fAcc->clusters.size() != 0) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 3);
        
        // Down-going
        if (fTof->betaH < 0.) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 4);

        // Charge
        //if (fTof->Qall < 0.8 || fTof->Qall > 1.3) continue;
        //if (track.QIn < 0.8 || track.QIn > 1.3) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 5);

        // TOF
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 6);
        
        if ((fTof->numOfExtCls[0]+fTof->numOfExtCls[1]) > 0 || 
            (fTof->numOfExtCls[2]+fTof->numOfExtCls[3]) > 1) continue; 
        hCut->fillH2D(fG4mc->primPart.mom, 7);
        
        if (fTof->numOfInTimeCls > 4) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 8);

        // No Interaction
        if (fG4mc->primVtx.status && fG4mc->primVtx.coo[2] > -100) continue;

        // REC hit
        HitTRKInfo * rec[9]; std::fill_n(rec, 9, nullptr);
        for (auto&& hit : track.hits) { rec[hit.layJ-1] = &hit; }

        // MC hit
        HitTRKMCInfo * mch[9]; std::fill_n(mch, 9, nullptr);
        for (auto&& hit : fG4mc->primPart.hits) { mch[hit.layJ-1] = &hit; }
        
        SegPARTMCInfo * mcs[9]; std::fill_n(mcs, 9, nullptr);
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec==0) mcs[seg.lay-1] = &seg; }
        
        SegPARTMCInfo * mtf[4]; std::fill_n(mtf, 4, nullptr);
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec==1) mtf[seg.lay-1] = &seg; }

        for (Int_t it = 2; it < 8; ++it) {
            if (!rec[it] || !mch[it] || !mcs[it]) continue;
            Double_t eta = mass/mch[it]->mom;
            Double_t res[2] = { rec[it]->coo[0] - mch[it]->coo[0], rec[it]->coo[1] - mch[it]->coo[1] };
            Short_t  ntp[2] = { rec[it]->nsr[0], rec[it]->nsr[1] };
            
            constexpr Double_t CM2UM = 1.0e4;
            if (ntp[0]!=0) hMrx->fillH2D(mch[it]->mom, CM2UM * res[0]);
            if (ntp[1]!=0) hMry->fillH2D(mch[it]->mom, CM2UM * res[1]);
            if (mch[it]->mom > 50.0) {
                if (ntp[0]!=0) hMrxNN->fillH1D(CM2UM * res[0]);
                if (ntp[0]==1) hMrxN1->fillH1D(CM2UM * res[0]);
                if (ntp[0]==2) hMrxN2->fillH1D(CM2UM * res[0]);
                if (ntp[0]>=3) hMrxN3->fillH1D(CM2UM * res[0]);
                
                if (ntp[1]!=0) hMryNN->fillH1D(CM2UM * res[1]);
                if (ntp[1]==1) hMryN1->fillH1D(CM2UM * res[1]);
                if (ntp[1]==2) hMryN2->fillH1D(CM2UM * res[1]);
                if (ntp[1]==3) hMryN3->fillH1D(CM2UM * res[1]);
                if (ntp[1]>=4) hMryN4->fillH1D(CM2UM * res[1]);
            }
            if (ntp[0]!=0 && ntp[1]!=0) hMedep->fillH2D(eta, mch[it]->edep*1.0e3*std::fabs(mcs[it]->dir[2]));
            
            if (ntp[0]!=0 && rec[it]->adc[0]>0) hMadcx->fillH2D(eta, rec[it]->adc[0]);
            if (ntp[1]!=0 && rec[it]->adc[1]>0) hMadcy->fillH2D(eta, rec[it]->adc[1]);
        }
            
        for (Int_t it = 0; it < 4; ++it) {
            if (!mtf[it] || fTof->Q[it]<=0) continue;
            Double_t eta = mass/mtf[it]->mom;
            Double_t dz = fTof->coo[it][2] - mtf[it]->coo[2];
            Double_t tx = mtf[it]->dir[0] / mtf[it]->dir[2];
            Double_t ty = mtf[it]->dir[1] / mtf[it]->dir[2];
            Double_t rx = fTof->coo[it][0] - (mtf[it]->coo[0] + tx * dz);
            Double_t ry = fTof->coo[it][1] - (mtf[it]->coo[1] + ty * dz);
            hMTFrx[it]->fillH2D(mtf[it]->mom, rx);
            hMTFry[it]->fillH2D(mtf[it]->mom, ry);
            hMTFadc->fillH2D(eta, fTof->Q[it]);
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

    return 0;
}
