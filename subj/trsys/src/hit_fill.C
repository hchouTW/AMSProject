#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

//#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/18Mar23/src/ClassDef.h"
#include "/ams_home/hchou/AMSCore/prod/18Mar23/src/ClassDef.h"

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
    TFile * ofle = new TFile(Form("%s/hit_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 100, 0.5, 2000., AxisScale::kLog);
    
    Double_t mass = PartInfo(PartType::Proton).mass();
    Axis AXeta("1/GammaBeta [1]", AXmom.nbin(), mass/AXmom.max(), mass/AXmom.min(), AxisScale::kLog);

    Double_t lbta = 1.0/std::sqrt(1.0+AXeta.max()*AXeta.max());
    Double_t ubta = 1.0/std::sqrt(1.0+AXeta.min()*AXeta.min());
    Axis AXbta("Beta [1]", AXmom.nbin(), lbta, ubta, AxisScale::kLog);

    // Cut
    Axis AXcut("Cut", 9, 0., 9.);
    Hist* hCut = Hist::New("hCut", HistAxis(AXmom, AXcut));
    Hist* hEvt = Hist::New("hEvt", HistAxis(AXmom, AXcut));
/*  
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
    
    Axis AXTKadc("TKadc", 800, 0., 8.);
    Hist* hTKadcx = Hist::New("hTKadcx", HistAxis(AXeta, AXTKadc));
    Hist* hTKadcy = Hist::New("hTKadcy", HistAxis(AXeta, AXTKadc));
    
    Axis AXTFadc("TFadc", 800, 0., 4.);
    Hist* hTFadc = Hist::New("hTFadc", HistAxis(AXeta, AXTFadc));
    
    Axis AXTDadc("TDadc", 800, 0., 4000.);
    Hist* hTDadc = Hist::New("hTDadc", HistAxis(AXeta, AXTDadc));
    
    Axis AXTDavg("TDavg", 400, 0., 1500.);
    Hist* hTDavg = Hist::New("hTDavg", HistAxis(AXeta, AXTDavg));
*/
    Axis AXTFtme("TFtme", 800, -20, 20);
    Hist* hTFtme = Hist::New("hTFtme", HistAxis(AXeta, AXTFtme));

    Long64_t printRate = static_cast<Long64_t>(0.05*dst->GetEntries());
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
/*
        // REC hit
        HitTRKInfo * rec[9]; std::fill_n(rec, 9, nullptr);
        for (auto&& hit : track.hits) { rec[hit.layJ-1] = &hit; }

        // MC hit
        HitTRKMCInfo * mch[9]; std::fill_n(mch, 9, nullptr);
        for (auto&& hit : fG4mc->primPart.hits) { mch[hit.layJ-1] = &hit; }
        
        SegPARTMCInfo * mcs[9]; std::fill_n(mcs, 9, nullptr);
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec==0) mcs[seg.lay] = &seg; }
        
        SegPARTMCInfo * mtf[4]; std::fill_n(mtf, 4, nullptr);
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec==1) mtf[seg.lay] = &seg; }
    
        for (Int_t it = 2; it < 8; ++it) {
            if (!rec[it] || !mch[it] || !mcs[it]) continue;
            Double_t eta = (mass/mch[it]->mom);
            if (eta < AXeta.min() || eta > AXeta.max()) continue;
            eta = AXeta.center(AXeta.find(eta), AxisScale::kLog);
            Double_t bta = 1.0/std::sqrt(1.0+eta*eta);
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
            if (ntp[0]!=0 && rec[it]->adc[0]>0) hTKadcx->fillH2D(eta, rec[it]->adc[0]*rec[it]->adc[0]*bta*bta);
            if (ntp[1]!=0 && rec[it]->adc[1]>0) hTKadcy->fillH2D(eta, rec[it]->adc[1]*rec[it]->adc[1]*bta*bta);
        }
            
        for (Int_t it = 0; it < 4; ++it) {
            if (!mtf[it] || fTof->Q[it]<=0) continue;
            Double_t eta = std::sqrt(1.0/fTof->mcBeta[it]/fTof->mcBeta[it]-1);
            if (eta < AXeta.min() || eta > AXeta.max()) continue;
            eta = AXeta.center(AXeta.find(eta), AxisScale::kLog);
            Double_t bta = 1.0/std::sqrt(1.0+eta*eta);
            hTFadc->fillH2D(eta, fTof->Q[it]*fTof->Q[it]*bta*bta);
        }

        std::vector<Double_t> vdedx;
        for (auto&& hit : fTrd->hits[0]) {
            if (hit.len <= 0 || hit.amp <= 0) continue;
            if (hit.len < 0.3) continue;
            Double_t dedx = hit.amp/hit.len;
            hTDadc->fillH2D(mass/fG4mc->primPart.mom, dedx);
            vdedx.push_back(dedx);
        }
        if (vdedx.size() >= 5) {
            std::sort(vdedx.begin(), vdedx.end());
            Double_t avg = std::accumulate(vdedx.begin()+1, vdedx.end()-1, 0) / (vdedx.size()-2);
            hTDavg->fillH2D(mass/fG4mc->primPart.mom, avg);
        }
*/
        SegPARTMCInfo* mcsTOF[4] = { nullptr };
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec == 1) mcsTOF[seg.lay] = &seg; }
        for (Int_t sl = 0; sl < 2; ++sl) {
            if (mcsTOF[sl*2+0] ==nullptr || mcsTOF[sl*2+1] == nullptr) continue;
            if (fTof->mcBeta[sl*2+0] < 0.001 || fTof->mcBeta[sl*2+1] < 0.001) continue;
            TrackSys::SVecD<3> vlen( (mcsTOF[sl*2+0]->coo[0] - mcsTOF[sl*2+1]->coo[0]), (mcsTOF[sl*2+0]->coo[1] - mcsTOF[sl*2+1]->coo[1]), (mcsTOF[sl*2+0]->coo[2] - mcsTOF[sl*2+1]->coo[2]) );
            Double_t rat = (fTof->coo[sl*2+0][2] - fTof->coo[sl*2+1][2]) / (mcsTOF[sl*2+0]->coo[2] - mcsTOF[sl*2+1]->coo[2]);
            Double_t len = rat * TrackSys::LA::Mag(vlen);
            Double_t tme = 0.5 * len * (1.0/fTof->mcBeta[sl*2+0] + 1.0/fTof->mcBeta[sl*2+1]);
            Double_t mes = 2.99792458e+01 * (fTof->T[sl*2+1] - fTof->T[sl*2+0]);
            Double_t dlt = mes - tme;
            hTFtme->fillH2D(mass/(mcsTOF[sl*2+0]->mom + mcsTOF[sl*2+1]->mom), dlt);
            //CERR("TME %14.8f MES %14.8f DLT %14.8f\n", tme, mes, dlt);
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
