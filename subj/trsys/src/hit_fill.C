//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/ams_home/hchou/AMSCore/prod/18Feb27/src/ClassDef.h"

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

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
    
    Axis AXmom("Momentum [GeV]", 100, 0.35, 4000., AxisScale::kLog);
    
    Double_t mass = 0.938272297;
    Axis AXeta("1/GammaBeta [1]", 100, mass/4000., mass/0.35, AxisScale::kLog);

    // MC
    Axis AXcut("Cut", 8, 0., 8.);
    Hist* hCut = Hist::New("hCut", HistAxis(AXmom, AXcut));
    Hist* hEvt = Hist::New("hEvt", HistAxis(AXmom, AXcut));

    // Coo
    Axis AXres("res [#mum]", 800, -200., 200.);
    Hist* hMrx = Hist::New("hMrx", HistAxis(AXmom, AXres));
    Hist* hMry = Hist::New("hMry", HistAxis(AXmom, AXres));
    
    Axis AXadc("ADC", 1200, 0., 400.);
    Hist* hMax = Hist::New("hMax", HistAxis(AXeta, AXadc));
    Hist* hMay = Hist::New("hMay", HistAxis(AXeta, AXadc));
    
    Axis AXedep("Edep", 1600, 0., 1.0);
    Hist* hMedep = Hist::New("hMedep", HistAxis(AXeta, AXedep));
    
    Axis AXloc("loc [1]", 40, -0.5, 0.5);
    Hist* hMrxNN = Hist::New("hMrxNN", HistAxis(AXres, "Events/Bin"));
    Hist* hMrxN1 = Hist::New("hMrxN1", HistAxis(AXres, "Events/Bin"));
    Hist* hMrxN2 = Hist::New("hMrxN2", HistAxis(AXres, "Events/Bin"));
    Hist* hMrxN3 = Hist::New("hMrxN3", HistAxis(AXres, "Events/Bin"));
    
    Hist* hMryNN = Hist::New("hMryNN", HistAxis(AXres, "Events/Bin"));
    Hist* hMryN1 = Hist::New("hMryN1", HistAxis(AXres, "Events/Bin"));
    Hist* hMryN2 = Hist::New("hMryN2", HistAxis(AXres, "Events/Bin"));
    Hist* hMryN3 = Hist::New("hMryN3", HistAxis(AXres, "Events/Bin"));
    Hist* hMryN4 = Hist::New("hMryN4", HistAxis(AXres, "Events/Bin"));

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
        if (fTof->Qall < 0.8 || fTof->Qall > 1.3) continue;
        if (track.QIn < 0.8 || track.QIn > 1.3) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 5);

        // TOF
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        hCut->fillH2D(fG4mc->primPart.mom, 6);
        
        if (fTof->numOfInTimeCls > 4) continue;
        if ((fTof->extClsN[0]+fTof->extClsN[1]) > 0 || 
            (fTof->extClsN[2]+fTof->extClsN[3]) > 1) continue; 
        hCut->fillH2D(fG4mc->primPart.mom, 7);

        // No Interaction
        if (fG4mc->primVtx.status && fG4mc->primVtx.coo[2] > -100) continue;

        // REC hit
        HitTRKInfo * rec[9]; std::fill_n(rec, 9, nullptr);
        for (auto&& hit : track.hits) { rec[hit.layJ-1] = &hit; }

        // MC hit
        HitTRKMCInfo * mch[9]; std::fill_n(mch, 9, nullptr);
        for (auto&& hit : fG4mc->primPart.hits) { mch[hit.layJ-1] = &hit; }
        
        SegPARTMCInfo * mcs[9]; std::fill_n(mcs, 9, nullptr);
        for (auto&& seg : fG4mc->primPart.segs) { mcs[seg.lay-1] = &seg; }

        for (Int_t it = 2; it < 8; ++it) {
            if (!rec[it] || !mch[it] || !mcs[it]) continue;
            Short_t  ntp[2] = { rec[it]->nsr[0], rec[it]->nsr[1] };
            Double_t loc[2] = { rec[it]->loc[0] - std::lrint(rec[it]->loc[0]), rec[it]->loc[1] - std::lrint(rec[it]->loc[1]) };
            Double_t res[2] = { rec[it]->coo[0] - mch[it]->coo[0], rec[it]->coo[1] - mch[it]->coo[1] };
            Double_t adc[2] = { rec[it]->sig[0], rec[it]->sig[1] };
            
            Double_t dsdz = std::fabs(mcs[it]->dir[2]);

            constexpr Double_t CM2UM = 1.0e4;
            if (ntp[0]!=0) hMrx->fillH2D(mch[it]->mom, CM2UM * res[0]);
            if (ntp[1]!=0) hMry->fillH2D(mch[it]->mom, CM2UM * res[1]);
            
            if (ntp[0]!=0) hMax->fillH2D(mass/mch[it]->mom, adc[0]*dsdz);
            if (ntp[1]!=0) hMay->fillH2D(mass/mch[it]->mom, adc[1]*dsdz);
            
            if (ntp[0]!=0 && ntp[1]!=0) hMedep->fillH2D(mass/mch[it]->mom, mch[it]->edep*1.0e3*dsdz);

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
