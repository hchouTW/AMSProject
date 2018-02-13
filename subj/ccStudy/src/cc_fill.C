//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/ams_home/hchou/AMSCore/prod/18Feb05/src/ClassDef.h"
#include "/ams_home/hchou/AMSCore/prod/18Feb05/src/ClassDef.C"

using namespace std;
using namespace MGROOT;
using namespace TrackSys;

int main(int argc, char * argv[]) {
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
    PartType type = PartType::Proton;
    PhyArg::SetOpt(true, true);
    
    TFile * ofle = new TFile(Form("%s/cc_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 20, 1., 4000., AxisScale::kLog);
    
    Axis AXrig("Rigidity [GV]", 20, 1., 4000., AxisScale::kLog);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);

    Axis AXres("Residual [cm * (30/Momentum+1)^-1]", 12000, -0.4, 0.4);
    Hist * htk[9] = { nullptr };
    Hist * htkh[9] = { nullptr };
    Hist * htkCut[9] = { nullptr };
    Hist * htkhCut[9] = { nullptr };
    Hist * hMCxtk[9] = { nullptr };
    Hist * hMCxtkh[9] = { nullptr };
    Hist * hMCytk[9] = { nullptr };
    Hist * hMCytkh[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il) {
        htk[il] = Hist::New(STR_FMT("htk%d", il+1), HistAxis(AXrig, "Events/Bin"));
        htkh[il] = Hist::New(STR_FMT("htkh%d", il+1), HistAxis(AXrig, AXres, "Events/Bin"));
        htkCut[il] = Hist::New(STR_FMT("htkCut%d", il+1), HistAxis(AXrig, "Events/Bin"));
        htkhCut[il] = Hist::New(STR_FMT("htkhCut%d", il+1), HistAxis(AXrig, AXres, "Events/Bin"));
        hMCxtk[il]  = Hist::New(STR_FMT("hMCxtk%d", il+1), HistAxis(AXmom, "Events/Bin"));
        hMCxtkh[il] = Hist::New(STR_FMT("hMCxtkh%d", il+1), HistAxis(AXmom, AXres, "Events/Bin"));
        hMCytk[il]  = Hist::New(STR_FMT("hMCytk%d", il+1), HistAxis(AXmom, "Events/Bin"));
        hMCytkh[il] = Hist::New(STR_FMT("hMCytkh%d", il+1), HistAxis(AXmom, AXres, "Events/Bin"));
    }
    
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
        if ((fTof->extClsN[0]+fTof->extClsN[1]) > 1 || 
            (fTof->extClsN[2]+fTof->extClsN[3]) > 2) continue; 

        // No Interaction
        Int_t IntType = 0;
        if (fG4mc->primVtx.status) {
            if (fG4mc->primVtx.coo[2] > -55.) IntType = 1;
        }

        // Choutko
        if (!track.status[0][0]) continue;

        Float_t mc_mom = fG4mc->primPart.mom;
        if (mc_mom < AXmom.min()) continue;

        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t bincen  = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
        Double_t wgtpow = std::pow(mc_mom, -1.7);
        
        for (auto&& hit : track.hits) {
            for (auto&& mchit : fG4mc->primPart.hits) {
                if (mchit.layJ != hit.layJ) continue;
            }
            if (hit.side != 3) continue;
            Double_t resx = (track.stateLJ[0][0][hit.layJ-1][0] - hit.coo[0]) / (30./bincen+1.0);
            Double_t resy = (track.stateLJ[0][0][hit.layJ-1][1] - hit.coo[1]) / (30./bincen+1.0);

            if (track.rigidity[0][0] > 0) htk[hit.layJ-1]->fillH1D(track.rigidity[0][0], wgtpow);
            if (track.rigidity[0][0] > 0) htkh[hit.layJ-1]->fillH2D(track.rigidity[0][0], resx, wgtpow);
            if (track.rigidity[0][0] > 0 && track.chisq[0][0][1] < 8) htkCut[hit.layJ-1]->fillH1D(track.rigidity[0][0], wgtpow);
            if (track.rigidity[0][0] > 0 && track.chisq[0][0][1] < 8) htkhCut[hit.layJ-1]->fillH2D(track.rigidity[0][0], resx, wgtpow);
            if (track.chisq[0][0][1] < 8) hMCxtk[hit.layJ-1]->fillH1D(mc_mom);
            if (track.chisq[0][0][1] < 8) hMCxtkh[hit.layJ-1]->fillH2D(mc_mom, resx);
            if (track.chisq[0][0][0] < 8) hMCytk[hit.layJ-1]->fillH1D(mc_mom);
            if (track.chisq[0][0][0] < 8) hMCytkh[hit.layJ-1]->fillH2D(mc_mom, resy);
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
