#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/18Mar23/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18Mar23/src/ClassDef.h"

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();
    
    TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    
    //TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");

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
    //PhyArg::SetOpt(true, false);
    PhyArg::SetOpt(true, true);
    
    TFile * ofle = new TFile(Form("%s/time_fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 200, 0.5, 4000., AxisScale::kLog);
    
    Double_t mass = 0.938272297;
    Axis AXeta("1/GammaBeta [1]", AXmom.nbin(), mass/AXmom.max(), mass/AXmom.min(), AxisScale::kLog);
    
    Double_t lbta = 1.0/std::sqrt(1.0+AXeta.max()*AXeta.max());
    Double_t ubta = 1.0/std::sqrt(1.0+AXeta.min()*AXeta.min());
    Axis AXbta("Beta [1]", AXmom.nbin(), lbta, ubta, AxisScale::kLog);
    
    // Prop
    Axis AXtme("Time [cm]", 800, -30, 30);
   
    Hist* hTme = Hist::New("hTme", HistAxis(AXeta, AXtme)); // (TH2D) 

    Long64_t printRate = static_cast<Long64_t>(0.04*dst->GetEntries());
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
        if ((fTof->numOfExtCls[0]+fTof->numOfExtCls[1]) > 0 || 
            (fTof->numOfExtCls[2]+fTof->numOfExtCls[3]) > 1) continue; 
        
        // No Interaction
        if (fG4mc->primVtx.status && fG4mc->primVtx.coo[2] > -100.) continue;

        // MC hit
        SegPARTMCInfo* mcs[4]; std::fill_n(mcs, 4, nullptr);
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec == 1) mcs[seg.lay] = &seg; }

        Bool_t check = true;
        for (Int_t il = 0; il < 4; ++il)
            if (!mcs[il]) { check = false; break; }
        if (!check) continue;

        for (Int_t sl = 0; sl < 2; ++sl) {
            SegPARTMCInfo* mcsU = mcs[2*sl+0];
            SegPARTMCInfo* mcsL = mcs[2*sl+1];
            
            PhySt part(PartType::Proton);
            part.set_state_with_cos(
                mcsU->coo[0], mcsU->coo[1], mcsU->coo[2],
                mcsU->dir[0], mcsU->dir[1], mcsU->dir[2]
            );
            part.set_mom(mcsU->mom);
            
            Double_t mc_eta = std::sqrt(1.0/fTof->mcBeta[2*sl+0]/fTof->mcBeta[2*sl+0] - 1);
            
            PhySt ppst(part); // Particle Status
            PropMgnt::PropToZ(fTof->coo[2*sl+0][2], ppst); // Propagate to Z with magnetic field
            if (!MGNumc::Valid(mc_eta) || MGNumc::EqualToZero(mc_eta)) mc_eta = mass / mcsU->mom;
            part.set_eta(mc_eta);
            
            Double_t tmeU = ppst.time();
            PropMgnt::PropToZ(fTof->coo[2*sl+1][2], ppst); // Propagate to Z with magnetic field
            Double_t tmeL = ppst.time();

            Double_t genTme = (tmeL - tmeU);
            Double_t recTme = (fTof->T[2*sl+1] - fTof->T[2*sl+0]) * HitStTOF::TRANS_NS_TO_CM;
            
            mc_eta = std::sqrt(1.0/fTof->mcBeta[2*sl+0]/fTof->mcBeta[2*sl+1] - 1.0);

            //COUT("L%d TME %14.8f %14.8f DIF %14.8f\n", il, genTme, recTme, recTme-genTme);
            hTme->fillH2D(mc_eta, recTme-genTme);
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
