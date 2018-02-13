//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

//#include "/ams_home/hchou/AMSCore/prod/18Jan29/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18Jan29/src/ClassDef.C"

#include "/ams_home/hchou/AMSCore/prod/18Feb05/src/ClassDef.h"
#include "/ams_home/hchou/AMSCore/prod/18Feb05/src/ClassDef.C"

//#include "/ams_home/hchou/AMSCore/prod/18Jan29/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/18Jan29/src/ClassDef.C"

//#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/17Dec23/src/ClassDef.h"
//#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/17Dec23/src/ClassDef.C"

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
    G4MC * fG4mc = (opt.type() == "MC" ) ? new G4MC : nullptr;
    RTI  * fRti  = (opt.type() == "ISS") ? new RTI  : nullptr;
    TRG  * fTrg  = new TRG ;
    TOF  * fTof  = new TOF ;
    ACC  * fAcc  = new ACC ;
    TRK  * fTrk  = new TRK ;
    TRD  * fTrd  = new TRD ;
    RICH * fRich = new RICH;
    ECAL * fEcal = new ECAL;

    dst->SetBranchAddress("list", &fList);
    if (opt.type() == "MC")
        dst->SetBranchAddress("g4mc", &fG4mc);
    if (opt.type() == "ISS")
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
    //PartType type = PartType::Electron;
    PhyArg::SetOpt(true, true);
    Bool_t optL1 = false;
    Bool_t optL9 = false;
    
    TFile * ofle = new TFile(Form("%s/fit_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 100, 0.5, 4000., AxisScale::kLog);
    
    Axis AXrig("Rigidity [GV]", 100, 0.5, 4000., AxisScale::kLog);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);

    // Time
    Axis AXtme("Time [ms]", 400, 0., 2.);
    Hist * hCKtme = Hist::New("hCKtme", HistAxis(AXmom, AXtme));
    Hist * hCNtme = Hist::New("hCNtme", HistAxis(AXmom, AXtme));
    Hist * hKFtme = Hist::New("hKFtme", HistAxis(AXmom, AXtme));
    Hist * hHCtme = Hist::New("hHCtme", HistAxis(AXmom, AXtme));
    
    // Fit Eff
    Hist * hevt = Hist::New("hevt", HistAxis(AXmom, "Events/Bin"));
    Hist * hCKnum = Hist::New("hCKnum", HistAxis(AXmom, "Events/Bin"));
    Hist * hCNnum = Hist::New("hCNnum", HistAxis(AXmom, "Events/Bin"));
    Hist * hKFnum = Hist::New("hKFnum", HistAxis(AXmom, "Events/Bin"));
    Hist * hHCnum = Hist::New("hHCnum", HistAxis(AXmom, "Events/Bin"));

    // Fit Res
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 2000, -10.0, 10.0);
    Hist * hCKRrso = Hist::New("hCKRrso", HistAxis(AXmom, AXRrso));
    Hist * hCNRrso = Hist::New("hCNRrso", HistAxis(AXmom, AXRrso));
    Hist * hKFRrso = Hist::New("hKFRrso", HistAxis(AXmom, AXRrso));
    Hist * hHCRrso = Hist::New("hHCRrso", HistAxis(AXmom, AXRrso));
    
    Hist * hCKRrsoI0 = Hist::New("hCKRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI0 = Hist::New("hCNRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hKFRrsoI0 = Hist::New("hKFRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI0 = Hist::New("hHCRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI1 = Hist::New("hCKRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI1 = Hist::New("hCNRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hKFRrsoI1 = Hist::New("hKFRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI1 = Hist::New("hHCRrsoI1", HistAxis(AXmom, AXRrso));
    
    Axis AXRchi("Log-Chi-square [1]", 800, -8.0, 8.0);
    Hist * hCKRchix = Hist::New("hCKRchix", HistAxis(AXmom, AXRchi));
    Hist * hCNRchix = Hist::New("hCNRchix", HistAxis(AXmom, AXRchi));
    Hist * hKFRchix = Hist::New("hKFRchix", HistAxis(AXmom, AXRchi));
    Hist * hHCRchix = Hist::New("hHCRchix", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRchixI0 = Hist::New("hCKRchixI0", HistAxis(AXmom, AXRchi));
    Hist * hCNRchixI0 = Hist::New("hCNRchixI0", HistAxis(AXmom, AXRchi));
    Hist * hKFRchixI0 = Hist::New("hKFRchixI0", HistAxis(AXmom, AXRchi));
    Hist * hHCRchixI0 = Hist::New("hHCRchixI0", HistAxis(AXmom, AXRchi));
    Hist * hCKRchixI1 = Hist::New("hCKRchixI1", HistAxis(AXmom, AXRchi));
    Hist * hCNRchixI1 = Hist::New("hCNRchixI1", HistAxis(AXmom, AXRchi));
    Hist * hKFRchixI1 = Hist::New("hKFRchixI1", HistAxis(AXmom, AXRchi));
    Hist * hHCRchixI1 = Hist::New("hHCRchixI1", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRchiy = Hist::New("hCKRchiy", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiy = Hist::New("hCNRchiy", HistAxis(AXmom, AXRchi));
    Hist * hKFRchiy = Hist::New("hKFRchiy", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiy = Hist::New("hHCRchiy", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRchiyI0 = Hist::New("hCKRchiyI0", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiyI0 = Hist::New("hCNRchiyI0", HistAxis(AXmom, AXRchi));
    Hist * hKFRchiyI0 = Hist::New("hKFRchiyI0", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiyI0 = Hist::New("hHCRchiyI0", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiyI1 = Hist::New("hCKRchiyI1", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiyI1 = Hist::New("hCNRchiyI1", HistAxis(AXmom, AXRchi));
    Hist * hKFRchiyI1 = Hist::New("hKFRchiyI1", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiyI1 = Hist::New("hHCRchiyI1", HistAxis(AXmom, AXRchi));
    
    Hist * hCKMCflux = Hist::New("hCKMCflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hCNMCflux = Hist::New("hCNMCflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hKFMCflux = Hist::New("hKFMCflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hHCMCflux = Hist::New("hHCMCflux", HistAxis(AXrig, "Events/Bin"));
    
    Hist * hCKIRflux = Hist::New("hCKIRflux", HistAxis(AXirig, AXRchi));
    Hist * hCNIRflux = Hist::New("hCNIRflux", HistAxis(AXirig, AXRchi));
    Hist * hKFIRflux = Hist::New("hKFIRflux", HistAxis(AXirig, AXRchi));
    Hist * hHCIRflux = Hist::New("hHCIRflux", HistAxis(AXirig, AXRchi));
    
    Hist * hCKRflux = Hist::New("hCKRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hCNRflux = Hist::New("hCNRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hKFRflux = Hist::New("hKFRflux", HistAxis(AXrig, "Events/Bin"));
    Hist * hHCRflux = Hist::New("hHCRflux", HistAxis(AXrig, "Events/Bin"));
    
    Hist * hCKMCfluxCut = Hist::New("hCKMCfluxCut", HistAxis(AXrig, "Events/Bin"));
    Hist * hCNMCfluxCut = Hist::New("hCNMCfluxCut", HistAxis(AXrig, "Events/Bin"));
    Hist * hKFMCfluxCut = Hist::New("hKFMCfluxCut", HistAxis(AXrig, "Events/Bin"));
    Hist * hHCMCfluxCut = Hist::New("hHCMCfluxCut", HistAxis(AXrig, "Events/Bin"));
    
    Hist * hCKRfluxCut = Hist::New("hCKRfluxCut", HistAxis(AXrig, "Events/Bin"));
    Hist * hCNRfluxCut = Hist::New("hCNRfluxCut", HistAxis(AXrig, "Events/Bin"));
    Hist * hKFRfluxCut = Hist::New("hKFRfluxCut", HistAxis(AXrig, "Events/Bin"));
    Hist * hHCRfluxCut = Hist::New("hHCRfluxCut", HistAxis(AXrig, "Events/Bin"));
    
    Axis AXRchi2("Log-Chi-square [1]", 2000, -8.0, 8.0);
    Hist * hCKPflux = Hist::New("hCKPflux", HistAxis(AXrig, AXRchi2));
    Hist * hCNPflux = Hist::New("hCNPflux", HistAxis(AXrig, AXRchi2));
    Hist * hKFPflux = Hist::New("hKFPflux", HistAxis(AXrig, AXRchi2));
    Hist * hHCPflux = Hist::New("hHCPflux", HistAxis(AXrig, AXRchi2));
    
    Hist * hCKNflux = Hist::New("hCKNflux", HistAxis(AXrig, AXRchi2));
    Hist * hCNNflux = Hist::New("hCNNflux", HistAxis(AXrig, AXRchi2));
    Hist * hKFNflux = Hist::New("hKFNflux", HistAxis(AXrig, AXRchi2));
    Hist * hHCNflux = Hist::New("hHCNflux", HistAxis(AXrig, AXRchi2));

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
            
        Bool_t hasMCL1 = false;
        Bool_t hasMCL9 = false;
        for (auto&& mchit : fG4mc->primPart.hits) {
            if (mchit.layJ == 1) hasMCL1 = true;
            if (mchit.layJ == 9) hasMCL9 = true;
        }
        
        Bool_t hasL1 = false;
        Bool_t hasL9 = false;
        TrFitPar fitPar(type);
        for (auto&& hit : track.hits) {
            HitSt mhit(hit.side%2==1, hit.side/2==1);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
            mhit.set_err(hit.nsr[0], hit.nsr[1], type);
          
            if (hit.layJ >= 2 && hit.layJ <= 8) fitPar.addHit(mhit);
            else {
                if (optL1 && hit.layJ == 1) { hasL1 = true; fitPar.addHit(mhit); }
                if (optL9 && hit.layJ == 9) { hasL9 = true; fitPar.addHit(mhit); }
            }
        }
        Short_t cutNHit = 4 + optL1 + optL9;
        if (fitPar.numOfHit() <= cutNHit) continue;

        if (optL1 && !(hasL1 && hasMCL1)) continue;
        if (optL9 && !(hasL9 && hasMCL9)) continue;
        Short_t patt = (optL1 + optL9 * 2);

        SegPARTMCInfo* topmc = nullptr;
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec == 0 && seg.lay >= 2-optL1) { topmc = &seg; break; } }
        if (topmc == nullptr) continue;

        Double_t mc_mom  = topmc->mom;
        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t bincen  = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
      
        //-------------------------------------//
        MGClock::HrsStopwatch sw; sw.start();
        //SimpleTrFit tr(fitPar);
        PhyTrFit tr(fitPar);
        sw.stop();
        Bool_t hc_succ = tr.status();
        //Bool_t hc_succ = false;
        Double_t hc_irig = tr.part().irig();
        hHCtme->fillH2D(mc_mom, sw.time()*1.0e3);
        //-------------------------------------//
        
        Bool_t ck_succ = track.status[0][patt];
        Bool_t cn_succ = track.status[1][patt];
        Bool_t kf_succ = track.status[2][patt];
        //Bool_t hc_succ = track.status[3][patt];

        hevt->fillH1D(mc_mom);
        if (ck_succ) hCKnum->fillH1D(mc_mom);
        if (cn_succ) hCNnum->fillH1D(mc_mom);
        if (kf_succ) hKFnum->fillH1D(mc_mom);
        if (hc_succ) hHCnum->fillH1D(mc_mom);
        
        if (ck_succ) hCKtme->fillH2D(mc_mom, track.cpuTime[0][patt]);
        if (cn_succ) hCNtme->fillH2D(mc_mom, track.cpuTime[1][patt]*0.1);
        if (kf_succ) hKFtme->fillH2D(mc_mom, track.cpuTime[2][patt]*0.01);
        if (hc_succ) hHCtme->fillH2D(mc_mom, track.cpuTime[3][patt]);

        Double_t ck_irig = (ck_succ ? MGMath::ONE/track.rigidity[0][patt] : 0.);
        Double_t cn_irig = (cn_succ ? MGMath::ONE/track.rigidity[1][patt] : 0.);
        Double_t kf_irig = (kf_succ ? MGMath::ONE/track.rigidity[2][patt] : 0.);
        //Double_t hc_irig = (hc_succ ? MGMath::ONE/track.rigidity[3][patt] : 0.);

        Double_t ck_rig = (ck_succ ? 1.0/ck_irig : 0.);
        Double_t cn_rig = (cn_succ ? 1.0/cn_irig : 0.);
        Double_t kf_rig = (kf_succ ? 1.0/kf_irig : 0.);
        Double_t hc_rig = (hc_succ ? 1.0/hc_irig : 0.);
        
        Double_t ck_chix = (ck_succ ? std::log(track.chisq[0][patt][0]) : 0.); 
        Double_t cn_chix = (cn_succ ? std::log(track.chisq[1][patt][0]) : 0.); 
        Double_t kf_chix = (kf_succ ? std::log(track.chisq[2][patt][0]) : 0.); 
        //Double_t hc_chix = (hc_succ ? std::log(track.chisq[3][patt][0]) : 0.); 
        Double_t hc_chix = (hc_succ ? std::log(tr.nchix())              : 0.); 
        
        Double_t ck_chiy = (ck_succ ? std::log(track.chisq[0][patt][1]) : 0.); 
        Double_t cn_chiy = (cn_succ ? std::log(track.chisq[1][patt][1]) : 0.); 
        Double_t kf_chiy = (kf_succ ? std::log(track.chisq[2][patt][1]) : 0.); 
        //Double_t hc_chiy = (hc_succ ? std::log(track.chisq[3][patt][1]) : 0.); 
        Double_t hc_chiy = (hc_succ ? std::log(tr.nchiy())              : 0.); 

        if (ck_succ) hCKRrso->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_succ) hCNRrso->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (kf_succ) hKFRrso->fillH2D(mc_mom, bincen * (kf_irig - mc_irig));
        if (hc_succ) hHCRrso->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_succ && IntType == 0) hCKRrsoI0->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_succ && IntType == 0) hCNRrsoI0->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (kf_succ && IntType == 0) hKFRrsoI0->fillH2D(mc_mom, bincen * (kf_irig - mc_irig));
        if (hc_succ && IntType == 0) hHCRrsoI0->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_succ && IntType == 1) hCKRrsoI1->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_succ && IntType == 1) hCNRrsoI1->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (kf_succ && IntType == 1) hKFRrsoI1->fillH2D(mc_mom, bincen * (kf_irig - mc_irig));
        if (hc_succ && IntType == 1) hHCRrsoI1->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_succ) hCKRchix->fillH2D(mc_mom, ck_chix);
        if (cn_succ) hCNRchix->fillH2D(mc_mom, cn_chix);
        if (kf_succ) hKFRchix->fillH2D(mc_mom, kf_chix);
        if (hc_succ) hHCRchix->fillH2D(mc_mom, hc_chix);
        
        if (ck_succ && IntType == 0) hCKRchixI0->fillH2D(mc_mom, ck_chix);
        if (cn_succ && IntType == 0) hCNRchixI0->fillH2D(mc_mom, cn_chix);
        if (kf_succ && IntType == 0) hKFRchixI0->fillH2D(mc_mom, kf_chix);
        if (hc_succ && IntType == 0) hHCRchixI0->fillH2D(mc_mom, hc_chix);
        if (ck_succ && IntType == 1) hCKRchixI1->fillH2D(mc_mom, ck_chix);
        if (cn_succ && IntType == 1) hCNRchixI1->fillH2D(mc_mom, cn_chix);
        if (kf_succ && IntType == 1) hKFRchixI1->fillH2D(mc_mom, kf_chix);
        if (hc_succ && IntType == 1) hHCRchixI1->fillH2D(mc_mom, hc_chix);
        
        if (ck_succ) hCKRchiy->fillH2D(mc_mom, ck_chiy);
        if (cn_succ) hCNRchiy->fillH2D(mc_mom, cn_chiy);
        if (kf_succ) hKFRchiy->fillH2D(mc_mom, kf_chiy);
        if (hc_succ) hHCRchiy->fillH2D(mc_mom, hc_chiy);
        
        if (ck_succ && IntType == 0) hCKRchiyI0->fillH2D(mc_mom, ck_chiy);
        if (cn_succ && IntType == 0) hCNRchiyI0->fillH2D(mc_mom, cn_chiy);
        if (kf_succ && IntType == 0) hKFRchiyI0->fillH2D(mc_mom, kf_chiy);
        if (hc_succ && IntType == 0) hHCRchiyI0->fillH2D(mc_mom, hc_chiy);
        if (ck_succ && IntType == 1) hCKRchiyI1->fillH2D(mc_mom, ck_chiy);
        if (cn_succ && IntType == 1) hCNRchiyI1->fillH2D(mc_mom, cn_chiy);
        if (kf_succ && IntType == 1) hKFRchiyI1->fillH2D(mc_mom, kf_chiy);
        if (hc_succ && IntType == 1) hHCRchiyI1->fillH2D(mc_mom, hc_chiy);
        
        const Double_t initR = 30.;
        //const Double_t initR = 1000.;
        Double_t pow27 = std::pow((mc_mom/initR), -1.7);
        Double_t powf = pow27;
        
        if (mc_mom > initR && ck_succ) hCKMCflux->fillH1D(mc_mom, powf);
        if (mc_mom > initR && cn_succ) hCNMCflux->fillH1D(mc_mom, powf);
        if (mc_mom > initR && kf_succ) hKFMCflux->fillH1D(mc_mom, powf);
        if (mc_mom > initR && hc_succ) hHCMCflux->fillH1D(mc_mom, powf);
        
        if (mc_mom > initR && ck_succ) hCKIRflux->fillH2D(ck_irig, ck_chiy, powf);
        if (mc_mom > initR && cn_succ) hCNIRflux->fillH2D(cn_irig, cn_chiy, powf);
        if (mc_mom > initR && kf_succ) hKFIRflux->fillH2D(kf_irig, kf_chiy, powf);
        if (mc_mom > initR && hc_succ) hHCIRflux->fillH2D(hc_irig, hc_chiy, powf);
        
        if (mc_mom > initR && ck_succ && ck_rig > 0) hCKRflux->fillH1D(ck_rig, powf);
        if (mc_mom > initR && cn_succ && cn_rig > 0) hCNRflux->fillH1D(cn_rig, powf);
        if (mc_mom > initR && kf_succ && kf_rig > 0) hKFRflux->fillH1D(kf_rig, powf);
        if (mc_mom > initR && hc_succ && hc_rig > 0) hHCRflux->fillH1D(hc_rig, powf);
        
        if (mc_mom > initR && ck_succ && ck_chiy < 1.90) hCKMCfluxCut->fillH1D(mc_mom, powf);
        if (mc_mom > initR && cn_succ && cn_chiy < 1.90) hCNMCfluxCut->fillH1D(mc_mom, powf);
        if (mc_mom > initR && kf_succ && kf_chiy < 1.90) hKFMCfluxCut->fillH1D(mc_mom, powf);
        if (mc_mom > initR && hc_succ && hc_chiy < 1.15) hHCMCfluxCut->fillH1D(mc_mom, powf);
        
        if (mc_mom > initR && ck_succ && ck_rig > 0 && ck_chiy < 1.90) hCKRfluxCut->fillH1D(ck_rig, powf);
        if (mc_mom > initR && cn_succ && cn_rig > 0 && cn_chiy < 1.90) hCNRfluxCut->fillH1D(cn_rig, powf);
        if (mc_mom > initR && kf_succ && kf_rig > 0 && kf_chiy < 1.90) hKFRfluxCut->fillH1D(kf_rig, powf);
        if (mc_mom > initR && hc_succ && hc_rig > 0 && hc_chiy < 1.15) hHCRfluxCut->fillH1D(hc_rig, powf);

        if (mc_mom > initR && ck_succ && ck_rig > 0) hCKPflux->fillH2D(ck_rig, ck_chiy, powf);
        if (mc_mom > initR && cn_succ && cn_rig > 0) hCNPflux->fillH2D(cn_rig, cn_chiy, powf);
        if (mc_mom > initR && kf_succ && kf_rig > 0) hKFPflux->fillH2D(kf_rig, kf_chiy, powf);
        if (mc_mom > initR && hc_succ && hc_rig > 0) hHCPflux->fillH2D(hc_rig, hc_chiy, powf);
        
        if (mc_mom > initR && ck_succ && ck_rig < 0) hCKNflux->fillH2D(-ck_rig, ck_chiy, powf);
        if (mc_mom > initR && cn_succ && cn_rig < 0) hCNNflux->fillH2D(-cn_rig, cn_chiy, powf);
        if (mc_mom > initR && kf_succ && kf_rig < 0) hKFNflux->fillH2D(-kf_rig, kf_chiy, powf);
        if (mc_mom > initR && hc_succ && hc_rig < 0) hHCNflux->fillH2D(-hc_rig, hc_chiy, powf);
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
