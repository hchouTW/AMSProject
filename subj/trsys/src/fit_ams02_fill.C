//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

//#include "/ams_home/hchou/AMSCore/prod/17Dec23/src/ClassDef.h"
//#include "/ams_home/hchou/AMSCore/prod/17Dec23/src/ClassDef.C"

#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/17Dec23/src/ClassDef.h"
#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/17Dec23/src/ClassDef.C"

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
    //Axis AXmom("Momentum [GeV]", 40, 0.25, 200., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 40, 200., 4000., AxisScale::kLog);
    Axis AXimm("1/Rigidity [1/GV]", AXmom, 1, true);

    // Eloss Frac
    Axis AXfrac("Eloss Frac [1]", 50, 0., 1.);
    Hist * hMCfrac = Hist::New("hMCfrac", HistAxis(AXmom, AXfrac));

    // Fit
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 2000, -10.0, 10.0);
    Hist * hCKRrso = Hist::New("hCKRrso", HistAxis(AXmom, AXRrso));
    Hist * hCNRrso = Hist::New("hCNRrso", HistAxis(AXmom, AXRrso));
    Hist * hHCRrso = Hist::New("hHCRrso", HistAxis(AXmom, AXRrso));
    
    Hist * hCKRrsoI0 = Hist::New("hCKRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI0 = Hist::New("hCNRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI0 = Hist::New("hHCRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI1 = Hist::New("hCKRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI1 = Hist::New("hCNRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI1 = Hist::New("hHCRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI2 = Hist::New("hCKRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI2 = Hist::New("hCNRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI2 = Hist::New("hHCRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI3 = Hist::New("hCKRrsoI3", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI3 = Hist::New("hCNRrsoI3", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI3 = Hist::New("hHCRrsoI3", HistAxis(AXmom, AXRrso));
    
    Axis AXRchi("Log-Chi-square [1]", 800, -8.0, 8.0);
    Hist * hCKRchix = Hist::New("hCKRchix", HistAxis(AXmom, AXRchi));
    Hist * hCNRchix = Hist::New("hCNRchix", HistAxis(AXmom, AXRchi));
    Hist * hHCRchix = Hist::New("hHCRchix", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRchixI0 = Hist::New("hCKRchixI0", HistAxis(AXmom, AXRchi));
    Hist * hCNRchixI0 = Hist::New("hCNRchixI0", HistAxis(AXmom, AXRchi));
    Hist * hHCRchixI0 = Hist::New("hHCRchixI0", HistAxis(AXmom, AXRchi));
    Hist * hCKRchixI1 = Hist::New("hCKRchixI1", HistAxis(AXmom, AXRchi));
    Hist * hCNRchixI1 = Hist::New("hCNRchixI1", HistAxis(AXmom, AXRchi));
    Hist * hHCRchixI1 = Hist::New("hHCRchixI1", HistAxis(AXmom, AXRchi));
    Hist * hCKRchixI2 = Hist::New("hCKRchixI2", HistAxis(AXmom, AXRchi));
    Hist * hCNRchixI2 = Hist::New("hCNRchixI2", HistAxis(AXmom, AXRchi));
    Hist * hHCRchixI2 = Hist::New("hHCRchixI2", HistAxis(AXmom, AXRchi));
    Hist * hCKRchixI3 = Hist::New("hCKRchixI3", HistAxis(AXmom, AXRchi));
    Hist * hCNRchixI3 = Hist::New("hCNRchixI3", HistAxis(AXmom, AXRchi));
    Hist * hHCRchixI3 = Hist::New("hHCRchixI3", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRchiy = Hist::New("hCKRchiy", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiy = Hist::New("hCNRchiy", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiy = Hist::New("hHCRchiy", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRchiyI0 = Hist::New("hCKRchiyI0", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiyI0 = Hist::New("hCNRchiyI0", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiyI0 = Hist::New("hHCRchiyI0", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiyI1 = Hist::New("hCKRchiyI1", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiyI1 = Hist::New("hCNRchiyI1", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiyI1 = Hist::New("hHCRchiyI1", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiyI2 = Hist::New("hCKRchiyI2", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiyI2 = Hist::New("hCNRchiyI2", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiyI2 = Hist::New("hHCRchiyI2", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiyI3 = Hist::New("hCKRchiyI3", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiyI3 = Hist::New("hCNRchiyI3", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiyI3 = Hist::New("hHCRchiyI3", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRrsoCut = Hist::New("hCKRrsoCut", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCut = Hist::New("hCNRrsoCut", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCut = Hist::New("hHCRrsoCut", HistAxis(AXmom, AXRrso));
    
    Hist * hCKRrsoCutI0 = Hist::New("hCKRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI0 = Hist::New("hCNRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI0 = Hist::New("hHCRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI1 = Hist::New("hCKRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI1 = Hist::New("hCNRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI1 = Hist::New("hHCRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI2 = Hist::New("hCKRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI2 = Hist::New("hCNRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI2 = Hist::New("hHCRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI3 = Hist::New("hCKRrsoCutI3", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI3 = Hist::New("hCNRrsoCutI3", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI3 = Hist::New("hHCRrsoCutI3", HistAxis(AXmom, AXRrso));
    
    Hist * hCKflux = Hist::New("hCKflux", HistAxis(AXimm));
    Hist * hCNflux = Hist::New("hCNflux", HistAxis(AXimm));
    Hist * hHCflux = Hist::New("hHCflux", HistAxis(AXimm));
    Hist * hCKflux2 = Hist::New("hCKflux2", HistAxis(AXimm));
    Hist * hCNflux2 = Hist::New("hCNflux2", HistAxis(AXimm));
    Hist * hHCflux2 = Hist::New("hHCflux2", HistAxis(AXimm));

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
            if      (std::fabs(fG4mc->primVtx.coo[2]) < 55.) IntType = 1;
            else if (fG4mc->primVtx.coo[2] > 55.)            IntType = 2;
            else if (fG4mc->primVtx.coo[2] > -80.)           IntType = 3;
            hMCfrac->fillH2D(fG4mc->primPart.mom, fG4mc->primVtx.ke.at(0)/fG4mc->primPart.ke);
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
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec == 0) { topmc = &seg; break; } }
        if (topmc == nullptr) continue;

        Double_t mc_mom = topmc->mom;
        Double_t bincen = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
       
        //-------------------------------------//
        SimpleTrFit tr(fitPar);
        Bool_t hc_succ = tr.status();
        //Bool_t hc_succ = false;
        Double_t hc_irig = tr.part().irig();
        //-------------------------------------//
        Bool_t ck_succ = track.status[0][patt];
        Bool_t cn_succ = track.status[1][patt];

        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t ck_irig = (ck_succ ? MGMath::ONE/track.rigidity[0][patt] : 0.);
        Double_t cn_irig = (cn_succ ? MGMath::ONE/track.rigidity[1][patt] : 0.);
        
        Double_t ck_chix = (ck_succ ? std::log(track.chisq[0][patt][0]) : 0.); 
        Double_t cn_chix = (cn_succ ? std::log(track.chisq[1][patt][0]) : 0.); 
        Double_t hc_chix = (hc_succ ? std::log(tr.nchix())              : 0.); 
        
        Double_t ck_chiy = (ck_succ ? std::log(track.chisq[0][patt][1]) : 0.); 
        Double_t cn_chiy = (cn_succ ? std::log(track.chisq[1][patt][1]) : 0.); 
        Double_t hc_chiy = (hc_succ ? std::log(tr.nchiy())              : 0.); 

        Double_t ck_cut = (ck_succ ? (ck_chiy < 2.0) : false); // 96%
        Double_t cn_cut = (cn_succ ? (cn_chiy < 2.0) : false); // 96%
        Double_t hc_cut = (hc_succ ? (hc_chiy < 1.0) : false); // 96%

        if (ck_succ) hCKRrso->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_succ) hCNRrso->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_succ) hHCRrso->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_succ && IntType == 0) hCKRrsoI0->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_succ && IntType == 0) hCNRrsoI0->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_succ && IntType == 0) hHCRrsoI0->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_succ && IntType == 1) hCKRrsoI1->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_succ && IntType == 1) hCNRrsoI1->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_succ && IntType == 1) hHCRrsoI1->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_succ && IntType == 2) hCKRrsoI2->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_succ && IntType == 2) hCNRrsoI2->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_succ && IntType == 2) hHCRrsoI2->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_succ && IntType == 3) hCKRrsoI3->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_succ && IntType == 3) hCNRrsoI3->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_succ && IntType == 3) hHCRrsoI3->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_succ) hCKRchix->fillH2D(mc_mom, ck_chix);
        if (cn_succ) hCNRchix->fillH2D(mc_mom, cn_chix);
        if (hc_succ) hHCRchix->fillH2D(mc_mom, hc_chix);
        
        if (ck_succ && IntType == 0) hCKRchixI0->fillH2D(mc_mom, ck_chix);
        if (cn_succ && IntType == 0) hCNRchixI0->fillH2D(mc_mom, cn_chix);
        if (hc_succ && IntType == 0) hHCRchixI0->fillH2D(mc_mom, hc_chix);
        if (ck_succ && IntType == 1) hCKRchixI1->fillH2D(mc_mom, ck_chix);
        if (cn_succ && IntType == 1) hCNRchixI1->fillH2D(mc_mom, cn_chix);
        if (hc_succ && IntType == 1) hHCRchixI1->fillH2D(mc_mom, hc_chix);
        if (ck_succ && IntType == 2) hCKRchixI2->fillH2D(mc_mom, ck_chix);
        if (cn_succ && IntType == 2) hCNRchixI2->fillH2D(mc_mom, cn_chix);
        if (hc_succ && IntType == 2) hHCRchixI2->fillH2D(mc_mom, hc_chix);
        if (ck_succ && IntType == 3) hCKRchixI3->fillH2D(mc_mom, ck_chix);
        if (cn_succ && IntType == 3) hCNRchixI3->fillH2D(mc_mom, cn_chix);
        if (hc_succ && IntType == 3) hHCRchixI3->fillH2D(mc_mom, hc_chix);
        
        if (ck_succ) hCKRchiy->fillH2D(mc_mom, ck_chiy);
        if (cn_succ) hCNRchiy->fillH2D(mc_mom, cn_chiy);
        if (hc_succ) hHCRchiy->fillH2D(mc_mom, hc_chiy);
        
        if (ck_succ && IntType == 0) hCKRchiyI0->fillH2D(mc_mom, ck_chiy);
        if (cn_succ && IntType == 0) hCNRchiyI0->fillH2D(mc_mom, cn_chiy);
        if (hc_succ && IntType == 0) hHCRchiyI0->fillH2D(mc_mom, hc_chiy);
        if (ck_succ && IntType == 1) hCKRchiyI1->fillH2D(mc_mom, ck_chiy);
        if (cn_succ && IntType == 1) hCNRchiyI1->fillH2D(mc_mom, cn_chiy);
        if (hc_succ && IntType == 1) hHCRchiyI1->fillH2D(mc_mom, hc_chiy);
        if (ck_succ && IntType == 2) hCKRchiyI2->fillH2D(mc_mom, ck_chiy);
        if (cn_succ && IntType == 2) hCNRchiyI2->fillH2D(mc_mom, cn_chiy);
        if (hc_succ && IntType == 2) hHCRchiyI2->fillH2D(mc_mom, hc_chiy);
        if (ck_succ && IntType == 3) hCKRchiyI3->fillH2D(mc_mom, ck_chiy);
        if (cn_succ && IntType == 3) hCNRchiyI3->fillH2D(mc_mom, cn_chiy);
        if (hc_succ && IntType == 3) hHCRchiyI3->fillH2D(mc_mom, hc_chiy);
        
        if (ck_cut) hCKRrsoCut->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut) hCNRrsoCut->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut) hHCRrsoCut->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_cut && IntType == 0) hCKRrsoCutI0->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut && IntType == 0) hCNRrsoCutI0->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut && IntType == 0) hHCRrsoCutI0->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_cut && IntType == 1) hCKRrsoCutI1->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut && IntType == 1) hCNRrsoCutI1->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut && IntType == 1) hHCRrsoCutI1->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_cut && IntType == 2) hCKRrsoCutI2->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut && IntType == 2) hCNRrsoCutI2->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut && IntType == 2) hHCRrsoCutI2->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_cut && IntType == 3) hCKRrsoCutI3->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut && IntType == 3) hCNRrsoCutI3->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut && IntType == 3) hHCRrsoCutI3->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        Double_t pow27 = std::pow((mc_mom/20.), -1.7);
        Double_t pow30 = std::pow((mc_mom/200.), -2.0);
        Double_t powf = pow30;
        if (mc_mom > 20. && ck_succ) hCKflux->fillH1D(ck_irig, powf);
        if (mc_mom > 20. && cn_succ) hCNflux->fillH1D(cn_irig, powf);
        if (mc_mom > 20. && hc_succ) hHCflux->fillH1D(hc_irig, powf);
        
        if (mc_mom > 20. && ck_cut) hCKflux2->fillH1D(ck_irig, powf);
        if (mc_mom > 20. && cn_cut) hCNflux2->fillH1D(cn_irig, powf);
        if (mc_mom > 20. && hc_cut) hHCflux2->fillH1D(hc_irig, powf);
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
