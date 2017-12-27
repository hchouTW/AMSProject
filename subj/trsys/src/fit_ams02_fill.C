//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/data3/hchou/AMSCore/prod/17Dec23/src/ClassDef.h"
#include "/data3/hchou/AMSCore/prod/17Dec23/src/ClassDef.C"

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
    PhyArg::SetOpt(false, true);
    Bool_t optL1 = false;
    Bool_t optL9 = false;
    
    TFile * ofle = new TFile(Form("%s/fit_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 100, 0.5, 4000., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 10, 20., 100., AxisScale::kLog);
    Axis AXimm("1/Momentum [1/GeV]", AXmom, 1, true);

    // Fit
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 2000, -10.0, 10.0);
    Hist * hCKRrso = Hist::New("hCKRrso", "hCKRrso", HistAxis(AXmom, AXRrso));
    Hist * hCNRrso = Hist::New("hCNRrso", "hCNRrso", HistAxis(AXmom, AXRrso));
    Hist * hHCRrso = Hist::New("hHCRrso", "hHCRrso", HistAxis(AXmom, AXRrso));
    
    Hist * hCKRrsoI0 = Hist::New("hCKRrsoI0", "hCKRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI0 = Hist::New("hCNRrsoI0", "hCNRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI0 = Hist::New("hHCRrsoI0", "hHCRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI1 = Hist::New("hCKRrsoI1", "hCKRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI1 = Hist::New("hCNRrsoI1", "hCNRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI1 = Hist::New("hHCRrsoI1", "hHCRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI2 = Hist::New("hCKRrsoI2", "hCKRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI2 = Hist::New("hCNRrsoI2", "hCNRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI2 = Hist::New("hHCRrsoI2", "hHCRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI3 = Hist::New("hCKRrsoI3", "hCKRrsoI3", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI3 = Hist::New("hCNRrsoI3", "hCNRrsoI3", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI3 = Hist::New("hHCRrsoI3", "hHCRrsoI3", HistAxis(AXmom, AXRrso));
    
    Axis AXRchi("Log-Chi-square [1]", 800, -8.0, 8.0);
    Hist * hCKRchi = Hist::New("hCKRchi", "hCKRchi", HistAxis(AXmom, AXRchi));
    Hist * hCNRchi = Hist::New("hCNRchi", "hCNRchi", HistAxis(AXmom, AXRchi));
    Hist * hHCRchi = Hist::New("hHCRchi", "hHCRchi", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRchiI0 = Hist::New("hCKRchiI0", "hCKRchiI0", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiI0 = Hist::New("hCNRchiI0", "hCNRchiI0", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiI0 = Hist::New("hHCRchiI0", "hHCRchiI0", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiI1 = Hist::New("hCKRchiI1", "hCKRchiI1", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiI1 = Hist::New("hCNRchiI1", "hCNRchiI1", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiI1 = Hist::New("hHCRchiI1", "hHCRchiI1", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiI2 = Hist::New("hCKRchiI2", "hCKRchiI2", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiI2 = Hist::New("hCNRchiI2", "hCNRchiI2", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiI2 = Hist::New("hHCRchiI2", "hHCRchiI2", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiI3 = Hist::New("hCKRchiI3", "hCKRchiI3", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiI3 = Hist::New("hCNRchiI3", "hCNRchiI3", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiI3 = Hist::New("hHCRchiI3", "hHCRchiI3", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRrsoCut = Hist::New("hCKRrsoCut", "hCKRrsoCut", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCut = Hist::New("hCNRrsoCut", "hCNRrsoCut", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCut = Hist::New("hHCRrsoCut", "hHCRrsoCut", HistAxis(AXmom, AXRrso));
    
    Hist * hCKRrsoCutI0 = Hist::New("hCKRrsoCutI0", "hCKRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI0 = Hist::New("hCNRrsoCutI0", "hCNRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI0 = Hist::New("hHCRrsoCutI0", "hHCRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI1 = Hist::New("hCKRrsoCutI1", "hCKRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI1 = Hist::New("hCNRrsoCutI1", "hCNRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI1 = Hist::New("hHCRrsoCutI1", "hHCRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI2 = Hist::New("hCKRrsoCutI2", "hCKRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI2 = Hist::New("hCNRrsoCutI2", "hCNRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI2 = Hist::New("hHCRrsoCutI2", "hHCRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI3 = Hist::New("hCKRrsoCutI3", "hCKRrsoCutI3", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI3 = Hist::New("hCNRrsoCutI3", "hCNRrsoCutI3", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI3 = Hist::New("hHCRrsoCutI3", "hHCRrsoCutI3", HistAxis(AXmom, AXRrso));
    
    Hist * hCKflux = Hist::New("hCKflux", "hCKflux", HistAxis(AXimm));
    Hist * hCNflux = Hist::New("hCNflux", "hCNflux", HistAxis(AXimm));
    Hist * hHCflux = Hist::New("hHCflux", "hHCflux", HistAxis(AXimm));
    Hist * hCKflux2 = Hist::New("hCKflux2", "hCKflux2", HistAxis(AXimm));
    Hist * hCNflux2 = Hist::New("hCNflux2", "hCNflux2", HistAxis(AXimm));
    Hist * hHCflux2 = Hist::New("hHCflux2", "hHCflux2", HistAxis(AXimm));

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
        }
            
        Bool_t hasMCL1 = false;
        Bool_t hasMCL9 = false;
        for (auto&& mchit : fG4mc->primPart.hits) {
            if (mchit.layJ == 1) hasMCL1 = true;
            if (mchit.layJ == 9) hasMCL9 = true;
        }
        
        Bool_t hasL1 = false;
        Bool_t hasL9 = false;
        std::vector<HitSt> mhits;
        for (auto&& hit : track.hits) {
            HitSt mhit(hit.side%2==1, hit.side/2==1);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
            mhit.set_err(hit.nsr[0], hit.nsr[1]);
          
            if (hit.layJ >= 2 && hit.layJ <= 8) mhits.push_back(mhit);
            else {
                if (optL1 && hit.layJ == 1) { hasL1 = true; mhits.push_back(mhit); }
                if (optL9 && hit.layJ == 9) { hasL9 = true; mhits.push_back(mhit); }
            }
        }
        Short_t cutNHit = 4 + optL1 + optL9;
        if (mhits.size() <= cutNHit) continue;

        if (optL1 && !(hasL1 && hasMCL1)) continue;
        if (optL9 && !(hasL9 && hasMCL9)) continue;
        Short_t patt = (optL1 + optL9 * 2);

        SegPARTMCInfo* topmc = nullptr;
        for (auto&& seg : fG4mc->primPart.segs) { if (seg.dec == 0) { topmc = &seg; break; } }
        if (topmc == nullptr) continue;

        Double_t mc_mom = topmc->mom;
        Double_t bincen = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
       
        //-------------------------------------//
        PhyTr tr(mhits);
        Bool_t hc_succ = tr.fit();
        //Bool_t hc_succ = false;
        Double_t hc_irig = tr.part().irig();
        //-------------------------------------//
        Bool_t ck_succ = track.status[0][patt];
        Bool_t cn_succ = track.status[1][patt];

        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t ck_irig = (ck_succ ? MGMath::ONE/track.rigidity[0][patt] : 0.);
        Double_t cn_irig = (cn_succ ? MGMath::ONE/track.rigidity[1][patt] : 0.);
        
        Double_t ck_lchi = (ck_succ ? std::log(track.chisq[0][patt][1]) : 0.); 
        Double_t cn_lchi = (cn_succ ? std::log(track.chisq[1][patt][1]) : 0.); 
        Double_t hc_lchi = (hc_succ ? std::log(tr.nchi())               : 0.); 

        Double_t ck_cut = (ck_succ ? (ck_lchi < 2.0) : false); // 96%
        Double_t cn_cut = (cn_succ ? (cn_lchi < 2.0) : false); // 96%
        Double_t hc_cut = (hc_succ ? (hc_lchi < 0.9) : false); // 96%

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
        
        if (ck_succ) hCKRchi->fillH2D(mc_mom, ck_lchi);
        if (cn_succ) hCNRchi->fillH2D(mc_mom, cn_lchi);
        if (hc_succ) hHCRchi->fillH2D(mc_mom, hc_lchi);
        
        if (ck_succ && IntType == 0) hCKRchiI0->fillH2D(mc_mom, ck_lchi);
        if (cn_succ && IntType == 0) hCNRchiI0->fillH2D(mc_mom, cn_lchi);
        if (hc_succ && IntType == 0) hHCRchiI0->fillH2D(mc_mom, hc_lchi);
        if (ck_succ && IntType == 1) hCKRchiI1->fillH2D(mc_mom, ck_lchi);
        if (cn_succ && IntType == 1) hCNRchiI1->fillH2D(mc_mom, cn_lchi);
        if (hc_succ && IntType == 1) hHCRchiI1->fillH2D(mc_mom, hc_lchi);
        if (ck_succ && IntType == 2) hCKRchiI2->fillH2D(mc_mom, ck_lchi);
        if (cn_succ && IntType == 2) hCNRchiI2->fillH2D(mc_mom, cn_lchi);
        if (hc_succ && IntType == 2) hHCRchiI2->fillH2D(mc_mom, hc_lchi);
        if (ck_succ && IntType == 3) hCKRchiI3->fillH2D(mc_mom, ck_lchi);
        if (cn_succ && IntType == 3) hCNRchiI3->fillH2D(mc_mom, cn_lchi);
        if (hc_succ && IntType == 3) hHCRchiI3->fillH2D(mc_mom, hc_lchi);
        
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
        
        Double_t pow27 = std::pow(20., 1.7) * std::pow(mc_mom, -1.7);
        if (mc_mom > 20. && ck_succ) hCKflux->fillH1D(ck_irig, pow27);
        if (mc_mom > 20. && cn_succ) hCNflux->fillH1D(cn_irig, pow27);
        if (mc_mom > 20. && hc_succ) hHCflux->fillH1D(hc_irig, pow27);
        
        if (mc_mom > 20. && ck_cut) hCKflux2->fillH1D(ck_irig, pow27);
        if (mc_mom > 20. && cn_cut) hCNflux2->fillH1D(cn_irig, pow27);
        if (mc_mom > 20. && hc_cut) hHCflux2->fillH1D(hc_irig, pow27);
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
