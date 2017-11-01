//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/data3/hchou/AMSCore/prod/17Oct30/src/ClassDef.h"
#include "/data3/hchou/AMSCore/prod/17Oct30/src/ClassDef.C"

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
    //if (opt.type() == "ISS")
    //    dst->SetBranchAddress("rti",  &fRti);
    //dst->SetBranchAddress("trg",  &fTrg);
    dst->SetBranchAddress("tof",  &fTof);
    //dst->SetBranchAddress("acc",  &fAcc);
    dst->SetBranchAddress("trk",  &fTrk);
    //dst->SetBranchAddress("trd",  &fTrd);
    //dst->SetBranchAddress("rich", &fRich);
    //dst->SetBranchAddress("ecal", &fEcal);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    PhyArg::SetOpt(true, true);
    Bool_t optL1 = false;
    Bool_t optL9 = false;
    
    TFile * ofle = new TFile(Form("%s/fit_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    //Axis AXmom("Momentum [GeV]", 200, 1.0, 800., AxisScale::kLog);
    Axis AXmom("Momentum [GeV]", 200, 0.5, 800., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 50, 20., 16000., AxisScale::kLog);

    // Hit
    Axis AXres("Residual [10^{-4} cm]", 800, -300., 300.);
    Hist * hXres = Hist::New("hXres", "hXres", HistAxis(AXmom, AXres));
    Hist * hYres = Hist::New("hYres", "hYres", HistAxis(AXmom, AXres));
    
    // Fit
    Axis AXrso("(1/Rm - 1/Rt) [1/GV]", 800, -1.5, 1.5);
    Hist * hCKrso = Hist::New("hCKrso", "hCKrso", HistAxis(AXmom, AXrso));
    Hist * hCNrso = Hist::New("hCNrso", "hCNrso", HistAxis(AXmom, AXrso));
    Hist * hHCrso = Hist::New("hHCrso", "hHCrso", HistAxis(AXmom, AXrso));

    Long64_t printRate = dst->GetEntries()/40;
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        //if (entry>(dst->GetEntries()/20)) break; // testcode
        //if (entry%10!=0) continue; // testcode
        dst->GetEntry(entry);
        
        Double_t mc_mom  = (fG4mc->primPart.mom);
        Double_t bincen = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
        
        if (fTof->betaH < 0.3 || fTof->betaH > 1.3) continue;
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        if (fTof->Qall < 0.8 || fTof->Qall > 1.8) continue;
        if (fTof->betaHPatt != 15) continue;
        
        if (fTrk->tracks.size() != 1) continue;
        TrackInfo& track = fTrk->tracks.at(0);
        
        Bool_t hasL1 = false;
        Bool_t hasL9 = false;
        std::vector<HitSt> mhits;
        for (auto&& hit : track.hits) {
            if (hit.side != 3) continue;
            if (hit.layJ == 1) hasL1 = true;
            if (hit.layJ == 9) hasL9 = true;
            Bool_t isInn = (hit.layJ >= 2 || hit.layJ <= 8);
            HitSt mhit(hit.coo[0], hit.coo[1], hit.coo[2]);
           
            if (isInn || (optL1 && hasL1) || (optL9 && hasL9)) {
                mhits.push_back(mhit);
            }

            for (auto&& mchit : fG4mc->primPart.hits) {
                if (mchit.layJ < 2 || mchit.layJ > 8) continue;
                if (mchit.layJ != hit.layJ) continue;
                Double_t resx = 1.0e4 * (hit.coo[0] - mchit.coo[0]);
                Double_t resy = 1.0e4 * (hit.coo[1] - mchit.coo[1]);
                hXres->fill(mc_mom, resx);
                hYres->fill(mc_mom, resy);
            }
        }
        Short_t cutNHit = 4 + optL1 + optL9;
        if (mhits.size() <= cutNHit) continue;

        if (optL1 && !hasL1) continue;
        if (optL9 && !hasL9) continue;

        Short_t patt = (optL1 + optL9 * 2);
        Double_t mc_irig = (fG4mc->primPart.chrg / fG4mc->primPart.mom);
        Double_t ck_irig = ((track.status[0][patt]) ? MGMath::ONE/track.rigidity[0][patt] : 0.);
        Double_t cn_irig = ((track.status[1][patt]) ? MGMath::ONE/track.rigidity[1][patt] : 0.);
        
        PhyTr tr(mhits);
        tr.fit();
        //tr.print();
        Double_t hc_irig = tr.part().irig();
        
        hCKrso->fill(mc_mom, bincen * (ck_irig - mc_irig));
        hCNrso->fill(mc_mom, bincen * (cn_irig - mc_irig));
        hHCrso->fill(mc_mom, bincen * (hc_irig - mc_irig));
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
