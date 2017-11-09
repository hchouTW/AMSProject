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
    Bool_t optL1 = true;
    Bool_t optL9 = false;
    
    TFile * ofle = new TFile(Form("%s/fit_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 40, 0.5, 10., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 40, 1.0, 800., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 100, 0.5, 800., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 50, 20., 16000., AxisScale::kLog);
    
    Axis AXeta("Eta [1]", 50, 0.01, 3.0, AxisScale::kLog);
    Axis AXbta("1/Bta [1]", 50,  1.0, 3.0, AxisScale::kLog);

    const Double_t cm2um = 10000.;
    const Double_t gev2mev = 1000.;

    // Hit
    Axis AXres("Residual [10^{-4} cm]", 800, -300., 300.);
    Hist * hXres = Hist::New("hXres", "hXres", HistAxis(AXmom, AXres));
    Hist * hYres = Hist::New("hYres", "hYres", HistAxis(AXmom, AXres));
    
    Axis AXedep("Edep [MeV]", 1600, 0.0, 1.0);
    Hist * hMGedep = Hist::New("hMGedep", "hMGedep", HistAxis(AXeta, AXedep));
    Hist * hMBedep = Hist::New("hMBedep", "hMBedep", HistAxis(AXbta, AXedep));
    
    // Fit
    Axis AXrso("(1/Rm - 1/Rt) [1/GV]", 800, -1.5, 1.5);
    Hist * hCKrso = Hist::New("hCKrso", "hCKrso", HistAxis(AXmom, AXrso));
    Hist * hCNrso = Hist::New("hCNrso", "hCNrso", HistAxis(AXmom, AXrso));
    Hist * hHCrso = Hist::New("hHCrso", "hHCrso", HistAxis(AXmom, AXrso));
    
    Axis AXresx("[um]", 800, -300., 300.);
    Hist * hCKresx = Hist::New("hCKresx", "hCKresx", HistAxis(AXmom, AXresx));
    Hist * hCNresx = Hist::New("hCNresx", "hCNresx", HistAxis(AXmom, AXresx));
    Hist * hHCresx = Hist::New("hHCresx", "hHCresx", HistAxis(AXmom, AXresx));
    
    Axis AXresy("[um]", 800, -300., 300.);
    Hist * hCKresy = Hist::New("hCKresy", "hCKresy", HistAxis(AXmom, AXresy));
    Hist * hCNresy = Hist::New("hCNresy", "hCNresy", HistAxis(AXmom, AXresy));
    Hist * hHCresy = Hist::New("hHCresy", "hHCresy", HistAxis(AXmom, AXresy));

    Long64_t printRate = dst->GetEntries()/40;
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        //if (entry>(dst->GetEntries()/20)) break; // testcode
        //if (entry%10!=0) continue; // testcode
        dst->GetEntry(entry);
        
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
            //if (hit.side != 3) continue;
            if (hit.layJ == 1) hasL1 = true;
            if (hit.layJ == 9) hasL9 = true;
            Bool_t isInn = (hit.layJ >= 2 && hit.layJ <= 8);
            HitSt mhit(hit.side%2==1, hit.side/2==1);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
          
            if (isInn) mhits.push_back(mhit);
            else {
                if (optL1 && hit.layJ == 1) mhits.push_back(mhit);
                if (optL9 && hit.layJ == 9) mhits.push_back(mhit);
            }
            
            Bool_t hasMCL1 = false;
            Bool_t hasMCL9 = false;
            for (auto&& mchit : fG4mc->primPart.hits) {
                if (mchit.layJ == 1) hasMCL1 = true;
                if (mchit.layJ == 9) hasMCL9 = true;
                
                if (mchit.layJ < 2 || mchit.layJ > 8) continue;
                if (mchit.layJ != hit.layJ) continue;
                Double_t resx = 1.0e4 * (hit.coo[0] - mchit.coo[0]);
                Double_t resy = 1.0e4 * (hit.coo[1] - mchit.coo[1]);
                hXres->fill(mchit.mom, resx);
                hYres->fill(mchit.mom, resy);

                Double_t eta = fG4mc->primPart.mass/mchit.mom;
                Double_t ibta = std::sqrt(1.+eta*eta);
                Double_t edep = gev2mev * mchit.edep * std::fabs(mchit.dir[2]); // [MeV]
                hMGedep->fill(eta,  edep);
                hMBedep->fill(ibta, edep);
            }
            if (hasL1 && !hasMCL1) hasL1 = false;
            if (hasL9 && !hasMCL9) hasL9 = false;
        }
        Short_t cutNHit = 4 + optL1 + optL9;
        if (mhits.size() <= cutNHit) continue;

        if (optL1 && !hasL1) continue;
        if (optL9 && !hasL9) continue;

        HitTRKMCInfo* topmc = (fG4mc->primPart.hits.size() == 0) ? nullptr : &fG4mc->primPart.hits.at(0);
        if (topmc == nullptr) continue;
        Int_t    mc_lay = topmc->layJ-1;
        Double_t mc_mom = topmc->mom;
        Double_t bincen = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);

        Short_t patt = (optL1 + optL9 * 2);
        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t ck_irig = ((track.status[0][patt]) ? MGMath::ONE/track.rigidity[0][patt] : 0.);
        Double_t cn_irig = ((track.status[1][patt]) ? MGMath::ONE/track.rigidity[1][patt] : 0.);
        
        PhyTr tr(mhits);
        //Bool_t succ = tr.fit();
        Bool_t succ = false;
        //tr.print();
        Double_t hc_irig = tr.part().irig();
        
        if (track.status[0][patt]) hCKrso->fill(mc_mom, bincen * (ck_irig - mc_irig));
        if (track.status[1][patt]) hCNrso->fill(mc_mom, bincen * (cn_irig - mc_irig));
        if (succ) hHCrso->fill(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (track.status[0][patt]) hCKresx->fill(mc_mom, cm2um * (track.stateLJ[0][patt][mc_lay][0] - topmc->coo[0]));
        if (track.status[1][patt]) hCNresx->fill(mc_mom, cm2um * (track.stateLJ[1][patt][mc_lay][0] - topmc->coo[0]));
        if (succ) hHCresx->fill(mc_mom, tr.part().cx() - topmc->coo[0]);
        
        if (track.status[0][patt]) hCKresy->fill(mc_mom, cm2um * (track.stateLJ[0][patt][mc_lay][1] - topmc->coo[1]));
        if (track.status[1][patt]) hCNresy->fill(mc_mom, cm2um * (track.stateLJ[1][patt][mc_lay][1] - topmc->coo[1]));
        if (succ) hHCresy->fill(mc_mom, cm2um * (tr.part().cy() - topmc->coo[1]));
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
