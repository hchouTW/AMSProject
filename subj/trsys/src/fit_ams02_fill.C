//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/data3/hchou/AMSCore/prod/17Nov24/src/ClassDef.h"
#include "/data3/hchou/AMSCore/prod/17Nov24/src/ClassDef.C"

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
    //PhyArg::SetOpt(false, true);
    //PhyArg::SetOpt(false, false);
    //PhyArg::SetOpt(true, false);
    PhyArg::SetOpt(true, true);
    Bool_t optL1 = false;
    Bool_t optL9 = false;
    
    TFile * ofle = new TFile(Form("%s/fit_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    //Axis AXmom("Momentum [GeV]", 40, 0.5, 100., AxisScale::kLog);
    Axis AXmom("Momentum [GeV]", 100, 0.5, 4000., AxisScale::kLog);
    Axis AXimm("1/Momentum [1/GeV]", AXmom, 1, true);
    
    Axis AXeta("Eta [1]", 50, 0.01, 3.0, AxisScale::kLog);
    Axis AXbta("1/Bta [1]", 50,  1.0, 3.0, AxisScale::kLog);

    const Double_t cm2um = 10000.;
    const Double_t gev2mev = 1000.;

    // Hit
    Axis AXres("Residual [10^{-4} cm]", 800, -300., 300.);
    Hist * hXres = Hist::New("hXres", "hXres", HistAxis(AXmom, AXres));
    Hist * hYres = Hist::New("hYres", "hYres", HistAxis(AXmom, AXres));
    Hist * hXres2 = Hist::New("hXres2", "hXres2", HistAxis(AXres));
    Hist * hYres2 = Hist::New("hYres2", "hYres2", HistAxis(AXres));
   
    Axis AXcos("Cos [1]", 20,  0.8, 1.0);
    Axis AXeac("Eta [1]", 40, -0.5, 0.5);
    Hist * hXresnt = Hist::New("hXresnt", "hXresnt", HistAxis(AXeac, AXres));
    Hist * hYresnt = Hist::New("hYresnt", "hYresnt", HistAxis(AXeac, AXres));
    Hist * hXresn1 = Hist::New("hXresn1", "hXresn1", HistAxis(AXeac, AXres));
    Hist * hYresn1 = Hist::New("hYresn1", "hYresn1", HistAxis(AXeac, AXres));
    Hist * hXresn2 = Hist::New("hXresn2", "hXresn2", HistAxis(AXeac, AXres));
    Hist * hYresn2 = Hist::New("hYresn2", "hYresn2", HistAxis(AXeac, AXres));
    Hist * hXresn3 = Hist::New("hXresn3", "hXresn3", HistAxis(AXeac, AXres));
    Hist * hYresn3 = Hist::New("hYresn3", "hYresn3", HistAxis(AXeac, AXres));
    Hist * hXresn4 = Hist::New("hXresn4", "hXresn4", HistAxis(AXeac, AXres));
    Hist * hYresn4 = Hist::New("hYresn4", "hYresn4", HistAxis(AXeac, AXres));
    
    Axis AXedep("Edep [MeV]", 1600, 0.0, 1.0);
    Hist * hMGedep = Hist::New("hMGedep", "hMGedep", HistAxis(AXeta, AXedep));
    Hist * hMBedep = Hist::New("hMBedep", "hMBedep", HistAxis(AXbta, AXedep));
    
    // Fit
    Hist * hCKflux = Hist::New("hCKflux", "hCKflux", HistAxis(AXimm));
    Hist * hCNflux = Hist::New("hCNflux", "hCNflux", HistAxis(AXimm));
    Hist * hHCflux = Hist::New("hHCflux", "hHCflux", HistAxis(AXimm));
    
    Axis AXrso("(1/Rm - 1/Rt) [1/GV]", 1600, -5.0, 5.0);
    Hist * hCKrso = Hist::New("hCKrso", "hCKrso", HistAxis(AXmom, AXrso));
    Hist * hCNrso = Hist::New("hCNrso", "hCNrso", HistAxis(AXmom, AXrso));
    Hist * hHCrso = Hist::New("hHCrso", "hHCrso", HistAxis(AXmom, AXrso));
    Hist * hCKrso2 = Hist::New("hCKrso2", "hCKrso2", HistAxis(AXrso));
    Hist * hCNrso2 = Hist::New("hCNrso2", "hCNrso2", HistAxis(AXrso));
    Hist * hHCrso2 = Hist::New("hHCrso2", "hHCrso2", HistAxis(AXrso));
    
    Axis AXresx("[um]", 800, -300., 300.);
    Hist * hCKresx = Hist::New("hCKresx", "hCKresx", HistAxis(AXmom, AXresx));
    Hist * hCNresx = Hist::New("hCNresx", "hCNresx", HistAxis(AXmom, AXresx));
    Hist * hHCresx = Hist::New("hHCresx", "hHCresx", HistAxis(AXmom, AXresx));
    
    Axis AXresy("[um]", 800, -300., 300.);
    Hist * hCKresy = Hist::New("hCKresy", "hCKresy", HistAxis(AXmom, AXresy));
    Hist * hCNresy = Hist::New("hCNresy", "hCNresy", HistAxis(AXmom, AXresy));
    Hist * hHCresy = Hist::New("hHCresy", "hHCresy", HistAxis(AXmom, AXresy));

    Long64_t printRate = dst->GetEntries();
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
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
            //if (hit.side != 3) continue; // testcode
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
                if (mchit.mom > 5.) hXres2->fill(resx);
                if (mchit.mom > 5.) hYres2->fill(resy);
               
                Double_t mcos = std::fabs(mchit.dir[2]);
                Double_t etax = (hit.loc[0] - std::lrint(hit.loc[0]));
                Double_t etay = (hit.loc[1] - std::lrint(hit.loc[1]));
                if (hit.side == 3) {
                    if (mchit.mom > 3. && hit.stripSigX.size() != 1) hXresnt->fill(etax, resx);
                    if (mchit.mom > 3. && hit.stripSigY.size() != 1) hYresnt->fill(etay, resy);
                    if (mchit.mom > 3. && hit.stripSigX.size() == 1) hXresn1->fill(etax, resx);
                    if (mchit.mom > 3. && hit.stripSigY.size() == 1) hYresn1->fill(etay, resy);
                    if (mchit.mom > 3. && hit.stripSigX.size() == 2) hXresn2->fill(etax, resx);
                    if (mchit.mom > 3. && hit.stripSigY.size() == 2) hYresn2->fill(etay, resy);
                    if (mchit.mom > 3. && hit.stripSigX.size() == 3) hXresn3->fill(etax, resx);
                    if (mchit.mom > 3. && hit.stripSigY.size() == 3) hYresn3->fill(etay, resy);
                    if (mchit.mom > 3. && hit.stripSigX.size() >= 4) hXresn4->fill(etax, resx);
                    if (mchit.mom > 3. && hit.stripSigY.size() >= 4) hYresn4->fill(etay, resy);
                }

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
        Bool_t succ = tr.fit();
        //Bool_t succ = false;
        //tr.print();
        Double_t hc_irig = tr.part().irig();
       
        Double_t pow27 = std::pow(100., 1.7) * std::pow(mc_mom, -1.7);
        if (mc_mom > 100. && track.status[0][patt]) hCKflux->fill(ck_irig, pow27);
        if (mc_mom > 100. && track.status[1][patt]) hCNflux->fill(cn_irig, pow27);
        if (mc_mom > 100. && succ) hHCflux->fill(hc_irig, pow27);
        
        if (track.status[0][patt]) hCKrso->fill(mc_mom, bincen * (ck_irig - mc_irig));
        if (track.status[1][patt]) hCNrso->fill(mc_mom, bincen * (cn_irig - mc_irig));
        if (succ) hHCrso->fill(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (mc_mom > 100.) {
            if (track.status[0][patt]) hCKrso2->fill(100. * (ck_irig - mc_irig));
            if (track.status[1][patt]) hCNrso2->fill(100. * (cn_irig - mc_irig));
            if (succ) hHCrso2->fill(100. * (hc_irig - mc_irig));
        }
        
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
