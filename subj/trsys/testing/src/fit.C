//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/data1/hchou/17SEP27/src/ClassDef.h"
#include "/data1/hchou/17SEP27/src/ClassDef.C"

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    TH1::AddDirectory(true);
    
    MGConfig::JobOpt opt(argc, argv);

    //MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree();
    MatFld&& mat = MatMgnt::Get(SVecD<3>(0, 0, 195.0), SVecD<3>(0, 0, 50.0));
    mat.print();

    ////MatFld&& matNaf = MatMgnt::Get(SVecD<3>(1, 1, -71.8), SVecD<3>(1, 1, -74.8));
    ////MatFld&& matAgl = MatMgnt::Get(SVecD<3>(20, 20, -71.8), SVecD<3>(20, 20, -74.8));
    //MatFld&& matNaf = MatMgnt::Get(SVecD<3>(1, 1, -73), SVecD<3>(1, 1, -76));
    //MatFld&& matAgl = MatMgnt::Get(SVecD<3>(20, 20, -73), SVecD<3>(20, 20, -76));
    //matNaf.print();
    //matAgl.print();

    //return false;


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
    //RICH * fRich = new RICH;
    //ECAL * fEcal = new ECAL;

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
    
    TFile * ofle = new TFile("fit.root", "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 10, 0.5, 10., AxisScale::kLog);

    // Hit
    Axis AXres("Residual [10^{-4} cm]", 800, -300., 300.);
    Hist * hXres = Hist::New("hXres", "hXres", HistAxis(AXmom, AXres));
    Hist * hYres = Hist::New("hYres", "hYres", HistAxis(AXmom, AXres));

    // Prop
    Axis AXcoo("Residual [cm * p#beta/Q^{2}]", 400, -3, 3);
    Hist * hMCXcoo = Hist::New("hMCXcoo", "hMCXcoo", HistAxis(AXmom, AXcoo));
    Hist * hMCYcoo = Hist::New("hMCYcoo", "hMCYcoo", HistAxis(AXmom, AXcoo));
    Hist * hTMXcoo = Hist::New("hTMXcoo", "hTMXcoo", HistAxis(AXmom, AXcoo));
    Hist * hTMYcoo = Hist::New("hTMYcoo", "hTMYcoo", HistAxis(AXmom, AXcoo));
    
    Axis AXagl("Residual [p#beta/Q^{2}]", 400, -0.05, 0.05);
    Hist * hMCXagl = Hist::New("hMCXagl", "hMCXagl", HistAxis(AXmom, AXagl));
    Hist * hMCYagl = Hist::New("hMCYagl", "hMCYagl", HistAxis(AXmom, AXagl));
    Hist * hTMXagl = Hist::New("hTMXagl", "hTMXagl", HistAxis(AXmom, AXagl));
    Hist * hTMYagl = Hist::New("hTMYagl", "hTMYagl", HistAxis(AXmom, AXagl));
    
    Axis AXels("Eloss [GeV * #beta^{2}/Q^{2}]", 400, 0.005, 0.05);
    Hist * hMCelsm = Hist::New("hMCelsm", "hMCelsm", HistAxis(AXmom, AXels));
    Hist * hTMelsm = Hist::New("hTMelsm", "hTMelsm", HistAxis(AXmom, AXels));
    
    Hist * hMCXcu = Hist::New("hMCXcu", "hMCXcu", HistAxis(AXcoo, AXagl));
    Hist * hMCYcu = Hist::New("hMCYcu", "hMCYcu", HistAxis(AXcoo, AXagl));
    Hist * hTMXcu = Hist::New("hTMXcu", "hTMXcu", HistAxis(AXcoo, AXagl));
    Hist * hTMYcu = Hist::New("hTMYcu", "hTMYcu", HistAxis(AXcoo, AXagl));

    // Fit
    Axis AXrso("(1/Rm - 1/Rt) [1/GV]", 800, -1.5, 1.5);
    Hist * hCKrso = Hist::New("hCKrso", "hCKrso", HistAxis(AXmom, AXrso));
    Hist * hCNrso = Hist::New("hCNrso", "hCNrso", HistAxis(AXmom, AXrso));
    Hist * hHYrso = Hist::New("hHYrso", "hHYrso", HistAxis(AXmom, AXrso));

    Long64_t printRate = dst->GetEntries()/50;
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        //if (entry > 100) break; // testcode
        dst->GetEntry(entry);
        
        Double_t mc_mom  = (fG4mc->primPart.mom);
        Double_t bincen = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
        //if (mc_mom > 0.8) continue;
        //COUT("=== MOM MC %14.8f ===\n", mc_mom);
       
        if (fTof->betaH < 0.3 || fTof->betaH > 1.3) continue;
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        if (fTof->Qall < 0.8 || fTof->Qall > 1.8) continue;
        if (fTof->betaHPatt != 15) continue;

        if (fTrk->tracks.size() != 1) continue;
        TrackInfo& track = fTrk->tracks.at(0);

        Short_t countMC = 0;
        for (auto&& hit : fG4mc->primPart.hits) {
            if (hit.layJ < 2 || hit.layJ > 8) continue;
            countMC++;
        }
        
        Short_t countMS = 0;
        for (auto&& hit : track.hits) {
            if (hit.layJ < 2 || hit.layJ > 8) continue;
            countMS++;
        }

        Short_t cutNHit = 4 + optL1 + optL9;
        if (countMC <= cutNHit || countMS <= cutNHit) continue;

        // Propagation Testing
        HitTRKMCInfo * mchitL1 = nullptr;
        HitTRKMCInfo * mchitL2 = nullptr;
        for (auto&& hit : fG4mc->primPart.hits) {
            if (hit.layJ == 1) mchitL1 = &hit;
            if (hit.layJ == 2) mchitL2 = &hit;
        }
        if (mchitL1 && mchitL2) {
            PhySt part(PartType::Proton);
            part.set_state(
                mchitL1->coo[0], mchitL1->coo[1], mchitL1->coo[2],
                mchitL1->dir[0], mchitL1->dir[1], mchitL1->dir[2]
            );
            part.set_mom(mchitL1->mom);
            //Double_t mc_mom  = mchitL1->mom;
            Double_t scl_mscat = (part.mom() * part.bta() / (part.chrg() * part.chrg()));
            Double_t scl_eloss = (part.bta() * part.bta() / (part.chrg() * part.chrg()));
            
            PhySt ppst(part);
            PropMgnt::PropToZ(mchitL2->coo[2], ppst);
            SVecD<3> refc = ppst.c();
            SVecD<3> refu = ppst.u();
            
            ppst = part;
            PropMgnt::PropToZWithMC(mchitL2->coo[2], ppst);
            Double_t mc_resc[2] = { mchitL2->coo[0] - refc(0), mchitL2->coo[1] - refc(1) };
            Double_t mc_resu[2] = { mchitL2->dir[0] - refu(0), mchitL2->dir[1] - refu(1) };
            Double_t mc_elsm    = (part.mom() - mchitL2->mom);
            Double_t tm_resc[2] = { ppst.cx() - refc(0), ppst.cy() - refc(1) };
            Double_t tm_resu[2] = { ppst.ux() - refu(0), ppst.uy() - refu(1) };
            Double_t tm_elsm    = (part.mom() - ppst.mom());
            hMCXcoo->fill(mc_mom, scl_mscat * mc_resc[0]);
            hMCYcoo->fill(mc_mom, scl_mscat * mc_resc[1]);
            hMCXagl->fill(mc_mom, scl_mscat * mc_resu[0]);
            hMCYagl->fill(mc_mom, scl_mscat * mc_resu[1]);
            hMCelsm->fill(mc_mom, scl_eloss * mc_elsm);
            hTMXcoo->fill(mc_mom, scl_mscat * tm_resc[0]);
            hTMYcoo->fill(mc_mom, scl_mscat * tm_resc[1]);
            hTMXagl->fill(mc_mom, scl_mscat * tm_resu[0]);
            hTMYagl->fill(mc_mom, scl_mscat * tm_resu[1]);
            hTMelsm->fill(mc_mom, scl_eloss * tm_elsm);
            
            hMCXcu->fill(scl_mscat * mc_resc[0], scl_mscat * mc_resu[0]);
            hMCYcu->fill(scl_mscat * mc_resc[1], scl_mscat * mc_resu[1]);
            hTMXcu->fill(scl_mscat * tm_resc[0], scl_mscat * tm_resu[0]);
            hTMYcu->fill(scl_mscat * tm_resc[1], scl_mscat * tm_resu[1]);
        }
        continue;





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
        if (mhits.size() <= cutNHit) continue;

        if (optL1 && !hasL1) continue;
        //if (optL9 && !hasL9) continue;

        Double_t mc_irig = (fG4mc->primPart.chrg / fG4mc->primPart.mom);
        Double_t ck_irig = ((track.status[0][0]) ? 1.0/track.rigidity[0][0] : 0.);
        Double_t cn_irig = ((track.status[1][0]) ? 1.0/track.rigidity[1][0] : 0.);

        // Fitting by H.Y.Chou
        //std::cout << Form("\n==== Entry %lld ====\n", entry);

        PhyTr tr(mhits);
        tr.fit();
        //tr.print();
        Double_t hy_irig = tr.part().irig();

        hCKrso->fill(mc_mom, bincen * (ck_irig - mc_irig));
        hCNrso->fill(mc_mom, bincen * (cn_irig - mc_irig));
        hHYrso->fill(mc_mom, bincen * (hy_irig - mc_irig));
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
    //if (fRich) { delete fRich; fRich = nullptr; }
    //if (fEcal) { delete fEcal; fEcal = nullptr; }

    return 0;
}
