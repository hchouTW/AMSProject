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
   
    //MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree();
    //
    //MatFld&& mf1 = MatMgnt::Get(SVecD<3>(0, 0, 175), SVecD<3>(0, 0, 50));
    //MatFld&& mf2 = MatMgnt::Get(SVecD<3>(0, 0, -30), SVecD<3>(0, 0, -100));
    //MatFld&& mf3 = MatMgnt::Get(SVecD<3>(0, 0, 50), SVecD<3>(0, 0, -50));
    //mf1.print();
    //mf2.print();
    //mf3.print();

    //return 0;

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
    //PhyArg::SetOpt(true, false);
    PhyArg::SetOpt(true, true);
    Int_t laySat = 4;
    Int_t layEnd = 5;
    
    TFile * ofle = new TFile(Form("%s/prop_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 100, 0.5, 4000., AxisScale::kLog);

    Axis AXnrl("Nrl", 800, 0.005, 0.02);
    Axis AXela("Ela", 800, 0.1, 0.3);
    Hist * hNrl = Hist::New("hNrl", "hNrl", HistAxis(AXmom, AXnrl)); // (TH2D) Number of Radiator Length [1]
    Hist * hEla = Hist::New("hEla", "hEla", HistAxis(AXmom, AXela)); // (TH2D) Electron Abundance [mol/cm^2]
    
    // Number of Radiation Length Map
    Axis AXcxy("Coo [cm]", 260, -65., 65.);
    Hist * hEvt = Hist::New("hEvt", "hEvt", HistAxis(AXcxy, AXcxy));
    Hist * hMap = Hist::New("hMap", "hMap", HistAxis(AXcxy, AXcxy));

    // Prop
    Axis AXcoo("Residual [cm * p#beta/Q * L^-1]", 2000, -0.70, 0.70);
    Hist * hMcx = Hist::New("hMcx", "hMcx", HistAxis(AXmom, AXcoo)); // (TH2D) MC: residual x
    Hist * hMcy = Hist::New("hMcy", "hMcy", HistAxis(AXmom, AXcoo)); // (TH2D) MC: residual y
    Hist * hTcx = Hist::New("hTcx", "hTcx", HistAxis(AXmom, AXcoo)); // (TH2D) ToyMC: residual x
    Hist * hTcy = Hist::New("hTcy", "hTcy", HistAxis(AXmom, AXcoo)); // (TH2D) ToyMC: residual y
    
    Axis AXagl("Residual [p#beta/Q]", 2000, -0.70, 0.70);
    Hist * hMux = Hist::New("hMux", "hMux", HistAxis(AXmom, AXagl)); // (TH2D) MC: cosine angle x
    Hist * hMuy = Hist::New("hMuy", "hMuy", HistAxis(AXmom, AXagl)); // (TH2D) MC: cosine angle y
    Hist * hTux = Hist::New("hTux", "hTux", HistAxis(AXmom, AXagl)); // (TH2D) ToyMC: cosine angle x
    Hist * hTuy = Hist::New("hTuy", "hTuy", HistAxis(AXmom, AXagl)); // (TH2D) ToyMC: cosine angle y
   
    Axis AXels("Eloss [MeV * #beta^{2}/Q^{2}]", 600, 0.2, 15);
    Hist * hMee = Hist::New("hMee", "hMee", HistAxis(AXmom, AXels)); // (TH2D) MC: kinetic energy difference
    Hist * hTee = Hist::New("hTee", "hTee", HistAxis(AXmom, AXels)); // (TH2D) ToyMC: kinetic energy difference
    
    Axis AXcoo2("Residual [cm * p#beta/Q * L^-1]", 600, -0.10, 0.10);
    Axis AXagl2("Residual [p#beta/Q]", 600, -0.10, 0.10);
    Hist * hMcux = Hist::New("hMcux", "hMcux", HistAxis(AXcoo2, AXagl2)); // (TH2D) MC: residual x vs. cosine angle x
    Hist * hMcuy = Hist::New("hMcuy", "hMcuy", HistAxis(AXcoo2, AXagl2)); // (TH2D) MC: residual y cosine angle y
    Hist * hTcux = Hist::New("hTcux", "hTcux", HistAxis(AXcoo2, AXagl2)); // (TH2D) ToyMC: residual x vs. cosine angle x
    Hist * hTcuy = Hist::New("hTcuy", "hTcuy", HistAxis(AXcoo2, AXagl2)); // (TH2D) ToyMC: residual y vs. cosine angle y

    Long64_t printRate = dst->GetEntries();
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        dst->GetEntry(entry);

        // MC hit from primary particle
        // HitTRKMCInfo from TrMCClusterR (GetGtrkID == AMSEventR->GetPrimaryMC->trkID)
        HitTRKMCInfo * mchitU = nullptr; // MC starting hit (SiTr-layerJ = laySat)
        HitTRKMCInfo * mchitL = nullptr; // MC ending   hit (SiTr-layerJ = layEnd)
        Double_t       kEngU  = 0; // kinetic energy
        Double_t       kEngL  = 0; // kinetic energy
        for (auto&& hit : fG4mc->primPart.hits) {
            Double_t radius = std::sqrt(hit.coo[0]*hit.coo[0] + hit.coo[1]*hit.coo[1]);
            Double_t maxxy  = std::max(std::fabs(hit.coo[0]), std::fabs(hit.coo[1]));
            Double_t keng   = std::sqrt(hit.mom * hit.mom + fG4mc->primPart.mass * fG4mc->primPart.mass) - fG4mc->primPart.mass;
            //if (maxxy > 16.) continue;
            //if (radius > 35.) continue;
            if (hit.layJ == laySat) { mchitU = &hit; kEngU = keng; }
            if (hit.layJ == layEnd) { mchitL = &hit; kEngL = keng; }
        }
        if (mchitU && mchitL) {
            PhySt part(PartType::Proton);
            part.set_state(
                mchitU->coo[0], mchitU->coo[1], mchitU->coo[2],
                mchitU->dir[0], mchitU->dir[1], mchitU->dir[2]
            );
            part.set_mom(mchitU->mom);
            Double_t mc_mom = part.mom();
            
            MatFld mfld;       // Material information
            PhySt ppst(part);  // Particle Status
            PropMgnt::PropToZ(mchitL->coo[2], ppst, &mfld); // Propagate to Z with magnetic field
            Double_t len = std::fabs(mchitL->coo[2]-mchitU->coo[2]); // Delta Z
            Double_t nrl = mfld.nrl();  // Number of Radiator Length [1]
            Double_t ela = mfld.ela();  // Electron Abundance [mol/cm^2]
            SVecD<3> refc = ppst.c();   // coord
            SVecD<3> refu = ppst.u();   // cosine angle
    
            const Double_t GeV2MeV = 1.0e3;
            Double_t scl_eloss = (part.bta() * part.bta()) / (part.chrg() * part.chrg()) / ela;       // normalized factor (energy loss)
            Double_t scl_mscat = (part.mom() * part.bta()) / std::fabs(part.chrg()) / std::sqrt(nrl); // normalized factor (multiple-scattering)
            
            hNrl->fill(mc_mom, nrl);
            hEla->fill(mc_mom, ela);
            
            ppst = part;
            PropMgnt::PropToZWithMC(mchitL->coo[2], ppst);
            Double_t mc_resc[2] = { mchitL->coo[0] - refc(0), mchitL->coo[1] - refc(1) }; // MC: residual xy [cm]
            Double_t mc_resu[2] = { mchitL->dir[0] - refu(0), mchitL->dir[1] - refu(1) }; // MC: cosine angle xy [1]
            Double_t mc_elsm    = GeV2MeV * (kEngU - kEngL);                              // MC: kinetic energy difference [GeV]
            Double_t tm_resc[2] = { ppst.cx() - refc(0), ppst.cy() - refc(1) };           // ToyMC: residual xy [cm]
            Double_t tm_resu[2] = { ppst.ux() - refu(0), ppst.uy() - refu(1) };           // ToyMC: cosine angle xy [1]
            Double_t tm_elsm    = GeV2MeV * (part.ke() - ppst.ke());                      // ToyMC: kinetic energy difference [GeV]

            hMcx->fill(mc_mom, scl_mscat * mc_resc[0] / len);
            hMcy->fill(mc_mom, scl_mscat * mc_resc[1] / len);
            hMux->fill(mc_mom, scl_mscat * mc_resu[0]);
            hMuy->fill(mc_mom, scl_mscat * mc_resu[1]);
            hMee->fill(mc_mom, scl_eloss * mc_elsm);
            hTcx->fill(mc_mom, scl_mscat * tm_resc[0] / len);
            hTcy->fill(mc_mom, scl_mscat * tm_resc[1] / len);
            hTux->fill(mc_mom, scl_mscat * tm_resu[0]);
            hTuy->fill(mc_mom, scl_mscat * tm_resu[1]);
            hTee->fill(mc_mom, scl_eloss * tm_elsm);
            
            if (mc_mom > 5.0 && mc_mom < 100.) hMcux->fill(scl_mscat * mc_resc[0] / len, scl_mscat * mc_resu[0]);
            if (mc_mom > 5.0 && mc_mom < 100.) hMcuy->fill(scl_mscat * mc_resc[1] / len, scl_mscat * mc_resu[1]);
            if (mc_mom > 5.0 && mc_mom < 100.) hTcux->fill(scl_mscat * tm_resc[0] / len, scl_mscat * tm_resu[0]);
            if (mc_mom > 5.0 && mc_mom < 100.) hTcuy->fill(scl_mscat * tm_resc[1] / len, scl_mscat * tm_resu[1]);

            Double_t cx = 0.5 * (mchitU->coo[0] + mchitL->coo[0]);
            Double_t cy = 0.5 * (mchitU->coo[1] + mchitL->coo[1]);
            hEvt->fill(cx, cy);
            hMap->fill(cx, cy, nrl);
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
