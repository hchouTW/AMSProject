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
   
    //MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree();
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
    PhyArg::SetOpt(true, true);
    Int_t layBeg = 1;
    Int_t layEnd = 2;
    
    TFile * ofle = new TFile(Form("%s/prop_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 40, 0.5, 10., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 200, 1.0, 800., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 200, 0.5, 800., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 50, 20., 16000., AxisScale::kLog);
    Axis AXcos("Cos [1]", 40, 0.9,  1.);

    Hist * hCos = Hist::New("hCos", "hCos", HistAxis(AXcos, AXcos));

    // Hit
    Axis AXres("Residual [10^{-4} cm]", 800, -300., 300.);
    Hist * hXres = Hist::New("hXres", "hXres", HistAxis(AXmom, AXres));
    Hist * hYres = Hist::New("hYres", "hYres", HistAxis(AXmom, AXres));
    
    Axis AXnrl("Nrl", 800, 0., 0.05);
    Axis AXelc("Elc", 800, 0., 1.00);
    Axis AXchg("Chg", 800, 15., 17.);
    Hist * hNrl = Hist::New("hNrl", "hNrl", HistAxis(AXnrl));
    Hist * hElc = Hist::New("hElc", "hElc", HistAxis(AXelc));
    Hist * hChg = Hist::New("hChg", "hChg", HistAxis(AXchg));
    Hist * hNrlElc = Hist::New("hNrlElc", "hNrlElc", HistAxis(AXnrl, AXelc));

    // Prop
    Axis AXcoo("Residual [cm * p#beta/Q^{2}]", 400, -4.0, 4.0);
    Hist * hMcx = Hist::New("hMcx", "hMcx", HistAxis(AXmom, AXcoo));
    Hist * hMcy = Hist::New("hMcy", "hMcy", HistAxis(AXmom, AXcoo));
    Hist * hTcx = Hist::New("hTcx", "hTcx", HistAxis(AXmom, AXcoo));
    Hist * hTcy = Hist::New("hTcy", "hTcy", HistAxis(AXmom, AXcoo));
    
    Axis AXagl("Residual [p#beta/Q^{2}]", 400, -0.05, 0.05);
    Hist * hMux = Hist::New("hMux", "hMux", HistAxis(AXmom, AXagl));
    Hist * hMuy = Hist::New("hMuy", "hMuy", HistAxis(AXmom, AXagl));
    Hist * hTux = Hist::New("hTux", "hTux", HistAxis(AXmom, AXagl));
    Hist * hTuy = Hist::New("hTuy", "hTuy", HistAxis(AXmom, AXagl));
   
    Axis AXels("Eloss [GeV * #beta^{2}/Q^{2}]", 400, 0.0015, 0.015);
    Hist * hMee = Hist::New("hMee", "hMee", HistAxis(AXmom, AXels));
    Hist * hTee = Hist::New("hTee", "hTee", HistAxis(AXmom, AXels));
    
    Hist * hMcux = Hist::New("hMcux", "hMcux", HistAxis(AXcoo, AXagl));
    Hist * hMcuy = Hist::New("hMcuy", "hMcuy", HistAxis(AXcoo, AXagl));
    Hist * hTcux = Hist::New("hTcux", "hTcux", HistAxis(AXcoo, AXagl));
    Hist * hTcuy = Hist::New("hTcuy", "hTcuy", HistAxis(AXcoo, AXagl));

    Long64_t printRate = dst->GetEntries()/40;
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        //if (entry>(dst->GetEntries()/20)) break; // testcode
        //if (entry%10!=0) continue; // testcode
        dst->GetEntry(entry);

        // Propagation Testing
        HitTRKMCInfo * mchitU = nullptr;
        HitTRKMCInfo * mchitL = nullptr;
        for (auto&& hit : fG4mc->primPart.hits) {
            Double_t radius = std::sqrt(hit.coo[0]*hit.coo[0] + hit.coo[1]*hit.coo[1]);
            Double_t maxxy  = std::max(std::fabs(hit.coo[0]), std::fabs(hit.coo[1]));
            Double_t cos    = std::fabs(hit.dir[2]);
            //if (maxxy > 30.) continue;
            //if (maxxy > 18.) continue;
            if (radius > 40.) continue;
            if (hit.layJ == layBeg) mchitU = &hit;
            if (hit.layJ == layEnd) mchitL = &hit;
        }
        if (mchitU && mchitL) {
            PhySt part(PartType::Proton);
            part.set_state(
                mchitU->coo[0], mchitU->coo[1], mchitU->coo[2],
                mchitU->dir[0], mchitU->dir[1], mchitU->dir[2]
            );
            part.set_mom(mchitU->mom);
            
            Double_t mc_mom = std::sqrt(mchitU->mom * mchitL->mom);
            Double_t mc_cos = std::sqrt(std::fabs(mchitU->dir[2] * mchitL->dir[2]));
            
            part.set_mom(mc_mom);
            Double_t scl_eloss = (part.bta() * part.bta()) / (part.chrg() * part.chrg());
            Double_t scl_mscat = (part.mom() * part.bta()) / std::fabs(part.chrg());
            part.set_mom(mchitU->mom);
            
            hCos->fill(std::fabs(mchitU->dir[2]), std::fabs(mchitL->dir[2]));

            PhySt ppst(part);
            
            MatFld mfld;
            PropMgnt::PropToZ(mchitL->coo[2], ppst, &mfld);
            
            //PropMgnt::PropToZ_AMSLibs(mchitL->coo[2], ppst);
            //MatFld&& mfld = MatMgnt::Get(SVecD<3>(mchitU->coo[0], mchitU->coo[1], mchitU->coo[2]), SVecD<3>(mchitL->coo[0], mchitL->coo[1], mchitL->coo[2]));
            
            SVecD<3> refc = ppst.c();
            SVecD<3> refu = ppst.u();
            
            Double_t nrl = mfld.num_rad_len();
            Double_t elc = mfld.elcloud_abundance();
           
            scl_mscat /= std::sqrt(nrl);
            scl_eloss /= elc;

            hNrl->fill(nrl);
            hElc->fill(elc);
            hChg->fill(elc/nrl);
            hNrlElc->fill(nrl, elc);
            
            ppst = part;
            PropMgnt::PropToZWithMC(mchitL->coo[2], ppst);
            Double_t mc_resc[2] = { mchitL->coo[0] - refc(0), mchitL->coo[1] - refc(1) };
            Double_t mc_resu[2] = { mchitL->dir[0] - refu(0), mchitL->dir[1] - refu(1) };
            Double_t mc_elsm    = (mchitU->mom - mchitL->mom);
            Double_t tm_resc[2] = { ppst.cx() - refc(0), ppst.cy() - refc(1) };
            Double_t tm_resu[2] = { ppst.ux() - refu(0), ppst.uy() - refu(1) };
            Double_t tm_elsm    = (mchitU->mom - ppst.mom());
            hMcx->fill(mc_mom, scl_mscat * mc_resc[0]);
            hMcy->fill(mc_mom, scl_mscat * mc_resc[1]);
            hMux->fill(mc_mom, scl_mscat * mc_resu[0]);
            hMuy->fill(mc_mom, scl_mscat * mc_resu[1]);
            hMee->fill(mc_mom, scl_eloss * mc_elsm);
            hTcx->fill(mc_mom, scl_mscat * tm_resc[0]);
            hTcy->fill(mc_mom, scl_mscat * tm_resc[1]);
            hTux->fill(mc_mom, scl_mscat * tm_resu[0]);
            hTuy->fill(mc_mom, scl_mscat * tm_resu[1]);
            hTee->fill(mc_mom, scl_eloss * tm_elsm);
            
            hMcux->fill(scl_mscat * mc_resc[0], scl_mscat * mc_resu[0]);
            hMcuy->fill(scl_mscat * mc_resc[1], scl_mscat * mc_resu[1]);
            hTcux->fill(scl_mscat * tm_resc[0], scl_mscat * tm_resu[0]);
            hTcuy->fill(scl_mscat * tm_resc[1], scl_mscat * tm_resu[1]);
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
