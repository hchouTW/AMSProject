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

    MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree();
    MatFld&& mat1 = MatMgnt::Get(SVecD<3>(0, 0, 195.0), SVecD<3>(0, 0, 50.0));
    MatFld&& mat2 = MatMgnt::Get(SVecD<3>(2, 2, -60.0), SVecD<3>(2, 2, -70.0));
    MatFld&& mat3 = MatMgnt::Get(SVecD<3>(0, 0, 58.0),  SVecD<3>(0, 0, 51.0));
    mat1.print();
    mat2.print();
    mat3.print();
    return 1;

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
    PhyArg::SetOpt(true, false);
    Int_t layBeg = 1;
    Int_t layEnd = 2;
    
    TFile * ofle = new TFile("prop_ams02.root", "RECREATE");
    
    //Axis AXmom("Momentum [GeV]", 40, 0.5, 10., AxisScale::kLog);
    Axis AXmom("Momentum [GeV]", 40, 1.0, 800., AxisScale::kLog);
    Axis AXcos("Cos [1]", 40, 0.9,  1.);

    Hist * hCos = Hist::New("hCos", "hCos", HistAxis(AXcos, AXcos));

    // Hit
    Axis AXres("Residual [10^{-4} cm]", 800, -300., 300.);
    Hist * hXres = Hist::New("hXres", "hXres", HistAxis(AXmom, AXres));
    Hist * hYres = Hist::New("hYres", "hYres", HistAxis(AXmom, AXres));

    // Prop
    Axis AXcoo("Residual [cm * p#beta/Q^{2}]", 400, -3, 3);
    Hist * hMcx = Hist::New("hMcx", "hMcx", HistAxis(AXmom, AXcoo));
    Hist * hMcy = Hist::New("hMcy", "hMcy", HistAxis(AXmom, AXcoo));
    Hist * hTcx = Hist::New("hTcx", "hTcx", HistAxis(AXmom, AXcoo));
    Hist * hTcy = Hist::New("hTcy", "hTcy", HistAxis(AXmom, AXcoo));
    
    Axis AXagl("Residual [p#beta/Q^{2}]", 400, -0.05, 0.05);
    Hist * hMux = Hist::New("hMux", "hMux", HistAxis(AXmom, AXagl));
    Hist * hMuy = Hist::New("hMuy", "hMuy", HistAxis(AXmom, AXagl));
    Hist * hTux = Hist::New("hTux", "hTux", HistAxis(AXmom, AXagl));
    Hist * hTuy = Hist::New("hTuy", "hTuy", HistAxis(AXmom, AXagl));
    
    Axis AXels("Eloss [GeV * #beta^{2}/Q^{2}]", 400, 0.005, 0.05);
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
        //if (entry > 100) break; // testcode
        dst->GetEntry(entry);
        
        //Double_t mc_mom  = (fG4mc->primPart.mom);
        //Double_t bincen = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
        //if (mc_mom > 0.8) continue;
        //COUT("=== MOM MC %14.8f ===\n", mc_mom);
       
        //if (fTof->betaH < 0.3 || fTof->betaH > 1.3) continue;
        //if (fTof->normChisqT > 10.) continue;
        //if (fTof->normChisqC > 10.) continue;
        //if (fTof->Qall < 0.8 || fTof->Qall > 1.8) continue;
        //if (fTof->betaHPatt != 15) continue;

        //if (fTrk->tracks.size() != 1) continue;
        //TrackInfo& track = fTrk->tracks.at(0);

        // Propagation Testing
        HitTRKMCInfo * mchitU = nullptr;
        HitTRKMCInfo * mchitL = nullptr;
        for (auto&& hit : fG4mc->primPart.hits) {
            Double_t radius = std::sqrt(hit.coo[0]*hit.coo[0] + hit.coo[1]*hit.coo[1]);
            Double_t cos = std::fabs(hit.dir[2]);
            if (radius > 30.) continue;
            //if (cos < 0.90) continue;
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
            Double_t mc_mom = mchitU->mom;
            Double_t mc_cos = std::fabs(mchitU->dir[2]);
            
            Double_t scl_eloss = mc_cos * (part.bta() * part.bta()) / (part.chrg() * part.chrg());
            Double_t scl_mscat = mc_cos * (part.mom() * part.bta()) / (part.chrg() * part.chrg());
            
            hCos->fill(std::fabs(mchitU->dir[2]), std::fabs(mchitL->dir[2]));
            
            PhySt ppst(part);
            PropMgnt::PropToZ(mchitL->coo[2], ppst);
            SVecD<3> refc = ppst.c();
            SVecD<3> refu = ppst.u();
            
            ppst = part;
            PropMgnt::PropToZWithMC(mchitL->coo[2], ppst);
            Double_t mc_resc[2] = { mchitL->coo[0] - refc(0), mchitL->coo[1] - refc(1) };
            Double_t mc_resu[2] = { mchitL->dir[0] - refu(0), mchitL->dir[1] - refu(1) };
            Double_t mc_elsm    = (part.mom() - mchitL->mom);
            Double_t tm_resc[2] = { ppst.cx() - refc(0), ppst.cy() - refc(1) };
            Double_t tm_resu[2] = { ppst.ux() - refu(0), ppst.uy() - refu(1) };
            Double_t tm_elsm    = (part.mom() - ppst.mom());
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

    TGraphErrors* gMcx_men = new TGraphErrors(); gMcx_men->SetNameTitle("gMcx_men", "");
    TGraphErrors* gMcx_sgm = new TGraphErrors(); gMcx_sgm->SetNameTitle("gMcx_sgm", "");
    TGraphErrors* gMcy_men = new TGraphErrors(); gMcy_men->SetNameTitle("gMcy_men", "");
    TGraphErrors* gMcy_sgm = new TGraphErrors(); gMcy_sgm->SetNameTitle("gMcy_sgm", "");
    TGraphErrors* gMux_men = new TGraphErrors(); gMux_men->SetNameTitle("gMux_men", "");
    TGraphErrors* gMux_sgm = new TGraphErrors(); gMux_sgm->SetNameTitle("gMux_sgm", "");
    TGraphErrors* gMuy_men = new TGraphErrors(); gMuy_men->SetNameTitle("gMuy_men", "");
    TGraphErrors* gMuy_sgm = new TGraphErrors(); gMuy_sgm->SetNameTitle("gMuy_sgm", "");
    TGraphErrors* gTcx_men = new TGraphErrors(); gTcx_men->SetNameTitle("gTcx_men", "");
    TGraphErrors* gTcx_sgm = new TGraphErrors(); gTcx_sgm->SetNameTitle("gTcx_sgm", "");
    TGraphErrors* gTcy_men = new TGraphErrors(); gTcy_men->SetNameTitle("gTcy_men", "");
    TGraphErrors* gTcy_sgm = new TGraphErrors(); gTcy_sgm->SetNameTitle("gTcy_sgm", "");
    TGraphErrors* gTux_men = new TGraphErrors(); gTux_men->SetNameTitle("gTux_men", "");
    TGraphErrors* gTux_sgm = new TGraphErrors(); gTux_sgm->SetNameTitle("gTux_sgm", "");
    TGraphErrors* gTuy_men = new TGraphErrors(); gTuy_men->SetNameTitle("gTuy_men", "");
    TGraphErrors* gTuy_sgm = new TGraphErrors(); gTuy_sgm->SetNameTitle("gTuy_sgm", "");
    
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    std::vector<Hist*> vhMcx = Hist::ProjectAll(HistProj::kY, hMcx);
    std::vector<Hist*> vhMcy = Hist::ProjectAll(HistProj::kY, hMcy);
    std::vector<Hist*> vhMux = Hist::ProjectAll(HistProj::kY, hMux);
    std::vector<Hist*> vhMuy = Hist::ProjectAll(HistProj::kY, hMuy);
    std::vector<Hist*> vhTcx = Hist::ProjectAll(HistProj::kY, hTcx);
    std::vector<Hist*> vhTcy = Hist::ProjectAll(HistProj::kY, hTcy);
    std::vector<Hist*> vhTux = Hist::ProjectAll(HistProj::kY, hTux);
    std::vector<Hist*> vhTuy = Hist::ProjectAll(HistProj::kY, hTuy);
    for (int it = 1; it <= AXmom.nbin(); ++it) {
        double mom = AXmom.center(it+1, AxisScale::kLog);
        PhySt part(PartType::Proton);
        part.set_mom(mom);
        Double_t bta = part.bta();
        Double_t val = mom/part.mass();
        //Double_t val = bta;
        
        gaus->SetParameters(1000, 0, (*vhMcx.at(it))()->GetRMS());
        (*vhMcx.at(it))()->Fit(gaus, "q0", "");
        (*vhMcx.at(it))()->Fit(gaus, "q0", "");
        gMcx_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gMcx_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gMcx_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gMcx_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhMcy.at(it))()->GetRMS());
        (*vhMcy.at(it))()->Fit(gaus, "q0", "");
        (*vhMcy.at(it))()->Fit(gaus, "q0", "");
        gMcy_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gMcy_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gMcy_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gMcy_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhMux.at(it))()->GetRMS());
        (*vhMux.at(it))()->Fit(gaus, "q0", "");
        (*vhMux.at(it))()->Fit(gaus, "q0", "");
        gMux_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gMux_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gMux_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gMux_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhMuy.at(it))()->GetRMS());
        (*vhMuy.at(it))()->Fit(gaus, "q0", "");
        (*vhMuy.at(it))()->Fit(gaus, "q0", "");
        gMuy_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gMuy_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gMuy_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gMuy_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhTcx.at(it))()->GetRMS());
        (*vhTcx.at(it))()->Fit(gaus, "q0", "");
        (*vhTcx.at(it))()->Fit(gaus, "q0", "");
        gTcx_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gTcx_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gTcx_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gTcx_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhTcy.at(it))()->GetRMS());
        (*vhTcy.at(it))()->Fit(gaus, "q0", "");
        (*vhTcy.at(it))()->Fit(gaus, "q0", "");
        gTcy_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gTcy_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gTcy_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gTcy_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhTux.at(it))()->GetRMS());
        (*vhTux.at(it))()->Fit(gaus, "q0", "");
        (*vhTux.at(it))()->Fit(gaus, "q0", "");
        gTux_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gTux_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gTux_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gTux_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhTuy.at(it))()->GetRMS());
        (*vhTuy.at(it))()->Fit(gaus, "q0", "");
        (*vhTuy.at(it))()->Fit(gaus, "q0", "");
        gTuy_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gTuy_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gTuy_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gTuy_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
    } 

    gMcx_men->Write();
    gMcx_sgm->Write();
    gMcy_men->Write();
    gMcy_sgm->Write();
    gMux_men->Write();
    gMux_sgm->Write();
    gMuy_men->Write();
    gMuy_sgm->Write();
    gTcx_men->Write();
    gTcx_sgm->Write();
    gTcy_men->Write();
    gTcy_sgm->Write();
    gTux_men->Write();
    gTux_sgm->Write();
    gTuy_men->Write();
    gTuy_sgm->Write();

    TGraphErrors* gMee_kpa = new TGraphErrors(); gMee_kpa->SetNameTitle("gMee_kpa",  "");
    TGraphErrors* gMee_mpv = new TGraphErrors(); gMee_mpv->SetNameTitle("gMee_mpv",  "");
    TGraphErrors* gMee_sgm = new TGraphErrors(); gMee_sgm->SetNameTitle("gMee_sgm",  "");
    TGraphErrors* gTee_kpa = new TGraphErrors(); gTee_kpa->SetNameTitle("gTee_kpa",  "");
    TGraphErrors* gTee_mpv = new TGraphErrors(); gTee_mpv->SetNameTitle("gTee_mpv",  "");
    TGraphErrors* gTee_sgm = new TGraphErrors(); gTee_sgm->SetNameTitle("gTee_sgm",  "");

    TF1 * feloss = new TF1("feloss", "[0] * TMath::Power( ([2]/x)/[1]/[1], ([2]/x)/[1]/[1] ) / TMath::Gamma( ([2]/x)/[1]/[1] ) * TMath::Exp(-(([2]/x)/[1]/[1]) * ((x-[2])/[3] + TMath::Exp(-(x-[2])/[3])) )");

    std::vector<Hist*> vhMee = Hist::ProjectAll(HistProj::kY, hMee);
    std::vector<Hist*> vhTee = Hist::ProjectAll(HistProj::kY, hTee);
    for (int it = 1; it <= AXmom.nbin(); ++it) {
        double mom = AXmom.center(it+1, AxisScale::kLog);
        PhySt part(PartType::Proton);
        part.set_mom(mom);
        Double_t bta = part.bta();
        Double_t val = mom/part.mass();
        //Double_t val = bta;
        
        double mpeak = (*vhMee.at(it))()->GetXaxis()->GetBinCenter((*vhMee.at(it))()->GetMaximumBin());
        feloss->SetParameters(1000., 1.0, mpeak, 0.1*mpeak, 1.0);
        (*vhMee.at(it))()->Fit(feloss, "q0", "", 0.8*mpeak, 3*mpeak);
        (*vhMee.at(it))()->Fit(feloss, "q0", "", 0.8*mpeak, 3*mpeak);
        (*vhMee.at(it))()->Fit(feloss, "q0", "", 0.8*mpeak, 3*mpeak);
        gMee_kpa->SetPoint     (it-1, val, feloss->GetParameter(1));
        gMee_kpa->SetPointError(it-1,  0., feloss->GetParError(1));
        gMee_mpv->SetPoint     (it-1, val, feloss->GetParameter(2));
        gMee_mpv->SetPointError(it-1,  0., feloss->GetParError(2));
        gMee_sgm->SetPoint     (it-1, val, feloss->GetParameter(3));
        gMee_sgm->SetPointError(it-1,  0., feloss->GetParError(3));
        
        double tpeak = (*vhTee.at(it))()->GetXaxis()->GetBinCenter((*vhTee.at(it))()->GetMaximumBin());
        feloss->SetParameters(1000., 1.0, tpeak, 0.1*tpeak, 1.0);
        (*vhTee.at(it))()->Fit(feloss, "q0", "", 0.8*tpeak, 3*tpeak);
        (*vhTee.at(it))()->Fit(feloss, "q0", "", 0.8*tpeak, 3*tpeak);
        (*vhTee.at(it))()->Fit(feloss, "q0", "", 0.8*tpeak, 3*tpeak);
        gTee_kpa->SetPoint     (it-1, val, feloss->GetParameter(1));
        gTee_kpa->SetPointError(it-1,  0., feloss->GetParError(1));
        gTee_mpv->SetPoint     (it-1, val, feloss->GetParameter(2));
        gTee_mpv->SetPointError(it-1,  0., feloss->GetParError(2));
        gTee_sgm->SetPoint     (it-1, val, feloss->GetParameter(3));
        gTee_sgm->SetPointError(it-1,  0., feloss->GetParError(3));
    }

    gMee_kpa->Write();
    gMee_mpv->Write();
    gMee_sgm->Write();
    gTee_kpa->Write();
    gTee_mpv->Write();
    gTee_sgm->Write();
    
    TGraphErrors* gMTcx = new TGraphErrors(); gMTcx->SetNameTitle("gMTcx", "");
    TGraphErrors* gMTcy = new TGraphErrors(); gMTcy->SetNameTitle("gMTcy", "");
    TGraphErrors* gMTux = new TGraphErrors(); gMTux->SetNameTitle("gMTux", "");
    TGraphErrors* gMTuy = new TGraphErrors(); gMTuy->SetNameTitle("gMTuy", "");
    TGraphErrors* gMTkpa = new TGraphErrors(); gMTkpa->SetNameTitle("gMTkpa", "");
    TGraphErrors* gMTmpv = new TGraphErrors(); gMTmpv->SetNameTitle("gMTmpv", "");
    TGraphErrors* gMTsgm = new TGraphErrors(); gMTsgm->SetNameTitle("gMTsgm", "");
    for (int it = 0; it < AXmom.nbin(); ++it) {
        gMTcx->SetPoint(it, gMcx_sgm->GetX()[it], gMcx_sgm->GetY()[it]/gTcx_sgm->GetY()[it]);
        gMTcy->SetPoint(it, gMcy_sgm->GetX()[it], gMcy_sgm->GetY()[it]/gTcy_sgm->GetY()[it]);
        gMTux->SetPoint(it, gMux_sgm->GetX()[it], gMux_sgm->GetY()[it]/gTux_sgm->GetY()[it]);
        gMTuy->SetPoint(it, gMuy_sgm->GetX()[it], gMuy_sgm->GetY()[it]/gTuy_sgm->GetY()[it]);
        gMTkpa->SetPoint(it, gMee_kpa->GetX()[it], gMee_kpa->GetY()[it]/gTee_kpa->GetY()[it]);
        gMTmpv->SetPoint(it, gMee_mpv->GetX()[it], gMee_mpv->GetY()[it]/gTee_mpv->GetY()[it]);
        gMTsgm->SetPoint(it, gMee_sgm->GetX()[it], gMee_sgm->GetY()[it]/gTee_sgm->GetY()[it]);
    }
    
    gMTcx ->Write(); 
    gMTcy ->Write(); 
    gMTux ->Write(); 
    gMTuy ->Write(); 
    gMTkpa->Write();
    gMTmpv->Write();
    gMTsgm->Write();

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
