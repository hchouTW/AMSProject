#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
//#define __HAS_AMS_OFFICE_LIBS__

#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();
   
    //MatGeoBoxTestProp::CreateMatGeoBox();
    //return 1;

    MGConfig::JobOpt opt(argc, argv);
       
    TChain * chain = new TChain("Ntuple");
    for (auto&& file : opt.flist()) chain->Add(file.c_str());

    std::vector<double> * cx = nullptr;
    std::vector<double> * cy = nullptr;
    std::vector<double> * cz = nullptr;
    std::vector<double> * mx = nullptr;
    std::vector<double> * my = nullptr;
    std::vector<double> * mz = nullptr;

    chain->SetBranchAddress("pos_x", &cx);
    chain->SetBranchAddress("pos_y", &cy);
    chain->SetBranchAddress("pos_z", &cz);
    chain->SetBranchAddress("mom_x", &mx);
    chain->SetBranchAddress("mom_y", &my);
    chain->SetBranchAddress("mom_z", &mz);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //PhyArg::SetOpt(true, false);
    PhyArg::SetOpt(true, true);
    Int_t layBeg = 0;
    Int_t layEnd = 2;
    
    TFile * ofle = new TFile(Form("%s/prop_smc_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");

    // Prop
    Axis AXcoo("Residual [cm * p#beta/Q^{2} * L^-1]", 400, -0.08, 0.08);
    Hist * hMcx = Hist::New("hMcx", "hMcx", HistAxis(AXcoo));
    Hist * hMcy = Hist::New("hMcy", "hMcy", HistAxis(AXcoo));
    Hist * hTcx = Hist::New("hTcx", "hTcx", HistAxis(AXcoo));
    Hist * hTcy = Hist::New("hTcy", "hTcy", HistAxis(AXcoo));
    
    Axis AXagl("Residual [p#beta/Q^{2}]", 400, -0.08, 0.08);
    Hist * hMux = Hist::New("hMux", "hMux", HistAxis(AXagl));
    Hist * hMuy = Hist::New("hMuy", "hMuy", HistAxis(AXagl));
    Hist * hTux = Hist::New("hTux", "hTux", HistAxis(AXagl));
    Hist * hTuy = Hist::New("hTuy", "hTuy", HistAxis(AXagl));
    
    Axis AXrel("Rel", 400, -0.02*0.02, 0.04*0.04);
    Hist * hMrx = Hist::New("hMrx", "hMrx", HistAxis(AXrel));
    Hist * hMry = Hist::New("hMry", "hMry", HistAxis(AXrel));
    Hist * hTrx = Hist::New("hTrx", "hTrx", HistAxis(AXrel));
    Hist * hTry = Hist::New("hTry", "hTry", HistAxis(AXrel));
   
    Axis AXels("Eloss [GeV * #beta^{2}/Q^{2}]", 400, 0.0010, 0.015);
    Hist * hMee = Hist::New("hMee", "hMee", HistAxis(AXels));
    Hist * hTee = Hist::New("hTee", "hTee", HistAxis(AXels));
    
    Hist * hMcux = Hist::New("hMcux", "hMcux", HistAxis(AXcoo, AXagl));
    Hist * hMcuy = Hist::New("hMcuy", "hMcuy", HistAxis(AXcoo, AXagl));
    Hist * hTcux = Hist::New("hTcux", "hTcux", HistAxis(AXcoo, AXagl));
    Hist * hTcuy = Hist::New("hTcuy", "hTcuy", HistAxis(AXcoo, AXagl));
    
    Axis AXagl2("Residual [p#beta/Q^{2}]", 200, -0.06, 0.06);
    Axis AXrel2("Rel", 400, -0.006, 0.006);
    Hist * hMurx = Hist::New("hMurx", "hMurx", HistAxis(AXagl2, AXrel2));
    Hist * hMury = Hist::New("hMury", "hMury", HistAxis(AXagl2, AXrel2));
    Hist * hTurx = Hist::New("hTurx", "hTurx", HistAxis(AXagl2, AXrel2));
    Hist * hTury = Hist::New("hTury", "hTury", HistAxis(AXagl2, AXrel2));
    
    Long64_t printRate = chain->GetEntries()/40;
    std::cout << Form("\n==== Totally Entries %lld ====\n", chain->GetEntries());
    for (Long64_t entry = 0; entry < chain->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, chain->GetEntries());
        chain->GetEntry(entry);
        //if (entry > 100000) break; // testcode
        if (cz->size() <= layBeg) continue;
        if (cz->size() <= layEnd) continue;
        if (mz->size() <= layBeg) continue;
        if (mz->size() <= layEnd) continue;
        
        double icoo[3] = { cx->at(layBeg), cy->at(layBeg), cz->at(layBeg) };
        double ivec[3] = { mx->at(layBeg), my->at(layBeg), mz->at(layBeg) };
        double imom    = std::sqrt(ivec[0]*ivec[0] + ivec[1]*ivec[1] + ivec[2]*ivec[2]);
        double idir[3] = { ivec[0]/imom, ivec[1]/imom, ivec[2]/imom };
        
        double jcoo[3] = { cx->at(layEnd), cy->at(layEnd), cz->at(layEnd) };
        double jvec[3] = { mx->at(layEnd), my->at(layEnd), mz->at(layEnd) };
        double jmom    = std::sqrt(jvec[0]*jvec[0] + jvec[1]*jvec[1] + jvec[2]*jvec[2]);
        double jdir[3] = { jvec[0]/jmom, jvec[1]/jmom, jvec[2]/jmom };

        // Prop
        PhySt part(PartType::Proton);
        part.set_state(icoo[0], icoo[1], icoo[2], ivec[0], ivec[1], ivec[2]);

        MatFld mfld;
        PhySt pst = part;
        PropMgnt::PropToZ(jcoo[2], pst, &mfld);
        double rcoo[3] = { pst.cx(), pst.cy(), pst.cz() };
        double rdir[3] = { pst.ux(), pst.uy(), pst.uz() };
        double rmom    = pst.mom();

        pst = part;
        PropMgnt::PropToZWithMC(jcoo[2], pst);
        double tcoo[3] = { pst.cx(), pst.cy(), pst.cz() };
        double tdir[3] = { pst.ux(), pst.uy(), pst.uz() };
        double tmom    = pst.mom();

        Double_t len = std::fabs(jcoo[2]-icoo[2]);
        Double_t nrl = mfld.num_rad_len();
        Double_t ela = mfld.elcloud_abundance();
        
        Double_t scl_eloss = (part.bta() * part.bta()) / (part.chrg() * part.chrg());
        Double_t scl_mscat = (part.mom() * part.bta()) / std::fabs(part.chrg());
        scl_mscat /= std::sqrt(nrl);
        scl_eloss /= ela;
        
        Double_t mc_resc[2] = { jcoo[0] - rcoo[0], jcoo[1] - rcoo[1] };
        Double_t mc_resu[2] = { jdir[0] - rdir[0], jdir[1] - rdir[1] };
        Double_t mc_elsm    = (imom - jmom);
        Double_t tm_resc[2] = { tcoo[0] - rcoo[0], tcoo[1] - rcoo[1] };
        Double_t tm_resu[2] = { tdir[0] - rdir[0], tdir[1] - rdir[1] };
        Double_t tm_elsm    = (imom - tmom);
        hMcx->fill(scl_mscat * mc_resc[0] / len);
        hMcy->fill(scl_mscat * mc_resc[1] / len);
        hMux->fill(scl_mscat * mc_resu[0]);
        hMuy->fill(scl_mscat * mc_resu[1]);
        hMrx->fill(scl_mscat * mc_resu[0] * scl_mscat * mc_resc[0] / len);
        hMry->fill(scl_mscat * mc_resu[1] * scl_mscat * mc_resc[1] / len);
        hMee->fill(scl_eloss * mc_elsm);
        hTcx->fill(scl_mscat * tm_resc[0] / len);
        hTcy->fill(scl_mscat * tm_resc[1] / len);
        hTux->fill(scl_mscat * tm_resu[0]);
        hTuy->fill(scl_mscat * tm_resu[1]);
        hTrx->fill(scl_mscat * tm_resu[0] * scl_mscat * tm_resc[0] / len);
        hTry->fill(scl_mscat * tm_resu[1] * scl_mscat * tm_resc[1] / len);
        hTee->fill(scl_eloss * tm_elsm);
        
        hMcux->fill(scl_mscat * mc_resc[0] / len, scl_mscat * mc_resu[0]);
        hMcuy->fill(scl_mscat * mc_resc[1] / len, scl_mscat * mc_resu[1]);
        hTcux->fill(scl_mscat * tm_resc[0] / len, scl_mscat * tm_resu[0]);
        hTcuy->fill(scl_mscat * tm_resc[1] / len, scl_mscat * tm_resu[1]);
        
        hMurx->fill(scl_mscat * mc_resu[0], scl_mscat * (mc_resc[0] / (len * (1.-mfld.loc())) - mc_resu[0]));
        hMury->fill(scl_mscat * mc_resu[1], scl_mscat * (mc_resc[1] / (len * (1.-mfld.loc())) - mc_resu[1]));
        hTurx->fill(scl_mscat * tm_resu[0], scl_mscat * (tm_resc[0] / (len * (1.-mfld.loc())) - tm_resu[0]));
        hTury->fill(scl_mscat * tm_resu[1], scl_mscat * (tm_resc[1] / (len * (1.-mfld.loc())) - tm_resu[1]));
    }

    Canvas cvsMTcx(STR_FMT("cvsMTcx%02d", opt.gi()));
    cvsMTcx.create();
    hMcx->style(Fill(), Line(kBlue), Marker(kBlue));
    hTcx->style(Fill(), Line(kRed) , Marker(kRed) );
    Hist::Collect("hMTcx", "hMTcx;Residual [cm * p#beta/Q^{2} * L^-1]", HistList({hMcx, hTcx}) )->Draw("nostack");
    cvsMTcx.write();
    
    Canvas cvsMTux(STR_FMT("cvsMTux%02d", opt.gi()));
    cvsMTux.create();
    hMux->style(Fill(), Line(kBlue), Marker(kBlue));
    hTux->style(Fill(), Line(kRed) , Marker(kRed) );
    Hist::Collect("hMTux", "hMTux;Residual [p#beta/Q^{2}]", HistList({hMux, hTux}) )->Draw("nostack");
    cvsMTux.write();
    
    Canvas cvsMTee(STR_FMT("cvsMTee%02d", opt.gi()));
    cvsMTee.create();
    hMee->style(Fill(), Line(kBlue), Marker(kBlue));
    hTee->style(Fill(), Line(kRed) , Marker(kRed) );
    Hist::Collect("hMTee", "hMTee;Eloss [GeV * #beta^{2}/Q^{2}]", HistList({hMee, hTee}) )->Draw("nostack");
    cvsMTee.write();
    
    ofle->Write();
    ofle->Close();

    return 0;
}
