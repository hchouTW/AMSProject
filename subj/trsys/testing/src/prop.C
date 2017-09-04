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
    //TH1::AddDirectory(true);

    std::string ifile = "/data3/hchou/AMSData/Lorentz/proton_0.8GeV_990k.root";
    MGConfig::JobOpt opt(argc, argv);
    //for (auto&& str : opt.flist()) chain->Add(str.c_str());

    TChain * chain = new TChain("Ntuple");
    chain->Add(ifile.c_str());

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

    if (chain->GetEntries() == 0)
        MGSys::ShowErrorAndExit("No entry. Exiting...");

    double refm = 0;
    for (Long64_t it = 0; it < chain->GetEntries(); ++it) {
        chain->GetEntry(it);
        if (mz->size() == 0) continue;
        double mom = 0.001 * std::sqrt(mx->at(0)*mx->at(0) + my->at(0)*my->at(0) + mz->at(0)*mz->at(0));
        refm = mom;
        break;
    }

    Int_t  nbinc = 400;
    double binsc[2] = { -0.3, 0.3 };
    TH1D * hMcx = new TH1D("hMcx", "hMcx", nbinc, binsc[0], binsc[1]);
    TH1D * hMcy = new TH1D("hMcy", "hMcy", nbinc, binsc[0], binsc[1]);
    TH1D * hTcx = new TH1D("hTcx", "hTcx", nbinc, binsc[0], binsc[1]);
    TH1D * hTcy = new TH1D("hTcy", "hTcy", nbinc, binsc[0], binsc[1]);

    Int_t  nbinu = 400;
    double binsu[2] = { -0.015, 0.015 };
    TH1D * hMux = new TH1D("hMux", "hMux", nbinu, binsu[0], binsu[1]);
    TH1D * hMuy = new TH1D("hMuy", "hMuy", nbinu, binsu[0], binsu[1]);
    TH1D * hTux = new TH1D("hTux", "hTux", nbinu, binsu[0], binsu[1]);
    TH1D * hTuy = new TH1D("hTuy", "hTuy", nbinu, binsu[0], binsu[1]);

    Int_t  nbinm = 400;
    double binsm[2] = { 0.0, 0.001 };
    TH1D * hMee = new TH1D("hMee", "hMee", nbinm, binsm[0], binsm[1]);
    TH1D * hTee = new TH1D("hTee", "hTee", nbinm, binsm[0], binsm[1]);

    const double m2c = 0.1;
    const double m2g = 1.0e-3;
    for (Long64_t it = 0; it < chain->GetEntries(); ++it) {
        int    idi = 0;
        int    idj = 1;
        
        chain->GetEntry(it++);
        if (cz->size() != 2) continue;
        double icoo[3] = { m2c * cx->at(idi), m2c * cy->at(idi), m2c * cz->at(idi) };
        double jcoo[3] = { m2c * cx->at(idj), m2c * cy->at(idj), m2c * cz->at(idj) };
        
        chain->GetEntry(it);
        if (mz->size() != 2) continue;
        double ivec[3] = { m2g * mx->at(idi), m2g * my->at(idi), m2g * mz->at(idi) };
        double jvec[3] = { m2g * mx->at(idj), m2g * my->at(idj), m2g * mz->at(idj) };
        
        
        double imom    = std::sqrt(ivec[0]*ivec[0] + ivec[1]*ivec[1] + ivec[2]*ivec[2]);
        double idir[3] = { ivec[0]/imom, ivec[1]/imom, ivec[2]/imom };
        
        double jmom    = std::sqrt(jvec[0]*jvec[0] + jvec[1]*jvec[1] + jvec[2]*jvec[2]);
        double jdir[3] = { jvec[0]/jmom, jvec[1]/jmom, jvec[2]/jmom };
       
        // Prop
        PhySt part(PartType::Proton);
        part.set_state(icoo[0], icoo[1], icoo[2], ivec[0], ivec[1], ivec[2]);

        PhySt pst = part;
        PropMgnt::PropToZ(jcoo[2], pst);
        double rcoo[3] = { pst.cx(), pst.cy(), pst.cz() };
        double rdir[3] = { pst.ux(), pst.uy(), pst.uz() };
        double rmom    = pst.mom();

        pst = part;
        PropMgnt::PropToZWithMC(jcoo[2], pst, MatArg(true, true));
        double tcoo[3] = { pst.cx(), pst.cy(), pst.cz() };
        double tdir[3] = { pst.ux(), pst.uy(), pst.uz() };
        double tmom    = pst.mom();

        hMcx->Fill( jcoo[0]-rcoo[0] );
        hMcy->Fill( jcoo[1]-rcoo[1] );
        hTcx->Fill( tcoo[0]-rcoo[0] );
        hTcy->Fill( tcoo[1]-rcoo[1] );
        
        hMux->Fill( jdir[0]-rdir[0] );
        hMuy->Fill( jdir[1]-rdir[1] );
        hTux->Fill( tdir[0]-rdir[0] );
        hTuy->Fill( tdir[1]-rdir[1] );

        hMee->Fill( imom-jmom );
        hTee->Fill( imom-tmom );
    }

    TFile * ofile = new TFile("out.root", "RECREATE");
    hMcx->Write();
    hMcy->Write();
    hTcx->Write();
    hTcy->Write();
    hMux->Write();
    hMuy->Write();
    hTux->Write();
    hTuy->Write();
    hMee->Write();
    hTee->Write();
    ofile->Write();
    ofile->Close();

    return 0;
}
