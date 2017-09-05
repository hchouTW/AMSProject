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
    TH1::AddDirectory(true);

    const int NL = 3;
    const int IDi = 0;
    const int IDj = 1;

    PartType part_type = PartType::Proton;
    PhySt part_info(part_type);
    TFile * ofle = new TFile("proton.root", "RECREATE");

    MGConfig::JobOpt opt(argc, argv);
    std::vector<double> momlst;
    std::vector<double> sclmscatlst;
    std::vector<double> sclenglslst;
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        COUT("Current InFile : %s\n", opt.flist(ifle).c_str());
        TChain * chain = new TChain("Ntuple");
        chain->Add(opt.flist(ifle).c_str());

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

        if (chain->GetEntries() == 0) continue;
        
        double refm = 0;
        double refci[3];
        double refcj[3];
        for (Long64_t it = 0; it < chain->GetEntries(); ++it) {
            chain->GetEntry(it);
            if (cz->size() != NL) continue;
            if (mz->size() != NL) continue;
            refm = std::sqrt(mx->at(IDi)*mx->at(IDi) + my->at(IDi)*my->at(IDi) + mz->at(IDi)*mz->at(IDi));
            refci[0] = cx->at(IDi);
            refci[1] = cy->at(IDi);
            refci[2] = cz->at(IDi);
            refcj[0] = cx->at(IDj);
            refcj[1] = cy->at(IDj);
            refcj[2] = cz->at(IDj);
            break;
        }
        momlst.push_back(refm);

        part_info.set_mom(refm);
        double scl_mscat = (part_info.mom() * part_info.bta() / (part_info.part().chrg() * part_info.part().chrg()));
        double scl_engls = (part_info.bta() * part_info.bta() / (part_info.part().chrg() * part_info.part().chrg()));
        sclmscatlst.push_back(scl_mscat);
        sclenglslst.push_back(scl_engls);

        ofle->cd();

        Int_t  nbinc = 800;
        double binsc[2] = { -0.2, 0.2 };
        TH1D * hMcx = new TH1D(Form("hMcx%02d", ifle), Form("hMcx%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D * hMcy = new TH1D(Form("hMcy%02d", ifle), Form("hMcy%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D * hTcx = new TH1D(Form("hTcx%02d", ifle), Form("hTcx%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D * hTcy = new TH1D(Form("hTcy%02d", ifle), Form("hTcy%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);

        Int_t  nbinu = 800;
        double binsu[2] = { -0.04, 0.04 };
        TH1D * hMux = new TH1D(Form("hMux%02d", ifle), Form("hMux%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D * hMuy = new TH1D(Form("hMuy%02d", ifle), Form("hMuy%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D * hTux = new TH1D(Form("hTux%02d", ifle), Form("hTux%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D * hTuy = new TH1D(Form("hTuy%02d", ifle), Form("hTuy%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);

        Int_t  nbinm = 1600;
        double binsm[2] = { 0.005, 0.06 };
        TH1D * hMee = new TH1D(Form("hMee%02d", ifle), Form("hMee%02d;Energy Loss [GeV * #beta^{2}/Q^{2}];Events/Bin", ifle), nbinm, binsm[0], binsm[1]);
        TH1D * hTee = new TH1D(Form("hTee%02d", ifle), Form("hTee%02d;Energy Loss [GeV * #beta^{2}/Q^{2}];Events/Bin", ifle), nbinm, binsm[0], binsm[1]);

        for (Long64_t it = 0; it < chain->GetEntries(); ++it) {
            chain->GetEntry(it);
            if (cz->size() != NL) continue;
            if (mz->size() != NL) continue;
            
            double icoo[3] = { cx->at(IDi), cy->at(IDi), cz->at(IDi) };
            double ivec[3] = { mx->at(IDi), my->at(IDi), mz->at(IDi) };
            double imom    = std::sqrt(ivec[0]*ivec[0] + ivec[1]*ivec[1] + ivec[2]*ivec[2]);
            double idir[3] = { ivec[0]/imom, ivec[1]/imom, ivec[2]/imom };
            
            double jcoo[3] = { cx->at(IDj), cy->at(IDj), cz->at(IDj) };
            double jvec[3] = { mx->at(IDj), my->at(IDj), mz->at(IDj) };
            double jmom    = std::sqrt(jvec[0]*jvec[0] + jvec[1]*jvec[1] + jvec[2]*jvec[2]);
            double jdir[3] = { jvec[0]/jmom, jvec[1]/jmom, jvec[2]/jmom };
           
            // Prop
            PhySt part(part_type);
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

            hMcx->Fill( scl_mscat * (jcoo[0]-rcoo[0]) );
            hMcy->Fill( scl_mscat * (jcoo[1]-rcoo[1]) );
            hTcx->Fill( scl_mscat * (tcoo[0]-rcoo[0]) );
            hTcy->Fill( scl_mscat * (tcoo[1]-rcoo[1]) );
            
            hMux->Fill( scl_mscat * (jdir[0]-rdir[0]) );
            hMuy->Fill( scl_mscat * (jdir[1]-rdir[1]) );
            hTux->Fill( scl_mscat * (tdir[0]-rdir[0]) );
            hTuy->Fill( scl_mscat * (tdir[1]-rdir[1]) );

            hMee->Fill( scl_engls * (imom-jmom) );
            hTee->Fill( scl_engls * (imom-tmom) );
        }
    }
   
    // Coord X
    TGraphErrors * gMcx_mn = new TGraphErrors();
    TGraphErrors * gMcx_sg = new TGraphErrors();
    gMcx_mn->SetNameTitle("gMcx_mn", "");
    gMcx_sg->SetNameTitle("gMcx_sg", "");
    gMcx_mn->GetXaxis()->SetTitle("Momentum [GeV]");
    gMcx_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gMcx_mn->GetYaxis()->SetTitle("Mean [cm * p#beta/Q^{2}]");
    gMcx_sg->GetYaxis()->SetTitle("Sigma [cm * p#beta/Q^{2}]");
    
    TGraphErrors * gTcx_mn = new TGraphErrors();
    TGraphErrors * gTcx_sg = new TGraphErrors();
    gTcx_mn->SetNameTitle("gTcx_mn", "");
    gTcx_sg->SetNameTitle("gTcx_sg", "");
    gTcx_mn->GetXaxis()->SetTitle("Momentum [GeV]");
    gTcx_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gTcx_mn->GetYaxis()->SetTitle("Mean [cm * p#beta/Q^{2}]");
    gTcx_sg->GetYaxis()->SetTitle("Sigma [cm * p#beta/Q^{2}]");
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double scl_mscat = 1.0 / sclmscatlst.at(ifle);
        TF1 * func = nullptr;

        scl_mscat = 1.0;
        
        TH1D * hMcx = (TH1D*) ofle->Get(Form("hMcx%02d", ifle));
        hMcx->Fit("gaus", "q0", "");
        func = hMcx->GetFunction("gaus");
        gMcx_mn->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(1));
        gMcx_mn->SetPointError(ifle,  0., scl_mscat * func->GetParError(1));
        gMcx_sg->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(2));
        gMcx_sg->SetPointError(ifle,  0., scl_mscat * func->GetParError(2));
        
        TH1D * hTcx = (TH1D*) ofle->Get(Form("hTcx%02d", ifle));
        hTcx->Fit("gaus", "q0", "");
        func = hTcx->GetFunction("gaus");
        gTcx_mn->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(1));
        gTcx_mn->SetPointError(ifle,  0., scl_mscat * func->GetParError(1));
        gTcx_sg->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(2));
        gTcx_sg->SetPointError(ifle,  0., scl_mscat * func->GetParError(2));
    }
    
    gMcx_mn->Write();
    gMcx_sg->Write();
    
    gTcx_mn->Write();
    gTcx_sg->Write();
    
    // Energy Loss
    TGraphErrors * gMee_pk = new TGraphErrors();
    TGraphErrors * gMee_sg = new TGraphErrors();
    gMee_pk->SetNameTitle("gMee_pk", "");
    gMee_sg->SetNameTitle("gMee_sg", "");
    gMee_pk->GetXaxis()->SetTitle("Momentum [GeV]");
    gMee_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gMee_pk->GetYaxis()->SetTitle("Peak [GeV * #beta^{2}/Q^{2}]");
    gMee_sg->GetYaxis()->SetTitle("Sigma [GeV * #beta^{2}/Q^{2}]");

    TGraphErrors * gTee_pk = new TGraphErrors();
    TGraphErrors * gTee_sg = new TGraphErrors();
    gTee_pk->SetNameTitle("gTee_pk", "");
    gTee_sg->SetNameTitle("gTee_sg", "");
    gTee_pk->GetXaxis()->SetTitle("Momentum [GeV]");
    gTee_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gTee_pk->GetYaxis()->SetTitle("Peak [GeV * #beta^{2}/Q^{2}]");
    gTee_sg->GetYaxis()->SetTitle("Sigma [GeV * #beta^{2}/Q^{2}]");
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double scl_engls = 1.0 / sclenglslst.at(ifle);
        TF1 * func = nullptr;
       
        scl_engls = 1.0;

        TH1D * hMee = (TH1D*) ofle->Get(Form("hMee%02d", ifle));
        hMee->Fit("landau", "q0", "");
        func = hMee->GetFunction("landau");
        gMee_pk->SetPoint     (ifle, mom, scl_engls * func->GetParameter(1));
        gMee_pk->SetPointError(ifle,  0., scl_engls * func->GetParError(1));
        gMee_sg->SetPoint     (ifle, mom, scl_engls * func->GetParameter(2));
        gMee_sg->SetPointError(ifle,  0., scl_engls * func->GetParError(2));
        
        TH1D * hTee = (TH1D*) ofle->Get(Form("hTee%02d", ifle));
        hTee->Fit("landau", "q0", "");
        func = hTee->GetFunction("landau");
        gTee_pk->SetPoint     (ifle, mom, scl_engls * func->GetParameter(1));
        gTee_pk->SetPointError(ifle,  0., scl_engls * func->GetParError(1));
        gTee_sg->SetPoint     (ifle, mom, scl_engls * func->GetParameter(2));
        gTee_sg->SetPointError(ifle,  0., scl_engls * func->GetParError(2));
    }    
    
    gMee_pk->Write();
    gMee_sg->Write();

    gTee_pk->Write();
    gTee_sg->Write();
    
    ofle->Write();
    ofle->Close();

    return 0;
}
