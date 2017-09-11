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

    const int NL = 2;
    const int IDi = 0;
    const int IDj = 1;

    const double sclc = 1.0;
    const double sclu = 1.0;
    const double scle = 1.0;

    PartType part_type = PartType::Proton;
    PhySt part_info(part_type);
    TFile * ofle = new TFile("proton.root", "RECREATE");
        
    TH1D* hMu = new TH1D("hMu", "hMu;Residual [p#beta/Q^{2}];Events/Bin", 1600, -0.06 * sclu, 0.06 * sclu);
    TH1D* hTu = new TH1D("hTu", "hTu;Residual [p#beta/Q^{2}];Events/Bin", 1600, -0.06 * sclu, 0.06 * sclu);

    MGConfig::JobOpt opt(argc, argv);
    std::vector<double> momlst;
    std::vector<double> sclmscatlst;
    std::vector<double> sclelosslst;
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
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

        COUT("Current InFile : %s  Entries %ld\n", opt.flist(ifle).c_str(), chain->GetEntries());
        if (chain->GetEntries() == 0) continue;
        
        double refm = 0;
        double refe = 0;
        double refci[3];
        double refcj[3];
        for (Long64_t it = 0; it < chain->GetEntries(); ++it) {
            chain->GetEntry(it);
            if (cz->size() < NL) continue;
            if (mz->size() < NL) continue;
            refm = std::sqrt(mx->at(IDi)*mx->at(IDi) + my->at(IDi)*my->at(IDi) + mz->at(IDi)*mz->at(IDi));
            refe = std::sqrt(mx->at(IDj)*mx->at(IDj) + my->at(IDj)*my->at(IDj) + mz->at(IDj)*mz->at(IDj));
            refci[0] = cx->at(IDi);
            refci[1] = cy->at(IDi);
            refci[2] = cz->at(IDi);
            refcj[0] = cx->at(IDj);
            refcj[1] = cy->at(IDj);
            refcj[2] = cz->at(IDj);
        
            SVecD<3> vcoo(refci[0], refci[1], refci[2]);
            SVecD<3> wcoo(refcj[0], refcj[1], refcj[2]);
            MatFld&& mfld = MatMgnt::Get(vcoo, wcoo);
            mfld.print();
            
            MatFld&& mfld2 = MatMgnt::Get(SVecD<3>(0, 0, 58.3), SVecD<3>(0, 0, 57.7));
            mfld2.print();
            break;
        }
        momlst.push_back(refm);
        if (MGNumc::EqualToZero(refm)) continue;

        part_info.set_mom(refm);
        double scl_mscat = (part_info.mom() * part_info.bta() / (part_info.part().chrg() * part_info.part().chrg()));
        double scl_eloss = (part_info.bta() * part_info.bta() / (part_info.part().chrg() * part_info.part().chrg()));
        sclmscatlst.push_back(scl_mscat);
        sclelosslst.push_back(scl_eloss);
        
        std::cerr << Form("============= MOM %14.4f (%14.4f)  GB %14.4f  LOSS %8.4f ===========\n", refm, refe, part_info.gmbta(), 100.*(refm-refe)/refm );

        ofle->cd();

        Int_t  nbinc = 1600;
        double binsc[2] = { -0.25*sclc, 0.25*sclc };
        TH1D* hMcx = new TH1D(Form("hMcx%02d", ifle), Form("hMcx%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D* hMcy = new TH1D(Form("hMcy%02d", ifle), Form("hMcy%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D* hTcx = new TH1D(Form("hTcx%02d", ifle), Form("hTcx%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D* hTcy = new TH1D(Form("hTcy%02d", ifle), Form("hTcy%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);

        Int_t  nbinu = 1600;
        double binsu[2] = { -0.06*sclu, 0.06*sclu };
        TH1D* hMux = new TH1D(Form("hMux%02d", ifle), Form("hMux%02d;Residual [p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D* hMuy = new TH1D(Form("hMuy%02d", ifle), Form("hMuy%02d;Residual [p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D* hTux = new TH1D(Form("hTux%02d", ifle), Form("hTux%02d;Residual [p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D* hTuy = new TH1D(Form("hTuy%02d", ifle), Form("hTuy%02d;Residual [p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        
        Int_t  nbinm = 1600;
        double binsm[2] = { 0.0005*scle, 0.006*scle };
        TH1D* hMee = new TH1D(Form("hMee%02d", ifle), Form("hMee%02d;Energy Loss [GeV * #beta^{2}/Q^{2}];Events/Bin", ifle), nbinm, binsm[0], binsm[1]);
        TH1D* hTee = new TH1D(Form("hTee%02d", ifle), Form("hTee%02d;Energy Loss [GeV * #beta^{2}/Q^{2}];Events/Bin", ifle), nbinm, binsm[0], binsm[1]);

        for (Long64_t it = 0; it < chain->GetEntries(); ++it) {
            if (it > 100000) break; // testcode
            chain->GetEntry(it);
            if (cz->size() < NL) continue;
            if (mz->size() < NL) continue;
            
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

            hMee->Fill( scl_eloss * (imom-jmom) );
            hTee->Fill( scl_eloss * (imom-tmom) );

            if (imom > 8.01) hMu->Fill( scl_mscat * (jdir[0]-rdir[0]) );
            if (imom > 8.01) hTu->Fill( scl_mscat * (tdir[0]-rdir[0]) );
        }
    }
 
    // Coord X
    TGraphErrors* gMcx_mn = new TGraphErrors();
    TGraphErrors* gMcx_sg = new TGraphErrors();
    gMcx_mn->SetNameTitle("gMcx_mn", "");
    gMcx_sg->SetNameTitle("gMcx_sg", "");
    gMcx_mn->GetXaxis()->SetTitle("Momentum [GeV]");
    gMcx_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gMcx_mn->GetYaxis()->SetTitle("Mean [cm * p#beta/Q^{2}]");
    gMcx_sg->GetYaxis()->SetTitle("Sigma [cm * p#beta/Q^{2}]");
    
    TGraphErrors* gTcx_mn = new TGraphErrors();
    TGraphErrors* gTcx_sg = new TGraphErrors();
    gTcx_mn->SetNameTitle("gTcx_mn", "");
    gTcx_sg->SetNameTitle("gTcx_sg", "");
    gTcx_mn->GetXaxis()->SetTitle("Momentum [GeV]");
    gTcx_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gTcx_mn->GetYaxis()->SetTitle("Mean [cm * p#beta/Q^{2}]");
    gTcx_sg->GetYaxis()->SetTitle("Sigma [cm * p#beta/Q^{2}]");
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double gmbta = part_info.gmbta();
        TF1 * func = nullptr;
        double scl_mscat = 1.0 / sclmscatlst.at(ifle);
        scl_mscat = 1.0;
        
        TH1D* hMcx = (TH1D*) ofle->Get(Form("hMcx%02d", ifle));
        hMcx->Fit("gaus", "q0", "");
        func = hMcx->GetFunction("gaus");
        if (func) {
            gMcx_mn->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(1));
            gMcx_mn->SetPointError(ifle,  0., scl_mscat * func->GetParError(1));
            gMcx_sg->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(2));
            gMcx_sg->SetPointError(ifle,  0., scl_mscat * func->GetParError(2));
        }
        
        TH1D* hTcx = (TH1D*) ofle->Get(Form("hTcx%02d", ifle));
        hTcx->Fit("gaus", "q0", "");
        func = hTcx->GetFunction("gaus");
        if (func) {
            gTcx_mn->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(1));
            gTcx_mn->SetPointError(ifle,  0., scl_mscat * func->GetParError(1));
            gTcx_sg->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(2));
            gTcx_sg->SetPointError(ifle,  0., scl_mscat * func->GetParError(2));
        }
    }
    
    gMcx_mn->Write();
    gMcx_sg->Write();
    
    gTcx_mn->Write();
    gTcx_sg->Write();
   
    
    // Momentum X
    TGraphErrors* gMux_mn = new TGraphErrors();
    TGraphErrors* gMux_sg = new TGraphErrors();
    gMux_mn->SetNameTitle("gMux_mn", "");
    gMux_sg->SetNameTitle("gMux_sg", "");
    gMux_mn->GetXaxis()->SetTitle("Momentum [GeV]");
    gMux_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gMux_mn->GetYaxis()->SetTitle("Mean [cm * p#beta/Q^{2}]");
    gMux_sg->GetYaxis()->SetTitle("Sigma [cm * p#beta/Q^{2}]");
    
    TGraphErrors* gTux_mn = new TGraphErrors();
    TGraphErrors* gTux_sg = new TGraphErrors();
    gTux_mn->SetNameTitle("gTux_mn", "");
    gTux_sg->SetNameTitle("gTux_sg", "");
    gTux_mn->GetXaxis()->SetTitle("Momentum [GeV]");
    gTux_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gTux_mn->GetYaxis()->SetTitle("Mean [cm * p#beta/Q^{2}]");
    gTux_sg->GetYaxis()->SetTitle("Sigma [cm * p#beta/Q^{2}]");
    
    TGraphErrors* gMTux_sg = new TGraphErrors();
    gMTux_sg->SetNameTitle("gMTux_sg", "");
    gMTux_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gMTux_sg->GetYaxis()->SetTitle("Sigma/Sigma [1]");
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double gmbta = part_info.gmbta();
        TF1 * func = nullptr;
        double scl_mscat = 1.0 / sclmscatlst.at(ifle);
        scl_mscat = 1.0;
        
        TH1D* hMux = (TH1D*) ofle->Get(Form("hMux%02d", ifle));
        hMux->Fit("gaus", "q0", "", -0.005, 0.005);
        func = hMux->GetFunction("gaus");
        if (func) { 
            gMux_mn->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(1));
            gMux_mn->SetPointError(ifle,  0., scl_mscat * func->GetParError(1));
            gMux_sg->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(2));
            gMux_sg->SetPointError(ifle,  0., scl_mscat * func->GetParError(2));
        }
        Double_t sgM = (func) ? func->GetParameter(2) : 1;
        
        TH1D* hTux = (TH1D*) ofle->Get(Form("hTux%02d", ifle));
        hTux->Fit("gaus", "q0", "", -0.005, 0.005);
        func = hTux->GetFunction("gaus");
        if (func) {
            gTux_mn->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(1));
            gTux_mn->SetPointError(ifle,  0., scl_mscat * func->GetParError(1));
            gTux_sg->SetPoint     (ifle, mom, scl_mscat * func->GetParameter(2));
            gTux_sg->SetPointError(ifle,  0., scl_mscat * func->GetParError(2));
        }
        Double_t sgT = (func) ? func->GetParameter(2) : 1;
        
        gMTux_sg->SetPoint(ifle, gmbta, sgM/sgT);
    }
    
    gMux_mn->Write();
    gMux_sg->Write();
    
    gTux_mn->Write();
    gTux_sg->Write();
   
    gMTux_sg->Write();


    // Energy Loss
    TGraphErrors* gMee_pk = new TGraphErrors();
    TGraphErrors* gMee_sg = new TGraphErrors();
    gMee_pk->SetNameTitle("gMee_pk", "");
    gMee_sg->SetNameTitle("gMee_sg", "");
    gMee_pk->GetXaxis()->SetTitle("Momentum [GeV]");
    gMee_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gMee_pk->GetYaxis()->SetTitle("Peak [GeV * #beta^{2}/Q^{2}]");
    gMee_sg->GetYaxis()->SetTitle("Sigma [GeV * #beta^{2}/Q^{2}]");

    TGraphErrors* gTee_pk = new TGraphErrors();
    TGraphErrors* gTee_sg = new TGraphErrors();

    gTee_pk->SetNameTitle("gTee_pk", "");
    gTee_sg->SetNameTitle("gTee_sg", "");
    gTee_pk->GetXaxis()->SetTitle("Momentum [GeV]");
    gTee_sg->GetXaxis()->SetTitle("Momentum [GeV]");
    gTee_pk->GetYaxis()->SetTitle("Peak [GeV * #beta^{2}/Q^{2}]");
    gTee_sg->GetYaxis()->SetTitle("Sigma [GeV * #beta^{2}/Q^{2}]");
    
    TGraphErrors * gpkMT = new TGraphErrors();
    gpkMT->SetNameTitle("gpkMT", "");
    TGraphErrors * gpkTM = new TGraphErrors();
    gpkTM->SetNameTitle("gpkTM", "");
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double gmbta = part_info.gmbta();
        TF1 * func = nullptr;
        double scl_eloss = 1.0 / sclelosslst.at(ifle);
        scl_eloss = 1.0;

        TH1D* hMee = (TH1D*) ofle->Get(Form("hMee%02d", ifle));
        hMee->Fit("landau", "q0", "");
        func = hMee->GetFunction("landau");
        if (func) {
            gMee_pk->SetPoint     (ifle, mom, scl_eloss * func->GetParameter(1));
            gMee_pk->SetPointError(ifle,  0., scl_eloss * func->GetParError(1));
            gMee_sg->SetPoint     (ifle, mom, scl_eloss * func->GetParameter(2));
            gMee_sg->SetPointError(ifle,  0., scl_eloss * func->GetParError(2));
        }
        
        TH1D* hTee = (TH1D*) ofle->Get(Form("hTee%02d", ifle));
        hTee->Fit("landau", "q0", "");
        func = hTee->GetFunction("landau");
        if (func) {
            gTee_pk->SetPoint     (ifle, mom, scl_eloss * func->GetParameter(1));
            gTee_pk->SetPointError(ifle,  0., scl_eloss * func->GetParError(1));
            gTee_sg->SetPoint     (ifle, mom, scl_eloss * func->GetParameter(2));
            gTee_sg->SetPointError(ifle,  0., scl_eloss * func->GetParError(2));
        }

        double pkm = hMee->GetXaxis()->GetBinCenter(hMee->GetMaximumBin());
        double pkt = hTee->GetXaxis()->GetBinCenter(hTee->GetMaximumBin());
        gpkMT->SetPoint(ifle, gmbta, pkm/pkt);
        gpkTM->SetPoint(ifle, gmbta, pkt/pkm);
    }    
    
    gMee_pk->Write();
    gMee_sg->Write();

    gTee_pk->Write();
    gTee_sg->Write();

    gpkMT->Write();
    gpkTM->Write();

    ofle->Write();
    ofle->Close();

    return 0;
}
