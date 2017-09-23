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
    PhyArg::SetOpt(true, true);

    const int NL = 3;
    const int IDi = 0;
    const int IDj = 1;

    Int_t  nbinc = 1600;
    Int_t  nbinu = 1600;
    Int_t  nbinm = 1600;
    const double sclc = 1.0;
    const double sclu = 1.0;
    const double scle = 1.0;

    PartType part_type = PartType::Proton;
    PhySt part_info(part_type);
    TFile * ofle = new TFile("proton.root", "RECREATE");
    
    TH1D* hMc = new TH1D("hMc", "hMc;Residual [cm * p#beta/Q^{2}];Events/Bin", nbinc, -0.2 * sclc, 0.2 * sclc);
    TH1D* hTc = new TH1D("hTc", "hTc;Residual [cm * p#beta/Q^{2}];Events/Bin", nbinc, -0.2 * sclc, 0.2 * sclc);
        
    TH1D* hMu = new TH1D("hMu", "hMu;Residual [p#beta/Q^{2}];Events/Bin", nbinu, -0.04 * sclu, 0.04 * sclu);
    TH1D* hTu = new TH1D("hTu", "hTu;Residual [p#beta/Q^{2}];Events/Bin", nbinu, -0.04 * sclu, 0.04 * sclu);
    
    TH2D* hMcu = new TH2D("hMcu", "hMcu;Residual [cm * p#beta/Q^{2}];Residual [p#beta/Q^{2}];Events/Bin", nbinc, -0.2 * sclc, 0.2 * sclc, nbinu, -0.04 * sclu, 0.04 * sclu);
    TH2D* hTcu = new TH2D("hTcu", "hTcu;Residual [cm * p#beta/Q^{2}];Residual [p#beta/Q^{2}];Events/Bin", nbinc, -0.2 * sclc, 0.2 * sclc, nbinu, -0.04 * sclu, 0.04 * sclu);
    
    TH1D* hMe = new TH1D("hMe", "hMe;Energy Loss [GeV * #beta^{2}/Q^{2}];Events/Bin", nbinm, 0.0005 * scle, 0.007 * scle);
    TH1D* hTe = new TH1D("hTe", "hTe;Energy Loss [GeV * #beta^{2}/Q^{2}];Events/Bin", nbinm, 0.0005 * scle, 0.007 * scle);

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

        
        double binsc[2] = { -0.2*sclc, 0.2*sclc };
        TH1D* hMcx = new TH1D(Form("hMcx%02d", ifle), Form("hMcx%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D* hMcy = new TH1D(Form("hMcy%02d", ifle), Form("hMcy%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D* hTcx = new TH1D(Form("hTcx%02d", ifle), Form("hTcx%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);
        TH1D* hTcy = new TH1D(Form("hTcy%02d", ifle), Form("hTcy%02d;Residual [cm * p#beta/Q^{2}];Events/Bin", ifle), nbinc, binsc[0], binsc[1]);

        double binsu[2] = { -0.04*sclu, 0.04*sclu };
        TH1D* hMux = new TH1D(Form("hMux%02d", ifle), Form("hMux%02d;Residual [p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D* hMuy = new TH1D(Form("hMuy%02d", ifle), Form("hMuy%02d;Residual [p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D* hTux = new TH1D(Form("hTux%02d", ifle), Form("hTux%02d;Residual [p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        TH1D* hTuy = new TH1D(Form("hTuy%02d", ifle), Form("hTuy%02d;Residual [p#beta/Q^{2}];Events/Bin", ifle), nbinu, binsu[0], binsu[1]);
        
        double binsm[2] = { 0.0005*scle, 0.007*scle };
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
            //PropMgnt::Prop(std::fabs(jcoo[2]-icoo[2]), pst);
            double rcoo[3] = { pst.cx(), pst.cy(), pst.cz() };
            double rdir[3] = { pst.ux(), pst.uy(), pst.uz() };
            double rmom    = pst.mom();

            pst = part;
            PropMgnt::PropToZWithMC(jcoo[2], pst);
            //PropMgnt::PropWithMC(std::fabs(jcoo[2]-icoo[2]), pst);
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
            
            if (imom > 8.01) hMc->Fill( scl_mscat * (jcoo[0]-rcoo[0]) );
            if (imom > 8.01) hTc->Fill( scl_mscat * (tcoo[0]-rcoo[0]) );

            if (imom > 8.01) hMu->Fill( scl_mscat * (jdir[0]-rdir[0]) );
            if (imom > 8.01) hTu->Fill( scl_mscat * (tdir[0]-rdir[0]) );
            
            if (imom > 8.01) hMcu->Fill( scl_mscat * (jcoo[0]-rcoo[0]), scl_mscat * (jdir[0]-rdir[0]) );
            if (imom > 8.01) hTcu->Fill( scl_mscat * (tcoo[0]-rcoo[0]), scl_mscat * (tdir[0]-rdir[0]) );
            
            if (imom > 8.01) hMe->Fill( scl_eloss * (imom-jmom) );
            if (imom > 8.01) hTe->Fill( scl_eloss * (imom-tmom) );

            //std::cout << Form("imom %14.8f ELOSS %14.8f %14.8f\n", imom, scl_eloss * (imom-jmom), scl_eloss * (imom-tmom));
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
    TF1 * sgaus = new TF1("sgaus", "gaus", -0.1, 0.1);
    TF1 * mgaus = new TF1("mgaus", "gaus", -0.1, 0.1);
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double gmbta = part_info.gmbta();
        double bta = part_info.bta();
        double val = bta;
        
        TH1D* hMcx = (TH1D*) ofle->Get(Form("hMcx%02d", ifle));
        sgaus->SetParameters(1000, 0, hMcx->GetRMS());
        hMcx->Fit(sgaus, "q0", "");
        hMcx->Fit(sgaus, "q0", "");
        gMcx_men->SetPoint     (ifle, val, sgaus->GetParameter(1));
        gMcx_men->SetPointError(ifle,  0., sgaus->GetParError(1));
        gMcx_sgm->SetPoint     (ifle, val, sgaus->GetParameter(2));
        gMcx_sgm->SetPointError(ifle,  0., sgaus->GetParError(2));
        
        TH1D* hMcy = (TH1D*) ofle->Get(Form("hMcy%02d", ifle));
        sgaus->SetParameters(1000, 0, hMcy->GetRMS());
        hMcy->Fit(sgaus, "q0", "");
        hMcy->Fit(sgaus, "q0", "");
        gMcy_men->SetPoint     (ifle, val, sgaus->GetParameter(1));
        gMcy_men->SetPointError(ifle,  0., sgaus->GetParError(1));
        gMcy_sgm->SetPoint     (ifle, val, sgaus->GetParameter(2));
        gMcy_sgm->SetPointError(ifle,  0., sgaus->GetParError(2));
        
        TH1D* hMux = (TH1D*) ofle->Get(Form("hMux%02d", ifle));
        sgaus->SetParameters(1000, 0, hMux->GetRMS());
        hMux->Fit(sgaus, "q0", "");
        hMux->Fit(sgaus, "q0", "");
        gMux_men->SetPoint     (ifle, val, sgaus->GetParameter(1));
        gMux_men->SetPointError(ifle,  0., sgaus->GetParError(1));
        gMux_sgm->SetPoint     (ifle, val, sgaus->GetParameter(2));
        gMux_sgm->SetPointError(ifle,  0., sgaus->GetParError(2));
        
        TH1D* hMuy = (TH1D*) ofle->Get(Form("hMuy%02d", ifle));
        sgaus->SetParameters(1000, 0, hMuy->GetRMS());
        hMuy->Fit(sgaus, "q0", "");
        hMuy->Fit(sgaus, "q0", "");
        gMuy_men->SetPoint     (ifle, val, sgaus->GetParameter(1));
        gMuy_men->SetPointError(ifle,  0., sgaus->GetParError(1));
        gMuy_sgm->SetPoint     (ifle, val, sgaus->GetParameter(2));
        gMuy_sgm->SetPointError(ifle,  0., sgaus->GetParError(2));
    }

    gMcx_men->Write();
    gMcx_sgm->Write();
    gMcy_men->Write();
    gMcy_sgm->Write();
    gMux_men->Write();
    gMux_sgm->Write();
    gMuy_men->Write();
    gMuy_sgm->Write();


    TGraphErrors* gTcx_men = new TGraphErrors(); gTcx_men->SetNameTitle("gTcx_men", "");
    TGraphErrors* gTcx_sgm = new TGraphErrors(); gTcx_sgm->SetNameTitle("gTcx_sgm", "");
    TGraphErrors* gTcy_men = new TGraphErrors(); gTcy_men->SetNameTitle("gTcy_men", "");
    TGraphErrors* gTcy_sgm = new TGraphErrors(); gTcy_sgm->SetNameTitle("gTcy_sgm", "");
    TGraphErrors* gTux_men = new TGraphErrors(); gTux_men->SetNameTitle("gTux_men", "");
    TGraphErrors* gTux_sgm = new TGraphErrors(); gTux_sgm->SetNameTitle("gTux_sgm", "");
    TGraphErrors* gTuy_men = new TGraphErrors(); gTuy_men->SetNameTitle("gTuy_men", "");
    TGraphErrors* gTuy_sgm = new TGraphErrors(); gTuy_sgm->SetNameTitle("gTuy_sgm", "");
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double gmbta = part_info.gmbta();
        double bta = part_info.bta();
        double val = bta;
        
        TH1D* hTcx = (TH1D*) ofle->Get(Form("hTcx%02d", ifle));
        sgaus->SetParameters(1000, 0, hTcx->GetRMS());
        hTcx->Fit(sgaus, "q0", "");
        hTcx->Fit(sgaus, "q0", "");
        gTcx_men->SetPoint     (ifle, val, sgaus->GetParameter(1));
        gTcx_men->SetPointError(ifle,  0., sgaus->GetParError(1));
        gTcx_sgm->SetPoint     (ifle, val, sgaus->GetParameter(2));
        gTcx_sgm->SetPointError(ifle,  0., sgaus->GetParError(2));
        
        TH1D* hTcy = (TH1D*) ofle->Get(Form("hTcy%02d", ifle));
        sgaus->SetParameters(1000, 0, hTcy->GetRMS());
        hTcy->Fit(sgaus, "q0", "");
        hTcy->Fit(sgaus, "q0", "");
        gTcy_men->SetPoint     (ifle, val, sgaus->GetParameter(1));
        gTcy_men->SetPointError(ifle,  0., sgaus->GetParError(1));
        gTcy_sgm->SetPoint     (ifle, val, sgaus->GetParameter(2));
        gTcy_sgm->SetPointError(ifle,  0., sgaus->GetParError(2));
        
        TH1D* hTux = (TH1D*) ofle->Get(Form("hTux%02d", ifle));
        sgaus->SetParameters(1000, 0, hTux->GetRMS());
        hTux->Fit(sgaus, "q0", "");
        hTux->Fit(sgaus, "q0", "");
        gTux_men->SetPoint     (ifle, val, sgaus->GetParameter(1));
        gTux_men->SetPointError(ifle,  0., sgaus->GetParError(1));
        gTux_sgm->SetPoint     (ifle, val, sgaus->GetParameter(2));
        gTux_sgm->SetPointError(ifle,  0., sgaus->GetParError(2));
        
        TH1D* hTuy = (TH1D*) ofle->Get(Form("hTuy%02d", ifle));
        sgaus->SetParameters(1000, 0, hTuy->GetRMS());
        hTuy->Fit(sgaus, "q0", "");
        hTuy->Fit(sgaus, "q0", "");
        gTuy_men->SetPoint     (ifle, val, sgaus->GetParameter(1));
        gTuy_men->SetPointError(ifle,  0., sgaus->GetParError(1));
        gTuy_sgm->SetPoint     (ifle, val, sgaus->GetParameter(2));
        gTuy_sgm->SetPointError(ifle,  0., sgaus->GetParError(2));
    }

    gTcx_men->Write();
    gTcx_sgm->Write();
    gTcy_men->Write();
    gTcy_sgm->Write();
    gTux_men->Write();
    gTux_sgm->Write();
    gTuy_men->Write();
    gTuy_sgm->Write();

    TGraphErrors* gMTcx_sgm  = new TGraphErrors(); gMTcx_sgm ->SetNameTitle("gMTcx_sgm",  "");
    TGraphErrors* gMTcy_sgm  = new TGraphErrors(); gMTcy_sgm ->SetNameTitle("gMTcy_sgm",  "");
    TGraphErrors* gMTux_sgm  = new TGraphErrors(); gMTux_sgm ->SetNameTitle("gMTux_sgm",  "");
    TGraphErrors* gMTuy_sgm  = new TGraphErrors(); gMTuy_sgm ->SetNameTitle("gMTuy_sgm",  "");
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        gMTcx_sgm->SetPoint     (ifle, gMcx_sgm->GetX()[ifle], gMcx_sgm->GetY()[ifle]/gTcx_sgm->GetY()[ifle]);
        gMTcx_sgm->SetPointError(ifle, 0., 0.);
        gMTcy_sgm->SetPoint     (ifle, gMcy_sgm->GetX()[ifle], gMcy_sgm->GetY()[ifle]/gTcy_sgm->GetY()[ifle]);
        gMTcy_sgm->SetPointError(ifle, 0., 0.);
        gMTux_sgm->SetPoint     (ifle, gMux_sgm->GetX()[ifle], gMux_sgm->GetY()[ifle]/gTux_sgm->GetY()[ifle]);
        gMTux_sgm->SetPointError(ifle, 0., 0.);
        gMTuy_sgm->SetPoint     (ifle, gMuy_sgm->GetX()[ifle], gMuy_sgm->GetY()[ifle]/gTuy_sgm->GetY()[ifle]);
        gMTuy_sgm->SetPointError(ifle, 0., 0.);
    }
   
    gMTcx_sgm->Write();
    gMTcy_sgm->Write();
    gMTux_sgm->Write();
    gMTuy_sgm->Write();

    // Energy Loss
    TF1 * feloss = new TF1("feloss", "[0] * TMath::Power( ([2]/x)/[1]/[1], ([2]/x)/[1]/[1] ) / TMath::Gamma( ([2]/x)/[1]/[1] ) * TMath::Exp(-(([2]/x)/[1]/[1]) * ((x-[2])/[3] + TMath::Exp(-(x-[2])/[3])) )");
    feloss->SetParameters(1000., 1.0, 0.001, 0.0002); 
    
    TGraphErrors* gMee_peak = new TGraphErrors(); gMee_peak->SetNameTitle("gMee_peak", "");
    TGraphErrors* gMee_kpa  = new TGraphErrors(); gMee_kpa ->SetNameTitle("gMee_kpa",  "");
    TGraphErrors* gMee_mpv  = new TGraphErrors(); gMee_mpv ->SetNameTitle("gMee_mpv",  "");
    TGraphErrors* gMee_sgm  = new TGraphErrors(); gMee_sgm ->SetNameTitle("gMee_sgm",  "");
    TGraphErrors* gMee_mos  = new TGraphErrors(); gMee_mos ->SetNameTitle("gMee_mos",  "");
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double gmbta = part_info.gmbta();
        double bta = part_info.bta();
        double val = bta;

        TH1D* hMee = (TH1D*) ofle->Get(Form("hMee%02d", ifle));
        double peak = hMee->GetXaxis()->GetBinCenter(hMee->GetMaximumBin());
        feloss->SetParameters(1000., 1.0, peak, 0.05*peak);
        hMee->Fit(feloss, "q0", "");
        hMee->Fit(feloss, "q0", "");
        hMee->Fit(feloss, "q0", "");

        gMee_peak->SetPoint     (ifle, val, peak);
        gMee_peak->SetPointError(ifle,    0., 0.);
        gMee_kpa->SetPoint     (ifle, val, feloss->GetParameter(1));
        gMee_kpa->SetPointError(ifle,    0., feloss->GetParError(1));
        gMee_mpv->SetPoint     (ifle, val, feloss->GetParameter(2));
        gMee_mpv->SetPointError(ifle,    0., feloss->GetParError(2));
        gMee_sgm->SetPoint     (ifle, val, feloss->GetParameter(3));
        gMee_sgm->SetPointError(ifle,    0., feloss->GetParError(3));
        gMee_mos->SetPoint     (ifle, val, feloss->GetParameter(2)/feloss->GetParameter(3));
        gMee_mos->SetPointError(ifle,    0., 0.);
    }
    
    TGraphErrors* gTee_peak = new TGraphErrors(); gTee_peak->SetNameTitle("gTee_peak", "");
    TGraphErrors* gTee_kpa  = new TGraphErrors(); gTee_kpa ->SetNameTitle("gTee_kpa",  "");
    TGraphErrors* gTee_mpv  = new TGraphErrors(); gTee_mpv ->SetNameTitle("gTee_mpv",  "");
    TGraphErrors* gTee_sgm  = new TGraphErrors(); gTee_sgm ->SetNameTitle("gTee_sgm",  "");
    TGraphErrors* gTee_mos  = new TGraphErrors(); gTee_mos ->SetNameTitle("gTee_mos",  "");
    
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        double mom = momlst.at(ifle);
        part_info.set_mom(mom);
        double gmbta = part_info.gmbta();
        double bta = part_info.bta();
        double val = bta;

        TH1D* hTee = (TH1D*) ofle->Get(Form("hTee%02d", ifle));
        double peak = hTee->GetXaxis()->GetBinCenter(hTee->GetMaximumBin());
        feloss->SetParameters(1000., 1.0, peak, 0.05*peak);
        hTee->Fit(feloss, "q0", "");
        hTee->Fit(feloss, "q0", "");
        hTee->Fit(feloss, "q0", "");

        gTee_peak->SetPoint     (ifle, val, peak);
        gTee_peak->SetPointError(ifle,    0., 0.);
        gTee_kpa->SetPoint     (ifle, val, feloss->GetParameter(1));
        gTee_kpa->SetPointError(ifle,    0., feloss->GetParError(1));
        gTee_mpv->SetPoint     (ifle, val, feloss->GetParameter(2));
        gTee_mpv->SetPointError(ifle,    0., feloss->GetParError(2));
        gTee_sgm->SetPoint     (ifle, val, feloss->GetParameter(3));
        gTee_sgm->SetPointError(ifle,    0., feloss->GetParError(3));
        gTee_mos->SetPoint     (ifle, val, feloss->GetParameter(2)/feloss->GetParameter(3));
        gTee_mos->SetPointError(ifle,    0., 0.);
    }
    
    TGraphErrors* gMTee_peak = new TGraphErrors(); gMTee_peak->SetNameTitle("gMTee_peak", "");
    TGraphErrors* gMTee_kpa  = new TGraphErrors(); gMTee_kpa ->SetNameTitle("gMTee_kpa",  "");
    TGraphErrors* gMTee_mpv  = new TGraphErrors(); gMTee_mpv ->SetNameTitle("gMTee_mpv",  "");
    TGraphErrors* gMTee_sgm  = new TGraphErrors(); gMTee_sgm ->SetNameTitle("gMTee_sgm",  "");
    TGraphErrors* gMTee_mos  = new TGraphErrors(); gMTee_mos ->SetNameTitle("gMTee_mos",  "");
    for (int ifle = 0; ifle < opt.fsize(); ++ifle) {
        gMTee_peak->SetPoint     (ifle, gMee_peak->GetX()[ifle], gMee_peak->GetY()[ifle]/gTee_peak->GetY()[ifle]);
        gMTee_peak->SetPointError(ifle, 0., 0.);
        gMTee_kpa->SetPoint     (ifle, gMee_kpa->GetX()[ifle], gMee_kpa->GetY()[ifle]/gTee_kpa->GetY()[ifle]);
        gMTee_kpa->SetPointError(ifle, 0., 0.);
        gMTee_mpv->SetPoint     (ifle, gMee_mpv->GetX()[ifle], gMee_mpv->GetY()[ifle]/gTee_mpv->GetY()[ifle]);
        gMTee_mpv->SetPointError(ifle, 0., 0.);
        gMTee_sgm->SetPoint     (ifle, gMee_sgm->GetX()[ifle], gMee_sgm->GetY()[ifle]/gTee_sgm->GetY()[ifle]);
        gMTee_sgm->SetPointError(ifle, 0., 0.);
        gMTee_mos->SetPoint     (ifle, gMee_mos->GetX()[ifle], gMee_mos->GetY()[ifle]/gTee_mos->GetY()[ifle]);
        gMTee_mos->SetPointError(ifle, 0., 0.);
    }
    
    gMee_peak->Write();
    gMee_kpa ->Write();
    gMee_mpv ->Write();
    gMee_sgm ->Write();
    gMee_mos ->Write();
    
    gTee_peak->Write();
    gTee_kpa ->Write();
    gTee_mpv ->Write();
    gTee_sgm ->Write();
    gTee_mos ->Write();
    
    gMTee_peak->Write();
    gMTee_kpa ->Write();
    gMTee_mpv ->Write();
    gMTee_sgm ->Write();
    gMTee_mos ->Write();
        
    ofle->Write();
    ofle->Close();

    return 0;
}
