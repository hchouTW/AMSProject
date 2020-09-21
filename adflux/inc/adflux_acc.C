#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "HistFit1D.h"
#include "HistFit1D.C"

static constexpr double Mproton    = 0.938272297;
static constexpr double Mdeuterium = 1.876123915;

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "66";
    
    short eaxis = 1; // 1, rig   0, ken
    std::string sube   = eaxis ? "rig" : "ken";
    std::string etitle = eaxis ? "Rigidity [GV]" : "Kinetic Energy per Nucleon [GeV/n]";

    TF1 *func_pol0 = new TF1("pol0", "[0]", 0.5, 100.0);
    func_pol0->SetParameter(0, 0.0);

    ////////////////////
    //==== INPUT =====//
    ////////////////////
    std::string path = "/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15";
    TFile* fiss = TFile::Open(Form("%s/iss%s/YiMdst.root", path.c_str(), subv.c_str()));
    TFile* ftme = TFile::Open("out/lvtme/RKlvtme.root");
    
    const double acc_geom_factor = (3.9 * 3.9 * TMath::Pi()); // (~47.78 m^2 sr)

    UInt_t cntev = 0;
        
    UInt_t cntad = 0;
    TFile* fmcad = TFile::Open(Form("%s/mcad%s/YiMdst.root", path.c_str(), subv.c_str()));
    TTree* tmcad = (TTree*)fmcad->Get("runlist");
    tmcad->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcad->GetEntries(); ++it) { tmcad->GetEntry(it); cntad+=cntev; }
    double density_fact_ad = static_cast<double>(cntad) / (std::log(1000.0) - std::log(0.2));

    //////////////////
    //==== TOF =====//
    //////////////////
    TH1D* hmcadTF_cnt = (TH1D*)fmcad->Get(Form("hTF%s_cnt_MC", sube.c_str()));
    Axis AXTFe(etitle, hmcadTF_cnt->GetXaxis());
   
    // TOF - Width
    Hist* hTFde = Hist::New("hTFde", HistAxis(AXTFe, eaxis?"Width [GV]":"Width [GeV/n]"));
    for (int ib = 1; ib <= AXTFe.nbin(); ++ib) {
        (*hTFde)()->SetBinContent(ib, AXTFe.width(ib));
        (*hTFde)()->SetBinError  (ib, 0.0);
    }

    // TOF - Accp
    Hist* hTFaccp = Hist::New("hTFaccp", HistAxis(AXTFe, "Effective Acceptance [m^{2} sr]"));
    for (int ib = 1; ib <= AXTFe.nbin(); ++ib) {
        double moml = AXTFe(ib - 1);
        double momu = AXTFe(ib);
        if (!eaxis) {
            moml = Mdeuterium * std::sqrt((1.0 + moml / Mproton) * (1.0 + moml / Mproton) - 1.0);
            momu = Mdeuterium * std::sqrt((1.0 + momu / Mproton) * (1.0 + momu / Mproton) - 1.0);
        }
        double dlogM = std::log(momu) - std::log(moml);
        double gcnt  = density_fact_ad * dlogM;
        double rerr  = hmcadTF_cnt->GetBinError(ib) / hmcadTF_cnt->GetBinContent(ib);
        
        double acc_val = acc_geom_factor * (hmcadTF_cnt->GetBinContent(ib) / gcnt);
        double acc_err = acc_val * rerr;
        
        if (!std::isfinite(acc_val) || acc_val < 0.0) continue;
        if (!std::isfinite(acc_err) || acc_err < 0.0) continue;
        (*hTFaccp)()->SetBinContent(ib, acc_val);
        (*hTFaccp)()->SetBinError  (ib, acc_err);
    }

    // TOF - Meff
    TH1D* hmcadTF_effn = (TH1D*)fmcad->Get(Form("hTF%s_effn_MC", sube.c_str()));
    TH1D* hmcadTF_effd = (TH1D*)fmcad->Get(Form("hTF%s_effd_MC", sube.c_str()));
    
    Hist* hTFmeff = Hist::New("hTFmeff", HistAxis(AXTFe, "Mass Selection Efficiency"));
    for (int ib = 1; ib <= AXTFe.nbin(); ++ib) {
        double meff = hmcadTF_effn->GetBinContent(ib) / hmcadTF_effd->GetBinContent(ib);
        double rerr = hmcadTF_effn->GetBinError  (ib) / hmcadTF_effn->GetBinContent(ib);
        (*hTFmeff)()->SetBinContent(ib, meff);
        (*hTFmeff)()->SetBinError  (ib, meff * rerr);
    }
    (*hTFmeff)()->Fit(func_pol0, "q0", "");
    double TFmeff = func_pol0->GetParameter(0);

    // TOF - Trig
    TH1D* hissTF_trgn = (TH1D*)fiss->Get(Form("hTF%s_trgn", sube.c_str()));
    TH1D* hissTF_trgd = (TH1D*)fiss->Get(Form("hTF%s_trgd", sube.c_str()));

    Hist* hTFtrig = Hist::New("hTFtrig", HistAxis(AXTFe, "Trigger Efficiency"));
    for (int ib = 1; ib <= AXTFe.nbin(); ++ib) {
        double trig = hissTF_trgn->GetBinContent(ib) / hissTF_trgd->GetBinContent(ib);
        double rerr = hissTF_trgn->GetBinError  (ib) / hissTF_trgn->GetBinContent(ib);
        if (!std::isfinite(trig)) { trig = 1.0; rerr = 0.0; }
        (*hTFtrig)()->SetBinContent(ib, trig);
        (*hTFtrig)()->SetBinError  (ib, trig * rerr);
    }
    (*hTFtrig)()->Fit(func_pol0, "q0", "");
    double TFtrig = func_pol0->GetParameter(0);

    // TOF - Livetime
    Hist* hTFlv = Hist::New("hTFlv", (TH1D*)(ftme->Get(Form("hTF%s_expt", sube.c_str()))));
    (*hTFlv)()->GetXaxis()->SetTitle(etitle.c_str());
    (*hTFlv)()->GetYaxis()->SetTitle("Exposure Time [s]");

    Hist* hTFcorr = Hist::New("hTFcorr", HistAxis(AXTFe, eaxis?"[(GV) m^{2} sr s]":"[(GeV/n) m^{2} sr s]"));
    for (int ib = 1; ib <= AXTFe.nbin(); ++ib) {
        double de   = (*hTFde  )()->GetBinContent(ib);
        double lv   = (*hTFlv  )()->GetBinContent(ib);
        double accp = (*hTFaccp)()->GetBinContent(ib);
        double meff = (*hTFmeff)()->GetBinContent(ib);
        double trig = (*hTFtrig)()->GetBinContent(ib);
        //double corr = de * lv * accp * meff * trig;
        double corr = de * lv * accp * TFmeff * TFtrig;

        (*hTFcorr)()->SetBinContent(ib, corr);
        (*hTFcorr)()->SetBinError  (ib, 0);
    }




    ///////////////////
    //==== RICH =====//
    ///////////////////
    TH1D* hmcadRH_cnt = (TH1D*)fmcad->Get(Form("hRH%s_cnt_MC", sube.c_str()));
    Axis AXRHe(etitle, hmcadRH_cnt->GetXaxis());
   
    // RICH - Width
    Hist* hRHde = Hist::New("hRHde", HistAxis(AXRHe, eaxis?"Width [GV]":"Width [GeV/n]"));
    for (int ib = 1; ib <= AXRHe.nbin(); ++ib) {
        (*hRHde)()->SetBinContent(ib, AXRHe.width(ib));
        (*hRHde)()->SetBinError  (ib, 0.0);
    }

    // RICH - Accp
    Hist* hRHaccp = Hist::New("hRHaccp", HistAxis(AXRHe, "Effective Acceptance [m^{2} sr]"));
    for (int ib = 1; ib <= AXRHe.nbin(); ++ib) {
        double moml = AXRHe(ib - 1);
        double momu = AXRHe(ib);
        if (!eaxis) {
            moml = Mdeuterium * std::sqrt((1.0 + moml / Mproton) * (1.0 + moml / Mproton) - 1.0);
            momu = Mdeuterium * std::sqrt((1.0 + momu / Mproton) * (1.0 + momu / Mproton) - 1.0);
        }
        double dlogM = std::log(momu) - std::log(moml);
        double gcnt  = density_fact_ad * dlogM;
        double rerr  = hmcadRH_cnt->GetBinError(ib) / hmcadRH_cnt->GetBinContent(ib);
        
        double acc_val = acc_geom_factor * (hmcadRH_cnt->GetBinContent(ib) / gcnt);
        double acc_err = acc_val * rerr;
        
        if (!std::isfinite(acc_val) || acc_val < 0.0) continue;
        if (!std::isfinite(acc_err) || acc_err < 0.0) continue;
        (*hRHaccp)()->SetBinContent(ib, acc_val);
        (*hRHaccp)()->SetBinError  (ib, acc_err);
    }

    // RICH - Meff
    TH1D* hmcadRH_effn = (TH1D*)fmcad->Get(Form("hRH%s_effn_MC", sube.c_str()));
    TH1D* hmcadRH_effd = (TH1D*)fmcad->Get(Form("hRH%s_effd_MC", sube.c_str()));
    
    Hist* hRHmeff = Hist::New("hRHmeff", HistAxis(AXRHe, "Mass Selection Efficiency"));
    for (int ib = 1; ib <= AXRHe.nbin(); ++ib) {
        double meff = hmcadRH_effn->GetBinContent(ib) / hmcadRH_effd->GetBinContent(ib);
        double rerr = hmcadRH_effn->GetBinError  (ib) / hmcadRH_effn->GetBinContent(ib);
        (*hRHmeff)()->SetBinContent(ib, meff);
        (*hRHmeff)()->SetBinError  (ib, meff * rerr);
    }
    (*hRHmeff)()->Fit(func_pol0, "q0", "");
    double RHmeff = func_pol0->GetParameter(0);

    // RICH - Trig
    TH1D* hissRH_trgn = (TH1D*)fiss->Get(Form("hRH%s_trgn", sube.c_str()));
    TH1D* hissRH_trgd = (TH1D*)fiss->Get(Form("hRH%s_trgd", sube.c_str()));

    Hist* hRHtrig = Hist::New("hRHtrig", HistAxis(AXRHe, "Trigger Efficiency"));
    for (int ib = 1; ib <= AXRHe.nbin(); ++ib) {
        double trig = hissRH_trgn->GetBinContent(ib) / hissRH_trgd->GetBinContent(ib);
        double rerr = hissRH_trgn->GetBinError  (ib) / hissRH_trgn->GetBinContent(ib);
        if (!std::isfinite(trig)) { trig = 1.0; rerr = 0.0; }
        (*hRHtrig)()->SetBinContent(ib, trig);
        (*hRHtrig)()->SetBinError  (ib, trig * rerr);
    }
    (*hRHtrig)()->Fit(func_pol0, "q0", "");
    double RHtrig = func_pol0->GetParameter(0);

    // RICH - Livetime
    Hist* hRHlv = Hist::New("hRHlv", (TH1D*)(ftme->Get(Form("hRH%s_expt", sube.c_str()))));
    (*hRHlv)()->GetXaxis()->SetTitle(etitle.c_str());
    (*hRHlv)()->GetYaxis()->SetTitle("Exposure Time [s]");

    Hist* hRHcorr = Hist::New("hRHcorr", HistAxis(AXRHe, eaxis?"[(GV) m^{2} sr s]":"[(GeV/n) m^{2} sr s]"));
    for (int ib = 1; ib <= AXRHe.nbin(); ++ib) {
        double de   = (*hRHde  )()->GetBinContent(ib);
        double lv   = (*hRHlv  )()->GetBinContent(ib);
        double accp = (*hRHaccp)()->GetBinContent(ib);
        double meff = (*hRHmeff)()->GetBinContent(ib);
        double trig = (*hRHtrig)()->GetBinContent(ib);
        //double corr = de * lv * accp * meff * trig;
        double corr = de * lv * accp * RHmeff * RHtrig;

        (*hRHcorr)()->SetBinContent(ib, corr);
        (*hRHcorr)()->SetBinError  (ib, 0);
    }



    /////////////////////
    //==== OUTPUT =====//
    /////////////////////
    PdfEditor editor(Window(WindowSize::kWideSliceLR), "adflux_acc", "out");
	
    hTFde  ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTFlv  ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTFaccp->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTFmeff->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTFtrig->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTFcorr->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    hRHde  ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRHlv  ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRHaccp->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRHmeff->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRHtrig->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRHcorr->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    std::cerr << Form("\n<TF>  meff %6.3f trig %6.3f\n\n", TFmeff, TFtrig);
    std::cerr << Form("\n<RH>  meff %6.3f trig %6.3f\n\n", RHmeff, RHtrig);

    Hist* hcanvas_TFtrig = Hist::New("hcanvas_TFtrig",
            HistAxis(hTFtrig->xaxis(), Axis("", 1000, 0.0, 1.05)));
    (*hcanvas_TFtrig)()->GetXaxis()->SetTitle(etitle.c_str());
    (*hcanvas_TFtrig)()->GetYaxis()->SetTitle("Trigger Efficiency");
    
    Hist* hcanvas_TFmeff = Hist::New("hcanvas_TFmeff",
            HistAxis(hTFmeff->xaxis(), Axis("", 1000, 0.0, 1.05)));
    (*hcanvas_TFmeff)()->GetXaxis()->SetTitle(etitle.c_str());
    (*hcanvas_TFmeff)()->GetYaxis()->SetTitle("Mass Selection Efficiency");

    editor.create("Low Velocity Region", 2, 2);
    editor.cd(1, PadAxis(0, 0));
    hTFlv->draw("hist");
    editor.cd(2, PadAxis(0, 0));
    hTFaccp->draw("hist");
    editor.cd(3, PadAxis(0, 0));
    hcanvas_TFtrig->draw();
    hTFtrig->draw("pe same");
    editor.cd(4, PadAxis(0, 0));
    hcanvas_TFmeff->draw();
    hTFmeff->draw("pe same");
    editor.save();
    
    editor.create("Low Velocity Region", 2, 1);
    editor.cd(1, PadAxis(0, 0));
    hTFlv->draw("hist");
    editor.cd(2, PadAxis(0, 0));
    hTFaccp->draw("hist");
    editor.save();
    
    editor.create("Low Velocity Region (Livetime)");
    editor.cd(0, PadAxis(0, 0));
    hTFlv->draw("hist");
    editor.save();
    
    editor.create("Low Velocity Region (Acceptance)");
    editor.cd(0, PadAxis(0, 0));
    hTFaccp->draw("hist");
    editor.save();
    
    Hist* hcanvas_RHtrig = Hist::New("hcanvas_RHtrig",
            HistAxis(hRHtrig->xaxis(), Axis("", 1000, 0.0, 1.05)));
    (*hcanvas_RHtrig)()->GetXaxis()->SetTitle(etitle.c_str());
    (*hcanvas_RHtrig)()->GetYaxis()->SetTitle("Trigger Efficiency");
    
    Hist* hcanvas_RHmeff = Hist::New("hcanvas_RHmeff",
            HistAxis(hRHmeff->xaxis(), Axis("", 1000, 0.0, 1.05)));
    (*hcanvas_RHmeff)()->GetXaxis()->SetTitle(etitle.c_str());
    (*hcanvas_RHmeff)()->GetYaxis()->SetTitle("Mass Selection Efficiency");
    
    editor.create("High Velocity Region", 2, 2);
    editor.cd(1, PadAxis(0, 0));
    hRHlv->draw("hist");
    editor.cd(2, PadAxis(0, 0));
    hRHaccp->draw("hist");
    editor.cd(3, PadAxis(0, 0));
    hcanvas_RHtrig->draw();
    hRHtrig->draw("pe same");
    editor.cd(4, PadAxis(0, 0));
    hcanvas_RHmeff->draw();
    hRHmeff->draw("pe same");
    editor.save();

    editor.create("High Velocity Region", 2, 1);
    editor.cd(1, PadAxis(0, 0));
    hRHlv->draw("hist");
    editor.cd(2, PadAxis(0, 0));
    hRHaccp->draw("hist");
    editor.save();
    
    editor.create("High Velocity Region (Livetime)");
    editor.cd(0, PadAxis(0, 0));
    hRHlv->draw("hist");
    editor.save();
    
    editor.create("High Velocity Region (Acceptance)");
    editor.cd(0, PadAxis(0, 0));
    hRHaccp->draw("hist");
    editor.save();
    
    editor.close();

    TFile* ofle = new TFile("out/adflux_acc.root", "RECREATE");
    ofle->cd();

    hTFde  ->write();
    hTFlv  ->write();
    hTFaccp->write();
    hTFmeff->write();
    hTFtrig->write();
    hTFcorr->write();
    
    hRHde  ->write();
    hRHlv  ->write();
    hRHaccp->write();
    hRHmeff->write();
    hRHtrig->write();
    hRHcorr->write();

    ofle->Write();
    ofle->Close();

    return 1;
}
