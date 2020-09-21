#include <CPPLibs.h>
#include <ROOTLibs.h>

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    const double TWO_PI = 2.0 * 3.1415927;

    TFile* ftme_psel = TFile::Open("others/ElectronPositron_prl.root");
    TFile* ftme_appr = TFile::Open("out/apflux_tme.root");

    Hist* hTR_psel = Hist::New("hTR_psel", (TH2D*)(ftme_psel->Get("hh_Ratio_StatErr")));
    const Axis& AXtme_psel = hTR_psel->xaxis();
    const Axis& AXrig_psel = hTR_psel->yaxis();
    
    std::vector<double> BINke_psel;
    for (int i = 0; i <= AXrig_psel.nbin(); ++i) {
        double ke = std::hypot(AXrig_psel(i), 0.000510999) - 0.000510999;
        BINke_psel.push_back(ke);
    }
    Axis AXke_psel("KE [GeV]", BINke_psel);
    
    Hist* hTR_appr = Hist::New("hTR_appr", (TH2D*)(ftme_appr->Get("hTAppStat")));
    const Axis& AXtme_appr = hTR_appr->xaxis();
    const Axis& AXrig_appr = hTR_appr->yaxis();

    std::vector<double> BINke_appr;
    for (int i = 0; i <= AXrig_appr.nbin(); ++i) {
        double ke = std::hypot(AXrig_appr(i), 0.938272297) - 0.938272297;
        BINke_appr.push_back(ke);
    }
    Axis AXke_appr("KE [GeV]", BINke_appr);


    const double SAT_UTIME = AXtme_appr.min();
    const double END_UTIME = AXtme_appr.max();
    const double DAY  = 86400.0;
    const double YEAR = 86400.0 * 365.25;

    const double REV_UTIME = 1372636800.0;
    double POLARITY_REVERSAL = (REV_UTIME - SAT_UTIME) / YEAR;
   
    const double DEF_AMP    = 0.2;
    const double DEF_REVTME = 0.0;
    const double DEF_PERIOD = 8.0;
    
    PdfEditor editor(Window(), "phys_tme", "out");

    // Positron / Electron
    TF1* tfmen_psel = new TF1("tfmen_psel", "[0]", SAT_UTIME, END_UTIME);
    TF1* tfunc_psel = new TF1("tfunc_psel", Form("[0] / (1.0 + [1] * cos(%14.8f * (((x - %f) / %f) - [2]) / [3]))", TWO_PI, REV_UTIME, YEAR), SAT_UTIME, END_UTIME);
    tfunc_psel->SetNpx(10000);
    tfunc_psel->SetParLimits(1,  0.0, 0.8);
    tfunc_psel->SetParLimits(2, -2.5, 2.5);
    tfunc_psel->SetParLimits(3,  5.0, 15.0);
    tfunc_psel->SetLineColor(kBlue);

    Hist* hTpsel_amp    = Hist::New("hTpsel_amp",    HistAxis(AXrig_psel));
    Hist* hTpsel_revtme = Hist::New("hTpsel_revtme", HistAxis(AXrig_psel));
    Hist* hTpsel_reftme = Hist::New("hTpsel_reftme", HistAxis(AXrig_psel));
    Hist* hTpsel_period = Hist::New("hTpsel_period", HistAxis(AXrig_psel));
    
    Hist* hTpselKE_amp    = Hist::New("hTpselKE_amp",    HistAxis(AXke_psel));
    Hist* hTpselKE_revtme = Hist::New("hTpselKE_revtme", HistAxis(AXke_psel));
    Hist* hTpselKE_reftme = Hist::New("hTpselKE_reftme", HistAxis(AXke_psel));
    Hist* hTpselKE_period = Hist::New("hTpselKE_period", HistAxis(AXke_psel));
    
    std::vector<Hist*> vhR_psel = Hist::ProjectAll(HistProj::kX, hTR_psel);
    for (int ir = 1; ir <= AXrig_psel.nbin(); ++ir) {
        Hist* hpsel     = vhR_psel.at(ir);
        TH1D* TH1Dhpsel = (TH1D*)((*hpsel)());

        tfmen_psel->SetParameter(0, TH1Dhpsel->GetBinContent(1));
        (*hpsel)()->Fit(tfmen_psel, "q0", "", AXtme_psel.min(), AXtme_psel.max());
        (*hpsel)()->Fit(tfmen_psel, "q0", "", AXtme_psel.min(), AXtme_psel.max());

        tfunc_psel->SetParameters(tfmen_psel->GetParameter(0), DEF_AMP, DEF_REVTME, DEF_PERIOD);
        (*hpsel)()->Fit(tfunc_psel, "q0", "", AXtme_psel.min(), AXtme_psel.max());
        (*hpsel)()->Fit(tfunc_psel, "q0", "", AXtme_psel.min(), AXtme_psel.max());
        (*hpsel)()->Fit(tfunc_psel, "q0", "", AXtme_psel.min(), AXtme_psel.max());

        double reftme_val = tfunc_psel->GetParameter(2) + 0.25 * tfunc_psel->GetParameter(3);
        double reftme_err = std::hypot(tfunc_psel->GetParError (2), 0.25 * tfunc_psel->GetParError (3));

        ((TH1D*)((*hTpsel_amp)()))->SetBinContent(ir, tfunc_psel->GetParameter(1));
        ((TH1D*)((*hTpsel_amp)()))->SetBinError  (ir, tfunc_psel->GetParError (1));
        
        ((TH1D*)((*hTpsel_revtme)()))->SetBinContent(ir, tfunc_psel->GetParameter(2));
        ((TH1D*)((*hTpsel_revtme)()))->SetBinError  (ir, tfunc_psel->GetParError (2));
       
        ((TH1D*)((*hTpsel_reftme)()))->SetBinContent(ir, reftme_val);
        ((TH1D*)((*hTpsel_reftme)()))->SetBinError  (ir, reftme_err);
        
        ((TH1D*)((*hTpsel_period)()))->SetBinContent(ir, tfunc_psel->GetParameter(3));
        ((TH1D*)((*hTpsel_period)()))->SetBinError  (ir, tfunc_psel->GetParError (3));
		
        ((TH1D*)((*hTpselKE_amp)()))->SetBinContent(ir, tfunc_psel->GetParameter(1));
        ((TH1D*)((*hTpselKE_amp)()))->SetBinError  (ir, tfunc_psel->GetParError (1));
        
        ((TH1D*)((*hTpselKE_revtme)()))->SetBinContent(ir, tfunc_psel->GetParameter(2));
        ((TH1D*)((*hTpselKE_revtme)()))->SetBinError  (ir, tfunc_psel->GetParError (2));
        
        ((TH1D*)((*hTpselKE_reftme)()))->SetBinContent(ir, reftme_val);
        ((TH1D*)((*hTpselKE_reftme)()))->SetBinError  (ir, reftme_err);
        
        ((TH1D*)((*hTpselKE_period)()))->SetBinContent(ir, tfunc_psel->GetParameter(3));
        ((TH1D*)((*hTpselKE_period)()))->SetBinError  (ir, tfunc_psel->GetParError (3));
		
        hpsel->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
        (*hpsel)()->GetYaxis()->SetTitle("e^{+}/e^{-} flux ratio");
        (*hpsel)()->GetXaxis()->SetTitle("");
        (*hpsel)()->GetXaxis()->SetTimeDisplay(true);
        (*hpsel)()->GetXaxis()->SetTimeOffset(0, "GMT");
        (*hpsel)()->GetXaxis()->SetTimeFormat("%Y/%m/%d");

        editor.create();
        TH1Dhpsel->Draw("pe");
        tfunc_psel->Draw("l same");
        TitleDraw(Form("%.2f - %.2f [GV]", AXrig_psel(ir-1), AXrig_psel(ir)));
        editor.save();
    }
    
    hTpsel_amp   ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hTpsel_revtme->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hTpsel_reftme->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hTpsel_period->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    
    hTpselKE_amp   ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hTpselKE_revtme->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hTpselKE_reftme->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hTpselKE_period->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    

    // Antiproton / Proton
    TF1* tfmen_appr = new TF1("tfmen_appr", "[0]", SAT_UTIME, END_UTIME);
    TF1* tfunc_appr = new TF1("tfunc_appr", Form("[0] * (1.0 + [1] * cos(%14.8f * (((x - %f) / %f) - [2]) / [3]))", TWO_PI, REV_UTIME, YEAR), SAT_UTIME, END_UTIME);
    tfunc_appr->SetNpx(10000);
    tfunc_appr->SetParLimits(1,  0.0, 0.8);
    tfunc_appr->SetParLimits(2, -2.5, 2.5);
    tfunc_appr->SetParLimits(3,  5.0, 15.0);
    tfunc_appr->SetLineColor(kBlue);

    Hist* hTappr_amp    = Hist::New("hTappr_amp",    HistAxis(AXrig_appr));
    Hist* hTappr_revtme = Hist::New("hTappr_revtme", HistAxis(AXrig_appr));
    Hist* hTappr_reftme = Hist::New("hTappr_reftme", HistAxis(AXrig_appr));
    Hist* hTappr_period = Hist::New("hTappr_period", HistAxis(AXrig_appr));
    
    Hist* hTapprKE_amp    = Hist::New("hTapprKE_amp",    HistAxis(AXke_appr));
    Hist* hTapprKE_revtme = Hist::New("hTapprKE_revtme", HistAxis(AXke_appr));
    Hist* hTapprKE_reftme = Hist::New("hTapprKE_reftme", HistAxis(AXke_appr));
    Hist* hTapprKE_period = Hist::New("hTapprKE_period", HistAxis(AXke_appr));
    
    std::vector<Hist*> vhR_appr = Hist::ProjectAll(HistProj::kX, hTR_appr);
    for (int ir = 1; ir <= AXrig_appr.nbin(); ++ir) {
        Hist* happr     = vhR_appr.at(ir);
        TH1D* TH1Dhappr = (TH1D*)((*happr)());
        
        tfmen_appr->SetParameter(0, TH1Dhappr->GetBinContent(1));
        (*happr)()->Fit(tfmen_appr, "q0", "", AXtme_appr.min(), AXtme_appr.max());
        (*happr)()->Fit(tfmen_appr, "q0", "", AXtme_appr.min(), AXtme_appr.max());
        
        tfunc_appr->SetParameters(tfmen_appr->GetParameter(0), DEF_AMP, DEF_REVTME, DEF_PERIOD);
        (*happr)()->Fit(tfunc_appr, "q0", "", AXtme_appr.min(), AXtme_appr.max());
        (*happr)()->Fit(tfunc_appr, "q0", "", AXtme_appr.min(), AXtme_appr.max());
        (*happr)()->Fit(tfunc_appr, "q0", "", AXtme_appr.min(), AXtme_appr.max());
        
        double reftme_val = tfunc_appr->GetParameter(2) + 0.25 * tfunc_appr->GetParameter(3);
        double reftme_err = std::hypot(tfunc_appr->GetParError (2), 0.25 * tfunc_appr->GetParError (3));

        ((TH1D*)((*hTappr_amp)()))->SetBinContent(ir, tfunc_appr->GetParameter(1));
        ((TH1D*)((*hTappr_amp)()))->SetBinError  (ir, tfunc_appr->GetParError (1));
        
        ((TH1D*)((*hTappr_revtme)()))->SetBinContent(ir, tfunc_appr->GetParameter(2));
        ((TH1D*)((*hTappr_revtme)()))->SetBinError  (ir, tfunc_appr->GetParError (2));
        
        ((TH1D*)((*hTappr_reftme)()))->SetBinContent(ir, reftme_val);
        ((TH1D*)((*hTappr_reftme)()))->SetBinError  (ir, reftme_err);
        
        ((TH1D*)((*hTappr_period)()))->SetBinContent(ir, tfunc_appr->GetParameter(3));
        ((TH1D*)((*hTappr_period)()))->SetBinError  (ir, tfunc_appr->GetParError (3));
        
        ((TH1D*)((*hTapprKE_amp)()))->SetBinContent(ir, tfunc_appr->GetParameter(1));
        ((TH1D*)((*hTapprKE_amp)()))->SetBinError  (ir, tfunc_appr->GetParError (1));
        
        ((TH1D*)((*hTapprKE_revtme)()))->SetBinContent(ir, tfunc_appr->GetParameter(2));
        ((TH1D*)((*hTapprKE_revtme)()))->SetBinError  (ir, tfunc_appr->GetParError (2));
        
        ((TH1D*)((*hTapprKE_reftme)()))->SetBinContent(ir, reftme_val);
        ((TH1D*)((*hTapprKE_reftme)()))->SetBinError  (ir, reftme_err);
        
        ((TH1D*)((*hTapprKE_period)()))->SetBinContent(ir, tfunc_appr->GetParameter(3));
        ((TH1D*)((*hTapprKE_period)()))->SetBinError  (ir, tfunc_appr->GetParError (3));
        
        happr->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
        (*happr)()->GetYaxis()->SetTitle("#bar{p}/p flux ratio");
        (*happr)()->GetXaxis()->SetTitle("");
        (*happr)()->GetXaxis()->SetTimeDisplay(true);
        (*happr)()->GetXaxis()->SetTimeOffset(0, "GMT");
        (*happr)()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
        
        editor.create();
        TH1Dhappr->Draw("pe");
        tfunc_appr->Draw("l same");
        TitleDraw(Form("%.2f - %.2f [GV]", AXrig_appr(ir-1), AXrig_appr(ir)));
        editor.save();
    }
        
    hTappr_amp   ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTappr_revtme->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTappr_reftme->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTappr_period->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    hTapprKE_amp   ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTapprKE_revtme->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTapprKE_reftme->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hTapprKE_period->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    // Summary
    Axis AXamp("Amp", 2000, -0.02, 0.6);
    Hist* hAmp = Hist::New("hAmp", HistAxis(AXrig_appr, AXamp));
    (*hAmp)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hAmp)()->GetXaxis()->SetMoreLogLabels();
    (*hAmp)()->GetXaxis()->CenterTitle();
    (*hAmp)()->GetYaxis()->CenterTitle();
    (*hAmp)()->GetXaxis()->SetRangeUser(1.0, 25.0);
    
    Hist* hKEAmp = Hist::New("hKEAmp", HistAxis(AXke_appr, AXamp));
    (*hKEAmp)()->GetXaxis()->SetTitle("KE [GeV]");
    (*hKEAmp)()->GetXaxis()->SetMoreLogLabels();
    (*hKEAmp)()->GetXaxis()->CenterTitle();
    (*hKEAmp)()->GetYaxis()->CenterTitle();
    (*hKEAmp)()->GetXaxis()->SetRangeUser(0.4, 25.0);
    
    editor.create();
    editor.cd(1, PadAxis(1));
    hAmp->draw();
    hTpsel_amp->draw("pe same");
    hTappr_amp->draw("pe same");
    Legend leg_amp("", TextStyle(kBlack, 30, 43), PadWindow(0.65, 0.90, 0.75, 0.90));
    leg_amp()->SetHeader();
    leg_amp()->AddEntry((*hTpsel_amp)(), "e^{+}/e^{-} flux ratio", "lp");
    leg_amp()->AddEntry((*hTappr_amp)(), "#bar{p}/p flux ratio", "lp");
    leg_amp()->SetFillColor(0);
    leg_amp.draw();
    editor.save();
    
    editor.create();
    editor.cd(1, PadAxis(1));
    hKEAmp->draw();
    hTpselKE_amp->draw("pe same");
    hTapprKE_amp->draw("pe same");
    leg_amp.draw();
    editor.save();
    
    Axis AXrevtme("t_{#phi=0}-t_{rev} (years)", 2000, -4, 2);
    Hist* hRevtme = Hist::New("hRevtme", HistAxis(AXrig_appr, AXrevtme));
    (*hRevtme)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hRevtme)()->GetXaxis()->SetMoreLogLabels();
    (*hRevtme)()->GetXaxis()->CenterTitle();
    (*hRevtme)()->GetYaxis()->CenterTitle();
    (*hRevtme)()->GetXaxis()->SetRangeUser(1.0, 4.5);
    (*hRevtme)()->GetXaxis()->SetRangeUser(1.0, 7.0);
    
    Hist* hKERevtme = Hist::New("hKERevtme", HistAxis(AXke_appr, AXrevtme));
    (*hKERevtme)()->GetXaxis()->SetTitle("KE [GeV]");
    (*hKERevtme)()->GetXaxis()->SetMoreLogLabels();
    (*hKERevtme)()->GetXaxis()->CenterTitle();
    (*hKERevtme)()->GetYaxis()->CenterTitle();
    (*hKERevtme)()->GetXaxis()->SetRangeUser(0.4, 4.0);
    
    editor.create();
    editor.cd(1, PadAxis(1));
    hRevtme->draw();
    hTpsel_revtme->draw("pe same");
    hTappr_revtme->draw("pe same");
    Legend leg_revtme("", TextStyle(kBlack, 30, 43), PadWindow(0.65, 0.90, 0.75, 0.90));
    leg_revtme()->SetHeader();
    leg_revtme()->AddEntry((*hTpsel_revtme)(), "e^{+}/e^{-} flux ratio", "lp");
    leg_revtme()->AddEntry((*hTappr_revtme)(), "#bar{p}/p flux ratio", "lp");
    leg_revtme()->SetFillColor(0);
    leg_revtme.draw();
    editor.save();
    
    editor.create();
    editor.cd(1, PadAxis(1));
    hKERevtme->draw();
    hTpselKE_revtme->draw("pe same");
    hTapprKE_revtme->draw("pe same");
    leg_revtme.draw();
    editor.save();
    
    Axis AXreftme("t_{#phi=#pi/2}-t_{ref} (years)", 2000, 0, 5);
    Hist* hReftme = Hist::New("hReftme", HistAxis(AXrig_appr, AXreftme));
    (*hReftme)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hReftme)()->GetXaxis()->SetMoreLogLabels();
    (*hReftme)()->GetXaxis()->CenterTitle();
    (*hReftme)()->GetYaxis()->CenterTitle();
    (*hReftme)()->GetXaxis()->SetRangeUser(1.0, 4.5);
    (*hReftme)()->GetXaxis()->SetRangeUser(1.0, 7.0);
    
    Hist* hKEReftme = Hist::New("hKEReftme", HistAxis(AXke_appr, AXreftme));
    (*hKEReftme)()->GetXaxis()->SetTitle("KE [GeV]");
    (*hKEReftme)()->GetXaxis()->SetMoreLogLabels();
    (*hKEReftme)()->GetXaxis()->CenterTitle();
    (*hKEReftme)()->GetYaxis()->CenterTitle();
    (*hKEReftme)()->GetXaxis()->SetRangeUser(0.4, 4.0);
    
    editor.create();
    editor.cd(1, PadAxis(1));
    hReftme->draw();
    hTpsel_reftme->draw("pe same");
    hTappr_reftme->draw("pe same");
    Legend leg_reftme("", TextStyle(kBlack, 30, 43), PadWindow(0.65, 0.90, 0.75, 0.90));
    leg_reftme()->SetHeader();
    leg_reftme()->AddEntry((*hTpsel_reftme)(), "e^{+}/e^{-} flux ratio", "lp");
    leg_reftme()->AddEntry((*hTappr_reftme)(), "#bar{p}/p flux ratio", "lp");
    leg_reftme()->SetFillColor(0);
    leg_reftme.draw();
    editor.save();
    
    Axis AXperiod("T (years)", 2000, 4, 14);
    Hist* hPeriod = Hist::New("hPeriod", HistAxis(AXrig_appr, AXperiod));
    (*hPeriod)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hPeriod)()->GetXaxis()->SetMoreLogLabels();
    (*hPeriod)()->GetXaxis()->CenterTitle();
    (*hPeriod)()->GetYaxis()->CenterTitle();
    (*hPeriod)()->GetXaxis()->SetRangeUser(1.0, 4.5);
    (*hPeriod)()->GetXaxis()->SetRangeUser(1.0, 7.0);
    
    Hist* hKEPeriod = Hist::New("hKEPeriod", HistAxis(AXke_appr, AXperiod));
    (*hKEPeriod)()->GetXaxis()->SetTitle("KE [GeV]");
    (*hKEPeriod)()->GetXaxis()->SetMoreLogLabels();
    (*hKEPeriod)()->GetXaxis()->CenterTitle();
    (*hKEPeriod)()->GetYaxis()->CenterTitle();
    (*hKEPeriod)()->GetXaxis()->SetRangeUser(0.4, 4.0);
    
    editor.create();
    editor.cd(1, PadAxis(1));
    hPeriod->draw();
    hTpsel_period->draw("pe same");
    hTappr_period->draw("pe same");
    Legend leg_period("", TextStyle(kBlack, 30, 43), PadWindow(0.65, 0.90, 0.75, 0.90));
    leg_period()->SetHeader();
    leg_period()->AddEntry((*hTpsel_period)(), "e^{+}/e^{-} flux ratio", "lp");
    leg_period()->AddEntry((*hTappr_period)(), "#bar{p}/p flux ratio", "lp");
    leg_period()->SetFillColor(0);
    leg_period.draw();
    editor.save();
    
    editor.create();
    editor.cd(1, PadAxis(1));
    hKEPeriod->draw();
    hTpselKE_period->draw("pe same");
    hTapprKE_period->draw("pe same");
    leg_period.draw();
    editor.save();
    
    editor.close();
    
    TFile * ofle = new TFile("out/phys_tme.root", "RECREATE");
    ofle->cd();

    (*hTpsel_amp   )()->Write();
    (*hTpsel_revtme)()->Write();
    (*hTpsel_period)()->Write();
    
    (*hTappr_amp   )()->Write();
    (*hTappr_revtme)()->Write();
    (*hTappr_period)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
