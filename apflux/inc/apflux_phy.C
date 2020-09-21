#include <CPPLibs.h>
#include <ROOTLibs.h>

static constexpr double TWO_PI = 2.0 * 3.1415927;
static constexpr double HALF_PI = 0.5 * 3.1415927;

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subt = "3M";

    short fitType = 1; // 0, +Z/-Z  1, -Z/+Z  flux ratio

    const double YEAR2DAY = 1.0;
    //const double YEAR2DAY = 365.25;
    TFile* felps = TFile::Open(Form("out/psphys_tme_%s.root", fitType?"np":"pn"));
    TFile* fappr = TFile::Open(Form("out/apphys_tme%s_%s.root", subt.c_str(), fitType?"np":"pn"));
    
    // Logistic function
    Hist* hRIG_elpsL_amp    = Hist::New("hRIG_elpsL_amp"   , (TH1D*)(felps->Get("hRIG_L_amp"   )));
    Hist* hRIG_elpsL_delay  = Hist::New("hRIG_elpsL_delay" , (TH1D*)(felps->Get("hRIG_L_delay" )));
    Hist* hRIG_elpsL_period = Hist::New("hRIG_elpsL_period", (TH1D*)(felps->Get("hRIG_L_period")));

    Hist* hRIG_apprL_amp    = Hist::New("hRIG_apprL_amp"   , (TH1D*)(fappr->Get("hRIG_L_amp"   )));
    Hist* hRIG_apprL_delay  = Hist::New("hRIG_apprL_delay" , (TH1D*)(fappr->Get("hRIG_L_delay" )));
    Hist* hRIG_apprL_period = Hist::New("hRIG_apprL_period", (TH1D*)(fappr->Get("hRIG_L_period")));
    
    Hist* hKE_elpsL_amp    = Hist::New("hKE_elpsL_amp"   , (TH1D*)(felps->Get("hKE_L_amp"   )));
    Hist* hKE_elpsL_delay  = Hist::New("hKE_elpsL_delay" , (TH1D*)(felps->Get("hKE_L_delay" )));
    Hist* hKE_elpsL_period = Hist::New("hKE_elpsL_period", (TH1D*)(felps->Get("hKE_L_period")));

    Hist* hKE_apprL_amp    = Hist::New("hKE_apprL_amp"   , (TH1D*)(fappr->Get("hKE_L_amp"   )));
    Hist* hKE_apprL_delay  = Hist::New("hKE_apprL_delay" , (TH1D*)(fappr->Get("hKE_L_delay" )));
    Hist* hKE_apprL_period = Hist::New("hKE_apprL_period", (TH1D*)(fappr->Get("hKE_L_period")));
    
    hRIG_elpsL_amp   ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hRIG_elpsL_delay ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hRIG_elpsL_period->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hKE_elpsL_amp    ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hKE_elpsL_delay  ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hKE_elpsL_period ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    
    hRIG_apprL_amp   ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRIG_apprL_delay ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRIG_apprL_period->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_apprL_amp    ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_apprL_delay  ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_apprL_period ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
   
    // Harmonic function
    Hist* hRIG_elpsH_amp    = Hist::New("hRIG_elpsH_amp"   , (TH1D*)(felps->Get("hRIG_H_amp"   )));
    Hist* hRIG_elpsH_delay  = Hist::New("hRIG_elpsH_delay" , (TH1D*)(felps->Get("hRIG_H_delay" )));
    Hist* hRIG_elpsH_period = Hist::New("hRIG_elpsH_period", (TH1D*)(felps->Get("hRIG_H_period")));

    Hist* hRIG_apprH_amp    = Hist::New("hRIG_apprH_amp"   , (TH1D*)(fappr->Get("hRIG_H_amp"   )));
    Hist* hRIG_apprH_delay  = Hist::New("hRIG_apprH_delay" , (TH1D*)(fappr->Get("hRIG_H_delay" )));
    Hist* hRIG_apprH_period = Hist::New("hRIG_apprH_period", (TH1D*)(fappr->Get("hRIG_H_period")));
    
    Hist* hKE_elpsH_amp    = Hist::New("hKE_elpsH_amp"   , (TH1D*)(felps->Get("hKE_H_amp"   )));
    Hist* hKE_elpsH_delay  = Hist::New("hKE_elpsH_delay" , (TH1D*)(felps->Get("hKE_H_delay" )));
    Hist* hKE_elpsH_period = Hist::New("hKE_elpsH_period", (TH1D*)(felps->Get("hKE_H_period")));

    Hist* hKE_apprH_amp    = Hist::New("hKE_apprH_amp"   , (TH1D*)(fappr->Get("hKE_H_amp"   )));
    Hist* hKE_apprH_delay  = Hist::New("hKE_apprH_delay" , (TH1D*)(fappr->Get("hKE_H_delay" )));
    Hist* hKE_apprH_period = Hist::New("hKE_apprH_period", (TH1D*)(fappr->Get("hKE_H_period")));
    
    hRIG_elpsH_amp   ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hRIG_elpsH_delay ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hRIG_elpsH_period->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hKE_elpsH_amp    ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hKE_elpsH_delay  ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hKE_elpsH_period ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    
    hRIG_apprH_amp   ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRIG_apprH_delay ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRIG_apprH_period->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_apprH_amp    ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_apprH_delay  ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_apprH_period ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    Hist* hRIG_Lshift = Hist::New("hRIG_Lshift", HistAxis(hRIG_apprL_amp->xaxis(), "#Deltat_{1/2}"));
    Hist* hRIG_Hshift = Hist::New("hRIG_Hshift", HistAxis(hRIG_apprL_amp->xaxis(), "#Deltat_{1/2}"));
    hRIG_Lshift->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(fitType?MarkerShape::kCircle:MarkerShape::kSquare)));
    hRIG_Hshift->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(fitType?MarkerShape::kCircle:MarkerShape::kSquare)));

    //PdfEditor editor(Window(WindowSize::kSliceLR), Form("apflux_phy%s", subt.c_str()), "out");
    PdfEditor editor(Window(WindowSize::kWideSliceLR), Form("apflux_phy%s_%s", subt.c_str(), fitType?"np":"pn"), "out");
    //PdfEditor editor(Window(600, 300), Form("apflux_phy%s_%s", subt.c_str(), fitType?"np":"pn"), "out");

    // Axis
    Axis AXrig("|Rigidity| [GV]", 5000, hRIG_apprL_amp->xaxis().min(), hRIG_apprL_amp->xaxis().max(), AxisScale::kLog);
    Axis AXke("Kinetic Energy [GeV]", 5000, hKE_apprL_amp->xaxis().min(), hKE_apprL_amp->xaxis().max(), AxisScale::kLog);
   
    // PRL
    TF1* baseline = new TF1("baseline", "[0]", 0.1, 100);
    baseline->SetNpx(10000);
    baseline->SetLineColor(kBlack);
    baseline->SetLineStyle(3);
    baseline->SetLineWidth(2.0);
    baseline->SetParameter(0, 0.0);
    
    TF1* elps_prl_delay = new TF1("elps_prl_delay", "[0]*pow(x, [1])", 0.1, 100);
    elps_prl_delay->SetNpx(10000);
    elps_prl_delay->SetLineColor(kGreen+2);
    elps_prl_delay->SetLineStyle(1);
    elps_prl_delay->SetLineWidth(2.0);
    elps_prl_delay->SetParameters(580.0, -0.33);
    elps_prl_delay->SetParameters(580.0 / 365.25, -0.33); // YEAR2DAY
    
    TF1* elps_prl_period = new TF1("elps_prl_period", "[0]", 0.1, 100);
    elps_prl_period->SetNpx(10000);
    elps_prl_period->SetLineColor(kGreen+2);
    elps_prl_period->SetLineStyle(1);
    elps_prl_period->SetLineWidth(2.0);
    elps_prl_period->SetParameter(0, 830.0);
    elps_prl_period->SetParameter(0, 830.0 / 365.25); // YEAR2DAY
    
    // Delay
    TF1* Lelps_delay = new TF1("Lelps_delay", "[0]*pow(x, [1])", 0.1, 100);
    Lelps_delay->SetNpx(10000);
    Lelps_delay->SetLineColor(kGreen+2);
    Lelps_delay->SetLineStyle(1);
    Lelps_delay->SetLineWidth(2.0);
    
    TF1* Lappr_delay = new TF1("Lappr_delay", "[0]*pow(x, [1])", 0.1, 100);
    Lappr_delay->SetNpx(10000);
    Lappr_delay->SetLineColor(kRed);
    Lappr_delay->SetLineStyle(1);
    Lappr_delay->SetLineWidth(2.0);

    Lelps_delay->SetParameters(1.0, -0.3);
    (*hRIG_elpsL_delay)()->Fit(Lelps_delay, "q0", "", 1.0, 7.5);
    std::array<double, 4> elpsL_delay({
            Lelps_delay->GetParameter(0), Lelps_delay->GetParError(0),
            Lelps_delay->GetParameter(1), Lelps_delay->GetParError(1)});
    
    Lappr_delay->SetParameters(1.0, -0.3);
    (*hRIG_apprL_delay)()->Fit(Lappr_delay, "q0", "", 1.0, 7.5);
    std::array<double, 4> apprL_delay({
            Lappr_delay->GetParameter(0), Lappr_delay->GetParError(0),
            Lappr_delay->GetParameter(1), Lappr_delay->GetParError(1)});

    TF1* Helps_delay = new TF1("Helps_delay", "[0]*pow(x, [1])", 0.1, 100);
    Helps_delay->SetNpx(10000);
    Helps_delay->SetLineColor(kGreen+2);
    Helps_delay->SetLineStyle(1);
    Helps_delay->SetLineWidth(2.0);
    
    TF1* Happr_delay = new TF1("Happr_delay", "[0]*pow(x, [1])", 0.1, 100);
    Happr_delay->SetNpx(10000);
    Happr_delay->SetLineColor(kRed);
    Happr_delay->SetLineStyle(1);
    Happr_delay->SetLineWidth(2.0);
    
    Helps_delay->SetParameters(1.0, -0.3, 0.0);
    (*hRIG_elpsH_delay)()->Fit(Helps_delay, "q0", "", 1.0, 7.5);
    std::array<double, 4> elpsH_delay({ 
            Helps_delay->GetParameter(0), Helps_delay->GetParError(0),
            Helps_delay->GetParameter(1), Helps_delay->GetParError(1) });
    
    Happr_delay->SetParameters(1.0, -0.3, 0.0);
    (*hRIG_apprH_delay)()->Fit(Happr_delay, "q0", "", 1.0, 7.5);
    std::array<double, 4> apprH_delay({
            Happr_delay->GetParameter(0), Happr_delay->GetParError(0),
            Happr_delay->GetParameter(1), Happr_delay->GetParError(1) });
   
    // Period
    TF1* Lelps_period = new TF1("Lelps_period", "[0]", 0.1, 100);
    Lelps_period->SetNpx(10000);
    Lelps_period->SetLineColor(kGreen+2);
    Lelps_period->SetLineStyle(1);
    Lelps_period->SetLineWidth(2.0);
    
    TF1* Lappr_period = new TF1("Lappr_period", "[0]", 0.1, 100);
    Lappr_period->SetNpx(10000);
    Lappr_period->SetLineColor(kRed);
    Lappr_period->SetLineStyle(1);
    Lappr_period->SetLineWidth(2.0);
    
    (*hRIG_elpsL_period)()->Fit(Lelps_period, "q0", "", 1.0, 7.5);
    std::array<double, 2> elpsL_period({ Lelps_period->GetParameter(0), Lelps_period->GetParError(0) });
    
    (*hRIG_apprL_period)()->Fit(Lappr_period, "q0", "", 1.0, 7.5);
    std::array<double, 2> apprL_period({ Lappr_period->GetParameter(0), Lappr_period->GetParError(0) });
    
    TF1* Helps_period = new TF1("Helps_period", "[0]", 0.1, 100);
    Helps_period->SetNpx(10000);
    Helps_period->SetLineColor(kGreen+2);
    Helps_period->SetLineStyle(1);
    Helps_period->SetLineWidth(2.0);
    
    TF1* Happr_period = new TF1("Happr_period", "[0]", 0.1, 100);
    Happr_period->SetNpx(10000);
    Happr_period->SetLineColor(kRed);
    Happr_period->SetLineStyle(1);
    Happr_period->SetLineWidth(2.0);
   
    (*hRIG_elpsH_period)()->Fit(Helps_period, "q0", "", 1.0, 7.5);
    std::array<double, 2> elpsH_period({ Helps_period->GetParameter(0), Helps_period->GetParError(0) });
    
    (*hRIG_apprH_period)()->Fit(Happr_period, "q0", "", 1.0, 7.5);
    std::array<double, 2> apprH_period({ Happr_period->GetParameter(0), Happr_period->GetParError(0) });

    for (int ib = 1; ib <= hRIG_apprL_delay->xaxis().nbin(); ++ib) {
        double rig = hRIG_apprL_delay->xaxis().center(ib, AxisScale::kLog);
        double val = (*hRIG_apprL_delay)()->GetBinContent(ib);
        double err = (*hRIG_apprL_delay)()->GetBinError  (ib);
        (*hRIG_Lshift)()->SetBinContent(ib, val - Lelps_delay->Eval(rig));
        (*hRIG_Lshift)()->SetBinError  (ib, std::hypot(err, Lelps_delay->GetParError(0)));
    }
    
    for (int ib = 1; ib <= hRIG_apprH_delay->xaxis().nbin(); ++ib) {
        double rig = hRIG_apprH_delay->xaxis().center(ib, AxisScale::kLog);
        double val = (*hRIG_apprH_delay)()->GetBinContent(ib);
        double err = (*hRIG_apprH_delay)()->GetBinError  (ib);
        (*hRIG_Hshift)()->SetBinContent(ib, val - Helps_delay->Eval(rig));
        (*hRIG_Hshift)()->SetBinError  (ib, std::hypot(err, Helps_delay->GetParError(0)));
    }

    // ====================
    // ===== Logistic =====
    // ====================
    editor.create();
    TextDraw("Logistic Function", TextStyle(kBlack, 0.1), TextAlign(0.5, 0.5));
    editor.save();

    // Logistic <Amp>
    //Axis AX_Lamp("C", 2000, -0.3, 1.6);
    Axis AX_Lamp("L_{C}^{-Z/+Z}", 2000, -0.3, 1.6);
    Hist* hRIG_Lamp = Hist::New("hRIG_Lamp", HistAxis(AXrig, AX_Lamp));
    (*hRIG_Lamp)()->GetXaxis()->SetMoreLogLabels();
    (*hRIG_Lamp)()->GetXaxis()->CenterTitle();
    (*hRIG_Lamp)()->GetYaxis()->CenterTitle();
    (*hRIG_Lamp)()->GetXaxis()->SetRangeUser(1.01, 45.0);
    Hist* hKE_Lamp = Hist::New("hKE_Lamp", HistAxis(AXke, AX_Lamp));
    (*hKE_Lamp)()->GetXaxis()->SetMoreLogLabels();
    (*hKE_Lamp)()->GetXaxis()->CenterTitle();
    (*hKE_Lamp)()->GetYaxis()->CenterTitle();
    (*hKE_Lamp)()->GetXaxis()->SetRangeUser(0.5, 45.0);
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIG_Lamp)()->Draw();
    baseline->SetRange(0.1, 100);
    baseline->Draw("l same");
    (*hRIG_elpsL_amp)()->Draw("pe same");
    (*hRIG_apprL_amp)()->Draw("pe same");
    Legend legRIG_Lamp("", TextStyle(kBlack, 30, 43), PadWindow(0.3, 0.6, 0.75, 0.90));
    //legRIG_Lamp()->SetHeader("Logistic");
    legRIG_Lamp()->AddEntry((*hRIG_elpsL_amp)(), "e^{-}/e^{+} flux ratio", "lp");
    legRIG_Lamp()->AddEntry((*hRIG_apprL_amp)(), "#bar{p}/p flux ratio", "lp");
    legRIG_Lamp()->SetFillColor(0);
    legRIG_Lamp.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hKE_Lamp)()->Draw();
    baseline->SetRange(0.1, 100);
    baseline->Draw("l same");
    (*hKE_elpsL_amp)()->Draw("pe same");
    (*hKE_apprL_amp)()->Draw("pe same");
    Legend legKE_Lamp("", TextStyle(kBlack, 30, 43), PadWindow(0.3, 0.6, 0.75, 0.90));
    //legKE_Lamp()->SetHeader("Logistic");
    legKE_Lamp()->AddEntry((*hKE_elpsL_amp)(), "e^{-}/e^{+} flux ratio", "lp");
    legKE_Lamp()->AddEntry((*hKE_apprL_amp)(), "#bar{p}/p flux ratio", "lp");
    legKE_Lamp()->SetFillColor(0);
    legKE_Lamp.draw();
    editor.save();

    // Logistic <Delay>
    //Axis AX_Ldelay("t_{1/2} - t_{rev} (years)", 2000, 0.0 * YEAR2DAY, 5.0 * YEAR2DAY);
    Axis AX_Ldelay("L_{S}^{-Z/+Z} (years)", 2000, 0.0 * YEAR2DAY, 5.0 * YEAR2DAY);
    Hist* hRIG_Ldelay = Hist::New("hRIG_Ldelay", HistAxis(AXrig, AX_Ldelay));
    (*hRIG_Ldelay)()->GetXaxis()->SetMoreLogLabels();
    (*hRIG_Ldelay)()->GetXaxis()->CenterTitle();
    (*hRIG_Ldelay)()->GetYaxis()->CenterTitle();
    (*hRIG_Ldelay)()->GetXaxis()->SetRangeUser(1.01, 7.5);

    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIG_Ldelay)()->Draw();
    //powline_elps_prl->SetRange(0.1, 100.0);
    //powline_elps_prl->Draw("l same");
    Lelps_delay->Draw("l same");
    //Lappr_delay->Draw("l same");
    (*hRIG_elpsL_delay)()->Draw("pe same");
    (*hRIG_apprL_delay)()->Draw("pe same");
    Legend legRIG_Ldelay("", TextStyle(kBlack, 30, 43), PadWindow(0.3, 0.6, 0.70, 0.90));
    //legRIG_Ldelay()->SetHeader("Logistic");
    legRIG_Ldelay()->AddEntry((*hRIG_elpsL_delay)(), "e^{-}/e^{+} flux ratio", "lp");
    legRIG_Ldelay()->AddEntry((*hRIG_apprL_delay)(), "#bar{p}/p flux ratio", "lp");
    legRIG_Ldelay()->AddEntry((TObject*)0, Form("L_{S}^{e^{-}/e^{+}} #approx (%4.2f#pm%4.2f)  R^{%4.2f#pm%4.2f}", elpsL_delay[0], elpsL_delay[1], elpsL_delay[2], elpsL_delay[3]), "");
    legRIG_Ldelay()->SetFillColor(0);
    legRIG_Ldelay.draw();
    editor.save();
    
    // Logistic <Period>
    //Axis AX_Lperiod("#Delta t (years)", 2000, 0.0 * YEAR2DAY, 7.0 * YEAR2DAY);
    Axis AX_Lperiod("L_{T}^{-Z/+Z} (years)", 2000, 0.0 * YEAR2DAY, 7.0 * YEAR2DAY);
    Hist* hRIG_Lperiod = Hist::New("hRIG_Lperiod", HistAxis(AXrig, AX_Lperiod));
    (*hRIG_Lperiod)()->GetXaxis()->SetMoreLogLabels();
    (*hRIG_Lperiod)()->GetXaxis()->CenterTitle();
    (*hRIG_Lperiod)()->GetYaxis()->CenterTitle();
    (*hRIG_Lperiod)()->GetXaxis()->SetRangeUser(1.01, 7.5);
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIG_Lperiod)()->Draw();
    //Lelps_period_prl->SetRange(0.1, 100.0);
    //Lelps_period_prl->Draw("l same");
    //Lelps_period->Draw("l same");
    //Lappr_period->Draw("l same");
    (*hRIG_elpsL_period)()->Draw("pe same");
    (*hRIG_apprL_period)()->Draw("pe same");
    Legend legRIG_Lperiod("", TextStyle(kBlack, 30, 43), PadWindow(0.15, 0.50, 0.60, 0.90));
    //legRIG_Lperiod()->SetHeader("Logistic");
    legRIG_Lperiod()->AddEntry((*hRIG_elpsL_period)(), "e^{-}/e^{+} flux ratio", "lp");
    //legRIG_Lperiod()->AddEntry((TObject*)0, Form("#Delta t_{e^{-}/e^{+}} = %.2f #pm %.2f (years)", elpsL_period[0], elpsL_period[1]), "");
    legRIG_Lperiod()->AddEntry((TObject*)0, Form("L_{T}^{e^{-}/e^{+}} = %.2f #pm %.2f (years)", elpsL_period[0], elpsL_period[1]), "");
    legRIG_Lperiod()->AddEntry((*hRIG_apprL_period)(), "#bar{p}/p flux ratio", "lp");
    //legRIG_Lperiod()->AddEntry((TObject*)0, Form("#Delta t_{#bar{p}/p} = %.2f #pm %.2f (years)", apprL_period[0], apprL_period[1]), "");
    legRIG_Lperiod()->AddEntry((TObject*)0, Form("L_{T}^{#bar{p}/p} = %.2f #pm %.2f (years)", apprL_period[0], apprL_period[1]), "");
    legRIG_Lperiod()->SetFillColor(0);
    legRIG_Lperiod.draw();
    editor.save();
    
    // ====================
    // ===== Harmonic =====
    // ====================
    editor.create();
    TextDraw("Harmonic Function", TextStyle(kBlack, 0.1), TextAlign(0.5, 0.5));
    editor.save();
    
    // Harmonic <Amp>
    //Axis AX_Hamp("C", 2000, -0.1, 0.6);
    Axis AX_Hamp("H_{C}^{-Z/+Z}", 2000, -0.1, 0.6);
    Hist* hRIG_Hamp = Hist::New("hRIG_Hamp", HistAxis(AXrig, AX_Hamp));
    (*hRIG_Hamp)()->GetXaxis()->SetMoreLogLabels();
    (*hRIG_Hamp)()->GetXaxis()->CenterTitle();
    (*hRIG_Hamp)()->GetYaxis()->CenterTitle();
    (*hRIG_Hamp)()->GetXaxis()->SetRangeUser(1.01, 45.0);
    Hist* hKE_Hamp = Hist::New("hKE_Hamp", HistAxis(AXke, AX_Hamp));
    (*hKE_Hamp)()->GetXaxis()->SetMoreLogLabels();
    (*hKE_Hamp)()->GetXaxis()->CenterTitle();
    (*hKE_Hamp)()->GetYaxis()->CenterTitle();
    (*hKE_Hamp)()->GetXaxis()->SetRangeUser(0.5, 45.0);
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIG_Hamp)()->Draw();
    baseline->SetRange(0.1, 100);
    baseline->Draw("l same");
    (*hRIG_elpsH_amp)()->Draw("pe same");
    (*hRIG_apprH_amp)()->Draw("pe same");
    Legend legRIG_Hamp("", TextStyle(kBlack, 30, 43), PadWindow(0.3, 0.6, 0.75, 0.90));
    //legRIG_Hamp()->SetHeader("Harmonic");
    legRIG_Hamp()->AddEntry((*hRIG_elpsH_amp)(), "e^{-}/e^{+} flux ratio", "lp");
    legRIG_Hamp()->AddEntry((*hRIG_apprH_amp)(), "#bar{p}/p flux ratio", "lp");
    legRIG_Hamp()->SetFillColor(0);
    legRIG_Hamp.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hKE_Hamp)()->Draw();
    baseline->SetRange(0.1, 100);
    baseline->Draw("l same");
    (*hKE_elpsH_amp)()->Draw("pe same");
    (*hKE_apprH_amp)()->Draw("pe same");
    Legend legKE_Hamp("", TextStyle(kBlack, 30, 43), PadWindow(0.3, 0.6, 0.75, 0.90));
    //legKE_Hamp()->SetHeader("Harmonic");
    legKE_Hamp()->AddEntry((*hKE_elpsH_amp)(), "e^{-}/e^{+} flux ratio", "lp");
    legKE_Hamp()->AddEntry((*hKE_apprH_amp)(), "#bar{p}/p flux ratio", "lp");
    legKE_Hamp()->SetFillColor(0);
    legKE_Hamp.draw();
    editor.save();

    // Harmonic <Delay>
    //Axis AX_Hdelay("t_{1/2} - t_{rev} (years)", 2000, 0.0 * YEAR2DAY, 5.0 * YEAR2DAY);
    Axis AX_Hdelay("H_{S}^{-Z/+Z} (years)", 2000, 0.0 * YEAR2DAY, 5.0 * YEAR2DAY);
    Hist* hRIG_Hdelay = Hist::New("hRIG_Hdelay", HistAxis(AXrig, AX_Hdelay));
    (*hRIG_Hdelay)()->GetXaxis()->SetMoreLogLabels();
    (*hRIG_Hdelay)()->GetXaxis()->CenterTitle();
    (*hRIG_Hdelay)()->GetYaxis()->CenterTitle();
    (*hRIG_Hdelay)()->GetXaxis()->SetRangeUser(1.01, 7.5);
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIG_Hdelay)()->Draw();
    Helps_delay->Draw("l same");
    //Happr_delay->Draw("l same");
    (*hRIG_elpsH_delay)()->Draw("pe same");
    (*hRIG_apprH_delay)()->Draw("pe same");
    Legend legRIG_Hdelay("", TextStyle(kBlack, 30, 43), PadWindow(0.3, 0.6, 0.70, 0.90));
    //legRIG_Hdelay()->SetHeader("Harmonic");
    legRIG_Hdelay()->AddEntry((*hRIG_elpsH_delay)(), "e^{-}/e^{+} flux ratio", "lp");
    legRIG_Hdelay()->AddEntry((*hRIG_apprH_delay)(), "#bar{p}/p flux ratio", "lp");
    legRIG_Hdelay()->AddEntry((TObject*)0, Form("H_{S}^{e^{-}/e^{+}} #approx (%4.2f#pm%4.2f)  R^{%4.2f#pm%4.2f}", elpsH_delay[0], elpsH_delay[1], elpsH_delay[2], elpsH_delay[3]), "");
    legRIG_Hdelay()->SetFillColor(0);
    legRIG_Hdelay.draw();
    editor.save();
    
    // Harmonic <Period>
    //Axis AX_Hperiod("#Delta t (years)", 2000, 4.0 * YEAR2DAY, 17.0 * YEAR2DAY);
    Axis AX_Hperiod("H_{T}^{-Z/+Z} (years)", 2000, 4.0 * YEAR2DAY, 17.0 * YEAR2DAY);
    Hist* hRIG_Hperiod = Hist::New("hRIG_Hperiod", HistAxis(AXrig, AX_Hperiod));
    (*hRIG_Hperiod)()->GetXaxis()->SetMoreLogLabels();
    (*hRIG_Hperiod)()->GetXaxis()->CenterTitle();
    (*hRIG_Hperiod)()->GetYaxis()->CenterTitle();
    (*hRIG_Hperiod)()->GetXaxis()->SetRangeUser(1.01, 7.5);
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIG_Hperiod)()->Draw();
    //Helps_period->Draw("l same");
    //Happr_period->Draw("l same");
    (*hRIG_elpsH_period)()->Draw("pe same");
    (*hRIG_apprH_period)()->Draw("pe same");
    Legend legRIG_Hperiod("", TextStyle(kBlack, 30, 43), PadWindow(0.15, 0.50, 0.60, 0.90));
    //legRIG_Hperiod()->SetHeader("Harmonic");
    legRIG_Hperiod()->AddEntry((*hRIG_elpsH_period)(), "e^{-}/e^{+} flux ratio", "lp");
    //legRIG_Hperiod()->AddEntry((TObject*)0, Form("#Delta t_{e^{-}/e^{+}} = %.2f #pm %.2f (years)", elpsH_period[0], elpsH_period[1]), "");
    legRIG_Hperiod()->AddEntry((TObject*)0, Form("H_{T}^{e^{-}/e^{+}} = %.2f #pm %.2f (years)", elpsH_period[0], elpsH_period[1]), "");
    legRIG_Hperiod()->AddEntry((*hRIG_apprH_period)(), "#bar{p}/p flux ratio", "lp");
    //legRIG_Hperiod()->AddEntry((TObject*)0, Form("#Delta t_{#bar{p}/p} = %.2f #pm %.2f (years)", apprH_period[0], apprH_period[1]), "");
    legRIG_Hperiod()->AddEntry((TObject*)0, Form("H_{T}^{#bar{p}/p} = %.2f #pm %.2f (years)", apprH_period[0], apprH_period[1]), "");
    legRIG_Hperiod()->SetFillColor(0);
    legRIG_Hperiod.draw();
    editor.save();
    
    // Logistic and Harmonic
    editor.create();
    TextDraw("Logistic & Harmonic", TextStyle(kBlack, 0.1), TextAlign(0.5, 0.5));
    editor.save();
    
    //Axis AX_shift("#Deltat_{1/2} (years)", 2000, -1.0 * YEAR2DAY, 3.5 * YEAR2DAY);
    Axis AX_shift("#DeltaS (years)", 2000, -1.0 * YEAR2DAY, 3.5 * YEAR2DAY);
    Hist* hRIG_shift = Hist::New("hRIG_shift", HistAxis(AXrig, AX_shift));
    (*hRIG_shift)()->GetXaxis()->SetMoreLogLabels();
    (*hRIG_shift)()->GetXaxis()->CenterTitle();
    (*hRIG_shift)()->GetYaxis()->CenterTitle();
    (*hRIG_shift)()->GetXaxis()->SetRangeUser(1.01, 7.5);
        
    TGraphErrors* gRIG_Lshift = new TGraphErrors((*hRIG_Lshift)());
    TGraphErrors* gRIG_Hshift = new TGraphErrors((*hRIG_Hshift)());
    for (int ip = 0; ip < gRIG_Lshift->GetN(); ++ip) gRIG_Lshift->GetX()[ip] *= 0.995;
    for (int ip = 0; ip < gRIG_Hshift->GetN(); ++ip) gRIG_Hshift->GetX()[ip] *= 1.005;
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIG_shift)()->Draw();
    baseline->SetRange(0.1, 100);
    baseline->SetLineColor(kGreen+2);
    baseline->Draw("l same");
    gRIG_Lshift->Draw("pe same");
    gRIG_Hshift->Draw("pe same");
    Legend legRIG_shift("", TextStyle(kBlack, 30, 43), PadWindow(0.30, 0.70, 0.75, 0.90));
    legRIG_shift()->AddEntry(gRIG_Lshift, "Logistic (L_{S}^{#bar{p}/p} - L_{S}^{e^{-}/e^{+}})", "lp");
    legRIG_shift()->AddEntry(gRIG_Hshift, "Harmonic (H_{S}^{#bar{p}/p} - H_{S}^{e^{-}/e^{+}})", "lp");
    legRIG_shift()->SetFillColor(0);
    legRIG_shift.draw();
    editor.save();
    
    editor.close();
    
    TFile * ofle = new TFile(Form("out/apflux_phy%s_%s.root", subt.c_str(), fitType?"np":"pn"), "RECREATE");
    ofle->cd();

    hRIG_Lshift->write();
    hRIG_Hshift->write();

    ofle->Write();
    ofle->Close();
   
    if (1) {
        TFile * ifle_np = TFile::Open(Form("out/apflux_phy%s_np.root", subt.c_str()));
        TGraphErrors* gRIG_Lshift_np = new TGraphErrors((TH1D*)(ifle_np->Get("hRIG_Lshift")));
        TGraphErrors* gRIG_Hshift_np = new TGraphErrors((TH1D*)(ifle_np->Get("hRIG_Hshift")));
        ifle_np->Close();

        for (int ip = 0; ip < gRIG_Lshift_np->GetN(); ++ip) gRIG_Lshift_np->GetX()[ip] *= 0.9775;
        for (int ip = 0; ip < gRIG_Hshift_np->GetN(); ++ip) gRIG_Hshift_np->GetX()[ip] *= 0.9925;

        TFile * ifle_pn = TFile::Open(Form("out/apflux_phy%s_pn.root", subt.c_str()));
        TGraphErrors* gRIG_Lshift_pn = new TGraphErrors((TH1D*)(ifle_pn->Get("hRIG_Lshift")));
        TGraphErrors* gRIG_Hshift_pn = new TGraphErrors((TH1D*)(ifle_pn->Get("hRIG_Hshift")));
        ifle_pn->Close();
        
        for (int ip = 0; ip < gRIG_Lshift_pn->GetN(); ++ip) gRIG_Lshift_pn->GetX()[ip] *= 1.0075;
        for (int ip = 0; ip < gRIG_Hshift_pn->GetN(); ++ip) gRIG_Hshift_pn->GetX()[ip] *= 1.0225;
        
        PdfEditor editor_merge(Window(WindowSize::kWideSliceLR), Form("apflux_phy%s", subt.c_str()), "out");
        editor_merge.create();
        editor_merge.cd(0, PadAxis(1, 0));
        (*hRIG_shift)()->Draw();
        baseline->SetRange(0.1, 100);
        baseline->SetLineColor(kGreen+2);
        baseline->Draw("l same");
        gRIG_Lshift_np->Draw("pe same");
        gRIG_Hshift_np->Draw("pe same");
        gRIG_Lshift_pn->Draw("pe same");
        gRIG_Hshift_pn->Draw("pe same");
        Legend legRIG_shift_merge("", TextStyle(kBlack, 30, 43), PadWindow(0.50, 0.85, 0.60, 0.90));
        legRIG_shift_merge()->AddEntry(gRIG_Lshift_np, "Logistic (L_{S}^{#bar{p}/p} - L_{S}^{e^{-}/e^{+}})", "lp");
        legRIG_shift_merge()->AddEntry(gRIG_Hshift_np, "Harmonic (H_{S}^{#bar{p}/p} - H_{S}^{e^{-}/e^{+}})", "lp");
        legRIG_shift_merge()->AddEntry(gRIG_Lshift_pn, "Logistic (L_{S}^{p/#bar{p}} - L_{S}^{e^{+}/e^{-}})", "lp");
        legRIG_shift_merge()->AddEntry(gRIG_Hshift_pn, "Harmonic (H_{S}^{p/#bar{p}} - H_{S}^{e^{+}/e^{-}})", "lp");
        legRIG_shift_merge()->SetFillColor(0);
        legRIG_shift_merge.draw();
        editor_merge.save();
        
        editor_merge.close();
    }

    return 1;
}
