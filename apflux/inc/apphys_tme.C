#include <CPPLibs.h>
#include <ROOTLibs.h>

static constexpr double TWO_PI = 2.0 * 3.1415927;
static constexpr double HALF_PI = 0.5 * 3.1415927;


static constexpr double FitRegLw = 1307750400;
static constexpr double FitRegUp = 1573689600;

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subt = "3M";
    
    short fitType = 0; // 0, +Z/-Z  1, -Z/+Z  flux ratio
    
    // Scale 10^5   \\ testcode
    //const double SCALE = 1.0;
    const double SCALE = 1.0e+5;
    
    //PdfEditor editor(Window(WindowSize::kWideSliceLR), Form("apphys_tme%s", subt.c_str()), "out");
    PdfEditor editor(Window(600, 350), Form("apphys_tme%s_%s", subt.c_str(), fitType?"np":"pn"), "out");

    TFile* ftme = TFile::Open(Form("out/apflux_tme%s.root", subt.c_str()));
    Hist* hTRappr = Hist::New("hTRappr", (TH2D*)(ftme->Get("hTAppTotl")));
    hTRappr->scale(SCALE);
    
    Axis AXtme = hTRappr->xaxis(); AXtme.set_title("Date");
    Axis AXrig = hTRappr->yaxis(); AXrig.set_title("|Rigidity| [GV]");

    std::vector<double> BINke;
    for (int i = 0; i <= AXrig.nbin(); ++i) {
        double ke = 0.938272297 * (std::hypot(AXrig(i), 1.0) - 1.0);
        BINke.push_back(ke);
    }
    Axis AXke("Kinetic Energy [GeV]", BINke);

    Hist* hFLCFlux = Hist::New("hFLCFlux", HistAxis(AXrig));
    Hist* hMIDFlux = Hist::New("hMIDFlux", HistAxis(AXrig));
    Hist* hMAXFlux = Hist::New("hMAXFlux", HistAxis(AXrig));
    Hist* hMINFlux = Hist::New("hMINFlux", HistAxis(AXrig));
    for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
        for (int it = 1; it <= AXtme.nbin(); ++it) {
            double val = (*hTRappr)()->GetBinContent(it, ir);
            double err = (*hTRappr)()->GetBinError  (it, ir);
            if (val <= 0.0) continue;
            bool update_maxflux = ((*hMAXFlux)()->GetBinContent(ir) > 0) ? (val > (*hMAXFlux)()->GetBinContent(ir)) : true;
            bool update_minflux = ((*hMINFlux)()->GetBinContent(ir) > 0) ? (val < (*hMINFlux)()->GetBinContent(ir)) : true;
            if (update_maxflux) { (*hMAXFlux)()->SetBinContent(ir, val); (*hMAXFlux)()->SetBinError(ir, err); }
            if (update_minflux) { (*hMINFlux)()->SetBinContent(ir, val); (*hMINFlux)()->SetBinError(ir, err); }
        }
        double midval = 0.5 * ((*hMAXFlux)()->GetBinContent(ir) + (*hMINFlux)()->GetBinContent(ir));
        double miderr = 0.5 * ((*hMAXFlux)()->GetBinContent(ir) - (*hMINFlux)()->GetBinContent(ir));
        (*hMIDFlux)()->SetBinContent(ir, midval);
        (*hMIDFlux)()->SetBinError  (ir, miderr);
        (*hFLCFlux)()->SetBinContent(ir, miderr / midval);
        (*hFLCFlux)()->SetBinError  (ir, 1.0e-6);
    }
    hFLCFlux->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hMIDFlux->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hMAXFlux->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hMINFlux->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));

    TF1* func_apflux = new TF1("func_apflux", "[0] * pow((x/45.0), [1]) / (1.0 + pow(x/[2], [3]))", 0.5, 1000.0);
    func_apflux->SetParameters(1.97737e-04 * SCALE, -4.30173e-02, 5.34225e+00, -2.09848e+00);
    (*hMIDFlux)()->Fit(func_apflux, "q0", "");
    (*hMIDFlux)()->Fit(func_apflux, "q0", "");
    (*hMIDFlux)()->Fit(func_apflux, "q0", "");

    TF1* func_apfluc = new TF1("func_apfluc", "[0] * pow(x, [1]) + [2]", 0.5, 1000.0);
    func_apfluc->SetParameters(5.92881e-01, -1.41457e+00, 2.63987e-02);
    (*hFLCFlux)()->Fit(func_apfluc, "q0", "");
    (*hFLCFlux)()->Fit(func_apfluc, "q0", "");
    (*hFLCFlux)()->Fit(func_apfluc, "q0", "");

    Hist* hAPFlucFlux = Hist::New("hAPFlucFlux", HistAxis(Axis(AXrig.title(), 10000, AXrig.min(), AXrig.max(), AxisScale::kLog)));
    for (int ir = 1; ir <= hAPFlucFlux->xaxis().nbin(); ++ir) {
        double cen = hAPFlucFlux->xaxis().center(ir, AxisScale::kLog);
        (*hAPFlucFlux)()->SetBinContent(ir, func_apflux->Eval(cen));
        (*hAPFlucFlux)()->SetBinError  (ir, func_apflux->Eval(cen) * func_apfluc->Eval(cen));
    }
    
    TGraphErrors* gAPFlucFlux = new TGraphErrors((*hAPFlucFlux)());
    gAPFlucFlux->SetMarkerColor(kRed);
    gAPFlucFlux->SetLineColor(kRed);
    gAPFlucFlux->SetFillColor(kRed);
    gAPFlucFlux->SetFillStyle(3002);
    
    Hist* hAppTotl = Hist::New("hAppTotl", (TH1D*)(TFile::Open("out/apflux_flx.root")->Get("hAppTotl")));
    hAppTotl->style(Line(kRed, 0, 3), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hAppTotl->scale(SCALE);
    
    TGraphErrors* gAppOffical = (TGraphErrors*) (TFile::Open("others/20160406MIT.root")->Get("gpbarp"));
    gAppOffical->SetLineColor(kBlue);
    gAppOffical->SetMarkerColor(kBlue);
    gAppOffical->SetLineWidth(3);

    Hist* hcvsAppFR = Hist::New("hcvsAppFR", HistAxis(AXrig, Axis("#bar{p}/p Flux Ratio [10^{5}]", 2000, 2.0e-6 * SCALE, 3.0e-4 * SCALE)));
    (*hcvsAppFR)()->GetXaxis()->SetMoreLogLabels();
    (*hcvsAppFR)()->GetXaxis()->CenterTitle();
    (*hcvsAppFR)()->GetYaxis()->CenterTitle();
    //(*hcvsAppFR)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    //(*hcvsAppFR)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
    (*hcvsAppFR)()->GetXaxis()->SetRangeUser(1.0, 25.0);

    editor.create();
    editor.cd(0, PadAxis(1, 1));
    (*hcvsAppFR)()->Draw();
    gAPFlucFlux->Draw("3");
    gAppOffical->Draw("pe same");
    (*hAppTotl)()->Draw("pe same");
    Legend leg_apflux_comparison("", TextStyle(kBlack, 50, 43), PadWindow(0.30, 0.90, 0.20, 0.40));
    leg_apflux_comparison()->AddEntry(gAppOffical, "PRL _{(PhysRevLett.117.091103)}", "lp");
    leg_apflux_comparison()->AddEntry((*hAppTotl)(), "Result", "lp");
    leg_apflux_comparison()->SetFillColor(0);
    leg_apflux_comparison.draw();
    editor.save();

    // Fitting
    const double REV_UTIME = 1372636800.0; // 2013, July, 01
    const double SAT_UTIME = AXtme.min();
    const double END_UTIME = AXtme.max();
    const double DAY  = 86400.0;
    const double YEAR = 86400.0 * 365.25;
    //const double YEAR2DAY = 365.25;
    const double YEAR2DAY = 1.0;
   
    // linear
    TF1* tfuncM = new TF1("tfuncM", "[0]", SAT_UTIME, END_UTIME);
    
    // Logistic function
    TF1* tfuncL = nullptr;
    if (fitType==0) tfuncL = new TF1("tfuncL", Form("[0] / (1.0 + [1] / ( 1.0 + exp(-((x - %f) / %f - [2]) / ([3] / 4.39)) ))", REV_UTIME, YEAR), SAT_UTIME, END_UTIME);
    if (fitType==1) tfuncL = new TF1("tfuncL", Form("[0] * (1.0 + [1] / ( 1.0 + exp(((x - %f) / %f - [2]) / ([3] / 4.39)) ))", REV_UTIME, YEAR), SAT_UTIME, END_UTIME);
    
    tfuncL->SetNpx(10000);
    tfuncL->SetParLimits(1, -0.1, 4.0);
    tfuncL->SetParLimits(2,  0.0, 4.0);
    tfuncL->SetParLimits(3,  0.5, 4.0);
    tfuncL->SetLineColor(kBlue);
    
    const double DEF_L_AMP    = 0.5;
    const double DEF_L_DELAY  = 2.0;
    const double DEF_L_PERIOD = 2.0;
    
    Hist* hRIG_L_amp    = Hist::New("hRIG_L_amp",    HistAxis(AXrig));
    Hist* hRIG_L_delay  = Hist::New("hRIG_L_delay",  HistAxis(AXrig));
    Hist* hRIG_L_period = Hist::New("hRIG_L_period", HistAxis(AXrig));
    
    Hist* hKE_L_amp    = Hist::New("hKE_L_amp",    HistAxis(AXke));
    Hist* hKE_L_delay  = Hist::New("hKE_L_delay",  HistAxis(AXke));
    Hist* hKE_L_period = Hist::New("hKE_L_period", HistAxis(AXke));
    
    // harmonic function
    TF1* tfuncH = nullptr;
    if (fitType==0) tfuncH = new TF1("tfuncH", Form("[0] / (1.0 - [1] * cos(%14.8f + %14.8f * (((x - %f) / %f) - [2]) / [3]))", HALF_PI, TWO_PI, REV_UTIME, YEAR), SAT_UTIME, END_UTIME);
    if (fitType==1) tfuncH = new TF1("tfuncH", Form("[0] * (1.0 + [1] * cos(%14.8f + %14.8f * (((x - %f) / %f) - [2]) / [3]))", HALF_PI, TWO_PI, REV_UTIME, YEAR), SAT_UTIME, END_UTIME);

    tfuncH->SetNpx(10000);
    tfuncH->SetParLimits(1,  0.0, 0.8);
    tfuncH->SetParLimits(2,  0.0, 4.0);
    tfuncH->SetParLimits(3,  6.0, 15.0);
    tfuncH->SetLineColor(kBlue);
    
    const double DEF_H_AMP    =  0.2;
    const double DEF_H_DELAY  =  0.2;
    const double DEF_H_PERIOD =  9.0;

    Hist* hRIG_H_amp    = Hist::New("hRIG_H_amp",    HistAxis(AXrig));
    Hist* hRIG_H_delay  = Hist::New("hRIG_H_delay",  HistAxis(AXrig));
    Hist* hRIG_H_period = Hist::New("hRIG_H_period", HistAxis(AXrig));
    
    Hist* hKE_H_amp    = Hist::New("hKE_H_amp",    HistAxis(AXke));
    Hist* hKE_H_delay  = Hist::New("hKE_H_delay",  HistAxis(AXke));
    Hist* hKE_H_period = Hist::New("hKE_H_period", HistAxis(AXke));


    Hist* hPOLTcanvas = Hist::New("hPOLTcanvas", HistAxis(Axis("", 100, 1305417600, 1573689600), Axis("", 1000, 0.0, 5.0e-4 * SCALE)));
    (*hPOLTcanvas)()->GetXaxis()->SetTitle("");
    (*hPOLTcanvas)()->GetXaxis()->SetTimeDisplay(true);
    (*hPOLTcanvas)()->GetXaxis()->SetTimeOffset(0, "GMT");
    (*hPOLTcanvas)()->GetXaxis()->SetTimeFormat("%y/%b");
    (*hPOLTcanvas)()->GetXaxis()->SetLabelSize((*hPOLTcanvas)()->GetXaxis()->GetLabelSize() * 1.8);
    (*hPOLTcanvas)()->GetYaxis()->SetLabelSize((*hPOLTcanvas)()->GetYaxis()->GetLabelSize() * 2.2);
    
    std::vector<Hist*> vhRfr = Hist::ProjectAll(HistProj::kX, hTRappr);
    for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        Hist* hfr  = vhRfr.at(ir);

        tfuncM->SetParameter(0, (*hfr)()->GetBinContent(1));
        (*hfr)()->Fit(tfuncM, "q0", "", FitRegLw, FitRegUp);
        (*hfr)()->Fit(tfuncM, "q0", "", FitRegLw, FitRegUp);
        
        // Logistic function
        tfuncL->SetParameters(tfuncM->GetParameter(0), DEF_L_AMP, DEF_L_DELAY, DEF_L_PERIOD);
        (*hfr)()->Fit(tfuncL, "q0", "", FitRegLw, FitRegUp);
        (*hfr)()->Fit(tfuncL, "q0", "", FitRegLw, FitRegUp);

        (*hRIG_L_amp)()->SetBinContent(ir, tfuncL->GetParameter(1));
        (*hRIG_L_amp)()->SetBinError  (ir, tfuncL->GetParError (1));
        
        (*hRIG_L_delay)()->SetBinContent(ir, YEAR2DAY * tfuncL->GetParameter(2));
        (*hRIG_L_delay)()->SetBinError  (ir, YEAR2DAY * tfuncL->GetParError (2));
        
        (*hRIG_L_period)()->SetBinContent(ir, YEAR2DAY * tfuncL->GetParameter(3));
        (*hRIG_L_period)()->SetBinError  (ir, YEAR2DAY * tfuncL->GetParError (3));
        
        (*hKE_L_amp)()->SetBinContent(ir, tfuncL->GetParameter(1));
        (*hKE_L_amp)()->SetBinError  (ir, tfuncL->GetParError (1));
        
        (*hKE_L_delay)()->SetBinContent(ir, YEAR2DAY * tfuncL->GetParameter(2));
        (*hKE_L_delay)()->SetBinError  (ir, YEAR2DAY * tfuncL->GetParError (2));
        
        (*hKE_L_period)()->SetBinContent(ir, YEAR2DAY * tfuncL->GetParameter(3));
        (*hKE_L_period)()->SetBinError  (ir, YEAR2DAY * tfuncL->GetParError (3));

        // Harmonic function
        tfuncH->SetParameters(tfuncM->GetParameter(0), DEF_H_AMP, DEF_H_DELAY, DEF_H_PERIOD);
        (*hfr)()->Fit(tfuncH, "q0", "", FitRegLw, FitRegUp);
        (*hfr)()->Fit(tfuncH, "q0", "", FitRegLw, FitRegUp);

        (*hRIG_H_amp)()->SetBinContent(ir, tfuncH->GetParameter(1));
        (*hRIG_H_amp)()->SetBinError  (ir, tfuncH->GetParError (1));
        
        (*hRIG_H_delay)()->SetBinContent(ir, YEAR2DAY * tfuncH->GetParameter(2));
        (*hRIG_H_delay)()->SetBinError  (ir, YEAR2DAY * tfuncH->GetParError (2));
       
        (*hRIG_H_period)()->SetBinContent(ir, YEAR2DAY * tfuncH->GetParameter(3));
        (*hRIG_H_period)()->SetBinError  (ir, YEAR2DAY * tfuncH->GetParError (3));
       
        (*hKE_H_amp)()->SetBinContent(ir, tfuncH->GetParameter(1));
        (*hKE_H_amp)()->SetBinError  (ir, tfuncH->GetParError (1));
        
        (*hKE_H_delay)()->SetBinContent(ir, YEAR2DAY * tfuncH->GetParameter(2));
        (*hKE_H_delay)()->SetBinError  (ir, YEAR2DAY * tfuncH->GetParError (2));
        
        (*hKE_H_period)()->SetBinContent(ir, YEAR2DAY * tfuncH->GetParameter(3));
        (*hKE_H_period)()->SetBinError  (ir, YEAR2DAY * tfuncH->GetParError (3));

        hfr->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
        (*hfr)()->GetYaxis()->SetTitle("#bar{p}/p flux ratio");
        (*hfr)()->GetXaxis()->SetTitle("");
        (*hfr)()->GetXaxis()->SetTimeDisplay(true);
        (*hfr)()->GetXaxis()->SetTimeOffset(0, "GMT");
        (*hfr)()->GetXaxis()->SetTimeFormat("%y/%b");
        
        double POLTcen = tfuncL->GetParameter(0) * 0.5 * (1.0 + 1.0 / (1.0 + tfuncL->GetParameter(1)));
        //double POLTlw  = POLTcen * (1.0 - 0.5);
        //double POLTup  = POLTcen * (1.0 + 0.5);
        //if (ir == 1 || ir == 2) POLTlw  = POLTcen * (1.0 - 1.0);
        //if (ir == 1 || ir == 2) POLTup  = POLTcen * (1.0 + 1.5);
        double POLTamp = 0.3 * (1.0 + 2.0 / std::sqrt(rig) + 3.0 / rig);
        //double POLTamp = 0.3 * (1.0 + 5.0 / rig);
        double POLTlw = POLTcen * (1.0 - POLTamp);
        double POLTup = POLTcen * (1.0 + POLTamp);
        if (POLTlw < 0.0) POLTlw = 0.0;

        editor.create("", PadMargin(0.1, 0.1, 0.08, 0.02));
        editor.cd(0, PadAxis(0, 0));
        (*hPOLTcanvas)()->GetYaxis()->SetRangeUser(POLTlw, POLTup);
        (*hPOLTcanvas)()->Draw();
        (*hfr)()->Draw("pe same");
        //tfuncL->Draw("l same");
        tfuncH->Draw("l same");
        TextDraw(Form("%.2f - %.2f [GV]", AXrig(ir-1), AXrig(ir)), TextStyle(kBlack, 0.1), TextAlign(0.7, 0.8));
        editor.save();
    }

    hRIG_L_amp   ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRIG_L_delay ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRIG_L_period->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    hKE_L_amp   ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_L_delay ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_L_period->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    hRIG_H_amp   ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRIG_H_delay ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hRIG_H_period->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    hKE_H_amp   ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_H_delay ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hKE_H_period->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    Axis AXamp("C", 2000, -0.2, 1.5);
    Hist* hRIGamp = Hist::New("hRIGamp", HistAxis(AXrig, AXamp));
    (*hRIGamp)()->GetXaxis()->SetMoreLogLabels();
    (*hRIGamp)()->GetXaxis()->CenterTitle();
    (*hRIGamp)()->GetYaxis()->CenterTitle();
    (*hRIGamp)()->GetXaxis()->SetRangeUser(1.01, 45.0);
    Hist* hKEamp = Hist::New("hKEamp", HistAxis(AXke, AXamp));
    (*hKEamp)()->GetXaxis()->SetMoreLogLabels();
    (*hKEamp)()->GetXaxis()->CenterTitle();
    (*hKEamp)()->GetYaxis()->CenterTitle();
    (*hKEamp)()->GetXaxis()->SetRangeUser(1.01, 45.0);
    
    Axis AXdelay("t_{1/2} - t_{rev} (years)", 2000, -4.0 * YEAR2DAY, 4.0 * YEAR2DAY);
    Hist* hRIGdelay = Hist::New("hRIGdelay", HistAxis(AXrig, AXdelay));
    (*hRIGdelay)()->GetXaxis()->SetMoreLogLabels();
    (*hRIGdelay)()->GetXaxis()->CenterTitle();
    (*hRIGdelay)()->GetYaxis()->CenterTitle();
    (*hRIGdelay)()->GetXaxis()->SetRangeUser(1.01, 7.0);
    Hist* hKEdelay = Hist::New("hKEdelay", HistAxis(AXke, AXdelay));
    (*hKEdelay)()->GetXaxis()->SetMoreLogLabels();
    (*hKEdelay)()->GetXaxis()->CenterTitle();
    (*hKEdelay)()->GetYaxis()->CenterTitle();
    (*hKEdelay)()->GetXaxis()->SetRangeUser(1.01, 7.0);
    
    Axis AXperiod("#Delta t (years)", 2000, 0.0 * YEAR2DAY, 15.0 * YEAR2DAY);
    Hist* hRIGperiod = Hist::New("hRIGperiod", HistAxis(AXrig, AXperiod));
    (*hRIGperiod)()->GetXaxis()->SetMoreLogLabels();
    (*hRIGperiod)()->GetXaxis()->CenterTitle();
    (*hRIGperiod)()->GetYaxis()->CenterTitle();
    (*hRIGperiod)()->GetXaxis()->SetRangeUser(1.01, 7.0);
    Hist* hKEperiod = Hist::New("hKEperiod", HistAxis(AXke, AXperiod));
    (*hKEperiod)()->GetXaxis()->SetMoreLogLabels();
    (*hKEperiod)()->GetXaxis()->CenterTitle();
    (*hKEperiod)()->GetYaxis()->CenterTitle();
    (*hKEperiod)()->GetXaxis()->SetRangeUser(1.01, 7.0);
        
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIGamp)()->Draw();
    (*hRIG_L_amp)()->Draw("pe same");
    editor.save();
    
    editor.create();
    (*hRIGdelay)()->Draw();
    (*hRIG_L_delay)()->Draw("pe same");
    editor.save();
    
    editor.create();
    (*hRIGperiod)()->Draw();
    (*hRIG_L_period)()->Draw("pe same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hRIGamp)()->Draw();
    (*hRIG_H_amp)()->Draw("pe same");
    editor.save();
    
    editor.create();
    (*hRIGdelay)()->Draw();
    (*hRIG_H_delay)()->Draw("pe same");
    editor.save();
    
    editor.create();
    (*hRIGperiod)()->Draw();
    (*hRIG_H_period)()->Draw("pe same");
    editor.save();
    
    TFile * ofle = new TFile(Form("out/apphys_tme%s_%s.root", subt.c_str(), fitType?"np":"pn"), "RECREATE");
    ofle->cd();
    
    (*hFLCFlux)()->Write();
    (*hMIDFlux)()->Write();
    (*hMAXFlux)()->Write();
    (*hMINFlux)()->Write();
    (*hAPFlucFlux)()->Write();

    (*hRIG_L_amp   )()->Write();
    (*hRIG_L_delay )()->Write();
    (*hRIG_L_period)()->Write();
    
    (*hKE_L_amp   )()->Write();
    (*hKE_L_delay )()->Write();
    (*hKE_L_period)()->Write();
    
    (*hRIG_H_amp   )()->Write();
    (*hRIG_H_delay )()->Write();
    (*hRIG_H_period)()->Write();
    
    (*hKE_H_amp   )()->Write();
    (*hKE_H_delay )()->Write();
    (*hKE_H_period)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
