#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "HistFit1D.h"
#include "HistFit1D.C"

void stdfmt(TH1* hist) {
    if (hist == nullptr) return;
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleFont(43);
    hist->GetXaxis()->SetTitleSize(20);
    hist->GetXaxis()->SetTitleOffset(2.5);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelSize(15);
    hist->GetYaxis()->SetTitleFont(43);
    hist->GetYaxis()->SetTitleSize(20);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelSize(15);
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subt = "3M";
    
    TFile* file_flx_CF = TFile::Open("out/apflux_flx_CF.root");
    TH1D* happ_stat_CF = (TH1D*) file_flx_CF->Get("hAppStat");
    
    TFile* file_flx  = TFile::Open("out/apflux_flx.root");
    TH1D*  hcnt_pr   = (TH1D*) file_flx->Get("hNcntPr");
    TH1D*  hcnt_ap   = (TH1D*) file_flx->Get("hNcntAp");
    TH1D*  happ_stat = (TH1D*) file_flx->Get("hAppStat");
    TH1D*  happ_totl = (TH1D*) file_flx->Get("hAppTotl");
    TH1D*  herr_stat = (TH1D*) file_flx->Get("hErrStat");
    TH1D*  herr_syst = (TH1D*) file_flx->Get("hErrSyst");
    TH1D*  herr_ncnt = (TH1D*) file_flx->Get("hErrNcnt");
    TH1D*  herr_accp = (TH1D*) file_flx->Get("hErrAccp");
    TH1D*  herr_rfnc = (TH1D*) file_flx->Get("hErrRfnc");
    
    TFile* file_Tflx  = TFile::Open(Form("out/apflux_tme%s.root", subt.c_str()));
    TH1D*  hTcnt_pr   = (TH1D*) file_Tflx->Get("hTNcntPr");
    TH1D*  hTcnt_ap   = (TH1D*) file_Tflx->Get("hTNcntAp");
    TH1D*  hTapp_stat = (TH1D*) file_Tflx->Get("hTAppStat");
    TH1D*  hTapp_totl = (TH1D*) file_Tflx->Get("hTAppTotl");
    TH1D*  hTerr_stat = (TH1D*) file_Tflx->Get("hTErrStat");
    TH1D*  hTerr_syst = (TH1D*) file_Tflx->Get("hTErrSyst");
    TH1D*  hTerr_ncnt = (TH1D*) file_Tflx->Get("hTErrNcnt");
    TH1D*  hTerr_accp = (TH1D*) file_Tflx->Get("hTErrAccp");
    TH1D*  hTerr_rfnc = (TH1D*) file_Tflx->Get("hTErrRfnc");
    
    const Axis AXrig("|Rigidity| [GV]", happ_totl->GetXaxis());
    
    const Axis AXTtme("Date", hTapp_totl->GetXaxis());
    const Axis AXTrig("|Rigidity| [GV]", hTapp_totl->GetYaxis());
    
    Hist* hErrCfftRaw = Hist::New("hErrCfftRaw", HistAxis(AXrig, "Relative Error (%)"));
    Hist* hErrCfftBD  = Hist::New("hErrCfftBD" , HistAxis(Axis("", 1000, 0.1, 1000, AxisScale::kLog), "Relative Error (%)"));
    Hist* hErrCfft    = Hist::New("hErrCfft"   , HistAxis(AXrig, "Relative Error (%)"));
		
    hErrCfftRaw->style(Line(kBlack, 0, 2), Marker(kBlack, MarkerStyle(MarkerShape::kCircle )));
    hErrCfftBD->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
    hErrCfft->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
    
    PdfEditor editor(Window(WindowSize::kWideSliceLR), "apflux_rlt", "out");

    // Geomagnetic Rigidity Cutoff
    for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
        double rat = happ_stat->GetBinContent(ir) / happ_stat_CF->GetBinContent(ir);
        double err = happ_stat_CF->GetBinError(ir) / happ_stat_CF->GetBinContent(ir);
        (*hErrCfftRaw)()->SetBinContent(ir, rat);
        (*hErrCfftRaw)()->SetBinError  (ir, err);
    }

    TF1* errfunc_cfft = new TF1("errfunc_cfft", "[0]*pow(x,[1])*exp([2]*x) + 1.0", 0.1, 1000);
    errfunc_cfft->SetParameters(0.05, 0.1, -0.1);
    errfunc_cfft->SetLineColor(kRed);
    errfunc_cfft->SetNpx(100000);

    (*hErrCfftRaw)()->Fit(errfunc_cfft, "q0", "");
    (*hErrCfftRaw)()->Fit(errfunc_cfft, "q0", "");
    std::cerr << Form("VAL %14.8f %14.8f %14.8f\n", errfunc_cfft->GetParameter(0), errfunc_cfft->GetParameter(1), errfunc_cfft->GetParameter(2));
    std::cerr << Form("ERR %14.8f %14.8f %14.8f\n", errfunc_cfft->GetParError(0), errfunc_cfft->GetParError(1), errfunc_cfft->GetParError(2));

    errfunc_cfft->SetParameter(0, errfunc_cfft->GetParameter(0) + 1.5 * errfunc_cfft->GetParError(0));
    for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
        (*hErrCfft  )()->SetBinContent(ir, errfunc_cfft->Eval(AXrig.center(ir, AxisScale::kLog))-1.0);
    }
    
    for (int ib = 1; ib <= hErrCfftBD->xaxis().nbin(); ++ib) {
        (*hErrCfftBD)()->SetBinContent(ib, 1.0);
        (*hErrCfftBD)()->SetBinError  (ib, errfunc_cfft->Eval(hErrCfftBD->xaxis().center(ib, AxisScale::kLog))-1.0);
    }
    TGraphErrors* gcfft = new TGraphErrors((*hErrCfftBD)());
    gcfft->SetMarkerColor(kRed);
    gcfft->SetLineColor(kRed);
    gcfft->SetFillColor(kRed);
    gcfft->SetFillStyle(3002);
    
    TF1* errfunc_cfft_U = new TF1("errfunc_cfft_U", "[0]*pow(x,[1])*exp([2]*x) + 1.0", 0.1, 1000);
    TF1* errfunc_cfft_L = new TF1("errfunc_cfft_L", "[0]*pow(x,[1])*exp([2]*x) + 1.0", 0.1, 1000);
    errfunc_cfft_U->SetParameters(+errfunc_cfft->GetParameter(0), errfunc_cfft->GetParameter(1), errfunc_cfft->GetParameter(2));
    errfunc_cfft_L->SetParameters(-errfunc_cfft->GetParameter(0), errfunc_cfft->GetParameter(1), errfunc_cfft->GetParameter(2));
    errfunc_cfft_U->SetLineColor(kRed);
    errfunc_cfft_L->SetLineColor(kRed);
    errfunc_cfft_U->SetNpx(100000);
    errfunc_cfft_L->SetNpx(100000);
    
    Hist* hcvsErrCfft = Hist::New("hcvsErrCfft", HistAxis(AXrig, Axis("", 100, 0.80, 1.20)));
    (*hcvsErrCfft)()->GetXaxis()->SetMoreLogLabels();
    (*hcvsErrCfft)()->GetXaxis()->CenterTitle();
    (*hcvsErrCfft)()->GetYaxis()->CenterTitle();
    (*hcvsErrCfft)()->GetXaxis()->SetRangeUser(1.0, 50.0);	
    (*hcvsErrCfft)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hcvsErrCfft)()->GetYaxis()->SetTitle("(#bar{p}/p)_{100%} / (#bar{p}/p)_{140%}");
    editor.create();
    editor.cd(1, PadAxis(1, 0));
    (*hcvsErrCfft)()->Draw();
    gcfft->Draw("3");
    errfunc_cfft_U->Draw("l same");
    errfunc_cfft_L->Draw("l same");
    (*hErrCfftRaw)()->Draw("pe same");
    Legend leg_cfft("", TextStyle(kBlack, 40, 43), PadWindow(0.55, 0.85, 0.65, 0.85));
    leg_cfft()->AddEntry((*hErrCfftRaw)(), "Data", "lp");
    leg_cfft()->AddEntry((*hErrCfftBD)(), "Systematic Error", "l");
    leg_cfft()->SetFillColor(0);
    leg_cfft.draw();
    editor.save();

    std::cerr << std::endl;
    std::cerr << Form("GeoCutoff PARAMS %14.8f %14.8f %14.8f\n", errfunc_cfft->GetParameter(0), errfunc_cfft->GetParameter(1), errfunc_cfft->GetParameter(2));
    std::cerr << std::endl;
   
    // Power Law Fit
    TF1* locpow = new TF1("locpow", "[0]*pow(x/[2], [1])", 0.1, 1000);
    locpow->SetParameters(1.0e-4, 0.01, 1.0);
    locpow->SetLineColor(kRed);
    
    std::vector<double> binFitR;
    for (int ir = 3; ir <= AXrig.nbin(); ir+=5) {
        binFitR.push_back(AXrig(ir));
    }
    Axis AXfitR("|Rigidity| [GV]", binFitR);

    TGraphErrors* gFlxPowC = new TGraphErrors();
    gFlxPowC->SetNameTitle("gFlxPowC", "gFlxPowC");

    TGraphErrors* gFlxPowK = new TGraphErrors();
    gFlxPowK->SetNameTitle("gFlxPowK", "gFlxPowK");

    for (int ir = 1; ir <= AXfitR.nbin(); ++ir) {
        double wgtr = 0.0;
        double wgte = 0.0;
        for (int jr = (ir-1)*5 + 4; jr <= ir*5 + 3; ++jr) {
           double cenr = AXrig.center(jr, AxisScale::kLog);
           double val  = happ_totl->GetBinContent(jr);
           double err  = happ_totl->GetBinError(jr) / val;
           wgtr += cenr * (1.0/err/err);
           wgte += (1.0/err/err);
        }
        double rig = (wgtr / wgte);

        double bdlw = AXfitR(ir-1);
        double bdup = AXfitR(ir);
        locpow->SetParameters(happ_totl->GetBinContent(ir * 5 + 3), 0.0, rig);
        locpow->FixParameter(2, rig);
        happ_totl->Fit(locpow, "q0", "", bdlw, bdup);
        happ_totl->Fit(locpow, "q0", "", bdlw, bdup);
        happ_totl->Fit(locpow, "q0", "", bdlw, bdup);
        gFlxPowC->SetPoint     (ir, rig, 1.0e+4 * locpow->GetParameter(0));
        gFlxPowC->SetPointError(ir, 0.0, 1.0e+4 * locpow->GetParError(0));
        gFlxPowK->SetPoint     (ir, rig, locpow->GetParameter(1));
        gFlxPowK->SetPointError(ir, 0.0, locpow->GetParError(1));
    }
    gFlxPowC->SetMarkerColor(kRed);
    gFlxPowC->SetLineColor(kRed);
    gFlxPowK->SetMarkerColor(kRed);
    gFlxPowK->SetLineColor(kRed);
    
    Hist* hcvsFlxPow = Hist::New("hcvsFlxPow", HistAxis(AXfitR, Axis("", 1000, -0.5, 2.5)));
    (*hcvsFlxPow)()->GetXaxis()->SetMoreLogLabels();
    (*hcvsFlxPow)()->GetXaxis()->CenterTitle();
    (*hcvsFlxPow)()->GetYaxis()->CenterTitle();
    (*hcvsFlxPow)()->GetXaxis()->SetTitle("|Rigidity| [GV]");

    TF1* baseline = new TF1("baseline", "[0]", 0.0, 1000.0);
    baseline->SetParameter(0, 0.0);
    baseline->SetLineWidth(1);
    baseline->SetLineStyle(3);

    editor.create();
    editor.cd(1, PadAxis(1, 0));
    (*hcvsFlxPow)()->GetYaxis()->SetRangeUser(0.0, 2.5);	
    (*hcvsFlxPow)()->GetYaxis()->SetTitle("C [10^{-4}]");
    (*hcvsFlxPow)()->Draw();
    baseline->SetParameter(0, 2.0);
    baseline->Draw("l same");
    gFlxPowC->Draw("pe");
    editor.save();
    
    editor.create();
    editor.cd(1, PadAxis(1, 0));
    (*hcvsFlxPow)()->GetYaxis()->SetRangeUser(-0.5, 2.5);
    (*hcvsFlxPow)()->GetYaxis()->SetTitle("k");
    (*hcvsFlxPow)()->Draw();
    baseline->SetParameter(0, 0.0);
    baseline->Draw("l same");
    gFlxPowK->Draw("pe");
    editor.save();
    
    // Table of Results (pbar/p)
    std::cerr << Form("==== The Antiproton-to-Proton Flux Ratio ====\n");
    for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
        double rigl = AXrig(ir-1);
        double rigu = AXrig(ir);
        double cntpr    = std::rint(hcnt_pr->GetBinContent(ir)) * 1.0e-6;
        double cntap    = std::rint(hcnt_ap->GetBinContent(ir));
        double app_vals = happ_totl->GetBinContent(ir) * 1.0e+4;

        int optpow = (app_vals < 1.0) ? (app_vals < 0.1 ? 2 : 1) : 0;
        if (optpow == 1) app_vals *= 1.0e+1;
        if (optpow == 2) app_vals *= 1.0e+2;
        
        int powidx = -4;
        if (optpow == 1) powidx = -5;
        if (optpow == 2) powidx = -6;

        double err_stat = herr_stat->GetBinContent(ir) * (app_vals / 100.0);
        double err_syst = herr_syst->GetBinContent(ir) * (app_vals / 100.0);
        double err_ncnt = herr_ncnt->GetBinContent(ir) * (app_vals / 100.0);
        double err_accp = herr_accp->GetBinContent(ir) * (app_vals / 100.0);
        double err_rfnc = herr_rfnc->GetBinContent(ir) * (app_vals / 100.0);

        //std::cout << Form("RIG %6.2f - %6.2f  CNT %6.2fM %5.0f  APP %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", rigl, rigu, cntpr, cntap, app_vals, err_stat, err_syst, err_ncnt, err_accp, err_rfnc);
        std::cout << Form("%6.2f - %6.2f & %6.2fM & %5.0f & (%5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f) & $\\times10^{%d}$\\\\\n", rigl, rigu, cntpr, cntap, app_vals, err_stat, err_syst, err_ncnt, err_accp, err_rfnc, powidx);
    }
    

    std::vector<int> exclude_times = std::vector<int>({ 16 });
    const int BRN_INIT = 2427; // BartelsRotationNumber
    std::cerr << Form("\n\n==== The Time Variation of Antiproton-to-Proton Flux Ratio ====\n");
    for (int it = 1; it <= AXTtme.nbin(); ++it) {
        int BRNsat = BRN_INIT + (it - 1) * 3;
        int BRNend = BRN_INIT + (it + 0) * 3 - 1;
        //std::cerr << Form("\n\n Time %d  Bartels Rotation Number %d %d\n", it, BRNsat, BRNend);
        
        bool skip = false;
        for (auto&& itime : exclude_times) { if (it == itime) { skip = true; break; } }
        if (skip) continue;
        
        for (int ir = 1; ir <= AXTrig.nbin(); ++ir) {
            double rigl = AXTrig(ir-1);
            double rigu = AXTrig(ir);
            double cntpr    = std::rint(hTcnt_pr->GetBinContent(it, ir)) * 1.0e-6;
            double cntap    = std::rint(hTcnt_ap->GetBinContent(it, ir));
            double app_vals = hTapp_totl->GetBinContent(it, ir) * 1.0e+4;

            int optpow = (app_vals < 1.0) ? (app_vals < 0.1 ? 2 : 1) : 0;
            if (optpow == 1) app_vals *= 1.0e+1;
            if (optpow == 2) app_vals *= 1.0e+2;
            
            int powidx = -4;
            if (optpow == 1) powidx = -5;
            if (optpow == 2) powidx = -6;

            double err_stat = hTerr_stat->GetBinContent(it, ir) * (app_vals / 100.0);
            double err_syst = hTerr_syst->GetBinContent(it, ir) * (app_vals / 100.0);
            double err_ncnt = hTerr_ncnt->GetBinContent(it, ir) * (app_vals / 100.0);
            double err_accp = hTerr_accp->GetBinContent(it, ir) * (app_vals / 100.0);
            double err_rfnc = hTerr_rfnc->GetBinContent(it, ir) * (app_vals / 100.0);

            //std::cout << Form("RIG %6.2f - %6.2f  CNT %6.2fM %5.0f  APP %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", rigl, rigu, cntpr, cntap, app_vals, err_stat, err_syst, err_ncnt, err_accp, err_rfnc);
            //std::cout << Form("%4d - %4d & %5.2f - %5.2f & %6.2fM & %5.0f & (%5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f) & $\\times10^{%d}$\\\\\n", BRNsat, BRNend, rigl, rigu, cntpr, cntap, app_vals, err_stat, err_syst, err_ncnt, err_accp, err_rfnc, powidx);
            std::cout << Form("%4d - %4d & %5.2f - %5.2f & %6.2fM & %5.0f & (%5.3f & %5.3f & %5.3f) & $\\times10^{%d}$\\\\\n", BRNsat, BRNend, rigl, rigu, cntpr, cntap, app_vals, err_stat, err_syst, powidx);
        }
        std::cerr << Form("\\hline\n");
    }

    // Function Fit on pbar/p flux ratio
    TF1* appfunc = new TF1("appfunc", "[0] * pow(x/45.0, [1]) / (1.0 + pow(x/[2], [1]-[3]))", 1, 500);
    appfunc->SetParameters(1.89044e-02, 2.10260e+00, 5.45452e+00, -6.30314e-02);
    appfunc->SetLineColor(kBlue);
    appfunc->SetNpx(10000);

    happ_totl->Fit(appfunc, "q0", "", 1.0, 500.0);
    happ_totl->Fit(appfunc, "q0", "", 1.0, 500.0);
    happ_totl->Fit(appfunc, "q0", "", 1.0, 500.0);
    
    TF1* appfuncfluc = new TF1("appfuncU", "[0] * pow(x/45.0, [1]) / (1.0 + pow(x/[2], [1]-[3]))", 1, 500);
    appfuncfluc->SetParameters(appfunc->GetParameter(0), appfunc->GetParameter(1), appfunc->GetParameter(2), appfunc->GetParameter(3));
    appfuncfluc->SetLineColor(kGreen+2);
    appfuncfluc->SetNpx(10000);

    TGraphErrors* gr_appfunc = new TGraphErrors();
    gr_appfunc->SetNameTitle("appfunc", "appfunc");
    gr_appfunc->SetLineColor(kYellow+1);
    gr_appfunc->SetMarkerColor(kYellow+1);
   
    TRandom3 rndm(0);
    Axis AXfunc("|Rigidity| [GV]", 200, 1.0, 1000.0, AxisScale::kLog);
    for (int ip = 0; ip < AXfunc.nbin(); ++ip) {
        double rig = AXfunc.center(ip+1, AxisScale::kLog);
        std::array<double, 4000> test;
        for (int itest = 0; itest < 4000; ++itest) {
            double par0 = appfunc->GetParameter(0) + 0.0 * rndm.Gaus() * appfunc->GetParError(0);
            double par1 = appfunc->GetParameter(1) + 0.0 * rndm.Gaus() * appfunc->GetParError(1);
            double par2 = appfunc->GetParameter(2) + 0.0 * rndm.Gaus() * appfunc->GetParError(2);
            double par3 = appfunc->GetParameter(3) + rndm.Gaus() * appfunc->GetParError(3);
            appfuncfluc->SetParameters(par0, par1, par2, par3);
            test.at(itest) = (appfuncfluc->Eval(rig) - appfunc->Eval(rig));
            test.at(itest) = test.at(itest) * test.at(itest);
        }
        double error = std::sqrt(std::accumulate(test.begin(), test.end(), 0.0) / test.size());
    
        gr_appfunc->SetPoint(ip, rig, appfunc->Eval(rig));
        gr_appfunc->SetPointError(ip, 0.0, error);
    }

    editor.create();
    editor.cd(1, PadAxis(1, 0));
    happ_totl->Draw("pe");
    appfunc->Draw("l same");
    gr_appfunc->Draw("le same");
    editor.save();

    editor.close();

    TFile * ofle = new TFile("out/apflux_rlt.root", "RECREATE");
    ofle->cd();

    (*hErrCfftRaw)()->Write();
    (*hErrCfftBD)()->Write();
    (*hErrCfft)()->Write();

    gFlxPowC->Write();
    gFlxPowK->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
