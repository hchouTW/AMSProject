#include <CPPLibs.h>
#include <ROOTLibs.h>

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
    
    Hist* hTCrrR = Hist::New("hTAppStat", (TH1*)TFile::Open("out/apflux_tme.root")->Get("hTAppStat"));
    std::vector<Hist*> vhTCrrR = Hist::ProjectAll(HistProj::kX, Hist::Head("hTAppStat"));
    
    for (int it = 1; it < vhTCrrR.size()-1; ++it) {
        //(*vhTCrrR.at(it))()->GetXaxis()->SetTitle("Date");
        (*vhTCrrR.at(it))()->GetXaxis()->SetTitle("");
        (*vhTCrrR.at(it))()->GetXaxis()->SetTimeDisplay(true);
        (*vhTCrrR.at(it))()->GetXaxis()->SetTimeOffset(0, "GMT");
        (*vhTCrrR.at(it))()->GetXaxis()->SetTimeFormat("%Y/%m");
        (*vhTCrrR.at(it))()->GetXaxis()->SetLabelFont(43);
        //(*vhTCrrR.at(it))()->GetXaxis()->SetLabelSize(35);
        (*vhTCrrR.at(it))()->GetXaxis()->SetLabelSize(70);
        (*vhTCrrR.at(it))()->GetXaxis()->SetTitleFont(43);
        (*vhTCrrR.at(it))()->GetXaxis()->SetTitleSize(100);
        (*vhTCrrR.at(it))()->GetXaxis()->SetTitleOffset(1.5);
        (*vhTCrrR.at(it))()->GetXaxis()->CenterTitle();
        //(*vhTCrrR.at(it))()->GetYaxis()->SetTitle("#bar{p}/p");
        (*vhTCrrR.at(it))()->GetYaxis()->SetTitle("");
        (*vhTCrrR.at(it))()->GetYaxis()->SetLabelFont(43);
        //(*vhTCrrR.at(it))()->GetYaxis()->SetLabelSize(35);
        (*vhTCrrR.at(it))()->GetYaxis()->SetLabelSize(70);
        (*vhTCrrR.at(it))()->GetYaxis()->SetTitleFont(43);
        (*vhTCrrR.at(it))()->GetYaxis()->SetTitleSize(100);
        (*vhTCrrR.at(it))()->GetYaxis()->SetTitleOffset(1.5);
        (*vhTCrrR.at(it))()->GetYaxis()->SetNoExponent();
        (*vhTCrrR.at(it))()->GetYaxis()->CenterTitle();
        vhTCrrR.at(it)->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 1.5));
        (*vhTCrrR.at(it))()->Scale(1.0e+5);
    }
    
    const Axis& AXTtme = Hist::Head("hTAppStat")->xaxis();
    const Axis& AXTrig = Hist::Head("hTAppStat")->yaxis();

    Hist* hTapp = Hist::New("hTapp", HistAxis(AXTtme, AXTrig));
    for (int it = 1; it <= AXTtme.nbin(); ++it) {
    for (int jt = 1; jt <= AXTrig.nbin(); ++jt) {
        double val = (*hTCrrR)()->GetBinContent(it, jt);
        double err = (*hTCrrR)()->GetBinError(it, jt);
        if (!std::isfinite(val)) continue;
        if (!std::isfinite(err)) continue;

        (*hTapp)()->SetBinContent(it, jt, val);
        (*hTapp)()->SetBinError  (it, jt, err);
    }}
    std::vector<Hist*> vhTapp = Hist::ProjectAll(HistProj::kX, Hist::Head("hTapp"));
   
    Hist* hTapp_dt = Hist::New("hTapp_dt", HistAxis(AXTrig));
    Hist* hTapp_st = Hist::New("hTapp_st", HistAxis(AXTrig));
    Hist* hTapp_cc = Hist::New("hTapp_cc", HistAxis(AXTrig));
    Hist* hTapp_cc_bl = Hist::New("hTapp_cc_baseline", HistAxis(AXTrig));
    (*hTapp_dt)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hTapp_dt)()->GetYaxis()->SetTitle("#Delta t (years)");
    (*hTapp_st)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hTapp_st)()->GetYaxis()->SetTitle("t_{1/2}-t_{rev} (years)");
    (*hTapp_cc)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hTapp_cc)()->GetYaxis()->SetTitle("C");

    PdfEditor editor2(Window(WindowSize::kMovieMR), "apflux_time2", "out", PdfEditor::Type::kPDF, PdfEditor::Mode::kMerge);

    //TF1* tfunc = new TF1("tfunc", Form("[0]*(1+[1]/( 1 + exp(-((x-%f)/86400. - [2]) / ([3] / 4.39)) ))", AXTtme.min()), AXTtme.min(), AXTtme.max());

    std::vector<TF1*> tfuncs(2+AXTrig.nbin(), nullptr);
    for (int it = 1; it <= AXTrig.nbin(); ++it) {
        //tfuncs.at(it) = new TF1(Form("tfunc%03d", it), Form("[0]*(1+[1]/( 1 + exp(((x-%f)/86400. - [2]) / ([3] / 4.39)) ))", AXTtme.min()), AXTtme.min(), AXTtme.max());
        //tfuncs.at(it) = new TF1(Form("tfunc%03d", it), Form("[0]/(1+[1]/( 1 + exp(-((x-%f)/86400. - [2]) / ([3] / 4.39)) ))", AXTtme.min()), AXTtme.min(), AXTtme.max());
        tfuncs.at(it) = new TF1(Form("tfunc%03d", it), Form("[0] * (1.0 + [1] * cos(2.0 * 3.141596 * ((x-%f)/(365.25*86400.) - [2]) / ([3])))", AXTtme.min()), AXTtme.min(), AXTtme.max());
        //tfuncs.at(it) = new TF1(Form("tfunc%03d", it), Form("[0] * (1.0 + [1] * cos(2.0 * 3.141596 * ((x-%f)/86400. - [2]) / (700 / 4.39)))", AXTtme.min()), AXTtme.min(), AXTtme.max());
        //tfuncs.at(it) = new TF1(Form("tfunc%03d", it), Form("[0]*(1-[1]/( 1 + exp(-((x-%f)/86400. - [2]) / ([3] / 4.39)) ))", AXTtme.min()), AXTtme.min(), AXTtme.max());
        TF1* tfunc = tfuncs.at(it);
        Hist* happ = vhTapp.at(it);
        tfunc->SetParameters((*happ)()->GetBinContent(1), 0.5, 2, 8);
        tfunc->SetParLimits(1, 0., 0.8);
        tfunc->SetParLimits(2, 0, 6);
        tfunc->SetParLimits(3, 4, 14);
        (*happ)()->Fit(tfunc, "q0", "", AXTtme.min(), AXTtme.max());
        (*happ)()->Fit(tfunc, "q0", "", AXTtme.min(), AXTtme.max());
        (*happ)()->Fit(tfunc, "q0", "", AXTtme.min(), AXTtme.max());
        //double dtime = tfunc->GetParameter(2) - (1372636800-AXTtme.min()) / 86400.;
        double dtime = tfunc->GetParameter(2);
        std::cerr << Form("R%2d PAR %14.8f %14.8f %14.8f %14.8f DTIME %14.8f\n", it, tfunc->GetParameter(0), tfunc->GetParameter(1), tfunc->GetParameter(2), tfunc->GetParameter(3), dtime);
        
        (*hTapp_dt)()->SetBinContent(it, tfunc->GetParameter(3));
        (*hTapp_dt)()->SetBinError  (it, tfunc->GetParError(3));

        (*hTapp_st)()->SetBinContent(it, dtime);
        (*hTapp_st)()->SetBinError  (it, tfunc->GetParError(2));

        (*hTapp_cc)()->SetBinContent(it, tfunc->GetParameter(1));
        (*hTapp_cc)()->SetBinError  (it, tfunc->GetParError(1));

        (*happ)()->GetYaxis()->SetTitle("#bar{p}/p");
        (*happ)()->GetXaxis()->SetTitle("");
        (*happ)()->GetXaxis()->SetTimeDisplay(true);
        (*happ)()->GetXaxis()->SetTimeOffset(0, "GMT");
        (*happ)()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
        happ->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 1.5));
        tfunc->SetLineColor(kBlue);

        editor2.create();
        happ->draw("pe");
        tfunc->Draw("l same");
    
        Legend leg("", TextStyle(kRed, 30, 43), PadWindow(0.6, 0.8, 0.1, 0.15));
        leg()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(it-1), AXTrig()(it)), "");
        leg.draw();

        editor2.save();
    } 
    
    Axis AXbins_dt("#Delta t (years)", 4000, 0, 15);
    Axis AXbins_st("t_{1/2}-t_{rev} (years)", 2000, 0, 6);
    Hist* hTapp2D_dt = Hist::New("hTapp2D_dt", HistAxis(AXTrig, AXbins_dt));
    Hist* hTapp2D_st = Hist::New("hTapp2D_st", HistAxis(AXTrig, AXbins_st));

    editor2.create("", PadMargin(0.1, 0.2, 0.1, 0.05));
    editor2.cd(1, PadAxis(1));
    (*hTapp2D_dt)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hTapp2D_dt)()->GetXaxis()->SetMoreLogLabels();
    (*hTapp2D_dt)()->GetXaxis()->CenterTitle();
    (*hTapp2D_dt)()->GetYaxis()->CenterTitle();
    (*hTapp2D_dt)()->GetXaxis()->SetRangeUser(1.0, 4.5);
    (*hTapp2D_dt)()->GetXaxis()->SetLabelSize(2 * (*hTapp_dt)()->GetXaxis()->GetLabelSize());
    (*hTapp2D_dt)()->GetYaxis()->SetLabelSize(2 * (*hTapp_dt)()->GetYaxis()->GetLabelSize());
    (*hTapp2D_dt)()->GetXaxis()->SetTitleSize(1.6 * (*hTapp_dt)()->GetXaxis()->GetTitleSize());
    (*hTapp2D_dt)()->GetYaxis()->SetTitleSize(1.6 * (*hTapp_dt)()->GetYaxis()->GetTitleSize());
    (*hTapp2D_dt)()->GetYaxis()->SetTitleOffset(0.55 * (*hTapp_dt)()->GetYaxis()->GetTitleOffset());
    (*hTapp2D_dt)()->GetYaxis()->SetRangeUser(0, 6000);
    hTapp2D_dt->draw();
    hTapp_dt->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 3));
    hTapp_dt->draw("pe same");
    editor2.save();

    editor2.create("", PadMargin(0.1, 0.2, 0.1, 0.05));
    editor2.cd(1, PadAxis(1));
    (*hTapp2D_st)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hTapp2D_st)()->GetXaxis()->SetMoreLogLabels();
    (*hTapp2D_st)()->GetXaxis()->CenterTitle();
    (*hTapp2D_st)()->GetYaxis()->CenterTitle();
    (*hTapp2D_st)()->GetXaxis()->SetRangeUser(1.0, 4.5);
    (*hTapp2D_st)()->GetXaxis()->SetLabelSize(2 * (*hTapp_st)()->GetXaxis()->GetLabelSize());
    (*hTapp2D_st)()->GetYaxis()->SetLabelSize(2 * (*hTapp_st)()->GetYaxis()->GetLabelSize());
    (*hTapp2D_st)()->GetXaxis()->SetTitleSize(1.6 * (*hTapp_st)()->GetXaxis()->GetTitleSize());
    (*hTapp2D_st)()->GetYaxis()->SetTitleSize(1.6 * (*hTapp_st)()->GetYaxis()->GetTitleSize());
    (*hTapp2D_st)()->GetYaxis()->SetTitleOffset(0.55 * (*hTapp_st)()->GetYaxis()->GetTitleOffset());
    (*hTapp2D_st)()->GetYaxis()->SetRangeUser(-6000, 6000);
    hTapp2D_st->draw();
    hTapp_st->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 3));
    hTapp_st->draw("pe same");
    editor2.save();
    
    editor2.create("", PadMargin(0.1, 0.2, 0.1, 0.05));
    editor2.cd(1, PadAxis(1));
    (*hTapp_cc)()->GetXaxis()->SetMoreLogLabels();
    (*hTapp_cc)()->GetXaxis()->CenterTitle();
    (*hTapp_cc)()->GetYaxis()->CenterTitle();
    (*hTapp_cc)()->GetXaxis()->SetRangeUser(1.0, 25);
    (*hTapp_cc)()->GetXaxis()->SetLabelSize(2 * (*hTapp_cc)()->GetXaxis()->GetLabelSize());
    (*hTapp_cc)()->GetYaxis()->SetLabelSize(2 * (*hTapp_cc)()->GetYaxis()->GetLabelSize());
    (*hTapp_cc)()->GetXaxis()->SetTitleSize(1.6 * (*hTapp_cc)()->GetXaxis()->GetTitleSize());
    (*hTapp_cc)()->GetYaxis()->SetTitleSize(1.6 * (*hTapp_cc)()->GetYaxis()->GetTitleSize());
    (*hTapp_cc)()->GetYaxis()->SetTitleOffset(0.55 * (*hTapp_cc)()->GetYaxis()->GetTitleOffset());
    hTapp_cc->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 3));
    //(*hTapp_cc)()->SetMininum(0);
    hTapp_cc->draw("pe");
    hTapp_cc_bl->draw("hist same");
    editor2.save();
    
    editor2.close();
    
    PdfEditor editor(Window(WindowSize::kA4Vertical), "apflux_time", "out", PdfEditor::Type::kPDF, PdfEditor::Mode::kMerge);

    TVirtualPad* pad = nullptr;
    editor.create("", 1, 4, PadMargin(0.03, 0.08, 0.1, 0.05));

    pad = editor.cd(1);
    //(*vhTCrrR.at(1))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(1)->draw("pe");
    tfuncs.at(1)->SetParameter(0, tfuncs.at(1)->GetParameter(0)*1.0e5);
    tfuncs.at(1)->Draw("l same");
    Legend leg1("", TextStyle(kRed, 75, 43), PadWindow(0.15, 0.4, 0.1, 0.25));
    leg1()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(0), AXTrig()(1)), "");
    leg1.draw();

    pad = editor.cd(2);
    (*vhTCrrR.at(2))()->GetYaxis()->SetTitle("#bar{p}/p #times 10^{5}");
    //(*vhTCrrR.at(2))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(2)->draw("pe");
    tfuncs.at(2)->SetParameter(0, tfuncs.at(2)->GetParameter(0)*1.0e5);
    tfuncs.at(2)->Draw("l same");
    Legend leg2("", TextStyle(kRed, 75, 43), PadWindow(0.15, 0.4, 0.1, 0.25));
    leg2()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(1), AXTrig()(2)), "");
    leg2.draw();
    
    pad = editor.cd(3);
    //(*vhTCrrR.at(3))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(3)->draw("pe");
    tfuncs.at(3)->SetParameter(0, tfuncs.at(3)->GetParameter(0)*1.0e5);
    tfuncs.at(3)->Draw("l same");
    Legend leg3("", TextStyle(kRed, 75, 43), PadWindow(0.15, 0.4, 0.1, 0.25));
    leg3()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(2), AXTrig()(3)), "");
    leg3.draw();
    
    pad = editor.cd(4);
    //(*vhTCrrR.at(4))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(5)->draw("pe");
    tfuncs.at(5)->SetParameter(0, tfuncs.at(5)->GetParameter(0)*1.0e5);
    tfuncs.at(5)->Draw("l same");
    Legend leg4("", TextStyle(kRed, 75, 43), PadWindow(0.15, 0.4, 0.1, 0.25));
    leg4()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(4), AXTrig()(5)), "");
    leg4.draw();
    
    editor.save();
    editor.close();


    TFile * ofle = new TFile("out/apflux_time.root", "RECREATE");
    ofle->cd();

    (*hTCrrR)()->Write();
    for (int it = 1; it < vhTCrrR.size()-1; ++it) {
        (*vhTCrrR.at(it))()->Write();
    }

    (*hTapp)()->Write();
    for (int it = 1; it < vhTapp.size()-1; ++it) {
        (*vhTapp.at(it))()->Write();
    }

    ofle->Write();
    ofle->Close();

    return 1;
}
