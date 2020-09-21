#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "DataFit2D.h"
#include "DataFit2D.C"

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
    
    Hist* hTCrrR = Hist::New("hTCrrR", (TH1*)TFile::Open("out_Dec04/antip_time_v3.root")->Get("hTCrrR"));
    std::vector<Hist*> vhTCrrR = Hist::ProjectAll(HistProj::kX, Hist::Head("hTCrrR"));
    
    for (int it = 1; it < vhTCrrR.size()-1; ++it) {
        //(*vhTCrrR.at(it))()->GetXaxis()->SetTitle("Date");
        (*vhTCrrR.at(it))()->GetXaxis()->SetTitle("");
        (*vhTCrrR.at(it))()->GetXaxis()->SetTimeDisplay(true);
        (*vhTCrrR.at(it))()->GetXaxis()->SetTimeOffset(0, "GMT");
        (*vhTCrrR.at(it))()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
        (*vhTCrrR.at(it))()->GetXaxis()->SetLabelFont(43);
        (*vhTCrrR.at(it))()->GetXaxis()->SetLabelSize(35);
        (*vhTCrrR.at(it))()->GetXaxis()->SetTitleFont(43);
        (*vhTCrrR.at(it))()->GetXaxis()->SetTitleSize(100);
        (*vhTCrrR.at(it))()->GetXaxis()->SetTitleOffset(1.5);
        (*vhTCrrR.at(it))()->GetXaxis()->CenterTitle();
        //(*vhTCrrR.at(it))()->GetYaxis()->SetTitle("#bar{p}/p");
        (*vhTCrrR.at(it))()->GetYaxis()->SetTitle("");
        (*vhTCrrR.at(it))()->GetYaxis()->SetLabelFont(43);
        (*vhTCrrR.at(it))()->GetYaxis()->SetLabelSize(35);
        (*vhTCrrR.at(it))()->GetYaxis()->SetTitleFont(43);
        (*vhTCrrR.at(it))()->GetYaxis()->SetTitleSize(100);
        (*vhTCrrR.at(it))()->GetYaxis()->SetTitleOffset(1.5);
        (*vhTCrrR.at(it))()->GetYaxis()->SetNoExponent();
        (*vhTCrrR.at(it))()->GetYaxis()->CenterTitle();
        vhTCrrR.at(it)->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 1.5));
        (*vhTCrrR.at(it))()->Scale(1.0e+5);
    }
    
    const Axis& AXTtme = Hist::Head("hTCrrR")->xaxis();
    const Axis& AXTrig = Hist::Head("hTCrrR")->yaxis();

    PdfEditor editor(Window(WindowSize::kA4Vertical), "antip_time", "out_Dec04", PdfEditor::Type::kPDF, PdfEditor::Mode::kMergeAndAlone);
    //PdfEditor editor(Window(WindowSize::kA4kHorizon), "antip_time", "out_Dec04", PdfEditor::Type::kPDF, PdfEditor::Mode::kMergeAndAlone);

    TVirtualPad* pad = nullptr;
    editor.create("", 1, 8, PadMargin(0.02, 0.07, 0.1, 0.05));

    pad = editor.cd(1);
    //(*vhTCrrR.at(1))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(1)->draw("pe");
    Legend leg1("", TextStyle(kRed, 50, 43), PadWindow(0.2, 0.6, 0.1, 0.25));
    leg1()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(0), AXTrig()(1)), "");
    leg1.draw();

    pad = editor.cd(2);
    //(*vhTCrrR.at(2))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(2)->draw("pe");
    Legend leg2("", TextStyle(kRed, 50, 43), PadWindow(0.2, 0.6, 0.1, 0.25));
    leg2()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(1), AXTrig()(2)), "");
    leg2.draw();
    
    pad = editor.cd(3);
    //(*vhTCrrR.at(3))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(3)->draw("pe");
    Legend leg3("", TextStyle(kRed, 50, 43), PadWindow(0.2, 0.6, 0.1, 0.25));
    leg3()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(2), AXTrig()(3)), "");
    leg3.draw();
    
    pad = editor.cd(4);
    //(*vhTCrrR.at(4))()->GetXaxis()->SetLabelSize(0);
    (*vhTCrrR.at(4))()->GetYaxis()->SetTitle("#bar{p}/p #times 10^{5}");
    vhTCrrR.at(4)->draw("pe");
    Legend leg4("", TextStyle(kRed, 50, 43), PadWindow(0.2, 0.6, 0.1, 0.25));
    leg4()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(3), AXTrig()(4)), "");
    leg4.draw();
    
    pad = editor.cd(5);
    //(*vhTCrrR.at(6))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(6)->draw("pe");
    Legend leg5("", TextStyle(kRed, 50, 43), PadWindow(0.2, 0.6, 0.1, 0.25));
    leg5()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(5), AXTrig()(6)), "");
    leg5.draw();
    
    pad = editor.cd(6);
    //(*vhTCrrR.at(8))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(8)->draw("pe");
    Legend leg6("", TextStyle(kRed, 50, 43), PadWindow(0.2, 0.6, 0.1, 0.25));
    leg6()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(7), AXTrig()(8)), "");
    leg6.draw();
    
    pad = editor.cd(7);
    //(*vhTCrrR.at(11))()->GetXaxis()->SetLabelSize(0);
    vhTCrrR.at(11)->draw("pe");
    Legend leg7("", TextStyle(kRed, 50, 43), PadWindow(0.2, 0.6, 0.1, 0.25));
    leg7()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(10), AXTrig()(11)), "");
    leg7.draw();
    
    pad = editor.cd(8);
    //(*vhTCrrR.at(15))()->GetXaxis()->SetLabelSize(35);
    vhTCrrR.at(15)->draw("pe");
    Legend leg8("", TextStyle(kRed, 50, 43), PadWindow(0.2, 0.6, 0.1, 0.25));
    leg8()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(14), AXTrig()(15)), "");
    leg8.draw();
    editor.save();
    editor.close();
    
    Hist* hTpap = Hist::New("hTpap", HistAxis(AXTtme, AXTrig));
    for (int it = 1; it <= AXTtme.nbin(); ++it) {
    for (int jt = 1; jt <= AXTrig.nbin(); ++jt) {
        double val = 1.0 / (*hTCrrR)()->GetBinContent(it, jt);
        double err = (*hTCrrR)()->GetBinError(it, jt) / (*hTCrrR)()->GetBinContent(it, jt);
        if (!std::isfinite(val)) continue;
        if (!std::isfinite(err)) continue;

        (*hTpap)()->SetBinContent(it, jt, val);
        (*hTpap)()->SetBinError  (it, jt, val * err);
    }}
    std::vector<Hist*> vhTpap = Hist::ProjectAll(HistProj::kX, Hist::Head("hTpap"));
   
    Hist* hTpap_dt = Hist::New("hTpap_dt", HistAxis(AXTrig));
    Hist* hTpap_st = Hist::New("hTpap_st", HistAxis(AXTrig));
    Hist* hTpap_cc = Hist::New("hTpap_cc", HistAxis(AXTrig));
    Hist* hTpap_cc_bl = Hist::New("hTpap_cc_baseline", HistAxis(AXTrig));
    (*hTpap_dt)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hTpap_dt)()->GetYaxis()->SetTitle("#Delta t (days)");
    (*hTpap_st)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hTpap_st)()->GetYaxis()->SetTitle("t_{1/2}-t_{rev} (days)");
    (*hTpap_cc)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hTpap_cc)()->GetYaxis()->SetTitle("C");

    PdfEditor editor2(Window(), "antip_time_out", "out_Dec04", PdfEditor::Type::kPDF, PdfEditor::Mode::kMerge);

    TF1* tfunc = new TF1("tfunc", Form("[0]*(1+[1]/( 1 + exp(-((x-%f)/86400. - [2]) / ([3] / 4.39)) ))", AXTtme.min()), AXTtme.min(), AXTtme.max());
    for (int it = 1; it <= AXTrig.nbin(); ++it) {
        Hist* hpap = vhTpap.at(it);
        tfunc->SetParameters((*hpap)()->GetBinContent(1), 0.5, 1400, 700);
        tfunc->SetParLimits(1, -0.2, 2.0);
        tfunc->SetParLimits(2, 700, 2500);
        tfunc->SetParLimits(3,  50, 4000);
        (*hpap)()->Fit(tfunc, "q0", "", AXTtme.min(), AXTtme.max());
        double dtime = tfunc->GetParameter(2) - (1372636800-1305417600) / 86400.;
        std::cerr << Form("R%2d PAR %14.8f %14.8f %14.8f %14.8f DTIME %14.8f\n", it, tfunc->GetParameter(0), tfunc->GetParameter(1), tfunc->GetParameter(2), tfunc->GetParameter(3), dtime);
        
        (*hTpap_dt)()->SetBinContent(it, tfunc->GetParameter(3));
        (*hTpap_dt)()->SetBinError  (it, tfunc->GetParError(3));

        (*hTpap_st)()->SetBinContent(it, dtime);
        (*hTpap_st)()->SetBinError  (it, tfunc->GetParError(2));

        (*hTpap_cc)()->SetBinContent(it, tfunc->GetParameter(1));
        (*hTpap_cc)()->SetBinError  (it, tfunc->GetParError(1));

        (*hpap)()->GetYaxis()->SetTitle("p/#bar{p}");
        (*hpap)()->GetXaxis()->SetTitle("");
        (*hpap)()->GetXaxis()->SetTimeDisplay(true);
        (*hpap)()->GetXaxis()->SetTimeOffset(0, "GMT");
        (*hpap)()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
        hpap->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 1.5));
        tfunc->SetLineColor(kBlue);

        editor2.create();
        hpap->draw("pe");
        tfunc->Draw("l same");
    
        Legend leg("", TextStyle(kRed, 30, 43), PadWindow(0.6, 0.8, 0.1, 0.15));
        leg()->AddEntry((TObject*)0, Form("[%5.2f ~ %5.2f] GV", AXTrig()(it-1), AXTrig()(it)), "");
        leg.draw();

        editor2.save();
    } 

    editor2.create();
    editor2.cd(1, PadAxis(1));
    (*hTpap_dt)()->GetXaxis()->SetMoreLogLabels();
    (*hTpap_dt)()->GetXaxis()->SetRangeUser(1.0, 5.5);
    hTpap_dt->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 1.5));
    hTpap_dt->draw("pe");
    editor2.save();

    editor2.create();
    editor2.cd(1, PadAxis(1));
    (*hTpap_st)()->GetXaxis()->SetMoreLogLabels();
    (*hTpap_st)()->GetXaxis()->SetRangeUser(1.0, 5.5);
    hTpap_st->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 1.5));
    hTpap_st->draw("pe");
    editor2.save();
    
    editor2.create();
    editor2.cd(1, PadAxis(1));
    (*hTpap_cc)()->GetXaxis()->SetMoreLogLabels();
    hTpap_cc->style(Line(kRed, 0, 1), Marker(kRed, MarkerStyle(MarkerShape::kCircle), 1.5));
    hTpap_cc->draw("pe");
    hTpap_cc_bl->draw("hist same");
    editor2.save();
    
    editor2.close();

    TFile * ofle = new TFile("out_Dec04/antip_time_out.root", "RECREATE");
    ofle->cd();

    (*hTCrrR)()->Write();
    for (int it = 1; it < vhTCrrR.size()-1; ++it) {
        (*vhTCrrR.at(it))()->Write();
    }

    (*hTpap)()->Write();
    for (int it = 1; it < vhTpap.size()-1; ++it) {
        (*vhTpap.at(it))()->Write();
    }

    ofle->Write();
    ofle->Close();

    return 1;
}
