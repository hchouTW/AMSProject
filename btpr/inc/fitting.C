#include <CPPLibs.h>
#include <ROOTLibs.h>
#include <TrSys.h>

double LinkSMG3func(double *x, double *p) {
    return TrSys::MultiGausFunc(x[0], 0.0, std::vector<std::array<double, 2>>{ {p[0], p[1]}, {p[2], p[3]}, {p[4], p[5]} });
}
double LinkSMG4func(double *x, double *p) {
    return TrSys::MultiGausFunc(x[0], 0.0, std::vector<std::array<double, 2>>{ {p[0], p[1]}, {p[2], p[3]}, {p[4], p[5]}, {p[6], p[7]} });
}
double LinkALGfunc(double *x, double *p) {
    return (p[4] * TrSys::LandauGausFunc(x[0], p[0], p[1], p[2], p[3]));
}

int main() {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    
    Hist::Load("YiMdst.root", "/eos/ams/user/h/hchou/AMSData/subj/btpr/20Jun11/btpr21");
    Hist* hEtd  = Hist::Head("hEtd");
    Hist* hEtf  = Hist::Head("hEtf");
    Hist* hEtk  = Hist::Head("hEtk");
    Hist* hEtkx = Hist::Head("hEtkx");
    Hist* hEtky = Hist::Head("hEtky");
    Hist* hEtkxy = Hist::Head("hEtkxy");
    
    hEtd ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hEtf ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hEtk ->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hEtkx->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hEtky->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hEtkxy->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));

    TF1* ALGfunc = new TF1("ALGfunc", LinkALGfunc, 0.0, 100.0, 5);
    ALGfunc->SetLineColor(kBlue);
    ALGfunc->SetNpx(100000);

    ALGfunc->SetParameter(0, 0.5);
    ALGfunc->SetParameter(1, 1.0);
    ALGfunc->SetParameter(2, 0.5);
    ALGfunc->SetParameter(3, 0.1);
    ALGfunc->SetParameter(4, 10000.0);
    ALGfunc->SetParLimits(0, 0.0, 1.0);
    ALGfunc->SetParLimits(1, 0.0, 5.0);
    ALGfunc->SetParLimits(2, 0.01, 10.0);
    ALGfunc->SetParLimits(3, 0.00, 10.0);

    (*hEtf)()->Fit(ALGfunc, "q0", "");
    (*hEtf)()->Fit(ALGfunc, "q0", "");
    (*hEtf)()->Fit(ALGfunc, "q0", "");

    std::cerr << Form("PARAM[0]  %14.8f %14.8f\n", ALGfunc->GetParameter(0), ALGfunc->GetParError(0)); 
    std::cerr << Form("PARAM[1]  %14.8f %14.8f\n", ALGfunc->GetParameter(1), ALGfunc->GetParError(1)); 
    std::cerr << Form("PARAM[2]  %14.8f %14.8f\n", ALGfunc->GetParameter(2), ALGfunc->GetParError(2)); 
    std::cerr << Form("PARAM[3]  %14.8f %14.8f\n", ALGfunc->GetParameter(3), ALGfunc->GetParError(3)); 

    PdfEditor editor(Window(), "fitting", "out");
    //PdfEditor editor(Window(WindowSize::kWideSliceLR), "fitting", "out");

    editor.create();
    editor.cd(0, PadAxis(0, 1));
    (*hEtf)()->Draw("hist");
    ALGfunc->Draw("l same");
    editor.save();
    
    TF1* Gfunc = new TF1("Gfunc", "gaus", 0.0, 5.0);
    Gfunc->SetParameters(1.0, 0.0, 1.0);

    Hist::Head("hRrso_ckIn")    ->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hRrso_ckL1")    ->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hRrso_ckL9")    ->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hRrso_ckFs")    ->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hRrso_hcIn_new")->style(Line(kRed, 0, 2),  Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hRrso_hcL1_new")->style(Line(kRed, 0, 2),  Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hRrso_hcL9_new")->style(Line(kRed, 0, 2),  Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hRrso_hcFs_new")->style(Line(kRed, 0, 2),  Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    (*Hist::Head("hRrso_ckIn"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_ckIn"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_ckIn"))()->Fit(Gfunc, "q0", "");
    double RSgmCKin = Gfunc->GetParameter(2);
    
    (*Hist::Head("hRrso_ckL1"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_ckL1"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_ckL1"))()->Fit(Gfunc, "q0", "");
    double RSgmCKl1 = Gfunc->GetParameter(2);
    
    (*Hist::Head("hRrso_ckL9"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_ckL9"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_ckL9"))()->Fit(Gfunc, "q0", "");
    double RSgmCKl9 = Gfunc->GetParameter(2);
    
    (*Hist::Head("hRrso_ckFs"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_ckFs"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_ckFs"))()->Fit(Gfunc, "q0", "");
    double RSgmCKfs = Gfunc->GetParameter(2);
    
    (*Hist::Head("hRrso_hcIn_new"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_hcIn_new"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_hcIn_new"))()->Fit(Gfunc, "q0", "");
    double RSgmHCin = Gfunc->GetParameter(2);
    
    (*Hist::Head("hRrso_hcL1_new"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_hcL1_new"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_hcL1_new"))()->Fit(Gfunc, "q0", "");
    double RSgmHCl1 = Gfunc->GetParameter(2);
    
    (*Hist::Head("hRrso_hcL9_new"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_hcL9_new"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_hcL9_new"))()->Fit(Gfunc, "q0", "");
    double RSgmHCl9 = Gfunc->GetParameter(2);
    
    (*Hist::Head("hRrso_hcFs_new"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_hcFs_new"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hRrso_hcFs_new"))()->Fit(Gfunc, "q0", "");
    double RSgmHCfs = Gfunc->GetParameter(2);
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    (*Hist::Head("hRrso_hcIn_new"))()->GetYaxis()->SetTitle("Events / Bin");
    Hist::Head("hRrso_hcIn_new")->draw("hist");
    Hist::Head("hRrso_ckIn"    )->draw("hist same");
    Legend leg_table_rin("", TextStyle(kBlack, 25, 43), PadWindow(0.65, 0.85, 0.70, 0.90));
    leg_table_rin()->AddEntry((*Hist::Head("hRrso_hcIn_new"))(), Form("New (#sigma=%5.3f)", RSgmHCin), "lp");
    leg_table_rin()->AddEntry((*Hist::Head("hRrso_ckIn"    ))(), Form("Official (#sigma=%5.3f)", RSgmCKin)     , "lp");
    leg_table_rin()->AddEntry((TObject*)0, Form("Improvement ~%5.1f%", 100.0 * (1.0 - RSgmHCin/RSgmCKin)), "");
    leg_table_rin()->SetFillColor(0);
    leg_table_rin.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    (*Hist::Head("hRrso_hcL1_new"))()->GetYaxis()->SetTitle("Events / Bin");
    Hist::Head("hRrso_hcL1_new")->draw("hist");
    Hist::Head("hRrso_ckL1"    )->draw("hist same");
    Legend leg_table_rl1("", TextStyle(kBlack, 25, 43), PadWindow(0.65, 0.85, 0.70, 0.90));
    leg_table_rl1()->AddEntry((*Hist::Head("hRrso_hcL1_new"))(), Form("New (#sigma=%5.3f)", RSgmHCl1), "lp");
    leg_table_rl1()->AddEntry((*Hist::Head("hRrso_ckL1"    ))(), Form("Official (#sigma=%5.3f)", RSgmCKl1)     , "lp");
    leg_table_rl1()->AddEntry((TObject*)0, Form("Improvement ~%5.1f%", 100.0 * (1.0 - RSgmHCl1/RSgmCKl1)), "");
    leg_table_rl1()->SetFillColor(0);
    leg_table_rl1.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    (*Hist::Head("hRrso_hcL9_new"))()->GetYaxis()->SetTitle("Events / Bin");
    Hist::Head("hRrso_hcL9_new")->draw("hist");
    Hist::Head("hRrso_ckL9"    )->draw("hist same");
    Legend leg_table_rl9("", TextStyle(kBlack, 25, 43), PadWindow(0.65, 0.85, 0.70, 0.90));
    leg_table_rl9()->AddEntry((*Hist::Head("hRrso_hcL9_new"))(), Form("New (#sigma=%5.3f)", RSgmHCl9), "lp");
    leg_table_rl9()->AddEntry((*Hist::Head("hRrso_ckL9"    ))(), Form("Official (#sigma=%5.3f)", RSgmCKl9)     , "lp");
    leg_table_rl9()->AddEntry((TObject*)0, Form("Improvement ~%5.1f%", 100.0 * (1.0 - RSgmHCl9/RSgmCKl9)), "");
    leg_table_rl9()->SetFillColor(0);
    leg_table_rl9.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    (*Hist::Head("hRrso_hcFs_new"))()->GetYaxis()->SetTitle("Events / Bin");
    Hist::Head("hRrso_hcFs_new")->draw("hist");
    Hist::Head("hRrso_ckFs"    )->draw("hist same");
    Legend leg_table_rfs("", TextStyle(kBlack, 25, 43), PadWindow(0.65, 0.85, 0.70, 0.90));
    leg_table_rfs()->AddEntry((*Hist::Head("hRrso_hcFs_new"))(), Form("New (#sigma=%5.3f)", RSgmHCfs), "lp");
    leg_table_rfs()->AddEntry((*Hist::Head("hRrso_ckFs"    ))(), Form("Official (#sigma=%5.3f)", RSgmCKfs)     , "lp");
    leg_table_rfs()->AddEntry((TObject*)0, Form("Improvement ~%5.1f%", 100.0 * (1.0 - RSgmHCfs/RSgmCKfs)), "");
    leg_table_rfs()->SetFillColor(0);
    leg_table_rfs.draw();
    editor.save();
    
    Gfunc->SetParameters(1.0, 1.0, 0.1);
    Hist::Head("hQtk_alg")->style(Line(kRed, 0, 2),  Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hQtk_off")->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kCircle)));
    (*Hist::Head("hQtk_alg"))()->GetYaxis()->SetTitle("Events / Bin");
    (*Hist::Head("hQtk_alg"))()->GetXaxis()->SetRangeUser(0.5, 2.5);

    (*Hist::Head("hQtk_alg"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hQtk_alg"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hQtk_alg"))()->Fit(Gfunc, "q0", "");
    double QtkSgmAlg = Gfunc->GetParameter(2);
    
    (*Hist::Head("hQtk_off"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hQtk_off"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hQtk_off"))()->Fit(Gfunc, "q0", "");
    double QtkSgmOff = Gfunc->GetParameter(2);

    Hist::Head("hQtf_alg")->style(Line(kRed, 0, 2),  Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    Hist::Head("hQtf_off")->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kCircle)));
    (*Hist::Head("hQtf_alg"))()->GetYaxis()->SetTitle("Events / Bin");
    (*Hist::Head("hQtf_alg"))()->GetXaxis()->SetRangeUser(0.5, 3.0);
    
    (*Hist::Head("hQtf_alg"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hQtf_alg"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hQtf_alg"))()->Fit(Gfunc, "q0", "");
    double QtfSgmAlg = Gfunc->GetParameter(2);
    
    (*Hist::Head("hQtf_off"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hQtf_off"))()->Fit(Gfunc, "q0", "");
    (*Hist::Head("hQtf_off"))()->Fit(Gfunc, "q0", "");
    double QtfSgmOff = Gfunc->GetParameter(2);
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    Hist::Head("hQtk_alg")->draw("hist");
    Hist::Head("hQtk_off")->draw("hist same");
    Legend leg_table_tk("", TextStyle(kBlack, 30, 43), PadWindow(0.50, 0.85, 0.70, 0.90));
    leg_table_tk()->AddEntry((*Hist::Head("hQtk_alg"))(), Form("New (#sigma=%5.3f)", QtkSgmAlg), "lp");
    leg_table_tk()->AddEntry((*Hist::Head("hQtk_off"))(), Form("Official (#sigma=%5.3f)", QtkSgmOff)     , "lp");
    leg_table_tk()->AddEntry((TObject*)0, Form("Improvement ~%5.1f%", 100.0 * (1.0 - QtkSgmAlg/QtkSgmOff)), "");
    leg_table_tk()->SetFillColor(0);
    leg_table_tk.draw();
    editor.save();

    editor.create();
    editor.cd(0, PadAxis(0, 1));
    Hist::Head("hQtf_alg")->draw("hist");
    Hist::Head("hQtf_off")->draw("hist same");
    Legend leg_table_tf("", TextStyle(kBlack, 30, 43), PadWindow(0.50, 0.85, 0.70, 0.90));
    leg_table_tf()->AddEntry((*Hist::Head("hQtf_alg"))(), Form("New (#sigma=%5.3f)", QtfSgmAlg), "lp");
    leg_table_tf()->AddEntry((*Hist::Head("hQtf_off"))(), Form("Official (#sigma=%5.3f)", QtfSgmOff)     , "lp");
    leg_table_tf()->AddEntry((TObject*)0, Form("Improvement ~%5.1f%", 100.0 * (1.0 - QtfSgmAlg/QtfSgmOff)), "");
    leg_table_tf()->SetFillColor(0);
    leg_table_tf.draw();
    editor.save();
    
    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    (*Hist::Head("hQtk_comp"))()->GetXaxis()->SetRangeUser(0.5, 2.5);
    (*Hist::Head("hQtk_comp"))()->GetYaxis()->SetRangeUser(0.5, 2.5);
    (*Hist::Head("hQtk_comp"))()->GetXaxis()->SetTitle("TRK Q [New]");
    (*Hist::Head("hQtk_comp"))()->GetYaxis()->SetTitle("TRK Q [Official]");
    Hist::Head("hQtk_comp")->draw("colz");
    editor.save();
    
    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    (*Hist::Head("hQtf_comp"))()->GetXaxis()->SetRangeUser(0.5, 3.0);
    (*Hist::Head("hQtf_comp"))()->GetYaxis()->SetRangeUser(0.5, 3.0);
    (*Hist::Head("hQtf_comp"))()->GetXaxis()->SetTitle("TOF Q [New]");
    (*Hist::Head("hQtf_comp"))()->GetYaxis()->SetTitle("TOF Q [Official]");
    Hist::Head("hQtf_comp")->draw("colz");
    editor.save();
    
    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    (*Hist::Head("hQtk_lchi"))()->GetXaxis()->SetRangeUser( 0.5, 2.5);
    (*Hist::Head("hQtk_lchi"))()->GetYaxis()->SetRangeUser(-4.0, 5.0);
    Hist::Head("hQtk_lchi")->draw("colz");
    editor.save();

    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    (*Hist::Head("hQtf_lchi"))()->GetXaxis()->SetRangeUser( 0.5, 3.0);
    (*Hist::Head("hQtf_lchi"))()->GetYaxis()->SetRangeUser(-4.0, 5.0);
    Hist::Head("hQtf_lchi")->draw("colz");
    editor.save();

    editor.close();


    return 1;
}
