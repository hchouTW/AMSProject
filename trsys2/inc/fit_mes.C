#include <CPPLibs.h>
#include <ROOTLibs.h>
#include <TrSys.h>

void stdfmt(TH1* hist) {
    if (hist == nullptr) return;
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleFont(43);
    hist->GetXaxis()->SetTitleSize(20);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelSize(15);
    hist->GetYaxis()->SetTitleFont(43);
    hist->GetYaxis()->SetTitleSize(20);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelSize(15);
}

double mgaus3(double* x, double* par) {
    return par[7] * TrSys::MultiGausFunc(x[0], par[6], { {par[0], par[1]}, {par[2], par[3]}, {par[4], par[5]}} );
}
double mgaus4(double* x, double* par) {
    return par[9] * TrSys::MultiGausFunc(x[0], par[8], { {par[0], par[1]}, {par[2], par[3]}, {par[4], par[5]}, {par[6], par[7]}} );
}
double landgaus(double* x, double* par) {
    return par[4] * TrSys::LandauGausFunc(x[0], par[0], par[1], par[2], par[3]);
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    TH1D* hresx_ref = (TH1D*) TFile::Open("out/doc_mes/mes_resx.root")->Get("hresx");
    TH1D* hresy_ref = (TH1D*) TFile::Open("out/doc_mes/mes_resy.root")->Get("hresy");
    TH1D* htfql_ref = (TH1D*) TFile::Open("out/doc_mes/mes_tfql.root")->Get("htfql");
    TH1D* htfqm_ref = (TH1D*) TFile::Open("out/doc_mes/mes_tfqm.root")->Get("htfqm");
    TH1D* htfqh_ref = (TH1D*) TFile::Open("out/doc_mes/mes_tfqh.root")->Get("htfqh");

    hresx_ref->Sumw2();
    hresy_ref->Sumw2();
    htfql_ref->Sumw2();
    htfqm_ref->Sumw2();
    htfqh_ref->Sumw2();

    hresx_ref->Scale(1.0/hresx_ref->Integral());
    hresy_ref->Scale(1.0/hresy_ref->Integral());
    htfql_ref->Scale(1.0/htfql_ref->Integral());
    htfqm_ref->Scale(1.0/htfqm_ref->Integral());
    htfqh_ref->Scale(1.0/htfqh_ref->Integral());

    Hist* hresx = Hist::New("resx", hresx_ref);
    Hist* hresy = Hist::New("resy", hresy_ref);
    Hist* htfql = Hist::New("tfql", htfql_ref);
    Hist* htfqm = Hist::New("tfqm", htfqm_ref);
    Hist* htfqh = Hist::New("tfqh", htfqh_ref);
    stdfmt((*hresx)());
    stdfmt((*hresy)());
    stdfmt((*htfql)());
    stdfmt((*htfqm)());
    stdfmt((*htfqh)());

    (*hresx)()->SetLineColor(kRed);
    (*hresx)()->SetMarkerColor(kRed);
    (*hresy)()->SetLineColor(kRed);
    (*hresy)()->SetMarkerColor(kRed);
    (*htfql)()->SetLineColor(kRed);
    (*htfql)()->SetMarkerColor(kRed);
    (*htfqm)()->SetLineColor(kRed);
    (*htfqm)()->SetMarkerColor(kRed);
    (*htfqh)()->SetLineColor(kRed);
    (*htfqh)()->SetMarkerColor(kRed);

    PdfEditor editor(Window(), "fit_mes", "out/doc_mes", PdfEditor::Type::kEPS, PdfEditor::Mode::kMergeAndAlone);
    
    TF1* mgausx = new TF1("mgausx", mgaus3, -1, 1, 8);
    mgausx->SetLineColor(kBlue);
    mgausx->SetNpx(10000);

    mgausx->SetParameters(1, 0.001, 0.3, 0.003, 0.1, 0.010, 0.0, 1000);
    mgausx->SetParLimits(0, 0.0001, 1.0);
    mgausx->SetParLimits(1, 0.0001, 0.1);
    mgausx->SetParLimits(2, 0.0001, 1.0);
    mgausx->SetParLimits(3, 0.0001, 0.1);
    mgausx->SetParLimits(4, 0.0001, 1.0);
    mgausx->SetParLimits(5, 0.0001, 0.1);
    mgausx->SetParLimits(6, -0.01,  0.01);
    hresx_ref->Fit(mgausx, "q0", "");
    hresx_ref->Fit(mgausx, "q0", "");
    hresx_ref->Fit(mgausx, "q0", "");
    hresx_ref->Fit(mgausx, "", "");

    TF1* mgausy = new TF1("mgausy", mgaus4, -1, 1, 10);
    mgausy->SetLineColor(kBlue);
    mgausy->SetNpx(10000);

    mgausy->SetParameters(1, 0.001, 0.3, 0.003, 0.1, 0.010, 0.02, 0.05, 0.0, 1000);
    mgausy->SetParLimits(0, 0.0001, 1.0);
    mgausy->SetParLimits(1, 0.0001, 0.1);
    mgausy->SetParLimits(2, 0.0001, 1.0);
    mgausy->SetParLimits(3, 0.0001, 0.1);
    mgausy->SetParLimits(4, 0.0001, 1.0);
    mgausy->SetParLimits(5, 0.0001, 0.1);
    mgausy->SetParLimits(6, 0.0001, 1.0);
    mgausy->SetParLimits(7, 0.0001, 0.1);
    mgausy->SetParLimits(8, -0.01,  0.01);
    hresy_ref->Fit(mgausy, "q0", "");
    hresy_ref->Fit(mgausy, "q0", "");
    hresy_ref->Fit(mgausy, "q0", "");
    hresy_ref->Fit(mgausy, "", "");

    TF1* lgausl = new TF1("lgausl", landgaus, 0, 10, 5);
    lgausl->SetLineColor(kBlue);
    lgausl->SetNpx(10000);

    lgausl->SetParameters(0.5, 2.1, 0.3, 0.1, 1);
    lgausl->SetParLimits(0, 0.0, 1.0);
    lgausl->SetParLimits(1, 0.01, 10.0);
    lgausl->SetParLimits(2, 0.01, 10.0);
    lgausl->SetParLimits(3, 0.01, 10.0);
    htfql_ref->Fit(lgausl, "q0", "");
    htfql_ref->Fit(lgausl, "q0", "");
    htfql_ref->Fit(lgausl, "q0", "");
    htfql_ref->Fit(lgausl, "", "");

    TF1* lgausm = new TF1("lgausm", landgaus, 0, 10, 5);
    lgausm->SetLineColor(kBlue);
    lgausm->SetNpx(10000);

    lgausm->SetParameters(0.10, 1.4, 0.3, 0.1, 1);
    lgausm->SetParLimits(0, 0.0, 1.0);
    lgausm->SetParLimits(1, 0.01, 10.0);
    lgausm->SetParLimits(2, 0.01, 10.0);
    lgausm->SetParLimits(3, 0.01, 10.0);
    htfqm_ref->Fit(lgausm, "q0", "");
    htfqm_ref->Fit(lgausm, "q0", "");
    htfqm_ref->Fit(lgausm, "q0", "");
    htfqm_ref->Fit(lgausm, "", "");

    TF1* lgaush = new TF1("lgaush", landgaus, 0, 10, 5);
    lgaush->SetLineColor(kBlue);
    lgaush->SetNpx(10000);

    lgaush->SetParameters(0.05, 1.3, 0.3, 0.1, 1);
    lgaush->SetParLimits(0, 0.0, 1.0);
    lgaush->SetParLimits(1, 0.01, 10.0);
    lgaush->SetParLimits(2, 0.01, 10.0);
    lgaush->SetParLimits(3, 0.01, 10.0);
    htfqh_ref->Fit(lgaush, "q0", "");
    htfqh_ref->Fit(lgaush, "q0", "");
    htfqh_ref->Fit(lgaush, "q0", "");
    htfqh_ref->Fit(lgaush, "", "");

    editor.create("");
    editor.cd(0, PadAxis(0, 1));
    (*hresx)()->GetXaxis()->SetTitle("Residual in X [cm]");
    (*hresx)()->GetYaxis()->SetTitle("Probability");
    (*hresx)()->Draw("pe");
    mgausx->Draw("l same");
    Legend leg_resx("", TextStyle(kBlack, 20, 43), PadWindow(0.70, 0.90, 0.75, 0.90));
    leg_resx()->AddEntry((TObject*)0, "Tracker", "");
    leg_resx()->AddEntry((*hresx)(), "Data", "lp");
    leg_resx()->AddEntry(mgausx, "SMG", "l");
    leg_resx.draw();
    editor.save();

    editor.create("");
    editor.cd(0, PadAxis(0, 1));
    (*hresy)()->GetXaxis()->SetTitle("Residual in Y [cm]");
    (*hresy)()->GetYaxis()->SetTitle("Probability");
    (*hresy)()->Draw("pe");
    mgausy->Draw("l same");
    Legend leg_resy("", TextStyle(kBlack, 20, 43), PadWindow(0.70, 0.90, 0.75, 0.90));
    leg_resy()->AddEntry((TObject*)0, "Tracker", "");
    leg_resy()->AddEntry((*hresy)(), "Data", "lp");
    leg_resy()->AddEntry(mgausy, "SMG", "l");
    leg_resy.draw();
    editor.save();

    editor.create("");
    editor.cd(0, PadAxis(0, 1));
    (*htfql)()->GetXaxis()->SetTitle("dE/dX [MIP]");
    (*htfql)()->GetYaxis()->SetTitle("Probability");
    (*htfql)()->Draw("pe");
    lgausl->Draw("l same");
    Legend leg_tfql("", TextStyle(kBlack, 20, 43), PadWindow(0.70, 0.90, 0.75, 0.90));
    leg_tfql()->AddEntry((TObject*)0, "TOF", "");
    leg_tfql()->AddEntry((*htfql)(), "Data", "lp");
    leg_tfql()->AddEntry(lgausl, "ALG", "l");
    leg_tfql.draw();
    editor.save();
    
    editor.create("");
    editor.cd(0, PadAxis(0, 1));
    (*htfqm)()->GetXaxis()->SetTitle("dE/dX [MIP]");
    (*htfqm)()->GetYaxis()->SetTitle("Probability");
    (*htfqm)()->Draw("pe");
    lgausm->Draw("l same");
    Legend leg_tfqm("", TextStyle(kBlack, 20, 43), PadWindow(0.70, 0.90, 0.75, 0.90));
    leg_tfqm()->AddEntry((TObject*)0, "TOF", "");
    leg_tfqm()->AddEntry((*htfqm)(), "Data", "lp");
    leg_tfqm()->AddEntry(lgausm, "ALG", "l");
    leg_tfqm.draw();
    editor.save();

    editor.create("");
    editor.cd(0, PadAxis(0, 1));
    (*htfqh)()->GetXaxis()->SetTitle("dE/dX [MIP]");
    (*htfqh)()->GetYaxis()->SetTitle("Probability");
    (*htfqh)()->Draw("pe");
    lgaush->Draw("l same");
    Legend leg_tfqh("", TextStyle(kBlack, 20, 43), PadWindow(0.70, 0.90, 0.75, 0.90));
    leg_tfqh()->AddEntry((TObject*)0, "TOF", "");
    leg_tfqh()->AddEntry((*htfqh)(), "Data", "lp");
    leg_tfqh()->AddEntry(lgaush, "ALG", "l");
    leg_tfqh.draw();
    editor.save();

    editor.close();
    
    TFile * ofle = new TFile("out/doc_mes/mes.root", "RECREATE");
    ofle->cd();

    hresx->write();
    hresy->write();
    htfql->write();
    htfqm->write();
    htfqh->write();

    ofle->Write();
    ofle->Close();

    return 1;
}


/*
FCN=1113.36 FROM MIGRAD    STATUS=CONVERGED     111 CALLS         112 TOTAL
EDM=2.8215e-07    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           9.99196e-01   1.00365e-01   1.93522e-03   2.02760e-03
2  p1           1.88412e-03   1.03164e-05   2.53366e-06  -5.75100e-01
3  p2           3.89323e-01   1.78710e-02   6.03515e-05  -7.57883e-02
4  p3           3.32112e-03   3.25429e-05   3.07783e-06  -9.98096e-01
5  p4           2.71260e-02   1.73197e-03   1.96203e-05  -6.71387e-02
6  p5           6.98928e-03   8.04874e-05   1.05542e-05  -2.11808e-01
7  p6          -6.52651e-06   2.12891e-06   3.47000e-06  -1.87515e-01
8  p7           8.38696e-03   1.13485e-05   1.11949e-07   6.50503e+00


FCN=2549.16 FROM MIGRAD    STATUS=CONVERGED     164 CALLS         165 TOTAL
EDM=1.94554e-08    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           4.78682e-01   7.37566e-02   1.18522e-04   8.35804e-04
2  p1           1.37797e-03   1.83575e-05   3.11005e-06   1.47638e-01
3  p2           8.21875e-02   1.35245e-02   3.17199e-05   1.57474e-02
4  p3           2.72244e-03   3.56429e-05   6.16318e-06  -3.75867e-02
5  p4           2.85380e-03   4.95330e-04   1.72952e-05   2.66679e-03
6  p5           6.83937e-03   1.30370e-04   3.10273e-05  -4.23147e-03
7  p6           9.26431e-01   1.26574e-01   2.38865e-04  -1.59470e-03
8  p7           7.50895e-04   5.50307e-06   3.05995e-06   4.77938e-02
9  p8          -1.12822e-04   9.20708e-07   2.26516e-06  -3.81866e-02
10  p9           1.85680e-02   2.54554e-05   3.49628e-07   1.36540e+00


FCN=561.809 FROM MIGRAD    STATUS=CONVERGED      65 CALLS          66 TOTAL
EDM=1.43674e-07    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           4.07524e-01   2.21759e-02   6.69885e-05   1.56047e-03
2  p1           2.14736e+00   2.93272e-03   2.11130e-06   4.05238e-01
3  p2           2.40172e-01   5.45350e-03   3.15483e-06   2.77808e-01
4  p3           1.18389e-01   7.39729e-03   1.27997e-05   6.54317e-02
5  p4           1.10024e-02   1.65592e-04   3.33370e-07   2.61687e+00


FCN=559.187 FROM MIGRAD    STATUS=CONVERGED      59 CALLS          60 TOTAL
EDM=2.14883e-14    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           2.19427e-02   9.77297e-04  -0.00000e+00   2.56466e-05
2  p1           1.39837e+00   6.78995e-04  -0.00000e+00   1.18909e-05
3  p2           1.01428e-01   1.33337e-03  -0.00000e+00  -5.42759e-06
4  p3           1.13978e-01   1.36265e-03   0.00000e+00  -1.55312e-05
5  p4           1.83862e-02   1.91705e-04   0.00000e+00  -6.37058e-04


FCN=506.369 FROM MIGRAD    STATUS=CONVERGED      64 CALLS          65 TOTAL
EDM=5.00253e-10    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           4.82158e-03   1.85709e-04   1.48051e-05   9.32746e-04
2  p1           1.12452e+00   5.88167e-04   1.76892e-06   9.96169e-02
3  p2           8.38657e-02   7.49960e-04   2.40292e-06   6.49245e-02
4  p3           9.36490e-02   9.53060e-04   6.54752e-06   2.03025e-02
5  p4           2.04397e-02   1.67304e-04   6.19620e-07   2.52229e-01
*/
