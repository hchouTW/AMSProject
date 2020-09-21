#include <CPPLibs.h>
#include <ROOTLibs.h>
#include <TrSys.h>

//#include "inc/SMGFit.h"
//#include "inc/ALGFit.h"

double SMG2func(double *x, double *p) {
    return TrSys::MultiGausFunc(x[0], 0.0, std::vector<std::array<double, 2>>{ {p[0], p[1]}, {p[2], p[3]} });
}
double SMG4func(double *x, double *p) {
    return TrSys::MultiGausFunc(x[0], 0.0, std::vector<std::array<double, 2>>{ {p[0], p[1]}, {p[2], p[3]}, {p[4], p[5]}, {p[6], p[7]} });
}
double ALGfunc(double *x, double *p) {
    return TrSys::LandauGausFunc(x[0], p[0], p[1], p[2], p[3]);
}

int main() {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
  
    const long MIN_NUM = 3;
    const long MAX_NUM = 25;
    Axis AXNUM("N_{Hit}", MAX_NUM-MIN_NUM+1, MIN_NUM-0.5, MAX_NUM+0.5);
   
    TF1* gsfunc = new TF1("gsfunc", "gaus");

    Hist::Load("nogaus_hist.root", "/afs/cern.ch/user/h/hchou/AMSProject/nogaus/out");
    
    PdfEditor editor(Window(), "nogaus_fit", "out");
    //PdfEditor editor(Window(WindowSize::kWideSliceLR), "nogaus_fit", "out");
   
    std::cerr << Form("==== SMG : LX ====\n");
    const std::array<double, 6> SMGLXpars({
        9.99196e-01, 1.88412e-03,
        3.89323e-01, 3.32112e-03,
        2.71260e-02, 6.98928e-03
    });
    const std::array<double, 8> SMGLYpars({
        9.26431e-01, 7.50895e-04, 
        4.78682e-01, 1.37797e-03,
        8.21875e-02, 2.72244e-03,
        2.85380e-03, 6.83937e-03
    });
    Hist* hSMGLXstdM = Hist::New("hSMGLXstdM", HistAxis(AXNUM));
    Hist* hSMGLXstdS = Hist::New("hSMGLXstdS", HistAxis(AXNUM));
    Hist* hSMGLXsmgM = Hist::New("hSMGLXsmgM", HistAxis(AXNUM));
    Hist* hSMGLXsmgS = Hist::New("hSMGLXsmgS", HistAxis(AXNUM));
    Hist* hSMGLXrsoR = Hist::New("hSMGLXrsoR", HistAxis(AXNUM));
        
    hSMGLXstdM->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kSquare)));
    hSMGLXstdS->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kSquare)));
	hSMGLXsmgM->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle)));
	hSMGLXsmgS->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle)));
	hSMGLXrsoR->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle)));
    
    std::vector<Hist*> vhSMGLXstd = Hist::ProjectAll(HistProj::kY, Hist::Head("hSMGLYstd"));
    std::vector<Hist*> vhSMGLXsmg = Hist::ProjectAll(HistProj::kY, Hist::Head("hSMGLYsmg"));
    //std::vector<Hist*> vhSMGLXstd = Hist::ProjectAll(HistProj::kY, Hist::Head("hSMGLXstd"));
    //std::vector<Hist*> vhSMGLXsmg = Hist::ProjectAll(HistProj::kY, Hist::Head("hSMGLXsmg"));

    for (int in = 1; in <= AXNUM.nbin(); ++in) {
        Hist* hstd = vhSMGLXstd.at(in);
        Hist* hsmg = vhSMGLXsmg.at(in);
		
        hstd->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kSquare)));
		hsmg->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle)));

        gsfunc->SetParameters(1.0, 0.0, 0.0002);
        (*hstd)()->Fit(gsfunc, "q0", "");
        (*hstd)()->Fit(gsfunc, "q0", "");
        (*hSMGLXstdM)()->SetBinContent(in, 1.0e+4 * gsfunc->GetParameter(1));
        (*hSMGLXstdM)()->SetBinError  (in, 1.0e+4 * gsfunc->GetParError (1));
        (*hSMGLXstdS)()->SetBinContent(in, 1.0e+4 * gsfunc->GetParameter(2));
        (*hSMGLXstdS)()->SetBinError  (in, 1.0e+4 * gsfunc->GetParError (2));
        
        gsfunc->SetParameters(1.0, 0.0, 0.0002);
        (*hsmg)()->Fit(gsfunc, "q0", "");
        (*hsmg)()->Fit(gsfunc, "q0", "");
        (*hSMGLXsmgM)()->SetBinContent(in, 1.0e+4 * gsfunc->GetParameter(1));
        (*hSMGLXsmgM)()->SetBinError  (in, 1.0e+4 * gsfunc->GetParError (1));
        (*hSMGLXsmgS)()->SetBinContent(in, 1.0e+4 * gsfunc->GetParameter(2));
        (*hSMGLXsmgS)()->SetBinError  (in, 1.0e+4 * gsfunc->GetParError (2));
       
        (*hSMGLXrsoR)()->SetBinContent(in, (*hSMGLXsmgS)()->GetBinContent(in) / (*hSMGLXstdS)()->GetBinContent(in));

        editor.create();
        THStack* hstack = Hist::Collect(Form("hsSMGLX%02ld\n", MIN_NUM+in-1), HistList({hstd, hsmg}));
        hstack->Draw("nostack hist");
        hstack->GetHistogram()->SetLineColor(0);
        hstack->GetHistogram()->SetMarkerColor(0);
        hstack->GetXaxis()->SetTitle("Value");
        hstack->GetYaxis()->SetTitle("Events/Bin");
        hstack->GetXaxis()->SetRangeUser(
            (*hstd)()->GetXaxis()->GetXmin() / std::sqrt(MIN_NUM+in-1),
            (*hstd)()->GetXaxis()->GetXmax() / std::sqrt(MIN_NUM+in-1)
        );
        hstack->Draw("nostack hist");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table()->SetHeader(Form("N_{Hit} = %ld", MIN_NUM+in-1));
        leg_table()->AddEntry((*hstd)(), Form("Gaussian #sigma=%5.2f", (*hSMGLXstdS)()->GetBinContent(in)), "l");
        leg_table()->AddEntry((*hsmg)(), Form("SMG-Func #sigma=%5.2f", (*hSMGLXsmgS)()->GetBinContent(in)), "l");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
    }

    THStack* hsSMGLXM = Hist::Collect("hsSMGLXmen", HistList({ hSMGLXstdM, hSMGLXsmgM }));
    editor.create();
    hsSMGLXM->Draw("nostack pe");
    hsSMGLXM->GetHistogram()->SetLineColor(0);
    hsSMGLXM->GetHistogram()->SetMarkerColor(0);
    hsSMGLXM->GetXaxis()->SetTitle("N_{Hit}");
    hsSMGLXM->GetYaxis()->SetTitle("Mean [10^{-4}]");
    hsSMGLXM->Draw("nostack pe");
    Legend leg_SMGLXM("", TextStyle(kBlack, 30, 43), PadWindow(0.60, 0.85, 0.15, 0.35));
    leg_SMGLXM()->SetHeader("Mean");
    leg_SMGLXM()->AddEntry((*hSMGLXstdM)(), "Gaussian", "l");
    leg_SMGLXM()->AddEntry((*hSMGLXsmgM)(), "SMG-Func", "l");
    leg_SMGLXM()->SetFillColor(0);
    leg_SMGLXM.draw();
    editor.save();

    THStack* hsSMGLXS = Hist::Collect("hsSMGLXmen", HistList({ hSMGLXstdS, hSMGLXsmgS }));
    editor.create();
    hsSMGLXS->Draw("nostack pe");
    hsSMGLXS->GetHistogram()->SetLineColor(0);
    hsSMGLXS->GetHistogram()->SetMarkerColor(0);
    hsSMGLXS->GetXaxis()->SetTitle("N_{Hit}");
    hsSMGLXS->GetYaxis()->SetTitle("Resolution [10^{-4}]");
    hsSMGLXS->Draw("nostack pe");
    Legend leg_SMGLXS("", TextStyle(kBlack, 30, 43), PadWindow(0.60, 0.85, 0.65, 0.85));
    leg_SMGLXS()->SetHeader("Resolution");
    leg_SMGLXS()->AddEntry((*hSMGLXstdM)(), "Gaussian", "l");
    leg_SMGLXS()->AddEntry((*hSMGLXsmgM)(), "SMG-Func", "l");
    leg_SMGLXS()->SetFillColor(0);
    leg_SMGLXS.draw();
    editor.save();

    editor.create();
    (*hSMGLXrsoR)()->GetYaxis()->SetTitle("#sigma_{SMG-Func} / #sigma_{Gaussian}");
    (*hSMGLXrsoR)()->Draw("pe");
    editor.save();

    
    
    
    
    
    
    
    
    
    
    std::cerr << Form("==== ALG : EH ====\n");
    const std::array<double, 4> ALGELpars({ 
        4.07524e-01, 2.14736e+00, 2.40172e-01, 1.18389e-01
    });
    const std::array<double, 4> ALGEMpars({ 
        2.19427e-02, 1.39837e+00, 1.01428e-01, 1.13978e-01
    });
    const std::array<double, 4> ALGEHpars({ 
        4.82158e-03, 1.12452e+00, 8.38657e-02, 9.36490e-02  
    });
    Hist* hALGEHstdM = Hist::New("hALGEHstdM", HistAxis(AXNUM));
    Hist* hALGEHstdS = Hist::New("hALGEHstdS", HistAxis(AXNUM));
    Hist* hALGEHalgM = Hist::New("hALGEHalgM", HistAxis(AXNUM));
    Hist* hALGEHalgS = Hist::New("hALGEHalgS", HistAxis(AXNUM));
    Hist* hALGEHrsoR = Hist::New("hALGEHrsoR", HistAxis(AXNUM));
        
    hALGEHstdM->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kSquare)));
    hALGEHstdS->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kSquare)));
	hALGEHalgM->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle)));
	hALGEHalgS->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle)));
	hALGEHrsoR->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle)));
    
    //std::vector<Hist*> vhALGEHstd = Hist::ProjectAll(HistProj::kY, Hist::Head("hALGELstd"));
    //std::vector<Hist*> vhALGEHalg = Hist::ProjectAll(HistProj::kY, Hist::Head("hALGELalg"));
    //std::vector<Hist*> vhALGEHstd = Hist::ProjectAll(HistProj::kY, Hist::Head("hALGEMstd"));
    //std::vector<Hist*> vhALGEHalg = Hist::ProjectAll(HistProj::kY, Hist::Head("hALGEMalg"));
    std::vector<Hist*> vhALGEHstd = Hist::ProjectAll(HistProj::kY, Hist::Head("hALGEHstd"));
    std::vector<Hist*> vhALGEHalg = Hist::ProjectAll(HistProj::kY, Hist::Head("hALGEHalg"));


    ////double ALGpeak = ALGELpars[1] + TrSys::ConvGaus::GetLGS(ALGELpars[0], ALGELpars[3]/ALGELpars[2]) * ALGELpars[2];
    ////double ALGpeak = ALGEMpars[1] + TrSys::ConvGaus::GetLGS(ALGEMpars[0], ALGEMpars[3]/ALGEMpars[2]) * ALGEMpars[2];
    //double ALGpeak = ALGEHpars[1] + TrSys::ConvGaus::GetLGS(ALGEHpars[0], ALGEHpars[3]/ALGEHpars[2]) * ALGEHpars[2];
    //double ALGpeak = ALGELpars[1];
    //double ALGpeak = ALGEMpars[1];
    double ALGpeak = ALGEHpars[1];
    std::cout << Form("\n==== ALG peak %14.8f ====\n", ALGpeak);

    for (int in = 1; in <= AXNUM.nbin(); ++in) {
        Hist* hstd = vhALGEHstd.at(in);
        Hist* halg = vhALGEHalg.at(in);
		
        hstd->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kSquare)));
		halg->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle)));

        gsfunc->SetParameters(1.0, 0.0, 0.0002);
        (*hstd)()->Fit(gsfunc, "q0", "");
        (*hstd)()->Fit(gsfunc, "q0", "");
        (*hALGEHstdM)()->SetBinContent(in, gsfunc->GetParameter(1) - ALGpeak);
        (*hALGEHstdM)()->SetBinError  (in, gsfunc->GetParError (1));
        (*hALGEHstdS)()->SetBinContent(in, gsfunc->GetParameter(2));
        (*hALGEHstdS)()->SetBinError  (in, gsfunc->GetParError (2));
        
        gsfunc->SetParameters(1.0, 0.0, 0.0002);
        (*halg)()->Fit(gsfunc, "q0", "");
        (*halg)()->Fit(gsfunc, "q0", "");
        (*hALGEHalgM)()->SetBinContent(in, gsfunc->GetParameter(1) - ALGpeak);
        (*hALGEHalgM)()->SetBinError  (in, gsfunc->GetParError (1));
        (*hALGEHalgS)()->SetBinContent(in, gsfunc->GetParameter(2));
        (*hALGEHalgS)()->SetBinError  (in, gsfunc->GetParError (2));
       
        (*hALGEHrsoR)()->SetBinContent(in, (*hALGEHalgS)()->GetBinContent(in) / (*hALGEHstdS)()->GetBinContent(in));

        editor.create();
        THStack* hstack = Hist::Collect(Form("hsALGEH%02ld\n", MIN_NUM+in-1), HistList({hstd, halg}));
        hstack->Draw("nostack hist");
        hstack->GetHistogram()->SetLineColor(0);
        hstack->GetHistogram()->SetMarkerColor(0);
        hstack->GetXaxis()->SetTitle("Value");
        hstack->GetYaxis()->SetTitle("Events/Bin");
        hstack->GetXaxis()->SetRangeUser(
            ALGpeak + ((*hstd)()->GetXaxis()->GetXmin() - ALGpeak) / std::sqrt(MIN_NUM+in-1),
            ALGpeak + ((*hstd)()->GetXaxis()->GetXmax() - ALGpeak) / std::sqrt(MIN_NUM+in-1)
        );
        hstack->Draw("nostack hist");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.58, 0.85, 0.65, 0.85));
        leg_table()->SetHeader(Form("N_{Hit} = %ld", MIN_NUM+in-1));
        leg_table()->AddEntry((*hstd)(), Form("Gaussian #sigma=%f", (*hALGEHstdS)()->GetBinContent(in)), "l");
        leg_table()->AddEntry((*halg)(), Form("ALG-Func #sigma=%f", (*hALGEHalgS)()->GetBinContent(in)), "l");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
    }

    THStack* hsALGEHM = Hist::Collect("hsALGEHmen", HistList({ hALGEHstdM, hALGEHalgM }));
    editor.create();
    hsALGEHM->Draw("nostack pe");
    hsALGEHM->GetHistogram()->SetLineColor(0);
    hsALGEHM->GetHistogram()->SetMarkerColor(0);
    hsALGEHM->GetXaxis()->SetTitle("N_{Hit}");
    hsALGEHM->GetYaxis()->SetTitle("Mean - #DeltaE_{MPV}");
    hsALGEHM->Draw("nostack pe");
    Legend leg_ALGEHM("", TextStyle(kBlack, 30, 43), PadWindow(0.60, 0.85, 0.40, 0.60));
    leg_ALGEHM()->SetHeader(Form("Mean - #DeltaE_{MPV}(%4.2f)", ALGpeak));
    leg_ALGEHM()->AddEntry((*hALGEHstdM)(), "Gaussian", "l");
    leg_ALGEHM()->AddEntry((*hALGEHalgM)(), "ALG-Func", "l");
    leg_ALGEHM()->SetFillColor(0);
    leg_ALGEHM.draw();
    editor.save();

    THStack* hsALGEHS = Hist::Collect("hsALGEHmen", HistList({ hALGEHstdS, hALGEHalgS }));
    editor.create();
    hsALGEHS->Draw("nostack pe");
    hsALGEHS->GetHistogram()->SetLineColor(0);
    hsALGEHS->GetHistogram()->SetMarkerColor(0);
    hsALGEHS->GetXaxis()->SetTitle("N_{Hit}");
    hsALGEHS->GetYaxis()->SetTitle("Resolution");
    hsALGEHS->Draw("nostack pe");
    Legend leg_ALGEHS("", TextStyle(kBlack, 30, 43), PadWindow(0.60, 0.85, 0.65, 0.85));
    leg_ALGEHS()->SetHeader("Resolution");
    leg_ALGEHS()->AddEntry((*hALGEHstdM)(), "Gaussian", "l");
    leg_ALGEHS()->AddEntry((*hALGEHalgM)(), "ALG-Func", "l");
    leg_ALGEHS()->SetFillColor(0);
    leg_ALGEHS.draw();
    editor.save();

    editor.create();
    (*hALGEHrsoR)()->GetYaxis()->SetTitle("#sigma_{ALG-Func} / #sigma_{Gaussian}");
    (*hALGEHrsoR)()->Draw("pe");
    editor.save();

    editor.close();

    TFile * ofle = new TFile("out/nogaus_fit.root", "RECREATE");
    ofle->cd();
    
    (*hSMGLXstdM)()->Write();
    (*hSMGLXstdS)()->Write();
    (*hSMGLXsmgM)()->Write();
    (*hSMGLXsmgS)()->Write();
    (*hSMGLXrsoR)()->Write();
    
    (*hALGEHstdM)()->Write();
    (*hALGEHstdS)()->Write();
    (*hALGEHalgM)()->Write();
    (*hALGEHalgS)()->Write();
    (*hALGEHrsoR)()->Write();

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
