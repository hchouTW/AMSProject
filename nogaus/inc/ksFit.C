#include <CPPLibs.h>
#include <ROOTLibs.h>
#include <TrSys.h>

#include "inc/SMGFit.h"
#include "inc/ALGFit.h"

double ALGfunc(double *x, double *p) {
    return TrSys::LandauGausFunc(x[0], p[0], p[1], p[2], p[3]);
}

int main() {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    
    Hist::Load("ksmap.root", "/afs/cern.ch/user/h/hchou/AMSProject/nogaus/out");
    std::vector<Hist*> vhksMap = Hist::ProjectAll(HistProj::kY, Hist::Head("hksMap"));
    Axis AXk = Hist::Head("hksMap")->xaxis();
    Axis AXs = Hist::Head("hksMap")->yaxis();
    AXk.set_title("#kappa");
    AXs.set_title("#sigma_{#nu}(#kappa)");

    Hist* hP0 = Hist::New("hP0", HistAxis(AXk));
    Hist* hP1 = Hist::New("hP1", HistAxis(AXk));
    Hist* hP2 = Hist::New("hP2", HistAxis(AXk));
    Hist* hP3 = Hist::New("hP3", HistAxis(AXk));

    hP0->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hP1->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hP2->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hP3->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));

    const std::array<double, 4> SMGpars({ 
        1.9, 0.25, 1.7, 1.0
    });
    //TF1* SMGfunc = new TF1("SMGfunc", "[0] * pow((x/[1]), [2]) / (1.0 + pow((x/[1]), [2]))", 0.0, 100.0);
    TF1* SMGfunc = new TF1("SMGfunc", "[0] * pow([1] * pow((x), [2]) / (1.0 + [1] * pow((x), [2])), [3])", 0.0, 100.0);
    for (unsigned long it = 0; it < SMGpars.size(); ++it) SMGfunc->SetParameter(it, SMGpars[it]);
    SMGfunc->SetLineColor(kRed);
    SMGfunc->SetNpx(10000);
 
    SMGfunc->SetParLimits(0, 0.0, 10.0);
    SMGfunc->SetParLimits(1, 0.0, 5.0);
    SMGfunc->SetParLimits(2, 0.0, 5.0);
    SMGfunc->FixParameter(3, 1.0);

    PdfEditor editor(Window(), "ksfit", "out");
   
    std::vector<double> vP0;
    std::vector<double> vP1;
    std::vector<double> vP2;
    std::vector<double> vP3;
    for (int ik = 1; ik <= AXk.nbin(); ++ik) {
        double ka = AXk.center(ik);

        Hist* hist = vhksMap.at(ik);
        hist->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kCircle)));
        (*hist)()->GetXaxis()->SetTitle("#sigma_{#nu}(#kappa)");
        //(*hist)()->GetYaxis()->SetTitle("#mathcal{z}_{0}(#kappa, #sigma_{#nu})");
        for (int is = 1; is <= AXs.nbin(); ++is) (*hist)()->SetBinError(is, 1.0e-6);
    
        //SMGfunc->FixParameter(1, 0.247339 + 3.57837 * k1 + (-7.07388) * k2 + 7.10421 * k3 + (-3.14195) * k4);
        
        (*hist)()->Fit(SMGfunc, "q0", "");
        (*hist)()->Fit(SMGfunc, "q0", "");
        (*hist)()->Fit(SMGfunc, "q0", "");
        (*hist)()->Fit(SMGfunc, "q0", "");
        (*hist)()->Fit(SMGfunc, "q0", "");

        (*hP0)()->SetBinContent(ik, SMGfunc->GetParameter(0));
        (*hP0)()->SetBinError  (ik, SMGfunc->GetParError (0));
        
        (*hP1)()->SetBinContent(ik, SMGfunc->GetParameter(1));
        (*hP1)()->SetBinError  (ik, SMGfunc->GetParError (1));
        
        (*hP2)()->SetBinContent(ik, SMGfunc->GetParameter(2));
        (*hP2)()->SetBinError  (ik, SMGfunc->GetParError (2));
        
        (*hP3)()->SetBinContent(ik, SMGfunc->GetParameter(3));
        (*hP3)()->SetBinError  (ik, SMGfunc->GetParError (3));
        
        vP0.push_back((*hP0)()->GetBinContent(ik));
        vP1.push_back((*hP1)()->GetBinContent(ik));
        vP2.push_back((*hP2)()->GetBinContent(ik));
        vP3.push_back((*hP3)()->GetBinContent(ik));

        std::cerr << Form("#kappa %6.4f  PARS %14.8f %14.8f %14.8f %14.8f\n", ka, vP0.back(), vP1.back(), vP2.back(), vP3.back());
        editor.create();
        editor.cd(1, PadAxis(1, 0));
        (*hist)()->Draw("pe");
        SMGfunc->Draw("l same");
        Legend leg_table("", TextStyle(kBlack, 40, 43), PadWindow(0.65, 0.85, 0.15, 0.35));
        leg_table()->SetHeader(Form("#kappa = %5.3f", ka));
        leg_table.draw();
        editor.save();
    }

    //editor.close();
   
    //int count = 1;
    //std::cerr << std::endl;
    //std::cerr << std::endl;
    //for (auto&& x : vP2) {
    //    std::cerr << Form("%7.5f, ", x);
    //    if (count%10 == 0) std::cerr << std::endl;
    //    count++;
    //}
    //std::cerr << std::endl;
    //std::cerr << std::endl;

    editor.create();
    editor.cd(1, PadAxis(0, 0));
    (*hP0)()->Draw("l");
    editor.save();
    
    editor.create();
    editor.cd(1, PadAxis(0, 0));
    (*hP1)()->Draw("l");
    editor.save();
    
    editor.create();
    editor.cd(1, PadAxis(0, 0));
    (*hP2)()->Draw("l");
    editor.save();
    
    editor.create();
    editor.cd(1, PadAxis(0, 0));
    (*hP3)()->Draw("l");
    editor.save();

    editor.close();

    TFile * ofle = new TFile("out/ksfit.root", "RECREATE");

    (*hP0)()->Write();
    (*hP1)()->Write();
    (*hP2)()->Write();
    (*hP3)()->Write();

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
