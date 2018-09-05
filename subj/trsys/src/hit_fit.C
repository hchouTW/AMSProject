#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
    
static TF1* flg = new TF1("flg", "[0] * TMath::Exp([1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) + (1-[1]) * TMath::Log(TMath::Landau(1.17741002*(x-[2])/[3]-2.22782980e-01)/1.80655634e-01))");

Double_t flgcov(Double_t* x, Double_t* par) {
    flg->SetParameters(1.0, par[1], par[2], par[3]);

    // control constants
    Int_t    np = 50;
    Double_t sc = 2.5;

    Double_t xv  = x[0];
    Double_t xlw = xv - sc * par[4];
    Double_t xup = xv + sc * par[4];
    Double_t stp = (xup - xlw) / static_cast<Double_t>(np);
    
    Double_t sum = 0.;
    Double_t wgt = 0.;
    for (Int_t it = 0; it <= np; ++it) {
        Double_t xx = xlw + stp * static_cast<Double_t>(it);
        Double_t gs = TMath::Gaus(xv, xx, par[4]);
        sum += gs * flg->Eval(xx);
        wgt += gs;
    }

    Double_t val = par[0] * (sum / wgt);
    return val;
}

/*
TF1* flggm = new TF1("flggm", "[0] * (1.0/sqrt(2.0*TMath::Pi())/[3]) * TMath::Exp((1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3])) + [4] * (pow([6], [5]) / TMath::Gamma([5])) * TMath::Power(x,[5]-1) * TMath::Exp(-[6]*x) * 0.5 * (TMath::Erf((x-[7])/[8])+1)", 0.1, 100);

TF1* flg = new TF1("flg", "[0] * (1.0/sqrt(2.0*TMath::Pi())/[3]) * TMath::Exp((1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]))");
TF1* fgm = new TF1("fgm", "[0] * (pow([2], [1]) / TMath::Gamma([1])) * TMath::Power(x,[1]-1) * TMath::Exp(-[2]*x) * 0.5 * (TMath::Erf((x-[3])/[4])+1)");
*/

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
 
    int ibin = std::atoi(argv[1]);

    Hist::Load("hit_fill.root", "/ams_home/hchou/AMSProject/subj/trsys/dat");

    // Fit
    //Hist* hAdc = Hist::Head("hTKadcx");
    //Hist* hAdc = Hist::Head("hTKadcy");
    Hist* hAdc = Hist::Head("hTFadc");
    //Hist* hAdc = Hist::Head("hTDavg");
    //Hist* hAdc = Hist::Head("hTDex2");
    std::vector<Hist*>&& vhAdc = Hist::ProjectAll(HistProj::kY, hAdc);

    const Axis& AXeta = hAdc->xaxis();
    const Axis& AXadc = hAdc->yaxis();
    
    TFile * ofle = new TFile(Form("/ams_home/hchou/AMSProject/subj/trsys/dat/hit_fit%04d.root", ibin), "RECREATE");
    ofle->cd();

    //Double_t pkpa = 0.001;
    //std::vector<Double_t> pmpv({ 1.14489e+01, 1.39737e+00, 7.35383e-03, 1.23680e+00, 2.73847e+00, 7.93824e-02 });
    //std::vector<Double_t> psgm({ 5.17642e-02, 1.11762e+01, 6.77707e+00, 6.52921e-01, 1.48594e-03, 1.24828e+00 });
    
    //std::vector<Double_t> pkpa({ 2.17486e-04, 4.27812e+00 });
    //std::vector<Double_t> pmpv({ 9.62033e+00, 1.26248e+00, 8.36500e-03, 1.24422e+00, 2.22255e+00, 8.06270e-02 });
    //std::vector<Double_t> psgm({ 2.13993e-01, 2.78925e+00, 1.65428e+00, 6.55129e-01, 1.17441e-01, 3.61303e-01 });

    //TF1* fkpa = new TF1("fkpa", "[0] + (1-[0]) * 0.5 * (1 + TMath::Erf(x-[1]))");
    //fkpa->SetParameters(2.17486e-04, 4.27812e+00);
    //
    //TF1* fmpv = new TF1("fmpv", "[0] * (1+x*x)^[3] * ([1] - [2]*(1+x*x)^(-[3]) - TMath::Log([4]+x^[5]))");
    //fmpv->SetParameters(9.62033e+00, 1.26248e+00, 8.36500e-03, 1.24422e+00, 2.22255e+00, 8.06270e-02);
    //
    //TF1* fsgm = new TF1("fsgm", "[0] * (1+x*x)^[3] * ([1] - [2]*(1+x*x)^(-[3]) - TMath::Log([4]+x^[5]))");
    //fsgm->SetParameters(2.13993e-01, 2.78925e+00, 1.65428e+00, 6.55129e-01, 1.17441e-01, 3.61303e-01);
    
    //TF1* fkpa = new TF1("fkpa", "(1 - [0]) + [0] * 0.5 * (1.0 + TMath::Erf([1] * TMath::Log(1+[3]*(x*x)) + [2]))");
    //fkpa->SetParameters(1.0, 1.0, 0.3, 1.0);
    TF1* fkpa = new TF1("fkpa", "0.5 * (1.0 + TMath::Erf([0] * TMath::Log(1+[1]*(1+x*x)^[2]) - [3]))");
    //fkpa->SetParameters(3.63156e+04, 7.49003e-07, 1.67437e+00, 1.99507e+00); // TKadcx
    //fkpa->SetParameters(1.81137e+04, 9.99578e-07, 1.92633e+00, 1.98827e+00); // TKadcy
    fkpa->SetParameters(7.23786e-02, 8.94961e+01, 2.54620e+01, 2.42002e+00); // TFadc

    Hist* hAdcK = Hist::New("hAdcK", HistAxis(AXeta, "Kappa"));
    Hist* hAdcM = Hist::New("hAdcM", HistAxis(AXeta, "Mpv"));
    Hist* hAdcS = Hist::New("hAdcS", HistAxis(AXeta, "Sigma"));
    Hist* hAdcF = Hist::New("hAdcF", HistAxis(AXeta, "Fluc"));

    TF1* func = new TF1("func", flgcov, 0, 10, 5);
    func->SetNpx(1000);
    //for (int it = 30; it <= AXeta.nbin(); ++it) {
    for (int it = ibin; it <= ibin; ++it) {
        Double_t eta = AXeta.center(it, AxisScale::kLog);
        COUT("Process ITER %d\n", it);
        Double_t kpa = 0.01;
        Double_t mpv = (*vhAdc.at(it))()->GetBinCenter((*vhAdc.at(it))()->GetMaximumBin());
        Double_t rms = (*vhAdc.at(it))()->GetRMS();
        func->SetParameters(1000, kpa, mpv, 0.3443236, 0.0772006);
        func->SetParLimits(1, 0.0, 1.0);
        func->SetParLimits(2, 0.0, 5.0*mpv);
        func->SetParLimits(3, 0.0, 20.0*rms);
        func->SetParLimits(4, 0.0, 20.0*rms);
        func->FixParameter(1, fkpa->Eval(eta));
        //func->FixParameter(4, 0.284416); // TKadcx
        //func->FixParameter(4, 0.175361); // TKadcy
        func->FixParameter(4, 0.0772006); // TFadc
       
        (*vhAdc.at(it))()->Fit(func, "q0", "", mpv-5*rms, mpv+10*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));

        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhAdc.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+15*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));

        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhAdc.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+15*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhAdc.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+15*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhAdc.at(it))()->Fit(func, "q0", "", mpv-4*rms, mpv+4*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
    
    
        (*hAdcK)()->SetBinContent(it, func->GetParameter(1));
        (*hAdcK)()->SetBinError  (it, func->GetParError(1));
        (*hAdcM)()->SetBinContent(it, func->GetParameter(2));
        (*hAdcM)()->SetBinError  (it, func->GetParError(2));
        (*hAdcS)()->SetBinContent(it, func->GetParameter(3));
        (*hAdcS)()->SetBinError  (it, func->GetParError(3));
        (*hAdcF)()->SetBinContent(it, func->GetParameter(4));
        (*hAdcF)()->SetBinError  (it, func->GetParError(4));

        (*vhAdc.at(it))()->Write();
        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXadc));
        for (int jt = 1; jt <= AXadc.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXadc.center(jt)));
        }

        (vhAdc.at(it))->style(Fill(), Line(kBlue), Marker(kBlue));
        tmpl->style(Fill(), Line(kRed), Marker(kRed));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhAdc.at(it), tmpl }));
        ch->Write();
    }
 
    hAdcK->Write();
    hAdcM->Write();
    hAdcS->Write();
    hAdcF->Write();

/*    
    std::vector<Double_t> pbta({ 1.83481e-01, 5.07634e-01, -5.30224e+00 });
    TF1* fbta = new TF1("fbta", "[0] + (1-[0]) * 0.5 * (1 + TMath::Erf([1]*(TMath::Log(x)-[2])))");
    fbta->SetParameters(1.83481e-01, 5.07634e-01, -5.30224e+00);

    std::vector<Double_t> pwgt({ 6.87250e-01, 8.49158e-01, -6.50955e+00 });
    TF1* fwgt = new TF1("fwgt", "[0] * 0.5 * TMath::Erfc([1]*(TMath::Log(x)-[2]))");
    fwgt->SetParameters(6.87250e-01, 8.49158e-01, -6.50955e+00 );

    Hist* hAdcW  = Hist::New("hAdcW",  HistAxis(AXeta, "Wgt"));
    Hist* hAdcW2  = Hist::New("hAdcW2",  HistAxis(AXeta, "Wgt"));

    Hist* hAdcA = Hist::New("hAdcA", HistAxis(AXeta, "Alpha"));
    Hist* hAdcB = Hist::New("hAdcB", HistAxis(AXeta, "Beta"));
    Hist* hAdcM = Hist::New("hAdcM", HistAxis(AXeta, "ErfM"));
    Hist* hAdcS = Hist::New("hAdcS", HistAxis(AXeta, "ErfS"));

    TF1* func = flggm;
    func->SetNpx(100000);
    for (int it = 1; it <= AXeta.nbin() && it <= 143; ++it) {
        Double_t eta = AXeta.center(it, AxisScale::kLog);
        COUT("Process ITER %d\n", it);
        
        for (Int_t i = 0; i <= 8; ++i) func->ReleaseParameter(i);
        flggm->SetParameters(1.84683e+04, 2.53358e-04, 2.31822e+00, 6.54225e-01, 8.39894e+03, 2.18565e+00, 2.14286e-01, 6.60662e+00, 1.77739e+00);
        
        func->FixParameter(1, fkpa->Eval(eta));
        func->FixParameter(2, fmpv->Eval(eta));
        func->FixParameter(3, fsgm->Eval(eta));
        
        //func->SetParLimits(5, 0.0, 20.0);
        //func->SetParLimits(6, 0.0, 20.0);
        //func->SetParLimits(7, 0.0, 10.0);
        //func->SetParLimits(8, 0.0,  5.0);
        
        func->FixParameter(5, 2.11427);
        func->FixParameter(6, fbta->Eval(eta));
        
        func->FixParameter(7, 6.52658434);
        func->FixParameter(8, 1.79160092);
        
        func->SetParLimits(4, 0.0, 109290.0*1000.0);
        func->SetParameter(0, 1.0e+4);
        func->SetParameter(4, 1.0e+3);

        flg->SetParameters(func->GetParameter(0), func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));
        fgm->SetParameters(func->GetParameter(4), func->GetParameter(5), func->GetParameter(6), func->GetParameter(7), func->GetParameter(8));
        
        (*vhAdc.at(it))()->Fit(func, "q0", "");
        (*vhAdc.at(it))()->Fit(func, "q0", "", 0.5, 30.);
        
        CERR("FIT == LG KPA %14.8f MPV %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));
        CERR("FIT == GM ALP %14.8f BTA %14.8f M %14.8f S %14.8f\n", func->GetParameter(5), func->GetParameter(6), func->GetParameter(7), func->GetParameter(8));
        
        flg->SetParameters(func->GetParameter(0), func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));
        fgm->SetParameters(func->GetParameter(4), func->GetParameter(5), func->GetParameter(6), func->GetParameter(7), func->GetParameter(8));
    
        (*hAdcA)()->SetBinContent(it, func->GetParameter(5));
        (*hAdcA)()->SetBinError  (it, func->GetParError(5));
        (*hAdcB)()->SetBinContent(it, func->GetParameter(6));
        (*hAdcB)()->SetBinError  (it, func->GetParError(6));
        (*hAdcM)()->SetBinContent(it, func->GetParameter(7));
        (*hAdcM)()->SetBinError  (it, func->GetParError(7));
        (*hAdcS)()->SetBinContent(it, func->GetParameter(8));
        (*hAdcS)()->SetBinError  (it, func->GetParError(8));
       
        (*vhAdc.at(it))()->Write();
        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXadc));
        Hist* tmpl_lg = Hist::New(Form("tmpl%03d_lg", it), HistAxis(AXadc));
        Hist* tmpl_gm = Hist::New(Form("tmpl%03d_gm", it), HistAxis(AXadc));
        for (int jt = 1; jt <= AXadc.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXadc.center(jt)));
            (*tmpl_lg)()->SetBinContent(jt, flg->Eval(AXadc.center(jt)));
            (*tmpl_gm)()->SetBinContent(jt, fgm->Eval(AXadc.center(jt)));
        }

        (vhAdc.at(it))->style(Fill(), Line(kBlue), Marker(kBlue));
        tmpl->style(Fill(), Line(kRed), Marker(kRed));
        tmpl_lg->style(Fill(), Line(kYellow+2), Marker(kYellow+2));
        tmpl_gm->style(Fill(), Line(kGreen+2), Marker(kGreen+2));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhAdc.at(it), tmpl, tmpl_lg, tmpl_gm }));
        ch->Write();

        Double_t wgtL = flg->Integral(0.01, 1000.);
        Double_t wgtG = fgm->Integral(0.01, 1000.);

        Double_t wgt = (wgtG / (wgtL + wgtG));
        (*hAdcW)()->SetBinContent(it, wgt);
       
        Double_t factL = func->GetParameter(0);
        Double_t factG = func->GetParameter(4);
        Double_t fact  = factG / (factL + factG);
        (*hAdcW2)()->SetBinContent(it, fact);

        CERR("FIT == LG %14.8f GM %14.8f WGT %14.8f Fact %14.8f\n", wgtL, wgtG, wgt, fact);
    }
 
    hAdcW->Write();

    hAdcA->Write();
    hAdcB->Write();
    hAdcS->Write();
    hAdcM->Write();
*/    
    ofle->Write();
    ofle->Close();

    return 1;
}
