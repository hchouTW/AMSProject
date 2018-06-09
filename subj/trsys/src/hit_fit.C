#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
    

static TF1* flg = new TF1("flg", "[0] * TMath::Exp( (1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) )");

Double_t flgcov(Double_t* x, Double_t* par) {
    flg->SetParameters(1.0, par[1], par[2], par[3]);

    // control constants
    Int_t    np = 40.0;
    Double_t sc =  3.0;

    Double_t xv  = x[0];
    Double_t xlw = xv - sc * par[4];
    Double_t xup = xv + sc * par[4];
    Double_t stp = (xup - xlw) / static_cast<Double_t>(np);
    
    Double_t sum = 0.;
    Double_t wgt = 0.;
    for (Int_t it = 0; it <= np; ++it) {
        Double_t xx = xlw + stp * static_cast<Double_t>(it);
        Double_t gs = TMath::Gaus(xv, xx, par[4]);
        sum += gs * ((xx <= 0) ? 0 : flg->Eval(xx));
        wgt += gs;
    }

    Double_t val = par[0] * (sum / wgt);
    return val;
}

/*
TF1* flggm = new TF1("flggm", "[0] * TMath::Exp((1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3])) * (TMath::Erf((x-[4])/[5])+1) + [6] * TMath::Power(x,[7]) * TMath::Exp(-[8]*x) * (TMath::Erf((x-[9])/[10])+1)");
//flggm->SetParameters(1.86244e+04, 8.01663e-03, 2.46953e+00, 6.79641e-01, 8.69772e-01, 2.90538e-01, 1.10869e+02, 3.18041e+00, 4.10684e-01, 6.50848e+00, 1.86977e+00);
//flggm->SetParLimits(1, 0.0, 1.0);
//flggm->SetNpx(10000);

TF1* flg = new TF1("flg", "[0] * TMath::Exp((1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3])) * (TMath::Erf((x-[4])/[5])+1)");
TF1* fgm = new TF1("fgm", "[0] * TMath::Power(x,[1]) * TMath::Exp(-[2]*x) * (TMath::Erf((x-[3])/[4])+1)");
*/
int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
  
    Hist::Load("hit_fill.root", "dat");

    // Fit
    //Hist* hAdc = Hist::Head("hTKadcx");
    //Hist* hAdc = Hist::Head("hTKadcy");
    Hist* hAdc = Hist::Head("hTFadc");
    //Hist* hAdc = Hist::Head("hTDavg");
    //Hist* hAdc = Hist::Head("hTDex");
    std::vector<Hist*>&& vhAdc = Hist::ProjectAll(HistProj::kY, hAdc);

    const Axis& AXeta = hAdc->xaxis();
    const Axis& AXadc = hAdc->yaxis();
    
    TFile * ofle = new TFile("hit_fit.root", "RECREATE");
    ofle->cd();
    
    Hist* hAdcK = Hist::New("hAdcK", HistAxis(AXeta, "Kappa"));
    Hist* hAdcM = Hist::New("hAdcM", HistAxis(AXeta, "Mpv"));
    Hist* hAdcS = Hist::New("hAdcS", HistAxis(AXeta, "Sigma"));
    Hist* hAdcF = Hist::New("hAdcF", HistAxis(AXeta, "Fluc"));

    TF1* func = new TF1("func", flgcov, 0, 10, 5);
    //for (int it = 1; it <= AXeta.nbin(); ++it) {
    for (int it = 40; it <= AXeta.nbin(); ++it) {
        Double_t eta = AXeta.center(it, AxisScale::kLog);
        COUT("Process ITER %d\n", it);
        Double_t kpa = 0.001;
        Double_t mpv = (*vhAdc.at(it))()->GetBinCenter((*vhAdc.at(it))()->GetMaximumBin());
        Double_t rms = (*vhAdc.at(it))()->GetRMS();
        //func->SetParameters(1000, kpa, mpv, rms, 0.3*rms);
        func->SetParameters(1000, kpa, mpv, 0.25, 0.01);
        func->SetParLimits(1, 0.0, 1.0);
        func->SetParLimits(2, 0.0, 5.0*mpv);
        func->SetParLimits(3, 0.0, 10.0*rms);
        func->SetParLimits(4, 0.0, 10.0*rms);
        func->FixParameter(4, 0.0824851); // TFadc
        //func->FixParameter(4, 0.266246); // TKadcx
        //func->FixParameter(4, 0.152518); // TKadcy
        //func->FixParameter(4, 0.17); // TDavg
        
        Double_t tmpsgm = rms;
        (*vhAdc.at(it))()->Fit(func, "q0", "", mpv-5*tmpsgm, mpv+8*tmpsgm);
        tmpsgm = std::sqrt(func->GetParameter(3) * func->GetParameter(3) + func->GetParameter(4) * func->GetParameter(4));

        (*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(2)-5*tmpsgm, func->GetParameter(2)+30.*tmpsgm); // TRK
        
        //tmpsgm = std::sqrt(func->GetParameter(3) * func->GetParameter(3) + func->GetParameter(4) * func->GetParameter(4));

        //(*vhAdc.at(it))()->Fit(func, "q0", ""); // TRK
        //(*vhAdc.at(it))()->Fit(func, "q0", ""); // TRK
        
        //(*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(2)-3*tmpsgm, func->GetParameter(2)+4.*tmpsgm); // TOF
        //(*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(2)-3*tmpsgm, func->GetParameter(2)+3.*tmpsgm); // TOF
        
    
        (*hAdcK)()->SetBinContent(it, func->GetParameter(1));
        (*hAdcK)()->SetBinError  (it, func->GetParError(1));
        (*hAdcM)()->SetBinContent(it, func->GetParameter(2));
        (*hAdcM)()->SetBinError  (it, func->GetParError(2));
        (*hAdcS)()->SetBinContent(it, func->GetParameter(3));
        (*hAdcS)()->SetBinError  (it, func->GetParError(3));
        (*hAdcF)()->SetBinContent(it, func->GetParameter(4));
        (*hAdcF)()->SetBinError  (it, func->GetParError(4));

        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));

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
    Hist* hAdcW  = Hist::New("hAdcW",  HistAxis(AXeta, "Wgt"));

    Hist* hAdcI_K  = Hist::New("hAdcI_K",  HistAxis(AXeta, "Kappa"));
    Hist* hAdcI_M  = Hist::New("hAdcI_M",  HistAxis(AXeta, "Mpv"));
    Hist* hAdcI_S  = Hist::New("hAdcI_S",  HistAxis(AXeta, "Sigma"));
    Hist* hAdcI_EM = Hist::New("hAdcI_EM", HistAxis(AXeta, "ErfM"));
    Hist* hAdcI_ES = Hist::New("hAdcI_ES", HistAxis(AXeta, "ErfS"));
    
    Hist* hAdcG_A  = Hist::New("hAdcG_A",  HistAxis(AXeta, "Alpha"));
    Hist* hAdcG_B  = Hist::New("hAdcG_B",  HistAxis(AXeta, "Beta"));
    Hist* hAdcG_EM = Hist::New("hAdcG_EM", HistAxis(AXeta, "ErfM"));
    Hist* hAdcG_ES = Hist::New("hAdcG_ES", HistAxis(AXeta, "ErfS"));

    flggm->SetNpx(10000);
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        Double_t eta = AXeta.center(it, AxisScale::kLog);
        COUT("Process ITER %d\n", it);
        Double_t kpa = 0.001;
        Double_t mpv = (*vhAdc.at(it))()->GetBinCenter((*vhAdc.at(it))()->GetMaximumBin());
        Double_t rms = (*vhAdc.at(it))()->GetRMS();
        
        for (Int_t i = 0; i <= 10; ++i) flggm->ReleaseParameter(i);
        flggm->SetParameters(1.64104e+04, 4.00868e-03, 2.34098e+00, 6.90771e-01, 5.68705e-01, 4.70127e-01, 1.05594e+02, 3.01722e+00, 3.55486e-01, 6.19620e+00, 1.50715e+00);
        
        flggm->SetParLimits(1,  0.0, 1.0);
        flggm->SetParLimits(2,  0.0, 30.0);
        flggm->SetParLimits(3,  0.0, 10.0);
        flggm->SetParLimits(4,  0.0, 5.0);
        flggm->SetParLimits(5,  0.0, 5.0);
        
        flggm->SetParLimits(7,  0.01, 30.0);
        flggm->SetParLimits(8,  0.01, 30.0);
        flggm->SetParLimits(9,  0.0, 10.0);
        flggm->SetParLimits(10, 0.0, 10.0);
    
        flggm->FixParameter(4,  5.68705e-01);
        flggm->FixParameter(5,  4.70127e-01);
        flggm->FixParameter(9,  6.19620e+00);
        flggm->FixParameter(10, 1.50715e+00);
        
        flg->SetParameters(flggm->GetParameter(0), flggm->GetParameter(1), flggm->GetParameter(2), flggm->GetParameter(3), flggm->GetParameter(4), flggm->GetParameter(5));
        fgm->SetParameters(flggm->GetParameter(6), flggm->GetParameter(7), flggm->GetParameter(8), flggm->GetParameter(9), flggm->GetParameter(10));
        
        Double_t wgtLv = flg->Integral(0.01, 100.);
        Double_t wgtGv = fgm->Integral(0.01, 100.);
        Double_t wgtv  = (wgtGv / (wgtLv + wgtGv));
        Double_t pred  = 1.99891e-01 * TMath::Erfc((std::log(eta) + 7.10465e+00) / 1.10334e+00);
        flggm->SetParameter(6, (pred / wgtv) * flggm->GetParameter(6));

        (*vhAdc.at(it))()->Fit(flggm, "q0", "");
        flggm->SetParameter(0, std::fabs(flggm->GetParameter(0)));
        flggm->SetParameter(6, std::fabs(flggm->GetParameter(6)));

        (*vhAdc.at(it))()->Fit(flggm, "q0", "");
        flggm->SetParameter(0, std::fabs(flggm->GetParameter(0)));
        flggm->SetParameter(6, std::fabs(flggm->GetParameter(6)));
        
        (*vhAdc.at(it))()->Fit(flggm, "q0", "", 0.5, 40.);
        flggm->SetParameter(0, std::fabs(flggm->GetParameter(0)));
        flggm->SetParameter(6, std::fabs(flggm->GetParameter(6)));
        
        (*vhAdc.at(it))()->Fit(flggm, "q0", "", 0.5, 40.);
        
        CERR("FIT == LG KPA %14.8f MPV %14.8f SMG %14.8f  EM %14.8f ES %14.8f\n", flggm->GetParameter(1), flggm->GetParameter(2), flggm->GetParameter(3), flggm->GetParameter(4), flggm->GetParameter(5));
        CERR("FIT == GM ALP %14.8f BTA %14.8f EM %14.8f ES %14.8f\n", flggm->GetParameter(7), flggm->GetParameter(8), flggm->GetParameter(9), flggm->GetParameter(10));
    
        (*hAdcI_K)()->SetBinContent(it, flggm->GetParameter(1));
        (*hAdcI_K)()->SetBinError  (it, flggm->GetParError(1));
        (*hAdcI_M)()->SetBinContent(it, flggm->GetParameter(2));
        (*hAdcI_M)()->SetBinError  (it, flggm->GetParError(2));
        (*hAdcI_S)()->SetBinContent(it, flggm->GetParameter(3));
        (*hAdcI_S)()->SetBinError  (it, flggm->GetParError(3));
        (*hAdcI_EM)()->SetBinContent(it, flggm->GetParameter(4));
        (*hAdcI_EM)()->SetBinError  (it, flggm->GetParError(4));
        (*hAdcI_ES)()->SetBinContent(it, flggm->GetParameter(5));
        (*hAdcI_ES)()->SetBinError  (it, flggm->GetParError(5));

        (*hAdcG_A)()->SetBinContent(it, flggm->GetParameter(7));
        (*hAdcG_A)()->SetBinError  (it, flggm->GetParError(7));
        (*hAdcG_B)()->SetBinContent(it, flggm->GetParameter(8));
        (*hAdcG_B)()->SetBinError  (it, flggm->GetParError(8));
        (*hAdcG_EM)()->SetBinContent(it, flggm->GetParameter(9));
        (*hAdcG_EM)()->SetBinError  (it, flggm->GetParError(9));
        (*hAdcG_ES)()->SetBinContent(it, flggm->GetParameter(10));
        (*hAdcG_ES)()->SetBinError  (it, flggm->GetParError(10));
       
        flg->SetParameters(flggm->GetParameter(0), flggm->GetParameter(1), flggm->GetParameter(2), flggm->GetParameter(3), flggm->GetParameter(4), flggm->GetParameter(5));
        fgm->SetParameters(flggm->GetParameter(6), flggm->GetParameter(7), flggm->GetParameter(8), flggm->GetParameter(9), flggm->GetParameter(10));

        (*vhAdc.at(it))()->Write();
        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXadc));
        Hist* tmpl_lg = Hist::New(Form("tmpl%03d_lg", it), HistAxis(AXadc));
        Hist* tmpl_gm = Hist::New(Form("tmpl%03d_gm", it), HistAxis(AXadc));
        for (int jt = 1; jt <= AXadc.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, flggm->Eval(AXadc.center(jt)));
            (*tmpl_lg)()->SetBinContent(jt, flg->Eval(AXadc.center(jt)));
            (*tmpl_gm)()->SetBinContent(jt, fgm->Eval(AXadc.center(jt)));
        }

        (vhAdc.at(it))->style(Fill(), Line(kBlue), Marker(kBlue));
        tmpl->style(Fill(), Line(kRed), Marker(kRed));
        tmpl_lg->style(Fill(), Line(kYellow+2), Marker(kYellow+2));
        tmpl_gm->style(Fill(), Line(kGreen+2), Marker(kGreen+2));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhAdc.at(it), tmpl, tmpl_lg, tmpl_gm }));
        ch->Write();

        Double_t wgtL = flg->Integral(0.01, 100.);
        Double_t wgtG = fgm->Integral(0.01, 100.);

        Double_t wgt = (wgtG / (wgtL + wgtG));
        (*hAdcW)()->SetBinContent(it, wgt);
        
        CERR("FIT == LG %14.8f GM %14.8f WGT %14.8f\n", wgtL, wgtG, wgt);
    }
 
    hAdcW->Write();

    hAdcI_K->Write();
    hAdcI_M->Write();
    hAdcI_S->Write();
    hAdcI_ES->Write();
    hAdcI_EM->Write();
    
    hAdcG_A->Write();
    hAdcG_B->Write();
    hAdcG_ES->Write();
    hAdcG_EM->Write();
*/
    ofle->Write();
    ofle->Close();

    return 1;
}
