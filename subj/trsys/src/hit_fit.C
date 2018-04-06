#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
    
static TF1* flg = new TF1("flg", "[0] * TMath::Exp( (1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) )");

Double_t flgcov(Double_t* x, Double_t* par) {
    flg->SetParameters(1.0, par[1], par[2], par[3]);

    // control constants
    Int_t    np = 50.0;
    Double_t sc =  3.5;

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


int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
  
    std::string iopath = "/ams_home/hchou/AMSProject/subj/trsys/dat";
    //std::string iopath = "/afs/cern.ch/work/h/hchou/AMSData/test18";
    Int_t idx = std::atoi(argv[1]);

    //Hist::Load("hit_fill.root", "dat");
    Hist::Load("hit_fill.root", iopath);

    // Fit
    Hist* hAdc = Hist::Head("hTKadcy");
    std::vector<Hist*>&& vhAdc = Hist::ProjectAll(HistProj::kY, hAdc);

    const Axis& AXeta = hAdc->xaxis();
    const Axis& AXadc = hAdc->yaxis();
    
    TFile * ofle = new TFile(Form("%s/hit_fit%04d.root", iopath.c_str(), idx), "RECREATE");
    ofle->cd();
    
    Hist* hAdcK = Hist::New("hAdcK", HistAxis(AXeta, "Kappa"));
    Hist* hAdcM = Hist::New("hAdcM", HistAxis(AXeta, "Mpv"));
    Hist* hAdcS = Hist::New("hAdcS", HistAxis(AXeta, "Sigma"));
    Hist* hAdcF = Hist::New("hAdcF", HistAxis(AXeta, "Fluc"));

    TF1* func = new TF1("func", flgcov, 0, 10, 5);
    //for (int it = 1; it <= AXeta.nbin(); ++it) {
    for (int it = idx; it <= idx && it <= AXeta.nbin(); ++it) {
        Double_t eta = AXeta.center(it, AxisScale::kLog);
        Double_t bta = 1.0/std::sqrt(1.0+eta*eta);
        COUT("Process ITER %d\n", it);
        Double_t kpa = 0.02;
        Double_t mpv = (*vhAdc.at(it))()->GetBinCenter((*vhAdc.at(it))()->GetMaximumBin());
        Double_t rms = (*vhAdc.at(it))()->GetRMS();
        func->SetParameters(1000, kpa, mpv, rms, 0.3*rms);
        func->SetParLimits(1, 0.0, 1.0);
        func->SetParLimits(2, 0.0, 5.0*mpv);
        func->SetParLimits(3, 0.0, 10.0*rms);
        func->SetParLimits(4, 0.0, 10.0*rms);
        //func->FixParameter(4, (bta*bta)*0.0829427); // TFadc
        //func->FixParameter(1, (1.0-0.5*TMath::Erfc(5.07474e-01*TMath::Log(1.0+4.16139e+00*eta)-1.96169e+00))); // TKadcx
        //func->FixParameter(4, (bta*bta)*0.230493); // TKadcx
        //func->FixParameter(4, 0.166633); // TKadcy
        
        Double_t tmpsgm = rms;
        (*vhAdc.at(it))()->Fit(func, "q0", "", mpv-2*tmpsgm, mpv+2*tmpsgm);
        tmpsgm = std::sqrt(func->GetParameter(3) * func->GetParameter(3) + func->GetParameter(4) * func->GetParameter(4));
        //(*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(2)-3*tmpsgm, func->GetParameter(2)+4.*tmpsgm); // TOF
        (*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(2)-3*tmpsgm, func->GetParameter(2)+5.*tmpsgm); // TRK
        tmpsgm = std::sqrt(func->GetParameter(3) * func->GetParameter(3) + func->GetParameter(4) * func->GetParameter(4));
        //(*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(2)-3*tmpsgm, func->GetParameter(2)+3.*tmpsgm); // TOF
        (*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(2)-3*tmpsgm, func->GetParameter(2)+5.*tmpsgm); // TRK
    
        (*hAdcK)()->SetBinContent(it, func->GetParameter(1));
        (*hAdcK)()->SetBinError  (it, func->GetParError(1));
        (*hAdcM)()->SetBinContent(it, (1.0/bta/bta)*func->GetParameter(2));
        (*hAdcM)()->SetBinError  (it, (1.0/bta/bta)*func->GetParError(2));
        (*hAdcS)()->SetBinContent(it, (1.0/bta/bta)*func->GetParameter(3));
        (*hAdcS)()->SetBinError  (it, (1.0/bta/bta)*func->GetParError(3));
        (*hAdcF)()->SetBinContent(it, (1.0/bta/bta)*func->GetParameter(4));
        (*hAdcF)()->SetBinError  (it, (1.0/bta/bta)*func->GetParError(4));

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

    ofle->Write();
    ofle->Close();

    return 1;
}
