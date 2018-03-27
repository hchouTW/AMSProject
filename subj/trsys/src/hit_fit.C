#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
    
static TF1* flg = new TF1("flg", "[0] * TMath::Exp( (1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) )");

Double_t flgcov(Double_t* x, Double_t* par) {
    flg->SetParameters(1.0, par[1], par[2], par[3]);

    // control constants
    Int_t    np = 60.0;
    Double_t sc =  4.0;

    Double_t xlw = x[0] - sc * par[4];
    Double_t xup = x[0] + sc * par[4];
    Double_t stp = (xup - xlw) / static_cast<Double_t>(np);
    
    Double_t sum = 0.;
    Double_t wgt = 0.;
    for (Int_t it = 0; it <= np; ++it) {
        Double_t xx = xlw + stp * static_cast<Double_t>(it);
        Double_t gs = TMath::Gaus(x[0], xx, par[4]);
        sum += gs * flg->Eval(xx);
        wgt += gs;
    }

    Double_t val = par[0] * (sum / wgt);
    return val;
}


int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
  
    std::string iopath = "/ams_home/hchou/AMSProject/subj/trsys/dat";
    //std::string iopath = "/afs/cern.ch/work/h/hchou/AMSData/test18";
    Int_t idx = std::atoi(argv[1]);

    //Hist::Load("hit_fill.root", "dat");
    Hist::Load("hit_fill.root", iopath);

    // Fit
    Hist* hAdc = Hist::Head("hMTFadc");
    std::vector<Hist*> vhAdc = Hist::ProjectAll(HistProj::kY, hAdc);

    const Axis& AXeta = hAdc->xaxis();
    const Axis& AXadc = hAdc->yaxis();
    
    TFile * ofle = new TFile(Form("%s/hit_fit%04d.root", iopath.c_str(), idx), "RECREATE");
    ofle->cd();

    Hist::AddDirectory();
    
    Hist* hAdcK = Hist::New("hAdcK", HistAxis(AXeta, "Kappa"));
    Hist* hAdcM = Hist::New("hAdcM", HistAxis(AXeta, "Mpv"));
    Hist* hAdcS = Hist::New("hAdcS", HistAxis(AXeta, "Sigma"));
    Hist* hAdcF = Hist::New("hAdcF", HistAxis(AXeta, "Fluc"));

    TF1* fcov = new TF1("fcov", flgcov, 0, 100, 5);
    TF1* func = fcov;
    func->SetParameters(1.0, 0.1, 0.0, 1.0, 1.0);
    //for (int it = 1; it <= AXeta.nbin(); ++it) {
    for (int it = idx; it <= idx && it <= AXeta.nbin(); ++it) {
        COUT("Process ITER %d\n", it);
        Double_t kpa = 0.1;
        Double_t mpv = (*vhAdc.at(it))()->GetBinCenter((*vhAdc.at(it))()->GetMaximumBin());
        Double_t rms = (*vhAdc.at(it))()->GetRMS();
        func->SetParameters(1000, kpa, mpv, rms, 0.3*rms);
        func->SetParLimits(1, 0.0, 1.0);
        func->SetParLimits(2, 0.0, 10.0*mpv);
        func->SetParLimits(3, 0.0, 10.0*rms);
        func->SetParLimits(4, 0.0, 10.0*rms);
        //func->FixParameter(4, 0.178087); // TKadcx
        //func->FixParameter(4, 0.166633); // TKadcy
        
        (*vhAdc.at(it))()->Fit(func, "q0", "");
        (*vhAdc.at(it))()->Fit(func, "q0", "");
        (*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(1)-4*func->GetParameter(2), func->GetParameter(1)+6.*func->GetParameter(2));
        (*vhAdc.at(it))()->Fit(func, "q0", "", func->GetParameter(1)-4*func->GetParameter(2), func->GetParameter(1)+6.*func->GetParameter(2));
    
        (*hAdcK)()->SetBinContent(it, func->GetParameter(1));
        (*hAdcK)()->SetBinError  (it, func->GetParError(1));
        (*hAdcM)()->SetBinContent(it, func->GetParameter(2));
        (*hAdcM)()->SetBinError  (it, func->GetParError(2));
        (*hAdcS)()->SetBinContent(it, func->GetParameter(3));
        (*hAdcS)()->SetBinError  (it, func->GetParError(3));
        (*hAdcF)()->SetBinContent(it, func->GetParameter(4));
        (*hAdcF)()->SetBinError  (it, func->GetParError(4));

        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXadc));
        for (int jt = 1; jt <= AXadc.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXadc.center(jt)));
        }
        (*vhAdc.at(it))()->Write();
    }
  
    ofle->Write();
    ofle->Close();

    return 0;
}
