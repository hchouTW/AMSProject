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
    if (par[4] <= 1.0e-08) val = par[0] * flg->Eval(xv);
    return val;
}


int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
 
    int ibin = std::atoi(argv[1]);

    Hist::Load("hit_fill.root", "/ams_home/hchou/AMSProject/subj/trsys/dat");

    // Fit
    Hist* hSqrQ = Hist::Head("hTKqx");
    //Hist* hSqrQ = Hist::Head("hTKqy");
    //Hist* hSqrQ = Hist::Head("hTFq");
    std::vector<Hist*>&& vhSqrQ = Hist::ProjectAll(HistProj::kY, hSqrQ);

    const Axis& AXigb = hSqrQ->xaxis();
    const Axis& AXsqr = hSqrQ->yaxis();
    
    TFile * ofle = new TFile(Form("/ams_home/hchou/AMSProject/subj/trsys/dat/hitq_fit%04d.root", ibin), "RECREATE");
    ofle->cd();

    TF1* fkpa = new TF1("fkpa", "0.5 * (1.0 + TMath::Erf([0] * TMath::Log(1+[1]*(1+x*x)^[2]) - [3]))");
    //fkpa->SetParameters(3.83633e+02, 1.22278e-03, 6.08278e-01, 2.28694e+00); // PR: TKqx
    //fkpa->SetParameters(9.94116e+00, 1.91792e+00, 1.20427e-01, 1.27033e+01); // PR: TKqy
    //fkpa->SetParameters(2.94768e+01, 2.21358e+01, 6.06181e-02, 9.45926e+01); // PR: TFq
    fkpa->SetParameters(4.09744e+04, 2.01209e-03, 5.41564e-03, 8.34717e+01); // HE4: TKqx
    //fkpa->SetParameters(3.80213e+00, 5.81308e+03, 4.34709e-01, 3.47169e+01); // HE4: TKqy
    //fkpa->SetParameters(3.17895e+00, 1.64932e+04, 5.05553e-01, 3.19234e+01); // HE4: TFq

    Hist* hSqrQK = Hist::New("hSqrQK", HistAxis(AXigb, "Kappa"));
    Hist* hSqrQM = Hist::New("hSqrQM", HistAxis(AXigb, "Mpv"));
    Hist* hSqrQS = Hist::New("hSqrQS", HistAxis(AXigb, "Sigma"));
    Hist* hSqrQF = Hist::New("hSqrQF", HistAxis(AXigb, "Fluc"));

    TF1* func = new TF1("func", flgcov, 0, 10, 5);
    func->SetNpx(1000);
    //for (int it = 1; it <= AXigb.nbin(); ++it) {
    for (int it = ibin; it <= ibin; ++it) {
        COUT("Process ITER %d\n", it);
        Double_t igb = AXigb.center(it, AxisScale::kLog);
        Double_t kpa = 0.01;
        Double_t mpv = (*vhSqrQ.at(it))()->GetBinCenter((*vhSqrQ.at(it))()->GetMaximumBin());
        Double_t rms = (*vhSqrQ.at(it))()->GetRMS();
        func->SetParameters(1000, kpa, mpv, 0.3043236, 0.372006);
        func->SetParLimits(1, 0.0, 1.0);
        func->SetParLimits(2, 0.0, 5.0*mpv);
        func->SetParLimits(3, 0.0, 20.0*rms);
        func->SetParLimits(4, 0.0, 20.0*rms);
        func->FixParameter(1, fkpa->Eval(igb));
        //func->FixParameter(4, 0.135); // PR: TKqx
        //func->FixParameter(4, 0.15);  // PR: TKqy
        //func->FixParameter(4, 0.082); // PR: TFq
        func->FixParameter(4, 0.08); // HE4: TKqx
        //func->FixParameter(4, 0.49);  // HE4: TKqy
        //func->FixParameter(4, 0.0); // HE4: TFq
    
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-5*rms, mpv+10*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));

        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-4*rms, mpv+12*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));

        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+10*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+8*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+6*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
    
    
        (*hSqrQK)()->SetBinContent(it, func->GetParameter(1));
        (*hSqrQK)()->SetBinError  (it, func->GetParError(1));
        (*hSqrQM)()->SetBinContent(it, func->GetParameter(2));
        (*hSqrQM)()->SetBinError  (it, func->GetParError(2));
        (*hSqrQS)()->SetBinContent(it, func->GetParameter(3));
        (*hSqrQS)()->SetBinError  (it, func->GetParError(3));
        (*hSqrQF)()->SetBinContent(it, func->GetParameter(4));
        (*hSqrQF)()->SetBinError  (it, func->GetParError(4));

        (*vhSqrQ.at(it))()->Write();
        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXsqr));
        for (int jt = 1; jt <= AXsqr.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXsqr.center(jt)));
        }

        (vhSqrQ.at(it))->style(Fill(), Line(kBlue), Marker(kBlue));
        tmpl->style(Fill(), Line(kRed), Marker(kRed));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhSqrQ.at(it), tmpl }));
        ch->Write();
    }
 
    hSqrQK->Write();
    hSqrQM->Write();
    hSqrQS->Write();
    hSqrQF->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
