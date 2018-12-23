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

    // Gamma
    //Double_t gval = TMath::Abs(par[5]) * (TMath::Power(xv,TMath::Abs(par[6])-1) * TMath::Exp(-TMath::Abs(par[7])*xv)) * (1.0+TMath::Erf(par[8]*xv-par[9]));
    //Double_t gval = TMath::Abs(par[5]) * (TMath::Power(xv,TMath::Abs(par[6])-1) * TMath::Exp(-TMath::Abs(par[7])*xv));
    //val += gval;
    ///////////

    return val;
}


int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();

    if (argc != 2) return false;
    int ibin = std::atoi(argv[1]);

    Hist::Load("hit_fill.root", "/ams_home/hchou/AMSProject/subj/trsys/dat");

    // Fit
    Hist* hSqrQ = Hist::Head("hTKqxy");
    //Hist* hSqrQ = Hist::Head("hTFq");
    //Hist* hSqrQ = Hist::Head("hTDq");
    std::vector<Hist*>&& vhSqrQ = Hist::ProjectAll(HistProj::kY, hSqrQ);

    const Axis& AXigb = hSqrQ->xaxis();
    const Axis& AXsqr = hSqrQ->yaxis();
    
    TFile * ofle = new TFile(Form("/ams_home/hchou/AMSProject/subj/trsys/dat/hitq_fit%04d.root", ibin), "RECREATE");
    //TFile * ofle = new TFile("/ams_home/hchou/AMSProject/subj/trsys/dat/hitq_fit.root", "RECREATE");
    ofle->cd();

    TF1* fkpa = new TF1("fkpa", "0.5 * (1.0 + TMath::Erf([0] * TMath::Log(1+[1]*(1+x*x)^[2]) - [3]))");
    fkpa->SetParameters(2.16299e+00, 4.12329e+03, 2.57464e-01, 1.96535e+01); // PR: TKqxy
    //fkpa->SetParameters(2.94768e+01, 2.21358e+01, 6.06181e-02, 9.45926e+01); // PR: TFq
    //fkpa->SetParameters(1.18528e-01, 1.43470e+04, 1.37285e+01, 2.95944e+00); // HE4: TKqxy
    //fkpa->SetParameters(3.17895e+00, 1.64932e+04, 5.05553e-01, 3.19234e+01); // HE4: TFq
    
    //TF1* fmpv = new TF1("fmpv", "[0] * (1+x*x)^[3] * ([1] + [2]*(1+x*x)^(-[3]) - TMath::Log([4]+(x*x)^[5]))");
    //fmpv->SetParameters(1.60568e-01, 6.77999e+00, -3.45644e+00, 1.01775e+00, 4.73970e-04, 6.58263e-01); // PR: TDq
    //
    //TF1* fsgm = new TF1("fsgm", "[0] * (1+x*x)^[3] * ([1] + [2]*(1+x*x)^(-[3]) - TMath::Log([4]+(x*x)^[5]))");
    //fsgm->SetParameters(9.97392e-02, 8.32891e+02, -8.29976e+02, 7.38042e-03, 1.57741e-03, 4.73758e-01); // PR: TDq

    Hist* hSqrQK = Hist::New("hSqrQK", HistAxis(AXigb, "Kappa"));
    Hist* hSqrQM = Hist::New("hSqrQM", HistAxis(AXigb, "Mpv"));
    Hist* hSqrQS = Hist::New("hSqrQS", HistAxis(AXigb, "Sigma"));
    Hist* hSqrQF = Hist::New("hSqrQF", HistAxis(AXigb, "Fluc"));
    
    Hist* hSqrQW = Hist::New("hSqrQW", HistAxis(AXigb, "Wgt"));
    Hist* hSqrQA = Hist::New("hSqrQA", HistAxis(AXigb, "Alp"));
    Hist* hSqrQB = Hist::New("hSqrQB", HistAxis(AXigb, "Bta"));
    
    Hist* hSqrQC = Hist::New("hSqrQC", HistAxis(AXigb, "C"));
    Hist* hSqrQD = Hist::New("hSqrQD", HistAxis(AXigb, "D"));

    TF1* func = new TF1("func", flgcov, 0, 10, 5);
    //TF1* func = new TF1("func", flgcov, 0, 10, 10); // GM
    func->SetNpx(1000);
    //for (int it = 1; it <= AXigb.nbin(); ++it) {
    for (int it = ibin; it <= ibin; ++it) {
        COUT("Process ITER %d\n", it);
        Double_t igb = AXigb.center(it, AxisScale::kLog);
        Double_t kpa = 0.01;
        Double_t mpv = (*vhSqrQ.at(it))()->GetBinCenter((*vhSqrQ.at(it))()->GetMaximumBin());
        Double_t rms = (*vhSqrQ.at(it))()->GetRMS();
        func->SetParameters(1000, kpa, mpv, 0.3043236, 0.372006,           2.44706e+01, 1.8, 0.2, 0.6, 4.0);
        func->SetParLimits(1, 0.0, 1.0);
        func->SetParLimits(2, 0.0, 5.0*mpv);
        func->SetParLimits(3, 0.0, 20.0*rms);
        func->SetParLimits(4, 0.0, 20.0*rms);
        func->FixParameter(1, fkpa->Eval(igb));
        func->FixParameter(4, 0.0911489);  // PR: TKqxy
        //func->FixParameter(4, 0.082); // PR: TFq
        //func->FixParameter(4, 0.340068); // HE4: TKqxy
        //func->FixParameter(4, 0.0); // HE4: TFq
        
        //Double_t fluc = 4.96717e-01 * (1.0 + TMath::Erf(2.78973e+00 * TMath::Log(igb) - 6.87224e-01)); // PR: TDq
        //func->FixParameter(1, 0.0);  // PR: TDq
        ////func->FixParameter(2, fmpv->Eval(igb));  // PR: TDq
        ////func->FixParameter(3, fsgm->Eval(igb));  // PR: TDq
        //func->FixParameter(4, fluc); // PR: TDq
        
       // func->FixParameter(6, 1.62271e+00); // PR: TDq
       // func->FixParameter(7, 1.74546e-01); // PR: TDq
       // func->FixParameter(8, 5.24428e-01); // PR: TDq
       // func->FixParameter(9, 3.33551e+00); // PR: TDq
       // 
       // (*vhSqrQ.at(it))()->Fit(func, "q0", "", 1.5, 30.0);
       // CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f WGT %14.8f ALP %14.8f BTA %14.8f ES %14.8f EM %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4), func->GetParameter(5), func->GetParameter(6), func->GetParameter(7), func->GetParameter(8), func->GetParameter(9));
       //
        //func->ReleaseParameter(8);
        //func->ReleaseParameter(9);
        //(*vhSqrQ.at(it))()->Fit(func, "q0", "");
        //CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f WGT %14.8f ALP %14.8f BTA %14.8f ES %14.8f EM %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4), func->GetParameter(5), func->GetParameter(6), func->GetParameter(7), func->GetParameter(8), func->GetParameter(9));
   
         
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-5*rms, mpv+8*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));

        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+7*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));

        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+6*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+5*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+4*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+4*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
       
        // Only for TRD
        //mpv = func->GetParameter(2);
        //rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        //(*vhSqrQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+3*rms);
        //CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        ///////////////
        
    
        (*hSqrQK)()->SetBinContent(it, func->GetParameter(1));
        (*hSqrQK)()->SetBinError  (it, func->GetParError(1));
        (*hSqrQM)()->SetBinContent(it, func->GetParameter(2));
        (*hSqrQM)()->SetBinError  (it, func->GetParError(2));
        (*hSqrQS)()->SetBinContent(it, func->GetParameter(3));
        (*hSqrQS)()->SetBinError  (it, func->GetParError(3));
        (*hSqrQF)()->SetBinContent(it, func->GetParameter(4));
        (*hSqrQF)()->SetBinError  (it, func->GetParError(4));
        
        (*hSqrQW)()->SetBinContent(it, func->GetParameter(5));
        (*hSqrQW)()->SetBinError  (it, func->GetParError(5));
        (*hSqrQA)()->SetBinContent(it, func->GetParameter(6));
        (*hSqrQA)()->SetBinError  (it, func->GetParError(6));
        (*hSqrQB)()->SetBinContent(it, func->GetParameter(7));
        (*hSqrQB)()->SetBinError  (it, func->GetParError(7));
        
        (*hSqrQC)()->SetBinContent(it, func->GetParameter(8));
        (*hSqrQC)()->SetBinError  (it, func->GetParError(8));
        (*hSqrQD)()->SetBinContent(it, func->GetParameter(9));
        (*hSqrQD)()->SetBinError  (it, func->GetParError(9));

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
    
    hSqrQW->Write();
    hSqrQA->Write();
    hSqrQB->Write();
    
    hSqrQC->Write();
    hSqrQD->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
