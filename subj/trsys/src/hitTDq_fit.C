#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

Double_t funclg (Double_t* x, Double_t* par) { return par[0] * TrackSys::LandauGaus::Func(x[0], par[1], par[2], par[3]); }
Double_t funclgs (Double_t* x, Double_t* par) { return par[0] * (1.0-par[1]) * TrackSys::LandauGaus::Func(x[0], par[2], par[3], par[4]) + (par[0] * par[1]) * TrackSys::SftLandauGaus::Func(x[0], par[5], par[6], par[7], par[8]); }

Double_t funckpa(Double_t* x, Double_t* par) { long double ibta = std::hypot(1.0, x[0]); return TrackSys::IonTrEloss::FuncKpa(ibta, std::array<long double, 8>({par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]})); }
Double_t funcmpv(Double_t* x, Double_t* par) { long double ibta = std::hypot(1.0, x[0]); return TrackSys::IonTrEloss::FuncMpv(ibta, std::array<long double, 8>({par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]})); }
Double_t funcsgm(Double_t* x, Double_t* par) { long double ibta = std::hypot(1.0, x[0]); return TrackSys::IonTrEloss::FuncSgm(ibta, std::array<long double, 8>({par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]})); }

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();

    Hist::Load("hit_fill.root", "dat");

    // Fit
    Hist* hQ = Hist::Head("hTDqb");
    std::vector<Hist*>&& vhQ = Hist::ProjectAll(HistProj::kY, hQ);

    const Axis& AXigb = hQ->xaxis();
    const Axis& AXq   = hQ->yaxis();
    
    TFile* ofle = new TFile("hitTDq_fit.root", "RECREATE");
    ofle->cd();
    
    Hist* hQW   = Hist::New("hQW", HistAxis(AXigb, "Weight"));

    Hist* hQK   = Hist::New("hQK", HistAxis(AXigb, "Kappa"));
    Hist* hQM   = Hist::New("hQM", HistAxis(AXigb, "Mpv"));
    Hist* hQS   = Hist::New("hQS", HistAxis(AXigb, "Sigma"));
    
    Hist* hQSK  = Hist::New("hQSK",  HistAxis(AXigb, "S-Kappa"));
    Hist* hQSM  = Hist::New("hQSM",  HistAxis(AXigb, "S-Mpv"));
    Hist* hQSS  = Hist::New("hQSS",  HistAxis(AXigb, "S-Sigma"));
    Hist* hQSSS = Hist::New("hQSSS", HistAxis(AXigb, "S-Shift"));
    
    TF1* fkpa = new TF1("fkpa", funckpa, 0,   1, 8);
    //TF1* fmpv = new TF1("fmpv", funcmpv, 0, 100, 8);
    //TF1* fsgm = new TF1("fsgm", funcsgm, 0, 100, 8);
   
    fkpa->SetParameters(2.66114e-01, 2.09585e+00, 8.08032e-03, 0.0, 0.0, 0.0, 0.0, 0.0);
    //fmpv->SetParameters(8.44593e-01, 1.14637e+00, 4.44571e-01, 3.36409e-01, 1.50348e+00, 1.08929e+00, 1.50062e-01, 1.58739e+00);
    //fsgm->SetParameters(2.92792e-01, 1.04848e+00, 4.01022e-01, 1.64130e-01, 1.67975e+00, 1.25329e-01, 5.05266e-01, 2.24477e+00);
    
    TF1* fmpv = new TF1("fmpv", "[0] + [1] * (1+x*x)  - [2] * TMath::Log([3]+(x*x))");
    TF1* fsgm = new TF1("fsgm", "[0] + [1] * (1+x*x)  - [2] * TMath::Log([3]+(x*x))");
    fmpv->SetParameters(-2.80726e-01, 1.09615e+00, 9.12295e-02, 4.72636e-06);
    fsgm->SetParameters(2.64028e-03, 3.06854e-01, 2.36350e-02, 9.70293e-06);

    //TF1* func = new TF1("func", funclg, 0, 10, 4);
    TF1* func = new TF1("func", funclgs, 0, 10, 9);
    func->SetNpx(100000);

    for (int it = 1; it <= AXigb.nbin(); ++it) {
        COUT("Process IBIN %d/%d\n", it, AXigb.nbin());
        Double_t igb = AXigb.center(it, AxisScale::kLog);
        Double_t mpv = (*vhQ.at(it))()->GetBinCenter((*vhQ.at(it))()->GetMaximumBin());
        Double_t rms = (*vhQ.at(it))()->GetRMS();
        
        func->SetParameters(1000, 0.0, fkpa->Eval(igb), fmpv->Eval(igb), fsgm->Eval(igb), 5.56454e-02, 7.50501e+00, 1.41060e+00, 6.29173e+00);

        func->SetParLimits(0, 0.0, 10000000.);
        func->SetParLimits(1, 0.0, 1.0);
        
        func->SetParLimits(2, 0.0, 1.0);
        func->SetParLimits(3, 0.0, 100.0*mpv);
        func->SetParLimits(4, 0.0, 100.0*rms);
        
        func->SetParLimits(5, 0.0, 1.0);
        func->SetParLimits(6, 0.0, 100.0*mpv);
        func->SetParLimits(7, 0.0, 100.0*rms);
        func->SetParLimits(8, -5.0, 10.0);
        
        func->FixParameter(2, fkpa->Eval(igb));
        func->FixParameter(3, fmpv->Eval(igb));
        func->FixParameter(4, fsgm->Eval(igb));
        
        func->FixParameter(1, 0.1);
        func->FixParameter(5, 9.45619e-01 * 0.5 * (1.0 + TMath::Erf(1.40853e+03 * igb - 2.50636e+00)) + (1.0 - 9.45619e-01) );
        func->FixParameter(6, 7.8);
        func->FixParameter(7, 1.6);
        func->FixParameter(8, 1.0);
        
        Bool_t onlyIon = (igb > 2.0e-03);
        (*vhQ.at(it))()->Fit(func, "q0", "");
        
        if (!onlyIon) {
            func->ReleaseParameter(1);
            func->SetParLimits(1, 0.0, 1.0);
        }
        else func->FixParameter(1, 0.0);
        func->ReleaseParameter(6);
        func->ReleaseParameter(7);
        (*vhQ.at(it))()->Fit(func, "q0", "");
        
        func->ReleaseParameter(5);
        func->SetParLimits(5, 0.0, 1.0);
        (*vhQ.at(it))()->Fit(func, "q0", "");
        
        func->ReleaseParameter(3);
        func->ReleaseParameter(4);
        (*vhQ.at(it))()->Fit(func, "q0", "");
        
        func->ReleaseParameter(2);
        func->SetParLimits(2, 0.0, 1.0);
        (*vhQ.at(it))()->Fit(func, "q0", "");
       
        if (!onlyIon) {
        func->ReleaseParameter(8);
        func->SetParLimits(8, -5.0, 10.0);
        (*vhQ.at(it))()->Fit(func, "q0", "");
        }

        (*vhQ.at(it))()->Fit(func, "q0", "");
        
        
        //Double_t igb1000sqr = (igb * 1.0e+03) * (igb * 1.0e+03);
        ////func->FixParameter(5, 4.20);
        ////func->FixParameter(6, 0.42);
        //func->FixParameter(5, (1.16098e+00 * igb1000sqr + 3.01886e-01) * (1.16098e+00 * igb1000sqr + 3.01886e-01) + 3.96549e+00);
        //func->FixParameter(6, (3.14549e-01 * igb1000sqr + 4.85015e-01) * (3.14549e-01 * igb1000sqr + 4.85015e-01) + 1.35373e-01);
        //func->FixParameter(7, 0.60);
        //func->FixParameter(8, 2.33509e-01 * 0.5 * (1.0 + TMath::Erf(8.92451e-01 * std::log(igb * igb) + 1.34014e+01)) + (1.0 - 2.33509e-01) );
      
        //Bool_t onlyIon = (igb > 3.0e-03);
        //if (onlyIon) func->FixParameter(4, 0.0);
        //else         func->ReleaseParameter(4);

        //func->SetParLimits(0, 0.0, 10000000.);
        //func->SetParLimits(4, 0.0, 10000000.);
    
        //Double_t fact = (2.12846e-01 - 1.22219e-01 * igb * 1.0e+03);
        //if (fact <= 0) fact = 0;
        //func->FixParameter(4, fact);

        //(*vhQ.at(it))()->Fit(func, "q0", "");
        //(*vhQ.at(it))()->Fit(func, "q0", "");
        //if (!onlyIon) func->ReleaseParameter(5);
        //if (!onlyIon) func->ReleaseParameter(6);
        
        //(*vhQ.at(it))()->Fit(func, "q0", "");
        //(*vhQ.at(it))()->Fit(func, "q0", "");
        //if (!onlyIon) func->ReleaseParameter(7);
        //if (!onlyIon) func->ReleaseParameter(8);
        
        //(*vhQ.at(it))()->Fit(func, "q0", "");
        //(*vhQ.at(it))()->Fit(func, "q0", "");
        //func->ReleaseParameter(1);
        //func->ReleaseParameter(2);
        //func->ReleaseParameter(3);
        
        //(*vhQ.at(it))()->Fit(func, "q0", "");
        //(*vhQ.at(it))()->Fit(func, "q0", "");

        /*
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-5*rms, mpv+8*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));

        mpv = func->GetParameter(2);
        rms = func->GetParameter(3);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+7*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));

        mpv = func->GetParameter(2);
        rms = func->GetParameter(3);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+6*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));
        
        mpv = func->GetParameter(2);
        rms = func->GetParameter(3);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+5*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));
        
        mpv = func->GetParameter(2);
        rms = func->GetParameter(3);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+4*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));
        
        mpv = func->GetParameter(2);
        rms = func->GetParameter(3);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+4*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3));
        */
        
        (*hQW)()->SetBinContent(it, func->GetParameter(1));
        (*hQW)()->SetBinError  (it, func->GetParError(1));

        (*hQK)()->SetBinContent(it, func->GetParameter(2));
        (*hQK)()->SetBinError  (it, func->GetParError(2));
        (*hQM)()->SetBinContent(it, func->GetParameter(3));
        (*hQM)()->SetBinError  (it, func->GetParError(3));
        (*hQS)()->SetBinContent(it, func->GetParameter(4));
        (*hQS)()->SetBinError  (it, func->GetParError(4));
        
        (*hQSK)()->SetBinContent(it, func->GetParameter(5));
        (*hQSK)()->SetBinError  (it, func->GetParError(5));
        (*hQSM)()->SetBinContent(it, func->GetParameter(6));
        (*hQSM)()->SetBinError  (it, func->GetParError(6));
        (*hQSS)()->SetBinContent(it, func->GetParameter(7));
        (*hQSS)()->SetBinError  (it, func->GetParError(7));
        (*hQSSS)()->SetBinContent(it, func->GetParameter(8));
        (*hQSSS)()->SetBinError  (it, func->GetParError(8));
        
        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXq));
        for (int jt = 1; jt <= AXq.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXq.center(jt)));
        }

        (vhQ.at(it))->style(Line(kBlue), Marker(kBlue));
        tmpl->style(Line(kRed), Marker(kRed));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhQ.at(it), tmpl }));
        ch->Write();
    }

    hQW->write();
    
    hQK->write();
    hQM->write();
    hQS->write();

    hQSK->write();
    hQSM->write();
    hQSS->write();
    hQSSS->write();

    ofle->Write();
    ofle->Close();

    return 1;
}
