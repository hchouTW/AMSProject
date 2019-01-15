#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

Double_t funclg (Double_t* x, Double_t* par) { return par[0] * TrackSys::LandauGaus::Func(x[0], par[1], par[2], par[3], par[4]); }

Double_t funckpa(Double_t* x, Double_t* par) { long double ibta = std::hypot(1.0, x[0]); return TrackSys::IonTrEloss::FuncKpa(ibta, std::array<long double, 7>({par[0], par[1], par[2], par[3], par[4], par[5], par[6]})); }
Double_t funcmpv(Double_t* x, Double_t* par) { long double ibta = std::hypot(1.0, x[0]); return TrackSys::IonTrEloss::FuncMpv(ibta, std::array<long double, 10>({par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9]})); }

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();

    //Hist::Load("PR_hit_fill.root", "dat");
    Hist::Load("HE4v3_hit_fill.root", "dat");

    // Fit
    Hist* hQ = Hist::Head("hTDq");
    //Hist* hQ = Hist::Head("hTDs");
    std::vector<Hist*>&& vhQ = Hist::ProjectAll(HistProj::kY, hQ);

    const Axis& AXigb = hQ->xaxis();
    const Axis& AXq   = hQ->yaxis();
    
    TFile* ofle = new TFile("hitq_fit.root", "RECREATE");
    ofle->cd();
    
    Hist* hQK = Hist::New("hQK", HistAxis(AXigb, "Kappa"));
    Hist* hQM = Hist::New("hQM", HistAxis(AXigb, "Mpv"));
    Hist* hQS = Hist::New("hQS", HistAxis(AXigb, "Sigma"));
    Hist* hQF = Hist::New("hQF", HistAxis(AXigb, "Fluc"));
    Hist* hQMOD = Hist::New("hQMOD", HistAxis(AXigb, "Mode"));
    
    TF1* fkpa = new TF1("fkpa", funckpa, 0,   1, 7);
    TF1* fmpv = new TF1("fmpv", funcmpv, 0, 100, 10);
    TF1* fsgm = new TF1("fsgm", funcmpv, 0, 100, 10);
    TF1* fmod = new TF1("fmod", funcmpv, 0, 100, 10);
    
    //TF1* fkpa = new TF1("fkpa", "0.5*(1+TMath::Erf([0]*log(x*x)-[1])) + [2]*0.5*TMath::Erfc([3]*log(x*x)+[4]) + [5]*0.5*TMath::Erfc([6]*log(x*x)+[7])");
    //TF1* fmpv = new TF1("fmpv", "[0] + [1] * (1+x*x) - [2] * TMath::Log((x*x)) + [3] * 0.5 * TMath::Erfc([4] * log(x*x) + [5])");
    //TF1* fsgm = new TF1("fsgm", "[0] + [1] * (1+x*x) - [2] * TMath::Log((x*x)) + [3] * 0.5 * TMath::Erfc([4] * log(x*x) + [5])");
    //TF1* fmod = new TF1("fmod", "[0] + [1] * (1+x*x) - [2] * TMath::Log((x*x)) + [3] * 0.5 * TMath::Erfc([4] * log(x*x) + [5])");

    // Proton (TDq)
    //fkpa->SetParameters(4.32812e-01, 1.17304e+00, 6.69675e-02, 1.57961e-01, 1.17835e+00, 3.33984e-01, 5.29717e-01, 7.31641e+00);
    //fmpv->SetParameters(-1.62672e-01, 1.74195e+00, 1.36595e-01, 4.56864e+00, 4.88916e-01, 7.08228e+00);
    //fsgm->SetParameters(-1.84424e-02, 3.21822e-01, 1.82253e-02, 1.30338e+00, 4.80543e-01, 6.78145e+00);
    //fmod->SetParameters(-7.85233e-02, 1.71828e+00, 1.34960e-01, 4.53189e+00, 4.91589e-01, 7.12373e+00);
    //Double_t fluc = 0.2;
    //Bool_t isq = true;
    
    // Proton (TDs)
    //fkpa->SetParameters(0.00000e+00, 5.56791e+00, 3.01595e-01, 3.37129e-01, 4.17765e+00, 0.0, 0.0, 0.0);
    //fmpv->SetParameters(1.10478e-01, 8.66312e-01, 5.61390e-02, 4.35985e+00, 4.58503e-01, 6.43397e+00);
    //fsgm->SetParameters(5.68283e-02, 3.28905e-01, 1.53413e-02, 1.05370e+00, 4.08763e-01, 5.45840e+00);
    //fmod->SetParameters(1.12279e-01, 8.64621e-01, 5.61388e-02, 4.35974e+00, 4.58504e-01, 6.43308e+00);
    //Double_t fluc = 0.0;
    //Bool_t isq = false;
   
    // He4 (TDq)
    //fkpa->SetParameters(5.02636e-01, 6.29782e-01, 3.96850e+00, 2.29894e-02, 1.54215e+00, 4.31925e-01, 5.88217e-01, 7.87117e+00);
    //fmpv->SetParameters(5.17292e-01, 6.72552e+00, 4.93536e-01, 1.73039e+01, 4.75069e-01, 6.86345e+00);
    //fsgm->SetParameters(2.55399e-01, 7.67649e-01, 3.17313e-02, 2.83073e+00, 4.64441e-01, 6.52800e+00);
    //fmod->SetParameters(8.63297e-01, 6.62279e+00, 4.89280e-01, 1.71202e+01, 4.78403e-01, 6.91607e+00);
    ////Double_t fluc = 0.75;
    //Double_t fluc = 0.0;
    //Bool_t isq = true;
    
    // He4 (TDs)
    //fkpa->SetParameters(7.47084e-01, 1.17180e+00, 4.69442e-01, 3.83918e-01, 5.18781e+00, 0.0, 0.0, 0.0);
    //fmpv->SetParameters(1.49998e+00, 2.02300e+00, 5.81187e-02, 7.36798e+00, 4.20838e-01, 5.92585e+00);
    //fsgm->SetParameters(8.24817e-01, 5.10777e-01, 0.00000e+00, 1.69321e+00, 4.01223e-01, 5.71100e+00);
    //fmod->SetParameters(1.50404e+00, 2.01905e+00, 5.81185e-02, 7.36787e+00, 4.20838e-01, 5.92503e+00);
    //Double_t fluc = 0.0;
    //Bool_t isq = false;
    
    // He4 (TDq)
    fkpa->SetParameters(4.64173e-01, 5.76225e-01, 1.95796e-01, 1.24523e-02, 2.52173e-01, 6.21123e-01, 7.99860e+00);
    fmpv->SetParameters(-4.36915e+00, 1.14534e+01, 7.03697e-01, 5.27722e-01, 4.75415e+00, 7.34352e-01, 9.79897e+00, 1.08747e+01, 5.18116e-01, 7.82150e+00);
    fsgm->SetParameters(-8.58152e-01, 2.10421e+00, 4.02927e-01, 2.70797e-02, 6.10239e-01, 7.70214e-01, 1.00874e+01, 1.69976e+00, 4.55862e-01, 6.77782e+00);
    fmod->SetParameters(-4.36915e+00, 1.14534e+01, 7.03697e-01, 5.27722e-01, 4.75415e+00, 7.34352e-01, 9.79897e+00, 1.08747e+01, 5.18116e-01, 7.82150e+00);
    Double_t fluc = 0.0;
   
    TF1* func = new TF1("func", funclg, 0, 40, 5);
    func->SetNpx(200000);
    for (int it = 1; it <= AXigb.nbin(); ++it) {
        COUT("Process IBIN %d/%d\n", it, AXigb.nbin());
        Double_t igb = AXigb.center(it, AxisScale::kLog);
        Double_t mpv = (*vhQ.at(it))()->GetBinCenter((*vhQ.at(it))()->GetMaximumBin());
        Double_t rms = (*vhQ.at(it))()->GetRMS();
        
        //func->SetParameters(1000, fkpa->Eval(igb), fmpv->Eval(igb), fsgm->Eval(igb), fluc);
        func->SetParameters(1000, 0.01, mpv, 0.5*rms, 0.07);

        func->SetParLimits(0, 0.0, 10000000.);
        func->SetParLimits(1, 0.0, 1.0);
        func->SetParLimits(2, 0.0, 100.0*mpv);
        func->SetParLimits(3, 0.0, 100.0*rms);
        func->SetParLimits(4, 0.0, 100.0*rms);
       
        func->FixParameter(1, fkpa->Eval(igb));
        func->FixParameter(2, fmpv->Eval(igb));
        func->FixParameter(3, fsgm->Eval(igb));
        func->FixParameter(4, fluc);
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        //(*vhQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+7*rms);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+4*rms);
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        //(*vhQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+6*rms);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-1.5*rms, mpv+3*rms);
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        //(*vhQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+6*rms);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-1.5*rms, mpv+2.5*rms);
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        //(*vhQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+6*rms);
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-1.5*rms, mpv+2.5*rms);
        CERR("FIT == KPA %14.8f MPV %14.8f SMG %14.8f FLUC %14.8f\n", func->GetParameter(1), func->GetParameter(2), func->GetParameter(3), func->GetParameter(4));
        
        (*hQK)()->SetBinContent(it, func->GetParameter(1));
        (*hQK)()->SetBinError  (it, func->GetParError(1));
        (*hQM)()->SetBinContent(it, func->GetParameter(2));
        (*hQM)()->SetBinError  (it, func->GetParError(2));
        (*hQS)()->SetBinContent(it, func->GetParameter(3));
        (*hQS)()->SetBinError  (it, func->GetParError(3));
        (*hQF)()->SetBinContent(it, func->GetParameter(4));
        (*hQF)()->SetBinError  (it, func->GetParError(4));
        
        (*hQMOD)()->SetBinContent(it, func->GetMaximumX());
        (*hQMOD)()->SetBinError  (it, 0.0001);

        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXq));
        for (int jt = 1; jt <= AXq.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXq.center(jt)));
        }

        (vhQ.at(it))->style(Line(kBlue), Marker(kBlue));
        tmpl->style(Line(kRed), Marker(kRed));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhQ.at(it), tmpl }));
        ch->Write();
    }

    hQK->write();
    hQM->write();
    hQS->write();
    hQF->write();
    hQMOD->write();

    ofle->Write();
    ofle->Close();

    return 1;
}
