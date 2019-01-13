#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

Double_t funclg (Double_t* x, Double_t* par) { return par[0] * TrackSys::LandauGaus::Func(x[0], par[1], par[2], par[3], par[4]); }

Double_t funckpa(Double_t* x, Double_t* par) { long double ibta = std::hypot(1.0, x[0]); return TrackSys::IonEloss::FuncKpa(ibta, std::array<long double, 5>({par[0], par[1], par[2], par[3], par[4]})); }
Double_t funcmpv(Double_t* x, Double_t* par) { long double ibta = std::hypot(1.0, x[0]); return TrackSys::IonEloss::FuncMpv(ibta, std::array<long double, 4>({par[0], par[1], par[2], par[3]})); }

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();

    //Hist::Load("PR_hit_fill.root", "dat");
    Hist::Load("HE4_hit_fill.root", "dat");

    // Fit
    //Hist* hQ = Hist::Head("hTKqxy");
    Hist* hQ = Hist::Head("hTFq");
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
    
    TF1* fkpa = new TF1("fkpa", funckpa, 0,   1, 5);
    TF1* fmpv = new TF1("fmpv", funcmpv, 0, 100, 4);
    TF1* fsgm = new TF1("fsgm", funcmpv, 0, 100, 4);
    TF1* fmod = new TF1("fmod", funcmpv, 0, 100, 4);
    
    //TF1* fkpa = new TF1("fkpa", "0.5*(1+TMath::Erf([0]*log(x*x)-[1])) + [2]*0.5*TMath::Erfc([3]*log(x*x)+[4])");
    //TF1* fmpv = new TF1("fmpv", "[0] + [1] * (1+x*x) - [2] * TMath::Log([3]+(x*x))");
    //TF1* fsgm = new TF1("fsgm", "[0] + [1] * (1+x*x) - [2] * TMath::Log([3]+(x*x))");
    //TF1* fmod = new TF1("fmod", "[0] + [1] * (1+x*x) - [2] * TMath::Log([3]+(x*x))");

    // Proton (TKq)
    //fkpa->SetParameters(3.08975e-01, 1.35977e+00, 7.18478e-03, 1.98220e-01, 8.06396e-01);
    //fmpv->SetParameters(5.28872e-02, 6.90213e-01, 5.71953e-02, 3.86236e-02);
    //fsgm->SetParameters(2.03321e-02, 7.09628e-02, 0.0, 0.0);
    //fmod->SetParameters(9.06801e-02, 6.83408e-01, 5.93934e-02, 4.14200e-02);
    //Double_t fluc = 0.089323;
   
    // He4 (TKq)
    //fkpa->SetParameters(5.37172e-01, 8.19362e-01, 5.63577e-03, 6.49063e-01, 1.87820e+00);
    //fmpv->SetParameters(5.52093e-02, 3.06309e+00, 1.70756e-01, 1.82024e-02);
    //fsgm->SetParameters(0.0, 3.27483e-01, 0.0, 0.0);
    //fmod->SetParameters(6.31500e-01, 2.76368e+00, 8.55058e-02, 4.55644e-03);
    //Double_t fluc = 0.174984;
    
    // Proton (TFq)
    //fkpa->SetParameters(8.77160e-01, 7.84314e-01, 5.45507e-03, 2.84523e+01, 1.61855e+01);
    //fmpv->SetParameters(6.63128e-02, 8.60609e-01, 1.57473e-02, 3.90757e-03);
    //fsgm->SetParameters(4.48021e-03, 7.18446e-02, 0.0, 0.0);
    //fmod->SetParameters(1.10112e-01, 8.47921e-01, 1.60048e-02, 4.16204e-03);
    //Double_t fluc = 0.0810673;
   
    // He4 (TFq)
    fkpa->SetParameters(4.75812e-01, 2.92518e-01, 8.86139e-02, 0.00000e+00, 4.92484e-02);
    fmpv->SetParameters(2.01919e+00, 2.32107e+00, 1.23154e-02, 4.92361e-05);
    fsgm->SetParameters(-4.90115e-02, 2.06891e-01, 9.39886e-02, 1.13678e-01);
    fmod->SetParameters(2.18181e+00, 2.22977e+00, 7.31254e-03, 3.37621e-06);
    Double_t fluc = 0.1588111;
   
    TF1* func = new TF1("func", funclg, 0, 10, 5);
    func->SetNpx(100000);

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
       
        //func->FixParameter(1, fkpa->Eval(igb));
        //func->FixParameter(2, fmpv->Eval(igb));
        //func->FixParameter(3, fsgm->Eval(igb));
        func->FixParameter(4, fluc);

        (*vhQ.at(it))()->Fit(func, "q0", "");
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+7*rms);
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+6*rms);
        
        mpv = func->GetParameter(2);
        rms = std::hypot(func->GetParameter(3), func->GetParameter(4));
        (*vhQ.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+6*rms);
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
