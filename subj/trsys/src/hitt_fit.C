#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
    
int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
 
    Hist::Load("hit_fill.root", "/ams_home/hchou/AMSProject/subj/trsys/dat");

    // Fit
    Hist* hTme = Hist::Head("hTFtme");
    std::vector<Hist*>&& vhTme = Hist::ProjectAll(HistProj::kY, hTme);

    const Axis& AXigb = hTme->xaxis();
    const Axis& AXtme = hTme->yaxis();
    
    TFile * ofle = new TFile("/ams_home/hchou/AMSProject/subj/trsys/dat/hitt_fit.root", "RECREATE");
    ofle->cd();

    Hist* hTmeM = Hist::New("hTmeM", HistAxis(AXigb, "Mean"));
    Hist* hTmeS = Hist::New("hTmeS", HistAxis(AXigb, "Sigma"));

    TF1* func = new TF1("func", "gaus");
    func->SetNpx(1000);
    for (int it = 1; it <= AXigb.nbin(); ++it) {
        COUT("Process ITER %d\n", it);
        Double_t igb = AXigb.center(it, AxisScale::kLog);
        Double_t mpv = (*vhTme.at(it))()->GetBinCenter((*vhTme.at(it))()->GetMaximumBin());
        Double_t rms = (*vhTme.at(it))()->GetRMS();
        func->SetParameters(1000, 0.0, 0.3);
    
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-4*rms, mpv+4*rms);
        CERR("FIT == MEN %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2));

        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+3*rms);
        CERR("FIT == MEN %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2));
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+2*rms);
        CERR("FIT == MEN %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2));
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+2*rms);
        CERR("FIT == MEN %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2));
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+2*rms);
        CERR("FIT == MEN %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2));

        (*hTmeM)()->SetBinContent(it, func->GetParameter(1) * TMath::Sqrt(0.5));
        (*hTmeM)()->SetBinError  (it, func->GetParError(1)  * TMath::Sqrt(0.5));
        (*hTmeS)()->SetBinContent(it, func->GetParameter(2) * TMath::Sqrt(0.5));
        (*hTmeS)()->SetBinError  (it, func->GetParError(2)  * TMath::Sqrt(0.5));

        (*vhTme.at(it))()->Write();
        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXtme));
        for (int jt = 1; jt <= AXtme.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXtme.center(jt)));
        }

        (vhTme.at(it))->style(Fill(), Line(kBlue), Marker(kBlue));
        tmpl->style(Fill(), Line(kRed), Marker(kRed));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhTme.at(it), tmpl }));
        ch->Write();
    }
 
    hTmeM->Write();
    hTmeS->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
