#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
    
int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
 
    Hist::Load("PR_hit_fill.root", "/ams_home/hchou/AMSProject/subj/trsys/dat");
    //Hist::Load("HE4_hit_fill.root", "/ams_home/hchou/AMSProject/subj/trsys/dat");

    // Fit
    Hist* hTme = Hist::Head("hTFtme");
    std::vector<Hist*>&& vhTme = Hist::ProjectAll(HistProj::kY, hTme);

    const Axis& AXigb = hTme->xaxis();
    const Axis& AXtme = hTme->yaxis();
    
    TFile * ofle = new TFile("hitt_fit.root", "RECREATE");
    ofle->cd();

    Hist* hTmeM = Hist::New("hTmeM", HistAxis(AXigb, "Mean"));
    Hist* hTmeS = Hist::New("hTmeS", HistAxis(AXigb, "Sigma"));
    Hist* hTmeInvS = Hist::New("hTmeInvS", HistAxis(AXigb, "1/Sigma"));

    TF1* func = new TF1("func", "gaus");
    func->SetNpx(1000);
    for (int it = 1; it <= AXigb.nbin(); ++it) {
        Double_t igb = AXigb.center(it, AxisScale::kLog);
        Double_t mpv = (*vhTme.at(it))()->GetBinCenter((*vhTme.at(it))()->GetMaximumBin());
        Double_t rms = (*vhTme.at(it))()->GetRMS();
        func->SetParameters(1000, 0.0, 0.3);
    
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-4*rms, mpv+4*rms);

        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+3*rms);
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+2*rms);
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-1.5*rms, mpv+1.5*rms);
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhTme.at(it))()->Fit(func, "q0", "", mpv-1.5*rms, mpv+1.5*rms);
        CERR("FIT == MEN %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2));

        (*hTmeM)()->SetBinContent(it, func->GetParameter(1) * TMath::Sqrt(0.5));
        (*hTmeM)()->SetBinError  (it, func->GetParError(1)  * TMath::Sqrt(0.5));
        (*hTmeS)()->SetBinContent(it, func->GetParameter(2) * TMath::Sqrt(0.5));
        (*hTmeS)()->SetBinError  (it, func->GetParError(2)  * TMath::Sqrt(0.5));
        (*hTmeInvS)()->SetBinContent(it, 1.0 / func->GetParameter(2) / TMath::Sqrt(0.5));
        (*hTmeInvS)()->SetBinError  (it, 0.0001);

        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXtme));
        for (int jt = 1; jt <= AXtme.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXtme.center(jt)));
        }

        (vhTme.at(it))->style(Line(kBlue), Marker(kBlue));
        tmpl->style(Line(kRed), Marker(kRed));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhTme.at(it), tmpl }));
        ch->Write();
    }
        
    for (int it = 1; it <= AXigb.nbin(); ++it) vhTme.at(it)->write();
 
    hTmeM->write();
    hTmeS->write();
    hTmeInvS->write();

    ofle->Write();
    ofle->Close();

    return 1;
}
