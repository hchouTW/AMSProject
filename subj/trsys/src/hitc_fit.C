#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
    
int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
 
    //Hist::Load("PR_hit_fill.root", "/ams_home/hchou/AMSProject/subj/trsys/dat");
    Hist::Load("HE4_hit_fill.root", "/ams_home/hchou/AMSProject/subj/trsys/dat");

    // Fit
    Hist* hCoo = Hist::Head("hMryInn");
    std::vector<Hist*>&& vhCoo = Hist::ProjectAll(HistProj::kY, hCoo);

    const Axis& AXigb = hCoo->xaxis();
    const Axis& AXtme = hCoo->yaxis();
    
    TFile * ofle = new TFile("hitc_fit.root", "RECREATE");
    ofle->cd();

    Hist* hCooM = Hist::New("hCooM", HistAxis(AXigb, "Mean"));
    Hist* hCooS = Hist::New("hCooS", HistAxis(AXigb, "Sigma"));
    Hist* hCooInvS = Hist::New("hCooInvS", HistAxis(AXigb, "1/Sigma"));

    TF1* func = new TF1("func", "gaus");
    func->SetNpx(1000);
    for (int it = 1; it <= AXigb.nbin(); ++it) {
        COUT("Process ITER %d\n", it);
        Double_t igb = AXigb.center(it, AxisScale::kLog);
        Double_t mpv = (*vhCoo.at(it))()->GetBinCenter((*vhCoo.at(it))()->GetMaximumBin());
        Double_t rms = (*vhCoo.at(it))()->GetRMS();
        func->SetParameters(1000, 0.0, 0.3);
    
        (*vhCoo.at(it))()->Fit(func, "q0", "", mpv-4*rms, mpv+4*rms);

        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhCoo.at(it))()->Fit(func, "q0", "", mpv-3*rms, mpv+3*rms);
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhCoo.at(it))()->Fit(func, "q0", "", mpv-2*rms, mpv+2*rms);
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhCoo.at(it))()->Fit(func, "q0", "", mpv-1.5*rms, mpv+1.5*rms);
        
        mpv = func->GetParameter(1);
        rms = func->GetParameter(2);
        (*vhCoo.at(it))()->Fit(func, "q0", "", mpv-1.5*rms, mpv+1.5*rms);
        CERR("FIT == MEN %14.8f SMG %14.8f\n", func->GetParameter(1), func->GetParameter(2));

        (*hCooM)()->SetBinContent(it, func->GetParameter(1));
        (*hCooM)()->SetBinError  (it, func->GetParError(1) );
        (*hCooS)()->SetBinContent(it, func->GetParameter(2));
        (*hCooS)()->SetBinError  (it, func->GetParError(2) );
        (*hCooInvS)()->SetBinContent(it, 1.0/func->GetParameter(2));
        (*hCooInvS)()->SetBinError  (it, 0.0001);

        Hist* tmpl = Hist::New(Form("tmpl%03d", it), HistAxis(AXtme));
        for (int jt = 1; jt <= AXtme.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, func->Eval(AXtme.center(jt)));
        }

        (vhCoo.at(it))->style(Line(kBlue), Marker(kBlue));
        tmpl->style(Line(kRed), Marker(kRed));
        THStack* ch = Hist::Collect(Form("ch%03d", it), HistList({ vhCoo.at(it), tmpl }));
        ch->Write();
    }
    
    for (int it = 1; it <= AXigb.nbin(); ++it) vhCoo.at(it)->write();
 
    hCooM->write();
    hCooS->write();
    hCooInvS->write();

    ofle->Write();
    ofle->Close();

    return 1;
}
