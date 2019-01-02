#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
   
    Hist::Load("fill.root", "dat");
    
    // Num
    Hist* hMCnum   = Hist::Head("hMCnum");
    Hist* hCKnum   = Hist::Head("hCKnum");
    Hist* hHCnum   = Hist::Head("hHCnum");
    
    Hist* hCKRDrso = Hist::Head("hCKRDrso");
    Hist* hHCRDrso = Hist::Head("hHCRDrso");

    // Tme
    Hist* hJFevt   = Hist::Head("hJFevt");
    Hist* hHCevt   = Hist::Head("hHCevt");
    Hist* hJFtme   = Hist::Head("hJFtme");
    Hist* hHCtme   = Hist::Head("hHCtme");

    Hist* hJFMrso = Hist::Head("hJFM");
    Hist* hHCMrso = Hist::Head("hHCM2");
    
    const Axis& AXrig = hCKRDrso->xaxis();
    
    TFile * ofle = new TFile("fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();
    
    COUT("rat\n");
    Hist* hCKrat = Hist::Calculate(HistArith::kDivide, "hCKrat", ";Momentum [GeV];Effienecy", hCKnum, hMCnum);
    Hist* hHCrat = Hist::Calculate(HistArith::kDivide, "hHCrat", ";Momentum [GeV];Effienecy", hHCnum, hMCnum);
    
    hCKnum->write();
    hHCnum->write();

    hCKrat->style(Line(kBlue), Marker(kBlue));
    hHCrat->style(Line(kRed), Marker(kRed));
    THStack* chRat = Hist::Collect("chrat", HistList({ hCKrat, hHCrat }));
    chRat->Write();
    
    COUT("tme\n");
    Hist* hJFtmeM = Hist::Calculate(HistArith::kDivide, "hJFtmeM", ";Rigidity [GV];Mean Time [msec]", hJFtme, hJFevt);
    Hist* hHCtmeM = Hist::Calculate(HistArith::kDivide, "hHCtmeM", ";Rigidity [GV];Mean Time [msec]", hHCtme, hHCevt);
    
    hJFtmeM->style(Line(kBlue), Marker(kBlue));
    hHCtmeM->style(Line(kRed), Marker(kRed));
    THStack* chTmeM = Hist::Collect("chTmeM", HistList({ hJFtmeM, hHCtmeM }));
    chTmeM->Write();
    
    const Double_t stable = 1.7;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);

    COUT("CKRDrso\n");
    Hist* hCKRDrsoM = Hist::New("hCKRDrsoM", HistAxis(AXrig, "Mean [1]"));
    Hist* hCKRDrsoS = Hist::New("hCKRDrsoS", HistAxis(AXrig, "Sigma [1]"));
    std::vector<Hist*> vhCKRDrso = Hist::ProjectAll(HistProj::kY, hCKRDrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhCKRDrso.at(it))()->GetBinCenter((*vhCKRDrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhCKRDrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhCKRDrso.at(it))()->Fit(gaus, "q0", "");
        (*vhCKRDrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhCKRDrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKRDrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKRDrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hCKRDrsoM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hCKRDrsoM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hCKRDrsoS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hCKRDrsoS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hCKRDrsoM->write();
    hCKRDrsoS->write();

    COUT("HCRDrso\n");
    Hist* hHCRDrsoM = Hist::New("hHCRDrsoM", HistAxis(AXrig, "Mean [1]"));
    Hist* hHCRDrsoS = Hist::New("hHCRDrsoS", HistAxis(AXrig, "Sigma [1]"));
    std::vector<Hist*> vhHCRDrso = Hist::ProjectAll(HistProj::kY, hHCRDrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhHCRDrso.at(it))()->GetBinCenter((*vhHCRDrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhHCRDrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhHCRDrso.at(it))()->Fit(gaus, "q0", "");
        (*vhHCRDrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhHCRDrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCRDrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCRDrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hHCRDrsoM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hHCRDrsoM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hHCRDrsoS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hHCRDrsoS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hHCRDrsoM->write();
    hHCRDrsoS->write();
    
    hCKRDrsoM->style(Line(kBlue), Marker(kBlue));
    hHCRDrsoM->style(Line(kRed), Marker(kRed));
    THStack* chRDrsoM = Hist::Collect("chRDrsoM", HistList({ hCKRDrsoM, hHCRDrsoM }));
    chRDrsoM->Write();
    
    hCKRDrsoS->style(Line(kBlue), Marker(kBlue));
    hHCRDrsoS->style(Line(kRed), Marker(kRed));
    THStack* chRDrsoS = Hist::Collect("chRDrsoS", HistList({ hCKRDrsoS, hHCRDrsoS }));
    chRDrsoS->Write();
    
    Hist* hHCCKRDrso = Hist::New("hHCCKRDrso", HistAxis(AXrig, "HYChou/Choutko Sigma Ratio [1]"));
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t cen = AXrig.center(it, AxisScale::kLog);
        hHCCKRDrso->fillH1D(cen, (*hHCRDrsoS)()->GetBinContent(it) / (*hCKRDrsoS)()->GetBinContent(it));
        (*hHCCKRDrso)()->SetBinError(it, 0.00001);
    }
    hHCCKRDrso->write();
    
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhCKRDrso.at(it)->style(Line(kBlue), Marker(kBlue));
        vhHCRDrso.at(it)->style(Line(kRed), Marker(kRed));
        THStack* cvhRDrso = Hist::Collect(Form("cvhRDrso%03d", it), HistList({ vhCKRDrso.at(it), vhHCRDrso.at(it) }));
        cvhRDrso->Write();
    }
    
    
    COUT("JFMrso\n");
    Hist* hJFMrsoM = Hist::New("hJFMrsoM", HistAxis(AXrig, "Mean [GeV]"));
    Hist* hJFMrsoS = Hist::New("hJFMrsoS", HistAxis(AXrig, "Sigma [GeV]"));
    std::vector<Hist*> vhJFMrso = Hist::ProjectAll(HistProj::kY, hJFMrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t max = (*vhJFMrso.at(it))()->GetBinCenter((*vhJFMrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhJFMrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhJFMrso.at(it))()->Fit(gaus, "q0", "");
        (*vhJFMrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhJFMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhJFMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhJFMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        //(*vhJFMrso.at(it))()->Scale(1.0/(*vhJFMrso.at(it))()->GetEntries());
    
        (*hJFMrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hJFMrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hJFMrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hJFMrsoS)()->SetBinError  (it, gaus->GetParError(2));
    } 
    hJFMrsoM->write();
    hJFMrsoS->write();

    COUT("HCMrso\n");
    Hist* hHCMrsoM = Hist::New("hHCMrsoM", HistAxis(AXrig, "Mean [GeV]"));
    Hist* hHCMrsoS = Hist::New("hHCMrsoS", HistAxis(AXrig, "Sigma [GeV]"));
    std::vector<Hist*> vhHCMrso = Hist::ProjectAll(HistProj::kY, hHCMrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t max = (*vhHCMrso.at(it))()->GetBinCenter((*vhHCMrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhHCMrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhHCMrso.at(it))()->Fit(gaus, "q0", "");
        (*vhHCMrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhHCMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hHCMrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hHCMrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hHCMrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hHCMrsoS)()->SetBinError  (it, gaus->GetParError(2));
    } 
    hHCMrsoM->write();
    hHCMrsoS->write();
    

    hJFMrsoM->style(Line(kBlue), Marker(kBlue));
    hHCMrsoM->style(Line(kRed), Marker(kRed));
    THStack* chMrsoM = Hist::Collect("chMrsoM", HistList({ hJFMrsoM, hHCMrsoM }));
    chMrsoM->Write();
    
    hJFMrsoS->style(Line(kBlue), Marker(kBlue));
    hHCMrsoS->style(Line(kRed), Marker(kRed));
    THStack* chMrsoS = Hist::Collect("chMrsoS", HistList({ hJFMrsoS, hHCMrsoS }));
    chMrsoS->Write();
    
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhJFMrso.at(it)->style(Line(kBlue), Marker(kBlue));
        vhHCMrso.at(it)->style(Line(kRed), Marker(kRed));
        THStack* cvhMrso = Hist::Collect(Form("cvhMrso%03d", it), HistList({ vhJFMrso.at(it), vhHCMrso.at(it) }));
        cvhMrso->Write();
    }
   
    ofle->Write();
    ofle->Close();

    return 0;
}
