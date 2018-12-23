#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
   
    Hist::Load("mu_fill.root", "dat");

    // Fit
    Hist* hCKR = Hist::Head("hCKRrso");
    Hist* hHCR = Hist::Head("hHCRrso");
    
    Hist* hCKM = Hist::Head("hCKMcut");
    //Hist* hCKM = Hist::Head("hHCM");
    Hist* hHCM = Hist::Head("hHCMcut");
   
    const Axis& AXmom = hHCM->xaxis();
    
    TFile * ofle = new TFile("mu_fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();

    //const Double_t stable = 1.7;
    const Double_t stable = 2.0;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);

    COUT("Rrso\n");
    Hist* hCKRrsoM = Hist::New("hCKRrsoM", HistAxis(AXmom, "(1/Rm - 1/Rt) Mean [1]"));
    Hist* hCKRrsoS = Hist::New("hCKRrsoS", HistAxis(AXmom, "(1/Rm - 1/Rt) Sigma [1]"));
    Hist* hHCRrsoM = Hist::New("hHCRrsoM", HistAxis(AXmom, "(1/Rm - 1/Rt) Mean [1]"));
    Hist* hHCRrsoS = Hist::New("hHCRrsoS", HistAxis(AXmom, "(1/Rm - 1/Rt) Sigma [1]"));
   
    std::vector<Hist*> vhCKR = Hist::ProjectAll(HistProj::kY, hCKR);
    std::vector<Hist*> vhHCR = Hist::ProjectAll(HistProj::kY, hHCR);

    for (int it = 1; it <= AXmom.nbin(); ++it) {
        Double_t sqrtcen = std::sqrt(AXmom.center(it, AxisScale::kLog));
        
        // Choutko
        Double_t CKRmax = (*vhCKR.at(it))()->GetBinCenter((*vhCKR.at(it))()->GetMaximumBin());
        Double_t CKRrms = (*vhCKR.at(it))()->GetRMS();
        gaus->SetParameters(1000, CKRmax, CKRrms);
        (*vhCKR.at(it))()->Fit(gaus, "q0", "");
        (*vhCKR.at(it))()->Fit(gaus, "q0", "", CKRmax-stable*CKRrms, CKRmax+stable*CKRrms);
        (*vhCKR.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKR.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKR.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKR.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hCKRrsoM)()->SetBinContent(it, sqrtcen * gaus->GetParameter(1));
        (*hCKRrsoM)()->SetBinError  (it, sqrtcen * gaus->GetParError(1));
        (*hCKRrsoS)()->SetBinContent(it, sqrtcen * gaus->GetParameter(2));
        (*hCKRrsoS)()->SetBinError  (it, sqrtcen * gaus->GetParError(2));
        
        // Hsin-Yi Chou
        Double_t HCRmax = (*vhHCR.at(it))()->GetBinCenter((*vhHCR.at(it))()->GetMaximumBin());
        Double_t HCRrms = (*vhHCR.at(it))()->GetRMS();
        gaus->SetParameters(1000, 0, HCRrms); // testcode
        //gaus->SetParameters(1000, HCRmax, HCRrms);
        (*vhHCR.at(it))()->Fit(gaus, "q0", "");
        (*vhHCR.at(it))()->Fit(gaus, "q0", "", HCRmax-stable*HCRrms, HCRmax+stable*HCRrms);
        (*vhHCR.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCR.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCR.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCR.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hHCRrsoM)()->SetBinContent(it, sqrtcen * gaus->GetParameter(1));
        (*hHCRrsoM)()->SetBinError  (it, sqrtcen * gaus->GetParError(1));
        (*hHCRrsoS)()->SetBinContent(it, sqrtcen * gaus->GetParameter(2));
        (*hHCRrsoS)()->SetBinError  (it, sqrtcen * gaus->GetParError(2));
    }

    for (int it = 1; it <= AXmom.nbin(); ++it) {
        vhCKR.at(it)->style(Fill(), Line(kBlue), Marker(kBlue));
        vhHCR.at(it)->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chRrso = Hist::Collect(Form("chRrso%04d", it), HistList({ vhCKR.at(it), vhHCR.at(it) }));
        chRrso->Write();
    }
        
    hCKRrsoM->write();
    hHCRrsoM->write();
    hCKRrsoS->write();
    hHCRrsoS->write();
   
    hCKRrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCRrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrsoM = Hist::Collect("chRrsoM", HistList({ hCKRrsoM, hHCRrsoM }));
    chRrsoM->Write();
    
    hCKRrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCRrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrsoS = Hist::Collect("chRrsoS", HistList({ hCKRrsoS, hHCRrsoS }));
    chRrsoS->Write();


    COUT("Mrso\n");
    Hist* hCKMrsoM = Hist::New("hCKMrsoM", HistAxis(AXmom, "Mass Mean [1]"));
    Hist* hCKMrsoS = Hist::New("hCKMrsoS", HistAxis(AXmom, "Mass Sigma [1]"));
    Hist* hHCMrsoM = Hist::New("hHCMrsoM", HistAxis(AXmom, "Mass Mean [1]"));
    Hist* hHCMrsoS = Hist::New("hHCMrsoS", HistAxis(AXmom, "Mass Sigma [1]"));
   
    std::vector<Hist*> vhCKM = Hist::ProjectAll(HistProj::kY, hCKM);
    std::vector<Hist*> vhHCM = Hist::ProjectAll(HistProj::kY, hHCM);

    for (int it = 1; it <= AXmom.nbin(); ++it) {
        // Choutko
        Double_t CKMmax = (*vhCKM.at(it))()->GetBinCenter((*vhCKM.at(it))()->GetMaximumBin());
        Double_t CKMrms = (*vhCKM.at(it))()->GetRMS();
        gaus->SetParameters(1000, CKMmax, CKMrms);
        (*vhCKM.at(it))()->Fit(gaus, "q0", "");
        (*vhCKM.at(it))()->Fit(gaus, "q0", "", CKMmax-stable*CKMrms, CKMmax+stable*CKMrms);
        (*vhCKM.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKM.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKM.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKM.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hCKMrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hCKMrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hCKMrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hCKMrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        // Hsin-Yi Chou
        Double_t HCMmax = (*vhHCM.at(it))()->GetBinCenter((*vhHCM.at(it))()->GetMaximumBin());
        Double_t HCMrms = (*vhHCM.at(it))()->GetRMS();
        gaus->SetParameters(1000, 0, HCMrms); // testcode
        //gaus->SetParameters(1000, HCMmax, HCMrms);
        (*vhHCM.at(it))()->Fit(gaus, "q0", "");
        (*vhHCM.at(it))()->Fit(gaus, "q0", "", HCMmax-stable*HCMrms, HCMmax+stable*HCMrms);
        (*vhHCM.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCM.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCM.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCM.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hHCMrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hHCMrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hHCMrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hHCMrsoS)()->SetBinError  (it, gaus->GetParError(2));
    }

    for (int it = 1; it <= AXmom.nbin(); ++it) {
        vhCKM.at(it)->style(Fill(), Line(kBlue), Marker(kBlue));
        vhHCM.at(it)->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chMrso = Hist::Collect(Form("chMrso%04d", it), HistList({ vhCKM.at(it), vhHCM.at(it) }));
        chMrso->Write();
    }
        
    hCKMrsoM->write();
    hHCMrsoM->write();
    hCKMrsoS->write();
    hHCMrsoS->write();
   
    hCKMrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCMrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMrsoM = Hist::Collect("chMrsoM", HistList({ hCKMrsoM, hHCMrsoM }));
    chMrsoM->Write();
    
    hCKMrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCMrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMrsoS = Hist::Collect("chMrsoS", HistList({ hCKMrsoS, hHCMrsoS }));
    chMrsoS->Write();

    ofle->Write();
    ofle->Close();

    return 0;
}
