#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
   
    Hist::Load("track_fill.root", "dat");

    // Fit
    Hist* hCKRrso = Hist::Head("hCKRrso");
    Hist* hHCRrso = Hist::Head("hHCRrso2");
    
    const Axis& AXmom = hCKRrso->xaxis();
    
    TFile * ofle = new TFile("track_fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();

    COUT("Rrso\n");
    Hist* hCKRrsoM = Hist::New("hCKRrsoM", HistAxis(AXmom));
    Hist* hCKRrsoS = Hist::New("hCKRrsoS", HistAxis(AXmom));
    Hist* hHCRrsoM = Hist::New("hHCRrsoM", HistAxis(AXmom));
    Hist* hHCRrsoS = Hist::New("hHCRrsoS", HistAxis(AXmom));
    
    Hist* hHCCKRrsoS = Hist::New("hHCCKRrsoS", HistAxis(AXmom, "HYChou/Choutko Sigma Ratio [1]"));
   
    //const Double_t stable = 1.7;
    const Double_t stable = 2.0;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    std::vector<Hist*> vhCKRrso = Hist::ProjectAll(HistProj::kY, hCKRrso);
    std::vector<Hist*> vhHCRrso = Hist::ProjectAll(HistProj::kY, hHCRrso);

    for (int it = 1; it <= AXmom.nbin(); ++it) {
        Double_t bincen = std::sqrt(AXmom.center(it, AxisScale::kLog));

        // Choutko
        Double_t CKRmax = (*vhCKRrso.at(it))()->GetBinCenter((*vhCKRrso.at(it))()->GetMaximumBin());
        Double_t CKRrms = (*vhCKRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, CKRmax, CKRrms);
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", CKRmax-stable*CKRrms, CKRmax+stable*CKRrms);
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hCKRrsoM)()->SetBinContent(it, bincen * gaus->GetParameter(1));
        (*hCKRrsoM)()->SetBinError  (it, bincen * gaus->GetParError(1));
        (*hCKRrsoS)()->SetBinContent(it, bincen * gaus->GetParameter(2));
        (*hCKRrsoS)()->SetBinError  (it, bincen * gaus->GetParError(2));
        
        // Hsin-Yi Chou
        Double_t HCRmax = (*vhHCRrso.at(it))()->GetBinCenter((*vhHCRrso.at(it))()->GetMaximumBin());
        Double_t HCRrms = (*vhHCRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, 0, HCRrms); // testcode
        //gaus->SetParameters(1000, HCRmax, HCRrms);
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", HCRmax-stable*HCRrms, HCRmax+stable*HCRrms);
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hHCRrsoM)()->SetBinContent(it, bincen * gaus->GetParameter(1));
        (*hHCRrsoM)()->SetBinError  (it, bincen * gaus->GetParError(1));
        (*hHCRrsoS)()->SetBinContent(it, bincen * gaus->GetParameter(2));
        (*hHCRrsoS)()->SetBinError  (it, bincen * gaus->GetParError(2));
        
        (*hHCCKRrsoS)()->SetBinContent(it, (*hHCRrsoS)()->GetBinContent(it)/(*hCKRrsoS)()->GetBinContent(it));
    } 
    
    for (int it = 1; it <= AXmom.nbin(); ++it) {
        vhCKRrso.at(it)->style(Line(kBlue), Marker(kBlue));
        vhHCRrso.at(it)->style(Line(kRed), Marker(kRed));
        THStack* cvRrso = Hist::Collect(Form("cvRrso%03d", it), HistList({ vhCKRrso.at(it), vhHCRrso.at(it) }));
        //vhCKRrso.at(it)->Write();
        //vhHCRrso.at(it)->Write();
        cvRrso->Write();
    }
  
    hCKRrsoM->write();
    hHCRrsoM->write();
    hCKRrsoS->write();
    hHCRrsoS->write();
    hHCCKRrsoS->write();
    
    hHCCKRrsoS->style(Line(kRed), Marker(kRed));

    hCKRrsoM->style(Line(kBlue), Marker(kBlue));
    hHCRrsoM->style(Line(kRed), Marker(kRed));
    THStack* chRrsoM = Hist::Collect("chRrsoM", HistList({ hCKRrsoM, hHCRrsoM }));
    chRrsoM->Write();
    
    hCKRrsoS->style(Line(kBlue), Marker(kBlue));
    hHCRrsoS->style(Line(kRed), Marker(kRed));
    THStack* chRrsoS = Hist::Collect("chRrsoS", HistList({ hCKRrsoS, hHCRrsoS }));
    chRrsoS->Write();

    ofle->Write();
    ofle->Close();

    return 0;
}
