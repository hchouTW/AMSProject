#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
   
    Hist::Load("bta_fill.root", "dat");

    // Fit
    Hist* hTFBrso = Hist::Head("hTFBrso");
    Hist* hHCBrso = Hist::Head("hHCTFBrso");
   
    const Axis& AXbta = hHCBrso->xaxis();
    
    TFile * ofle = new TFile("bta_fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();

    //const Double_t stable = 1.7;
    const Double_t stable = 2.0;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);

    COUT("Brso\n");
    Hist* hTFBrsoM = Hist::New("hTFBrsoM", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Mean [1]"));
    Hist* hTFBrsoS = Hist::New("hTFBrsoS", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Sigma [1]"));
    Hist* hHCBrsoM = Hist::New("hHCBrsoM", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Mean [1]"));
    Hist* hHCBrsoS = Hist::New("hHCBrsoS", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Sigma [1]"));
    
    Hist* hHCTFBrsoS = Hist::New("hHCTFBrsoS", HistAxis(AXbta, "HYChou/Choutko Sigma Ratio [1]"));
   
    std::vector<Hist*> vhTFBrso = Hist::ProjectAll(HistProj::kY, hTFBrso);
    std::vector<Hist*> vhHCBrso = Hist::ProjectAll(HistProj::kY, hHCBrso);

    for (int it = 1; it <= AXbta.nbin(); ++it) {
        // Choutko
        Double_t CKBmax = (*vhTFBrso.at(it))()->GetBinCenter((*vhTFBrso.at(it))()->GetMaximumBin());
        Double_t CKBrms = (*vhTFBrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, CKBmax, CKBrms);
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "");
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", CKBmax-stable*CKBrms, CKBmax+stable*CKBrms);
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hTFBrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hTFBrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hTFBrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hTFBrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        // Hsin-Yi Chou
        Double_t HCBmax = (*vhHCBrso.at(it))()->GetBinCenter((*vhHCBrso.at(it))()->GetMaximumBin());
        Double_t HCBrms = (*vhHCBrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, 0, HCBrms); // testcode
        //gaus->SetParameters(1000, HCBmax, HCBrms);
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "");
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", HCBmax-stable*HCBrms, HCBmax+stable*HCBrms);
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hHCBrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hHCBrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hHCBrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hHCBrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        (*hHCTFBrsoS)()->SetBinContent(it, (*hHCBrsoS)()->GetBinContent(it)/(*hTFBrsoS)()->GetBinContent(it));
    } 
    
    hTFBrsoM->write();
    hHCBrsoM->write();
    hTFBrsoS->write();
    hHCBrsoS->write();
    hHCTFBrsoS->write();
   
    hTFBrsoM->style(Fill(), Line(kGreen), Marker(kBlue));
    hHCBrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chBrsoM = Hist::Collect("chBrsoM", HistList({ hTFBrsoM, hHCBrsoM }));
    chBrsoM->Write();
    
    hHCTFBrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chBrsoS = Hist::Collect("chBrsoS", HistList({ hHCTFBrsoS }));
    chBrsoS->Write();

    ofle->Write();
    ofle->Close();

    return 0;
}
