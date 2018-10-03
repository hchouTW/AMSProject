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
    Hist* hKFRrso = Hist::Head("hKFRrso");
    Hist* hHCRrso = Hist::Head("hHCRrso");
    
    Hist* hCKBrso = Hist::Head("hCKBrso");
    Hist* hKFBrso = Hist::Head("hKFBrso");
    Hist* hHCBrso = Hist::Head("hHCBrso");
   
    const Axis& AXmom = hCKRrso->xaxis();
    const Axis& AXbta = hCKBrso->xaxis();
    
    TFile * ofle = new TFile("track_fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();

    COUT("Rrso\n");
    Hist* hCKRrsoM = Hist::New("hCKRrsoM", HistAxis(AXmom));
    Hist* hCKRrsoS = Hist::New("hCKRrsoS", HistAxis(AXmom));
    Hist* hKFRrsoM = Hist::New("hKFRrsoM", HistAxis(AXmom));
    Hist* hKFRrsoS = Hist::New("hKFRrsoS", HistAxis(AXmom));
    Hist* hHCRrsoM = Hist::New("hHCRrsoM", HistAxis(AXmom));
    Hist* hHCRrsoS = Hist::New("hHCRrsoS", HistAxis(AXmom));
    
    Hist* hKFCKRrsoS = Hist::New("hKFCKRrsoS", HistAxis(AXmom, "KF/Choutko Sigma Ratio [1]"));
    Hist* hHCCKRrsoS = Hist::New("hHCCKRrsoS", HistAxis(AXmom, "HYChou/Choutko Sigma Ratio [1]"));
   
    //const Double_t stable = 1.7;
    const Double_t stable = 2.0;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    std::vector<Hist*> vhCKRrso = Hist::ProjectAll(HistProj::kY, hCKRrso);
    std::vector<Hist*> vhKFRrso = Hist::ProjectAll(HistProj::kY, hKFRrso);
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
        
        // Kalman Filter
        Double_t KFRmax = (*vhKFRrso.at(it))()->GetBinCenter((*vhKFRrso.at(it))()->GetMaximumBin());
        Double_t KFRrms = (*vhKFRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, KFRmax, KFRrms);
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", KFRmax-stable*KFRrms, KFRmax+stable*KFRrms);
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hKFRrsoM)()->SetBinContent(it, bincen * gaus->GetParameter(1));
        (*hKFRrsoM)()->SetBinError  (it, bincen * gaus->GetParError(1));
        (*hKFRrsoS)()->SetBinContent(it, bincen * gaus->GetParameter(2));
        (*hKFRrsoS)()->SetBinError  (it, bincen * gaus->GetParError(2));
        
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
        
        (*hKFCKRrsoS)()->SetBinContent(it, (*hKFRrsoS)()->GetBinContent(it)/(*hCKRrsoS)()->GetBinContent(it));
        (*hHCCKRrsoS)()->SetBinContent(it, (*hHCRrsoS)()->GetBinContent(it)/(*hCKRrsoS)()->GetBinContent(it));
    } 
    
    for (int it = 1; it <= AXmom.nbin(); ++it) {
        vhCKRrso.at(it)->style(Fill(), Line(kGreen), Marker(kGreen));
        vhKFRrso.at(it)->style(Fill(), Line(kBlue), Marker(kBlue));
        vhHCRrso.at(it)->style(Fill(), Line(kRed), Marker(kRed));
        THStack* cvRrso = Hist::Collect(Form("cvRrso%03d", it), HistList({ vhCKRrso.at(it), vhKFRrso.at(it), vhHCRrso.at(it) }));
        //vhCKRrso.at(it)->Write();
        //vhKFRrso.at(it)->Write();
        //vhHCRrso.at(it)->Write();
        cvRrso->Write();
    }
  
    hCKRrsoM->write();
    hKFRrsoM->write();
    hHCRrsoM->write();
    hCKRrsoS->write();
    hKFRrsoS->write();
    hHCRrsoS->write();
    hKFCKRrsoS->write();
    hHCCKRrsoS->write();

    hCKRrsoM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFRrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCRrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrsoM = Hist::Collect("chRrsoM", HistList({ hCKRrsoM, hKFRrsoM, hHCRrsoM }));
    chRrsoM->Write();
    
    hKFCKRrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCCKRrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrsoS = Hist::Collect("chRrsoS", HistList({ hKFCKRrsoS, hHCCKRrsoS }));
    chRrsoS->Write();


    COUT("Brso\n");
    Hist* hCKBrsoM = Hist::New("hCKBrsoM", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Mean [1]"));
    Hist* hCKBrsoS = Hist::New("hCKBrsoS", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Sigma [1]"));
    Hist* hKFBrsoM = Hist::New("hKFBrsoM", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Mean [1]"));
    Hist* hKFBrsoS = Hist::New("hKFBrsoS", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Sigma [1]"));
    Hist* hHCBrsoM = Hist::New("hHCBrsoM", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Mean [1]"));
    Hist* hHCBrsoS = Hist::New("hHCBrsoS", HistAxis(AXbta, "(Beta_{rec}/Beta_{gen}-1) Sigma [1]"));
    
    Hist* hKFCKBrsoS = Hist::New("hKFCKBrsoS", HistAxis(AXbta, "KF/Choutko Sigma Ratio [1]"));
    Hist* hHCCKBrsoS = Hist::New("hHCCKBrsoS", HistAxis(AXbta, "HYChou/Choutko Sigma Ratio [1]"));
   
    std::vector<Hist*> vhCKBrso = Hist::ProjectAll(HistProj::kY, hCKBrso);
    std::vector<Hist*> vhKFBrso = Hist::ProjectAll(HistProj::kY, hKFBrso);
    std::vector<Hist*> vhHCBrso = Hist::ProjectAll(HistProj::kY, hHCBrso);

    for (int it = 1; it <= AXbta.nbin(); ++it) {
        // Choutko
        Double_t CKBmax = (*vhCKBrso.at(it))()->GetBinCenter((*vhCKBrso.at(it))()->GetMaximumBin());
        Double_t CKBrms = (*vhCKBrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, CKBmax, CKBrms);
        (*vhCKBrso.at(it))()->Fit(gaus, "q0", "");
        (*vhCKBrso.at(it))()->Fit(gaus, "q0", "", CKBmax-stable*CKBrms, CKBmax+stable*CKBrms);
        (*vhCKBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hCKBrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hCKBrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hCKBrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hCKBrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        // Kalman Filter
        Double_t KFBmax = (*vhKFBrso.at(it))()->GetBinCenter((*vhKFBrso.at(it))()->GetMaximumBin());
        Double_t KFBrms = (*vhKFBrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, KFBmax, KFBrms);
        (*vhKFBrso.at(it))()->Fit(gaus, "q0", "");
        (*vhKFBrso.at(it))()->Fit(gaus, "q0", "", KFBmax-stable*KFBrms, KFBmax+stable*KFBrms);
        (*vhKFBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhKFBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhKFBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhKFBrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hKFBrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hKFBrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hKFBrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hKFBrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
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
        
        (*hKFCKBrsoS)()->SetBinContent(it, (*hKFBrsoS)()->GetBinContent(it)/(*hCKBrsoS)()->GetBinContent(it));
        (*hHCCKBrsoS)()->SetBinContent(it, (*hHCBrsoS)()->GetBinContent(it)/(*hCKBrsoS)()->GetBinContent(it));
    } 
    
    hCKBrsoM->write();
    hKFBrsoM->write();
    hHCBrsoM->write();
    hCKBrsoS->write();
    hKFBrsoS->write();
    hHCBrsoS->write();
    hKFCKBrsoS->write();
    hHCCKBrsoS->write();
   
    hCKBrsoM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFBrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCBrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chBrsoM = Hist::Collect("chBrsoM", HistList({ hCKBrsoM, hKFBrsoM, hHCBrsoM }));
    chBrsoM->Write();
    
    hKFCKBrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCCKBrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chBrsoS = Hist::Collect("chBrsoS", HistList({ hKFCKBrsoS, hHCCKBrsoS }));
    chBrsoS->Write();

    ofle->Write();
    ofle->Close();

    return 0;
}
