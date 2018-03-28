#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
   
    //Hist::Load("track_fill.root", "dat");
    Hist::Load("track_fill.root", "/afs/cern.ch/work/h/hchou/AMSData/test23");

    // Fit
    Hist* hCKRrso = Hist::Head("hCKRrso");
    Hist* hKFRrso = Hist::Head("hKFRrso");
    Hist* hHCRrso = Hist::Head("hHCRrso");
    
    Hist* hCKBrso = Hist::Head("hCKBrso");
    Hist* hKFBrso = Hist::Head("hKFBrso");
    Hist* hHCBrso = Hist::Head("hHCBrso");
 
    std::vector<Hist*> hCKresx(9, nullptr);
    std::vector<Hist*> hKFresx(9, nullptr);
    std::vector<Hist*> hHCresx(9, nullptr);
    std::vector<Hist*> hCKresy(9, nullptr);
    std::vector<Hist*> hKFresy(9, nullptr);
    std::vector<Hist*> hHCresy(9, nullptr);
    std::vector<Hist*> hCKcosx(9, nullptr);
    std::vector<Hist*> hKFcosx(9, nullptr);
    std::vector<Hist*> hHCcosx(9, nullptr);
    std::vector<Hist*> hCKcosy(9, nullptr);
    std::vector<Hist*> hKFcosy(9, nullptr);
    std::vector<Hist*> hHCcosy(9, nullptr);

    for (Int_t it = 0; it < 9; ++it) {
        hCKresx.at(it) = Hist::Head(STR("hCKresxL%d", it+1));
        hKFresx.at(it) = Hist::Head(STR("hKFresxL%d", it+1));
        hHCresx.at(it) = Hist::Head(STR("hHCresxL%d", it+1));
        
        hCKresy.at(it) = Hist::Head(STR("hCKresyL%d", it+1));
        hKFresy.at(it) = Hist::Head(STR("hKFresyL%d", it+1));
        hHCresy.at(it) = Hist::Head(STR("hHCresyL%d", it+1));
        
        hCKcosx.at(it) = Hist::Head(STR("hCKcosxL%d", it+1));
        hKFcosx.at(it) = Hist::Head(STR("hKFcosxL%d", it+1));
        hHCcosx.at(it) = Hist::Head(STR("hHCcosxL%d", it+1));
        
        hCKcosy.at(it) = Hist::Head(STR("hCKcosyL%d", it+1));
        hKFcosy.at(it) = Hist::Head(STR("hKFcosyL%d", it+1));
        hHCcosy.at(it) = Hist::Head(STR("hHCcosyL%d", it+1));
    }
   
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
   
    const Double_t stable = 1.7;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    std::vector<Hist*> vhCKRrso = Hist::ProjectAll(HistProj::kY, hCKRrso);
    std::vector<Hist*> vhKFRrso = Hist::ProjectAll(HistProj::kY, hKFRrso);
    std::vector<Hist*> vhHCRrso = Hist::ProjectAll(HistProj::kY, hHCRrso);

    for (int it = 1; it <= AXmom.nbin(); ++it) {
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
    
        (*hCKRrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hCKRrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hCKRrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hCKRrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
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
    
        (*hKFRrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hKFRrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hKFRrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hKFRrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        // Hsin-Yi Chou
        Double_t HCRmax = (*vhHCRrso.at(it))()->GetBinCenter((*vhHCRrso.at(it))()->GetMaximumBin());
        Double_t HCRrms = (*vhHCRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, HCRmax, HCRrms);
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", HCRmax-stable*HCRrms, HCRmax+stable*HCRrms);
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hHCRrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hHCRrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hHCRrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hHCRrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        (*hKFCKRrsoS)()->SetBinContent(it, (*hKFRrsoS)()->GetBinContent(it)/(*hCKRrsoS)()->GetBinContent(it));
        (*hHCCKRrsoS)()->SetBinContent(it, (*hHCRrsoS)()->GetBinContent(it)/(*hCKRrsoS)()->GetBinContent(it));
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
    Hist* hCKBrsoM = Hist::New("hCKBrsoM", HistAxis(AXbta));
    Hist* hCKBrsoS = Hist::New("hCKBrsoS", HistAxis(AXbta));
    Hist* hKFBrsoM = Hist::New("hKFBrsoM", HistAxis(AXbta));
    Hist* hKFBrsoS = Hist::New("hKFBrsoS", HistAxis(AXbta));
    Hist* hHCBrsoM = Hist::New("hHCBrsoM", HistAxis(AXbta));
    Hist* hHCBrsoS = Hist::New("hHCBrsoS", HistAxis(AXbta));
    
    Hist* hKFCKBrsoS = Hist::New("hKFCKBrsoS", HistAxis(AXbta, "KF/Choutko Sigma Ratio [1]"));
    Hist* hHCCKBrsoS = Hist::New("hHCCKBrsoS", HistAxis(AXbta, "HYChou/Choutko Sigma Ratio [1]"));
   
    std::vector<Hist*> vhCKBrso = Hist::ProjectAll(HistProj::kY, hCKBrso);
    std::vector<Hist*> vhKFBrso = Hist::ProjectAll(HistProj::kY, hKFBrso);
    std::vector<Hist*> vhHCBrso = Hist::ProjectAll(HistProj::kY, hHCBrso);

    for (int it = 1; it <= AXbta.nbin(); ++it) {
        // Choutko
        Double_t CKBmax = (*vhCKBrso.at(it))()->GetBinCenter((*vhCKBrso.at(it))()->GetMaximumBin());
        //Double_t CKRmax = 0;
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
        //Double_t KFRmax = 0;
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
        //Double_t HCRmax = 0;
        Double_t HCBrms = (*vhHCBrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, HCBmax, HCBrms);
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

/*
    COUT("resx\n");
    for (Int_t lay = 0; lay < 9; ++lay) {
        Hist* hCKresxS = Hist::New(CSTR("hCKresxL%dS", lay+1), HistAxis(AXmom));
        Hist* hKFresxS = Hist::New(CSTR("hKFresxL%dS", lay+1), HistAxis(AXmom));
        Hist* hHCresxS = Hist::New(CSTR("hHCresxL%dS", lay+1), HistAxis(AXmom));
        
        Hist* hKFCKresxS = Hist::New(CSTR("hKFCKresxL%dS", lay+1), HistAxis(AXmom, "KF/Choutko Sigma Ratio [1]"));
        Hist* hHCCKresxS = Hist::New(CSTR("hHCCKresxL%dS", lay+1), HistAxis(AXmom, "HYChou/Choutko Sigma Ratio [1]"));
        
        std::vector<Hist*> vhCKresx = Hist::ProjectAll(HistProj::kY, hCKresx.at(lay));
        std::vector<Hist*> vhKFresx = Hist::ProjectAll(HistProj::kY, hKFresx.at(lay));
        std::vector<Hist*> vhHCresx = Hist::ProjectAll(HistProj::kY, hHCresx.at(lay));

        for (int it = 1; it <= AXmom.nbin(); ++it) {
            // Choutko
            Double_t CKRmax = (*vhCKresx.at(it))()->GetBinCenter((*vhCKresx.at(it))()->GetMaximumBin());
            Double_t CKRrms = (*vhCKresx.at(it))()->GetRMS();
            gaus->SetParameters(1000, CKRmax, CKRrms);
            (*vhCKresx.at(it))()->Fit(gaus, "q0", "");
            (*vhCKresx.at(it))()->Fit(gaus, "q0", "", CKRmax-stable*CKRrms, CKRmax+stable*CKRrms);
            (*vhCKresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhCKresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhCKresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hCKresxS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hCKresxS)()->SetBinError  (it, gaus->GetParError(2));
            
            // Kalman Filter
            Double_t KFRmax = (*vhKFresx.at(it))()->GetBinCenter((*vhKFresx.at(it))()->GetMaximumBin());
            Double_t KFRrms = (*vhKFresx.at(it))()->GetRMS();
            gaus->SetParameters(1000, KFRmax, KFRrms);
            (*vhKFresx.at(it))()->Fit(gaus, "q0", "");
            (*vhKFresx.at(it))()->Fit(gaus, "q0", "", KFRmax-stable*KFRrms, KFRmax+stable*KFRrms);
            (*vhKFresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhKFresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhKFresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hKFresxS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hKFresxS)()->SetBinError  (it, gaus->GetParError(2));
            
            // Hsin-Yi Chou
            Double_t HCRmax = (*vhHCresx.at(it))()->GetBinCenter((*vhHCresx.at(it))()->GetMaximumBin());
            Double_t HCRrms = (*vhHCresx.at(it))()->GetRMS();
            gaus->SetParameters(1000, HCRmax, HCRrms);
            (*vhHCresx.at(it))()->Fit(gaus, "q0", "");
            (*vhHCresx.at(it))()->Fit(gaus, "q0", "", HCRmax-stable*HCRrms, HCRmax+stable*HCRrms);
            (*vhHCresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhHCresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhHCresx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hHCresxS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hHCresxS)()->SetBinError  (it, gaus->GetParError(2));
            
            (*hKFCKresxS)()->SetBinContent(it, (*hKFresxS)()->GetBinContent(it)/(*hCKresxS)()->GetBinContent(it));
            (*hHCCKresxS)()->SetBinContent(it, (*hHCresxS)()->GetBinContent(it)/(*hCKresxS)()->GetBinContent(it));
        } 
   
        hKFCKresxS->style(Fill(), Line(kBlue), Marker(kBlue));
        hHCCKresxS->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chresxS = Hist::Collect(CSTR("chresxL%dS", lay+1), HistList({ hKFCKresxS, hHCCKresxS }));
        chresxS->Write();
    }


    COUT("resy\n");
    for (Int_t lay = 0; lay < 9; ++lay) {
        Hist* hCKresyS = Hist::New(CSTR("hCKresyL%dS", lay+1), HistAxis(AXmom));
        Hist* hKFresyS = Hist::New(CSTR("hKFresyL%dS", lay+1), HistAxis(AXmom));
        Hist* hHCresyS = Hist::New(CSTR("hHCresyL%dS", lay+1), HistAxis(AXmom));
        
        Hist* hKFCKresyS = Hist::New(CSTR("hKFCKresyL%dS", lay+1), HistAxis(AXmom, "KF/Choutko Sigma Ratio [1]"));
        Hist* hHCCKresyS = Hist::New(CSTR("hHCCKresyL%dS", lay+1), HistAxis(AXmom, "HYChou/Choutko Sigma Ratio [1]"));
        
        std::vector<Hist*> vhCKresy = Hist::ProjectAll(HistProj::kY, hCKresy.at(lay));
        std::vector<Hist*> vhKFresy = Hist::ProjectAll(HistProj::kY, hKFresy.at(lay));
        std::vector<Hist*> vhHCresy = Hist::ProjectAll(HistProj::kY, hHCresy.at(lay));

        for (int it = 1; it <= AXmom.nbin(); ++it) {
            // Choutko
            Double_t CKRmax = (*vhCKresy.at(it))()->GetBinCenter((*vhCKresy.at(it))()->GetMaximumBin());
            Double_t CKRrms = (*vhCKresy.at(it))()->GetRMS();
            gaus->SetParameters(1000, CKRmax, CKRrms);
            (*vhCKresy.at(it))()->Fit(gaus, "q0", "");
            (*vhCKresy.at(it))()->Fit(gaus, "q0", "", CKRmax-stable*CKRrms, CKRmax+stable*CKRrms);
            (*vhCKresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhCKresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhCKresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hCKresyS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hCKresyS)()->SetBinError  (it, gaus->GetParError(2));
            
            // Kalman Filter
            Double_t KFRmax = (*vhKFresy.at(it))()->GetBinCenter((*vhKFresy.at(it))()->GetMaximumBin());
            Double_t KFRrms = (*vhKFresy.at(it))()->GetRMS();
            gaus->SetParameters(1000, KFRmax, KFRrms);
            (*vhKFresy.at(it))()->Fit(gaus, "q0", "");
            (*vhKFresy.at(it))()->Fit(gaus, "q0", "", KFRmax-stable*KFRrms, KFRmax+stable*KFRrms);
            (*vhKFresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhKFresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhKFresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hKFresyS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hKFresyS)()->SetBinError  (it, gaus->GetParError(2));
            
            // Hsin-Yi Chou
            Double_t HCRmax = (*vhHCresy.at(it))()->GetBinCenter((*vhHCresy.at(it))()->GetMaximumBin());
            Double_t HCRrms = (*vhHCresy.at(it))()->GetRMS();
            gaus->SetParameters(1000, HCRmax, HCRrms);
            (*vhHCresy.at(it))()->Fit(gaus, "q0", "");
            (*vhHCresy.at(it))()->Fit(gaus, "q0", "", HCRmax-stable*HCRrms, HCRmax+stable*HCRrms);
            (*vhHCresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhHCresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhHCresy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hHCresyS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hHCresyS)()->SetBinError  (it, gaus->GetParError(2));
            
            (*hKFCKresyS)()->SetBinContent(it, (*hKFresyS)()->GetBinContent(it)/(*hCKresyS)()->GetBinContent(it));
            (*hHCCKresyS)()->SetBinContent(it, (*hHCresyS)()->GetBinContent(it)/(*hCKresyS)()->GetBinContent(it));
        } 
   
        hKFCKresyS->style(Fill(), Line(kBlue), Marker(kBlue));
        hHCCKresyS->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chresyS = Hist::Collect(CSTR("chresyL%dS", lay+1), HistList({ hKFCKresyS, hHCCKresyS }));
        chresyS->Write();
    }


    COUT("cosx\n");
    for (Int_t lay = 0; lay < 9; ++lay) {
        Hist* hCKcosxS = Hist::New(CSTR("hCKcosxL%dS", lay+1), HistAxis(AXmom));
        Hist* hKFcosxS = Hist::New(CSTR("hKFcosxL%dS", lay+1), HistAxis(AXmom));
        Hist* hHCcosxS = Hist::New(CSTR("hHCcosxL%dS", lay+1), HistAxis(AXmom));
        
        Hist* hKFCKcosxS = Hist::New(CSTR("hKFCKcosxL%dS", lay+1), HistAxis(AXmom, "KF/Choutko Sigma Ratio [1]"));
        Hist* hHCCKcosxS = Hist::New(CSTR("hHCCKcosxL%dS", lay+1), HistAxis(AXmom, "HYChou/Choutko Sigma Ratio [1]"));
        
        std::vector<Hist*> vhCKcosx = Hist::ProjectAll(HistProj::kY, hCKcosx.at(lay));
        std::vector<Hist*> vhKFcosx = Hist::ProjectAll(HistProj::kY, hKFcosx.at(lay));
        std::vector<Hist*> vhHCcosx = Hist::ProjectAll(HistProj::kY, hHCcosx.at(lay));

        for (int it = 1; it <= AXmom.nbin(); ++it) {
            // Choutko
            Double_t CKRmax = (*vhCKcosx.at(it))()->GetBinCenter((*vhCKcosx.at(it))()->GetMaximumBin());
            Double_t CKRrms = (*vhCKcosx.at(it))()->GetRMS();
            gaus->SetParameters(1000, CKRmax, CKRrms);
            (*vhCKcosx.at(it))()->Fit(gaus, "q0", "");
            (*vhCKcosx.at(it))()->Fit(gaus, "q0", "", CKRmax-stable*CKRrms, CKRmax+stable*CKRrms);
            (*vhCKcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhCKcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhCKcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hCKcosxS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hCKcosxS)()->SetBinError  (it, gaus->GetParError(2));
            
            // Kalman Filter
            Double_t KFRmax = (*vhKFcosx.at(it))()->GetBinCenter((*vhKFcosx.at(it))()->GetMaximumBin());
            Double_t KFRrms = (*vhKFcosx.at(it))()->GetRMS();
            gaus->SetParameters(1000, KFRmax, KFRrms);
            (*vhKFcosx.at(it))()->Fit(gaus, "q0", "");
            (*vhKFcosx.at(it))()->Fit(gaus, "q0", "", KFRmax-stable*KFRrms, KFRmax+stable*KFRrms);
            (*vhKFcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhKFcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhKFcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hKFcosxS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hKFcosxS)()->SetBinError  (it, gaus->GetParError(2));
            
            // Hsin-Yi Chou
            Double_t HCRmax = (*vhHCcosx.at(it))()->GetBinCenter((*vhHCcosx.at(it))()->GetMaximumBin());
            Double_t HCRrms = (*vhHCcosx.at(it))()->GetRMS();
            gaus->SetParameters(1000, HCRmax, HCRrms);
            (*vhHCcosx.at(it))()->Fit(gaus, "q0", "");
            (*vhHCcosx.at(it))()->Fit(gaus, "q0", "", HCRmax-stable*HCRrms, HCRmax+stable*HCRrms);
            (*vhHCcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhHCcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhHCcosx.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hHCcosxS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hHCcosxS)()->SetBinError  (it, gaus->GetParError(2));
            
            (*hKFCKcosxS)()->SetBinContent(it, (*hKFcosxS)()->GetBinContent(it)/(*hCKcosxS)()->GetBinContent(it));
            (*hHCCKcosxS)()->SetBinContent(it, (*hHCcosxS)()->GetBinContent(it)/(*hCKcosxS)()->GetBinContent(it));
        } 
   
        hKFCKcosxS->style(Fill(), Line(kBlue), Marker(kBlue));
        hHCCKcosxS->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chcosxS = Hist::Collect(CSTR("chcosxL%dS", lay+1), HistList({ hKFCKcosxS, hHCCKcosxS }));
        chcosxS->Write();
    }


    COUT("cosy\n");
    for (Int_t lay = 0; lay < 9; ++lay) {
        Hist* hCKcosyS = Hist::New(CSTR("hCKcosyL%dS", lay+1), HistAxis(AXmom));
        Hist* hKFcosyS = Hist::New(CSTR("hKFcosyL%dS", lay+1), HistAxis(AXmom));
        Hist* hHCcosyS = Hist::New(CSTR("hHCcosyL%dS", lay+1), HistAxis(AXmom));
        
        Hist* hKFCKcosyS = Hist::New(CSTR("hKFCKcosyL%dS", lay+1), HistAxis(AXmom, "KF/Choutko Sigma Ratio [1]"));
        Hist* hHCCKcosyS = Hist::New(CSTR("hHCCKcosyL%dS", lay+1), HistAxis(AXmom, "HYChou/Choutko Sigma Ratio [1]"));
        
        std::vector<Hist*> vhCKcosy = Hist::ProjectAll(HistProj::kY, hCKcosy.at(lay));
        std::vector<Hist*> vhKFcosy = Hist::ProjectAll(HistProj::kY, hKFcosy.at(lay));
        std::vector<Hist*> vhHCcosy = Hist::ProjectAll(HistProj::kY, hHCcosy.at(lay));

        for (int it = 1; it <= AXmom.nbin(); ++it) {
            // Choutko
            Double_t CKRmax = (*vhCKcosy.at(it))()->GetBinCenter((*vhCKcosy.at(it))()->GetMaximumBin());
            Double_t CKRrms = (*vhCKcosy.at(it))()->GetRMS();
            gaus->SetParameters(1000, CKRmax, CKRrms);
            (*vhCKcosy.at(it))()->Fit(gaus, "q0", "");
            (*vhCKcosy.at(it))()->Fit(gaus, "q0", "", CKRmax-stable*CKRrms, CKRmax+stable*CKRrms);
            (*vhCKcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhCKcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhCKcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hCKcosyS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hCKcosyS)()->SetBinError  (it, gaus->GetParError(2));
            
            // Kalman Filter
            Double_t KFRmax = (*vhKFcosy.at(it))()->GetBinCenter((*vhKFcosy.at(it))()->GetMaximumBin());
            Double_t KFRrms = (*vhKFcosy.at(it))()->GetRMS();
            gaus->SetParameters(1000, KFRmax, KFRrms);
            (*vhKFcosy.at(it))()->Fit(gaus, "q0", "");
            (*vhKFcosy.at(it))()->Fit(gaus, "q0", "", KFRmax-stable*KFRrms, KFRmax+stable*KFRrms);
            (*vhKFcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhKFcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhKFcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hKFcosyS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hKFcosyS)()->SetBinError  (it, gaus->GetParError(2));
            
            // Hsin-Yi Chou
            Double_t HCRmax = (*vhHCcosy.at(it))()->GetBinCenter((*vhHCcosy.at(it))()->GetMaximumBin());
            Double_t HCRrms = (*vhHCcosy.at(it))()->GetRMS();
            gaus->SetParameters(1000, HCRmax, HCRrms);
            (*vhHCcosy.at(it))()->Fit(gaus, "q0", "");
            (*vhHCcosy.at(it))()->Fit(gaus, "q0", "", HCRmax-stable*HCRrms, HCRmax+stable*HCRrms);
            (*vhHCcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhHCcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhHCcosy.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        
            (*hHCcosyS)()->SetBinContent(it, gaus->GetParameter(2));
            (*hHCcosyS)()->SetBinError  (it, gaus->GetParError(2));
            
            (*hKFCKcosyS)()->SetBinContent(it, (*hKFcosyS)()->GetBinContent(it)/(*hCKcosyS)()->GetBinContent(it));
            (*hHCCKcosyS)()->SetBinContent(it, (*hHCcosyS)()->GetBinContent(it)/(*hCKcosyS)()->GetBinContent(it));
        } 
   
        hKFCKcosyS->style(Fill(), Line(kBlue), Marker(kBlue));
        hHCCKcosyS->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chcosyS = Hist::Collect(CSTR("chcosyL%dS", lay+1), HistList({ hKFCKcosyS, hHCCKcosyS }));
        chcosyS->Write();
    }
*/
    ofle->Write();
    ofle->Close();

    return 0;
}
