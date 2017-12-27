//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__

#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
   
    Hist::Load("fit_ams02_fill.root", "dat");

    // Fit
    Hist * hCKRrso = Hist::Head("hCKRrso");
    Hist * hCNRrso = Hist::Head("hCNRrso");
    Hist * hHCRrso = Hist::Head("hHCRrso");
    
    Hist * hCKRrsoCut = Hist::Head("hCKRrsoCut");
    Hist * hCNRrsoCut = Hist::Head("hCNRrsoCut");
    Hist * hHCRrsoCut = Hist::Head("hHCRrsoCut");
    
    const Axis& AXmom = hCKRrso->xaxis();
    
    TFile * ofle = new TFile("fit_ams02_fit.root", "RECREATE");
    ofle->cd();
    
    Hist::AddDirectory();

    TGraphErrors* gCKRrso_men = new TGraphErrors(); gCKRrso_men->SetNameTitle("gCKRrso_men", "");
    TGraphErrors* gCKRrso_sgm = new TGraphErrors(); gCKRrso_sgm->SetNameTitle("gCKRrso_sgm", "");
    TGraphErrors* gCNRrso_men = new TGraphErrors(); gCNRrso_men->SetNameTitle("gCNRrso_men", "");
    TGraphErrors* gCNRrso_sgm = new TGraphErrors(); gCNRrso_sgm->SetNameTitle("gCNRrso_sgm", "");
    TGraphErrors* gHCRrso_men = new TGraphErrors(); gHCRrso_men->SetNameTitle("gHCRrso_men", "");
    TGraphErrors* gHCRrso_sgm = new TGraphErrors(); gHCRrso_sgm->SetNameTitle("gHCRrso_sgm", "");
    
    TGraphErrors* gHCCKRrso_sgm = new TGraphErrors(); gHCCKRrso_sgm->SetNameTitle("gHCCKRrso_sgm", "");
    gHCCKRrso_sgm->GetXaxis()->SetTitle("Momentum [GeV]");
    gHCCKRrso_sgm->GetYaxis()->SetTitle("HYChou/Choutko Sigma Ratio [1]");
    
    TGraphErrors* gHCCNRrso_sgm = new TGraphErrors(); gHCCNRrso_sgm->SetNameTitle("gHCCNRrso_sgm", "");
    gHCCNRrso_sgm->GetXaxis()->SetTitle("Momentum [GeV]");
    gHCCNRrso_sgm->GetYaxis()->SetTitle("HYChou/CN Sigma Ratio [1]");
   
    const Double_t stable = 1.5;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    std::vector<Hist*> vhCKRrso = Hist::ProjectAll(HistProj::kY, hCKRrso);
    std::vector<Hist*> vhCNRrso = Hist::ProjectAll(HistProj::kY, hCNRrso);
    std::vector<Hist*> vhHCRrso = Hist::ProjectAll(HistProj::kY, hHCRrso);
    //std::vector<Hist*> vhCKRrso = Hist::ProjectAll(HistProj::kY, hCKRrsoCut);
    //std::vector<Hist*> vhCNRrso = Hist::ProjectAll(HistProj::kY, hCNRrsoCut);
    //std::vector<Hist*> vhHCRrso = Hist::ProjectAll(HistProj::kY, hHCRrsoCut);
    for (int it = 1; it <= AXmom.nbin(); ++it) {
        double mom = AXmom.center(it+1, AxisScale::kLog);
        PhySt part(PartType::Proton);
        part.set_mom(mom);
        Double_t val = part.mom();
        
        Double_t CKRmax = (*vhCKRrso.at(it))()->GetBinCenter((*vhCKRrso.at(it))()->GetMaximumBin());
        Double_t CKRrms = (*vhCKRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, CKRmax, CKRrms);
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", CKRmax-stable*CKRrms, CKRmax+stable*CKRrms);
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        gCKRrso_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gCKRrso_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gCKRrso_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gCKRrso_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        Double_t CNRmax = (*vhCNRrso.at(it))()->GetBinCenter((*vhCNRrso.at(it))()->GetMaximumBin());
        Double_t CNRrms = (*vhCNRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, CNRmax, CNRrms);
        (*vhCNRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhCNRrso.at(it))()->Fit(gaus, "q0", "", CNRmax-stable*CNRrms, CNRmax+stable*CNRrms);
        (*vhCNRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCNRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCNRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        gCNRrso_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gCNRrso_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gCNRrso_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gCNRrso_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        Double_t HCRmax = (*vhHCRrso.at(it))()->GetBinCenter((*vhHCRrso.at(it))()->GetMaximumBin());
        Double_t HCRrms = (*vhHCRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, HCRmax, HCRrms);
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", HCRmax-stable*HCRrms, HCRmax+stable*HCRrms);
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        gHCRrso_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gHCRrso_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gHCRrso_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gHCRrso_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gHCCKRrso_sgm->SetPoint     (it-1, val, gHCRrso_sgm->GetY()[it-1]/gCKRrso_sgm->GetY()[it-1]);
        gHCCKRrso_sgm->SetPointError(it-1,  0., 0.);
        
        gHCCNRrso_sgm->SetPoint     (it-1, val, gHCRrso_sgm->GetY()[it-1]/gCNRrso_sgm->GetY()[it-1]);
        gHCCNRrso_sgm->SetPointError(it-1,  0., 0.);
    } 

    gCKRrso_men->Write();
    gCKRrso_sgm->Write();
    gCNRrso_men->Write();
    gCNRrso_sgm->Write();
    gHCRrso_men->Write();
    gHCRrso_sgm->Write();

    gHCCKRrso_sgm->Write();
    gHCCNRrso_sgm->Write();

    ofle->Write();
    ofle->Close();

    return 0;
}
