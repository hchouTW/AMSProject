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
    //Hist * hCNRrso = Hist::Head("hCNRrso");
    Hist * hKFRrso = Hist::Head("hKFRrso");
    Hist * hHCRrso = Hist::Head("hHCRrso");
    
    Hist * hCKPflux = Hist::Head("hCKPflux");
    //Hist * hCNPflux = Hist::Head("hCNPflux");
    Hist * hKFPflux = Hist::Head("hKFPflux");
    Hist * hHCPflux = Hist::Head("hHCPflux");
    
    Hist * hCKNflux = Hist::Head("hCKNflux");
    //Hist * hCNNflux = Hist::Head("hCNNflux");
    Hist * hKFNflux = Hist::Head("hKFNflux");
    Hist * hHCNflux = Hist::Head("hHCNflux");
    
    Hist * hCKRflux = Hist::Head("hCKRflux");
    //Hist * hCNRflux = Hist::Head("hCNRflux");
    Hist * hKFRflux = Hist::Head("hKFRflux");
    Hist * hHCRflux = Hist::Head("hHCRflux");
    
    Hist * hCKMCflux = Hist::Head("hCKMCflux");
    //Hist * hCNMCflux = Hist::Head("hCNMCflux");
    Hist * hKFMCflux = Hist::Head("hKFMCflux");
    Hist * hHCMCflux = Hist::Head("hHCMCflux");
    
    Hist * hCKRfluxCut = Hist::Head("hCKRfluxCut");
    //Hist * hCNRfluxCut = Hist::Head("hCNRfluxCut");
    Hist * hKFRfluxCut = Hist::Head("hKFRfluxCut");
    Hist * hHCRfluxCut = Hist::Head("hHCRfluxCut");
    
    Hist * hCKMCfluxCut = Hist::Head("hCKMCfluxCut");
    //Hist * hCNMCfluxCut = Hist::Head("hCNMCfluxCut");
    Hist * hKFMCfluxCut = Hist::Head("hKFMCfluxCut");
    Hist * hHCMCfluxCut = Hist::Head("hHCMCfluxCut");
    
    const Axis& AXmom = hCKRrso->xaxis();
    
    TFile * ofle = new TFile("fit_ams02_fit.root", "RECREATE");
    ofle->cd();
    
    Hist::AddDirectory();

    Hist* hCKRrsoM = Hist::New("hCKRrsoM", HistAxis(AXmom));
    Hist* hCKRrsoS = Hist::New("hCKRrsoS", HistAxis(AXmom));
    //Hist* hCNRrsoM = Hist::New("hCNRrsoM", HistAxis(AXmom));
    //Hist* hCNRrsoS = Hist::New("hCNRrsoS", HistAxis(AXmom));
    Hist* hKFRrsoM = Hist::New("hKFRrsoM", HistAxis(AXmom));
    Hist* hKFRrsoS = Hist::New("hKFRrsoS", HistAxis(AXmom));
    Hist* hHCRrsoM = Hist::New("hHCRrsoM", HistAxis(AXmom));
    Hist* hHCRrsoS = Hist::New("hHCRrsoS", HistAxis(AXmom));
    
    //Hist* hCNCKRrsoS = Hist::New("hCNCKRrsoS", HistAxis(AXmom, "CN/Choutko Sigma Ratio [1]"));
    Hist* hKFCKRrsoS = Hist::New("hKFCKRrsoS", HistAxis(AXmom, "KF/Choutko Sigma Ratio [1]"));
    Hist* hHCCKRrsoS = Hist::New("hHCCKRrsoS", HistAxis(AXmom, "HYChou/Choutko Sigma Ratio [1]"));
   
    const Double_t stable = 1.5;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    std::vector<Hist*> vhCKRrso = Hist::ProjectAll(HistProj::kY, hCKRrso);
    //std::vector<Hist*> vhCNRrso = Hist::ProjectAll(HistProj::kY, hCNRrso);
    std::vector<Hist*> vhKFRrso = Hist::ProjectAll(HistProj::kY, hKFRrso);
    std::vector<Hist*> vhHCRrso = Hist::ProjectAll(HistProj::kY, hHCRrso);

    for (int it = 1; it <= AXmom.nbin(); ++it) {
        double mom = AXmom.center(it, AxisScale::kLog);
        PhySt part(PartType::Proton);
        part.set_mom(mom);
        Double_t val = part.mom();
       
        // Choutko
        Double_t CKRmax = (*vhCKRrso.at(it))()->GetBinCenter((*vhCKRrso.at(it))()->GetMaximumBin());
        Double_t CKRrms = (*vhCKRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, CKRmax, CKRrms);
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", CKRmax-stable*CKRrms, CKRmax+stable*CKRrms);
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hCKRrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hCKRrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hCKRrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hCKRrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        // Chikanian
        //Double_t CNRmax = (*vhCNRrso.at(it))()->GetBinCenter((*vhCNRrso.at(it))()->GetMaximumBin());
        //Double_t CNRrms = (*vhCNRrso.at(it))()->GetRMS();
        //gaus->SetParameters(1000, CNRmax, CNRrms);
        //(*vhCNRrso.at(it))()->Fit(gaus, "q0", "");
        //(*vhCNRrso.at(it))()->Fit(gaus, "q0", "", CNRmax-stable*CNRrms, CNRmax+stable*CNRrms);
        //(*vhCNRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        //(*vhCNRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        //(*vhCNRrso.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        //(*hCNRrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        //(*hCNRrsoM)()->SetBinError  (it, gaus->GetParError(1));
        //(*hCNRrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        //(*hCNRrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        // Kalman Filter
        Double_t KFRmax = (*vhKFRrso.at(it))()->GetBinCenter((*vhKFRrso.at(it))()->GetMaximumBin());
        Double_t KFRrms = (*vhKFRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, KFRmax, KFRrms);
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", KFRmax-stable*KFRrms, KFRmax+stable*KFRrms);
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
    
        (*hHCRrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hHCRrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hHCRrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hHCRrsoS)()->SetBinError  (it, gaus->GetParError(2));
        
        //(*hCNCKRrsoS)()->SetBinContent(it, (*hCNRrsoS)()->GetBinContent(it)/(*hCKRrsoS)()->GetBinContent(it));
        (*hKFCKRrsoS)()->SetBinContent(it, (*hKFRrsoS)()->GetBinContent(it)/(*hCKRrsoS)()->GetBinContent(it));
        (*hHCCKRrsoS)()->SetBinContent(it, (*hHCRrsoS)()->GetBinContent(it)/(*hCKRrsoS)()->GetBinContent(it));
    } 
   
    hCKRrsoM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFRrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCRrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrsoM = Hist::Collect("chRrsoM", HistList({ hCKRrsoM, hKFRrsoM, hHCRrsoM }));
    chRrsoM->Write();
    
    hKFCKRrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCCKRrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrsoS = Hist::Collect("chRrsoS", HistList({ hKFCKRrsoS, hHCCKRrsoS }));
    chRrsoS->Write();

    const Axis& AXrig = hCKPflux->xaxis();
    
    std::vector<Hist*> vhCKPflux = Hist::ProjectAll(HistProj::kY, hCKPflux);
    std::vector<Hist*> vhKFPflux = Hist::ProjectAll(HistProj::kY, hKFPflux);
    std::vector<Hist*> vhHCPflux = Hist::ProjectAll(HistProj::kY, hHCPflux);
    
    std::vector<Hist*> vhCKNflux = Hist::ProjectAll(HistProj::kY, hCKNflux);
    std::vector<Hist*> vhKFNflux = Hist::ProjectAll(HistProj::kY, hKFNflux);
    std::vector<Hist*> vhHCNflux = Hist::ProjectAll(HistProj::kY, hHCNflux);
    
    Hist* hCKrat = Hist::New("hCKrat", AXrig);
    Hist* hKFrat = Hist::New("hKFrat", AXrig);
    Hist* hHCrat = Hist::New("hHCrat", AXrig);
    
    Hist* hCKratc = Hist::New("hCKratc", AXrig);
    Hist* hKFratc = Hist::New("hKFratc", AXrig);
    Hist* hHCratc = Hist::New("hHCratc", AXrig);
    
    Hist* hCKcc = Hist::New("hCKcc", AXrig);
    Hist* hKFcc = Hist::New("hKFcc", AXrig);
    Hist* hHCcc = Hist::New("hHCcc", AXrig);
    
    Hist* hCKcc99 = Hist::New("hCKcc99", AXrig);
    Hist* hCKcc97 = Hist::New("hCKcc97", AXrig);
    Hist* hCKcc95 = Hist::New("hCKcc95", AXrig);
    Hist* hKFcc99 = Hist::New("hKFcc99", AXrig);
    Hist* hKFcc97 = Hist::New("hKFcc97", AXrig);
    Hist* hKFcc95 = Hist::New("hKFcc95", AXrig);
    Hist* hHCcc99 = Hist::New("hHCcc99", AXrig);
    Hist* hHCcc97 = Hist::New("hHCcc97", AXrig);
    Hist* hHCcc95 = Hist::New("hHCcc95", AXrig);
    
    Hist* hCKnm = Hist::New("hCKnm", AXrig);
    Hist* hKFnm = Hist::New("hKFnm", AXrig);
    Hist* hHCnm = Hist::New("hHCnm", AXrig);
    
    Hist* hCKnm99 = Hist::New("hCKnm99", AXrig);
    Hist* hCKnm97 = Hist::New("hCKnm97", AXrig);
    Hist* hCKnm95 = Hist::New("hCKnm95", AXrig);
    Hist* hKFnm99 = Hist::New("hKFnm99", AXrig);
    Hist* hKFnm97 = Hist::New("hKFnm97", AXrig);
    Hist* hKFnm95 = Hist::New("hKFnm95", AXrig);
    Hist* hHCnm99 = Hist::New("hHCnm99", AXrig);
    Hist* hHCnm97 = Hist::New("hHCnm97", AXrig);
    Hist* hHCnm95 = Hist::New("hHCnm95", AXrig);

    Hist* hKFCKnm = Hist::New("hKFCKnm", AXrig);
    Hist* hHCCKnm = Hist::New("hHCCKnm", AXrig);
    
    Hist* hKFCKnm99 = Hist::New("hKFCKnm99", AXrig);
    Hist* hKFCKnm97 = Hist::New("hKFCKnm97", AXrig);
    Hist* hKFCKnm95 = Hist::New("hKFCKnm95", AXrig);
    Hist* hHCCKnm99 = Hist::New("hHCCKnm99", AXrig);
    Hist* hHCCKnm97 = Hist::New("hHCCKnm97", AXrig);
    Hist* hHCCKnm95 = Hist::New("hHCCKnm95", AXrig);

    std::vector<Int_t> vcutcc({ 99, 97, 94, 90, 85, 80, 75, 70 });
    std::vector<Hist*> vhCKvl(vcutcc.size(), nullptr);
    std::vector<Hist*> vhKFvl(vcutcc.size(), nullptr);
    std::vector<Hist*> vhHCvl(vcutcc.size(), nullptr);
    std::vector<Hist*> vhHCcc(vcutcc.size(), nullptr);
    for (Int_t it = 0; it < vcutcc.size(); ++it) {
        vhCKvl.at(it) = Hist::New(Form("hCKvlEFF%02d", vcutcc.at(it)), AXrig);
        vhKFvl.at(it) = Hist::New(Form("hKFvlEFF%02d", vcutcc.at(it)), AXrig);
        vhHCvl.at(it) = Hist::New(Form("hHCvlEFF%02d", vcutcc.at(it)), AXrig);
        vhHCcc.at(it) = Hist::New(Form("hHCccEFF%02d", vcutcc.at(it)), AXrig);
    }

    for (int it = 1; it <= AXrig.nbin(); ++it) {
        double rig = AXrig.center(it, AxisScale::kLog);
        if (rig < 60.) continue;

        Double_t CKnum = (*vhCKPflux.at(it))()->Integral(0, (*vhCKPflux.at(it))()->GetNbinsX()+1);
        TGraph grCKcc, grCKnm, grCKvl;
        for (Int_t ib = 1; ib <= (*vhCKPflux.at(it))()->GetNbinsX()+1; ++ib) {
            Double_t selP = (*vhCKPflux.at(it))()->Integral(0, ib);
            Double_t selN = (*vhCKNflux.at(it))()->Integral(0, ib);
            Double_t rtNP = (MGNumc::Compare(selP) > 0) ? (selN/selP) : 0;
            grCKcc.SetPoint(ib-1, selP/CKnum, rtNP);
            grCKnm.SetPoint(ib-1, selP/CKnum, selN);
            grCKvl.SetPoint(ib-1, selP/CKnum, (*vhCKPflux.at(it))()->GetXaxis()->GetBinCenter(ib));
        }
        (*hCKrat )()->SetBinContent(it, (*hCKRflux)()->GetBinContent(it)/(*hCKMCflux)()->GetBinContent(it));
        (*hCKratc)()->SetBinContent(it, (*hCKRfluxCut)()->GetBinContent(it)/(*hCKMCfluxCut)()->GetBinContent(it));
        (*hCKcc  )()->SetBinContent(it, grCKcc.Eval(1.00));
        (*hCKcc99)()->SetBinContent(it, grCKcc.Eval(0.99));
        (*hCKcc97)()->SetBinContent(it, grCKcc.Eval(0.97));
        (*hCKcc95)()->SetBinContent(it, grCKcc.Eval(0.95));
        (*hCKnm  )()->SetBinContent(it, grCKnm.Eval(1.00));
        (*hCKnm99)()->SetBinContent(it, grCKnm.Eval(0.99));
        (*hCKnm97)()->SetBinContent(it, grCKnm.Eval(0.97));
        (*hCKnm95)()->SetBinContent(it, grCKnm.Eval(0.95));
        
        Double_t KFnum = (*vhKFPflux.at(it))()->Integral(0, (*vhKFPflux.at(it))()->GetNbinsX()+1);
        TGraph grKFcc, grKFnm, grKFvl;
        for (Int_t ib = 1; ib <= (*vhKFPflux.at(it))()->GetNbinsX()+1; ++ib) {
            Double_t selP = (*vhKFPflux.at(it))()->Integral(0, ib);
            Double_t selN = (*vhKFNflux.at(it))()->Integral(0, ib);
            Double_t rtNP = (MGNumc::Compare(selP) > 0) ? (selN/selP) : 0;
            grKFcc.SetPoint(ib-1, selP/KFnum, rtNP);
            grKFnm.SetPoint(ib-1, selP/KFnum, selN);
            grKFvl.SetPoint(ib-1, selP/KFnum, (*vhKFPflux.at(it))()->GetXaxis()->GetBinCenter(ib));
        }
        (*hKFrat )()->SetBinContent(it, (*hKFRflux)()->GetBinContent(it)/(*hKFMCflux)()->GetBinContent(it));
        (*hKFratc)()->SetBinContent(it, (*hKFRfluxCut)()->GetBinContent(it)/(*hKFMCfluxCut)()->GetBinContent(it));
        (*hKFcc  )()->SetBinContent(it, grKFcc.Eval(1.00));
        (*hKFcc99)()->SetBinContent(it, grKFcc.Eval(0.99));
        (*hKFcc97)()->SetBinContent(it, grKFcc.Eval(0.97));
        (*hKFcc95)()->SetBinContent(it, grKFcc.Eval(0.95));
        (*hKFnm  )()->SetBinContent(it, grKFnm.Eval(1.00));
        (*hKFnm99)()->SetBinContent(it, grKFnm.Eval(0.99));
        (*hKFnm97)()->SetBinContent(it, grKFnm.Eval(0.97));
        (*hKFnm95)()->SetBinContent(it, grKFnm.Eval(0.95));
        
        Double_t HCnum = (*vhHCPflux.at(it))()->Integral(0, (*vhHCPflux.at(it))()->GetNbinsX()+1);
        TGraph grHCcc, grHCnm, grHCvl;
        for (Int_t ib = 1; ib <= (*vhHCPflux.at(it))()->GetNbinsX()+1; ++ib) {
            Double_t selP = (*vhHCPflux.at(it))()->Integral(0, ib);
            Double_t selN = (*vhHCNflux.at(it))()->Integral(0, ib);
            Double_t rtNP = (MGNumc::Compare(selP) > 0) ? (selN/selP) : 0;
            grHCcc.SetPoint(ib-1, selP/HCnum, rtNP);
            grHCnm.SetPoint(ib-1, selP/HCnum, selN);
            grHCvl.SetPoint(ib-1, selP/HCnum, (*vhHCPflux.at(it))()->GetXaxis()->GetBinCenter(ib));
        }
        (*hHCrat )()->SetBinContent(it, (*hHCRflux)()->GetBinContent(it)/(*hHCMCflux)()->GetBinContent(it));
        (*hHCratc)()->SetBinContent(it, (*hHCRfluxCut)()->GetBinContent(it)/(*hHCMCfluxCut)()->GetBinContent(it));
        (*hHCcc  )()->SetBinContent(it, grHCcc.Eval(1.00));
        (*hHCcc99)()->SetBinContent(it, grHCcc.Eval(0.99));
        (*hHCcc97)()->SetBinContent(it, grHCcc.Eval(0.97));
        (*hHCcc95)()->SetBinContent(it, grHCcc.Eval(0.95));
        (*hHCnm  )()->SetBinContent(it, grHCnm.Eval(1.00));
        (*hHCnm99)()->SetBinContent(it, grHCnm.Eval(0.99));
        (*hHCnm97)()->SetBinContent(it, grHCnm.Eval(0.97));
        (*hHCnm95)()->SetBinContent(it, grHCnm.Eval(0.95));
        
        (*hKFCKnm  )()->SetBinContent(it, (*hKFnm  )()->GetBinContent(it)/(*hCKnm  )()->GetBinContent(it));
        (*hKFCKnm99)()->SetBinContent(it, (*hKFnm99)()->GetBinContent(it)/(*hCKnm99)()->GetBinContent(it));
        (*hKFCKnm97)()->SetBinContent(it, (*hKFnm97)()->GetBinContent(it)/(*hCKnm97)()->GetBinContent(it));
        (*hKFCKnm95)()->SetBinContent(it, (*hKFnm95)()->GetBinContent(it)/(*hCKnm95)()->GetBinContent(it));
        (*hHCCKnm  )()->SetBinContent(it, (*hHCnm  )()->GetBinContent(it)/(*hCKnm  )()->GetBinContent(it));
        (*hHCCKnm99)()->SetBinContent(it, (*hHCnm99)()->GetBinContent(it)/(*hCKnm99)()->GetBinContent(it));
        (*hHCCKnm97)()->SetBinContent(it, (*hHCnm97)()->GetBinContent(it)/(*hCKnm97)()->GetBinContent(it));
        (*hHCCKnm95)()->SetBinContent(it, (*hHCnm95)()->GetBinContent(it)/(*hCKnm95)()->GetBinContent(it));

        Double_t refnm = grHCnm.Eval(1.0);
        for (Int_t jt = 0; jt < vcutcc.size(); ++jt) {
            (*vhCKvl.at(jt))()->SetBinContent(it, grCKvl.Eval(vcutcc.at(jt)*0.01));
            (*vhKFvl.at(jt))()->SetBinContent(it, grKFvl.Eval(vcutcc.at(jt)*0.01));
            (*vhHCvl.at(jt))()->SetBinContent(it, grHCvl.Eval(vcutcc.at(jt)*0.01));
            (*vhHCcc.at(jt))()->SetBinContent(it, grHCnm.Eval(vcutcc.at(jt)*0.01)/grHCnm.Eval(1.0));
        }
    }
    
    hCKrat->style(Fill(), Line(kGreen));
    hKFrat->style(Fill(), Line(kBlue));
    hHCrat->style(Fill(), Line(kRed));
    THStack* chrat = Hist::Collect("chrat", HistList({ hCKrat, hKFrat, hHCrat }));
    chrat->Write();
    
    hCKratc->style(Fill(), Line(kGreen));
    hKFratc->style(Fill(), Line(kBlue));
    hHCratc->style(Fill(), Line(kRed));
    THStack* chratc = Hist::Collect("chratc", HistList({ hCKratc, hKFratc, hHCratc }));
    chratc->Write();

    hKFCKnm  ->style(Fill(), Line(kBlue));
    hKFCKnm99->style(Fill(), Line(kBlue));
    hKFCKnm97->style(Fill(), Line(kBlue));
    hKFCKnm95->style(Fill(), Line(kBlue));
    hHCCKnm  ->style(Fill(), Line(kRed));
    hHCCKnm99->style(Fill(), Line(kRed));
    hHCCKnm97->style(Fill(), Line(kRed));
    hHCCKnm95->style(Fill(), Line(kRed));
    THStack* chnm   = Hist::Collect("chnm"  , HistList({ hKFCKnm  , hHCCKnm   }));
    THStack* chnm99 = Hist::Collect("chnm99", HistList({ hKFCKnm99, hHCCKnm99 }));
    THStack* chnm97 = Hist::Collect("chnm97", HistList({ hKFCKnm97, hHCCKnm97 }));
    THStack* chnm95 = Hist::Collect("chnm95", HistList({ hKFCKnm95, hHCCKnm95 }));
    chnm  ->Write();
    chnm99->Write();
    chnm97->Write();
    chnm95->Write();
    
    hCKcc  ->style(Fill(), Line(kGreen));
    hCKcc99->style(Fill(), Line(kGreen));
    hCKcc97->style(Fill(), Line(kGreen));
    hCKcc95->style(Fill(), Line(kGreen));
    hKFcc  ->style(Fill(), Line(kBlue));
    hKFcc99->style(Fill(), Line(kBlue));
    hKFcc97->style(Fill(), Line(kBlue));
    hKFcc95->style(Fill(), Line(kBlue));
    hHCcc  ->style(Fill(), Line(kRed));
    hHCcc99->style(Fill(), Line(kRed));
    hHCcc97->style(Fill(), Line(kRed));
    hHCcc95->style(Fill(), Line(kRed));
    
    THStack* chcc   = Hist::Collect("chcc"  , HistList({ hCKcc  , hKFcc  , hHCcc   }));
    THStack* chcc99 = Hist::Collect("chcc99", HistList({ hCKcc99, hKFcc99, hHCcc99 }));
    THStack* chcc97 = Hist::Collect("chcc97", HistList({ hCKcc97, hKFcc97, hHCcc97 }));
    THStack* chcc95 = Hist::Collect("chcc95", HistList({ hCKcc95, hKFcc95, hHCcc95 }));
    chcc  ->Write();
    chcc99->Write();
    chcc97->Write();
    chcc95->Write();
    
    for (Int_t it = 0; it < vhHCcc.size(); ++it) vhHCcc.at(it)->style(Fill(), Line(it+2));
    THStack* chHCcc = Hist::Collect("chHCcc", vhHCcc);
    chHCcc->Write();
    
    ofle->Write();
    ofle->Close();

    return 0;
}
