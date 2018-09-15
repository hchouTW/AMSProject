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
    Hist* hKFnum   = Hist::Head("hKFnum");
    Hist* hHCnum   = Hist::Head("hHCnum");

    // Tme
    Hist* hMCtme   = Hist::Head("hMCtme");
    Hist* hCKtme   = Hist::Head("hCKtme");
    Hist* hKFtme   = Hist::Head("hKFtme");
    Hist* hHCtme   = Hist::Head("hHCtme");

    // Fit
    Hist* hCKRrso   = Hist::Head("hCKRrso");
    Hist* hKFRrso   = Hist::Head("hKFRrso");
    Hist* hHCRrso   = Hist::Head("hHCRrso");

/*    
    Hist* hHCBrso = Hist::Head("hHCBrso");
    Hist* hTFBrso = Hist::Head("hTFBrso");
    Hist* hRHBrso = Hist::Head("hRHBrso");

    Hist* hCKMrso = Hist::Head("hCKMrso");
    Hist* hKFMrso = Hist::Head("hKFMrso");
    Hist* hHCMrso = Hist::Head("hHCMrso");
 */   
/*
    Hist* hHCMrsoCO = Hist::Head("hHCMrsoCO");
    Hist* hHCMrsoEO = Hist::Head("hHCMrsoEO");
    Hist* hHCMrsoSO = Hist::Head("hHCMrsoSO");
    
    Hist* hHCMrsoC[4] = { nullptr };
    Hist* hHCMrsoE[4] = { nullptr };
    Hist* hHCMrsoS[4] = { nullptr };
    for (Int_t is = 0; is < 4; ++is) {
        hHCMrsoC[is] = Hist::Head(Form("hHCMrsoC%d", is));
        hHCMrsoE[is] = Hist::Head(Form("hHCMrsoE%d", is));
        hHCMrsoS[is] = Hist::Head(Form("hHCMrsoS%d", is));
    }
  */  
    
    Hist* hCKcx = Hist::Head("hCKcx");
    Hist* hKFcx = Hist::Head("hKFcx");
    Hist* hHCcx = Hist::Head("hHCcx");
    
    Hist* hCKcy = Hist::Head("hCKcy");
    Hist* hKFcy = Hist::Head("hKFcy");
    Hist* hHCcy = Hist::Head("hHCcy");
    
    Hist* hCKux = Hist::Head("hCKux");
    Hist* hKFux = Hist::Head("hKFux");
    Hist* hHCux = Hist::Head("hHCux");
    
    Hist* hCKuy = Hist::Head("hCKuy");
    Hist* hKFuy = Hist::Head("hKFuy");
    Hist* hHCuy = Hist::Head("hHCuy");
 
    const Axis& AXrig = hCKRrso->xaxis();
    
    TFile * ofle = new TFile("fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();
    
    COUT("rat\n");
    Hist* hCKrat   = Hist::Calculate(HistArith::kDivide, "hCKrat", ";Momentum [GeV];Effienecy", hCKnum,   hMCnum);
    Hist* hKFrat   = Hist::Calculate(HistArith::kDivide, "hKFrat", ";Momentum [GeV];Effienecy", hKFnum,   hMCnum);
    Hist* hHCrat   = Hist::Calculate(HistArith::kDivide, "hHCrat", ";Momentum [GeV];Effienecy", hHCnum,   hMCnum);
    
    hCKrat  ->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFrat  ->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCrat  ->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRat = Hist::Collect("chrat", HistList({ hCKrat, hKFrat, hHCrat }));
    chRat->Write();
    
    COUT("tme\n");
    Hist* hCKtmeM = Hist::Calculate(HistArith::kDivide, "hCKtmeM", ";Momentum [GeV];Mean Time [msec]", hCKtme,   hMCtme);
    Hist* hKFtmeM = Hist::Calculate(HistArith::kDivide, "hKFtmeM", ";Momentum [GeV];Mean Time [msec]", hKFtme,   hMCtme);
    Hist* hHCtmeM = Hist::Calculate(HistArith::kDivide, "hHCtmeM", ";Momentum [GeV];Mean Time [msec]", hHCtme,   hMCtme);
    
    hCKtmeM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFtmeM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCtmeM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chTmeM = Hist::Collect("chTmeM", HistList({ hCKtmeM, hKFtmeM, hHCtmeM }));
    chTmeM->Write();
    
    const Double_t stable = 2.0;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);

    COUT("CKRrso\n");
    Hist* hCKRrsoM = Hist::New("hCKRrsoM", HistAxis(AXrig, "Mean [1]"));
    Hist* hCKRrsoS = Hist::New("hCKRrsoS", HistAxis(AXrig, "Sigma [1]"));
    std::vector<Hist*> vhCKRrso = Hist::ProjectAll(HistProj::kY, hCKRrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhCKRrso.at(it))()->GetBinCenter((*vhCKRrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhCKRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hCKRrsoM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hCKRrsoM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hCKRrsoS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hCKRrsoS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hCKRrsoM->write();
    hCKRrsoS->write();


    COUT("KFRrso\n");
    Hist* hKFRrsoM = Hist::New("hKFRrsoM", HistAxis(AXrig, "Mean [1]"));
    Hist* hKFRrsoS = Hist::New("hKFRrsoS", HistAxis(AXrig, "Sigma [1]"));
    std::vector<Hist*> vhKFRrso = Hist::ProjectAll(HistProj::kY, hKFRrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhKFRrso.at(it))()->GetBinCenter((*vhKFRrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhKFRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hKFRrsoM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hKFRrsoM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hKFRrsoS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hKFRrsoS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hKFRrsoM->write();
    hKFRrsoS->write();


    COUT("HCRrso\n");
    Hist* hHCRrsoM = Hist::New("hHCRrsoM", HistAxis(AXrig, "Mean [1]"));
    Hist* hHCRrsoS = Hist::New("hHCRrsoS", HistAxis(AXrig, "Sigma [1]"));
    std::vector<Hist*> vhHCRrso = Hist::ProjectAll(HistProj::kY, hHCRrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhHCRrso.at(it))()->GetBinCenter((*vhHCRrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhHCRrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "");
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCRrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hHCRrsoM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hHCRrsoM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hHCRrsoS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hHCRrsoS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hHCRrsoM->write();
    hHCRrsoS->write();
    
    hCKRrsoM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFRrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCRrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrsoM = Hist::Collect("chRrsoM", HistList({ hCKRrsoM, hKFRrsoM, hHCRrsoM }));
    chRrsoM->Write();
    
    hCKRrsoS->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFRrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCRrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrsoS = Hist::Collect("chRrsoS", HistList({ hCKRrsoS, hKFRrsoS, hHCRrsoS }));
    chRrsoS->Write();
    

    Hist* hKFCKRrso = Hist::New("hKFCKRrso", HistAxis(AXrig, "KF/Choutko Sigma Ratio [1]"));
    Hist* hHCCKRrso = Hist::New("hHCCKRrso", HistAxis(AXrig, "HYChou/Choutko Sigma Ratio [1]"));
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t cen = AXrig.center(it, AxisScale::kLog);
        hKFCKRrso->fillH1D(cen, (*hKFRrsoS)()->GetBinContent(it) / (*hCKRrsoS)()->GetBinContent(it));
        hHCCKRrso->fillH1D(cen, (*hHCRrsoS)()->GetBinContent(it) / (*hCKRrsoS)()->GetBinContent(it));
        (*hKFCKRrso)()->SetBinError(it, 0);
        (*hHCCKRrso)()->SetBinError(it, 0);
    }
    
    hKFCKRrso  ->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCCKRrso  ->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRrso = Hist::Collect("chRrso", HistList({ hKFCKRrso, hHCCKRrso }));
    chRrso->Write();
    
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhCKRrso.at(it)->style(Fill(), Line(kGreen), Marker(kGreen));
        vhKFRrso.at(it)->style(Fill(), Line(kBlue), Marker(kBlue));
        vhHCRrso.at(it)->style(Fill(), Line(kRed), Marker(kRed));
        THStack* cvhRrso = Hist::Collect(Form("cvhRrso%03d", it), HistList({ vhCKRrso.at(it), vhKFRrso.at(it), vhHCRrso.at(it) }));
        cvhRrso->Write();
    }

  /* 
    COUT("CKMrso\n");
    Hist* hCKMrsoM = Hist::New("hCKMrsoM", HistAxis(AXrig, "Mean [GeV]"));
    Hist* hCKMrsoS = Hist::New("hCKMrsoS", HistAxis(AXrig, "Sigma [GeV]"));
    std::vector<Hist*> vhCKMrso = Hist::ProjectAll(HistProj::kY, hCKMrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t max = (*vhCKMrso.at(it))()->GetBinCenter((*vhCKMrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhCKMrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhCKMrso.at(it))()->Fit(gaus, "q0", "");
        (*vhCKMrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhCKMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        //(*vhCKMrso.at(it))()->Scale(1.0/(*vhCKMrso.at(it))()->GetEntries());
    
        (*hCKMrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hCKMrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hCKMrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hCKMrsoS)()->SetBinError  (it, gaus->GetParError(2));
    } 
    hCKMrsoM->write();
    hCKMrsoS->write();


    COUT("KFMrso\n");
    Hist* hKFMrsoM = Hist::New("hKFMrsoM", HistAxis(AXrig, "Mean [GeV]"));
    Hist* hKFMrsoS = Hist::New("hKFMrsoS", HistAxis(AXrig, "Sigma [GeV]"));
    std::vector<Hist*> vhKFMrso = Hist::ProjectAll(HistProj::kY, hKFMrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t max = (*vhKFMrso.at(it))()->GetBinCenter((*vhKFMrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhKFMrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhKFMrso.at(it))()->Fit(gaus, "q0", "");
        (*vhKFMrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhKFMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFMrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        //(*vhKFMrso.at(it))()->Scale(1.0/(*vhKFMrso.at(it))()->GetEntries());
    
        (*hKFMrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hKFMrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hKFMrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hKFMrsoS)()->SetBinError  (it, gaus->GetParError(2));
    } 
    hKFMrsoM->write();
    hKFMrsoS->write();


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
    

    hCKMrsoM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFMrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCMrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMrsoM = Hist::Collect("chMrsoM", HistList({ hCKMrsoM, hKFMrsoM, hHCMrsoM }));
    chMrsoM->Write();
    
    hCKMrsoS->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFMrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCMrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMrsoS = Hist::Collect("chMrsoS", HistList({ hCKMrsoS, hKFMrsoS, hHCMrsoS }));
    chMrsoS->Write();
    
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhCKMrso.at(it)->style(Fill(), Line(kGreen), Marker(kGreen));
        vhKFMrso.at(it)->style(Fill(), Line(kBlue), Marker(kBlue));
        vhHCMrso.at(it)->style(Fill(), Line(kRed), Marker(kRed));
        THStack* cvhMrso = Hist::Collect(Form("cvhMrso%03d", it), HistList({ vhCKMrso.at(it), vhKFMrso.at(it), vhHCMrso.at(it) }));
        cvhMrso->Write();
    }
    */
   
/*
    Hist* hHCMrsoSM[4] = { nullptr };
    Hist* hHCMrsoSS[4] = { nullptr };
    std::vector<Hist*> vhHCMrsoS[4];
    for (Int_t is = 0; is < 4; ++is) {
        COUT("HCMrsoS%d\n", is);
        hHCMrsoSM[is] = Hist::New(Form("hHCMrsoSM%d", is), HistAxis(AXrig, "Mean [GeV]"));
        hHCMrsoSS[is] = Hist::New(Form("hHCMrsoSS%d", is), HistAxis(AXrig, "Sigma [GeV]"));
        vhHCMrsoS[is] = Hist::ProjectAll(HistProj::kY, hHCMrsoS[is]);
        for (int it = 1; it <= AXrig.nbin(); ++it) {
            Double_t max = (*(vhHCMrsoS[is]).at(it))()->GetBinCenter((*(vhHCMrsoS[is]).at(it))()->GetMaximumBin());
            Double_t rms = 0.5 * (*(vhHCMrsoS[is]).at(it))()->GetRMS();
            gaus->SetParameters(1000, max, rms);
            (*(vhHCMrsoS[is]).at(it))()->Fit(gaus, "q0", "");
            (*(vhHCMrsoS[is]).at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
            (*(vhHCMrsoS[is]).at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
            (*(vhHCMrsoS[is]).at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
            (*(vhHCMrsoS[is]).at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
            (vhHCMrsoS[is]).at(it)->style(Fill(), Line(kViolet+is), Marker(kViolet+is));
        
            (*(hHCMrsoSM[is]))()->SetBinContent(it, gaus->GetParameter(1));
            (*(hHCMrsoSM[is]))()->SetBinError  (it, gaus->GetParError(1));
            (*(hHCMrsoSS[is]))()->SetBinContent(it, gaus->GetParameter(2));
            (*(hHCMrsoSS[is]))()->SetBinError  (it, gaus->GetParError(2));
        } 
        hHCMrsoSM[is]->write();
        hHCMrsoSS[is]->write();
        hHCMrsoSM[is]->style(Fill(), Line(kViolet+is), Marker(kViolet+is));
        hHCMrsoSS[is]->style(Fill(), Line(kViolet+is), Marker(kViolet+is));
        
        (*(hHCMrsoC[is]))()->Divide((*hHCMrsoCO)());
        hHCMrsoC[is]->style(Fill(), Line(kViolet+is), Marker(kViolet+is));

        (*(hHCMrsoE[is]))()->Divide((*hHCMrsoEO)());
        hHCMrsoE[is]->style(Fill(), Line(kViolet+is), Marker(kViolet+is));
    }
    hHCMrsoEO->style(Fill(), Line(kRed), Marker(kRed));
    
    THStack* chSetMrsoC = Hist::Collect("chSetMrsoC", HistList({ hHCMrsoC[3], hHCMrsoC[2], hHCMrsoC[1] }));
    chSetMrsoC->Write();
    
    THStack* chSetMrsoE = Hist::Collect("chSetMrsoE", HistList({ hHCMrsoE[3], hHCMrsoE[2], hHCMrsoE[1] }));
    chSetMrsoE->Write();
    
    THStack* chSetMrsoM = Hist::Collect("chSetMrsoM", HistList({ hHCMrsoSM[3], hHCMrsoSM[2], hHCMrsoSM[1] }));
    chSetMrsoM->Write();
    
    THStack* chSetMrsoS = Hist::Collect("chSetMrsoS", HistList({ hHCMrsoSS[3], hHCMrsoSS[2], hHCMrsoSS[1] }));
    chSetMrsoS->Write();
    
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        (*vhHCMrso.at(it))()->Scale(1.0/(*vhHCMrso.at(it))()->GetEntries());
        (*vhHCMrsoS[0].at(it))()->Scale(1.0/(*vhHCMrsoS[0].at(it))()->GetEntries());
        (*vhHCMrsoS[1].at(it))()->Scale(1.0/(*vhHCMrsoS[1].at(it))()->GetEntries());
        (*vhHCMrsoS[2].at(it))()->Scale(1.0/(*vhHCMrsoS[2].at(it))()->GetEntries());
        (*vhHCMrsoS[3].at(it))()->Scale(1.0/(*vhHCMrsoS[3].at(it))()->GetEntries());

        THStack* cvhSetMrso = Hist::Collect(Form("cvhSetMrso%03d", it), HistList({ vhHCMrso.at(it), vhHCMrsoS[3].at(it), vhHCMrsoS[2].at(it), vhHCMrsoS[1].at(it) }));
        cvhSetMrso->Write();
    }
*/
    
   /* 
    COUT("HCBrso\n");
    Hist* hHCBrsoM = Hist::New("hHCBrsoM", HistAxis(AXrig, "Mean [1]"));
    Hist* hHCBrsoS = Hist::New("hHCBrsoS", HistAxis(AXrig, "Sigma [1]"));
    std::vector<Hist*> vhHCBrso = Hist::ProjectAll(HistProj::kY, hHCBrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t igb = (0.938272297 / AXrig.center(it, AxisScale::kLog));
        Double_t scl = std::sqrt(igb * igb + 1.0);
        Double_t max = (*vhHCBrso.at(it))()->GetBinCenter((*vhHCBrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhHCBrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "");
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCBrso.at(it))()->Scale(1.0/(*vhHCBrso.at(it))()->GetEntries());
    
        (*hHCBrsoM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hHCBrsoM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hHCBrsoS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hHCBrsoS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hHCBrsoM->write();
    hHCBrsoS->write();

    
    COUT("TFBrso\n");
    Hist* hTFBrsoM = Hist::New("hTFBrsoM", HistAxis(AXrig, "Mean [1]"));
    Hist* hTFBrsoS = Hist::New("hTFBrsoS", HistAxis(AXrig, "Sigma [1]"));
    std::vector<Hist*> vhTFBrso = Hist::ProjectAll(HistProj::kY, hTFBrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t igb = (0.938272297 / AXrig.center(it, AxisScale::kLog));
        Double_t scl = std::sqrt(igb * igb + 1.0);
        Double_t max = (*vhTFBrso.at(it))()->GetBinCenter((*vhTFBrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhTFBrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "");
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhTFBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhTFBrso.at(it))()->Scale(1.0/(*vhTFBrso.at(it))()->GetEntries());
    
        (*hTFBrsoM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hTFBrsoM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hTFBrsoS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hTFBrsoS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hTFBrsoM->write();
    hTFBrsoS->write();

    
    COUT("RHBrso\n");
    Hist* hRHBrsoM = Hist::New("hRHBrsoM", HistAxis(AXrig, "Mean [1]"));
    Hist* hRHBrsoS = Hist::New("hRHBrsoS", HistAxis(AXrig, "Sigma [1]"));
    std::vector<Hist*> vhRHBrso = Hist::ProjectAll(HistProj::kY, hRHBrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t igb = (0.938272297 / AXrig.center(it, AxisScale::kLog));
        Double_t scl = std::sqrt(igb * igb + 1.0);
        Double_t max = (*vhRHBrso.at(it))()->GetBinCenter((*vhRHBrso.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhRHBrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhRHBrso.at(it))()->Fit(gaus, "q0", "");
        (*vhRHBrso.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhRHBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhRHBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhRHBrso.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhRHBrso.at(it))()->Scale(1.0/(*vhRHBrso.at(it))()->GetEntries());
    
        (*hRHBrsoM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hRHBrsoM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hRHBrsoS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hRHBrsoS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hRHBrsoM->write();
    hRHBrsoS->write();

    
    hTFBrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCBrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chTFBrsoM = Hist::Collect("chTFBrsoM", HistList({ hTFBrsoM, hHCBrsoM }));
    chTFBrsoM->Write();
    
    hTFBrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCBrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chTFBrsoS = Hist::Collect("chTFBrsoS", HistList({ hTFBrsoS, hHCBrsoS }));
    chTFBrsoS->Write();
    
    hRHBrsoM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCBrsoM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRHBrsoM = Hist::Collect("chRHBrsoM", HistList({ hRHBrsoM, hHCBrsoM }));
    chRHBrsoM->Write();
    
    hRHBrsoS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCBrsoS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chRHBrsoS = Hist::Collect("chRHBrsoS", HistList({ hRHBrsoS, hHCBrsoS }));
    chRHBrsoS->Write();
    
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhTFBrso  .at(it)->style(Fill(), Line(kBlue), Marker(kBlue));
        vhHCBrso  .at(it)->style(Fill(), Line(kRed), Marker(kRed));
        THStack* cvhTFBrso = Hist::Collect(Form("cvhTFBrso%03d", it), HistList({ vhTFBrso.at(it), vhHCBrso.at(it) }));
        cvhTFBrso->Write();
    }
    
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhRHBrso  .at(it)->style(Fill(), Line(kBlue), Marker(kBlue));
        vhHCBrso  .at(it)->style(Fill(), Line(kRed), Marker(kRed));
        THStack* cvhRHBrso = Hist::Collect(Form("cvhRHBrso%03d", it), HistList({ vhRHBrso.at(it), vhHCBrso.at(it) }));
        cvhRHBrso->Write();
    }
*/    
    
    COUT("CKcx\n");
    Hist* hCKcxM = Hist::New("hCKcxM", HistAxis(AXrig, "Mean"));
    Hist* hCKcxS = Hist::New("hCKcxS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhCKcx = Hist::ProjectAll(HistProj::kY, hCKcx);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhCKcx.at(it))()->GetBinCenter((*vhCKcx.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhCKcx.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhCKcx.at(it))()->Fit(gaus, "q0", "");
        (*vhCKcx.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhCKcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hCKcxM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hCKcxM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hCKcxS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hCKcxS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hCKcxM->write();
    hCKcxS->write();
    

    COUT("KFcx\n");
    Hist* hKFcxM = Hist::New("hKFcxM", HistAxis(AXrig, "Mean"));
    Hist* hKFcxS = Hist::New("hKFcxS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhKFcx = Hist::ProjectAll(HistProj::kY, hKFcx);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhKFcx.at(it))()->GetBinCenter((*vhKFcx.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhKFcx.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhKFcx.at(it))()->Fit(gaus, "q0", "");
        (*vhKFcx.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhKFcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hKFcxM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hKFcxM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hKFcxS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hKFcxS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hKFcxM->write();
    hKFcxS->write();
    

    COUT("HCcx\n");
    Hist* hHCcxM = Hist::New("hHCcxM", HistAxis(AXrig, "Mean"));
    Hist* hHCcxS = Hist::New("hHCcxS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhHCcx = Hist::ProjectAll(HistProj::kY, hHCcx);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhHCcx.at(it))()->GetBinCenter((*vhHCcx.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhHCcx.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhHCcx.at(it))()->Fit(gaus, "q0", "");
        (*vhHCcx.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhHCcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCcx.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hHCcxM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hHCcxM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hHCcxS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hHCcxS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hHCcxM->write();
    hHCcxS->write();


    hCKcxM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFcxM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCcxM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chcxM = Hist::Collect("chcxM", HistList({ hCKcxM, hKFcxM, hHCcxM }));
    chcxM->Write();
    
    hCKcxS->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFcxS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCcxS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chcxS = Hist::Collect("chcxS", HistList({ hCKcxS, hKFcxS, hHCcxS }));
    chcxS->Write();


    COUT("CKcy\n");
    Hist* hCKcyM = Hist::New("hCKcyM", HistAxis(AXrig, "Mean"));
    Hist* hCKcyS = Hist::New("hCKcyS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhCKcy = Hist::ProjectAll(HistProj::kY, hCKcy);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhCKcy.at(it))()->GetBinCenter((*vhCKcy.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhCKcy.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhCKcy.at(it))()->Fit(gaus, "q0", "");
        (*vhCKcy.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhCKcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hCKcyM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hCKcyM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hCKcyS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hCKcyS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hCKcyM->write();
    hCKcyS->write();
    

    COUT("KFcy\n");
    Hist* hKFcyM = Hist::New("hKFcyM", HistAxis(AXrig, "Mean"));
    Hist* hKFcyS = Hist::New("hKFcyS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhKFcy = Hist::ProjectAll(HistProj::kY, hKFcy);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhKFcy.at(it))()->GetBinCenter((*vhKFcy.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhKFcy.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhKFcy.at(it))()->Fit(gaus, "q0", "");
        (*vhKFcy.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhKFcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hKFcyM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hKFcyM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hKFcyS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hKFcyS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hKFcyM->write();
    hKFcyS->write();
    

    COUT("HCcy\n");
    Hist* hHCcyM = Hist::New("hHCcyM", HistAxis(AXrig, "Mean"));
    Hist* hHCcyS = Hist::New("hHCcyS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhHCcy = Hist::ProjectAll(HistProj::kY, hHCcy);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhHCcy.at(it))()->GetBinCenter((*vhHCcy.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhHCcy.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhHCcy.at(it))()->Fit(gaus, "q0", "");
        (*vhHCcy.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhHCcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCcy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hHCcyM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hHCcyM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hHCcyS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hHCcyS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hHCcyM->write();
    hHCcyS->write();


    hCKcyM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFcyM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCcyM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chcyM = Hist::Collect("chcyM", HistList({ hCKcyM, hKFcyM, hHCcyM }));
    chcyM->Write();
    
    hCKcyS->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFcyS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCcyS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chcyS = Hist::Collect("chcyS", HistList({ hCKcyS, hKFcyS, hHCcyS }));
    chcyS->Write();
    
    
    COUT("CKux\n");
    Hist* hCKuxM = Hist::New("hCKuxM", HistAxis(AXrig, "Mean"));
    Hist* hCKuxS = Hist::New("hCKuxS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhCKux = Hist::ProjectAll(HistProj::kY, hCKux);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhCKux.at(it))()->GetBinCenter((*vhCKux.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhCKux.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhCKux.at(it))()->Fit(gaus, "q0", "");
        (*vhCKux.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhCKux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hCKuxM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hCKuxM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hCKuxS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hCKuxS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hCKuxM->write();
    hCKuxS->write();
    

    COUT("KFux\n");
    Hist* hKFuxM = Hist::New("hKFuxM", HistAxis(AXrig, "Mean"));
    Hist* hKFuxS = Hist::New("hKFuxS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhKFux = Hist::ProjectAll(HistProj::kY, hKFux);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhKFux.at(it))()->GetBinCenter((*vhKFux.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhKFux.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhKFux.at(it))()->Fit(gaus, "q0", "");
        (*vhKFux.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhKFux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hKFuxM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hKFuxM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hKFuxS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hKFuxS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hKFuxM->write();
    hKFuxS->write();
    

    COUT("HCux\n");
    Hist* hHCuxM = Hist::New("hHCuxM", HistAxis(AXrig, "Mean"));
    Hist* hHCuxS = Hist::New("hHCuxS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhHCux = Hist::ProjectAll(HistProj::kY, hHCux);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhHCux.at(it))()->GetBinCenter((*vhHCux.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhHCux.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhHCux.at(it))()->Fit(gaus, "q0", "");
        (*vhHCux.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhHCux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCux.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hHCuxM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hHCuxM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hHCuxS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hHCuxS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hHCuxM->write();
    hHCuxS->write();


    hCKuxM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFuxM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCuxM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chuxM = Hist::Collect("chuxM", HistList({ hCKuxM, hKFuxM, hHCuxM }));
    chuxM->Write();
    
    hCKuxS->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFuxS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCuxS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chuxS = Hist::Collect("chuxS", HistList({ hCKuxS, hKFuxS, hHCuxS }));
    chuxS->Write();


    COUT("CKuy\n");
    Hist* hCKuyM = Hist::New("hCKuyM", HistAxis(AXrig, "Mean"));
    Hist* hCKuyS = Hist::New("hCKuyS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhCKuy = Hist::ProjectAll(HistProj::kY, hCKuy);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhCKuy.at(it))()->GetBinCenter((*vhCKuy.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhCKuy.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhCKuy.at(it))()->Fit(gaus, "q0", "");
        (*vhCKuy.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhCKuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhCKuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hCKuyM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hCKuyM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hCKuyS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hCKuyS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hCKuyM->write();
    hCKuyS->write();
    

    COUT("KFuy\n");
    Hist* hKFuyM = Hist::New("hKFuyM", HistAxis(AXrig, "Mean"));
    Hist* hKFuyS = Hist::New("hKFuyS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhKFuy = Hist::ProjectAll(HistProj::kY, hKFuy);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhKFuy.at(it))()->GetBinCenter((*vhKFuy.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhKFuy.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhKFuy.at(it))()->Fit(gaus, "q0", "");
        (*vhKFuy.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhKFuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhKFuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hKFuyM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hKFuyM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hKFuyS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hKFuyS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hKFuyM->write();
    hKFuyS->write();
    

    COUT("HCuy\n");
    Hist* hHCuyM = Hist::New("hHCuyM", HistAxis(AXrig, "Mean"));
    Hist* hHCuyS = Hist::New("hHCuyS", HistAxis(AXrig, "Sigma"));
    std::vector<Hist*> vhHCuy = Hist::ProjectAll(HistProj::kY, hHCuy);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t scl = std::sqrt(1.0 / AXrig.center(it, AxisScale::kLog));
        Double_t max = (*vhHCuy.at(it))()->GetBinCenter((*vhHCuy.at(it))()->GetMaximumBin());
        Double_t rms = 0.5 * (*vhHCuy.at(it))()->GetRMS();
        gaus->SetParameters(1000, max, rms);
        (*vhHCuy.at(it))()->Fit(gaus, "q0", "");
        (*vhHCuy.at(it))()->Fit(gaus, "q0", "", max-stable*rms, max+stable*rms);
        (*vhHCuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
        (*vhHCuy.at(it))()->Fit(gaus, "q0", "", max-stable*gaus->GetParameter(2), max+stable*gaus->GetParameter(2));
    
        (*hHCuyM)()->SetBinContent(it, scl * gaus->GetParameter(1));
        (*hHCuyM)()->SetBinError  (it, scl * gaus->GetParError(1));
        (*hHCuyS)()->SetBinContent(it, scl * gaus->GetParameter(2));
        (*hHCuyS)()->SetBinError  (it, scl * gaus->GetParError(2));
    } 
    hHCuyM->write();
    hHCuyS->write();


    hCKuyM->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFuyM->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCuyM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chuyM = Hist::Collect("chuyM", HistList({ hCKuyM, hKFuyM, hHCuyM }));
    chuyM->Write();
    
    hCKuyS->style(Fill(), Line(kGreen), Marker(kGreen));
    hKFuyS->style(Fill(), Line(kBlue), Marker(kBlue));
    hHCuyS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chuyS = Hist::Collect("chuyS", HistList({ hCKuyS, hKFuyS, hHCuyS }));
    chuyS->Write();
   
   
    ofle->Write();
    ofle->Close();

    return 0;
}
