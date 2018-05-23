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
    //Hist::Load("fill.root", "/afs/cern.ch/work/h/hchou/AMSData/test58");

    // Fit
    Hist* hMrso = Hist::Head("hMpos");
    Hist* hM2rso = Hist::Head("hM2pos");
 
    const Axis& AXrig = hMrso->xaxis();
    
    TFile * ofle = new TFile("fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();
    
    const Double_t stable = 1.7;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);

    COUT("Mrso\n");
    Hist* hMrsoM = Hist::New("hMrsoM", HistAxis(AXrig, "Mass Mean [GeV]"));
    Hist* hMrsoS = Hist::New("hMrsoS", HistAxis(AXrig, "Mass Sigma [GeV]"));
    
    std::vector<Hist*> vhMrso = Hist::ProjectAll(HistProj::kY, hMrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t Mmax = (*vhMrso.at(it))()->GetBinCenter((*vhMrso.at(it))()->GetMaximumBin());
        Double_t Mrms = 0.5 * (*vhMrso.at(it))()->GetRMS();
        gaus->SetParameters(1000, Mmax, Mrms);
        (*vhMrso.at(it))()->Fit(gaus, "q0", "");
        (*vhMrso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*Mrms, Mmax+stable*Mrms);
        (*vhMrso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhMrso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhMrso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhMrso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
    
        (*hMrsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hMrsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hMrsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hMrsoS)()->SetBinError  (it, gaus->GetParError(2));
    } 
    
    hMrsoM->write();
    hMrsoS->write();
    
    COUT("M2rso\n");
    Hist* hM2rsoM = Hist::New("hM2rsoM", HistAxis(AXrig, "Mass Mean [GeV]"));
    Hist* hM2rsoS = Hist::New("hM2rsoS", HistAxis(AXrig, "Mass Sigma [GeV]"));
    
    std::vector<Hist*> vhM2rso = Hist::ProjectAll(HistProj::kY, hM2rso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t Mmax = (*vhM2rso.at(it))()->GetBinCenter((*vhM2rso.at(it))()->GetMaximumBin());
        Double_t Mrms = 0.5 * (*vhM2rso.at(it))()->GetRMS();
        gaus->SetParameters(1000, Mmax, Mrms);
        (*vhM2rso.at(it))()->Fit(gaus, "q0", "");
        (*vhM2rso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*Mrms, Mmax+stable*Mrms);
        (*vhM2rso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhM2rso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhM2rso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhM2rso.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
    
        (*hM2rsoM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hM2rsoM)()->SetBinError  (it, gaus->GetParError(1));
        (*hM2rsoS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hM2rsoS)()->SetBinError  (it, gaus->GetParError(2));
    } 
    
    hM2rsoM->write();
    hM2rsoS->write();

    ofle->Write();
    ofle->Close();

    return 0;
}
