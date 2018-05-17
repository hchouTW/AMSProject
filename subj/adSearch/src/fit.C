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
    Hist* hMres = Hist::Head("hMres");
 
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

    COUT("Mres\n");
    Hist* hMresM = Hist::New("hMresM", HistAxis(AXrig, "Mean"));
    Hist* hMresS = Hist::New("hMresS", HistAxis(AXrig, "Sigma"));
    
    std::vector<Hist*> vhMres = Hist::ProjectAll(HistProj::kY, hMres);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        Double_t Mmax = (*vhMres.at(it))()->GetBinCenter((*vhMres.at(it))()->GetMaximumBin());
        Double_t Mrms = 0.5 * (*vhMres.at(it))()->GetRMS();
        gaus->SetParameters(1000, Mmax, Mrms);
        (*vhMres.at(it))()->Fit(gaus, "q0", "");
        (*vhMres.at(it))()->Fit(gaus, "q0", "", Mmax-stable*Mrms, Mmax+stable*Mrms);
        (*vhMres.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhMres.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhMres.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
        (*vhMres.at(it))()->Fit(gaus, "q0", "", Mmax-stable*gaus->GetParameter(2), Mmax+stable*gaus->GetParameter(2));
    
        (*hMresM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hMresM)()->SetBinError  (it, gaus->GetParError(1));
        (*hMresS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hMresS)()->SetBinError  (it, gaus->GetParError(2));
    } 
    
    hMresM->write();
    hMresS->write();

    ofle->Write();
    ofle->Close();

    return 0;
}
