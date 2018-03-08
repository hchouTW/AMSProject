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
   
    Hist::Load("hit_fill.root", "dat");

    // Fit
    Hist* hMedep = Hist::Head("hMay");

    const Axis& AXeta = hMedep->xaxis();
    
    TFile * ofle = new TFile("hit_fit.root", "RECREATE");
    ofle->cd();
    
    Hist::AddDirectory();

    Hist* hMedepP = Hist::New("hMedepP", HistAxis(AXeta, "Peak"));
    Hist* hMedepM = Hist::New("hMedepM", HistAxis(AXeta, "Mean"));
    Hist* hMedepS = Hist::New("hMedepS", HistAxis(AXeta, "Sigma"));
    
    const Double_t stable = 1.5;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    std::vector<Hist*> vhMedep = Hist::ProjectAll(HistProj::kY, hMedep);

    for (int it = 1; it <= AXeta.nbin(); ++it) {
        Double_t Mmax = (*vhMedep.at(it))()->GetBinCenter((*vhMedep.at(it))()->GetMaximumBin());
        Double_t Mrms = (*vhMedep.at(it))()->GetRMS();
        gaus->SetParameters(1000, Mmax, Mrms);
        (*vhMedep.at(it))()->Fit(gaus, "q0", "");
        (*vhMedep.at(it))()->Fit(gaus, "q0", "", Mmax-stable*Mrms, Mmax+stable*Mrms);
        (*vhMedep.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhMedep.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
        (*vhMedep.at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
        (*hMedepP)()->SetBinContent(it, Mmax);
        (*hMedepP)()->SetBinError  (it, 1.0e-6);
        (*hMedepM)()->SetBinContent(it, gaus->GetParameter(1));
        (*hMedepM)()->SetBinError  (it, gaus->GetParError(1));
        (*hMedepS)()->SetBinContent(it, gaus->GetParameter(2));
        (*hMedepS)()->SetBinError  (it, gaus->GetParError(2));
    }
   
    ofle->Write();
    ofle->Close();

    return 0;
}
