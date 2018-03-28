#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>
    

static constexpr Double_t stable = 2.0;    
static TF1* fgaus = new TF1("fgaus", "gaus", -3.0, 3.0);

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
    
    //Hist::Load("time_fill.root", "dat");
    Hist::Load("time_fill.root", "/afs/cern.ch/work/h/hchou/AMSData/test22");

    // Prop
    Hist* hTme = Hist::Head("hTme");

    const Axis& AXeta = hTme->xaxis();
    const Axis& AXtme = hTme->yaxis();
    
    TFile * ofle = new TFile("time_fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();

    std::vector<Hist*> vhTme = Hist::ProjectAll(HistProj::kY, hTme);
    Hist* hTmeM = Hist::New("hTmeM", HistAxis(AXeta, "Mean"));
    Hist* hTmeS = Hist::New("hTmeS", HistAxis(AXeta, "Sigma"));
    for (int it = 35; it <= AXeta.nbin(); ++it) {
        COUT("Tme ITER %d\n", it);
        Double_t men = (*vhTme.at(it))()->GetBinCenter((*vhTme.at(it))()->GetMaximumBin());
        Double_t rms = (*vhTme.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhTme.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhTme.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTme.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTme.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hTmeM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hTmeM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hTmeS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hTmeS)()->SetBinError  (it, fgaus->GetParError(2));
    }

    hTmeM->Write();
    hTmeS->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
