#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

    constexpr Int_t niter = 5000;
    Double_t mpv = 0.0;
    Double_t sgm = 1.0;
    Double_t cut = 5000.0;

    Axis AXit("Iter", niter, 1, niter+1);
    Axis AXel("Eloss", 2000, -5, 200);
    Hist* hEloss = Hist::New("hEloss", HistAxis(AXit, AXel));
    
    for (Int_t ev = 0; ev < 100000; ++ev) {
        Double_t eloss = 0;
        Int_t iter = 1;
        while (iter <= niter) {
            Double_t el = mpv + sgm * MGROOT::Rndm::Landau();
            if (el > cut) continue;
            eloss += el;
            hEloss->fillH2D(iter, eloss/iter);
            iter++;
        }
    }

    TFile * ofle = new TFile(Form("/ams_home/hchou/AMSProject/subj/trsys/dat/landau%03ld.root", std::atol(argv[1])), "RECREATE");
    hEloss->Write();
    ofle->Close();

    return 1;
}
