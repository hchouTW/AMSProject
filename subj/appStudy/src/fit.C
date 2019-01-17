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
    
    Hist* hCKMrso    = Hist::Head("hCKmassC2");
    Hist* hHCMrso    = Hist::Head("hmassC2");
    Hist* hHCMrsoAll = Hist::Head("hmassC2All");
    
    Hist* hCKLMrso    = Hist::Head("hCKLmassC2");
    Hist* hHCLMrso    = Hist::Head("hLmassC2");
    Hist* hHCLMrsoAll = Hist::Head("hLmassC2All");
    
    Hist* hM1llx = Hist::Head("hPRM1llx");
    Hist* hM2llx = Hist::Head("hPRM2llx");
    Hist* hM1lly = Hist::Head("hPRM1lly");
    Hist* hM2lly = Hist::Head("hPRM2lly");
    
    const Axis& AXrig = hHCMrso->xaxis();
    
    TFile * ofle = new TFile("fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();
    const Double_t stable = 1.7;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);

    std::vector<Hist*> vhCKMrso       = Hist::ProjectAll(HistProj::kY, hCKMrso);
    std::vector<Hist*> vhHCMrso       = Hist::ProjectAll(HistProj::kY, hHCMrso);
    std::vector<Hist*> vhHCMrsoAll    = Hist::ProjectAll(HistProj::kY, hHCMrsoAll);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhCKMrso.at(it)->style(Line(kGreen+2, 0, 2), Marker(kGreen+2));
        vhHCMrso.at(it)->style(Line(kBlue, 0, 2), Marker(kBlue));
        vhHCMrsoAll.at(it)->style(Line(kRed, 0, 2), Marker(kRed));
        THStack* cvhMrso = Hist::Collect(Form("cvhMrso%03d", it), HistList({ vhCKMrso.at(it), vhHCMrso.at(it), vhHCMrsoAll.at(it) }));
        cvhMrso->SetMinimum(0.5);
        cvhMrso->Write();
    } 
    
    hCKLMrso->style(Line(kGreen+2, 0, 2), Marker(kGreen+2));
    hHCLMrso->style(Line(kBlue, 0, 2), Marker(kBlue));
    hHCLMrsoAll->style(Line(kRed, 0, 2), Marker(kRed));
    THStack* chLMrso = Hist::Collect("chLMrso", HistList({ hCKLMrso, hHCLMrso, hHCLMrsoAll }));
    chLMrso->SetMinimum(0.5);
    chLMrso->Write();
    
    hM1llx->style(Line(kBlue, 0, 2), Marker(kBlue));
    hM2llx->style(Line(kRed, 0, 2), Marker(kRed));
    THStack* chMllx = Hist::Collect("chMllx", HistList({ hM1llx, hM2llx }));
    chMllx->SetMinimum(0.5);
    chMllx->Write();
   
    hM1lly->style(Line(kBlue, 0, 2), Marker(kBlue));
    hM2lly->style(Line(kRed, 0, 2), Marker(kRed));
    THStack* chMlly = Hist::Collect("chMlly", HistList({ hM1lly, hM2lly }));
    chMlly->SetMinimum(0.5);
    chMlly->Write();
   
    ofle->Write();
    ofle->Close();

    return 0;
}
