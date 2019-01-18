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
    
    Hist* hMrso = Hist::Head("hmassC");
    Hist* hCKMrso = Hist::Head("hCKmassC");
    Hist* hHCMrso = Hist::Head("hHCmassC");
    
    Hist* hLMrso = Hist::Head("hLmassC");
    Hist* hCKLMrso = Hist::Head("hCKLmassC");
    Hist* hHCLMrso = Hist::Head("hHCLmassC");
    
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

    std::vector<Hist*> vhMrso   = Hist::ProjectAll(HistProj::kY, hMrso);
    std::vector<Hist*> vhCKMrso = Hist::ProjectAll(HistProj::kY, hCKMrso);
    std::vector<Hist*> vhHCMrso = Hist::ProjectAll(HistProj::kY, hHCMrso);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhMrso.at(it)->style(Line(kGreen+2, 0, 3), Marker(kGreen+2));
        vhCKMrso.at(it)->style(Line(kBlue, 0, 3), Marker(kBlue));
        vhHCMrso.at(it)->style(Line(kRed, 0, 3), Marker(kRed));
        THStack* cvhMrso = Hist::Collect(Form("cvhMrso%03d", it), HistList({ vhMrso.at(it), vhCKMrso.at(it), vhHCMrso.at(it) }));
        cvhMrso->SetMinimum(0.5);
        cvhMrso->Write();
    } 
    
    hLMrso->style(Line(kGreen+2, 0, 3), Marker(kGreen+2));
    hCKLMrso->style(Line(kBlue, 0, 3), Marker(kBlue));
    hHCLMrso->style(Line(kRed, 0, 3), Marker(kRed));
    THStack* chLMrso = Hist::Collect("chLMrso", HistList({ hLMrso, hCKLMrso, hHCLMrso }));
    chLMrso->SetMinimum(0.5);
    chLMrso->Write();
    
    hM1llx->style(Line(kBlue, 0, 3), Marker(kBlue));
    hM2llx->style(Line(kRed, 0, 3), Marker(kRed));
    THStack* chMllx = Hist::Collect("chMllx", HistList({ hM1llx, hM2llx }));
    chMllx->SetMinimum(0.5);
    chMllx->Write();
   
    hM1lly->style(Line(kBlue, 0, 3), Marker(kBlue));
    hM2lly->style(Line(kRed, 0, 3), Marker(kRed));
    THStack* chMlly = Hist::Collect("chMlly", HistList({ hM1lly, hM2lly }));
    chMlly->SetMinimum(0.5);
    chMlly->Write();
   
    ofle->Write();
    ofle->Close();

    return 0;
}
