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
    
    Hist* hOFFbta = Hist::Head("hNEWbetaCut");
    Hist* hNEWbta = Hist::Head("hNEWbetaCut3");
    
    Hist* hOFFmass = Hist::Head("hNEWmassCut");
    Hist* hNEWmass = Hist::Head("hNEWmassCut3");
    
    const Axis& AXrig = hOFFbta->xaxis();
    const Axis& AXbta = hOFFmass->xaxis();
    
    TFile * ofle = new TFile("fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();
    
    const Double_t stable = 1.7;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    
    std::vector<Hist*> vhOFFbta = Hist::ProjectAll(HistProj::kY, hOFFbta);
    std::vector<Hist*> vhNEWbta = Hist::ProjectAll(HistProj::kY, hNEWbta);
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        vhOFFbta.at(it)->style(Line(kBlue), Marker(kBlue));
        vhNEWbta.at(it)->style(Line(kRed), Marker(kRed));
        THStack* chbta = Hist::Collect(Form("chbta%03d", it), HistList({ vhOFFbta.at(it), vhNEWbta.at(it) }));
        chbta->SetTitle(Form("Rigidity (%6.2f ~ %6.2f)", AXrig()(it-1), AXrig()(it)));
        chbta->Write();
    }
    
    std::vector<Hist*> vhOFFmass = Hist::ProjectAll(HistProj::kY, hOFFmass);
    std::vector<Hist*> vhNEWmass = Hist::ProjectAll(HistProj::kY, hNEWmass);
    for (int it = 1; it <= AXbta.nbin(); ++it) {
        vhOFFmass.at(it)->style(Line(kBlue), Marker(kBlue));
        vhNEWmass.at(it)->style(Line(kRed), Marker(kRed));
        THStack* chmass = Hist::Collect(Form("chmass%03d", it), HistList({ vhOFFmass.at(it), vhNEWmass.at(it) }));
        chmass->SetTitle(Form("Beta (%6.3f ~ %6.3f)", AXbta()(it-1), AXbta()(it)));
        chmass->Write();
    }

    ofle->Write();
    ofle->Close();

    return 0;
}
