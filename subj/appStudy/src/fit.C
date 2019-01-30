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
    
    //Hist* hMrso = Hist::Head("hmassC");
    //Hist* hCKMrso = Hist::Head("hCKmassC");
    //Hist* hHCMrso = Hist::Head("hHCmassC");
    
    Hist* hAD_TF_CKMrso  = Hist::Head("hAD_TF_CKmass");
    Hist* hAD_TF_CKMrsoQ = Hist::Head("hAD_TF_CKmassQ");
    Hist* hAD_TF_HCMrso  = Hist::Head("hAD_TF_mass");
    Hist* hAD_TF_HCMrsoQ = Hist::Head("hAD_TF_massQ");
    
    Hist* hAD_RH_CKMrso  = Hist::Head("hAD_RH_CKmass");
    Hist* hAD_RH_CKMrsoQ = Hist::Head("hAD_RH_CKmassQ");
    Hist* hAD_RH_HCMrso  = Hist::Head("hAD_RH_mass");
    Hist* hAD_RH_HCMrsoQ = Hist::Head("hAD_RH_massQ");
    
    const Axis& AXrig = hAD_TF_HCMrso->xaxis();
    
    TFile * ofle = new TFile("fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();
    const Double_t stable = 1.7;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);

    //std::vector<Hist*> vhMrso   = Hist::ProjectAll(HistProj::kY, hMrso);
    //std::vector<Hist*> vhCKMrso = Hist::ProjectAll(HistProj::kY, hCKMrso);
    //std::vector<Hist*> vhHCMrso = Hist::ProjectAll(HistProj::kY, hHCMrso);
    //for (int it = 1; it <= AXrig.nbin(); ++it) {
    //    vhMrso.at(it)->style(Line(kGreen+2, 0, 3), Marker(kGreen+2));
    //    vhCKMrso.at(it)->style(Line(kBlue, 0, 3), Marker(kBlue));
    //    vhHCMrso.at(it)->style(Line(kRed, 0, 3), Marker(kRed));
    //    THStack* cvhMrso = Hist::Collect(Form("cvhMrso%03d", it), HistList({ vhMrso.at(it), vhCKMrso.at(it), vhHCMrso.at(it) }));
    //    cvhMrso->SetMinimum(0.5);
    //    cvhMrso->Write();
    //} 
    
    hAD_TF_CKMrso ->style(Line(kBlue+1, 0, 1), Marker(kBlue-1));
    hAD_TF_CKMrsoQ->style(Line(kBlue, 0, 3), Marker(kBlue));
    hAD_TF_HCMrso ->style(Line(kRed+1, 0, 1), Marker(kRed-1));
    hAD_TF_HCMrsoQ->style(Line(kRed, 0, 3), Marker(kRed));
    THStack* chAD_TF_Mrso = Hist::Collect("chAD_TF_Mrso", HistList({ hAD_TF_CKMrso, hAD_TF_HCMrso, hAD_TF_CKMrsoQ, hAD_TF_HCMrsoQ }));
    chAD_TF_Mrso->SetMinimum(0.5);
    chAD_TF_Mrso->Write();
   
    hAD_RH_CKMrso ->style(Line(kBlue+1, 0, 1), Marker(kBlue+1));
    hAD_RH_CKMrsoQ->style(Line(kBlue, 0, 3), Marker(kBlue));
    hAD_RH_HCMrso ->style(Line(kRed+1, 0, 1), Marker(kRed+1));
    hAD_RH_HCMrsoQ->style(Line(kRed, 0, 3), Marker(kRed));
    THStack* chAD_RH_Mrso = Hist::Collect("chAD_RH_Mrso", HistList({ hAD_RH_CKMrso, hAD_RH_HCMrso, hAD_RH_CKMrsoQ, hAD_RH_HCMrsoQ }));
    chAD_RH_Mrso->SetMinimum(0.5);
    chAD_RH_Mrso->Write();
   
    ofle->Write();
    ofle->Close();

    return 0;
}
