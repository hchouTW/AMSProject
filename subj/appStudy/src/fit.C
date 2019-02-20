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
    Hist* hAD_TF_HCMrsoQ1 = Hist::Head("hAD_TF_massQ1");
    Hist* hAD_TF_HCMrsoQ2 = Hist::Head("hAD_TF_massQ2");
    Hist* hAD_TF_HCMrsoQ3 = Hist::Head("hAD_TF_massQ3");
    Hist* hAD_TF_HCMrsoQ4 = Hist::Head("hAD_TF_massQ4");
    Hist* hAD_TF_HCMrsoQ5 = Hist::Head("hAD_TF_massQ5");
    Hist* hAD_TF_HCMrsoQ6 = Hist::Head("hAD_TF_massQ6");
    
    Hist* hAD_RH_CKMrso  = Hist::Head("hAD_RH_CKmass");
    Hist* hAD_RH_CKMrsoQ = Hist::Head("hAD_RH_CKmassQ");
    Hist* hAD_RH_CKMrsoQ1 = Hist::Head("hAD_RH_CKmassQ1");
    Hist* hAD_RH_CKMrsoQ2 = Hist::Head("hAD_RH_CKmassQ2");
    Hist* hAD_RH_CKMrsoQ3 = Hist::Head("hAD_RH_CKmassQ3");
    Hist* hAD_RH_HCMrso  = Hist::Head("hAD_RH_mass");
    Hist* hAD_RH_HCMrsoQ = Hist::Head("hAD_RH_massQ");
    Hist* hAD_RH_HCMrsoQ1 = Hist::Head("hAD_RH_massQ1");
    Hist* hAD_RH_HCMrsoQ2 = Hist::Head("hAD_RH_massQ2");
    Hist* hAD_RH_HCMrsoQ3 = Hist::Head("hAD_RH_massQ3");
    
    const Axis& AXrig = hAD_TF_HCMrso->xaxis();
    
    TFile * ofle = new TFile("fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();
    //const Double_t stable = 1.7;
    //TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);

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
    
    //hAD_TF_CKMrso ->style(Line(kBlue+1, 0, 1), Marker(kBlue-1));
    hAD_TF_CKMrsoQ->style(Line(kBlue, 0, 3), Marker(kBlue));
    //hAD_TF_HCMrso ->style(Line(kRed+1, 0, 1), Marker(kRed-1));
    hAD_TF_HCMrsoQ->style(Line(kRed, 0, 3), Marker(kRed));
    hAD_TF_HCMrsoQ1->style(Line(kGreen+0, 0, 3), Marker(kGreen+0));
    hAD_TF_HCMrsoQ2->style(Line(kGreen+1, 0, 3), Marker(kGreen+1));
    hAD_TF_HCMrsoQ3->style(Line(kGreen+2, 0, 3), Marker(kGreen+2));
    hAD_TF_HCMrsoQ4->style(Line(kGreen+3, 0, 3), Marker(kGreen+3));
    hAD_TF_HCMrsoQ5->style(Line(kGreen+4, 0, 3), Marker(kGreen+4));
    hAD_TF_HCMrsoQ6->style(Line(kGreen+5, 0, 3), Marker(kGreen+5));
    //THStack* chAD_TF_Mrso = Hist::Collect("chAD_TF_Mrso", HistList({ hAD_TF_CKMrso, hAD_TF_HCMrso, hAD_TF_CKMrsoQ, hAD_TF_HCMrsoQ }));
    THStack* chAD_TF_Mrso = Hist::Collect("chAD_TF_Mrso", HistList({ hAD_TF_CKMrsoQ, hAD_TF_HCMrsoQ, hAD_TF_HCMrsoQ1, hAD_TF_HCMrsoQ2, hAD_TF_HCMrsoQ3, hAD_TF_HCMrsoQ4, hAD_TF_HCMrsoQ5, hAD_TF_HCMrsoQ6 }));
    chAD_TF_Mrso->SetMinimum(0.5);
    chAD_TF_Mrso->Write();
   
    hAD_RH_CKMrso ->style(Line(kBlue+1, 0, 1), Marker(kBlue+1));
    hAD_RH_CKMrsoQ->style(Line(kBlue, 0, 3), Marker(kBlue));
    hAD_RH_HCMrso ->style(Line(kRed+1, 0, 1), Marker(kRed+1));
    hAD_RH_HCMrsoQ->style(Line(kRed, 0, 3), Marker(kRed));
    THStack* chAD_RH_Mrso = Hist::Collect("chAD_RH_Mrso", HistList({ hAD_RH_CKMrso, hAD_RH_HCMrso, hAD_RH_CKMrsoQ, hAD_RH_HCMrsoQ }));
    chAD_RH_Mrso->SetMinimum(0.5);
    chAD_RH_Mrso->Write();
    
    hAD_RH_CKMrsoQ1->style(Line(kBlue+1, 0, 2), Marker(kBlue));
    hAD_RH_CKMrsoQ2->style(Line(kBlue+2, 0, 2), Marker(kBlue));
    hAD_RH_CKMrsoQ3->style(Line(kBlue+3, 0, 2), Marker(kBlue));
    THStack* chAD_RH_MCK = Hist::Collect("chAD_RH_MCK", HistList({ hAD_RH_CKMrsoQ, hAD_RH_CKMrsoQ1, hAD_RH_CKMrsoQ2, hAD_RH_CKMrsoQ3 }));
    chAD_RH_MCK->SetMinimum(0.5);
    chAD_RH_MCK->Write();
    
    hAD_RH_HCMrsoQ1->style(Line(kRed+1, 0, 2), Marker(kRed));
    hAD_RH_HCMrsoQ2->style(Line(kRed+2, 0, 2), Marker(kRed));
    hAD_RH_HCMrsoQ3->style(Line(kRed+3, 0, 2), Marker(kRed));
    THStack* chAD_RH_MHC = Hist::Collect("chAD_RH_MHC", HistList({ hAD_RH_HCMrsoQ, hAD_RH_HCMrsoQ1, hAD_RH_HCMrsoQ2, hAD_RH_HCMrsoQ3 }));
    chAD_RH_MHC->SetMinimum(0.5);
    chAD_RH_MHC->Write();
   
    ofle->Write();
    ofle->Close();

    return 0;
}
