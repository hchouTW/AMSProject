#include <CPPLibs.h>
#include <ROOTLibs.h>

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    if (argc != 2) return 0;
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/trsys/%s", argv[1]));

    PdfEditor editor(Window(), "fit_eloss", "out/doc");

    TF1* func = new TF1("func", "gaus");

    const Axis& AXigb = Hist::Head("hEaMC")->xaxis();

    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<Hist*> vhMC_igb = Hist::ProjectAll(HistProj::kY, Hist::Head("hEaMC"));
    std::vector<Hist*> vhLL_igb = Hist::ProjectAll(HistProj::kY, Hist::Head("hEaLD"));
    std::vector<Hist*> vhBB_igb = Hist::ProjectAll(HistProj::kY, Hist::Head("hEaBB"));
    std::vector<Hist*> vhSM_igb = Hist::ProjectAll(HistProj::kY, Hist::Head("hEaSM"));
    
    Hist* hMCmen = Hist::New("hMCmen", HistAxis(AXigb));
    Hist* hLLmen = Hist::New("hLLmen", HistAxis(AXigb));
    Hist* hBBmen = Hist::New("hBBmen", HistAxis(AXigb));
    Hist* hSMmen = Hist::New("hSMmen", HistAxis(AXigb));
    
    Hist* hMCrso = Hist::New("hMCrso", HistAxis(AXigb));
    Hist* hLLrso = Hist::New("hLLrso", HistAxis(AXigb));
    Hist* hBBrso = Hist::New("hBBrso", HistAxis(AXigb));
    Hist* hSMrso = Hist::New("hSMrso", HistAxis(AXigb));
    
    Hist* hMCrelmen = Hist::New("hMCrelmen", HistAxis(AXigb));
    Hist* hLLrelmen = Hist::New("hLLrelmen", HistAxis(AXigb));
    Hist* hBBrelmen = Hist::New("hBBrelmen", HistAxis(AXigb));
    Hist* hSMrelmen = Hist::New("hSMrelmen", HistAxis(AXigb));

    for (int it = 1; it <= AXigb.nbin(); ++it) {
        double igb = AXigb.center(it, AxisScale::kLog);
        double scl = (1.0 + igb * igb);

        vhMC_igb.at(it)->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle )));
        vhLL_igb.at(it)->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kDiamond)));
        vhBB_igb.at(it)->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kSquare )));
        vhSM_igb.at(it)->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

        THStack* stack = Hist::Collect(Form("h_igb_%03d", it), HistList({ vhMC_igb.at(it), vhLL_igb.at(it), vhBB_igb.at(it), vhSM_igb.at(it) }));
        stack->SetTitle(Form("1/(#gamma#beta) (%6.4f ~ %6.4f)", AXigb()(it-1), AXigb()(it)));
        //vhstack.push_back(stack);
        //stack->Write();

        editor.create();
        editor()().SetLogx(0);
        stack->Draw("nostack hist");
        stack->GetHistogram()->SetLineColor(0);
        stack->GetHistogram()->SetMarkerColor(0);
        stack->GetHistogram()->GetXaxis()->SetTitle("Energy Loss [MeV]");
        stack->GetHistogram()->GetYaxis()->SetTitle("Events / Bin");
        stack->Draw("nostack hist");
        Legend leg("", PadWindow(0.60, 0.80, 0.65, 0.80));
        leg()->AddEntry((TObject*)nullptr, stack->GetTitle(), "");
        leg()->AddEntry((*vhMC_igb.at(it))(), "Monte Carlo", "lp");
        leg()->AddEntry((*vhLL_igb.at(it))(), "Landau", "lp");
        leg()->AddEntry((*vhBB_igb.at(it))(), "Bathe", "lp");
        leg()->AddEntry((*vhSM_igb.at(it))(), "Mixture", "lp");
        leg()->SetFillColor(0);
        leg.draw();
        editor.save();

        const double width = 1.0;

        (*vhMC_igb.at(it))()->Fit(func, "", "q0", 0, 400);
        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -8.0 * width * func->GetParameter(2) + func->GetParameter(1), 8.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -5.0 * width * func->GetParameter(2) + func->GetParameter(1), 5.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hMCrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hMCrso)()->SetBinError  (it, scl * func->GetParError(2));

        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hMCmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hMCmen)()->SetBinError  (it, scl * func->GetParError(1));
        
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", 0, 400);
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -8.0 * width * func->GetParameter(2) + func->GetParameter(1), 8.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -5.0 * width * func->GetParameter(2) + func->GetParameter(1), 5.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hLLrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hLLrso)()->SetBinError  (it, scl * func->GetParError(2));

        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhLL_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hLLmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hLLmen)()->SetBinError  (it, scl * func->GetParError(1));
        
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", 0, 400);
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -8.0 * width * func->GetParameter(2) + func->GetParameter(1), 8.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -5.0 * width * func->GetParameter(2) + func->GetParameter(1), 5.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hBBrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hBBrso)()->SetBinError  (it, scl * func->GetParError(2));

        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhBB_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hBBmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hBBmen)()->SetBinError  (it, scl * func->GetParError(1));
        
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", 0, 400);
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -8.0 * width * func->GetParameter(2) + func->GetParameter(1), 8.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -5.0 * width * func->GetParameter(2) + func->GetParameter(1), 5.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hSMrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hSMrso)()->SetBinError  (it, scl * func->GetParError(2));

        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hSMmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hSMmen)()->SetBinError  (it, scl * func->GetParError(1));
        
        double relmen_mc = (*hMCmen)()->GetBinContent(it) / (*hMCmen)()->GetBinContent(it);
        double relmen_ll = (*hLLmen)()->GetBinContent(it) / (*hMCmen)()->GetBinContent(it);
        double relmen_bb = (*hBBmen)()->GetBinContent(it) / (*hMCmen)()->GetBinContent(it);
        double relmen_sm = (*hSMmen)()->GetBinContent(it) / (*hMCmen)()->GetBinContent(it);

        (*hMCrelmen)()->SetBinContent(it, relmen_mc);
        (*hMCrelmen)()->SetBinError  (it, 0.0);
        
        (*hLLrelmen)()->SetBinContent(it, relmen_ll);
        (*hLLrelmen)()->SetBinError  (it, 0.0);
        
        (*hBBrelmen)()->SetBinContent(it, relmen_bb);
        (*hBBrelmen)()->SetBinError  (it, 0.0);

        (*hSMrelmen)()->SetBinContent(it, relmen_sm);
        (*hSMrelmen)()->SetBinError  (it, 0.0);
    }
    
    hMCmen->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle )));
    hLLmen->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kDiamond)));
    hBBmen->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kSquare )));
    hSMmen->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hMCrso->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle )));
    hLLrso->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kDiamond)));
    hBBrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kSquare )));
    hSMrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hMCrelmen->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle )));
    hLLrelmen->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kDiamond)));
    hBBrelmen->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kSquare )));
    hSMrelmen->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

    THStack* combined_men = Hist::Collect("h_combind_men", HistList({ hMCmen, hLLmen, hBBmen, hSMmen }));
    vhstack.push_back(combined_men);

    vhist.push_back(hMCmen);
    vhist.push_back(hLLmen);
    vhist.push_back(hBBmen);
    vhist.push_back(hSMmen);

    THStack* combined_rso = Hist::Collect("h_combind_rso", HistList({ hMCrso, hLLrso, hBBrso, hSMrso }));
    vhstack.push_back(combined_rso);
    
    vhist.push_back(hMCrso);
    vhist.push_back(hLLrso);
    vhist.push_back(hBBrso);
    vhist.push_back(hSMrso);
    
    THStack* combined_relmen = Hist::Collect("h_combind_relmen", HistList({ hMCrelmen, hLLrelmen, hBBrelmen, hSMrelmen }));
    vhstack.push_back(combined_relmen);
    
    vhist.push_back(hMCrelmen);
    vhist.push_back(hLLrelmen);
    vhist.push_back(hBBrelmen);
    vhist.push_back(hSMrelmen);
    
    editor.create();
    editor()().SetLogx();
    combined_men->Draw("nostack hist");
    combined_men->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_men->GetHistogram()->SetLineColor(0);
    combined_men->GetHistogram()->SetMarkerColor(0);
    combined_men->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
    combined_men->GetHistogram()->GetYaxis()->SetTitle("Energy Loss [MeV]");
    combined_men->Draw("nostack p");
    Legend leg_men("", PadWindow(0.60, 0.80, 0.55, 0.70));
    leg_men()->AddEntry((*hMCmen)(), "Monte Carlo", "lp");
    leg_men()->AddEntry((*hLLmen)(), "Landau", "lp");
    leg_men()->AddEntry((*hBBmen)(), "Bathe", "lp");
    leg_men()->AddEntry((*hSMmen)(), "Mixture", "lp");
    leg_men()->SetFillColor(0);
    leg_men.draw();
    editor.save();
    
    editor.create();
    editor()().SetLogx();
    combined_rso->Draw("nostack hist");
    combined_rso->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_rso->GetHistogram()->SetLineColor(0);
    combined_rso->GetHistogram()->SetMarkerColor(0);
    combined_rso->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
    combined_rso->GetHistogram()->GetYaxis()->SetTitle("Energy Loss [MeV]");
    combined_rso->Draw("nostack p");
    Legend leg_rso("", PadWindow(0.60, 0.80, 0.55, 0.70));
    leg_rso()->AddEntry((*hMCrso)(), "Monte Carlo", "lp");
    leg_rso()->AddEntry((*hLLrso)(), "Landau", "lp");
    leg_rso()->AddEntry((*hBBrso)(), "Bathe", "lp");
    leg_rso()->AddEntry((*hSMrso)(), "Mixture", "lp");
    leg_rso()->SetFillColor(0);
    leg_rso.draw();
    editor.save();
    
    editor.create();
    editor()().SetLogx();
    combined_relmen->Draw("nostack hist");
    combined_relmen->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_relmen->GetHistogram()->SetLineColor(0);
    combined_relmen->GetHistogram()->SetMarkerColor(0);
    combined_relmen->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
    combined_relmen->GetHistogram()->GetYaxis()->SetTitle("Relative Resolution");
    combined_relmen->Draw("nostack p");
    Legend leg_relmen("", PadWindow(0.60, 0.80, 0.55, 0.70));
    leg_relmen()->AddEntry((*hMCrelmen)(), "Monte Carlo", "lp");
    leg_relmen()->AddEntry((*hLLrelmen)(), "Landau", "lp");
    leg_relmen()->AddEntry((*hBBrelmen)(), "Bathe", "lp");
    leg_relmen()->AddEntry((*hSMrelmen)(), "Mixture", "lp");
    leg_relmen()->SetFillColor(0);
    leg_relmen.draw();
    editor.save();

    editor.close();

    TFile * ofle = new TFile("out/doc/fit_eloss.root", "RECREATE");
    ofle->cd();
   
    for (auto&& stack : vhstack) {
        stack->Write();
    }
  
    for (auto&& hist : vhist) {
        (*hist)()->Write();
    }

    for (auto&& hist : vhist) {
        TGraphErrors* gr = new TGraphErrors((*hist)());
        gr->SetName(Form("gr_%s", (*hist)()->GetName()));
        gr->Write();
    }

    ofle->Write();
    ofle->Close();

    return 1;
}
