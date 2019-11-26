#include <CPPLibs.h>
#include <ROOTLibs.h>

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    if (argc != 2) return 0;
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/trsys/%s", argv[1]));

    PdfEditor editor(Window(), "fit_mscat", "out/doc");

    TF1* func = new TF1("func", "gaus");

    const Axis& AXigb = Hist::Head("hUxMC")->xaxis();

    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<Hist*> vhMC_igb = Hist::ProjectAll(HistProj::kY, Hist::Head("hUxMC"));
    std::vector<Hist*> vhSM_igb = Hist::ProjectAll(HistProj::kY, Hist::Head("hUxSM"));
    
    Hist* hMCmen = Hist::New("hMCmen", HistAxis(AXigb));
    Hist* hSMmen = Hist::New("hSMmen", HistAxis(AXigb));
    
    Hist* hMCrso = Hist::New("hMCrso", HistAxis(AXigb));
    Hist* hSMrso = Hist::New("hSMrso", HistAxis(AXigb));
    
    Hist* hMCrelrso = Hist::New("hMCrelrso", HistAxis(AXigb));
    Hist* hSMrelrso = Hist::New("hSMrelrso", HistAxis(AXigb));

    for (int it = 1; it <= AXigb.nbin(); ++it) {
        double igb = AXigb.center(it, AxisScale::kLog);
        double scl = igb;

        vhMC_igb.at(it)->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle )));
        vhSM_igb.at(it)->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

        THStack* stack = Hist::Collect(Form("h_igb_%03d", it), HistList({ vhMC_igb.at(it), vhSM_igb.at(it) }));
        stack->SetTitle(Form("1/(#gamma#beta) (%6.4f ~ %6.4f)", AXigb()(it-1), AXigb()(it)));
        //vhstack.push_back(stack);

        editor.create();
        editor()().SetLogx(0);
        stack->Draw("nostack hist");
        stack->GetHistogram()->SetLineColor(0);
        stack->GetHistogram()->SetMarkerColor(0);
        stack->GetHistogram()->GetXaxis()->SetTitle("Cos-Angle Difference");
        stack->GetHistogram()->GetYaxis()->SetTitle("Events / Bin");
        stack->Draw("nostack hist");
        Legend leg("", PadWindow(0.60, 0.80, 0.65, 0.80));
        leg()->AddEntry((TObject*)nullptr, stack->GetTitle(), "");
        leg()->AddEntry((*vhMC_igb.at(it))(), "Monte Carlo", "lp");
        leg()->AddEntry((*vhSM_igb.at(it))(), "Mixture", "lp");
        leg()->SetFillColor(0);
        leg.draw();
        editor.save();

        const double width = 1.0;

        (*vhMC_igb.at(it))()->Fit(func, "", "q0", -1, 1);
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
        
        (*vhSM_igb.at(it))()->Fit(func, "", "q0", -1, 1);
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
        
        double relrso_mc = (*hMCrso)()->GetBinContent(it) / (*hMCrso)()->GetBinContent(it);
        double relrso_sm = (*hSMrso)()->GetBinContent(it) / (*hMCrso)()->GetBinContent(it);

        (*hMCrelrso)()->SetBinContent(it, relrso_mc);
        (*hMCrelrso)()->SetBinError  (it, 0.0);

        (*hSMrelrso)()->SetBinContent(it, relrso_sm);
        (*hSMrelrso)()->SetBinError  (it, 0.0);
    }
    
    hMCmen->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle )));
    hSMmen->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hMCrso->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle )));
    hSMrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hMCrelrso->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle )));
    hSMrelrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

    THStack* combined_men = Hist::Collect("h_combind_men", HistList({ hMCmen, hSMmen }));
    vhstack.push_back(combined_men);

    vhist.push_back(hMCmen);
    vhist.push_back(hSMmen);

    THStack* combined_rso = Hist::Collect("h_combind_rso", HistList({ hMCrso, hSMrso }));
    vhstack.push_back(combined_rso);
    
    vhist.push_back(hMCrso);
    vhist.push_back(hSMrso);
    
    THStack* combined_relrso = Hist::Collect("h_combind_relrso", HistList({ hMCrelrso, hSMrelrso }));
    vhstack.push_back(combined_relrso);
    
    vhist.push_back(hMCrelrso);
    vhist.push_back(hSMrelrso);
    
    editor.create();
    editor()().SetLogx();
    combined_men->Draw("nostack hist");
    combined_men->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_men->GetHistogram()->SetLineColor(0);
    combined_men->GetHistogram()->SetMarkerColor(0);
    combined_men->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
    combined_men->GetHistogram()->GetYaxis()->SetTitle("Cos-Angle Difference");
    combined_men->Draw("nostack p");
    Legend leg_men("", PadWindow(0.60, 0.80, 0.55, 0.70));
    leg_men()->AddEntry((*hMCmen)(), "Monte Carlo", "lp");
    leg_men()->AddEntry((*hSMmen)(), "ToyMC", "lp");
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
    combined_rso->GetHistogram()->GetYaxis()->SetTitle("Cos-Angle Difference");
    combined_rso->Draw("nostack p");
    Legend leg_rso("", PadWindow(0.60, 0.80, 0.55, 0.70));
    leg_rso()->AddEntry((*hMCrso)(), "Monte Carlo", "lp");
    leg_rso()->AddEntry((*hSMrso)(), "ToyMC", "lp");
    leg_rso()->SetFillColor(0);
    leg_rso.draw();
    editor.save();
    
    editor.create();
    editor()().SetLogx();
    combined_relrso->Draw("nostack hist");
    combined_relrso->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_relrso->GetHistogram()->SetLineColor(0);
    combined_relrso->GetHistogram()->SetMarkerColor(0);
    combined_relrso->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
    combined_relrso->GetHistogram()->GetYaxis()->SetTitle("Relative Resolution");
    combined_relrso->Draw("nostack p");
    Legend leg_relrso("", PadWindow(0.60, 0.80, 0.25, 0.40));
    leg_relrso()->AddEntry((*hMCrelrso)(), "Monte Carlo", "lp");
    leg_relrso()->AddEntry((*hSMrelrso)(), "ToyMC", "lp");
    leg_relrso()->SetFillColor(0);
    leg_relrso.draw();
    editor.save();
    
    editor.close();

    TFile * ofle = new TFile("out/doc/fit_mscat.root", "RECREATE");
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
