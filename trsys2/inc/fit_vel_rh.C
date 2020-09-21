#include <CPPLibs.h>
#include <ROOTLibs.h>

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    if (argc != 2) return 0;
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/trsys/%s", argv[1]));

    PdfEditor editor(Window(), "fit_vel_rh", "out/doc");

    TF1* func = new TF1("func", "gaus");

    const Axis& AXbta = Hist::Head("hRH_B_bta")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<double> rh_rng({ 0.96, 0.999 });

    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXbta));
    
    for (int it = 1; it <= AXbta.nbin(); ++it) {
        double tme_hc = (*Hist::Head("hHC_B_numt"))()->GetBinContent(it) / (*Hist::Head("hHC_B_dent"))()->GetBinContent(it);
        
        (*hHCtme)()->SetBinContent(it, tme_hc);
        (*hHCtme)()->SetBinError  (it, 0.0);
    }
        
    std::vector<Hist*> vhOF_bta = Hist::ProjectAll(HistProj::kY, Hist::Head("hRH_B_bta"));
    std::vector<Hist*> vhHC_bta = Hist::ProjectAll(HistProj::kY, Hist::Head("hHC_B_bta"));
    
    Hist* hOFeft = Hist::New("hOFeft", HistAxis(AXbta));
    Hist* hHCeft = Hist::New("hHCeft", HistAxis(AXbta));

    Hist* hOFmen = Hist::New("hOFmen", HistAxis(AXbta));
    Hist* hHCmen = Hist::New("hHCmen", HistAxis(AXbta));
    
    Hist* hOFrso = Hist::New("hOFrso", HistAxis(AXbta));
    Hist* hHCrso = Hist::New("hHCrso", HistAxis(AXbta));
    
    Hist* hOFrelrso = Hist::New("hOFrelrso", HistAxis(AXbta));
    Hist* hHCrelrso = Hist::New("hHCrelrso", HistAxis(AXbta));

    for (int it = 1; it <= AXbta.nbin(); ++it) {
        double bta = AXbta.center(it, AxisScale::kLog);
        double scl = 100.0;
        
        if (bta < rh_rng[0]) continue;
        if (bta > rh_rng[1]) continue;

        const Axis& AXdb = vhOF_bta.at(it)->xaxis();
        double db_width = AXdb.width(1);

        vhOF_bta.at(it)->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
        vhHC_bta.at(it)->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

        double nevt_of = (*vhOF_bta.at(it))()->Integral();
        double nevt_hc = (*vhHC_bta.at(it))()->Integral();

        double efft_of = nevt_of / nevt_of;
        double efft_hc = nevt_hc / nevt_of;
    
        (*hOFeft)()->SetBinContent(it, efft_of);
        (*hHCeft)()->SetBinContent(it, efft_hc);

        THStack* stack = Hist::Collect(Form("h_bta_%03d", it), HistList({ vhOF_bta.at(it), vhHC_bta.at(it) }));
        stack->SetTitle(Form("Velocity (%5.3f ~ %5.3f)", AXbta()(it-1), AXbta()(it)));
        //vhstack.push_back(stack);

        editor.create();
        stack->Draw("nostack hist");
        stack->GetHistogram()->SetLineColor(0);
        stack->GetHistogram()->SetMarkerColor(0);
        stack->GetHistogram()->GetXaxis()->SetTitle("#beta_{rec} / #beta_{gen} - 1");
        stack->GetHistogram()->GetYaxis()->SetTitle("Events / Bin");
        stack->Draw("nostack hist");
        Legend leg("", PadWindow(0.60, 0.80, 0.60, 0.80));
        leg()->AddEntry((TObject*)nullptr, stack->GetTitle(), "");
        leg()->AddEntry((*vhOF_bta.at(it))(), "Official Method", "lp");
        leg()->AddEntry((*vhHC_bta.at(it))(), "New", "lp");
        leg()->SetFillColor(0);
        leg.draw();
        editor.save();

        const double width = 1.0;

        (*vhOF_bta.at(it))()->Fit(func, "", "q0", -10, 10);
        (*vhOF_bta.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_bta.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_bta.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_bta.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hOFrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hOFrso)()->SetBinError  (it, scl * func->GetParError(2));
        
        (*vhOF_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hOFmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hOFmen)()->SetBinError  (it, scl * func->GetParError(1));

        (*vhHC_bta.at(it))()->Fit(func, "", "q0", -10, 10);
        (*vhHC_bta.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_bta.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_bta.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_bta.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hHCrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hHCrso)()->SetBinError  (it, scl * func->GetParError(2));
        
        (*vhHC_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hHCmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hHCmen)()->SetBinError  (it, scl * func->GetParError(1));

        double relrso_tt = (*hOFrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);
        double relrso_tq = (*hHCrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);

        (*hOFrelrso)()->SetBinContent(it, relrso_tt);
        (*hOFrelrso)()->SetBinError  (it, 0.0);
        
        (*hHCrelrso)()->SetBinContent(it, relrso_tq);
        (*hHCrelrso)()->SetBinError  (it, 0.0);
    }
    
    hHCtme->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hOFeft->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hHCeft->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hOFmen->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hHCmen->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

    hOFrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hHCrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

    hOFrelrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hHCrelrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    THStack* combined_tme = Hist::Collect("h_combind_tme", HistList({ hHCtme }));
    vhstack.push_back(combined_tme);
   
    vhist.push_back(hHCtme);

    THStack* combined_eft = Hist::Collect("h_combind_eft", HistList({ hOFeft, hHCeft }));
    vhstack.push_back(combined_eft);
    
    vhist.push_back(hOFeft);
    vhist.push_back(hHCeft);
    
    THStack* combined_men = Hist::Collect("h_combind_men", HistList({ hOFmen, hHCmen }));
    vhstack.push_back(combined_men);
    
    vhist.push_back(hOFmen);
    vhist.push_back(hHCmen);
    
    THStack* combined_rso = Hist::Collect("h_combind_rso", HistList({ hOFrso, hHCrso }));
    vhstack.push_back(combined_rso);
    
    vhist.push_back(hOFrso);
    vhist.push_back(hHCrso);
    
    THStack* combined_relrso = Hist::Collect("h_combind_relrso", HistList({ hOFrelrso, hHCrelrso }));
    vhstack.push_back(combined_relrso);

    vhist.push_back(hOFrelrso);
    vhist.push_back(hHCrelrso);

    editor.create();
    combined_tme->Draw("nostack hist");
    combined_tme->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_tme->GetHistogram()->SetLineColor(0);
    combined_tme->GetHistogram()->SetMarkerColor(0);
    combined_tme->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_tme->GetHistogram()->GetYaxis()->SetTitle("Computation Time per Event [ms]");
    combined_tme->GetHistogram()->GetXaxis()->SetRangeUser(rh_rng[0], rh_rng[1]);
    combined_tme->Draw("nostack p");
    Legend leg_tme("", PadWindow(0.60, 0.80, 0.30, 0.50));
    leg_tme()->AddEntry((*hHCtme)(), "New", "lp");
    leg_tme()->SetFillColor(0);
    leg_tme.draw();
    editor.save();
    
    editor.create();
    combined_eft->Draw("nostack hist");
    combined_eft->GetHistogram()->SetLineColor(0);
    combined_eft->GetHistogram()->SetMarkerColor(0);
    combined_eft->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_eft->GetHistogram()->GetYaxis()->SetTitle("Fitting Efficiency (Normalized by Official Method)");
    combined_eft->GetHistogram()->GetXaxis()->SetRangeUser(rh_rng[0], rh_rng[1]);
    combined_eft->Draw("nostack p");
    Legend leg_eft("", PadWindow(0.60, 0.80, 0.30, 0.40));
    leg_eft()->AddEntry((*hOFeft)(), "Official Method", "lp");
    leg_eft()->AddEntry((*hHCeft)(), "New", "lp");
    leg_eft()->SetFillColor(0);
    leg_eft.draw();
    editor.save();
    
    editor.create();
    combined_men->Draw("nostack hist");
    combined_men->GetHistogram()->SetLineColor(0);
    combined_men->GetHistogram()->SetMarkerColor(0);
    combined_men->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_men->GetHistogram()->GetYaxis()->SetTitle("Peak [%]");
    combined_men->GetHistogram()->GetXaxis()->SetRangeUser(rh_rng[0], rh_rng[1]);
    combined_men->Draw("nostack p");
    Legend leg_men("", PadWindow(0.60, 0.80, 0.20, 0.40));
    leg_men()->AddEntry((*hOFmen)(), "Official Method", "lp");
    leg_men()->AddEntry((*hHCmen)(), "New", "lp");
    leg_men()->SetFillColor(0);
    leg_men.draw();
    editor.save();
    
    editor.create();
    combined_rso->Draw("nostack hist");
    combined_rso->GetHistogram()->SetLineColor(0);
    combined_rso->GetHistogram()->SetMarkerColor(0);
    combined_rso->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution [%]");
    combined_rso->GetHistogram()->GetXaxis()->SetRangeUser(rh_rng[0], rh_rng[1]);
    combined_rso->Draw("nostack p");
    Legend leg_rso("", PadWindow(0.60, 0.80, 0.20, 0.40));
    leg_rso()->AddEntry((*hOFrso)(), "Official Method", "lp");
    leg_rso()->AddEntry((*hHCrso)(), "New", "lp");
    leg_rso()->SetFillColor(0);
    leg_rso.draw();
    editor.save();

    editor.create();
    combined_relrso->Draw("nostack hist");
    combined_relrso->GetHistogram()->SetLineColor(0);
    combined_relrso->GetHistogram()->SetMarkerColor(0);
    combined_relrso->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_relrso->GetHistogram()->GetYaxis()->SetTitle("Relative Resolution");
    combined_relrso->GetHistogram()->GetXaxis()->SetRangeUser(rh_rng[0], rh_rng[1]);
    combined_relrso->Draw("nostack p");
    Legend leg_relrso("", PadWindow(0.60, 0.80, 0.20, 0.40));
    leg_relrso()->AddEntry((*hOFrelrso)(), "Official Method", "lp");
    leg_relrso()->AddEntry((*hHCrelrso)(), "New", "lp");
    leg_relrso()->SetFillColor(0);
    leg_relrso.draw();
    editor.save();
    

    editor.close();

    TFile * ofle = new TFile("out/doc/fit_vel_rh.root", "RECREATE");
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
