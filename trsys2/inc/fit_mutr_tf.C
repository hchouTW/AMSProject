#include <CPPLibs.h>
#include <ROOTLibs.h>

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    if (argc != 2) return 0;
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/trsys/%s", argv[1]));
    double sqrm_ref = (std::string(argv[1]).find("pr") != std::string::npos) ? 0.938272297*0.938272297 : 3.727379240*3.727379240/4.0;

    PdfEditor editor(Window(), "fit_mutr_tf", "out/doc");

    TF1* func = new TF1("func", "gaus");

    const Axis& AXbta = Hist::Head("hTF_T_sqrm")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<double> tof_rng({ 0.55, 0.99 });

    Hist* hHCtme = Hist::New("hHCtme", HistAxis(AXbta));
    
    for (int it = 1; it <= AXbta.nbin(); ++it) {
        double tme_hc = (*Hist::Head("hHCmutr_TQ_numt"))()->GetBinContent(it) / (*Hist::Head("hHCmutr_TQ_dent"))()->GetBinContent(it);
        
        (*hHCtme)()->SetBinContent(it, tme_hc);
        (*hHCtme)()->SetBinError  (it, 0.0);
    }
        
    std::vector<Hist*> vhOF_sqrm = Hist::ProjectAll(HistProj::kY, Hist::Head("hTF_T_sqrm"));
    std::vector<Hist*> vhHC_sqrm = Hist::ProjectAll(HistProj::kY, Hist::Head("hHCmutr_TQ_sqrm"));
    
    Hist* hOFeft = Hist::New("hOFeft", HistAxis(AXbta));
    Hist* hHCeft = Hist::New("hHCeft", HistAxis(AXbta));

    Hist* hOFmen = Hist::New("hOFmen", HistAxis(AXbta));
    Hist* hHCmen = Hist::New("hHCmen", HistAxis(AXbta));
    
    Hist* hOFrso = Hist::New("hOFrso", HistAxis(AXbta));
    Hist* hHCrso = Hist::New("hHCrso", HistAxis(AXbta));
    
    Hist* hOFrelmen = Hist::New("hOFrelmen", HistAxis(AXbta));
    Hist* hHCrelmen = Hist::New("hHCrelmen", HistAxis(AXbta));
    
    Hist* hOFrelrso = Hist::New("hOFrelrso", HistAxis(AXbta));
    Hist* hHCrelrso = Hist::New("hHCrelrso", HistAxis(AXbta));

    for (int it = 1; it <= AXbta.nbin(); ++it) {
        double bta = AXbta.center(it, AxisScale::kLog);
        
        if (bta < tof_rng[0]) continue;
        if (bta > tof_rng[1]) continue;

        const Axis& AXdm = vhOF_sqrm.at(it)->xaxis();
        double dm_width = AXdm.width(1);

        vhOF_sqrm.at(it)->style(Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kDiamond )));
        vhHC_sqrm.at(it)->style(Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross)));

        double nevt_of = (*vhOF_sqrm.at(it))()->Integral();
        double nevt_hc = (*vhHC_sqrm.at(it))()->Integral();

        double efft_of = nevt_of / nevt_of;
        double efft_hc = nevt_hc / nevt_of;
        
        (*hOFeft)()->SetBinContent(it, efft_of);
        (*hHCeft)()->SetBinContent(it, efft_hc);

        THStack* stack = Hist::Collect(Form("h_sqrm_%03d", it), HistList({ vhOF_sqrm.at(it), vhHC_sqrm.at(it) }));
        stack->SetTitle(Form("Velocity (%4.2f ~ %4.2f)", AXbta()(it-1), AXbta()(it)));
        //vhstack.push_back(stack);

        editor.create();
        stack->Draw("nostack hist");
        stack->GetHistogram()->SetLineColor(0);
        stack->GetHistogram()->SetMarkerColor(0);
        stack->GetHistogram()->GetXaxis()->SetTitle("Mass^{2}/Z^{2} [(GeV/c^{2})^{2}]");
        stack->GetHistogram()->GetYaxis()->SetTitle("Events / Bin");
        stack->Draw("nostack hist");
        Legend leg("", PadWindow(0.60, 0.80, 0.65, 0.80));
        leg()->AddEntry((TObject*)nullptr, stack->GetTitle(), "");
        leg()->AddEntry((*vhOF_sqrm.at(it))(), "Official Method", "lp");
        leg()->AddEntry((*vhHC_sqrm.at(it))(), "New (Time, dE/dX)", "lp");
        leg()->SetFillColor(0);
        leg.draw();
        editor.save();

        const double width = 1.0;

        (*vhOF_sqrm.at(it))()->Fit(func, "", "q0", -30, 30);
        (*vhOF_sqrm.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_sqrm.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_sqrm.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_sqrm.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hOFrso)()->SetBinContent(it, func->GetParameter(2));
        (*hOFrso)()->SetBinError  (it, func->GetParError(2));

        (*vhOF_sqrm.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_sqrm.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_sqrm.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hOFmen)()->SetBinContent(it, func->GetParameter(1));
        (*hOFmen)()->SetBinError  (it, func->GetParError(1));

        (*vhHC_sqrm.at(it))()->Fit(func, "", "q0", -30, 30);
        (*vhHC_sqrm.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_sqrm.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_sqrm.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_sqrm.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hHCrso)()->SetBinContent(it, func->GetParameter(2));
        (*hHCrso)()->SetBinError  (it, func->GetParError(2));
        
        (*vhHC_sqrm.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_sqrm.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhHC_sqrm.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hHCmen)()->SetBinContent(it, func->GetParameter(1));
        (*hHCmen)()->SetBinError  (it, func->GetParError(1));
        
        double relmen_of = (*hOFmen)()->GetBinContent(it) / sqrm_ref;
        double relmen_hc = (*hHCmen)()->GetBinContent(it) / sqrm_ref;

        (*hOFrelmen)()->SetBinContent(it, relmen_of);
        (*hOFrelmen)()->SetBinError  (it, 0.0);

        (*hHCrelmen)()->SetBinContent(it, relmen_hc);
        (*hHCrelmen)()->SetBinError  (it, 0.0);
        
        double relrso_of = (*hOFrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);
        double relrso_hc = (*hHCrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);

        (*hOFrelrso)()->SetBinContent(it, relrso_of);
        (*hOFrelrso)()->SetBinError  (it, 0.0);

        (*hHCrelrso)()->SetBinContent(it, relrso_hc);
        (*hHCrelrso)()->SetBinError  (it, 0.0);
    }
    
    hHCtme->style(Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));
    
    hOFeft->style(Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kDiamond)));
    hHCeft->style(Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));
    
    hOFmen->style(Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kDiamond)));
    hHCmen->style(Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));

    hOFrso->style(Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kDiamond)));
    hHCrso->style(Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));
    
    hOFrelmen->style(Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kDiamond)));
    hHCrelmen->style(Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));

    hOFrelrso->style(Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kDiamond)));
    hHCrelrso->style(Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));
    
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
    
    THStack* combined_relmen = Hist::Collect("h_combind_relmen", HistList({ hOFrelmen, hHCrelmen }));
    vhstack.push_back(combined_relmen);
    
    vhist.push_back(hOFrelmen);
    vhist.push_back(hHCrelmen);

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
    combined_tme->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_tme->Draw("nostack p");
    Legend leg_tme("", PadWindow(0.60, 0.80, 0.70, 0.80));
    leg_tme()->AddEntry((*hHCtme)(), "New (Time, dE/dX)", "lp");
    leg_tme()->SetFillColor(0);
    leg_tme.draw();
    editor.save();
    
    editor.create();
    combined_eft->Draw("nostack hist");
    combined_eft->GetHistogram()->SetLineColor(0);
    combined_eft->GetHistogram()->SetMarkerColor(0);
    combined_eft->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_eft->GetHistogram()->GetYaxis()->SetTitle("Fitting Efficiency (Normalized by Official Method)");
    combined_eft->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_eft->Draw("nostack p");
    Legend leg_eft("", PadWindow(0.60, 0.80, 0.20, 0.40));
    leg_eft()->AddEntry((*hOFeft)(), "Official Method", "lp");
    leg_eft()->AddEntry((*hHCeft)(), "New (Time, dE/dX)", "lp");
    leg_eft()->SetFillColor(0);
    leg_eft.draw();
    editor.save();
    
    editor.create();
    combined_men->Draw("nostack hist");
    combined_men->GetHistogram()->SetLineColor(0);
    combined_men->GetHistogram()->SetMarkerColor(0);
    combined_men->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_men->GetHistogram()->GetYaxis()->SetTitle("Peak of Mass^{2}/Z^{2} [(GeV/c^{2})^{2}]");
    combined_men->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_men->Draw("nostack p");
    Legend leg_men("", PadWindow(0.20, 0.40, 0.20, 0.40));
    leg_men()->AddEntry((*hOFmen)(), "Official Method", "lp");
    leg_men()->AddEntry((*hHCmen)(), "New (Time, dE/dX)", "lp");
    leg_men()->SetFillColor(0);
    leg_men.draw();
    editor.save();
    
    editor.create();
    combined_rso->Draw("nostack hist");
    combined_rso->GetHistogram()->SetLineColor(0);
    combined_rso->GetHistogram()->SetMarkerColor(0);
    combined_rso->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Mass^{2}/Z^{2} [(GeV/c^{2})^{2}]");
    combined_rso->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_rso->Draw("nostack p");
    Legend leg_rso("", PadWindow(0.60, 0.80, 0.60, 0.80));
    leg_rso()->AddEntry((*hOFrso)(), "Official Method", "lp");
    leg_rso()->AddEntry((*hHCrso)(), "New (Time, dE/dX)", "lp");
    leg_rso()->SetFillColor(0);
    leg_rso.draw();
    editor.save();

    editor.create();
    combined_relmen->Draw("nostack hist");
    combined_relmen->GetHistogram()->SetLineColor(0);
    combined_relmen->GetHistogram()->SetMarkerColor(0);
    combined_relmen->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_relmen->GetHistogram()->GetYaxis()->SetTitle("Relative Mean");
    combined_relmen->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_relmen->Draw("nostack p");
    Legend leg_relmen("", PadWindow(0.60, 0.80, 0.20, 0.40));
    leg_relmen()->AddEntry((*hOFrelmen)(), "Official Method", "lp");
    leg_relmen()->AddEntry((*hHCrelmen)(), "New (Time, dE/dX)", "lp");
    leg_relmen()->SetFillColor(0);
    leg_relmen.draw();
    editor.save();
    
    editor.create();
    combined_relrso->Draw("nostack hist");
    combined_relrso->GetHistogram()->SetLineColor(0);
    combined_relrso->GetHistogram()->SetMarkerColor(0);
    combined_relrso->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_relrso->GetHistogram()->GetYaxis()->SetTitle("Relative Resolution");
    combined_relrso->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_relrso->Draw("nostack p");
    Legend leg_relrso("", PadWindow(0.60, 0.80, 0.20, 0.40));
    leg_relrso()->AddEntry((*hOFrelrso)(), "Official Method", "lp");
    leg_relrso()->AddEntry((*hHCrelrso)(), "New (Time, dE/dX)", "lp");
    leg_relrso()->SetFillColor(0);
    leg_relrso.draw();
    editor.save();

    editor.close();

    TFile * ofle = new TFile("out/doc/fit_mutr_tf.root", "RECREATE");
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
