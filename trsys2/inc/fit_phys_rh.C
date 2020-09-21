#include <CPPLibs.h>
#include <ROOTLibs.h>

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    if (argc != 2) return 0;
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/trsys/%s", argv[1]));
    int chrg = (std::string(argv[1]).find("pr") != std::string::npos) ? 1 : 2;

    PdfEditor editor(Window(), "fit_phys_rh", "out/doc");

    TF1* func = new TF1("func", "gaus");

    const Axis& AXrig = Hist::Head("hHC_inrh_rig")->xaxis();

    std::array<double, 2> trpt_rng({ 3.3 * chrg, 40.0 * chrg });

    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    Hist* hGMtme = Hist::New("hGMtme", HistAxis(AXrig));
    Hist* hPHtme = Hist::New("hPHtme", HistAxis(AXrig));
    
    for (int it = 1; it <= AXrig.nbin(); ++it) {
        double tme_gm = (*Hist::Head("hHC_inrh_numt"))()->GetBinContent(it) / (*Hist::Head("hHC_inrh_dent"))()->GetBinContent(it);
        double tme_ph = (*Hist::Head("hHCphys_B_numt"))()->GetBinContent(it) / (*Hist::Head("hHCphys_B_dent"))()->GetBinContent(it);
        
        (*hGMtme)()->SetBinContent(it, tme_gm);
        (*hGMtme)()->SetBinError  (it, 0.0);

        (*hPHtme)()->SetBinContent(it, tme_ph);
        (*hPHtme)()->SetBinError  (it, 0.0);
    }
        
    std::vector<Hist*> vhOF_rig = Hist::ProjectAll(HistProj::kY, Hist::Head("hCK_inrh_rig"));
    std::vector<Hist*> vhGM_rig = Hist::ProjectAll(HistProj::kY, Hist::Head("hHC_inrh_rig"));
    std::vector<Hist*> vhPH_rig = Hist::ProjectAll(HistProj::kY, Hist::Head("hHCphys_B_rig"));
    
    Hist* hOFeft = Hist::New("hOFeft", HistAxis(AXrig));
    Hist* hGMeft = Hist::New("hGMeft", HistAxis(AXrig));
    Hist* hPHeft = Hist::New("hPHeft", HistAxis(AXrig));

    Hist* hOFmen = Hist::New("hOFmen", HistAxis(AXrig));
    Hist* hGMmen = Hist::New("hGMmen", HistAxis(AXrig));
    Hist* hPHmen = Hist::New("hPHmen", HistAxis(AXrig));
    
    Hist* hOFrso = Hist::New("hOFrso", HistAxis(AXrig));
    Hist* hGMrso = Hist::New("hGMrso", HistAxis(AXrig));
    Hist* hPHrso = Hist::New("hPHrso", HistAxis(AXrig));
    
    Hist* hOFrelrso = Hist::New("hOFrelrso", HistAxis(AXrig));
    Hist* hGMrelrso = Hist::New("hGMrelrso", HistAxis(AXrig));
    Hist* hPHrelrso = Hist::New("hPHrelrso", HistAxis(AXrig));

    for (int it = 1; it <= AXrig.nbin(); ++it) {
        double mom = AXrig.center(it, AxisScale::kLog);
        double scl = std::sqrt(mom) * 100.0;

        if (mom < trpt_rng[0]) continue;
        if (mom > trpt_rng[1]) continue;

        const Axis& AXdr = vhGM_rig.at(it)->xaxis();
        double dr_width = AXdr.width(1);

        vhOF_rig.at(it)->style(Line(kBlack), Marker(kBlack, MarkerStyle(MarkerShape::kCircle)));
        vhGM_rig.at(it)->style(Line(kBlue),  Marker(kBlue, MarkerStyle(MarkerShape::kDiamond)));
        vhPH_rig.at(it)->style(Line(kRed),   Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));

        double nevt_of = (*vhOF_rig.at(it))()->Integral();
        double nevt_gm = (*vhGM_rig.at(it))()->Integral();
        double nevt_ph = (*vhPH_rig.at(it))()->Integral();

        double efft_of = nevt_gm / nevt_of;
        double efft_gm = nevt_gm / nevt_of;
        double efft_ph = nevt_ph / nevt_of;
        
        (*hOFeft)()->SetBinContent(it, efft_of);
        (*hGMeft)()->SetBinContent(it, efft_gm);
        (*hPHeft)()->SetBinContent(it, efft_ph);

        THStack* stack = Hist::Collect(Form("h_rig_%03d", it), HistList({ vhOF_rig.at(it), vhGM_rig.at(it), vhPH_rig.at(it) }));
        stack->SetTitle(Form("Rigidity (%6.2f ~ %6.2f)", AXrig()(it-1), AXrig()(it)));
        //vhstack.push_back(stack);
        //stack->Write();

        editor.create();
        editor()().SetLogx(0);
        stack->Draw("nostack hist");
        stack->GetHistogram()->SetLineColor(0);
        stack->GetHistogram()->SetMarkerColor(0);
        stack->GetHistogram()->GetXaxis()->SetTitle("R_{gen}^{1/2} * (1/R_{rec} - 1/R_{gen})");
        stack->GetHistogram()->GetYaxis()->SetTitle("Events / Bin");
        stack->Draw("nostack hist");
        Legend leg("<< Inner >>", PadWindow(0.60, 0.80, 0.65, 0.80));
        leg()->AddEntry((TObject*)nullptr, stack->GetTitle(), "");
        leg()->AddEntry((*vhOF_rig.at(it))(), "Official Fitting", "lp");
        leg()->AddEntry((*vhGM_rig.at(it))(), "Rigidity Fitting", "lp");
        leg()->AddEntry((*vhPH_rig.at(it))(), "Track Fitting", "lp");
        leg()->SetFillColor(0);
        leg.draw();
        editor.save();

        const double width = 1.0;

        (*vhOF_rig.at(it))()->Fit(func, "", "q0", -10, 10);
        (*vhOF_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hOFrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hOFrso)()->SetBinError  (it, scl * func->GetParError(2));

        (*vhOF_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhOF_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hOFmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hOFmen)()->SetBinError  (it, scl * func->GetParError(1));

        (*vhGM_rig.at(it))()->Fit(func, "", "q0", -10, 10);
        (*vhGM_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhGM_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhGM_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhGM_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hGMrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hGMrso)()->SetBinError  (it, scl * func->GetParError(2));

        (*vhGM_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhGM_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhGM_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hGMmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hGMmen)()->SetBinError  (it, scl * func->GetParError(1));

        (*vhPH_rig.at(it))()->Fit(func, "", "q0", -10, 10);
        (*vhPH_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhPH_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhPH_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhPH_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hPHrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hPHrso)()->SetBinError  (it, scl * func->GetParError(2));
        
        (*vhPH_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhPH_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhPH_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hPHmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hPHmen)()->SetBinError  (it, scl * func->GetParError(1));

        double relrso_of = (*hOFrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);
        double relrso_gm = (*hGMrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);
        double relrso_ph = (*hPHrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);

        (*hOFrelrso)()->SetBinContent(it, relrso_of);
        (*hOFrelrso)()->SetBinError  (it, 0.0);
        
        (*hGMrelrso)()->SetBinContent(it, relrso_gm);
        (*hGMrelrso)()->SetBinError  (it, 0.0);
        
        (*hPHrelrso)()->SetBinContent(it, relrso_ph);
        (*hPHrelrso)()->SetBinError  (it, 0.0);
    }
    
    hGMtme->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hPHtme->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hGMeft->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle)));
    hGMeft->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hPHeft->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hGMmen->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle)));
    hGMmen->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hPHmen->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

    hGMrso->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle)));
    hGMrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hPHrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

    hGMrelrso->style(Line(kBlack),   Marker(kBlack,   MarkerStyle(MarkerShape::kCircle)));
    hGMrelrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hPHrelrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    THStack* combined_tme = Hist::Collect("h_combind_tme", HistList({ hGMtme, hPHtme }));
    vhstack.push_back(combined_tme);
    
    vhist.push_back(hGMtme);
    vhist.push_back(hPHtme);

    THStack* combined_eft = Hist::Collect("h_combind_eft", HistList({ hOFeft, hGMeft, hPHeft }));
    vhstack.push_back(combined_eft);
    
    vhist.push_back(hOFeft);
    vhist.push_back(hGMeft);
    vhist.push_back(hPHeft);
    
    THStack* combined_men = Hist::Collect("h_combind_men", HistList({ hOFmen, hGMmen, hPHmen }));
    vhstack.push_back(combined_men);
    
    vhist.push_back(hOFmen);
    vhist.push_back(hGMmen);
    vhist.push_back(hPHmen);
    
    THStack* combined_rso = Hist::Collect("h_combind_rso", HistList({ hOFrso, hGMrso, hPHrso }));
    vhstack.push_back(combined_rso);
    
    vhist.push_back(hOFrso);
    vhist.push_back(hGMrso);
    vhist.push_back(hPHrso);
    
    THStack* combined_relrso = Hist::Collect("h_combind_relrso", HistList({ hOFrelrso, hGMrelrso, hPHrelrso }));
    vhstack.push_back(combined_relrso);
    
    vhist.push_back(hOFrelrso);
    vhist.push_back(hGMrelrso);
    vhist.push_back(hPHrelrso);

    editor.create();
    editor()().SetLogx();
    combined_tme->Draw("nostack hist");
    combined_tme->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_tme->GetHistogram()->SetLineColor(0);
    combined_tme->GetHistogram()->SetMarkerColor(0);
    combined_tme->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GeV/c]");
    combined_tme->GetHistogram()->GetYaxis()->SetTitle("Computation Time per Event [ms]");
    combined_tme->GetHistogram()->GetXaxis()->SetRangeUser(trpt_rng[0], trpt_rng[1]);
    combined_tme->Draw("nostack p");
    Legend leg_tme("<< Inner >>", PadWindow(0.60, 0.80, 0.65, 0.80));
    leg_tme()->AddEntry((*hGMtme)(), "Rigidity Fitting", "lp");
    leg_tme()->AddEntry((*hPHtme)(), "Track Fitting", "lp");
    leg_tme()->SetFillColor(0);
    leg_tme.draw();
    editor.save();
    
    editor.create();
    editor()().SetLogx();
    combined_eft->Draw("nostack hist");
    combined_eft->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_eft->GetHistogram()->SetLineColor(0);
    combined_eft->GetHistogram()->SetMarkerColor(0);
    combined_eft->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GeV/c]");
    combined_eft->GetHistogram()->GetYaxis()->SetTitle("Fitting Efficiency (Normalized by Rigidity Fitting)");
    combined_eft->GetHistogram()->GetXaxis()->SetRangeUser(trpt_rng[0], trpt_rng[1]);
    combined_eft->Draw("nostack p");
    Legend leg_eft("<< Inner >>", PadWindow(0.60, 0.80, 0.20, 0.35));
    leg_eft()->AddEntry((*hOFeft)(), "Official Fitting", "lp");
    leg_eft()->AddEntry((*hGMeft)(), "Rigidity Fitting", "lp");
    leg_eft()->AddEntry((*hPHeft)(), "Track Fitting", "lp");
    leg_eft()->SetFillColor(0);
    leg_eft.draw();
    editor.save();
    
    editor.create();
    editor()().SetLogx();
    combined_men->Draw("nostack hist");
    combined_men->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_men->GetHistogram()->SetLineColor(0);
    combined_men->GetHistogram()->SetMarkerColor(0);
    combined_men->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GeV/c]");
    combined_men->GetHistogram()->GetYaxis()->SetTitle("Peak [%]");
    combined_men->GetHistogram()->GetXaxis()->SetRangeUser(trpt_rng[0], trpt_rng[1]);
    combined_men->Draw("nostack p");
    Legend leg_men("<< Inner >>", PadWindow(0.60, 0.80, 0.20, 0.35));
    leg_men()->AddEntry((*hOFmen)(), "Official Fitting", "lp");
    leg_men()->AddEntry((*hGMmen)(), "Rigidity Fitting", "lp");
    leg_men()->AddEntry((*hPHmen)(), "Track Fitting", "lp");
    leg_men()->SetFillColor(0);
    leg_men.draw();
    editor.save();
    
    editor.create();
    editor()().SetLogx();
    combined_rso->Draw("nostack hist");
    combined_rso->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_rso->GetHistogram()->SetLineColor(0);
    combined_rso->GetHistogram()->SetMarkerColor(0);
    combined_rso->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GeV/c]");
    combined_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution [%]");
    combined_rso->GetHistogram()->GetXaxis()->SetRangeUser(trpt_rng[0], trpt_rng[1]);
    combined_rso->Draw("nostack p");
    Legend leg_rso("<< Inner >>", PadWindow(0.60, 0.80, 0.20, 0.35));
    leg_rso()->AddEntry((*hOFrso)(), "Official Fitting", "lp");
    leg_rso()->AddEntry((*hGMrso)(), "Rigidity Fitting", "lp");
    leg_rso()->AddEntry((*hPHrso)(), "Track Fitting", "lp");
    leg_rso()->SetFillColor(0);
    leg_rso.draw();
    editor.save();

    editor.create();
    editor()().SetLogx();
    combined_relrso->Draw("nostack hist");
    combined_relrso->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_relrso->GetHistogram()->SetLineColor(0);
    combined_relrso->GetHistogram()->SetMarkerColor(0);
    combined_relrso->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GeV/c]");
    combined_relrso->GetHistogram()->GetYaxis()->SetTitle("Relative Resolution");
    combined_relrso->GetHistogram()->GetXaxis()->SetRangeUser(trpt_rng[0], trpt_rng[1]);
    combined_relrso->Draw("nostack p");
    Legend leg_relrso("<< Inner >>", PadWindow(0.60, 0.80, 0.20, 0.35));
    leg_relrso()->AddEntry((*hOFrelrso)(), "Official Fitting", "lp");
    leg_relrso()->AddEntry((*hGMrelrso)(), "Rigidity Fitting", "lp");
    leg_relrso()->AddEntry((*hPHrelrso)(), "Track Fitting", "lp");
    leg_relrso()->SetFillColor(0);
    leg_relrso.draw();
    editor.save();
    
    editor.close();

    TFile * ofle = new TFile("out/doc/fit_phys_rh.root", "RECREATE");
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
