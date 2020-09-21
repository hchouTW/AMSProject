#include <CPPLibs.h>
#include <ROOTLibs.h>

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    if (argc != 2) return 0;
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/trsys/%s", argv[1]));

    PdfEditor editor(Window(), "fit_vel_tf", "out/doc");

    TF1* func = new TF1("func", "gaus");

    const Axis& AXbta = Hist::Head("hTF_T_bta")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<double> tof_rng({ 0.55, 0.99 });

    Hist* hTTtme = Hist::New("hTTtme", HistAxis(AXbta));
    Hist* hTQtme = Hist::New("hTQtme", HistAxis(AXbta));
    
    for (int it = 1; it <= AXbta.nbin(); ++it) {
        double tme_tt = (*Hist::Head("hHC_T_numt"))()->GetBinContent(it) / (*Hist::Head("hHC_T_dent"))()->GetBinContent(it);
        double tme_tq = (*Hist::Head("hHC_TQ_numt"))()->GetBinContent(it) / (*Hist::Head("hHC_TQ_dent"))()->GetBinContent(it);
        
        (*hTTtme)()->SetBinContent(it, tme_tt);
        (*hTTtme)()->SetBinError  (it, 0.0);
        (*hTQtme)()->SetBinContent(it, tme_tq);
        (*hTQtme)()->SetBinError  (it, 0.0);
    }
        
    std::vector<Hist*> vhOF_bta = Hist::ProjectAll(HistProj::kY, Hist::Head("hTF_T_bta"));
    std::vector<Hist*> vhTT_bta = Hist::ProjectAll(HistProj::kY, Hist::Head("hHC_T_bta"));
    std::vector<Hist*> vhTQ_bta = Hist::ProjectAll(HistProj::kY, Hist::Head("hHC_TQ_bta"));
    
    Hist* hOFeft = Hist::New("hOFeft", HistAxis(AXbta));
    Hist* hTTeft = Hist::New("hTTeft", HistAxis(AXbta));
    Hist* hTQeft = Hist::New("hTQeft", HistAxis(AXbta));

    Hist* hOFmen = Hist::New("hOFmen", HistAxis(AXbta));
    Hist* hTTmen = Hist::New("hTTmen", HistAxis(AXbta));
    Hist* hTQmen = Hist::New("hTQmen", HistAxis(AXbta));
    
    Hist* hOFrso = Hist::New("hOFrso", HistAxis(AXbta));
    Hist* hTTrso = Hist::New("hTTrso", HistAxis(AXbta));
    Hist* hTQrso = Hist::New("hTQrso", HistAxis(AXbta));
    
    Hist* hOFrelrso = Hist::New("hOFrelrso", HistAxis(AXbta));
    Hist* hTTrelrso = Hist::New("hTTrelrso", HistAxis(AXbta));
    Hist* hTQrelrso = Hist::New("hTQrelrso", HistAxis(AXbta));

    for (int it = 1; it <= AXbta.nbin(); ++it) {
        double bta = AXbta.center(it, AxisScale::kLinear);
        double scl = 100.0;
        
        if (bta < tof_rng[0]) continue;
        if (bta > tof_rng[1]) continue;

        const Axis& AXdb = vhTT_bta.at(it)->xaxis();
        double db_width = AXdb.width(1);

        vhOF_bta.at(it)->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
        vhTT_bta.at(it)->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
        vhTQ_bta.at(it)->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

        double nevt_of = (*vhOF_bta.at(it))()->Integral();
        double nevt_tt = (*vhTT_bta.at(it))()->Integral();
        double nevt_tq = (*vhTQ_bta.at(it))()->Integral();

        double efft_of = nevt_of / nevt_of;
        double efft_tt = nevt_tt / nevt_of;
        double efft_tq = nevt_tq / nevt_of;
        
        (*hOFeft)()->SetBinContent(it, efft_of);
        (*hTTeft)()->SetBinContent(it, efft_tt);
        (*hTQeft)()->SetBinContent(it, efft_tq);

        THStack* stack = Hist::Collect(Form("h_bta_%03d", it), HistList({ vhOF_bta.at(it), vhTT_bta.at(it), vhTQ_bta.at(it) }));
        stack->SetTitle(Form("Velocity (%4.2f ~ %4.2f)", AXbta()(it-1), AXbta()(it)));
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
        leg()->AddEntry((*vhTT_bta.at(it))(), "New (Time)", "lp");
        leg()->AddEntry((*vhTQ_bta.at(it))(), "New (Time, dE/dX)", "lp");
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

        (*vhTT_bta.at(it))()->Fit(func, "", "q0", -10, 10);
        (*vhTT_bta.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTT_bta.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTT_bta.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTT_bta.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hTTrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hTTrso)()->SetBinError  (it, scl * func->GetParError(2));
        
        (*vhTT_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTT_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTT_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hTTmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hTTmen)()->SetBinError  (it, scl * func->GetParError(1));

        (*vhTQ_bta.at(it))()->Fit(func, "", "q0", -10, 10);
        (*vhTQ_bta.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTQ_bta.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTQ_bta.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTQ_bta.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hTQrso)()->SetBinContent(it, scl * func->GetParameter(2));
        (*hTQrso)()->SetBinError  (it, scl * func->GetParError(2));
        
        (*vhTQ_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTQ_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*vhTQ_bta.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
        (*hTQmen)()->SetBinContent(it, scl * func->GetParameter(1));
        (*hTQmen)()->SetBinError  (it, scl * func->GetParError(1));

        double relrso_of = (*hOFrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);
        double relrso_tt = (*hTTrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);
        double relrso_tq = (*hTQrso)()->GetBinContent(it) / (*hOFrso)()->GetBinContent(it);

        (*hOFrelrso)()->SetBinContent(it, relrso_of);
        (*hOFrelrso)()->SetBinError  (it, 0.0);

        (*hTTrelrso)()->SetBinContent(it, relrso_tt);
        (*hTTrelrso)()->SetBinError  (it, 0.0);
        
        (*hTQrelrso)()->SetBinContent(it, relrso_tq);
        (*hTQrelrso)()->SetBinError  (it, 0.0);
    }
    
    hTTtme->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hTQtme->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hOFeft->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
    hTTeft->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hTQeft->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    hOFmen->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
    hTTmen->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hTQmen->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

    hOFrso->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
    hTTrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hTQrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

    hOFrelrso->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
    hTTrelrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
    hTQrelrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
    
    THStack* combined_tme = Hist::Collect("h_combind_tme", HistList({ hTTtme, hTQtme }));
    vhstack.push_back(combined_tme);

    vhist.push_back(hTTtme);
    vhist.push_back(hTQtme);

    THStack* combined_eft = Hist::Collect("h_combind_eft", HistList({ hOFeft, hTTeft, hTQeft }));
    vhstack.push_back(combined_eft);
    
    vhist.push_back(hOFeft);
    vhist.push_back(hTTeft);
    vhist.push_back(hTQeft);
    
    THStack* combined_men = Hist::Collect("h_combind_men", HistList({ hOFmen, hTTmen, hTQmen }));
    vhstack.push_back(combined_men);
    
    vhist.push_back(hOFmen);
    vhist.push_back(hTTmen);
    vhist.push_back(hTQmen);
    
    THStack* combined_rso = Hist::Collect("h_combind_rso", HistList({ hOFrso, hTTrso, hTQrso }));
    vhstack.push_back(combined_rso);
    
    vhist.push_back(hOFrso);
    vhist.push_back(hTTrso);
    vhist.push_back(hTQrso);
    
    THStack* combined_relrso = Hist::Collect("h_combind_relrso", HistList({ hOFrelrso, hTTrelrso, hTQrelrso }));
    vhstack.push_back(combined_relrso);
    
    vhist.push_back(hOFrelrso);
    vhist.push_back(hTTrelrso);
    vhist.push_back(hTQrelrso);

    editor.create();
    combined_tme->Draw("nostack hist");
    combined_tme->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    combined_tme->GetHistogram()->SetLineColor(0);
    combined_tme->GetHistogram()->SetMarkerColor(0);
    combined_tme->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_tme->GetHistogram()->GetYaxis()->SetTitle("Computation Time per Event [ms]");
    combined_tme->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_tme->Draw("nostack p");
    Legend leg_tme("", PadWindow(0.60, 0.80, 0.30, 0.50));
    leg_tme()->AddEntry((*hTTtme)(), "New (Time)", "lp");
    leg_tme()->AddEntry((*hTQtme)(), "New (Time, dE/dX)", "lp");
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
    leg_eft()->AddEntry((*hTTeft)(), "New (Time)", "lp");
    leg_eft()->AddEntry((*hTQeft)(), "New (Time, dE/dX)", "lp");
    leg_eft()->SetFillColor(0);
    leg_eft.draw();
    editor.save();
    
    editor.create();
    combined_men->Draw("nostack hist");
    combined_men->GetHistogram()->SetLineColor(0);
    combined_men->GetHistogram()->SetMarkerColor(0);
    combined_men->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_men->GetHistogram()->GetYaxis()->SetTitle("Peak [%]");
    combined_men->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_men->Draw("nostack p");
    Legend leg_men("", PadWindow(0.60, 0.80, 0.20, 0.40));
    leg_men()->AddEntry((*hOFmen)(), "Official Method", "lp");
    leg_men()->AddEntry((*hTTmen)(), "New (Time)", "lp");
    leg_men()->AddEntry((*hTQmen)(), "New (Time, dE/dX)", "lp");
    leg_men()->SetFillColor(0);
    leg_men.draw();
    editor.save();
    
    editor.create();
    combined_rso->Draw("nostack hist");
    combined_rso->GetHistogram()->SetLineColor(0);
    combined_rso->GetHistogram()->SetMarkerColor(0);
    combined_rso->GetHistogram()->GetXaxis()->SetTitle("Velocity");
    combined_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution [%]");
    combined_rso->GetHistogram()->GetXaxis()->SetRangeUser(tof_rng[0], tof_rng[1]);
    combined_rso->Draw("nostack p");
    Legend leg_rso("", PadWindow(0.60, 0.80, 0.20, 0.40));
    leg_rso()->AddEntry((*hOFrso)(), "Official Method", "lp");
    leg_rso()->AddEntry((*hTTrso)(), "New (Time)", "lp");
    leg_rso()->AddEntry((*hTQrso)(), "New (Time, dE/dX)", "lp");
    leg_rso()->SetFillColor(0);
    leg_rso.draw();
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
    leg_relrso()->AddEntry((*hTTrelrso)(), "New (Time)", "lp");
    leg_relrso()->AddEntry((*hTQrelrso)(), "New (Time, dE/dX)", "lp");
    leg_relrso()->SetFillColor(0);
    leg_relrso.draw();
    editor.save();
    

    editor.close();

    TFile * ofle = new TFile("out/doc/fit_vel_tf.root", "RECREATE");
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
