#include <CPPLibs.h>
#include <ROOTLibs.h>

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    if (argc != 2) return 0;
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/trsys/%s", argv[1]));
    int chrg = (std::string(argv[1]).find("pr") != std::string::npos) ? 1 : 2;

    PdfEditor editor(Window(), "fit_geom", "out/doc");

    TF1* func = new TF1("func", "gaus");

    const Axis& AXrig = Hist::Head("hCK_in_rig")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::array<std::string, 4> trpts_name({"<< Inner >>", "<< InnerL1 >>", "<< InnerL9 >>", "<< FullSpan >>"});
    std::array<std::string, 4> trpts({"in", "l1", "l9", "fs"});
    
    std::vector<std::vector<double>> trpts_rng({ { 0.6 * chrg, 400.0 * chrg }, { 0.6 * chrg, 1200.0 * chrg }, { 0.6 * chrg, 1800.0 * chrg }, { 0.6 * chrg, 2500.0 * chrg } });

    for (int ip = 0; ip < trpts.size(); ++ip) {
        Hist* hKFtme = Hist::New(Form("hKFtme_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hCKtme = Hist::New(Form("hCKtme_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hHCtme = Hist::New(Form("hHCtme_%s", trpts[ip].c_str()), HistAxis(AXrig));
        
        for (int it = 1; it <= AXrig.nbin(); ++it) {
            double tme_kf = (*Hist::Head(Form("hKF_%s_numt", trpts[ip].c_str())))()->GetBinContent(it) / (*Hist::Head(Form("hKF_%s_dent", trpts[ip].c_str())))()->GetBinContent(it);
            double tme_ck = (*Hist::Head(Form("hCK_%s_numt", trpts[ip].c_str())))()->GetBinContent(it) / (*Hist::Head(Form("hCK_%s_dent", trpts[ip].c_str())))()->GetBinContent(it);
            double tme_hc = (*Hist::Head(Form("hHC_%s_numt", trpts[ip].c_str())))()->GetBinContent(it) / (*Hist::Head(Form("hHC_%s_dent", trpts[ip].c_str())))()->GetBinContent(it);
            
            (*hKFtme)()->SetBinContent(it, tme_kf);
            (*hKFtme)()->SetBinError  (it, 0.0);
            
            (*hCKtme)()->SetBinContent(it, tme_ck);
            (*hCKtme)()->SetBinError  (it, 0.0);

            (*hHCtme)()->SetBinContent(it, tme_hc);
            (*hHCtme)()->SetBinError  (it, 0.0);
        }
            
        std::vector<Hist*> vhKF_rig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hKF_%s_rig", trpts[ip].c_str())));
        std::vector<Hist*> vhCK_rig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hCK_%s_rig", trpts[ip].c_str())));
        std::vector<Hist*> vhHC_rig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHC_%s_rig", trpts[ip].c_str())));
        
        Hist* hKFeft = Hist::New(Form("hKFeft_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hCKeft = Hist::New(Form("hCKeft_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hHCeft = Hist::New(Form("hHCeft_%s", trpts[ip].c_str()), HistAxis(AXrig));

        Hist* hKFmen = Hist::New(Form("hKFmen_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hCKmen = Hist::New(Form("hCKmen_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hHCmen = Hist::New(Form("hHCmen_%s", trpts[ip].c_str()), HistAxis(AXrig));
        
        Hist* hKFrso = Hist::New(Form("hKFrso_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hCKrso = Hist::New(Form("hCKrso_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hHCrso = Hist::New(Form("hHCrso_%s", trpts[ip].c_str()), HistAxis(AXrig));
        
        Hist* hKFrelrso = Hist::New(Form("hKFrelrso_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hCKrelrso = Hist::New(Form("hCKrelrso_%s", trpts[ip].c_str()), HistAxis(AXrig));
        Hist* hHCrelrso = Hist::New(Form("hHCrelrso_%s", trpts[ip].c_str()), HistAxis(AXrig));

        for (int it = 1; it <= AXrig.nbin(); ++it) {
            double mom = AXrig.center(it, AxisScale::kLog);
            double scl = std::sqrt(mom) * 100.0;
        
            if (mom < trpts_rng[ip][0]) continue;
            if (mom > trpts_rng[ip][1]) continue;

            const Axis& AXdr = vhCK_rig.at(it)->xaxis();
            double dr_width = AXdr.width(1);

            vhKF_rig.at(it)->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
            vhCK_rig.at(it)->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
            vhHC_rig.at(it)->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

            double nevt_kf = (*vhKF_rig.at(it))()->Integral();
            double nevt_ck = (*vhCK_rig.at(it))()->Integral();
            double nevt_hc = (*vhHC_rig.at(it))()->Integral();

            double efft_kf = nevt_kf / nevt_ck;
            double efft_ck = nevt_ck / nevt_ck;
            double efft_hc = nevt_hc / nevt_ck;
            
            (*hKFeft)()->SetBinContent(it, efft_kf);
            (*hCKeft)()->SetBinContent(it, efft_ck);
            (*hHCeft)()->SetBinContent(it, efft_hc);

            THStack* stack = Hist::Collect(Form("h_%s_rig_%03d", trpts[ip].c_str(), it), HistList({ vhKF_rig.at(it), vhCK_rig.at(it), vhHC_rig.at(it) }));
            stack->SetTitle(Form("Rigidity (%6.2f ~ %6.2f)", AXrig()(it-1), AXrig()(it)));
            //vhstack.push_back(stack);

            editor.create();
            editor()().SetLogx(0);
            stack->Draw("nostack hist");
            stack->GetHistogram()->SetLineColor(0);
            stack->GetHistogram()->SetMarkerColor(0);
            stack->GetHistogram()->GetXaxis()->SetTitle("R_{gen}^{1/2} * (1/R_{rec} - 1/R_{gen})");
            stack->GetHistogram()->GetYaxis()->SetTitle("Events / Bin");
            stack->GetHistogram()->SetTitle(trpts[ip].c_str());
            stack->Draw("nostack hist");
            Legend leg(trpts_name[ip].c_str(), PadWindow(0.60, 0.80, 0.60, 0.80));
            leg()->AddEntry((TObject*)nullptr, stack->GetTitle(), "");
            leg()->AddEntry((*vhKF_rig.at(it))(), "Kalman Fitter", "lp");
            leg()->AddEntry((*vhCK_rig.at(it))(), "Choutko", "lp");
            leg()->AddEntry((*vhHC_rig.at(it))(), "New", "lp");
            leg()->SetFillColor(0);
            leg.draw();
            editor.save();

            const double width = 1.0;

            (*vhKF_rig.at(it))()->Fit(func, "", "q0", -10, 10);
            (*vhKF_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhKF_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhKF_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhKF_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*hKFrso)()->SetBinContent(it, scl * func->GetParameter(2));
            (*hKFrso)()->SetBinError  (it, scl * func->GetParError(2));
            
            (*vhKF_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhKF_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhKF_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*hKFmen)()->SetBinContent(it, scl * func->GetParameter(1));
            (*hKFmen)()->SetBinError  (it, scl * func->GetParError(1));

            (*vhCK_rig.at(it))()->Fit(func, "", "q0", -10, 10);
            (*vhCK_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhCK_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhCK_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhCK_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*hCKrso)()->SetBinContent(it, scl * func->GetParameter(2));
            (*hCKrso)()->SetBinError  (it, scl * func->GetParError(2));
            
            (*vhCK_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhCK_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhCK_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*hCKmen)()->SetBinContent(it, scl * func->GetParameter(1));
            (*hCKmen)()->SetBinError  (it, scl * func->GetParError(1));

            (*vhHC_rig.at(it))()->Fit(func, "", "q0", -10, 10);
            (*vhHC_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhHC_rig.at(it))()->Fit(func, "", "q0", -3.0 * width * func->GetParameter(2) + func->GetParameter(1), 3.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhHC_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhHC_rig.at(it))()->Fit(func, "", "q0", -2.0 * width * func->GetParameter(2) + func->GetParameter(1), 2.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*hHCrso)()->SetBinContent(it, scl * func->GetParameter(2));
            (*hHCrso)()->SetBinError  (it, scl * func->GetParError(2));
            
            (*vhHC_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhHC_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*vhHC_rig.at(it))()->Fit(func, "", "q0", -1.0 * width * func->GetParameter(2) + func->GetParameter(1), 1.0 * width * func->GetParameter(2) + func->GetParameter(1));
            (*hHCmen)()->SetBinContent(it, scl * func->GetParameter(1));
            (*hHCmen)()->SetBinError  (it, scl * func->GetParError(1));

            double relrso_kf = (*hKFrso)()->GetBinContent(it) / (*hCKrso)()->GetBinContent(it);
            double relrso_ck = (*hCKrso)()->GetBinContent(it) / (*hCKrso)()->GetBinContent(it);
            double relrso_hc = (*hHCrso)()->GetBinContent(it) / (*hCKrso)()->GetBinContent(it);

            (*hKFrelrso)()->SetBinContent(it, relrso_kf);
            (*hKFrelrso)()->SetBinError  (it, 0.0);

            (*hCKrelrso)()->SetBinContent(it, relrso_ck);
            (*hCKrelrso)()->SetBinError  (it, 0.0);
            
            (*hHCrelrso)()->SetBinContent(it, relrso_hc);
            (*hHCrelrso)()->SetBinError  (it, 0.0);
        }
        
        hKFtme->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
        hCKtme->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
        hHCtme->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
        
        hKFeft->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
        hCKeft->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
        hHCeft->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
        
        hKFmen->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
        hCKmen->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
        hHCmen->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
            
        hKFrso->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
        hCKrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
        hHCrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));

        hKFrelrso->style(Line(kGreen+1), Marker(kGreen+1, MarkerStyle(MarkerShape::kCircle )));
        hCKrelrso->style(Line(kBlue),    Marker(kBlue,    MarkerStyle(MarkerShape::kDiamond)));
        hHCrelrso->style(Line(kRed),     Marker(kRed,     MarkerStyle(MarkerShape::kCross  )));
        
        THStack* combined_tme = Hist::Collect(Form("h_combind_tme_%s", trpts[ip].c_str()), HistList({ hKFtme, hCKtme, hHCtme }));
        vhstack.push_back(combined_tme);
       
        vhist.push_back(hKFtme);
        vhist.push_back(hCKtme);
        vhist.push_back(hHCtme);

        THStack* combined_eft = Hist::Collect(Form("h_combind_eft_%s", trpts[ip].c_str()), HistList({ hKFeft, hCKeft, hHCeft }));
        vhstack.push_back(combined_eft);
        
        vhist.push_back(hKFeft);
        vhist.push_back(hCKeft);
        vhist.push_back(hHCeft);
        
        THStack* combined_men = Hist::Collect(Form("h_combind_men_%s", trpts[ip].c_str()), HistList({ hKFmen, hCKmen, hHCmen }));
        vhstack.push_back(combined_men);
        
        vhist.push_back(hKFmen);
        vhist.push_back(hCKmen);
        vhist.push_back(hHCmen);
        
        THStack* combined_rso = Hist::Collect(Form("h_combind_rso_%s", trpts[ip].c_str()), HistList({ hKFrso, hCKrso, hHCrso }));
        vhstack.push_back(combined_rso);
        
        vhist.push_back(hKFrso);
        vhist.push_back(hCKrso);
        vhist.push_back(hHCrso);
        
        THStack* combined_relrso = Hist::Collect(Form("h_combind_relrso_%s", trpts[ip].c_str()), HistList({ hKFrelrso, hCKrelrso, hHCrelrso }));
        vhstack.push_back(combined_relrso);
        
        vhist.push_back(hKFrelrso);
        vhist.push_back(hCKrelrso);
        vhist.push_back(hHCrelrso);

        editor.create();
        editor()().SetLogx();
        combined_tme->Draw("nostack hist");
        combined_tme->GetHistogram()->GetXaxis()->SetMoreLogLabels();
        combined_tme->GetHistogram()->SetLineColor(0);
        combined_tme->GetHistogram()->SetMarkerColor(0);
        combined_tme->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GeV/c]");
        combined_tme->GetHistogram()->GetYaxis()->SetTitle("Computation Time per Event [ms]");
        combined_tme->GetHistogram()->SetTitle(trpts[ip].c_str());
        combined_tme->GetHistogram()->GetXaxis()->SetRangeUser(trpts_rng[ip][0], trpts_rng[ip][1]);
        combined_tme->Draw("nostack p");
        Legend leg_tme(trpts_name[ip].c_str(), PadWindow(0.60, 0.80, 0.60, 0.80));
        leg_tme()->AddEntry((*hKFtme)(), "Kalman Fitter", "lp");
        leg_tme()->AddEntry((*hCKtme)(), "Choutko", "lp");
        leg_tme()->AddEntry((*hHCtme)(), "New", "lp");
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
        combined_eft->GetHistogram()->GetYaxis()->SetTitle("Fitting Efficiency (Normalized by Choutko)");
        combined_eft->GetHistogram()->SetTitle(trpts[ip].c_str());
        combined_eft->GetHistogram()->GetXaxis()->SetRangeUser(trpts_rng[ip][0], trpts_rng[ip][1]);
        combined_eft->Draw("nostack p");
        Legend leg_eft(trpts_name[ip].c_str(), PadWindow(0.60, 0.80, 0.20, 0.40));
        leg_eft()->AddEntry((*hKFeft)(), "Kalman Fitter", "lp");
        leg_eft()->AddEntry((*hCKeft)(), "Choutko", "lp");
        leg_eft()->AddEntry((*hHCeft)(), "New", "lp");
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
        combined_men->GetHistogram()->SetTitle(trpts[ip].c_str());
        combined_men->GetHistogram()->GetXaxis()->SetRangeUser(trpts_rng[ip][0], trpts_rng[ip][1]);
        combined_men->Draw("nostack p");
        Legend leg_men(trpts_name[ip].c_str(), PadWindow(0.60, 0.80, 0.60, 0.80));
        leg_men()->AddEntry((*hKFmen)(), "Kalman Fitter", "lp");
        leg_men()->AddEntry((*hCKmen)(), "Choutko", "lp");
        leg_men()->AddEntry((*hHCmen)(), "New", "lp");
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
        combined_rso->GetHistogram()->SetTitle(trpts[ip].c_str());
        combined_rso->GetHistogram()->GetXaxis()->SetRangeUser(trpts_rng[ip][0], trpts_rng[ip][1]);
        combined_rso->SetMaximum(150.0);
        combined_rso->Draw("nostack p");
        Legend leg_rso(trpts_name[ip].c_str(), PadWindow(0.20, 0.40, 0.60, 0.80));
        leg_rso()->AddEntry((*hKFrso)(), "Kalman Fitter", "lp");
        leg_rso()->AddEntry((*hCKrso)(), "Choutko", "lp");
        leg_rso()->AddEntry((*hHCrso)(), "New", "lp");
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
        combined_relrso->GetHistogram()->SetTitle(trpts[ip].c_str());
        combined_relrso->GetHistogram()->GetXaxis()->SetRangeUser(trpts_rng[ip][0], trpts_rng[ip][1]);
        combined_relrso->Draw("nostack p");
        Legend leg_relrso(trpts_name[ip].c_str(), PadWindow(0.60, 0.80, 0.20, 0.40));
        leg_relrso()->AddEntry((*hKFrelrso)(), "Kalman Fitter", "lp");
        leg_relrso()->AddEntry((*hCKrelrso)(), "Choutko", "lp");
        leg_relrso()->AddEntry((*hHCrelrso)(), "New", "lp");
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        editor.save();
    }
    editor.close();
    
    TFile * ofle = new TFile("out/doc/fit_geom.root", "RECREATE");
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
