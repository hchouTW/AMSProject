#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "HistFit1D.h"
#include "HistFit1D.C"

void stdfmt(TH1* hist) {
    if (hist == nullptr) return;
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleFont(43);
    hist->GetXaxis()->SetTitleSize(20);
    hist->GetXaxis()->SetTitleOffset(2.5);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelSize(15);
    hist->GetYaxis()->SetTitleFont(43);
    hist->GetYaxis()->SetTitleSize(20);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelSize(15);
}

double counting_error_func(double rig) { // in GeV
    double eff_sel = 0.003;
    double err_cfr = 0.05 * std::pow((1.0 / rig), 2.7);
    double error = std::hypot(eff_sel, err_cfr);
    return error;
}

double acceptance_error_func(double rig) { // in GeV
    double err_accp = 2.0e-02 + 2.66521e-02 * std::pow(rig, -1.23130e-01) * std::exp(-2.91143e-02 * rig);
    if (!std::isfinite(err_accp)) err_accp = 2.0e-02;
    return err_accp;
}

double rigidity_scale_func(double low, double up) { // in GeV
    const double err_magmap = 0.00025; // magnetic field map for AMS inner tracker
    const double err_tmpcrr = 0.00100; // temperature correction to the magnetic field map
    const double shift = 1.0 / 26.0; // in 1/TeV
    const double powlaw = 2.7;
    double intu = 1.0 / (low * 0.001); // in 1/TeV
    double intl = 1.0 / ( up * 0.001); // in 1/TeV
    double Fpos = std::pow(intu + shift, powlaw - 1.0) - std::pow(intl + shift, powlaw - 1.0);
    double Fneg = std::pow(intu - shift, powlaw - 1.0) - std::pow(intl - shift, powlaw - 1.0);
    double Fcen = std::pow(intu        , powlaw - 1.0) - std::pow(intl        , powlaw - 1.0);
    double err_misalg = 0.5 * (Fpos - Fneg) / Fcen; // tracker misalignment
    double rscl = std::hypot(err_misalg, std::hypot(err_magmap, err_tmpcrr));
    if (!std::isfinite(rscl)) rscl = 0.0;
    return rscl;
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "45";
    
    TGraphErrors* offapp = (TGraphErrors*) (TFile::Open("others/20160406MIT.root")->Get("gpbarp"));

    TFile* file_acc = TFile::Open("out/apflux_acc.root");

    TFile* file_cc = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr_l1o9flux%s/YiMdst.root", subv.c_str()));
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/iss%s", subv.c_str()));

    const Axis& AXrig = Hist::Head("hLtfP_sqrm")->xaxis();
    
    Hist* hNchi = Hist::New("hNchi", HistAxis(AXrig, "#chi^{2}"));
    
    Hist* hNcntPr = Hist::New("hNcntPr", HistAxis(AXrig, "Number"));
    Hist* hNcntAp = Hist::New("hNcntAp", HistAxis(AXrig, "Number"));
                                    
    Hist* hErrNcnt = Hist::New("hErrNcnt", HistAxis(AXrig, "RelError of Count (%)"));
    Hist* hErrTmpl = Hist::New("hErrTmpl", HistAxis(AXrig, "RelError of Template (%)"));
    Hist* hErrRscl = Hist::New("hErrRscl", HistAxis(AXrig, "RelError of Rigidity Scale (%)"));
    Hist* hErrAccp = Hist::New("hErrAccp", HistAxis(AXrig, "RelError of Acceptance (%)"));
    Hist* hErrCC   = Hist::New("hErrCC",   HistAxis(AXrig, "RelError of CC (%)"));
    
    Hist* hErrStat = Hist::New("hErrStat", HistAxis(AXrig, "RelError of Stat (%)"));
    Hist* hErrSyst = Hist::New("hErrSyst", HistAxis(AXrig, "RelError of Syst (%)"));
    Hist* hErrTotl = Hist::New("hErrTotl", HistAxis(AXrig, "RelError of Totl (%)"));
    
    Hist* hAccCorr = Hist::New("hAccCorr", HistAxis(AXrig, "Acceptance Correction"));
    Hist* hAppStat = Hist::New("hAppStat", HistAxis(AXrig, "#bar{p}/p"));
    Hist* hAppTotl = Hist::New("hAppTotl", HistAxis(AXrig, "#bar{p}/p"));
   
    PdfEditor editor(Window(), "apflux_app", "out");

    bool sw_fluc = false; 
    bool sw_build_hist = true; 

    /////////////////
    ////         ////
    ////   LTF   ////
    ////         ////
    /////////////////
    TH1D* hLTF_acc = (TH1D*) file_acc->Get("hLtf_accp");
    std::vector<Hist*> vhLTF_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hLtfP_sqrm"));
    std::vector<Hist*> vhLTF_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hLtfN_sqrm"));
    std::vector<Hist*> vhLTF_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hLtfP_sqrm_pr"));
    std::vector<Hist*> vhLTF_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hLtfN_sqrm_el"));
    
    std::array<int, 2> RANGE_LTF({ 1, 11 });
    HistFit::Axis1D AX1D_LTF(
        "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
        "Events/Bin",
        (*Hist::Head("hLtfP_sqrm"))()->GetYaxis());

    for (int ir = RANGE_LTF[0]; ir <= RANGE_LTF[1]; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        TH1D*              hlink_smp  = (TH1D*)((*vhLTF_neg.at(ir))());
        std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhLTF_sig.at(ir))()), (TH1D*)((*vhLTF_bkg.at(ir))()) });

        std::string prefix = Form("hLTF_R%03d_", ir);
        HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhLTF_pos.at(ir))()), AX1D_LTF);
        HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_LTF, prefix, sw_fluc, sw_build_hist);
        if (!fit1D.status()) continue;
        double nchi = fit1D.nchi();

        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dsig = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dbkg = fit1D.wgt_tmps(1);

        double num_ref = h1Dref.data().sum();
        double num_smp = fit1D.nsmp();
        double num_sig = fit1D.wgts(0);
        double num_bkg = fit1D.wgts(1);

        double err_ncnt = counting_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_tmpl = fit1D.fluc(0).err / num_sig;
        double err_rscl = 2.0 * rigidity_scale_func(AXrig()(ir-1), AXrig()(ir));
        double err_accp = acceptance_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_cc   = 0.0;

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_tmpl*err_tmpl + err_rscl*err_rscl + err_accp*err_accp + err_cc*err_cc);
        double err_totl = std::hypot(err_stat, err_syst);

        double acc_corr     = hLTF_acc->Interpolate(rig);
        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
        
        std::cerr << Form("LTF STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
            fit1D.status(), ir, AXrig()(ir-1), AXrig()(ir),
            num_smp, num_sig, num_bkg,
            app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
            nchi);

        (*hNchi)()->SetBinContent(ir, nchi); 
        
        (*hNcntPr)()->SetBinContent(ir, num_ref); 
        (*hNcntAp)()->SetBinContent(ir, num_sig); 
        
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrRscl)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        (*hErrCC  )()->SetBinContent(ir, 100.0 * err_cc  );
        
        (*hErrStat)()->SetBinContent(ir, 100.0 * err_stat); 
        (*hErrSyst)()->SetBinContent(ir, 100.0 * err_syst);
        (*hErrTotl)()->SetBinContent(ir, 100.0 * err_totl);
        
        (*hAccCorr)()->SetBinContent(ir, acc_corr);

        (*hAppStat)()->SetBinContent(ir, app_val); 
        (*hAppStat)()->SetBinError  (ir, app_err_stat); 
        
        (*hAppTotl)()->SetBinContent(ir, app_val);
        (*hAppTotl)()->SetBinError  (ir, app_err_totl);

        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hsig = Hist::New(h1Dsig.get());
        Hist* hbkg = Hist::New(h1Dbkg.get());
        
        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
		hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
   
        THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
        hfit->Draw("nostack hist");

        Hist* hcanvas = Hist::New(
            Form("%scanvas", prefix.c_str()), 
            HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.3 * hfit->GetHistogram()->GetMaximum())));
        (*hcanvas)()->GetXaxis()->SetTitle(AX1D_LTF.name_x().c_str());
        (*hcanvas)()->GetYaxis()->SetTitle(AX1D_LTF.name_y().c_str());
        (*hcanvas)()->Draw();
        
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.1f - %.1f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
        Hist* hfluc = Hist::New(h1Dfluc.get());
        
        editor.create();
		hfluc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfluc->draw("hist");
        editor.save();
    }


    /////////////////
    ////         ////
    ////   LRH   ////
    ////         ////
    /////////////////
    TH1D* hLRH_acc = (TH1D*) file_acc->Get("hLrh_accp");
    std::vector<Hist*> vhLRH_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hLrhP_sqrm"));
    std::vector<Hist*> vhLRH_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hLrhN_sqrm"));
    std::vector<Hist*> vhLRH_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hLrhP_sqrm_pr"));
    std::vector<Hist*> vhLRH_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hLrhN_sqrm_el"));
    
    std::array<int, 2> RANGE_LRH({ 12, 25 });
    HistFit::Axis1D AX1D_LRH(
        "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
        "Events/Bin",
        (*Hist::Head("hLrhP_sqrm"))()->GetYaxis());

    for (int ir = RANGE_LRH[0]; ir <= RANGE_LRH[1]; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        TH1D*              hlink_smp  = (TH1D*)((*vhLRH_neg.at(ir))());
        std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhLRH_sig.at(ir))()), (TH1D*)((*vhLRH_bkg.at(ir))()) });

        std::string prefix = Form("hLRH_R%03d_", ir);
        HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhLRH_pos.at(ir))()), AX1D_LRH);
        HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_LRH, prefix, sw_fluc, sw_build_hist);
        if (!fit1D.status()) continue;
        double nchi = fit1D.nchi();

        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dsig = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dbkg = fit1D.wgt_tmps(1);

        double num_ref = h1Dref.data().sum();
        double num_smp = fit1D.nsmp();
        double num_sig = fit1D.wgts(0);
        double num_bkg = fit1D.wgts(1);

        double err_ncnt = counting_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_tmpl = fit1D.fluc(0).err / num_sig;
        double err_rscl = 2.0 * rigidity_scale_func(AXrig()(ir-1), AXrig()(ir));
        double err_accp = acceptance_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_cc   = 0.0;

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_tmpl*err_tmpl + err_rscl*err_rscl + err_accp*err_accp + err_cc*err_cc);
        double err_totl = std::hypot(err_stat, err_syst);

        double acc_corr     = hLRH_acc->Interpolate(rig);
        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
        
        std::cerr << Form("LRH STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
            fit1D.status(), ir, AXrig()(ir-1), AXrig()(ir),
            num_smp, num_sig, num_bkg,
            app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
            nchi);

        (*hNchi)()->SetBinContent(ir, nchi); 
        
        (*hNcntPr)()->SetBinContent(ir, num_ref); 
        (*hNcntAp)()->SetBinContent(ir, num_sig); 
        
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrRscl)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        (*hErrCC  )()->SetBinContent(ir, 100.0 * err_cc  );
        
        (*hErrStat)()->SetBinContent(ir, 100.0 * err_stat); 
        (*hErrSyst)()->SetBinContent(ir, 100.0 * err_syst);
        (*hErrTotl)()->SetBinContent(ir, 100.0 * err_totl);
        
        (*hAccCorr)()->SetBinContent(ir, acc_corr);

        (*hAppStat)()->SetBinContent(ir, app_val); 
        (*hAppStat)()->SetBinError  (ir, app_err_stat); 
        
        (*hAppTotl)()->SetBinContent(ir, app_val);
        (*hAppTotl)()->SetBinError  (ir, app_err_totl);

        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hsig = Hist::New(h1Dsig.get());
        Hist* hbkg = Hist::New(h1Dbkg.get());
        
        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
		hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
   
        THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
        hfit->Draw("nostack hist");

        Hist* hcanvas = Hist::New(
            Form("%scanvas", prefix.c_str()), 
            HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.3 * hfit->GetHistogram()->GetMaximum())));
        (*hcanvas)()->GetXaxis()->SetTitle(AX1D_LRH.name_x().c_str());
        (*hcanvas)()->GetYaxis()->SetTitle(AX1D_LRH.name_y().c_str());
        (*hcanvas)()->Draw();
        
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.1f - %.1f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
        Hist* hfluc = Hist::New(h1Dfluc.get());
        
        editor.create();
		hfluc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfluc->draw("hist");
        editor.save();
    }



    /////////////////
    ////         ////
    ////   IIN   ////
    ////         ////
    /////////////////
    TH1D* hIIN_acc = (TH1D*) file_acc->Get("hIin_accp");
    std::vector<Hist*> vhIIN_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hIinP_llr"));
    std::vector<Hist*> vhIIN_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hIinN_llr"));
    std::vector<Hist*> vhIIN_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hIinP_llr_pr"));
    std::vector<Hist*> vhIIN_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hIinN_llr_el"));
    
    std::array<int, 2> RANGE_IIN({ 26, 36 });
    HistFit::Axis1D AX1D_IIN(
        "TRD Estimator",
        "Events/Bin",
        (*Hist::Head("hIinP_llr"))()->GetYaxis());

    for (int ir = RANGE_IIN[0]; ir <= RANGE_IIN[1]; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        TH1D*              hlink_smp  = (TH1D*)((*vhIIN_neg.at(ir))());
        std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhIIN_sig.at(ir))()), (TH1D*)((*vhIIN_bkg.at(ir))()) });

        std::string prefix = Form("hIIN_R%03d_", ir);
        HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhIIN_pos.at(ir))()), AX1D_IIN);
        HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_IIN, prefix, sw_fluc, sw_build_hist);
        if (!fit1D.status()) continue;
        double nchi = fit1D.nchi();

        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dsig = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dbkg = fit1D.wgt_tmps(1);

        double num_ref = h1Dref.data().sum();
        double num_smp = fit1D.nsmp();
        double num_sig = fit1D.wgts(0);
        double num_bkg = fit1D.wgts(1);

        double err_ncnt = counting_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_tmpl = fit1D.fluc(0).err / num_sig;
        double err_rscl = 2.0 * rigidity_scale_func(AXrig()(ir-1), AXrig()(ir));
        double err_accp = acceptance_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_cc   = 0.0;

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_tmpl*err_tmpl + err_rscl*err_rscl + err_accp*err_accp + err_cc*err_cc);
        double err_totl = std::hypot(err_stat, err_syst);

        double acc_corr     = hIIN_acc->Interpolate(rig);
        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
        
        std::cerr << Form("IIN STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
            fit1D.status(), ir, AXrig()(ir-1), AXrig()(ir),
            num_smp, num_sig, num_bkg,
            app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
            nchi);

        (*hNchi)()->SetBinContent(ir, nchi); 
        
        (*hNcntPr)()->SetBinContent(ir, num_ref); 
        (*hNcntAp)()->SetBinContent(ir, num_sig); 
        
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrRscl)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        (*hErrCC  )()->SetBinContent(ir, 100.0 * err_cc  );
        
        (*hErrStat)()->SetBinContent(ir, 100.0 * err_stat); 
        (*hErrSyst)()->SetBinContent(ir, 100.0 * err_syst);
        (*hErrTotl)()->SetBinContent(ir, 100.0 * err_totl);
        
        (*hAccCorr)()->SetBinContent(ir, acc_corr);

        (*hAppStat)()->SetBinContent(ir, app_val); 
        (*hAppStat)()->SetBinError  (ir, app_err_stat); 
        
        (*hAppTotl)()->SetBinContent(ir, app_val);
        (*hAppTotl)()->SetBinError  (ir, app_err_totl);

        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hsig = Hist::New(h1Dsig.get());
        Hist* hbkg = Hist::New(h1Dbkg.get());
        
        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
		hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
   
        THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
        hfit->Draw("nostack hist");

        Hist* hcanvas = Hist::New(
            Form("%scanvas", prefix.c_str()), 
            HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.3 * hfit->GetHistogram()->GetMaximum())));
        (*hcanvas)()->GetXaxis()->SetTitle(AX1D_IIN.name_x().c_str());
        (*hcanvas)()->GetYaxis()->SetTitle(AX1D_IIN.name_y().c_str());
        (*hcanvas)()->Draw();
        
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.1f - %.1f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("e^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
        Hist* hfluc = Hist::New(h1Dfluc.get());
        
        editor.create();
		hfluc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfluc->draw("hist");
        editor.save();
    }




    /////////////////
    ////         ////
    ////   IEX   ////
    ////         ////
    /////////////////
    TH1D* hIEX_acc = (TH1D*) file_acc->Get("hIex_accp");
    std::vector<Hist*> vhIEX_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hIexP_llr"));
    std::vector<Hist*> vhIEX_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hIexN_llr"));
    std::vector<Hist*> vhIEX_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hIexP_llr_pr"));
    std::vector<Hist*> vhIEX_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hIexN_llr_el"));
    
    std::array<int, 2> RANGE_IEX({ 37, 45 });
    HistFit::Axis1D AX1D_IEX(
        "TRD Estimator",
        "Events/Bin",
        (*Hist::Head("hIexP_llr"))()->GetYaxis());

    for (int ir = RANGE_IEX[0]; ir <= RANGE_IEX[1]; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        TH1D*              hlink_smp  = (TH1D*)((*vhIEX_neg.at(ir))());
        std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhIEX_sig.at(ir))()), (TH1D*)((*vhIEX_bkg.at(ir))()) });

        std::string prefix = Form("hIEX_R%03d_", ir);
        HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhIEX_pos.at(ir))()), AX1D_IEX);
        HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_IEX, prefix, sw_fluc, sw_build_hist);
        if (!fit1D.status()) continue;
        double nchi = fit1D.nchi();

        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dsig = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dbkg = fit1D.wgt_tmps(1);

        double num_ref = h1Dref.data().sum();
        double num_smp = fit1D.nsmp();
        double num_sig = fit1D.wgts(0);
        double num_bkg = fit1D.wgts(1);

        double err_ncnt = counting_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_tmpl = fit1D.fluc(0).err / num_sig;
        double err_rscl = 2.0 * rigidity_scale_func(AXrig()(ir-1), AXrig()(ir));
        double err_accp = acceptance_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_cc   = 0.0;

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_tmpl*err_tmpl + err_rscl*err_rscl + err_accp*err_accp + err_cc*err_cc);
        double err_totl = std::hypot(err_stat, err_syst);

        double acc_corr     = hIEX_acc->Interpolate(rig);
        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
        
        std::cerr << Form("IEX STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
            fit1D.status(), ir, AXrig()(ir-1), AXrig()(ir),
            num_smp, num_sig, num_bkg,
            app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
            nchi);

        (*hNchi)()->SetBinContent(ir, nchi); 
        
        (*hNcntPr)()->SetBinContent(ir, num_ref); 
        (*hNcntAp)()->SetBinContent(ir, num_sig); 
        
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrRscl)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        (*hErrCC  )()->SetBinContent(ir, 100.0 * err_cc  );
        
        (*hErrStat)()->SetBinContent(ir, 100.0 * err_stat); 
        (*hErrSyst)()->SetBinContent(ir, 100.0 * err_syst);
        (*hErrTotl)()->SetBinContent(ir, 100.0 * err_totl);
        
        (*hAccCorr)()->SetBinContent(ir, acc_corr);

        (*hAppStat)()->SetBinContent(ir, app_val); 
        (*hAppStat)()->SetBinError  (ir, app_err_stat); 
        
        (*hAppTotl)()->SetBinContent(ir, app_val);
        (*hAppTotl)()->SetBinError  (ir, app_err_totl);

        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hsig = Hist::New(h1Dsig.get());
        Hist* hbkg = Hist::New(h1Dbkg.get());
        
        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
		hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
   
        THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
        hfit->Draw("nostack hist");

        Hist* hcanvas = Hist::New(
            Form("%scanvas", prefix.c_str()), 
            HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.3 * hfit->GetHistogram()->GetMaximum())));
        (*hcanvas)()->GetXaxis()->SetTitle(AX1D_IEX.name_x().c_str());
        (*hcanvas)()->GetYaxis()->SetTitle(AX1D_IEX.name_y().c_str());
        (*hcanvas)()->Draw();
        
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.1f - %.1f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("e^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
        Hist* hfluc = Hist::New(h1Dfluc.get());
        
        editor.create();
		hfluc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfluc->draw("hist");
        editor.save();
    }



    /////////////////
    ////         ////
    ////   HEX   ////
    ////         ////
    /////////////////
    TH1D* hHEX_acc = (TH1D*) file_acc->Get("hHex_accp");
    Hist* hHEX_pr = Hist::New("hHex_cc_pr", (TH2D*)((*Hist::Head("hHexP_cc"))()));
    Hist* hHEX_cc = Hist::New("hHex_cc_cc", (TH2D*)(file_cc->Get("hHexN_cc_MC_FLUX27")));
    std::vector<Hist*> vhHEX_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hHexP_cc"));
    std::vector<Hist*> vhHEX_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hHexN_cc"));
    std::vector<Hist*> vhHEX_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hHex_cc_pr"));
    std::vector<Hist*> vhHEX_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hHex_cc_cc"));
    
    std::array<int, 2> RANGE_HEX({ 46, 54 });
    HistFit::Axis1D AX1D_HEX(
        "CC Estimator",
        "Events/Bin",
        (*Hist::Head("hHexP_cc"))()->GetYaxis(), 5, 50);

    for (int ir = RANGE_HEX[0]; ir <= RANGE_HEX[1]; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        TH1D*              hlink_smp  = (TH1D*)((*vhHEX_neg.at(ir))());
        std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhHEX_sig.at(ir))()), (TH1D*)((*vhHEX_bkg.at(ir))()) });

        std::string prefix = Form("hHEX_R%03d_", ir);
        HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhHEX_pos.at(ir))()), AX1D_HEX);
        HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_HEX, prefix, sw_fluc, sw_build_hist);
        if (!fit1D.status()) continue;
        double nchi = fit1D.nchi();
        
        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dsig = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dbkg = fit1D.wgt_tmps(1);

        double num_ref = h1Dref.data().sum();
        double num_smp = fit1D.nsmp();
        double num_sig = fit1D.wgts(0);
        double num_bkg = fit1D.wgts(1);

        double err_ncnt = counting_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_tmpl = fit1D.fluc(0).err / num_sig;
        double err_rscl = 2.0 * rigidity_scale_func(AXrig()(ir-1), AXrig()(ir));
        double err_accp = acceptance_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_cc   = 0.0;

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_tmpl*err_tmpl + err_rscl*err_rscl + err_accp*err_accp + err_cc*err_cc);
        double err_totl = std::hypot(err_stat, err_syst);

        double acc_corr     = hHEX_acc->Interpolate(rig);
        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
        
        std::cerr << Form("HEX STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
            fit1D.status(), ir, AXrig()(ir-1), AXrig()(ir),
            num_smp, num_sig, num_bkg,
            app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
            nchi);

        (*hNchi)()->SetBinContent(ir, nchi); 
        
        (*hNcntPr)()->SetBinContent(ir, num_ref); 
        (*hNcntAp)()->SetBinContent(ir, num_sig); 
        
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrRscl)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        (*hErrCC  )()->SetBinContent(ir, 100.0 * err_cc  );
        
        (*hErrStat)()->SetBinContent(ir, 100.0 * err_stat); 
        (*hErrSyst)()->SetBinContent(ir, 100.0 * err_syst);
        (*hErrTotl)()->SetBinContent(ir, 100.0 * err_totl);
        
        (*hAccCorr)()->SetBinContent(ir, acc_corr);

        (*hAppStat)()->SetBinContent(ir, app_val); 
        (*hAppStat)()->SetBinError  (ir, app_err_stat); 
        
        (*hAppTotl)()->SetBinContent(ir, app_val);
        (*hAppTotl)()->SetBinError  (ir, app_err_totl);

        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hsig = Hist::New(h1Dsig.get());
        Hist* hbkg = Hist::New(h1Dbkg.get());
        
        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
		hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
   
        THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
        hfit->Draw("nostack hist");

        Hist* hcanvas = Hist::New(
            Form("%scanvas", prefix.c_str()), 
            HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.3 * hfit->GetHistogram()->GetMaximum())));
        (*hcanvas)()->GetXaxis()->SetTitle(AX1D_HEX.name_x().c_str());
        (*hcanvas)()->GetYaxis()->SetTitle(AX1D_HEX.name_y().c_str());
        (*hcanvas)()->Draw();
        
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.1f - %.1f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("p^{+} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
        Hist* hfluc = Hist::New(h1Dfluc.get());
        
        editor.create();
		hfluc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfluc->draw("hist");
        editor.save();
    }



    /////////////////
    ////         ////
    ////   HFS   ////
    ////         ////
    /////////////////
    TH1D* hHFS_acc = (TH1D*) file_acc->Get("hHfs_accp");
    Hist* hHFS_pr = Hist::New("hHfs_cc_pr", (TH2D*)((*Hist::Head("hHfsP_cc"))()));
    Hist* hHFS_cc = Hist::New("hHfs_cc_cc", (TH2D*)(file_cc->Get("hHfsN_cc_MC_FLUX27")));
    std::vector<Hist*> vhHFS_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hHfsP_cc"));
    std::vector<Hist*> vhHFS_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hHfsN_cc"));
    std::vector<Hist*> vhHFS_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hHfs_cc_pr"));
    std::vector<Hist*> vhHFS_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hHfs_cc_cc"));
    
    std::array<int, 2> RANGE_HFS({ 55, 58 });
    HistFit::Axis1D AX1D_HFS(
        "CC Estimator",
        "Events/Bin",
        (*Hist::Head("hHfsP_cc"))()->GetYaxis(), 5, 50);
   
    for (int ir = RANGE_HFS[0]; ir <= RANGE_HFS[1]; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        TH1D*              hlink_smp  = (TH1D*)((*vhHFS_neg.at(ir))());
        std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhHFS_sig.at(ir))()), (TH1D*)((*vhHFS_bkg.at(ir))()) });

        std::string prefix = Form("hHFS_R%03d_", ir);
        HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhHFS_pos.at(ir))()), AX1D_HFS);
        HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_HFS, prefix, sw_fluc, sw_build_hist);
        if (!fit1D.status()) continue;
        double nchi = fit1D.nchi();

        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dsig = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dbkg = fit1D.wgt_tmps(1);

        double num_ref = h1Dref.data().sum();
        double num_smp = fit1D.nsmp();
        double num_sig = fit1D.wgts(0);
        double num_bkg = fit1D.wgts(1);

        double err_ncnt = counting_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_tmpl = fit1D.fluc(0).err / num_sig;
        double err_rscl = 2.0 * rigidity_scale_func(AXrig()(ir-1), AXrig()(ir));
        double err_accp = acceptance_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_cc   = 0.0;

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_tmpl*err_tmpl + err_rscl*err_rscl + err_accp*err_accp + err_cc*err_cc);
        double err_totl = std::hypot(err_stat, err_syst);

        double acc_corr     = hHFS_acc->Interpolate(rig);
        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
        
        std::cerr << Form("HFS STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
            fit1D.status(), ir, AXrig()(ir-1), AXrig()(ir),
            num_smp, num_sig, num_bkg,
            app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
            nchi);

        (*hNchi)()->SetBinContent(ir, nchi); 
        
        (*hNcntPr)()->SetBinContent(ir, num_ref); 
        (*hNcntAp)()->SetBinContent(ir, num_sig); 
        
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrRscl)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        (*hErrCC  )()->SetBinContent(ir, 100.0 * err_cc  );
        
        (*hErrStat)()->SetBinContent(ir, 100.0 * err_stat); 
        (*hErrSyst)()->SetBinContent(ir, 100.0 * err_syst);
        (*hErrTotl)()->SetBinContent(ir, 100.0 * err_totl);
        
        (*hAccCorr)()->SetBinContent(ir, acc_corr);

        (*hAppStat)()->SetBinContent(ir, app_val); 
        (*hAppStat)()->SetBinError  (ir, app_err_stat); 
        
        (*hAppTotl)()->SetBinContent(ir, app_val);
        (*hAppTotl)()->SetBinError  (ir, app_err_totl);

        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hsig = Hist::New(h1Dsig.get());
        Hist* hbkg = Hist::New(h1Dbkg.get());
        
        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
		hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
   
        THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
        hfit->Draw("nostack hist");

        Hist* hcanvas = Hist::New(
            Form("%scanvas", prefix.c_str()), 
            HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.3 * hfit->GetHistogram()->GetMaximum())));
        (*hcanvas)()->GetXaxis()->SetTitle(AX1D_HFS.name_x().c_str());
        (*hcanvas)()->GetYaxis()->SetTitle(AX1D_HFS.name_y().c_str());
        (*hcanvas)()->Draw();
        
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.1f - %.1f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("p^{+} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
        Hist* hfluc = Hist::New(h1Dfluc.get());
        
        editor.create();
		hfluc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfluc->draw("hist");
        editor.save();
    }

    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hNchi)()->Draw("hist");
    editor.save();
    

    editor.create();
    editor.cd(0, PadAxis(1, 1));
    (*hNcntPr)()->Draw("hist");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    (*hNcntAp)()->Draw("hist");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hErrStat)()->Draw("hist");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hErrSyst)()->Draw("hist");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hAppStat)()->Draw("pe");
    offapp->Draw("pe same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    (*hAppStat)()->Draw("pe");
    offapp->Draw("pe same");
    editor.save();

    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hAppTotl)()->Draw("pe");
    offapp->Draw("pe same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    (*hAppTotl)()->Draw("pe");
    offapp->Draw("pe same");
    editor.save();


    editor.close();

    TFile * ofle = new TFile("out/apflux_app.root", "RECREATE");
    ofle->cd();

    (*hNchi)()->Write();
    
    (*hNcntPr)()->Write();
    (*hNcntAp)()->Write();
    
    (*hErrNcnt)()->Write();
    (*hErrTmpl)()->Write();
    (*hErrRscl)()->Write();
    (*hErrAccp)()->Write();
    (*hErrCC  )()->Write();
    
    (*hErrStat)()->Write();
    (*hErrSyst)()->Write();
    (*hErrTotl)()->Write();
    
    (*hAccCorr)()->Write();
    (*hAppStat)()->Write();
    (*hAppTotl)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
