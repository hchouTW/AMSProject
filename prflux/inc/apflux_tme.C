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
    std::string subt = "1M";

    TFile* file_acc = TFile::Open("out/apflux_acc.root");

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/iss%s", subv.c_str()));

    const Axis& AXtme = Hist::Head(Form("hT%sLtfP_sqrm", subt.c_str()))->xaxis();
    const Axis& AXrig = Hist::Head(Form("hT%sLtfP_sqrm", subt.c_str()))->yaxis();
    
    Hist* hTNchi = Hist::New("hTNchi", HistAxis(AXtme, AXrig, "#chi^{2}"));
    
    Hist* hTNcntPr = Hist::New("hTNcntPr", HistAxis(AXtme, AXrig, "Number"));
    Hist* hTNcntAp = Hist::New("hTNcntAp", HistAxis(AXtme, AXrig, "Number"));
    
    Hist* hTErrNcnt = Hist::New("hTErrNcnt", HistAxis(AXtme, AXrig, "RelError of Count (%)"));
    Hist* hTErrTmpl = Hist::New("hTErrTmpl", HistAxis(AXtme, AXrig, "RelError of Template (%)"));
    Hist* hTErrRscl = Hist::New("hTErrRscl", HistAxis(AXtme, AXrig, "RelError of Rigidity Scale (%)"));
    Hist* hTErrAccp = Hist::New("hTErrAccp", HistAxis(AXtme, AXrig, "RelError of Acceptance (%)"));
    Hist* hTErrCC   = Hist::New("hTErrCC",   HistAxis(AXtme, AXrig, "RelError of CC (%)"));
    
    Hist* hTErrStat = Hist::New("hTErrStat", HistAxis(AXtme, AXrig, "RelError of Stat (%)"));
    Hist* hTErrSyst = Hist::New("hTErrSyst", HistAxis(AXtme, AXrig, "RelError of Syst (%)"));
    Hist* hTErrTotl = Hist::New("hTErrTotl", HistAxis(AXtme, AXrig, "RelError of Totl (%)"));
    
    Hist* hTAccCorr = Hist::New("hTAccCorr", HistAxis(AXtme, AXrig, "Acceptance Correction"));
    Hist* hTAppStat = Hist::New("hTAppStat", HistAxis(AXtme, AXrig, "#bar{p}/p"));
    Hist* hTAppTotl = Hist::New("hTAppTotl", HistAxis(AXtme, AXrig, "#bar{p}/p"));
   
    std::vector<int> exclude_times;
    if (subt == "1M") exclude_times = std::vector<int>({ 46, 47, 98, 100 });

    /////////////////
    ////         ////
    ////   LTF   ////
    ////         ////
    /////////////////
    //PdfEditor editor_ltf(Window(), "apflux_tme_ltf", "out");
    
    TH1D* hLTF_acc = (TH1D*) file_acc->Get("hLtf_accp");
    std::vector<Hist*> vhTLTF_pos = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sLtfP_sqrm", subt.c_str())));
    std::vector<Hist*> vhTLTF_neg = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sLtfN_sqrm", subt.c_str())));
    std::vector<Hist*> vhTLTF_sig = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sLtfP_sqrm_pr", subt.c_str())));
    std::vector<Hist*> vhTLTF_bkg = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sLtfN_sqrm_el", subt.c_str())));
    
    std::vector<Hist*> vhRLTF_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hRLtfP_sqrm_pr"));
    std::vector<Hist*> vhRLTF_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hRLtfN_sqrm_el"));
    
    //std::array<int, 2> RANGE_LTF({ 1, 5 });
    std::array<int, 2> RANGE_LTF({ 1, 5 });
    HistFit::Axis1D AX1D_LTF(
        "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
        "Events/Bin",
        (*Hist::Head(Form("hT%sLtfP_sqrm", subt.c_str())))()->GetZaxis());

    for (int it = 1; it <= AXtme.nbin(); ++it) {
        bool skip = false;
        for (auto&& itime : exclude_times) { if (it == itime) { skip = true; break; } }
        if (skip) continue;

        double tme = AXtme.center(it, AxisScale::kLinear);
        std::vector<Hist*> vhRLTF_pos = Hist::ProjectAll(HistProj::kY, vhTLTF_pos.at(it));
        std::vector<Hist*> vhRLTF_neg = Hist::ProjectAll(HistProj::kY, vhTLTF_neg.at(it));
        //std::vector<Hist*> vhRLTF_sig = Hist::ProjectAll(HistProj::kY, vhTLTF_sig.at(it));
        //std::vector<Hist*> vhRLTF_bkg = Hist::ProjectAll(HistProj::kY, vhTLTF_bkg.at(it));

        for (int ir = RANGE_LTF[0]; ir <= RANGE_LTF[1]; ++ir) {
            double rig = AXrig.center(ir, AxisScale::kLog);
            TH1D*              hlink_smp  = (TH1D*)((*vhRLTF_neg.at(ir))());
            std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhRLTF_sig.at(ir))()), (TH1D*)((*vhRLTF_bkg.at(ir))()) });

            std::string prefix = Form("hLTF_T%03dR%03d_", it, ir);
            HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhRLTF_pos.at(ir))()), AX1D_LTF);
            HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_LTF, prefix, false, true);
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
            
            std::cerr << Form("LTF STATUS %d T%03d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
                fit1D.status(), it, ir, AXrig()(ir-1), AXrig()(ir),
                num_smp, num_sig, num_bkg,
                app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
                nchi);

            (*hTNchi)()->SetBinContent(it, ir, nchi); 
            
            (*hTNcntPr)()->SetBinContent(it, ir, num_ref); 
            (*hTNcntAp)()->SetBinContent(it, ir, num_sig); 
            
            (*hTErrNcnt)()->SetBinContent(it, ir, 100.0 * err_ncnt); 
            (*hTErrTmpl)()->SetBinContent(it, ir, 100.0 * err_tmpl);
            (*hTErrRscl)()->SetBinContent(it, ir, 100.0 * err_rscl);
            (*hTErrAccp)()->SetBinContent(it, ir, 100.0 * err_accp);
            (*hTErrCC  )()->SetBinContent(it, ir, 100.0 * err_cc  );
            
            (*hTErrStat)()->SetBinContent(it, ir, 100.0 * err_stat); 
            (*hTErrSyst)()->SetBinContent(it, ir, 100.0 * err_syst);
            (*hTErrTotl)()->SetBinContent(it, ir, 100.0 * err_totl);
            
            (*hTAccCorr)()->SetBinContent(it, ir, acc_corr); 
            
            (*hTAppStat)()->SetBinContent(it, ir, app_val); 
            (*hTAppStat)()->SetBinError  (it, ir, app_err_stat); 
            
            (*hTAppTotl)()->SetBinContent(it, ir, app_val);
            (*hTAppTotl)()->SetBinError  (it, ir, app_err_totl);

            Hist* hsmp = Hist::New(h1Dsmp.get());
            Hist* hsum = Hist::New(h1Dsum.get());
            Hist* hsig = Hist::New(h1Dsig.get());
            Hist* hbkg = Hist::New(h1Dbkg.get());
            
            hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
	    	hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
	    	hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
	    	hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
            
            /*
            editor_ltf.create();
   
            THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
            hfit->Draw("nostack hist");

            Hist* hcanvas = Hist::New(
                Form("%scanvas", prefix.c_str()), 
                HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.5 * hfit->GetHistogram()->GetMaximum())));
            (*hcanvas)()->GetXaxis()->SetTitle(AX1D_LTF.name_x().c_str());
            (*hcanvas)()->GetYaxis()->SetTitle(AX1D_LTF.name_y().c_str());
            (*hcanvas)()->Draw();
            
            hfit->Draw("nostack hist same");
            (*hsmp)()->Draw("pe same");
            
            Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
            leg_table()->SetHeader(Form("Time %3d Rigidity %.1f - %.1f [GV/c]", it, AXrig()(ir-1), AXrig()(ir)));
            leg_table()->AddEntry((*hsmp)(), "Data", "lp");
            leg_table()->AddEntry((*hsum)(), "Sum", "l");
            leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
            leg_table()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
            leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
            leg_table()->SetFillColor(0);
            leg_table.draw();
            
            editor_ltf.save();
            */
            
            Hist::Delete(hsmp);
            Hist::Delete(hsum);
            Hist::Delete(hsig);
            Hist::Delete(hbkg);

            fit1D.delete_hist();
        }
        
        for (auto&& hist : vhRLTF_pos) Hist::Delete(hist);
        for (auto&& hist : vhRLTF_neg) Hist::Delete(hist);
        //for (auto&& hist : vhRLTF_sig) Hist::Delete(hist);
        //for (auto&& hist : vhRLTF_bkg) Hist::Delete(hist);
    }
    for (auto&& hist : vhTLTF_pos) Hist::Delete(hist);
    for (auto&& hist : vhTLTF_neg) Hist::Delete(hist);
    for (auto&& hist : vhTLTF_sig) Hist::Delete(hist);
    for (auto&& hist : vhTLTF_bkg) Hist::Delete(hist);
        
    for (auto&& hist : vhRLTF_sig) Hist::Delete(hist);
    for (auto&& hist : vhRLTF_bkg) Hist::Delete(hist);

    //editor_ltf.close();

    /////////////////
    ////         ////
    ////   LRH   ////
    ////         ////
    /////////////////
    //PdfEditor editor_lrh(Window(), "apflux_tme_lrh", "out");
    
    TH1D* hLRH_acc = (TH1D*) file_acc->Get("hLrh_accp");
    std::vector<Hist*> vhTLRH_pos = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sLrhP_sqrm", subt.c_str())));
    std::vector<Hist*> vhTLRH_neg = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sLrhN_sqrm", subt.c_str())));
    std::vector<Hist*> vhTLRH_sig = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sLrhP_sqrm_pr", subt.c_str())));
    std::vector<Hist*> vhTLRH_bkg = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sLrhN_sqrm_el", subt.c_str())));
    
    std::vector<Hist*> vhRLRH_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hRLrhP_sqrm_pr"));
    std::vector<Hist*> vhRLRH_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hRLrhN_sqrm_el"));
    
    //std::array<int, 2> RANGE_LRH({ 6, 12 });
    std::array<int, 2> RANGE_LRH({ 6, 12 });
    HistFit::Axis1D AX1D_LRH(
        "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
        "Events/Bin",
        (*Hist::Head(Form("hT%sLrhP_sqrm", subt.c_str())))()->GetZaxis());

    for (int it = 1; it <= AXtme.nbin(); ++it) {
        bool skip = false;
        for (auto&& itime : exclude_times) { if (it == itime) { skip = true; break; } }
        if (skip) continue;

        double tme = AXtme.center(it, AxisScale::kLinear);
        std::vector<Hist*> vhRLRH_pos = Hist::ProjectAll(HistProj::kY, vhTLRH_pos.at(it));
        std::vector<Hist*> vhRLRH_neg = Hist::ProjectAll(HistProj::kY, vhTLRH_neg.at(it));
        //std::vector<Hist*> vhRLRH_sig = Hist::ProjectAll(HistProj::kY, vhTLRH_sig.at(it));
        //std::vector<Hist*> vhRLRH_bkg = Hist::ProjectAll(HistProj::kY, vhTLRH_bkg.at(it));

        for (int ir = RANGE_LRH[0]; ir <= RANGE_LRH[1]; ++ir) {
            double rig = AXrig.center(ir, AxisScale::kLog);
            TH1D*              hlink_smp  = (TH1D*)((*vhRLRH_neg.at(ir))());
            std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhRLRH_sig.at(ir))()), (TH1D*)((*vhRLRH_bkg.at(ir))()) });

            //if (ir == 5 || ir == 6 || ir == 7) {
            //    ((TH1D*)((*vhRLRH_bkg.at(ir))()))->Reset();
            //    hlink_tmps = std::vector<TH1D*>({ (TH1D*)((*vhRLRH_sig.at(ir))()),  (TH1D*)((*vhRLRH_bkg.at(ir))()) });
            //}

            std::string prefix = Form("hLRH_T%03dR%03d_", it, ir);
            HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhRLRH_pos.at(ir))()), AX1D_LRH);
            HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_LRH, prefix, false, true);
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
            
            std::cerr << Form("LRH STATUS %d T%03d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
                fit1D.status(), it, ir, AXrig()(ir-1), AXrig()(ir),
                num_smp, num_sig, num_bkg,
                app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
                nchi);

            (*hTNchi)()->SetBinContent(it, ir, nchi); 
            
            (*hTNcntPr)()->SetBinContent(it, ir, num_ref); 
            (*hTNcntAp)()->SetBinContent(it, ir, num_sig); 
            
            (*hTErrNcnt)()->SetBinContent(it, ir, 100.0 * err_ncnt); 
            (*hTErrTmpl)()->SetBinContent(it, ir, 100.0 * err_tmpl);
            (*hTErrRscl)()->SetBinContent(it, ir, 100.0 * err_rscl);
            (*hTErrAccp)()->SetBinContent(it, ir, 100.0 * err_accp);
            (*hTErrCC  )()->SetBinContent(it, ir, 100.0 * err_cc  );
            
            (*hTErrStat)()->SetBinContent(it, ir, 100.0 * err_stat); 
            (*hTErrSyst)()->SetBinContent(it, ir, 100.0 * err_syst);
            (*hTErrTotl)()->SetBinContent(it, ir, 100.0 * err_totl);
            
            (*hTAccCorr)()->SetBinContent(it, ir, acc_corr); 
            
            (*hTAppStat)()->SetBinContent(it, ir, app_val); 
            (*hTAppStat)()->SetBinError  (it, ir, app_err_stat); 
            
            (*hTAppTotl)()->SetBinContent(it, ir, app_val);
            (*hTAppTotl)()->SetBinError  (it, ir, app_err_totl);

            Hist* hsmp = Hist::New(h1Dsmp.get());
            Hist* hsum = Hist::New(h1Dsum.get());
            Hist* hsig = Hist::New(h1Dsig.get());
            Hist* hbkg = Hist::New(h1Dbkg.get());
            
            hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
	    	hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
	    	hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
	    	hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
            
            /*
            editor_lrh.create();
   
            THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
            hfit->Draw("nostack hist");

            Hist* hcanvas = Hist::New(
                Form("%scanvas", prefix.c_str()), 
                HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.5 * hfit->GetHistogram()->GetMaximum())));
            (*hcanvas)()->GetXaxis()->SetTitle(AX1D_LRH.name_x().c_str());
            (*hcanvas)()->GetYaxis()->SetTitle(AX1D_LRH.name_y().c_str());
            (*hcanvas)()->Draw();
            
            hfit->Draw("nostack hist same");
            (*hsmp)()->Draw("pe same");
            
            Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
            leg_table()->SetHeader(Form("Time %3d Rigidity %.1f - %.1f [GV/c]", it, AXrig()(ir-1), AXrig()(ir)));
            leg_table()->AddEntry((*hsmp)(), "Data", "lp");
            leg_table()->AddEntry((*hsum)(), "Sum", "l");
            leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
            leg_table()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
            leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
            leg_table()->SetFillColor(0);
            leg_table.draw();
            
            editor_lrh.save();
            */
            
            Hist::Delete(hsmp);
            Hist::Delete(hsum);
            Hist::Delete(hsig);
            Hist::Delete(hbkg);
            
            fit1D.delete_hist();
        }
        
        for (auto&& hist : vhRLRH_pos) Hist::Delete(hist);
        for (auto&& hist : vhRLRH_neg) Hist::Delete(hist);
        //for (auto&& hist : vhRLRH_sig) Hist::Delete(hist);
        //for (auto&& hist : vhRLRH_bkg) Hist::Delete(hist);
    }
    for (auto&& hist : vhTLRH_pos) Hist::Delete(hist);
    for (auto&& hist : vhTLRH_neg) Hist::Delete(hist);
    for (auto&& hist : vhTLRH_sig) Hist::Delete(hist);
    for (auto&& hist : vhTLRH_bkg) Hist::Delete(hist);
        
    for (auto&& hist : vhRLRH_sig) Hist::Delete(hist);
    for (auto&& hist : vhRLRH_bkg) Hist::Delete(hist);

    //editor_lrh.close();

    /////////////////
    ////         ////
    ////   IIN   ////
    ////         ////
    /////////////////
    //PdfEditor editor_iin(Window(), "apflux_tme_iin", "out");
    
    TH1D* hIIN_acc = (TH1D*) file_acc->Get("hIin_accp");
    std::vector<Hist*> vhTIIN_pos = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sIinP_llr", subt.c_str())));
    std::vector<Hist*> vhTIIN_neg = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sIinN_llr", subt.c_str())));
    std::vector<Hist*> vhTIIN_sig = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sIinP_llr_pr", subt.c_str())));
    std::vector<Hist*> vhTIIN_bkg = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sIinN_llr_el", subt.c_str())));
    
    std::vector<Hist*> vhRIIN_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hRIinP_llr_pr"));
    std::vector<Hist*> vhRIIN_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hRIinN_llr_el"));
    
    //std::array<int, 2> RANGE_IIN({ 13, 17 });
    std::array<int, 2> RANGE_IIN({ 13, 17 });
    HistFit::Axis1D AX1D_IIN(
        "TRD Estimator",
        "Events/Bin",
        (*Hist::Head(Form("hT%sIinP_llr", subt.c_str())))()->GetZaxis());

    for (int it = 1; it <= AXtme.nbin(); ++it) {
        bool skip = false;
        for (auto&& itime : exclude_times) { if (it == itime) { skip = true; break; } }
        if (skip) continue;

        double tme = AXtme.center(it, AxisScale::kLinear);
        std::vector<Hist*> vhRIIN_pos = Hist::ProjectAll(HistProj::kY, vhTIIN_pos.at(it));
        std::vector<Hist*> vhRIIN_neg = Hist::ProjectAll(HistProj::kY, vhTIIN_neg.at(it));
        //std::vector<Hist*> vhRIIN_sig = Hist::ProjectAll(HistProj::kY, vhTIIN_sig.at(it));
        //std::vector<Hist*> vhRIIN_bkg = Hist::ProjectAll(HistProj::kY, vhTIIN_bkg.at(it));

        for (int ir = RANGE_IIN[0]; ir <= RANGE_IIN[1]; ++ir) {
            double rig = AXrig.center(ir, AxisScale::kLog);
            TH1D*              hlink_smp  = (TH1D*)((*vhRIIN_neg.at(ir))());
            std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhRIIN_sig.at(ir))()), (TH1D*)((*vhRIIN_bkg.at(ir))()) });

            std::string prefix = Form("hIIN_T%03dR%03d_", it, ir);
            HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhRIIN_pos.at(ir))()), AX1D_IIN);
            HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_IIN, prefix, false, true);
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
            
            std::cerr << Form("IIN STATUS %d T%3d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
                fit1D.status(), it, ir, AXrig()(ir-1), AXrig()(ir),
                num_smp, num_sig, num_bkg,
                app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
                nchi);

            (*hTNchi)()->SetBinContent(it, ir, nchi); 
            
            (*hTNcntPr)()->SetBinContent(it, ir, num_ref); 
            (*hTNcntAp)()->SetBinContent(it, ir, num_sig); 
            
            (*hTErrNcnt)()->SetBinContent(it, ir, 100.0 * err_ncnt); 
            (*hTErrTmpl)()->SetBinContent(it, ir, 100.0 * err_tmpl);
            (*hTErrRscl)()->SetBinContent(it, ir, 100.0 * err_rscl);
            (*hTErrAccp)()->SetBinContent(it, ir, 100.0 * err_accp);
            (*hTErrCC  )()->SetBinContent(it, ir, 100.0 * err_cc  );
            
            (*hTErrStat)()->SetBinContent(it, ir, 100.0 * err_stat); 
            (*hTErrSyst)()->SetBinContent(it, ir, 100.0 * err_syst);
            (*hTErrTotl)()->SetBinContent(it, ir, 100.0 * err_totl);
            
            (*hTAccCorr)()->SetBinContent(it, ir, acc_corr); 
            
            (*hTAppStat)()->SetBinContent(it, ir, app_val); 
            (*hTAppStat)()->SetBinError  (it, ir, app_err_stat); 
            
            (*hTAppTotl)()->SetBinContent(it, ir, app_val);
            (*hTAppTotl)()->SetBinError  (it, ir, app_err_totl);

            Hist* hsmp = Hist::New(h1Dsmp.get());
            Hist* hsum = Hist::New(h1Dsum.get());
            Hist* hsig = Hist::New(h1Dsig.get());
            Hist* hbkg = Hist::New(h1Dbkg.get());
            
            hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
	    	hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
	    	hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
	    	hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
     
            /*
            editor_iin.create();
   
            THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
            hfit->Draw("nostack hist");

            Hist* hcanvas = Hist::New(
                Form("%scanvas", prefix.c_str()), 
                HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.5 * hfit->GetHistogram()->GetMaximum())));
            (*hcanvas)()->GetXaxis()->SetTitle(AX1D_IIN.name_x().c_str());
            (*hcanvas)()->GetYaxis()->SetTitle(AX1D_IIN.name_y().c_str());
            (*hcanvas)()->Draw();
            
            hfit->Draw("nostack hist same");
            (*hsmp)()->Draw("pe same");
            
            Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
            leg_table()->SetHeader(Form("Time %3d Rigidity %.1f - %.1f [GV/c]", it, AXrig()(ir-1), AXrig()(ir)));
            leg_table()->AddEntry((*hsmp)(), "Data", "lp");
            leg_table()->AddEntry((*hsum)(), "Sum", "l");
            leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
            leg_table()->AddEntry((*hbkg)(), Form("e^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
            leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
            leg_table()->SetFillColor(0);
            leg_table.draw();
            
            editor_iin.save();
            */

            Hist::Delete(hsmp);
            Hist::Delete(hsum);
            Hist::Delete(hsig);
            Hist::Delete(hbkg);
            
            fit1D.delete_hist();
        }
        
        for (auto&& hist : vhRIIN_pos) Hist::Delete(hist);
        for (auto&& hist : vhRIIN_neg) Hist::Delete(hist);
        //for (auto&& hist : vhRIIN_sig) Hist::Delete(hist);
        //for (auto&& hist : vhRIIN_bkg) Hist::Delete(hist);
    }
    for (auto&& hist : vhTIIN_pos) Hist::Delete(hist);
    for (auto&& hist : vhTIIN_neg) Hist::Delete(hist);
    for (auto&& hist : vhTIIN_sig) Hist::Delete(hist);
    for (auto&& hist : vhTIIN_bkg) Hist::Delete(hist);

    for (auto&& hist : vhRIIN_sig) Hist::Delete(hist);
    for (auto&& hist : vhRIIN_bkg) Hist::Delete(hist);
    
    //editor_iin.close();


    /////////////////
    ////         ////
    ////   IEX   ////
    ////         ////
    /////////////////
    //PdfEditor editor_iex(Window(), "apflux_tme_iex", "out");
    
    TH1D* hIEX_acc = (TH1D*) file_acc->Get("hIex_accp");
    std::vector<Hist*> vhTIEX_pos = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sIexP_llr", subt.c_str())));
    std::vector<Hist*> vhTIEX_neg = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sIexN_llr", subt.c_str())));
    std::vector<Hist*> vhTIEX_sig = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sIexP_llr_pr", subt.c_str())));
    std::vector<Hist*> vhTIEX_bkg = Hist::ProjectAll(HistProj::kZY, Hist::Head(Form("hT%sIexN_llr_el", subt.c_str())));
    
    std::vector<Hist*> vhRIEX_sig = Hist::ProjectAll(HistProj::kY, Hist::Head("hRIexP_llr_pr"));
    std::vector<Hist*> vhRIEX_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head("hRIexN_llr_el"));
    
    //std::array<int, 2> RANGE_IEX({ 18, 21 });
    std::array<int, 2> RANGE_IEX({ 18, 22 });
    HistFit::Axis1D AX1D_IEX(
        "TRD Estimator",
        "Events/Bin",
        (*Hist::Head(Form("hT%sIexP_llr", subt.c_str())))()->GetZaxis());

    for (int it = 1; it <= AXtme.nbin(); ++it) {
        bool skip = false;
        for (auto&& itime : exclude_times) { if (it == itime) { skip = true; break; } }
        if (skip) continue;

        double tme = AXtme.center(it, AxisScale::kLinear);
        std::vector<Hist*> vhRIEX_pos = Hist::ProjectAll(HistProj::kY, vhTIEX_pos.at(it));
        std::vector<Hist*> vhRIEX_neg = Hist::ProjectAll(HistProj::kY, vhTIEX_neg.at(it));
        //std::vector<Hist*> vhRIEX_sig = Hist::ProjectAll(HistProj::kY, vhTIEX_sig.at(it));
        //std::vector<Hist*> vhRIEX_bkg = Hist::ProjectAll(HistProj::kY, vhTIEX_bkg.at(it));

        for (int ir = RANGE_IEX[0]; ir <= RANGE_IEX[1]; ++ir) {
            double rig = AXrig.center(ir, AxisScale::kLog);
            TH1D*              hlink_smp  = (TH1D*)((*vhRIEX_neg.at(ir))());
            std::vector<TH1D*> hlink_tmps({ (TH1D*)((*vhRIEX_sig.at(ir))()), (TH1D*)((*vhRIEX_bkg.at(ir))()) });

            std::string prefix = Form("hIEX_T%03dR%03d_", it, ir);
            HistFit::Hist1D    h1Dref(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhRIEX_pos.at(ir))()), AX1D_IEX);
            HistFit::HistFit1D fit1D(hlink_smp, hlink_tmps, AX1D_IEX, prefix, false, true);
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
            
            std::cerr << Form("IEX STATUS %d T%3d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
                fit1D.status(), it, ir, AXrig()(ir-1), AXrig()(ir),
                num_smp, num_sig, num_bkg,
                app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
                nchi);

            (*hTNchi)()->SetBinContent(it, ir, nchi); 
            
            (*hTNcntPr)()->SetBinContent(it, ir, num_ref); 
            (*hTNcntAp)()->SetBinContent(it, ir, num_sig); 
            
            (*hTErrNcnt)()->SetBinContent(it, ir, 100.0 * err_ncnt); 
            (*hTErrTmpl)()->SetBinContent(it, ir, 100.0 * err_tmpl);
            (*hTErrRscl)()->SetBinContent(it, ir, 100.0 * err_rscl);
            (*hTErrAccp)()->SetBinContent(it, ir, 100.0 * err_accp);
            (*hTErrCC  )()->SetBinContent(it, ir, 100.0 * err_cc  );
            
            (*hTErrStat)()->SetBinContent(it, ir, 100.0 * err_stat); 
            (*hTErrSyst)()->SetBinContent(it, ir, 100.0 * err_syst);
            (*hTErrTotl)()->SetBinContent(it, ir, 100.0 * err_totl);
            
            (*hTAccCorr)()->SetBinContent(it, ir, acc_corr); 
            
            (*hTAppStat)()->SetBinContent(it, ir, app_val); 
            (*hTAppStat)()->SetBinError  (it, ir, app_err_stat); 
            
            (*hTAppTotl)()->SetBinContent(it, ir, app_val);
            (*hTAppTotl)()->SetBinError  (it, ir, app_err_totl);

            Hist* hsmp = Hist::New(h1Dsmp.get());
            Hist* hsum = Hist::New(h1Dsum.get());
            Hist* hsig = Hist::New(h1Dsig.get());
            Hist* hbkg = Hist::New(h1Dbkg.get());
            
            hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
	    	hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
	    	hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
	    	hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
     
            /*
            editor_iex.create();
   
            THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hbkg, hsig }));
            hfit->Draw("nostack hist");

            Hist* hcanvas = Hist::New(
                Form("%scanvas", prefix.c_str()), 
                HistAxis(hsum->xaxis(), Axis("", 10000, hfit->GetHistogram()->GetMinimum(), 1.5 * hfit->GetHistogram()->GetMaximum())));
            (*hcanvas)()->GetXaxis()->SetTitle(AX1D_IEX.name_x().c_str());
            (*hcanvas)()->GetYaxis()->SetTitle(AX1D_IEX.name_y().c_str());
            (*hcanvas)()->Draw();
            
            hfit->Draw("nostack hist same");
            (*hsmp)()->Draw("pe same");
            
            Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
            leg_table()->SetHeader(Form("Time %3d Rigidity %.1f - %.1f [GV/c]", it, AXrig()(ir-1), AXrig()(ir)));
            leg_table()->AddEntry((*hsmp)(), "Data", "lp");
            leg_table()->AddEntry((*hsum)(), "Sum", "l");
            leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
            leg_table()->AddEntry((*hbkg)(), Form("e^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
            leg_table()->AddEntry((TObject*)0, Form("#chi^{2} %.2f", nchi), "");
            leg_table()->SetFillColor(0);
            leg_table.draw();
            
            editor_iex.save();
            */

            Hist::Delete(hsmp);
            Hist::Delete(hsum);
            Hist::Delete(hsig);
            Hist::Delete(hbkg);
            
            fit1D.delete_hist();
        }
        
        for (auto&& hist : vhRIEX_pos) Hist::Delete(hist);
        for (auto&& hist : vhRIEX_neg) Hist::Delete(hist);
        //for (auto&& hist : vhRIEX_sig) Hist::Delete(hist);
        //for (auto&& hist : vhRIEX_bkg) Hist::Delete(hist);
    }
    for (auto&& hist : vhTIEX_pos) Hist::Delete(hist);
    for (auto&& hist : vhTIEX_neg) Hist::Delete(hist);
    for (auto&& hist : vhTIEX_sig) Hist::Delete(hist);
    for (auto&& hist : vhTIEX_bkg) Hist::Delete(hist);

    for (auto&& hist : vhRIEX_sig) Hist::Delete(hist);
    for (auto&& hist : vhRIEX_bkg) Hist::Delete(hist);
    
    //editor_iex.close();



    PdfEditor editor(Window(), "apflux_tme", "out");

    editor.create();
    editor.cd(0, PadAxis(0, 1, 0));
    (*hTNchi)()->Draw("colz");
    editor.save();

    editor.create();
    editor.cd(0, PadAxis(0, 1, 1));
    (*hTNcntPr)()->Draw("colz");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1, 1));
    (*hTNcntAp)()->Draw("colz");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1, 0));
    (*hTErrStat)()->Draw("colz");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1, 0));
    (*hTErrSyst)()->Draw("colz");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1, 0));
    (*hTAppStat)()->Draw("colz");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1, 1));
    (*hTAppStat)()->Draw("colz");
    editor.save();


    editor.close();

    TFile * ofle = new TFile("out/apflux_tme.root", "RECREATE");
    ofle->cd();

    (*hTNchi)()->Write();
    
    (*hTNcntPr)()->Write();
    (*hTNcntAp)()->Write();
    
    (*hTErrNcnt)()->Write();
    (*hTErrTmpl)()->Write();
    (*hTErrRscl)()->Write();
    (*hTErrAccp)()->Write();
    (*hTErrCC  )()->Write();
    
    (*hTErrStat)()->Write();
    (*hTErrSyst)()->Write();
    (*hTErrTotl)()->Write();
    
    (*hTAccCorr)()->Write();
    (*hTAppStat)()->Write();
    (*hTAppTotl)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
