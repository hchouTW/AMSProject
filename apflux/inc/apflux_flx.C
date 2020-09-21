#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "HistFit1D.h"
#include "HistFit1D.C"

double get_error_cutoff(double arig) { return (0.04227366 * std::pow(arig, -0.00379330) * std::exp(-0.07817582 * arig)); }
double get_error_selection(double arig) { return 0.01 * (0.5 + 1.5 / (1.0 + 200.0 / arig)); }

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

TF1* stdgaus = new TF1("stdgaus", "gaus", 0, 10.0);

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "08";
    std::string subg = "";
    //subg = "_CF";
    
    bool sw_fluc = true; 
    bool sw_build_hist = true; 
    
    TGraphErrors* offapp = (TGraphErrors*) (TFile::Open("others/20160406MIT.root")->Get("gpbarp"));
    offapp->SetLineColor(kBlue);
    offapp->SetMarkerColor(kBlue);

    TFile* file_acc = TFile::Open("out/apflux_acc.root");
    TH1D* haccp_smooth = (TH1D*) file_acc->Get("hist_accptance_factor_smooth");
    TH1D* haccp = (TH1D*) file_acc->Get("hist_accptance_factor");
    TH1D* hacce = (TH1D*) file_acc->Get("hist_accptance_error");
    TH1D* hrfnc = (TH1D*) file_acc->Get("hist_rigfunc_error");
    TH1D* hmccc = (TH1D*) file_acc->Get("hist_mccc");
    
    TFile* file_mc = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr_l1o9flux%s/YiMdst.root", subv.c_str()));

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/iss%s", subv.c_str()));

    const Axis& AXrig = Hist::Head("hLtfP_sqrm")->xaxis();
    
    Hist* hNchi = Hist::New("hNchi", HistAxis(AXrig, "#chi^{2}/NDF"));
    
    Hist* hNcntPr = Hist::New("hNcntPr", HistAxis(AXrig, "Events/Bin"));
    Hist* hNcntAp = Hist::New("hNcntAp", HistAxis(AXrig, "Events/Bin"));
    
    Hist* hAccCorr = Hist::New("hAccCorr", HistAxis(AXrig, "Effective Acceptance Correction Factor (A^{#bar{p}/p}_{CF})"));

    Hist* hErrCfft = Hist::New("hErrCfft", HistAxis(AXrig, "RelError of Cutoff (%)"));
    Hist* hErrTmpl = Hist::New("hErrTmpl", HistAxis(AXrig, "RelError of Template (%)"));
    Hist* hErrNcnt = Hist::New("hErrNcnt", HistAxis(AXrig, "RelError of Count (%)"));
    Hist* hErrRfnc = Hist::New("hErrRfnc", HistAxis(AXrig, "RelError of Rigidity (%)"));
    Hist* hErrAccp = Hist::New("hErrAccp", HistAxis(AXrig, "RelError of Acceptance (%)"));
    
    Hist* hErrStat = Hist::New("hErrStat", HistAxis(AXrig, "RelError of Stat (%)"));
    Hist* hErrSyst = Hist::New("hErrSyst", HistAxis(AXrig, "RelError of Syst (%)"));
    Hist* hErrTotl = Hist::New("hErrTotl", HistAxis(AXrig, "RelError of Totl (%)"));
    
    Hist* hAppStat = Hist::New("hAppStat", HistAxis(AXrig, "#bar{p}/p"));
    Hist* hAppTotl = Hist::New("hAppTotl", HistAxis(AXrig, "#bar{p}/p"));
    
    Hist* hAppCC   = Hist::New("hAppCC",   HistAxis(AXrig, "p_{CC}/p"));

    PdfEditor editor(Window(WindowSize::kWideSliceLR), Form("apflux_flx%s", subg.c_str()), "out");

    /////////////////
    ////         ////
    ////   LTF   ////
    ////         ////
    /////////////////
    std::vector<Hist*> vhLTF_pos = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hLtfP%s_sqrm", subg.c_str())));
    std::vector<Hist*> vhLTF_neg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hLtfN%s_sqrm", subg.c_str())));
    std::vector<Hist*> vhLTF_sig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hLtfP%s_sqrm_pr", subg.c_str())));
    std::vector<Hist*> vhLTF_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hLtfN%s_sqrm_el", subg.c_str())));
    
    std::array<int, 2> RANGE_LTF({ 1, 11 });

    for (int ir = RANGE_LTF[0]; ir <= RANGE_LTF[1]; ++ir) {
        HistFit::Axis1D AX1D_LTF(
            "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
            "Events/Bin",
            (*Hist::Head(Form("hLtfP%s_sqrm", subg.c_str())))()->GetYaxis());

        double rig = AXrig.center(ir, AxisScale::kLog);
        double acc_corr = haccp->Interpolate(rig);
        
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
        
        const int MAX_ii = 60;
        std::array<double, MAX_ii> hfitdv;
        Hist* hfitsqrm = Hist::New(Form("%sfitsqrm", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hLtfP%s_sqrm", subg.c_str())))()->GetYaxis())));
        Hist* hchisqrm = Hist::New(Form("%schisqrm", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hLtfP%s_sqrm", subg.c_str())))()->GetYaxis())));
        for (int ii = 1; ii <= MAX_ii; ++ii) {
            HistFit::Axis1D AX1D_LTFi(
                "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
                "Events/Bin",
                (*Hist::Head(Form("hLtfP%s_sqrm", subg.c_str())))()->GetYaxis(), ii, 100);
            HistFit::HistFit1D fit1Di(hlink_smp, hlink_tmps, AX1D_LTFi, prefix, false, false);
            HistFit::Hist1D    h1Drefi(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhLTF_pos.at(ir))()), AX1D_LTFi);
            (*hfitsqrm)()->SetBinContent(ii, acc_corr * fit1Di.wgts(0) / h1Drefi.data().sum());
            (*hfitsqrm)()->SetBinError  (ii, (*hfitsqrm)()->GetBinContent(ii) * fit1Di.errs(0) / fit1Di.wgts(0));
            (*hchisqrm)()->SetBinContent(ii, fit1Di.nchi());
            hfitdv.at(ii-1) = (*hfitsqrm)()->GetBinContent(ii) / (acc_corr * (num_sig / num_ref)) - 1.0;
            hfitdv.at(ii-1) = hfitdv.at(ii-1) * hfitdv.at(ii-1);
        }
        (*hfitsqrm)()->GetXaxis()->SetTitle("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        (*hchisqrm)()->GetXaxis()->SetTitle("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        (*hfitsqrm)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
        (*hchisqrm)()->GetYaxis()->SetTitle("#chi^{2}/NDF");
        (*hfitsqrm)()->GetXaxis()->SetRange(1, MAX_ii);
        (*hchisqrm)()->GetXaxis()->SetRange(1, MAX_ii);
        double err_cuts = std::sqrt(std::accumulate(hfitdv.begin(), hfitdv.end(), 0.0) / hfitdv.size());

        double err_cfft = get_error_cutoff(rig);
        double err_tmpl = std::hypot(fit1D.fluc(0).err / num_sig, get_error_selection(rig));
        double err_ncnt = std::hypot(err_cfft, err_tmpl);
        double err_rscl = 0.01 * hrfnc->Interpolate(rig);
        double err_accp = 0.01 * hacce->Interpolate(rig);

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_rscl*err_rscl + err_accp*err_accp);
        double err_totl = std::hypot(err_stat, err_syst);

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
        
        (*hErrCfft)()->SetBinContent(ir, 100.0 * err_cfft);
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrRfnc)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        
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
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2}/NDF %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
       
        if (sw_fluc) {
            const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
            Hist* hfluc = Hist::New(h1Dfluc.get()->GetName(), 
                                    HistAxis(Axis("#bar{p}/p Flux Ratio [10^{-4}]", 
                                                  h1Dfluc.get()->GetXaxis()->GetNbins(), 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmin() / num_ref, 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmax() / num_ref), "Events/Bin"));
            for (int ib = 1; ib <= hfluc->xaxis().nbin(); ++ib) {
                (*hfluc)()->SetBinContent(ib, h1Dfluc.get()->GetBinContent(ib));
                (*hfluc)()->SetBinError  (ib, h1Dfluc.get()->GetBinError  (ib));
            }
            stdgaus->SetParameters(1.0, app_val, app_val * err_tmpl);
            stdgaus->SetLineColor(kBlue);
            stdgaus->SetNpx(100000);
            (*hfluc)()->Fit(stdgaus, "q0", "", hfluc->xaxis().min(), hfluc->xaxis().max());

            editor.create();
		    hfluc->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
            hfluc->draw("pe");
            stdgaus->Draw("l same");
            Legend leg_fluc("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.65, 0.85));
            leg_fluc()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
            leg_fluc()->AddEntry((*hfluc)(), "Data", "lp");
            leg_fluc()->AddEntry(stdgaus, Form("Fit  #sigma_{tmpl.} = %.2f%", 100.0 * err_tmpl), "l");
            leg_fluc()->SetFillColor(0);
            leg_fluc.draw();
            editor.save();
        }
        
        editor.create();
		hfitsqrm->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfitsqrm->draw("pe");
        editor.save();
        
        editor.create();
		hchisqrm->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hchisqrm->draw("hist");
        editor.save();
    }


    /////////////////
    ////         ////
    ////   LRH   ////
    ////         ////
    /////////////////
    std::vector<Hist*> vhLRH_pos = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hLrhP%s_sqrm", subg.c_str())));
    std::vector<Hist*> vhLRH_neg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hLrhN%s_sqrm", subg.c_str())));
    std::vector<Hist*> vhLRH_sig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hLrhP%s_sqrm_pr", subg.c_str())));
    std::vector<Hist*> vhLRH_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hLrhN%s_sqrm_el", subg.c_str())));
    
    std::array<int, 2> RANGE_LRH({ 12, 25 });

    for (int ir = RANGE_LRH[0]; ir <= RANGE_LRH[1]; ++ir) {
        HistFit::Axis1D AX1D_LRH(
            "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
            "Events/Bin",
            (*Hist::Head(Form("hLrhP%s_sqrm", subg.c_str())))()->GetYaxis());

        double rig = AXrig.center(ir, AxisScale::kLog);
        double acc_corr = haccp->Interpolate(rig);
        
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
        
        const int MAX_ii = 60;
        std::array<double, MAX_ii> hfitdv;
        Hist* hfitsqrm = Hist::New(Form("%sfitsqrm", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hLrhP%s_sqrm", subg.c_str())))()->GetYaxis())));
        Hist* hchisqrm = Hist::New(Form("%schisqrm", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hLrhP%s_sqrm", subg.c_str())))()->GetYaxis())));
        for (int ii = 1; ii <= MAX_ii; ++ii) {
            HistFit::Axis1D AX1D_LRHi(
                "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
                "Events/Bin",
                (*Hist::Head(Form("hLrhP%s_sqrm", subg.c_str())))()->GetYaxis(), ii, 125);
            HistFit::HistFit1D fit1Di(hlink_smp, hlink_tmps, AX1D_LRHi, prefix, false, false);
            HistFit::Hist1D    h1Drefi(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhLRH_pos.at(ir))()), AX1D_LRHi);
            (*hfitsqrm)()->SetBinContent(ii, acc_corr * fit1Di.wgts(0) / h1Drefi.data().sum());
            (*hfitsqrm)()->SetBinError  (ii, (*hfitsqrm)()->GetBinContent(ii) * fit1Di.errs(0) / fit1Di.wgts(0));
            (*hchisqrm)()->SetBinContent(ii, fit1Di.nchi());
            hfitdv.at(ii-1) = (*hfitsqrm)()->GetBinContent(ii) / (acc_corr * (num_sig / num_ref)) - 1.0;
            hfitdv.at(ii-1) = hfitdv.at(ii-1) * hfitdv.at(ii-1);
        }
        (*hfitsqrm)()->GetXaxis()->SetTitle("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        (*hchisqrm)()->GetXaxis()->SetTitle("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        (*hfitsqrm)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
        (*hchisqrm)()->GetYaxis()->SetTitle("#chi^{2}/NDF");
        (*hfitsqrm)()->GetXaxis()->SetRange(1, MAX_ii);
        (*hchisqrm)()->GetXaxis()->SetRange(1, MAX_ii);
        double err_cuts = std::sqrt(std::accumulate(hfitdv.begin(), hfitdv.end(), 0.0) / hfitdv.size());

        double err_cfft = get_error_cutoff(rig);
        double err_tmpl = std::hypot(fit1D.fluc(0).err / num_sig, get_error_selection(rig));
        double err_ncnt = std::hypot(err_cfft, err_tmpl);
        double err_rscl = 0.01 * hrfnc->Interpolate(rig);
        double err_accp = 0.01 * hacce->Interpolate(rig);

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_rscl*err_rscl + err_accp*err_accp);
        double err_totl = std::hypot(err_stat, err_syst);

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
        
        (*hErrCfft)()->SetBinContent(ir, 100.0 * err_cfft);
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrRfnc)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        
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
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.58, 0.85, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [Gc]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2}/NDF %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        if (sw_fluc) {
            const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
            Hist* hfluc = Hist::New(h1Dfluc.get()->GetName(), 
                                    HistAxis(Axis("#bar{p}/p Flux Ratio [10^{-4}]", 
                                                  h1Dfluc.get()->GetXaxis()->GetNbins(), 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmin() / num_ref, 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmax() / num_ref), "Events/Bin"));
            for (int ib = 1; ib <= hfluc->xaxis().nbin(); ++ib) {
                (*hfluc)()->SetBinContent(ib, h1Dfluc.get()->GetBinContent(ib));
                (*hfluc)()->SetBinError  (ib, h1Dfluc.get()->GetBinError  (ib));
            }
            stdgaus->SetParameters(1.0, app_val, app_val * err_tmpl);
            stdgaus->SetLineColor(kBlue);
            stdgaus->SetNpx(100000);
            (*hfluc)()->Fit(stdgaus, "q0", "", hfluc->xaxis().min(), hfluc->xaxis().max());

            editor.create();
		    hfluc->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
            hfluc->draw("pe");
            stdgaus->Draw("l same");
            Legend leg_fluc("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.65, 0.85));
            leg_fluc()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
            leg_fluc()->AddEntry((*hfluc)(), "Data", "lp");
            leg_fluc()->AddEntry(stdgaus, Form("Fit  #sigma_{tmpl.} = %.2f%", 100.0 * err_tmpl), "l");
            leg_fluc()->SetFillColor(0);
            leg_fluc.draw();
            editor.save();
        }
        
        editor.create();
		hfitsqrm->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfitsqrm->draw("pe");
        editor.save();
        
        editor.create();
		hchisqrm->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hchisqrm->draw("hist");
        editor.save();
    }



    /////////////////
    ////         ////
    ////   IIN   ////
    ////         ////
    /////////////////
    std::vector<Hist*> vhIIN_pos = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hIinP%s_llr", subg.c_str())));
    std::vector<Hist*> vhIIN_neg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hIinN%s_llr", subg.c_str())));
    std::vector<Hist*> vhIIN_sig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hIinP%s_llr_pr", subg.c_str())));
    std::vector<Hist*> vhIIN_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hIinN%s_llr_el", subg.c_str())));
    
    std::array<int, 2> RANGE_IIN({ 26, 36 });

    for (int ir = RANGE_IIN[0]; ir <= RANGE_IIN[1]; ++ir) {
        HistFit::Axis1D AX1D_IIN(
            "TRD Estimator (#Lambda_{TRD})",
            "Events/Bin",
            (*Hist::Head(Form("hIinP%s_llr", subg.c_str())))()->GetYaxis());

        double rig = AXrig.center(ir, AxisScale::kLog);
        double acc_corr = haccp->Interpolate(rig);
        
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
        
        const int MAX_ii = 80;
        std::array<double, MAX_ii> hfitdv;
        Hist* hfitllr = Hist::New(Form("%sfitllr", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hIinP%s_llr", subg.c_str())))()->GetYaxis())));
        Hist* hchillr = Hist::New(Form("%schillr", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hIinP%s_llr", subg.c_str())))()->GetYaxis())));
        for (int ii = 1; ii <= MAX_ii; ++ii) {
            HistFit::Axis1D AX1D_IINi(
                "TRD Estimator (#Lambda_{TRD})",
                "Events/Bin",
                (*Hist::Head(Form("hIinP%s_llr", subg.c_str())))()->GetYaxis(), ii, 150);
            HistFit::HistFit1D fit1Di(hlink_smp, hlink_tmps, AX1D_IINi, prefix, false, false);
            HistFit::Hist1D    h1Drefi(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhIIN_pos.at(ir))()), AX1D_IINi);
            (*hfitllr)()->SetBinContent(ii, acc_corr * fit1Di.wgts(0) / h1Drefi.data().sum());
            (*hfitllr)()->SetBinError  (ii, (*hfitllr)()->GetBinContent(ii) * fit1Di.errs(0) / fit1Di.wgts(0));
            (*hchillr)()->SetBinContent(ii, fit1Di.nchi());
            hfitdv.at(ii-1) = (*hfitllr)()->GetBinContent(ii) / (acc_corr * (num_sig / num_ref)) - 1.0;
            hfitdv.at(ii-1) = hfitdv.at(ii-1) * hfitdv.at(ii-1);
        }
        (*hfitllr)()->GetXaxis()->SetTitle("TRD Estimator (#Lambda_{TRD})");
        (*hchillr)()->GetXaxis()->SetTitle("TRD Estimator (#Lambda_{TRD})");
        (*hfitllr)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
        (*hchillr)()->GetYaxis()->SetTitle("#chi^{2}/NDF");
        (*hfitllr)()->GetXaxis()->SetRange(1, MAX_ii);
        (*hchillr)()->GetXaxis()->SetRange(1, MAX_ii);
        double err_cuts = std::sqrt(std::accumulate(hfitdv.begin(), hfitdv.end(), 0.0) / hfitdv.size());

        double err_cfft = get_error_cutoff(rig);
        double err_tmpl = std::hypot(fit1D.fluc(0).err / num_sig, get_error_selection(rig));
        double err_ncnt = std::hypot(err_cfft, err_tmpl);
        double err_rscl = 0.01 * hrfnc->Interpolate(rig);
        double err_accp = 0.01 * hacce->Interpolate(rig);

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_rscl*err_rscl + err_accp*err_accp);
        double err_totl = std::hypot(err_stat, err_syst);

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
        
        (*hErrCfft)()->SetBinContent(ir, 100.0 * err_cfft);
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrRfnc)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        
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
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.58, 0.85, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("e^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2}/NDF %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        if (sw_fluc) {
            const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
            Hist* hfluc = Hist::New(h1Dfluc.get()->GetName(), 
                                    HistAxis(Axis("#bar{p}/p Flux Ratio [10^{-4}]", 
                                                  h1Dfluc.get()->GetXaxis()->GetNbins(), 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmin() / num_ref, 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmax() / num_ref), "Events/Bin"));
            for (int ib = 1; ib <= hfluc->xaxis().nbin(); ++ib) {
                (*hfluc)()->SetBinContent(ib, h1Dfluc.get()->GetBinContent(ib));
                (*hfluc)()->SetBinError  (ib, h1Dfluc.get()->GetBinError  (ib));
            }
            stdgaus->SetParameters(1.0, app_val, app_val * err_tmpl);
            stdgaus->SetLineColor(kBlue);
            stdgaus->SetNpx(100000);
            (*hfluc)()->Fit(stdgaus, "q0", "", hfluc->xaxis().min(), hfluc->xaxis().max());

            editor.create();
		    hfluc->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
            hfluc->draw("pe");
            stdgaus->Draw("l same");
            Legend leg_fluc("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.65, 0.85));
            leg_fluc()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
            leg_fluc()->AddEntry((*hfluc)(), "Data", "lp");
            leg_fluc()->AddEntry(stdgaus, Form("Fit  #sigma_{tmpl.} = %.2f%", 100.0 * err_tmpl), "l");
            leg_fluc()->SetFillColor(0);
            leg_fluc.draw();
            editor.save();
        }
        
        editor.create();
		hfitllr->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfitllr->draw("pe");
        editor.save();
        
        editor.create();
		hchillr->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hchillr->draw("hist");
        editor.save();
    }




    /////////////////
    ////         ////
    ////   IEX   ////
    ////         ////
    /////////////////
    std::vector<Hist*> vhIEX_pos = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hIexP%s_llr", subg.c_str())));
    std::vector<Hist*> vhIEX_neg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hIexN%s_llr", subg.c_str())));
    std::vector<Hist*> vhIEX_sig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hIexP%s_llr_pr", subg.c_str())));
    std::vector<Hist*> vhIEX_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hIexN%s_llr_el", subg.c_str())));
    
    std::array<int, 2> RANGE_IEX({ 37, 45 });

    for (int ir = RANGE_IEX[0]; ir <= RANGE_IEX[1]; ++ir) {
        HistFit::Axis1D AX1D_IEX(
            "TRD Estimator (#Lambda_{TRD})",
            "Events/Bin",
            (*Hist::Head(Form("hIexP%s_llr", subg.c_str())))()->GetYaxis());

        double rig = AXrig.center(ir, AxisScale::kLog);
        double acc_corr = haccp->Interpolate(rig);
        
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
        
        const int MAX_ii = 60;
        std::array<double, MAX_ii> hfitdv;
        Hist* hfitllr = Hist::New(Form("%sfitllr", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hIexP%s_llr", subg.c_str())))()->GetYaxis())));
        Hist* hchillr = Hist::New(Form("%schillr", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hIexP%s_llr", subg.c_str())))()->GetYaxis())));
        for (int ii = 1; ii <= MAX_ii; ++ii) {
            HistFit::Axis1D AX1D_IEXi(
                "TRD Estimator (#Lambda_{TRD})",
                "Events/Bin",
                (*Hist::Head(Form("hIexP%s_llr", subg.c_str())))()->GetYaxis(), ii, 100);
            HistFit::HistFit1D fit1Di(hlink_smp, hlink_tmps, AX1D_IEXi, prefix, false, false);
            HistFit::Hist1D    h1Drefi(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhIEX_pos.at(ir))()), AX1D_IEXi);
            (*hfitllr)()->SetBinContent(ii, acc_corr * fit1Di.wgts(0) / h1Drefi.data().sum());
            (*hfitllr)()->SetBinError  (ii, (*hfitllr)()->GetBinContent(ii) * fit1Di.errs(0) / fit1Di.wgts(0));
            (*hchillr)()->SetBinContent(ii, fit1Di.nchi());
            hfitdv.at(ii-1) = (*hfitllr)()->GetBinContent(ii) / (acc_corr * (num_sig / num_ref)) - 1.0;
            hfitdv.at(ii-1) = hfitdv.at(ii-1) * hfitdv.at(ii-1);
        }
        (*hfitllr)()->GetXaxis()->SetTitle("TRD Estimator (#Lambda_{TRD})");
        (*hchillr)()->GetXaxis()->SetTitle("TRD Estimator (#Lambda_{TRD})");
        (*hfitllr)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
        (*hchillr)()->GetYaxis()->SetTitle("#chi^{2}/NDF");
        (*hfitllr)()->GetXaxis()->SetRange(1, MAX_ii);
        (*hchillr)()->GetXaxis()->SetRange(1, MAX_ii);
        double err_cuts = std::sqrt(std::accumulate(hfitdv.begin(), hfitdv.end(), 0.0) / hfitdv.size());

        double err_cfft = get_error_cutoff(rig);
        double err_tmpl = std::hypot(fit1D.fluc(0).err / num_sig, get_error_selection(rig));
        double err_ncnt = std::hypot(err_cfft, err_tmpl);
        double err_rscl = 0.01 * hrfnc->Interpolate(rig);
        double err_accp = 0.01 * hacce->Interpolate(rig);

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_rscl*err_rscl + err_accp*err_accp);
        double err_totl = std::hypot(err_stat, err_syst);

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
        
        (*hErrCfft)()->SetBinContent(ir, 100.0 * err_cfft);
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrRfnc)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        
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
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.58, 0.85, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("e^{-} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2}/NDF %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        if (sw_fluc) {
            const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
            Hist* hfluc = Hist::New(h1Dfluc.get()->GetName(), 
                                    HistAxis(Axis("#bar{p}/p Flux Ratio [10^{-4}]", 
                                                  h1Dfluc.get()->GetXaxis()->GetNbins(), 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmin() / num_ref, 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmax() / num_ref), "Events/Bin"));
            for (int ib = 1; ib <= hfluc->xaxis().nbin(); ++ib) {
                (*hfluc)()->SetBinContent(ib, h1Dfluc.get()->GetBinContent(ib));
                (*hfluc)()->SetBinError  (ib, h1Dfluc.get()->GetBinError  (ib));
            }
            stdgaus->SetParameters(1.0, app_val, app_val * err_tmpl);
            stdgaus->SetLineColor(kBlue);
            stdgaus->SetNpx(100000);
            (*hfluc)()->Fit(stdgaus, "q0", "", hfluc->xaxis().min(), hfluc->xaxis().max());

            editor.create();
		    hfluc->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
            hfluc->draw("pe");
            stdgaus->Draw("l same");
            Legend leg_fluc("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.65, 0.85));
            leg_fluc()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
            leg_fluc()->AddEntry((*hfluc)(), "Data", "lp");
            leg_fluc()->AddEntry(stdgaus, Form("Fit  #sigma_{tmpl.} = %.2f%", 100.0 * err_tmpl), "l");
            leg_fluc()->SetFillColor(0);
            leg_fluc.draw();
            editor.save();
        }
        
        editor.create();
		hfitllr->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfitllr->draw("pe");
        editor.save();
        
        editor.create();
		hchillr->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hchillr->draw("hist");
        editor.save();
    }



    /////////////////
    ////         ////
    ////   HEX   ////
    ////         ////
    /////////////////
    Hist* hHEX_pr = Hist::New(Form("hHex%s_cc_pr", subg.c_str()), (TH2D*)((*Hist::Head(Form("hHexP%s_cc", subg.c_str())))()));
    Hist* hHEX_cc = Hist::New(Form("hHex%s_cc_cc", subg.c_str()), (TH2D*)(file_mc->Get(Form("hHexN_cc"))));
    std::vector<Hist*> vhHEX_pos = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHexP%s_cc", subg.c_str())));
    std::vector<Hist*> vhHEX_neg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHexN%s_cc", subg.c_str())));
    std::vector<Hist*> vhHEX_sig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHex%s_cc_pr", subg.c_str())));
    std::vector<Hist*> vhHEX_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHex%s_cc_cc", subg.c_str())));
    
    std::array<int, 2> RANGE_HEX({ 46, 55 });
    std::map<int, int> FITREG_HEX;
    FITREG_HEX[46] = 20;
    FITREG_HEX[47] = 10;
    FITREG_HEX[48] = 10;
    FITREG_HEX[49] = 20;
    FITREG_HEX[50] = 40;
    FITREG_HEX[51] = 30;
    FITREG_HEX[52] = 20;
    FITREG_HEX[53] = 10;
    FITREG_HEX[54] = 20;
    FITREG_HEX[55] = 10;

    for (int ir = RANGE_HEX[0]; ir <= RANGE_HEX[1]; ++ir) {
        HistFit::Axis1D AX1D_HEX(
            "CC Estimator (#Lambda_{CC})",
            "Events/Bin",
            (*Hist::Head(Form("hHexP%s_cc", subg.c_str())))()->GetYaxis(), FITREG_HEX[ir], 100);

        double rig = AXrig.center(ir, AxisScale::kLog);
        double acc_corr = haccp->Interpolate(rig);
        
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
        
        const int MAX_ii = 50+30;
        std::array<double, MAX_ii> hfitdv;
        Hist* hfitcc = Hist::New(Form("%sfitcc", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hHexP%s_cc", subg.c_str())))()->GetYaxis())));
        Hist* hchicc = Hist::New(Form("%schicc", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hHexP%s_cc", subg.c_str())))()->GetYaxis())));
        Hist* hflccc = Hist::New(Form("%sflccc", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hHexP%s_cc", subg.c_str())))()->GetYaxis())));
        for (int ii = 1; ii <= MAX_ii; ++ii) {
            HistFit::Axis1D AX1D_HEXi(
                "CC Estimator (#Lambda_{CC})",
                "Events/Bin",
                (*Hist::Head(Form("hHexP%s_cc", subg.c_str())))()->GetYaxis(), ii, 100);
            HistFit::HistFit1D fit1Di(hlink_smp, hlink_tmps, AX1D_HEXi, prefix, false, false);
            //HistFit::HistFit1D fit1Di(hlink_smp, hlink_tmps, AX1D_HEXi, prefix, false, false);
            HistFit::Hist1D    h1Drefi(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhHEX_pos.at(ir))()), AX1D_HEXi);
            (*hfitcc)()->SetBinContent(ii, acc_corr * fit1Di.wgts(0) / h1Drefi.data().sum());
            (*hfitcc)()->SetBinError  (ii, (*hfitcc)()->GetBinContent(ii) * fit1Di.errs(0) / fit1Di.wgts(0));
            (*hchicc)()->SetBinContent(ii, fit1Di.nchi());
            (*hflccc)()->SetBinContent(ii, 100.0 * (fit1Di.fluc(0).err / fit1Di.wgts(0)));
            hfitdv.at(ii-1) = (*hfitcc)()->GetBinContent(ii) / (acc_corr * (num_sig / num_ref)) - 1.0;
            hfitdv.at(ii-1) = hfitdv.at(ii-1) * hfitdv.at(ii-1);
        }
        (*hfitcc)()->GetXaxis()->SetTitle("CC Estimator (#Lambda_{CC})");
        (*hchicc)()->GetXaxis()->SetTitle("CC Estimator (#Lambda_{CC})");
        (*hflccc)()->GetXaxis()->SetTitle("CC Estimator (#Lambda_{CC})");
        (*hfitcc)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
        (*hchicc)()->GetYaxis()->SetTitle("#chi^{2}/NDF");
        (*hflccc)()->GetYaxis()->SetTitle("Relative Error (%)");
        (*hfitcc)()->GetXaxis()->SetRange(1, MAX_ii);
        (*hchicc)()->GetXaxis()->SetRange(1, MAX_ii);
        (*hflccc)()->GetXaxis()->SetRange(1, MAX_ii);
        double err_cuts = std::sqrt(std::accumulate(hfitdv.begin(), hfitdv.end(), 0.0) / hfitdv.size());

        double err_cfft = get_error_cutoff(rig);
        double err_tmpl = std::hypot(fit1D.fluc(0).err / num_sig, get_error_selection(rig));
        double err_ncnt = std::hypot(err_cfft, err_tmpl);
        double err_rscl = 0.01 * hrfnc->Interpolate(rig);
        double err_accp = 0.01 * hacce->Interpolate(rig);

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_rscl*err_rscl + err_accp*err_accp);
        double err_totl = std::hypot(err_stat, err_syst);

        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
        
        double app_cc_val = (fit1D.wgts(1) / num_ref);
        double app_cc_err = (fit1D.errs(1) / num_ref);
        
        std::cerr << Form("HEX STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
            fit1D.status(), ir, AXrig()(ir-1), AXrig()(ir),
            num_smp, num_sig, num_bkg,
            app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
            nchi);

        (*hNchi)()->SetBinContent(ir, nchi); 
        
        (*hNcntPr)()->SetBinContent(ir, num_ref); 
        (*hNcntAp)()->SetBinContent(ir, num_sig); 
        
        (*hErrCfft)()->SetBinContent(ir, 100.0 * err_cfft);
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrRfnc)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        
        (*hErrStat)()->SetBinContent(ir, 100.0 * err_stat); 
        (*hErrSyst)()->SetBinContent(ir, 100.0 * err_syst);
        (*hErrTotl)()->SetBinContent(ir, 100.0 * err_totl);
        
        (*hAccCorr)()->SetBinContent(ir, acc_corr);

        (*hAppStat)()->SetBinContent(ir, app_val); 
        (*hAppStat)()->SetBinError  (ir, app_err_stat); 
        
        (*hAppTotl)()->SetBinContent(ir, app_val);
        (*hAppTotl)()->SetBinError  (ir, app_err_totl);
        
        (*hAppCC)()->SetBinContent(ir, app_cc_val);
        (*hAppCC)()->SetBinError  (ir, app_cc_err);

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
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.35, 0.65, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("p^{+} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2}/NDF %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        if (sw_fluc) {
            const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
            Hist* hfluc = Hist::New(h1Dfluc.get()->GetName(), 
                                    HistAxis(Axis("#bar{p}/p Flux Ratio [10^{-4}]", 
                                                  h1Dfluc.get()->GetXaxis()->GetNbins(), 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmin() / num_ref, 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmax() / num_ref), "Events/Bin"));
            for (int ib = 1; ib <= hfluc->xaxis().nbin(); ++ib) {
                (*hfluc)()->SetBinContent(ib, h1Dfluc.get()->GetBinContent(ib));
                (*hfluc)()->SetBinError  (ib, h1Dfluc.get()->GetBinError  (ib));
            }
            stdgaus->SetParameters(1.0, app_val, app_val * err_tmpl);
            stdgaus->SetLineColor(kBlue);
            stdgaus->SetNpx(100000);
            (*hfluc)()->Fit(stdgaus, "q0", "", hfluc->xaxis().min(), hfluc->xaxis().max());

            editor.create();
		    hfluc->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
            hfluc->draw("pe");
            stdgaus->Draw("l same");
            Legend leg_fluc("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.65, 0.85));
            leg_fluc()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
            leg_fluc()->AddEntry((*hfluc)(), "Data", "lp");
            leg_fluc()->AddEntry(stdgaus, Form("Fit  #sigma_{tmpl.} = %.2f%", 100.0 * err_tmpl), "l");
            leg_fluc()->SetFillColor(0);
            leg_fluc.draw();
            editor.save();
        }
        
        editor.create();
		hfitcc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfitcc->draw("pe");
        editor.save();
        
        editor.create();
		hchicc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hchicc->draw("hist");
        editor.save();
        
        editor.create();
		hflccc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hflccc->draw("hist");
        editor.save();
    }



    /////////////////
    ////         ////
    ////   HFS   ////
    ////         ////
    /////////////////
    Hist* hHFS_pr = Hist::New(Form("hHfs%s_cc_pr", subg.c_str()), (TH2D*)((*Hist::Head(Form("hHfsP%s_cc", subg.c_str())))()));
    Hist* hHFS_cc = Hist::New(Form("hHfs%s_cc_cc", subg.c_str()), (TH2D*)(file_mc->Get(Form("hHfsN_cc"))));
    std::vector<Hist*> vhHFS_pos = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHfsP%s_cc", subg.c_str())));
    std::vector<Hist*> vhHFS_neg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHfsN%s_cc", subg.c_str())));
    std::vector<Hist*> vhHFS_sig = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHfs%s_cc_pr", subg.c_str())));
    std::vector<Hist*> vhHFS_bkg = Hist::ProjectAll(HistProj::kY, Hist::Head(Form("hHfs%s_cc_cc", subg.c_str())));
    
    std::array<int, 2> RANGE_HFS({ 56, 58 });
    std::map<int, int> FITREG_HFS;
    FITREG_HFS[56] = 10;
    FITREG_HFS[57] = 10;
    FITREG_HFS[58] = 10;
   
    for (int ir = RANGE_HFS[0]; ir <= RANGE_HFS[1]; ++ir) {
        HistFit::Axis1D AX1D_HFS(
            "CC Estimator (#Lambda_{CC})",
            "Events/Bin",
            (*Hist::Head(Form("hHfsP%s_cc", subg.c_str())))()->GetYaxis(), FITREG_HFS[ir], 50);
        
        double rig = AXrig.center(ir, AxisScale::kLog);
        double acc_corr = haccp->Interpolate(rig);
        
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
      
        const int MAX_ii = 25+15;
        std::array<double, MAX_ii> hfitdv;
        Hist* hfitcc = Hist::New(Form("%sfitcc", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hHfsP%s_cc", subg.c_str())))()->GetYaxis())));
        Hist* hchicc = Hist::New(Form("%schicc", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hHfsP%s_cc", subg.c_str())))()->GetYaxis())));
        Hist* hflccc = Hist::New(Form("%sflccc", prefix.c_str()), "", HistAxis(Axis((*Hist::Head(Form("hHfsP%s_cc", subg.c_str())))()->GetYaxis())));
        for (int ii = 1; ii <= MAX_ii; ++ii) {
            HistFit::Axis1D AX1D_HFSi(
                "CC Estimator (#Lambda_{CC})",
                "Events/Bin",
                (*Hist::Head(Form("hHfsP%s_cc", subg.c_str())))()->GetYaxis(), ii, 50);
            HistFit::HistFit1D fit1Di(hlink_smp, hlink_tmps, AX1D_HFSi, prefix, false, false);
            HistFit::Hist1D    h1Drefi(Form("%sREF", prefix.c_str()), "", (TH1D*)((*vhHFS_pos.at(ir))()), AX1D_HFSi);
            (*hfitcc)()->SetBinContent(ii, acc_corr * fit1Di.wgts(0) / h1Drefi.data().sum());
            (*hfitcc)()->SetBinError  (ii, (*hfitcc)()->GetBinContent(ii) * fit1Di.errs(0) / fit1Di.wgts(0));
            (*hchicc)()->SetBinContent(ii, fit1Di.nchi());
            (*hflccc)()->SetBinContent(ii, 100.0 * (fit1Di.fluc(0).err / fit1Di.wgts(0)));
            hfitdv.at(ii-1) = (*hfitcc)()->GetBinContent(ii) / (acc_corr * (num_sig / num_ref)) - 1.0;
            hfitdv.at(ii-1) = hfitdv.at(ii-1) * hfitdv.at(ii-1);
        }
        (*hfitcc)()->GetXaxis()->SetTitle("CC Estimator (#Lambda_{CC})");
        (*hchicc)()->GetXaxis()->SetTitle("CC Estimator (#Lambda_{CC})");
        (*hflccc)()->GetXaxis()->SetTitle("CC Estimator (#Lambda_{CC})");
        (*hfitcc)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
        (*hchicc)()->GetYaxis()->SetTitle("#chi^{2}/NDF");
        (*hflccc)()->GetYaxis()->SetTitle("Relative Error (%)");
        (*hfitcc)()->GetXaxis()->SetRange(1, MAX_ii);
        (*hchicc)()->GetXaxis()->SetRange(1, MAX_ii);
        (*hflccc)()->GetXaxis()->SetRange(1, MAX_ii);
        double err_cuts = std::sqrt(std::accumulate(hfitdv.begin(), hfitdv.end(), 0.0) / hfitdv.size());

        double err_cfft = get_error_cutoff(rig);
        double err_tmpl = std::hypot(fit1D.fluc(0).err / num_sig, get_error_selection(rig));
        double err_ncnt = std::hypot(err_cfft, err_tmpl);
        double err_rscl = 0.01 * hrfnc->Interpolate(rig);
        double err_accp = 0.01 * hacce->Interpolate(rig);

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_rscl*err_rscl + err_accp*err_accp);
        double err_totl = std::hypot(err_stat, err_syst);

        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
       
        double app_cc_val = (fit1D.wgts(1) / num_ref);
        double app_cc_err = (fit1D.errs(1) / num_ref);

        std::cerr << Form("HFS STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
            fit1D.status(), ir, AXrig()(ir-1), AXrig()(ir),
            num_smp, num_sig, num_bkg,
            app_val*1.0e+4, err_stat*100.0, err_syst*100.0,
            nchi);

        (*hNchi)()->SetBinContent(ir, nchi); 
        
        (*hNcntPr)()->SetBinContent(ir, num_ref); 
        (*hNcntAp)()->SetBinContent(ir, num_sig); 
        
        (*hErrCfft)()->SetBinContent(ir, 100.0 * err_cfft);
        (*hErrTmpl)()->SetBinContent(ir, 100.0 * err_tmpl);
        (*hErrNcnt)()->SetBinContent(ir, 100.0 * err_ncnt); 
        (*hErrRfnc)()->SetBinContent(ir, 100.0 * err_rscl);
        (*hErrAccp)()->SetBinContent(ir, 100.0 * err_accp);
        
        (*hErrStat)()->SetBinContent(ir, 100.0 * err_stat); 
        (*hErrSyst)()->SetBinContent(ir, 100.0 * err_syst);
        (*hErrTotl)()->SetBinContent(ir, 100.0 * err_totl);
        
        (*hAccCorr)()->SetBinContent(ir, acc_corr);

        (*hAppStat)()->SetBinContent(ir, app_val); 
        (*hAppStat)()->SetBinError  (ir, app_err_stat); 
        
        (*hAppTotl)()->SetBinContent(ir, app_val);
        (*hAppTotl)()->SetBinError  (ir, app_err_totl);
        
        (*hAppCC)()->SetBinContent(ir, app_cc_val);
        (*hAppCC)()->SetBinError  (ir, app_cc_err);

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
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.35, 0.65, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "l");
        leg_table()->AddEntry((*hsig)(), Form("p^{-} (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "l");
        leg_table()->AddEntry((*hbkg)(), Form("p^{+} (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "l");
        leg_table()->AddEntry((TObject*)0, Form("#chi^{2}/NDF %.2f", nchi), "");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
        
        if (sw_fluc) {
            const HistFit::Hist1D& h1Dfluc = fit1D.fluc_hists(0);
            Hist* hfluc = Hist::New(h1Dfluc.get()->GetName(), 
                                    HistAxis(Axis("#bar{p}/p Flux Ratio [10^{-4}]", 
                                                  h1Dfluc.get()->GetXaxis()->GetNbins(), 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmin() / num_ref, 
                                                  1.0e+4 * acc_corr * h1Dfluc.get()->GetXaxis()->GetXmax() / num_ref), "Events/Bin"));
            for (int ib = 1; ib <= hfluc->xaxis().nbin(); ++ib) {
                (*hfluc)()->SetBinContent(ib, h1Dfluc.get()->GetBinContent(ib));
                (*hfluc)()->SetBinError  (ib, h1Dfluc.get()->GetBinError  (ib));
            }
            stdgaus->SetParameters(1.0, app_val, app_val * err_tmpl);
            stdgaus->SetLineColor(kBlue);
            stdgaus->SetNpx(100000);
            (*hfluc)()->Fit(stdgaus, "q0", "", hfluc->xaxis().min(), hfluc->xaxis().max());

            editor.create();
		    hfluc->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
            hfluc->draw("pe");
            stdgaus->Draw("l same");
            Legend leg_fluc("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.65, 0.85));
            leg_fluc()->SetHeader(Form("Rigidity %.2f - %.2f [GV]", AXrig()(ir-1), AXrig()(ir)));
            leg_fluc()->AddEntry((*hfluc)(), "Data", "lp");
            leg_fluc()->AddEntry(stdgaus, Form("Fit  #sigma_{tmpl.} = %.2f%", 100.0 * err_tmpl), "l");
            leg_fluc()->SetFillColor(0);
            leg_fluc.draw();
            editor.save();
        }
        
        editor.create();
		hfitcc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hfitcc->draw("pe");
        editor.save();
        
        editor.create();
		hchicc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hchicc->draw("hist");
        editor.save();
        
        editor.create();
		hflccc->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        hflccc->draw("hist");
        editor.save();
    }


    // Summary (Chi2)
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hNchi->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
    (*hNchi)()->GetXaxis()->SetMoreLogLabels();
    (*hNchi)()->Draw("hist");
    editor.save();

    // Summary (Ncnt)
    Hist* hcvsAppNcnt = Hist::New("hcvsAppNcnt", HistAxis( AXrig, Axis("Events/Bin", 1000, 1.0, 5.0 * (*hNcntPr)()->GetBinContent((*hNcntPr)()->GetMaximumBin()), AxisScale::kLog) ));
    (*hcvsAppNcnt)()->GetXaxis()->SetMoreLogLabels();
    (*hcvsAppNcnt)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hcvsAppNcnt)()->GetYaxis()->SetTitle("Events/Bin");
    
    hNcntPr->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kSquare )));
    hNcntAp->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle )));
    double totl_numpr = 0.0;
    double totl_numap = 0.0;
    for (int ib = 1; ib <= AXrig.nbin(); ++ib) totl_numpr += (*hNcntPr)()->GetBinContent(ib);
    for (int ib = 1; ib <= AXrig.nbin(); ++ib) totl_numap += (*hNcntAp)()->GetBinContent(ib);
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    hcvsAppNcnt->draw();
    (*hNcntPr)()->Draw("hist same");
    (*hNcntAp)()->Draw("hist same");
    (*hNcntPr)()->Draw("p same");
    (*hNcntAp)()->Draw("p same");
    Legend leg_ncnt("", TextStyle(kBlack, 20, 43), PadWindow(0.30, 0.60, 0.15, 0.35));
    leg_ncnt()->AddEntry((*hNcntPr)(), Form("N^{p} (Total %.2f B)", totl_numpr * 1e-9), "lp");
    leg_ncnt()->AddEntry((*hNcntAp)(), Form("N^{#bar{p}} (Total %.2f M)", totl_numap * 1e-6), "lp");
    leg_ncnt()->SetFillColor(0);
    leg_ncnt.draw();
    editor.save();

    // Sumary (Error)
    Hist* hcvsError = Hist::New("hcvsError", HistAxis( AXrig, Axis("Events/Bin", 1000, 0.0, 40.0) ));
    (*hcvsError)()->GetXaxis()->SetMoreLogLabels();
    (*hcvsError)()->GetXaxis()->CenterTitle();
    (*hcvsError)()->GetYaxis()->CenterTitle();
    (*hcvsError)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hcvsError)()->GetYaxis()->SetTitle("Relative Error (%)");
    
    hErrStat->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle )));
    hErrSyst->style(Line(kRed , 0, 2), Marker(kRed , MarkerStyle(MarkerShape::kCircle )));

    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hcvsError)()->Draw();
    (*hErrStat)()->Draw("hist same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hcvsError)()->Draw();
    (*hErrSyst)()->Draw("hist same");
    editor.save();
    
    hAppCC->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
    (*hAppCC)()->GetXaxis()->SetRange(46, 58);
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    (*hAppCC)()->Draw("pe");
    editor.save();
    
    THStack* hstack_err_ncnt = new THStack("hstack_err_ncnt", "hstack_err_ncnt");
    hErrNcnt->style(Line(kRed    , 0, 4), Marker(kRed     , MarkerStyle(MarkerShape::kCircle)));
    hErrTmpl->style(Line(kBlue   , 0, 2), Marker(kBlue    , MarkerStyle(MarkerShape::kCircle)));
    hErrCfft->style(Line(kGreen+2, 0, 2), Marker(kGreen+2 , MarkerStyle(MarkerShape::kCircle)));
    hstack_err_ncnt->Add((*hErrNcnt)());
    hstack_err_ncnt->Add((*hErrCfft)());
    hstack_err_ncnt->Add((*hErrTmpl)());
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hstack_err_ncnt->Draw("nostack hist");
    hstack_err_ncnt->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    hstack_err_ncnt->GetHistogram()->GetXaxis()->CenterTitle();
    hstack_err_ncnt->GetHistogram()->GetYaxis()->CenterTitle();
    hstack_err_ncnt->GetHistogram()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    hstack_err_ncnt->GetHistogram()->GetYaxis()->SetTitle("Relative Error (%)");
    hstack_err_ncnt->Draw("nostack hist");
    Legend leg_err_ncnt("", TextStyle(kBlack, 20, 43), PadWindow(0.35, 0.65, 0.65, 0.85));
    leg_err_ncnt()->AddEntry((*hErrNcnt)(), "Total Systematic Errors", "l");
    leg_err_ncnt()->AddEntry((*hErrTmpl)(), "Shape of Templates", "l");
    leg_err_ncnt()->AddEntry((*hErrCfft)(), "Geomagnetic Rigidity Cutoff", "l");
    leg_err_ncnt()->SetFillColor(0);
    leg_err_ncnt.draw();
    editor.save();
   
    THStack* hstack_err_syst = new THStack("hstack_err_syst", "hstack_err_syst");
    hErrSyst->style(Line(kRed     , 0, 4), Marker(kRed     , MarkerStyle(MarkerShape::kCircle)));
    hErrNcnt->style(Line(kBlue    , 0, 2), Marker(kBlue    , MarkerStyle(MarkerShape::kCircle)));
    hErrAccp->style(Line(kGreen+2 , 0, 2), Marker(kGreen+2 , MarkerStyle(MarkerShape::kCircle)));
    hErrRfnc->style(Line(kYellow+1, 0, 2), Marker(kYellow+1, MarkerStyle(MarkerShape::kCircle)));
    hstack_err_syst->Add((*hErrSyst)());
    hstack_err_syst->Add((*hErrRfnc)());
    hstack_err_syst->Add((*hErrAccp)());
    hstack_err_syst->Add((*hErrNcnt)());
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    //(*hcvsError)()->Draw();
    hstack_err_syst->Draw("nostack hist");
    hstack_err_syst->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    hstack_err_syst->GetHistogram()->GetXaxis()->CenterTitle();
    hstack_err_syst->GetHistogram()->GetYaxis()->CenterTitle();
    hstack_err_syst->GetHistogram()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    hstack_err_syst->GetHistogram()->GetYaxis()->SetTitle("Relative Error (%)");
    hstack_err_syst->Draw("nostack hist");
    Legend leg_err_syst("", TextStyle(kBlack, 20, 43), PadWindow(0.35, 0.65, 0.65, 0.85));
    leg_err_syst()->AddEntry((*hErrSyst)(), "Total Systematic Errors", "l");
    leg_err_syst()->AddEntry((*hErrNcnt)(), "Event Selection", "l");
    leg_err_syst()->AddEntry((*hErrAccp)(), "Effective Acceptance", "l");
    leg_err_syst()->AddEntry((*hErrRfnc)(), "Rigidity Measurement", "l");
    leg_err_syst()->SetFillColor(0);
    leg_err_syst.draw();
    editor.save();
    
    THStack* hstack_err_totl = new THStack("hstack_err_totl", "hstack_err_totl");
    hErrTotl->style(Line(kRed    , 0, 4), Marker(kRed    , MarkerStyle(MarkerShape::kCircle)));
    hErrStat->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle)));
    hErrSyst->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));
    hstack_err_totl->Add((*hErrTotl)());
    hstack_err_totl->Add((*hErrSyst)());
    hstack_err_totl->Add((*hErrStat)());
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    //(*hcvsError)()->Draw();
    hstack_err_totl->Draw("nostack hist");
    hstack_err_totl->GetHistogram()->GetXaxis()->SetMoreLogLabels();
    hstack_err_totl->GetHistogram()->GetXaxis()->CenterTitle();
    hstack_err_totl->GetHistogram()->GetYaxis()->CenterTitle();
    hstack_err_totl->GetHistogram()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    hstack_err_totl->GetHistogram()->GetYaxis()->SetTitle("Relative Error (%)");
    hstack_err_totl->Draw("nostack hist");
    Legend leg_err_totl("", TextStyle(kBlack, 20, 43), PadWindow(0.35, 0.65, 0.65, 0.85));
    leg_err_totl()->AddEntry((*hErrTotl)(), "Total Errors", "l");
    leg_err_totl()->AddEntry((*hErrStat)(), "Statistical Errors", "l");
    leg_err_totl()->AddEntry((*hErrSyst)(), "Systematic Errors", "l");
    leg_err_totl.draw();
    editor.save();
    
    Hist* hcvsAppFR = Hist::New("hcvsAppFR", HistAxis(AXrig, Axis("#bar{p}/p Flux Ratio", 2000, 3.0e-6, 3.0e-4)));
    (*hcvsAppFR)()->GetXaxis()->SetMoreLogLabels();
    (*hcvsAppFR)()->GetXaxis()->CenterTitle();
    (*hcvsAppFR)()->GetYaxis()->CenterTitle();
    (*hcvsAppFR)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hcvsAppFR)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
	
    hAppStat->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
    hAppTotl->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
        
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hcvsAppFR->draw();
    (*hAppStat)()->Draw("pe same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    hcvsAppFR->draw();
    (*hAppStat)()->Draw("pe same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hcvsAppFR->draw();
    (*hAppTotl)()->Draw("pe same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    hcvsAppFR->draw();
    (*hAppTotl)()->Draw("pe same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hcvsAppFR->draw();
    offapp->Draw("pe same");
    (*hAppStat)()->Draw("pe same");
    Legend leg_appr_stat("", TextStyle(kBlack, 40, 43), PadWindow(0.60, 0.85, 0.20, 0.35));
    leg_appr_stat()->AddEntry(offapp, "PRL", "lp");
    leg_appr_stat()->AddEntry((*hAppStat)(), "Result", "lp");
    leg_appr_stat()->SetFillColor(0);
    leg_appr_stat.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    hcvsAppFR->draw();
    offapp->Draw("pe same");
    (*hAppStat)()->Draw("pe same");
    leg_appr_stat.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    hcvsAppFR->draw();
    offapp->Draw("pe same");
    (*hAppStat)()->Draw("pe same");
    leg_appr_stat.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hcvsAppFR->draw();
    offapp->Draw("pe same");
    (*hAppTotl)()->Draw("pe same");
    Legend leg_appr_totl("", TextStyle(kBlack, 40, 43), PadWindow(0.60, 0.85, 0.20, 0.35));
    leg_appr_totl()->AddEntry(offapp, "PRL", "lp");
    leg_appr_totl()->AddEntry((*hAppTotl)(), "Result", "lp");
    leg_appr_totl()->SetFillColor(0);
    leg_appr_totl.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    hcvsAppFR->draw();
    offapp->Draw("pe same");
    (*hAppTotl)()->Draw("pe same");
    leg_appr_totl.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    hcvsAppFR->draw();
    offapp->Draw("pe same");
    (*hAppTotl)()->Draw("pe same");
    leg_appr_totl.draw();
    editor.save();

    editor.close();

    TFile * ofle = new TFile(Form("out/apflux_flx%s.root", subg.c_str()), "RECREATE");
    ofle->cd();

    (*hNchi)()->Write();
    
    (*hNcntPr)()->Write();
    (*hNcntAp)()->Write();
    
    (*hErrCfft)()->Write();
    (*hErrTmpl)()->Write();
    (*hErrNcnt)()->Write();
    (*hErrRfnc)()->Write();
    (*hErrAccp)()->Write();
    
    (*hErrStat)()->Write();
    (*hErrSyst)()->Write();
    (*hErrTotl)()->Write();
   
    hstack_err_ncnt->Write();
    hstack_err_syst->Write();
    hstack_err_totl->Write();

    (*hAccCorr)()->Write();
    (*hAppStat)()->Write();
    (*hAppTotl)()->Write();
    
    (*hAppCC)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
