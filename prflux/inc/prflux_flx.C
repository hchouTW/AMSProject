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
    
    TGraphErrors* offapp = (TGraphErrors*) (TFile::Open("others/database_prflux.root")->Get("gr_exp1"));

    TFile* file_acc = TFile::Open("out/prflux_acc.root");
    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/iss%s", subv.c_str()));

    const Axis& AXrig = Hist::Head("hFlxP_cnt")->xaxis();
    
    Hist* hNchi = Hist::New("hNchi", HistAxis(AXrig, "#chi^{2}"));
    
    Hist* hNcntPr = Hist::New("hNcntPr", HistAxis(AXrig, "Number"));
                                    
    Hist* hErrNcnt = Hist::New("hErrNcnt", HistAxis(AXrig, "RelError of Count (%)"));
    Hist* hErrRscl = Hist::New("hErrRscl", HistAxis(AXrig, "RelError of Rigidity Scale (%)"));
    Hist* hErrAccp = Hist::New("hErrAccp", HistAxis(AXrig, "RelError of Acceptance (%)"));
    
    Hist* hErrStat = Hist::New("hErrStat", HistAxis(AXrig, "RelError of Stat (%)"));
    Hist* hErrSyst = Hist::New("hErrSyst", HistAxis(AXrig, "RelError of Syst (%)"));
    Hist* hErrTotl = Hist::New("hErrTotl", HistAxis(AXrig, "RelError of Totl (%)"));
    
    Hist* hAccCorr = Hist::New("hAccCorr", HistAxis(AXrig, "Acceptance Correction"));
    Hist* hAppStat = Hist::New("hAppStat", HistAxis(AXrig, "p"));
    Hist* hAppTotl = Hist::New("hAppTotl", HistAxis(AXrig, "p"));
   
    PdfEditor editor(Window(), "prflux_flx", "out");

    bool sw_fluc = false; 
    bool sw_build_hist = true; 

    /////////////////
    ////         ////
    ////   FLX   ////
    ////         ////
    /////////////////
    TH1D* hFLX_dR  = (TH1D*) file_acc->Get("hFlx_dR");
    TH1D* hFLX_acc = (TH1D*) file_acc->Get("hFlx_acc");
    TH1D* hFLX_lv  = (TH1D*) file_acc->Get("hFlx_lv");
    TH1D* hFLX_trg = (TH1D*) file_acc->Get("hFlx_trg");
   
    Hist* hFLX_cnt = Hist::Head("hFlxP_cnt");
    std::array<int, 2> RANGE_FLX({ 1, 58 });
    for (int ir = RANGE_FLX[0]; ir <= RANGE_FLX[1]; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);

        double ncnt_value = (*hFLX_cnt)()->GetBinContent(ir);
        double ncnt_rlerr = (*hFLX_cnt)()->GetBinError(ir) / (*hFLX_cnt)()->GetBinContent(ir);

        double dR  = hFLX_dR->GetBinContent(ir);
        double lv  = hFLX_lv->GetBinContent(ir);
        double trg = hFLX_trg->GetBinContent(ir);
        double acc = hFLX_acc->GetBinContent(ir);

        double rat_val = ncnt_value / (dR * lv);
        double rat_err = rat_val * ncnt_rlerr;

        double flx_val = rat_val / (trg * acc);
        double flx_err = rat_err / (trg * acc);

        std::cerr << Form("FLX R%3d (%.1f %.1f) SIG %8.1f PR %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f)\n",
            ir, AXrig()(ir-1), AXrig()(ir),
            ncnt_value, 
            flx_val*1.0e-4*std::pow(rig, 2.7), flx_err/flx_val*100.0, 0.0);



/*
        double err_ncnt = counting_error_func(AXrig.center(ir, AxisScale::kLog));
        double err_rscl = rigidity_scale_func(AXrig()(ir-1), AXrig()(ir));
        double err_accp = acceptance_error_func(AXrig.center(ir, AxisScale::kLog));

        double err_stat = fit1D.errs(0) / num_sig;
        double err_syst = std::sqrt(err_ncnt*err_ncnt + err_tmpl*err_tmpl + err_rscl*err_rscl + err_accp*err_accp + err_cc*err_cc);
        double err_totl = std::hypot(err_stat, err_syst);

        double acc_corr     = hFLX_acc->Interpolate(rig);
        double app_val      = acc_corr * (num_sig / num_ref);
        double app_err_stat = app_val * err_stat;
        double app_err_syst = app_val * err_syst;
        double app_err_totl = app_val * err_totl;
        
        std::cerr << Form("FLX STATUS %d R%3d (%.1f %.1f) SMP %8.1f (SIG %8.1f BKG %8.1f) APP %10.4f(e-4) RELERR(STAT %10.4f SYST %10.4f) NCHI %10.3f\n",
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
*/
        (*hAppStat)()->SetBinContent(ir, flx_val); 
        (*hAppStat)()->SetBinError  (ir, flx_err); 
        
        (*hAppTotl)()->SetBinContent(ir, flx_val * std::pow(rig, 2.7));
        (*hAppTotl)()->SetBinError  (ir, flx_err * std::pow(rig, 2.7));
    }

    TFile * ofle = new TFile("out/prflux_flx.root", "RECREATE");
    ofle->cd();

    (*hNchi)()->Write();
    
    (*hNcntPr)()->Write();
    
    (*hErrNcnt)()->Write();
    (*hErrRscl)()->Write();
    (*hErrAccp)()->Write();
    
    (*hErrStat)()->Write();
    (*hErrSyst)()->Write();
    (*hErrTotl)()->Write();
    
    (*hAccCorr)()->Write();
    (*hAppStat)()->Write();
    (*hAppTotl)()->Write();

    offapp->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
