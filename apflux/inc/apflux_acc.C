#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "HistFit1D.h"
#include "HistFit1D.C"

inline double eval_rigrso(double arig, int opt = 0) {
    double pars[3] = { 9.835e-2, 9.964e-3, 4.2183e-4 };
    //double rigrso = std::sqrt(pars[0] * pars[0] + pars[1] * pars[1] * arig + pars[2] * pars[2] * arig * arig) / (1.0 + 18.0 * std::exp(-2.0 * arig));
    double rigrso = std::sqrt(pars[0] * pars[0] + pars[1] * pars[1] * arig + pars[2] * pars[2] * arig * arig);
    if (opt) rigrso /= (1.0 + 18.0 * std::exp(-2.0 * arig));
    return rigrso;
}
inline double eval_prflux(double arig, double powidx = 0.0) { 
    double prflux = 1.27152e+04 * std::pow((arig / 45.0), -2.99752e-01 + powidx)
                                * std::pow(1.0 + std::pow((arig / 2.28781e+00), -1.56940e+00), -2.20494e+00)
                                * std::pow(1.0 + std::pow((arig / 3.63316e+02),  5.56393e+00),  2.29502e-02);
    return prflux;
}
inline double eval_apflux(double arig, double powidx = 0.0) { 
    double prflux = eval_prflux(arig, powidx);
    double appr   = 2.04533e-04 * std::pow((arig / 45.0), -1.60252e-01)
                                * std::pow(1.0 + std::pow((arig / 4.01379e+00), -1.49877e+00), -1.93532e+00);
    double apflux = prflux * appr;
    return apflux; 
}
inline double cross_section_fluc_func(double arig) {
    double fluc_crss = 1.0 + 1.33 * std::pow(arig, -1.23130e-01) * std::exp(-2.91143e-02 * arig);
    if (!std::isfinite(fluc_crss)) fluc_crss = 1.0;
    return fluc_crss;
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "08";

    constexpr int NSET = 6;
    std::array<std::string, NSET> name({ "hLtf", "hLrh", "hIin", "hIex", "hHex", "hHfs"});
    std::array<std::string, NSET> name_of_hcnt({ "hLtf_cnt_MC", "hLrh_cnt_MC", "hIin_cnt_MC", "hIex_cnt_MC", "hHex_cnt_MC", "hHfs_cnt_MC"});
    std::array<std::string, NSET> name_of_accp({ "hLtf_accp", "hLrh_accp", "hIin_accp", "hIex_accp", "hHex_accp", "hHfs_accp" });
    std::array<std::string, NSET> name_of_crss({ "hLtf_crss", "hLrh_crss", "hIin_crss", "hIex_crss", "hHex_crss", "hHfs_crss" });
    
    std::vector<double> vrig( {
         1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
         3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
         8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
        19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
        41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
        93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 330.00, 525.00 } );
    std::vector<std::array<double, 2>> range({ { 1.00, 3.64 }, { 3.64, 13.0 }, { 13.0, 31.1 }, { 31.1, 60.3 }, { 60.3, 211.0 }, { 211.0, 525.0 } });
    std::vector<std::array<int, 2>> range_bins({ { 1, 11 }, { 12, 25 }, { 26, 36 }, { 37, 45 }, { 46, 55 }, { 56, 58 } });
    Axis AXrig("|Rigidity| [GV]", vrig);

    std::array<TH1D*, NSET> hist_accp_raw;
    std::array<TH1D*, NSET> hist_crss_raw;
    TH1D* hist_accp_raw_data = new TH1D("haccp_raw", "haccp_raw;|Rigidity| [GV/c];Acceptance Correction", vrig.size()-1, vrig.data());
    TH1D* hist_crss_raw_data = new TH1D("hcrss_raw", "hcrss_raw;|Rigidity| [GV/c];Cross Section"        , vrig.size()-1, vrig.data());

    for (int is = 0; is < NSET; ++is) {
        std::cerr << Form("SET %2d:\n", is);
        UInt_t cntev = 0;
    
        // Acceptance Correction
        UInt_t cntpr = 0;
        TFile* fmcpr = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr%s/YiMdst.root", subv.c_str()));
        TH1D*  hmcpr = (TH1D*)fmcpr->Get(name_of_hcnt[is].c_str());
        TTree* tmcpr = (TTree*)fmcpr->Get("runlist");
        tmcpr->SetBranchAddress("event", &cntev);
        for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
        
        UInt_t cntap = 0;
        TFile* fmcap = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap%s/YiMdst.root", subv.c_str()));
        TH1D*  hmcap = (TH1D*)fmcap->Get(name_of_hcnt[is].c_str());
        TTree* tmcap = (TTree*)fmcap->Get("runlist");
        tmcap->SetBranchAddress("event", &cntev);
        for (int it = 0; it < tmcap->GetEntries(); ++it) { tmcap->GetEntry(it); cntap+=cntev; }

        TH1D* haccp_raw = new TH1D(Form("%s_raw", name_of_accp[is].c_str()), "", hmcpr->GetXaxis()->GetNbins(), hmcpr->GetXaxis()->GetXbins()->GetArray());
        for (int ib = 1; ib <= haccp_raw->GetXaxis()->GetNbins(); ++ib) {
            double errpr = hmcpr->GetBinError(ib) / hmcpr->GetBinContent(ib);
            double errap = hmcap->GetBinError(ib) / hmcap->GetBinContent(ib);
            double error = std::sqrt(errpr * errpr + errap * errap);
            double accp = (hmcap->GetBinContent(ib) / static_cast<double>(cntap)) / (hmcpr->GetBinContent(ib) / static_cast<double>(cntpr));
            if (!std::isfinite(accp) || accp <= 0.0) continue;
            haccp_raw->SetBinContent(ib, 1.0/accp);
            haccp_raw->SetBinError  (ib, 1.0/accp * error);

            if (ib < range_bins.at(is)[0] || ib > range_bins.at(is)[1]) continue;
            hist_accp_raw_data->SetBinContent(ib, haccp_raw->GetBinContent(ib));
            hist_accp_raw_data->SetBinError  (ib, haccp_raw->GetBinError  (ib));
        }
        haccp_raw->GetXaxis()->SetTitle(hist_accp_raw_data->GetXaxis()->GetTitle());
        haccp_raw->GetYaxis()->SetTitle(hist_accp_raw_data->GetYaxis()->GetTitle());
        hist_accp_raw[is] = haccp_raw;
       
        // Cross Section
        UInt_t cntapp = 0;
        TFile* fmcapp = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap_plus%s/YiMdst.root", subv.c_str()));
        TH1D*  hmcapp = (TH1D*)fmcapp->Get(name_of_hcnt[is].c_str());
        TTree* tmcapp = (TTree*)fmcapp->Get("runlist");
        tmcapp->SetBranchAddress("event", &cntev);
        for (int it = 0; it < tmcapp->GetEntries(); ++it) { tmcapp->GetEntry(it); cntapp+=cntev; }
        
        UInt_t cntapm = 0;
        TFile* fmcapm = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap_minus%s/YiMdst.root", subv.c_str()));
        TH1D*  hmcapm = (TH1D*)fmcapm->Get(name_of_hcnt[is].c_str());
        TTree* tmcapm = (TTree*)fmcapm->Get("runlist");
        tmcapm->SetBranchAddress("event", &cntev);
        for (int it = 0; it < tmcapm->GetEntries(); ++it) { tmcapm->GetEntry(it); cntapm+=cntev; }
        
        TH1D* hcrss_raw = new TH1D(Form("%s_raw", name_of_crss[is].c_str()), "", hmcapp->GetXaxis()->GetNbins(), hmcapp->GetXaxis()->GetXbins()->GetArray());
        for (int ib = 1; ib <= hcrss_raw->GetXaxis()->GetNbins(); ++ib) {
            double errapp = hmcapp->GetBinError(ib) / hmcapp->GetBinContent(ib);
            double errapm = hmcapm->GetBinError(ib) / hmcapm->GetBinContent(ib);
            double error = std::sqrt(errapp * errapp + errapm * errapm);
            double cross = (hmcapp->GetBinContent(ib) / static_cast<double>(cntapp)) / (hmcapm->GetBinContent(ib) / static_cast<double>(cntapm));
            if (!std::isfinite(cross) || cross <= 0.0) continue;
            hcrss_raw->SetBinContent(ib, cross - 1.0);
            hcrss_raw->SetBinError  (ib, cross * error);
            
            if (ib < range_bins.at(is)[0] || ib > range_bins.at(is)[1]) continue;
            hist_crss_raw_data->SetBinContent(ib, hcrss_raw->GetBinContent(ib));
            hist_crss_raw_data->SetBinError  (ib, hcrss_raw->GetBinError  (ib));
        }
        hcrss_raw->GetXaxis()->SetTitle(hist_crss_raw_data->GetXaxis()->GetTitle());
        hcrss_raw->GetYaxis()->SetTitle(hist_crss_raw_data->GetYaxis()->GetTitle());
        hist_crss_raw[is] = hcrss_raw;
    }
    
    Axis AXrig_smooth("|Rigidity| [GV]", 10000, 1.0, 525.0, AxisScale::kLog);
   
    //**** Acceptance Error ****//
    TH1D* hist_acce_raw_data = new TH1D("hacce_raw", "hacce_raw;|Rigidity| [GV/c];Acceptance Error", AXrig.nbin(), &(AXrig(0)));
    TH1D* hist_accp_smooth = new TH1D("haccp_smooth", "haccp;|Rigidity| [GV/c];Acceptance Correction Factor", AXrig_smooth.nbin(), &(AXrig_smooth(0)));
    TH1D* hist_acce_smooth = new TH1D("hacce_smooth", "hacce;|Rigidity| [GV/c];Error of Acceptance Correction Factor", AXrig_smooth.nbin(), &(AXrig_smooth(0)));
    {
        TF1* accp_func = new TF1("accp_func", "[0] + [1] * pow(x, -[2]) * (1.0 + [3] * erf([4] * (log(1.0 + x) - [5]))) * (1.0 + [6] * erf([7] * (log(1.0 + x) - [8]))) * (1.0 - [9] * ((x-[10])/[11]) / exp(abs(x-[10])/[11]))", 0.5, 10000);
        accp_func->SetParameters(1.06, 0.5, 1.2, 0.16, 8.0, 1.1); 
        accp_func->SetParameters(1.06, 0.5, 1.2, 0.16, 8.0, 1.1, 0.1, 8.0, 4.5); 
        accp_func->SetParameter(0, 1.06);
        accp_func->SetParameter(1, 0.5);
        accp_func->SetParameter(2, 1.2);
        accp_func->SetParameter(3, 0.16); 
        accp_func->SetParameter(4, 8.0);
        accp_func->SetParameter(5, 1.1);
        accp_func->SetParameter(6, 0.1);
        accp_func->SetParameter(7, 8.0);
        accp_func->SetParameter(8, 4.5); 
        accp_func->SetParameter(9, 0.2); 
        accp_func->SetParameter(10, 3.5); 
        accp_func->SetParameter(11, 1.0); 
        TH1D* haccp_raw_data_clone = (TH1D*) hist_accp_raw_data->Clone();
        haccp_raw_data_clone->Fit(accp_func, "q0", "", vrig.front(), vrig.back());
        haccp_raw_data_clone->Fit(accp_func, "q0", "", vrig.front(), vrig.back());
        haccp_raw_data_clone->Fit(accp_func, "q0", "", vrig.front(), vrig.back());
        haccp_raw_data_clone->Fit(accp_func, "q0", "", vrig.front(), vrig.back());
        haccp_raw_data_clone->Fit(accp_func, "q0", "", vrig.front(), vrig.back());

        std::vector<double> vd2;
        for (int ib = 1; ib <= hist_acce_raw_data->GetXaxis()->GetNbins(); ++ib) {
            hist_acce_raw_data->SetBinContent(ib, hist_accp_raw_data->GetBinContent(ib) / accp_func->Eval(hist_accp_raw_data->GetXaxis()->GetBinCenter(ib)) - 1.0);
            hist_acce_raw_data->SetBinError  (ib, hist_accp_raw_data->GetBinError  (ib) / accp_func->Eval(hist_accp_raw_data->GetXaxis()->GetBinCenter(ib)));
            vd2.push_back(hist_acce_raw_data->GetBinContent(ib) * hist_acce_raw_data->GetBinContent(ib));
        }
        double derr = std::sqrt(std::accumulate(vd2.begin(), vd2.end(), 0.0) / vd2.size());
        
        for (int ir = 1; ir <= AXrig_smooth.nbin(); ++ir) {
            double rig = AXrig_smooth.center(ir, AxisScale::kLog);
            hist_accp_smooth->SetBinContent(ir, accp_func->Eval(rig));
            hist_accp_smooth->SetBinError  (ir, accp_func->Eval(rig) * derr);
            hist_acce_smooth->SetBinContent(ir, derr);
        }
    }
   
    //**** Cross Section ****//
    TH1D* hist_crss_smooth = new TH1D("hcrss", "hcrss;|Rigidity| [GV/c];Cross Section", AXrig_smooth.nbin(), &(AXrig_smooth(0)));
    { 
        std::vector<double> vd2;
        for (int ib = 1; ib <= hist_crss_raw_data->GetXaxis()->GetNbins(); ++ib) {
            vd2.push_back(hist_crss_raw_data->GetBinContent(ib) * hist_crss_raw_data->GetBinContent(ib));
        }
        double derr = std::sqrt(std::accumulate(vd2.begin(), vd2.end(), 0.0) / vd2.size());
        
        for (int ir = 1; ir <= AXrig_smooth.nbin(); ++ir) {
            double err_crss = cross_section_fluc_func(AXrig_smooth.center(ir, AxisScale::kLog));
            err_crss *= 2.0; // testcode
            hist_crss_smooth->SetBinContent(ir, derr * err_crss);
        }
    }
    
    Hist* hist_summary_accptance_factor_smooth = Hist::New("hist_accptance_factor_smooth", HistAxis(AXrig_smooth));
    Hist* hist_summary_accptance_factor_raw    = Hist::New("hist_accptance_factor_raw"   , HistAxis(AXrig));
    Hist* hist_summary_accptance_factor        = Hist::New("hist_accptance_factor"       , HistAxis(AXrig));
    Hist* hist_summary_accptance_error_accp    = Hist::New("hist_accptance_error_accp"   , HistAxis(AXrig));
    Hist* hist_summary_accptance_error_crss    = Hist::New("hist_accptance_error_crss"   , HistAxis(AXrig));
    Hist* hist_summary_accptance_error         = Hist::New("hist_accptance_error"        , HistAxis(AXrig));
    for (int ir = 1; ir <= AXrig_smooth.nbin(); ++ir) {
        (*hist_summary_accptance_factor_smooth)()->SetBinContent(ir,  hist_accp_smooth->GetBinContent(ir));
        (*hist_summary_accptance_factor_smooth)()->SetBinError  (ir,  hist_accp_smooth->GetBinError  (ir));
    }
    for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        (*hist_summary_accptance_factor_raw)()->SetBinContent(ir, hist_accp_raw_data->GetBinContent(ir));
        (*hist_summary_accptance_factor_raw)()->SetBinError  (ir, hist_accp_raw_data->GetBinError  (ir));
        (*hist_summary_accptance_factor)()->SetBinContent(ir, hist_accp_smooth->Interpolate(rig));
        //(*hist_summary_accptance_factor)()->SetBinError  (ir, hist_accp_smooth->Interpolate(rig) * hist_acce_smooth->Interpolate(rig));
        (*hist_summary_accptance_error_accp)()->SetBinContent(ir, 100.0 * hist_acce_smooth->Interpolate(rig));
        (*hist_summary_accptance_error_crss)()->SetBinContent(ir, 100.0 * hist_crss_smooth->Interpolate(rig));

        double total_error = std::hypot((*hist_summary_accptance_error_accp)()->GetBinContent(ir), (*hist_summary_accptance_error_crss)()->GetBinContent(ir));
        (*hist_summary_accptance_error)()->SetBinContent(ir, total_error);
    }


    //**** Rigidity Scale and Migration Effect ****//
    TH1D* hist_mgrt = new TH1D("hmgrt", "hmgrt;|Rigidity| [GV];Migration Effect", AXrig.nbin(), &(AXrig(0)));
    TH1D* hist_mflc = new TH1D("hmflc", "hmflc;|Rigidity| [GV];Migration Effect", AXrig.nbin(), &(AXrig(0)));
    TH1D* hist_mcrr = new TH1D("hmcrr", "hmcrr;|Rigidity| [GV];Migration Effect", AXrig.nbin(), &(AXrig(0)));
    TH1D* hist_mcrrORG = new TH1D("hmcrrORG", "hmcrrORG;|Rigidity| [GV];Migration Effect", AXrig.nbin(), &(AXrig(0)));
    {
        Axis AXmom("", 10000, 0.3, 60000.0, AxisScale::kLog);
        TH1D* hpr_raw = new TH1D("hmgrt_pr_raw", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hap_raw = new TH1D("hmgrt_ap_raw", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hpr     = new TH1D("hmgrt_pr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hap     = new TH1D("hmgrt_ap", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hprORG = new TH1D("hmgrt_prORG", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hapORG = new TH1D("hmgrt_apORG", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hp10pr = new TH1D("hmgrt_p10pr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hp10ap = new TH1D("hmgrt_p10ap", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hm10pr = new TH1D("hmgrt_m10pr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hm10ap = new TH1D("hmgrt_m10ap", "", AXrig.nbin(), &(AXrig(0)));
        for (int im = 1; im <= AXmom.nbin(); ++im) {
            double mom = AXmom.center(im, AxisScale::kLog);
            double dm  = AXmom.width(im);
            double rsoORG = eval_rigrso(mom, 0);
            double rso = eval_rigrso(mom, 0);
            double pr  = eval_prflux(mom) * dm;
            double ap  = eval_apflux(mom) * dm;
            hpr_raw->SetBinContent(hpr_raw->FindBin(mom), hpr_raw->GetBinContent(hpr_raw->FindBin(mom)) + pr);
            hap_raw->SetBinContent(hap_raw->FindBin(mom), hap_raw->GetBinContent(hap_raw->FindBin(mom)) + ap);
            for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
                double lbdORG = (mom/AXrig(ir)   - 1.0) / (rsoORG);
                double ubdORG = (mom/AXrig(ir-1) - 1.0) / (rsoORG);
                double lbd = (mom/AXrig(ir)   - 1.0) / (rso);
                double ubd = (mom/AXrig(ir-1) - 1.0) / (rso);
                double lbdp10 = (mom/AXrig(ir)   - 1.0) / (0.9 * rso);
                double ubdp10 = (mom/AXrig(ir-1) - 1.0) / (0.9 * rso);
                double lbdm10 = (mom/AXrig(ir)   - 1.0) / (1.1 * rso);
                double ubdm10 = (mom/AXrig(ir-1) - 1.0) / (1.1 * rso);
                double wORG = 0.5 * (std::erf(ubdORG) - std::erf(lbdORG));
                double w    = 0.5 * (std::erf(ubd) - std::erf(lbd));
                double wp10 = 0.5 * (std::erf(ubdp10) - std::erf(lbdp10));
                double wm10 = 0.5 * (std::erf(ubdm10) - std::erf(lbdm10));
                hprORG->SetBinContent(ir, hprORG->GetBinContent(ir) + pr * wORG);
                hapORG->SetBinContent(ir, hapORG->GetBinContent(ir) + ap * wORG);
                hpr->SetBinContent(ir, hpr->GetBinContent(ir) + pr * w);
                hap->SetBinContent(ir, hap->GetBinContent(ir) + ap * w);
                hp10pr->SetBinContent(ir, hp10pr->GetBinContent(ir) + pr * wp10);
                hp10ap->SetBinContent(ir, hp10ap->GetBinContent(ir) + ap * wp10);
                hm10pr->SetBinContent(ir, hm10pr->GetBinContent(ir) + pr * wm10);
                hm10ap->SetBinContent(ir, hm10ap->GetBinContent(ir) + ap * wm10);
            }
        }
        TH1D* hmcrr = new TH1D("hmcrr_appr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hmcrrORG = new TH1D("hmcrrORG_appr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hp10appr = new TH1D("hmgrt_p10appr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hm10appr = new TH1D("hmgrt_m10appr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hpm10    = new TH1D("hmgrt_pm10"   , "", AXrig.nbin(), &(AXrig(0)));
        for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
            double p10appr = hp10ap->GetBinContent(ir) / hp10pr->GetBinContent(ir);
            double m10appr = hm10ap->GetBinContent(ir) / hm10pr->GetBinContent(ir);
            double diff    = (p10appr - m10appr) / (p10appr + m10appr);
            hp10appr->SetBinContent(ir, p10appr);
            hm10appr->SetBinContent(ir, m10appr);
            hpm10->SetBinContent(ir, diff);
            
            double prrat = hpr->GetBinContent(ir) / hpr_raw->GetBinContent(ir);
            double aprat = hap->GetBinContent(ir) / hap_raw->GetBinContent(ir);
            double mcrr = (prrat / aprat - 1.0);
            hmcrr->SetBinContent(ir, mcrr);
            
            double prratORG = hprORG->GetBinContent(ir) / hpr_raw->GetBinContent(ir);
            double apratORG = hapORG->GetBinContent(ir) / hap_raw->GetBinContent(ir);
            double mcrrORG = (prratORG / apratORG - 1.0);
            hmcrrORG->SetBinContent(ir, mcrrORG);
        }
        for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
            hist_mflc->SetBinContent(ir, hpm10->GetBinContent(ir));
            hist_mcrr->SetBinContent(ir, hmcrr->GetBinContent(ir));
            hist_mcrrORG->SetBinContent(ir, hmcrrORG->GetBinContent(ir));
           
            double mgrt = std::hypot(hist_mcrrORG->GetBinContent(ir) / 4.0, 0.0);
            hist_mgrt->SetBinContent(ir, mgrt);
        }
    }

    TH1D* hist_rscl = new TH1D("hrscl", "hrscl;|Rigidity| [GV/c];Rigidity Scale", AXrig.nbin(), &(AXrig(0)));
    {
        const double shift = 1.0 / 26000.0; // (1/26 TV)
        TH1D* hspappr = new TH1D("rscl_hspappr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hsmappr = new TH1D("rscl_hsmappr", "", AXrig.nbin(), &(AXrig(0)));
        TH1D* hpsm    = new TH1D("rscl_hpsm"   , "", AXrig.nbin(), &(AXrig(0)));
        for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
            double lbdsp = 1.0 / (1.0/AXrig(ir-1) + 1.0/26000.0);
            double ubdsp = 1.0 / (1.0/AXrig(ir)   + 1.0/26000.0);
            double lbdsm = 1.0 / (1.0/AXrig(ir-1) - 1.0/26000.0);
            double ubdsm = 1.0 / (1.0/AXrig(ir)   - 1.0/26000.0);
            double spappr = std::sqrt(eval_apflux(lbdsm) * eval_apflux(ubdsm)) / std::sqrt(eval_prflux(lbdsp) * eval_prflux(ubdsp));
            double smappr = std::sqrt(eval_apflux(lbdsp) * eval_apflux(ubdsp)) / std::sqrt(eval_prflux(lbdsm) * eval_prflux(ubdsm));
            double diff   = -2.0 * (spappr - smappr) / (spappr + smappr);
            if (diff <= 0.0) diff = 0.0;
            hspappr->SetBinContent(ir, spappr);
            hsmappr->SetBinContent(ir, smappr);
            hpsm->SetBinContent(ir, diff);
        }
        for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
            hist_rscl->SetBinContent(ir, hpsm->GetBinContent(ir));
        }
    }

    TH1D* hist_magm = new TH1D("hmagm", "hmagm;|Rigidity| [GV/c];Magnetic Field Map", AXrig.nbin(), &(AXrig(0)));
    TH1D* hist_tmpc = new TH1D("htmpc", "htmpc;|Rigidity| [GV/c];Temperature Correction", AXrig.nbin(), &(AXrig(0)));
    const double err_magm = 0.00025; // magnetic field map for AMS inner tracker
    const double err_tmpc = 0.00100; // temperature correction to the magnetic field map
    for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
        hist_magm->SetBinContent(ir, err_magm);
        hist_tmpc->SetBinContent(ir, err_tmpc);
    }
    
    Hist* hist_summary_rigfunc_error_magm = Hist::New("hist_rigfunc_error_magm", HistAxis(AXrig));
    Hist* hist_summary_rigfunc_error_tmpc = Hist::New("hist_rigfunc_error_tmpc", HistAxis(AXrig));
    Hist* hist_summary_rigfunc_error_mgrt = Hist::New("hist_rigfunc_error_mgrt", HistAxis(AXrig));
    Hist* hist_summary_rigfunc_error_rscl = Hist::New("hist_rigfunc_error_rscl", HistAxis(AXrig));
    Hist* hist_summary_rigfunc_error      = Hist::New("hist_rigfunc_error"     , HistAxis(AXrig));

    for (int ir = 1; ir <= AXrig.nbin(); ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        (*hist_summary_rigfunc_error_magm)()->SetBinContent(ir, 100.0 * hist_magm->GetBinContent(ir));
        (*hist_summary_rigfunc_error_tmpc)()->SetBinContent(ir, 100.0 * hist_tmpc->GetBinContent(ir));
        (*hist_summary_rigfunc_error_mgrt)()->SetBinContent(ir, 100.0 * hist_mgrt->GetBinContent(ir));
        (*hist_summary_rigfunc_error_rscl)()->SetBinContent(ir, 100.0 * hist_rscl->GetBinContent(ir));

        double total_error =
            std::hypot(std::hypot(std::hypot(
                (*hist_summary_rigfunc_error_magm)()->GetBinContent(ir),
                (*hist_summary_rigfunc_error_tmpc)()->GetBinContent(ir)),
                (*hist_summary_rigfunc_error_mgrt)()->GetBinContent(ir)),
                (*hist_summary_rigfunc_error_rscl)()->GetBinContent(ir));
        (*hist_summary_rigfunc_error)()->SetBinContent(ir, total_error);
    }

    TFile * ofle = new TFile("out/apflux_acc.root", "RECREATE");
    ofle->cd();

    hist_mcrr->Write();
    hist_mcrrORG->Write();
    
    (*hist_summary_accptance_factor_smooth)()->Write();
    (*hist_summary_accptance_factor_raw)()->Write();
    (*hist_summary_accptance_factor)()->Write();
    (*hist_summary_accptance_error_accp)()->Write();
    (*hist_summary_accptance_error_crss)()->Write();
    (*hist_summary_accptance_error)()->Write(); 
    
    (*hist_summary_rigfunc_error_magm)()->Write();
    (*hist_summary_rigfunc_error_tmpc)()->Write();
    (*hist_summary_rigfunc_error_mgrt)()->Write();
    (*hist_summary_rigfunc_error_rscl)()->Write();
    (*hist_summary_rigfunc_error)()->Write(); 

    ofle->Write();
    ofle->Close();

    PdfEditor editor(Window(WindowSize::kWideSliceLR), "apflux_acc", "out");
   
    THStack* hstack_error_accp = Hist::Collect("hstack_error_accp", HistList({ hist_summary_accptance_error, hist_summary_accptance_error_accp, hist_summary_accptance_error_crss }));
    hist_summary_accptance_error     ->style(Line(kRed     , 0, 4), Marker(kRed    , MarkerStyle(MarkerShape::kCircle)));
    hist_summary_accptance_error_crss->style(Line(kBlue    , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle)));
    hist_summary_accptance_error_accp->style(Line(kGreen+2 , 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle)));

    THStack* hstack_error_rfnc = Hist::Collect("hstack_error_rfnc", HistList({ hist_summary_rigfunc_error, hist_summary_rigfunc_error_tmpc, hist_summary_rigfunc_error_magm, hist_summary_rigfunc_error_mgrt, hist_summary_rigfunc_error_rscl }));
    hist_summary_rigfunc_error     ->style(Line(kRed      , 0, 4), Marker(kRed      , MarkerStyle(MarkerShape::kCircle)));
    hist_summary_rigfunc_error_rscl->style(Line(kBlue     , 0, 2), Marker(kBlue     , MarkerStyle(MarkerShape::kCircle)));
    hist_summary_rigfunc_error_mgrt->style(Line(kGreen+2  , 0, 2), Marker(kGreen+2  , MarkerStyle(MarkerShape::kCircle)));
    hist_summary_rigfunc_error_magm->style(Line(kMagenta+1, 0, 2), Marker(kMagenta+1, MarkerStyle(MarkerShape::kCircle)));
    hist_summary_rigfunc_error_tmpc->style(Line(kYellow+1 , 0, 2), Marker(kYellow+1 , MarkerStyle(MarkerShape::kCircle)));

    TGraphAsymmErrors* graph_summary_accptance_factor_smooth = new TGraphAsymmErrors((*hist_summary_accptance_factor_smooth)());
    graph_summary_accptance_factor_smooth->SetMarkerColor(kRed);
    graph_summary_accptance_factor_smooth->SetLineColor(kRed);
    graph_summary_accptance_factor_smooth->SetFillColor(kRed);
    graph_summary_accptance_factor_smooth->SetFillStyle(3002);
    
    hist_summary_accptance_factor_raw->style(Line(kBlack, 0, 2), Marker(kBlack, MarkerStyle(MarkerShape::kCircle )));
    hist_summary_accptance_factor_smooth->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
    hist_summary_accptance_factor->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    graph_summary_accptance_factor_smooth->Draw("a3");
    graph_summary_accptance_factor_smooth->GetXaxis()->SetMoreLogLabels();
    graph_summary_accptance_factor_smooth->GetXaxis()->SetTitle("|Rigidity| [GV]");
    graph_summary_accptance_factor_smooth->GetYaxis()->SetTitle("Effective Acceptance Correction Factor (A^{#bar{p}/p}_{CF})");
    graph_summary_accptance_factor_smooth->Draw("a3");
    (*hist_summary_accptance_factor_raw)()->Draw("pe same");
    (*hist_summary_accptance_factor)()->Draw("l same");
    Legend leg_accp("", TextStyle(kBlack, 40, 43), PadWindow(0.50, 0.90, 0.65, 0.85));
    leg_accp()->AddEntry((*hist_summary_accptance_factor_raw)()   , "Data", "lp");
    leg_accp()->AddEntry((*hist_summary_accptance_factor_smooth)(), "Parameterized", "l");
    leg_accp()->SetFillColor(0);
    leg_accp.draw();
    editor.save();

    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hstack_error_accp->Draw("nostack hist");
    hstack_error_accp->GetHistogram()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    hstack_error_accp->GetHistogram()->GetYaxis()->SetTitle("Relative Error (%)");
    hstack_error_accp->Draw("nostack hist");
    Legend leg_err_accp("", TextStyle(kBlack, 20, 43), PadWindow(0.50, 0.90, 0.70, 0.90));
    leg_err_accp()->AddEntry((*hist_summary_accptance_error     )(), "Total Systematic Errors", "l");
    leg_err_accp()->AddEntry((*hist_summary_accptance_error_accp)(), "Effective Acceptance Correction", "l");
    leg_err_accp()->AddEntry((*hist_summary_accptance_error_crss)(), "Interaction Cross Sections", "l");
    leg_err_accp()->SetFillColor(0);
    leg_err_accp.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hstack_error_rfnc->Draw("nostack hist");
    hstack_error_rfnc->GetHistogram()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    hstack_error_rfnc->GetHistogram()->GetYaxis()->SetTitle("Relative Error (%)");
    hstack_error_rfnc->Draw("nostack hist");
    Legend leg_err_rfnc("", TextStyle(kBlack, 20, 43), PadWindow(0.30, 0.70, 0.70, 0.90));
    leg_err_rfnc()->AddEntry((*hist_summary_rigfunc_error     )(), "Total Systematic Errors", "lp");
    leg_err_rfnc()->AddEntry((*hist_summary_rigfunc_error_rscl)(), "Rigidity Scale"  , "l");
    leg_err_rfnc()->AddEntry((*hist_summary_rigfunc_error_mgrt)(), "Migration Effect", "l");
    leg_err_rfnc()->AddEntry((*hist_summary_rigfunc_error_tmpc)(), "Temperature Correction", "l");
    leg_err_rfnc()->AddEntry((*hist_summary_rigfunc_error_magm)(), "Magnetic Field Map", "l");
    leg_err_rfnc()->SetFillColor(0);
    leg_err_rfnc.draw();
    editor.save();
    
    editor.close();

    return 1;
}
