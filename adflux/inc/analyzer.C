#ifndef __Analyzer_C__
#define __Analyzer_C__

#include "analyzer.h"

using namespace MGROOT;
        
static constexpr double Mproton    = 0.938272297;
static constexpr double Mdeuterium = 1.876123915;

static constexpr double MPlw = 0.8;
static constexpr double MPup = 1.3;

static constexpr double MDlw = 1.6;
static constexpr double MDup = 2.5;

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

void Analyzer::set_environment() {
    std::cout << Format("\n====  Set Environment ====\n");
    LOG(INFO) << Format("\n====  Set Environment ====\n");
}


void Analyzer::process_events() {
    if (!build_hist()) return;
    
    Stopwatch stopwatch;
    std::string statement_start;
    statement_start += Format("\n**----------------------**\n");
	statement_start += Format("**  Process Alls START  **\n");
	statement_start += Format("**----------------------**\n\n");
    std::cout << statement_start;
    LOG(INFO) << statement_start;
    
    std::cout << Format("\n----==== TF ====----\n");
    LOG(INFO) << Format("\n----==== TF ====----\n");
    process_data(varsTF.get_chain(), &Analyzer::process_data_tf); 
    
    std::cout << Format("\n----==== RH ====----\n");
    LOG(INFO) << Format("\n----==== RH ====----\n");
    process_data(varsRH.get_chain(), &Analyzer::process_data_rh); 
    
    std::string statement_end;
    statement_end += Format("\n**--------------------**\n");
	statement_end += Format("**  Process Alls END  **\n");
	statement_end += Format("**--------------------**\n\n");
    std::cout << statement_end;
    LOG(INFO) << statement_end;

    stopwatch.stop();
    std::string statement_time = stopwatch.str();
    std::cout << statement_time;
    LOG(INFO) << statement_time;
}


bool Analyzer::process_data(TChain* mdst, bool (Analyzer::*process)()) {
    if (mdst == nullptr || mdst->GetEntries() == 0) return false;
    if (file == nullptr) return false;
    Stopwatch stopwatch;

    std::string statement_start;
    statement_start += Format("\n**------------------------**\n");
	statement_start += Format("**  Process Events START  **\n");
	statement_start += Format("**------------------------**\n\n");
    std::cout << statement_start;
    LOG(INFO) << statement_start;

    long num_passed = 0;
	long num_processed = 0;
    long loop_entries = mdst->GetEntries();
	long print_rate =  loop_entries / 5;
	for (long ientry = 0; ientry < loop_entries; ++ientry) {
		if (num_processed%print_rate == 0) {
            stopwatch.stop();
            std::string statement_mid;
            statement_mid += Format("\nInfo :: %lf %\n"                              , 100. * float(num_processed) / float(loop_entries));
			statement_mid += Format("        Passed / Processed : %ld / %ld / %ld\n" , num_passed, num_processed, loop_entries);
			statement_mid += Format("        Passed Ratio       : %lf %\n"           , ((num_processed == 0) ? 0. : (100. * float(num_passed) / float(num_processed))));
			statement_mid += Format("        Real Time          : %9.2f (second)\n"  , stopwatch.time());
			statement_mid += Format("        Processed Rate     : %8.2f (Hz)\n"      , num_processed / stopwatch.time());
            std::cout << statement_mid;
            LOG(INFO) << statement_mid;
        }
        num_processed++;
    
        //if (ientry > 50000) break;
        mdst->GetEvent(ientry);

        if (!(this->*process)()) continue;

        num_passed++;
    }
    
    stopwatch.stop();
    std::string statement_fin;
    statement_fin += Format("\nInfo :: %lf %\n"                              , 100. * float(num_processed) / float(loop_entries));
	statement_fin += Format("        Passed / Processed : %ld / %ld / %ld\n" , num_passed, num_processed, loop_entries);
	statement_fin += Format("        Passed Ratio       : %lf %\n"           , ((num_processed == 0) ? 0. : (100. * float(num_passed) / float(num_processed))));
	statement_fin += Format("        Real Time          : %9.2f (second)\n"  , stopwatch.time());
	statement_fin += Format("        Processed Rate     : %8.2f (Hz)\n"      , num_processed / stopwatch.time());
    std::cout << statement_fin;
    LOG(INFO) << statement_fin;
    
    std::string statement_end;
    statement_end += Format("\n**----------------------**\n");
	statement_end += Format("**  Process Events END  **\n");
	statement_end += Format("**----------------------**\n\n");
    std::cout << statement_end;
    LOG(INFO) << statement_end;

    stopwatch.stop();
    std::string statement_time = stopwatch.str();
    std::cout << statement_time;
    LOG(INFO) << statement_time;

    return true;
}


bool Analyzer::build_hist() {
    Hist::New("hRIG_expt", (TH1D*)(TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/adflux/out/lvtme/RIGlvtme.root")->Get("hRIG_expt")));
    file->cd();

    // TOF
    Axis AXTFbta("", 30, 0.50, 0.80);
    std::vector<double> AXTFrig_bins;
    std::vector<double> AXTFken_bins;
    for (int ii = 0; ii <= AXTFbta.nbin(); ++ii) {
        double ibta = (1.0 / AXTFbta(ii));
        double rig_fact = 1.0 / std::sqrt((ibta - 1.0) * (ibta + 1.0));
        double rig = Mdeuterium * rig_fact;
        AXTFrig_bins.push_back(rig);

        double ken_fact = (std::hypot(1.0, rig_fact) - 1.0);
        double ken = ken_fact * Mproton;
        AXTFken_bins.push_back(ken);
    }
    Axis AXTFrig("Rigidity [GV]", AXTFrig_bins);
    Axis AXTFken("Kinetic Energy per Nucleon [GeV/n]", AXTFken_bins);
    
    // RICH
    Axis AXRHbta("", 20, 0.96, 0.98);
    std::vector<double> AXRHrig_bins;
    std::vector<double> AXRHken_bins;
    for (int ii = 0; ii <= AXRHbta.nbin(); ++ii) {
        double ibta = (1.0 / AXRHbta(ii));
        double rig_fact = 1.0 / std::sqrt((ibta - 1.0) * (ibta + 1.0));
        double rig = Mdeuterium * rig_fact;
        AXRHrig_bins.push_back(rig);

        double ken_fact = (std::hypot(1.0, rig_fact) - 1.0);
        double ken = ken_fact * Mproton;
        AXRHken_bins.push_back(ken);
    }
    Axis AXRHrig("Rigidity [GV]", AXRHrig_bins);
    Axis AXRHken("Kinetic Energy per Nucleon [GeV/n]", AXRHken_bins);

    
    Axis AXTFmass("Mass/Z [(GV/c^{2})]", 240, -3.0, 3.0);
    Axis AXTFllr("LLR", 160, 0.0, 1.6);
    Axis AXTFchiFR("Log(#chi^{2}/NDF)  [Free-Mass]", 100, -2.0, 8.0);
    Axis AXTFchiFX("Log(#chi^{2}/NDF)  [Fixed-Mass]", 100, -2.0, 8.0);
    
    Hist::New("hTF_llr" , HistAxis(AXTFmass, AXTFllr));

    Hist::New("hTF_mass_NO_PRE"     , HistAxis(AXTFmass));
    Hist::New("hTF_mass_NO_SEL_MUTR", HistAxis(AXTFmass));
    Hist::New("hTF_mass_NO_SEL_PHYS", HistAxis(AXTFmass));
    Hist::New("hTF_mass_NO"         , HistAxis(AXTFmass));
    Hist::New("hTF_mass_CF_PRE"     , HistAxis(AXTFmass));
    Hist::New("hTF_mass_CF_SEL_MUTR", HistAxis(AXTFmass));
    Hist::New("hTF_mass_CF_SEL_PHYS", HistAxis(AXTFmass));
    Hist::New("hTF_mass_CF"         , HistAxis(AXTFmass));
    
    Hist::New("hTF_mass_de", HistAxis(AXTFmass));
    Hist::New("hTF_mass_pr", HistAxis(AXTFmass));
    Hist::New("hTF_mass_ap", HistAxis(AXTFmass));
    Hist::New("hTF_mass_pp", HistAxis(AXTFmass));

    Hist::New("hTF_mutr_lx"  , HistAxis(AXTFmass, AXTFchiFR));
    Hist::New("hTF_mutr_ly"  , HistAxis(AXTFmass, AXTFchiFR));
    Hist::New("hTF_mutr_lb"  , HistAxis(AXTFmass, AXTFchiFR));
    Hist::New("hTF_mutr_lxyb", HistAxis(AXTFmass, AXTFchiFR));
    
    Hist::New("hTF_phys_lx"  , HistAxis(AXTFmass, AXTFchiFX));
    Hist::New("hTF_phys_ly"  , HistAxis(AXTFmass, AXTFchiFX));
    Hist::New("hTF_phys_lb"  , HistAxis(AXTFmass, AXTFchiFX));
    Hist::New("hTF_phys_lxyb", HistAxis(AXTFmass, AXTFchiFX));
    
    Hist::New("hTFrig_trgn", HistAxis(AXTFrig));
    Hist::New("hTFrig_trgd", HistAxis(AXTFrig));
    Hist::New("hTFken_trgn", HistAxis(AXTFken));
    Hist::New("hTFken_trgd", HistAxis(AXTFken));
    
    Hist::New("hTFrig_effn_MC", HistAxis(AXTFrig));
    Hist::New("hTFrig_effd_MC", HistAxis(AXTFrig));
    Hist::New("hTFken_effn_MC", HistAxis(AXTFken));
    Hist::New("hTFken_effd_MC", HistAxis(AXTFken));
    
    Hist::New("hTFrig_cnt"   , HistAxis(AXTFrig));
    Hist::New("hTFrig_cnt_MC", HistAxis(AXTFrig));
    Hist::New("hTFken_cnt"   , HistAxis(AXTFken));
    Hist::New("hTFken_cnt_MC", HistAxis(AXTFken));
    
    Axis AXRHmass("Mass/Z [(GV/c^{2})]", 240, -3.0, 3.0);
    Axis AXRHllr("LLR", 160, 0.0, 1.6);
    Axis AXRHchiFR("Log(#chi^{2}/NDF)  [Free-Mass]", 100, -2.0, 8.0);
    Axis AXRHchiFX("Log(#chi^{2}/NDF)  [Fixed-Mass]", 100, -2.0, 8.0);
    Axis AXRHbta2("#beta", 400, 0.96, 1.01);
    Axis AXRHq("Q", 400, 0.0, 4.00);
    
    Hist::New("hRH_llr" , HistAxis(AXRHmass, AXRHllr));
    Hist::New("hRH_mass_q" , HistAxis(AXRHmass, AXRHq));
    
    Hist::New("hRH_mass_NO_PRE"     , HistAxis(AXRHmass));
    Hist::New("hRH_mass_NO_SEL_MUTR", HistAxis(AXRHmass));
    Hist::New("hRH_mass_NO_SEL_PHYS", HistAxis(AXRHmass));
    Hist::New("hRH_mass_NO"         , HistAxis(AXRHmass));
    Hist::New("hRH_mass_CF_PRE"     , HistAxis(AXRHmass));
    Hist::New("hRH_mass_CF_SEL_MUTR", HistAxis(AXRHmass));
    Hist::New("hRH_mass_CF_SEL_PHYS", HistAxis(AXRHmass));
    Hist::New("hRH_mass_CF"         , HistAxis(AXRHmass));
    
    Hist::New("hRH_mass_NO_V2"      , HistAxis(AXRHmass));
    Hist::New("hRH_mass_CF_V2"      , HistAxis(AXRHmass));
    
    Hist::New("hRH_mass_de", HistAxis(AXRHmass));
    Hist::New("hRH_mass_pr", HistAxis(AXRHmass));
    Hist::New("hRH_mass_ap", HistAxis(AXRHmass));
    Hist::New("hRH_mass_pp", HistAxis(AXRHmass));
    
    Hist::New("hRH_mass_de_V2", HistAxis(AXRHmass));
    Hist::New("hRH_mass_pr_V2", HistAxis(AXRHmass));
    Hist::New("hRH_mass_ap_V2", HistAxis(AXRHmass));
    Hist::New("hRH_mass_pp_V2", HistAxis(AXRHmass));
    
    Hist::New("hRH_tf_mutr_lx"  , HistAxis(AXRHmass, AXRHchiFR));
    Hist::New("hRH_tf_mutr_ly"  , HistAxis(AXRHmass, AXRHchiFR));
    Hist::New("hRH_tf_mutr_lb"  , HistAxis(AXRHmass, AXRHchiFR));
    Hist::New("hRH_tf_mutr_lxyb", HistAxis(AXRHmass, AXRHchiFR));
    
    Hist::New("hRH_rh_mutr_lx"  , HistAxis(AXRHmass, AXRHchiFR));
    Hist::New("hRH_rh_mutr_ly"  , HistAxis(AXRHmass, AXRHchiFR));
    Hist::New("hRH_rh_mutr_lb"  , HistAxis(AXRHmass, AXRHchiFR));
    Hist::New("hRH_rh_mutr_lxyb", HistAxis(AXRHmass, AXRHchiFR));
    
    Hist::New("hRH_tf_phys_lx"  , HistAxis(AXRHmass, AXRHchiFX));
    Hist::New("hRH_tf_phys_ly"  , HistAxis(AXRHmass, AXRHchiFX));
    Hist::New("hRH_tf_phys_lb"  , HistAxis(AXRHmass, AXRHchiFX));
    Hist::New("hRH_tf_phys_lxyb", HistAxis(AXRHmass, AXRHchiFX));
    
    Hist::New("hRH_rh_phys_lx"  , HistAxis(AXRHmass, AXRHchiFX));
    Hist::New("hRH_rh_phys_ly"  , HistAxis(AXRHmass, AXRHchiFX));
    Hist::New("hRH_rh_phys_lb"  , HistAxis(AXRHmass, AXRHchiFX));
    Hist::New("hRH_rh_phys_lxyb", HistAxis(AXRHmass, AXRHchiFX));
    
    Hist::New("hRHbta"   , HistAxis(AXRHbta2));
    Hist::New("hRHbta_de", HistAxis(AXRHbta2));
    Hist::New("hRHbta_pr", HistAxis(AXRHbta2));
    Hist::New("hRHbta_ap", HistAxis(AXRHbta2));
    
    Hist::New("hRHbta2"   , HistAxis(AXRHbta2));
    Hist::New("hRHbta2_de", HistAxis(AXRHbta2));
    Hist::New("hRHbta2_pr", HistAxis(AXRHbta2));
    Hist::New("hRHbta2_ap", HistAxis(AXRHbta2));
    
    Hist::New("hRHrig_trgn", HistAxis(AXRHrig));
    Hist::New("hRHrig_trgd", HistAxis(AXRHrig));
    Hist::New("hRHken_trgn", HistAxis(AXRHken));
    Hist::New("hRHken_trgd", HistAxis(AXRHken));
    
    Hist::New("hRHrig_effn_MC", HistAxis(AXRHrig));
    Hist::New("hRHrig_effd_MC", HistAxis(AXRHrig));
    Hist::New("hRHken_effn_MC", HistAxis(AXRHken));
    Hist::New("hRHken_effd_MC", HistAxis(AXRHken));
    
    Hist::New("hRHrig_cnt"   , HistAxis(AXRHrig));
    Hist::New("hRHrig_cnt_MC", HistAxis(AXRHrig));
    Hist::New("hRHken_cnt"   , HistAxis(AXRHken));
    Hist::New("hRHken_cnt_MC", HistAxis(AXRHken));
    
    Hist::New("hRHrig_trgn_V2", HistAxis(AXRHrig));
    Hist::New("hRHrig_trgd_V2", HistAxis(AXRHrig));
    Hist::New("hRHken_trgn_V2", HistAxis(AXRHken));
    Hist::New("hRHken_trgd_V2", HistAxis(AXRHken));
    
    Hist::New("hRHrig_effn_MC_V2", HistAxis(AXRHrig));
    Hist::New("hRHrig_effd_MC_V2", HistAxis(AXRHrig));
    Hist::New("hRHken_effn_MC_V2", HistAxis(AXRHken));
    Hist::New("hRHken_effd_MC_V2", HistAxis(AXRHken));
    
    Hist::New("hRHrig_cnt_V2"   , HistAxis(AXRHrig));
    Hist::New("hRHrig_cnt_MC_V2", HistAxis(AXRHrig));
    Hist::New("hRHken_cnt_V2"   , HistAxis(AXRHken));
    Hist::New("hRHken_cnt_MC_V2", HistAxis(AXRHken));
    
    return true;
}


bool Analyzer::process_data_tf() {
    VarsTF* vars = &varsTF;
    
    short  sign = vars->sign;
    double bta  = vars->b_top_bta[1];
    double rig  = vars->b_top_rig[1];
    double arig = std::abs(rig);
    double cfr  = vars->mc ? 0.0 : (arig / vars->cfsec);
   /* 
    if (vars->nhx < 4) return false;
    if (vars->nhy < 5) return false;
    
    if (vars->ext && vars->extlx > 3.0) return false;
    if (vars->ext && vars->extly > 3.0) return false;

    if (!vars->veto) return false; 
    if (vars->ncls != 0) return false; 
    
    if (vars->nhinn > 1) return false;
    if (vars->nhout > 3) return false;
    
    if (vars->nhinn2 > 1) return false;
    if (vars->nhout2 > 3) return false;
    
    if ((vars->nvtxx[0]+vars->nvtxx[1]) > 10) return false;
    if ((vars->nvtxy[0]+vars->nvtxy[1]) > 15) return false;
    */
    double mutr_lx   = std::max(vars->mutr_lx[0], vars->mutr_lx[1]);
    double mutr_ly   = std::max(vars->mutr_ly[0], vars->mutr_ly[1]);
    double mutr_lb   = std::max(vars->mutr_lb[0], vars->mutr_lb[1]);
    double mutr_lxyb = std::max(std::max(mutr_lx, mutr_ly), mutr_lb);
    
    double phys_lx   = std::max(std::min(vars->a_phys_lx[0], vars->b_phys_lx[0]),  std::min(vars->a_phys_lx[1], vars->b_phys_lx[1]));
    double phys_ly   = std::max(std::min(vars->a_phys_ly[0], vars->b_phys_ly[0]),  std::min(vars->a_phys_ly[1], vars->b_phys_ly[1]));
    double phys_lb   = std::max(std::min(vars->a_phys_lb[0], vars->b_phys_lb[0]),  std::min(vars->a_phys_lb[1], vars->b_phys_lb[1]));
    double phys_lxyb = std::max(std::max(phys_lx, phys_ly), phys_lb);

    // select
    bool is_in_windows = ((vars->bta > 0.5 && vars->bta < 0.8) && vars->sqrm > 1.0e-6);
    //bool isovercf_iss = (vars->mc || cfr > 1.2);
    //bool isovercf_iss = (vars->mc || cfr < 0.8);
    bool isovercf_iss = true; // testcode
    
    double MCrig = (vars->mc) ? std::abs(vars->mc_rig) : 0.0;
    double MCken = (vars->mc) ? Mproton * (std::hypot(1.0, std::abs(vars->mc_rig / Mdeuterium)) - 1.0) : 0.0;
    double ken   = Mproton * (std::hypot(1.0, (arig / Mdeuterium)) - 1.0);
    double mass  = vars->sign * std::sqrt(vars->sqrm);
    
    TH1D* hRIG_expt = (TH1D*) (*Hist::Head("hRIG_expt"))();
    double wgtcf = vars->mc ? hRIG_expt->GetBinContent(hRIG_expt->FindBin(arig)) * 1.0e-7 : 1.0;
    double wgtde = vars->mc ? wgtcf * eval_prflux(MCrig, -1.7) / eval_prflux(1.0, -1.7) : 1.0;
    double wgtpr = vars->mc ? wgtcf * eval_prflux(MCrig, -1.7) / eval_prflux(1.0, -1.7) : 1.0;
    double wgtap = vars->mc ? wgtcf * eval_apflux(MCrig, -1.7) / eval_apflux(1.0, -1.7) : 1.0;
    bool isovercf_mc = (vars->mc && wgtcf > 0.0);
    
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_llr")->fillH2D(mass, (vars->llr / vars->bta), vars->wgt);
    //if ((vars->llr / vars->bta) < 0.70) return false;
    
    if (is_in_windows && vars->trg) Hist::Head("hTF_mass_NO_PRE")->fillH1D(mass, vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_mass_CF_PRE")->fillH1D(mass, vars->wgt);
    
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_mutr_lx")  ->fillH2D(mass, mutr_lx  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_mutr_ly")  ->fillH2D(mass, mutr_ly  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_mutr_lb")  ->fillH2D(mass, mutr_lb  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_mutr_lxyb")->fillH2D(mass, mutr_lxyb, vars->wgt);
    
    //if (mutr_lxyb > 1.75) return false;
    //if (mutr_lxyb > 2.5) return false;
    
    if (is_in_windows && vars->trg) Hist::Head("hTF_mass_NO_SEL_MUTR")->fillH1D(mass, vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_mass_CF_SEL_MUTR")->fillH1D(mass, vars->wgt);

    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_phys_lx")  ->fillH2D(mass, phys_lx  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_phys_ly")  ->fillH2D(mass, phys_ly  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_phys_lb")  ->fillH2D(mass, phys_lb  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_phys_lxyb")->fillH2D(mass, phys_lxyb, vars->wgt);

    //if (phys_lxyb > 1.75) return false;
    //if (phys_lxyb > 2.5) return false;
    
    if (is_in_windows && vars->trg) Hist::Head("hTF_mass_NO_SEL_PHYS")->fillH1D(mass, vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_mass_CF_SEL_PHYS")->fillH1D(mass, vars->wgt);
    
    if (is_in_windows && vars->trg) Hist::Head("hTF_mass_NO")->fillH1D(mass, vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hTF_mass_CF")->fillH1D(mass, vars->wgt);
    if (!isovercf_iss) return false;
    
    if (is_in_windows && mass > MDlw && mass < MDup && vars->trg) Hist::Head("hTFrig_trgn")->fillH1D(arig, vars->wgt);
    if (is_in_windows && mass > MDlw && mass < MDup && vars->trg) Hist::Head("hTFken_trgn")->fillH1D(ken , vars->wgt);
    if (is_in_windows && mass > MDlw && mass < MDup) Hist::Head("hTFrig_trgd")->fillH1D(arig, vars->wgt * (vars->trg?1.0:100.0));
    if (is_in_windows && mass > MDlw && mass < MDup) Hist::Head("hTFken_trgd")->fillH1D(ken , vars->wgt * (vars->trg?1.0:100.0));
    if (!vars->trg) return false;

    if (is_in_windows && isovercf_mc) Hist::Head("hTF_mass_de")->fillH1D( mass, wgtde * vars->wgt);
    if (is_in_windows && isovercf_mc) Hist::Head("hTF_mass_pr")->fillH1D( mass, wgtpr * vars->wgt);
    if (is_in_windows && isovercf_mc) Hist::Head("hTF_mass_ap")->fillH1D( mass, wgtap * vars->wgt);
    if (is_in_windows && isovercf_mc) Hist::Head("hTF_mass_pp")->fillH1D(-mass, wgtap * vars->wgt);
    
    if (vars->mc && mass < -MDlw && mass > -MDup) Hist::Head("hTFrig_effn_MC")->fillH1D(MCrig, vars->wgt);
    if (vars->mc && mass < -MDlw && mass > -MDup) Hist::Head("hTFken_effn_MC")->fillH1D(MCken, vars->wgt);
    if (vars->mc) Hist::Head("hTFrig_effd_MC")->fillH1D(MCrig, vars->wgt);
    if (vars->mc) Hist::Head("hTFken_effd_MC")->fillH1D(MCken, vars->wgt);

    Hist::Head("hTFrig_cnt")->fillH1D(arig, vars->wgt);
    Hist::Head("hTFken_cnt")->fillH1D(ken , vars->wgt);
    if (vars->mc) Hist::Head("hTFrig_cnt_MC")->fillH1D(MCrig, vars->wgt);
    if (vars->mc) Hist::Head("hTFken_cnt_MC")->fillH1D(MCken, vars->wgt);

    if (is_in_windows && (mass < -MDlw && mass > -MDup)) {
        adtreeTF->Fill();
    }

    return true;
}


bool Analyzer::process_data_rh() {
    VarsRH* vars = &varsRH;
    
    short  sign = vars->sign;
    double bta  = vars->b_top_bta[3];
    double rig  = vars->b_top_rig[3];
    double arig = std::abs(rig);
    double cfr  = vars->mc ? 0.0 : (arig / vars->cfsec);
    
    if (vars->nhx < 4) return false;
    if (vars->nhy < 5) return false;
    if (vars->L2 != 3) return false;
    
    if (vars->ext && vars->extlx > 3.0) return false;
    if (vars->ext && vars->extly > 3.0) return false;
  
    if (std::sqrt(vars->selfqq) > 2.5) return false;
    if (vars->nclsZ1 == 0) return false;
    if (vars->nclsZ1 >  2) return false;

    if (vars->nhinn > 1) return false;
    if (vars->nhout > 2) return false;
    
    if (vars->nhinn2 > 1) return false;
    if (vars->nhout2 > 2) return false;
    
    if ((vars->nvtxx[0]+vars->nvtxx[1]) > 3) return false;
    if ((vars->nvtxy[0]+vars->nvtxy[1]) > 8) return false;

    double tf_mutr_lx   = std::max(vars->mutr_lx[0], vars->mutr_lx[1]);
    double tf_mutr_ly   = std::max(vars->mutr_ly[0], vars->mutr_ly[1]);
    double tf_mutr_lb   = std::max(vars->mutr_lb[0], vars->mutr_lb[1]);
    double tf_mutr_lxyb = std::max(std::max(tf_mutr_lx, tf_mutr_ly), tf_mutr_lb);
    
    double rh_mutr_lx   = std::max(vars->mutr_lx[2], vars->mutr_lx[3]);
    double rh_mutr_ly   = std::max(vars->mutr_ly[2], vars->mutr_ly[3]);
    double rh_mutr_lb   = std::max(vars->mutr_lb[2], vars->mutr_lb[3]);
    double rh_mutr_lxyb = std::max(std::max(rh_mutr_lx, rh_mutr_ly), rh_mutr_lb);

    double tf_phys_lx   = std::max(std::min(vars->a_phys_lx[0], vars->b_phys_lx[0]),  std::min(vars->a_phys_lx[1], vars->b_phys_lx[1]));
    double tf_phys_ly   = std::max(std::min(vars->a_phys_ly[0], vars->b_phys_ly[0]),  std::min(vars->a_phys_ly[1], vars->b_phys_ly[1]));
    double tf_phys_lb   = std::max(std::min(vars->a_phys_lb[0], vars->b_phys_lb[0]),  std::min(vars->a_phys_lb[1], vars->b_phys_lb[1]));
    double tf_phys_lxyb = std::max(std::max(tf_phys_lx, tf_phys_ly), tf_phys_lb);
    
    double rh_phys_lx   = std::max(std::min(vars->a_phys_lx[2], vars->b_phys_lx[2]),  std::min(vars->a_phys_lx[3], vars->b_phys_lx[3]));
    double rh_phys_ly   = std::max(std::min(vars->a_phys_ly[2], vars->b_phys_ly[2]),  std::min(vars->a_phys_ly[3], vars->b_phys_ly[3]));
    double rh_phys_lb   = std::max(std::min(vars->a_phys_lb[2], vars->b_phys_lb[2]),  std::min(vars->a_phys_lb[3], vars->b_phys_lb[3]));
    double rh_phys_lxyb = std::max(std::max(rh_phys_lx, rh_phys_ly), rh_phys_lb);
  
    // select
    bool is_in_windows = ((vars->bta > 0.96 && vars->bta < 0.98) && vars->sqrm > 1.0e-6);
    //bool isovercf_iss = (vars->mc || cfr > 1.2);
    bool isovercf_iss = (vars->mc || cfr < 0.8);
    
    double mass  = vars->sign * std::sqrt(vars->sqrm);
    double ken   = Mproton * (std::hypot(1.0, (arig / Mdeuterium)) - 1.0);
    double MCrig = (vars->mc) ? std::abs(vars->mc_rig) : 0.0;
    double MCken = (vars->mc) ? Mproton * (std::hypot(1.0, std::abs(vars->mc_rig / Mdeuterium)) - 1.0) : 0.0;
   
    TH1D* hRIG_expt = (TH1D*) (*Hist::Head("hRIG_expt"))();
    double wgtcf = vars->mc ? hRIG_expt->GetBinContent(hRIG_expt->FindBin(arig)) * 1.0e-6 : 1.0;
    double wgtde = vars->mc ? wgtcf * eval_prflux(MCrig, -1.7) / eval_prflux(3.0, -1.7) : 1.0;
    double wgtpr = vars->mc ? wgtcf * eval_prflux(MCrig, -1.7) / eval_prflux(3.0, -1.7) : 1.0;
    double wgtap = vars->mc ? wgtcf * eval_apflux(MCrig, -1.7) / eval_apflux(3.0, -1.7) : 1.0;
    bool isovercf_mc = (vars->mc && wgtcf > 0.0);

    if (arig > 15 && isovercf_iss) Hist::Head("hRHbta")->fillH1D(vars->rhbta[1], vars->wgt);
    if (arig > 15 && isovercf_mc) Hist::Head("hRHbta_de")->fillH1D(vars->rhbta[1], wgtde * vars->wgt);
    if (arig > 15 && isovercf_mc) Hist::Head("hRHbta_pr")->fillH1D(vars->rhbta[1], wgtpr * vars->wgt);
    if (arig > 15 && isovercf_mc) Hist::Head("hRHbta_ap")->fillH1D(vars->rhbta[1], wgtap * vars->wgt);
    
    if (arig > 15 && isovercf_iss) Hist::Head("hRHbta2")->fillH1D(vars->rhbta[0], vars->wgt);
    if (arig > 15 && isovercf_mc) Hist::Head("hRHbta2_de")->fillH1D(vars->rhbta[0], wgtde * vars->wgt);
    if (arig > 15 && isovercf_mc) Hist::Head("hRHbta2_pr")->fillH1D(vars->rhbta[0], wgtpr * vars->wgt);
    if (arig > 15 && isovercf_mc) Hist::Head("hRHbta2_ap")->fillH1D(vars->rhbta[0], wgtap * vars->wgt);
    
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_llr")->fillH2D(mass, vars->llr, vars->wgt);
    if (is_in_windows && vars->llr < 0.75) return false;
    
    if (is_in_windows && vars->trg) Hist::Head("hRH_mass_NO_PRE")->fillH1D(mass, vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_mass_CF_PRE")->fillH1D(mass, vars->wgt);
    
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_tf_mutr_lx"  )->fillH2D(mass, tf_mutr_lx  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_tf_mutr_ly"  )->fillH2D(mass, tf_mutr_ly  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_tf_mutr_lb"  )->fillH2D(mass, tf_mutr_lb  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_tf_mutr_lxyb")->fillH2D(mass, tf_mutr_lxyb, vars->wgt);
    
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_rh_mutr_lx"  )->fillH2D(mass, rh_mutr_lx  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_rh_mutr_ly"  )->fillH2D(mass, rh_mutr_ly  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_rh_mutr_lb"  )->fillH2D(mass, rh_mutr_lb  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_rh_mutr_lxyb")->fillH2D(mass, rh_mutr_lxyb, vars->wgt);
    
    if (tf_mutr_lxyb > 2.00) return false;
    if (rh_mutr_lxyb > 1.75) return false;
    
    if (is_in_windows && vars->trg) Hist::Head("hRH_mass_NO_SEL_MUTR")->fillH1D(mass, vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_mass_CF_SEL_MUTR")->fillH1D(mass, vars->wgt);
    
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_tf_phys_lx"  )->fillH2D(mass, tf_phys_lx  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_tf_phys_ly"  )->fillH2D(mass, tf_phys_ly  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_tf_phys_lb"  )->fillH2D(mass, tf_phys_lb  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_tf_phys_lxyb")->fillH2D(mass, tf_phys_lxyb, vars->wgt);
    
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_rh_phys_lx"  )->fillH2D(mass, rh_phys_lx  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_rh_phys_ly"  )->fillH2D(mass, rh_phys_ly  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_rh_phys_lb"  )->fillH2D(mass, rh_phys_lb  , vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_rh_phys_lxyb")->fillH2D(mass, rh_phys_lxyb, vars->wgt);
    
    if (tf_phys_lxyb > 2.00) return false;
    if (rh_phys_lxyb > 1.75) return false;

    if (is_in_windows && vars->trg) Hist::Head("hRH_mass_NO_SEL_PHYS")->fillH1D(mass, vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_mass_CF_SEL_PHYS")->fillH1D(mass, vars->wgt);

    if (is_in_windows && vars->trg) Hist::Head("hRH_mass_NO")->fillH1D(mass, vars->wgt);
    if (is_in_windows && vars->trg && isovercf_iss) Hist::Head("hRH_mass_CF")->fillH1D(mass, vars->wgt);
    
    if (is_in_windows && rh_phys_lxyb < 1.5 && vars->trg) Hist::Head("hRH_mass_NO_V2")->fillH1D(mass, vars->wgt);
    if (is_in_windows && rh_phys_lxyb < 1.5 && vars->trg && isovercf_iss) Hist::Head("hRH_mass_CF_V2")->fillH1D(mass, vars->wgt);
    if (!isovercf_iss) return false;
    
    if (is_in_windows && mass > MDlw && mass < MDup && vars->trg) Hist::Head("hRHrig_trgn")->fillH1D(arig, vars->wgt);
    if (is_in_windows && mass > MDlw && mass < MDup && vars->trg) Hist::Head("hRHken_trgn")->fillH1D(ken , vars->wgt);
    if (is_in_windows && mass > MDlw && mass < MDup) Hist::Head("hRHrig_trgd")->fillH1D(arig, vars->wgt * (vars->trg?1.0:100.0));
    if (is_in_windows && mass > MDlw && mass < MDup) Hist::Head("hRHken_trgd")->fillH1D(ken , vars->wgt * (vars->trg?1.0:100.0));
   
    if (rh_phys_lxyb < 1.5) {
    if (is_in_windows && mass > MDlw && mass < MDup && vars->trg) Hist::Head("hRHrig_trgn_V2")->fillH1D(arig, vars->wgt);
    if (is_in_windows && mass > MDlw && mass < MDup && vars->trg) Hist::Head("hRHken_trgn_V2")->fillH1D(ken , vars->wgt);
    if (is_in_windows && mass > MDlw && mass < MDup) Hist::Head("hRHrig_trgd_V2")->fillH1D(arig, vars->wgt * (vars->trg?1.0:100.0));
    if (is_in_windows && mass > MDlw && mass < MDup) Hist::Head("hRHken_trgd_V2")->fillH1D(ken , vars->wgt * (vars->trg?1.0:100.0));
    }
    
    //if (!vars->trg) return false;
  

    if (is_in_windows) Hist::Head("hRH_mass_q")->fillH2D(mass, std::sqrt(vars->selfqq), vars->wgt);

    bool signout = (vars->mc && vars->mc_rig > 0 && mass < -MDlw);
    if (is_in_windows && isovercf_mc && !signout) Hist::Head("hRH_mass_de")->fillH1D( mass, wgtde * vars->wgt);
    if (is_in_windows && isovercf_mc && !signout) Hist::Head("hRH_mass_pr")->fillH1D( mass, wgtpr * vars->wgt);
    if (is_in_windows && isovercf_mc && !signout) Hist::Head("hRH_mass_ap")->fillH1D( mass, wgtap * vars->wgt);
    if (is_in_windows && isovercf_mc && !signout) Hist::Head("hRH_mass_pp")->fillH1D(-mass, wgtap * vars->wgt);
    
    if (vars->mc && mass < -MDlw && mass > -MDup) Hist::Head("hRHrig_effn_MC")->fillH1D(MCrig, vars->wgt);
    if (vars->mc && mass < -MDlw && mass > -MDup) Hist::Head("hRHken_effn_MC")->fillH1D(MCken, vars->wgt);
    if (vars->mc) Hist::Head("hRHrig_effd_MC")->fillH1D(MCrig, vars->wgt);
    if (vars->mc) Hist::Head("hRHken_effd_MC")->fillH1D(MCken, vars->wgt);

    Hist::Head("hRHrig_cnt")->fillH1D(arig, vars->wgt);
    Hist::Head("hRHken_cnt")->fillH1D(ken , vars->wgt);
    if (vars->mc) Hist::Head("hRHrig_cnt_MC")->fillH1D(MCrig, vars->wgt);
    if (vars->mc) Hist::Head("hRHken_cnt_MC")->fillH1D(MCken, vars->wgt);
    
    if (rh_phys_lxyb < 1.5) {
    if (is_in_windows && isovercf_mc && !signout) Hist::Head("hRH_mass_de_V2")->fillH1D( mass, wgtde * vars->wgt);
    if (is_in_windows && isovercf_mc && !signout) Hist::Head("hRH_mass_pr_V2")->fillH1D( mass, wgtpr * vars->wgt);
    if (is_in_windows && isovercf_mc && !signout) Hist::Head("hRH_mass_ap_V2")->fillH1D( mass, wgtap * vars->wgt);
    if (is_in_windows && isovercf_mc && !signout) Hist::Head("hRH_mass_pp_V2")->fillH1D(-mass, wgtap * vars->wgt);
    
    if (vars->mc && mass < -MDlw && mass > -MDup) Hist::Head("hRHrig_effn_MC_V2")->fillH1D(MCrig, vars->wgt);
    if (vars->mc && mass < -MDlw && mass > -MDup) Hist::Head("hRHken_effn_MC_V2")->fillH1D(MCken, vars->wgt);
    if (vars->mc) Hist::Head("hRHrig_effd_MC_V2")->fillH1D(MCrig, vars->wgt);
    if (vars->mc) Hist::Head("hRHken_effd_MC_V2")->fillH1D(MCken, vars->wgt);
    
    Hist::Head("hRHrig_cnt_V2")->fillH1D(arig, vars->wgt);
    Hist::Head("hRHken_cnt_V2")->fillH1D(ken , vars->wgt);
    if (vars->mc) Hist::Head("hRHrig_cnt_MC_V2")->fillH1D(MCrig, vars->wgt);
    if (vars->mc) Hist::Head("hRHken_cnt_MC_V2")->fillH1D(MCken, vars->wgt);
    }

    //if (is_in_windows && (mass < -MDlw && mass > -MDup) && rh_phys_lxyb < 1.5) {
    //    adtreeRH->Fill();
    //}

    bool cover =
        (vars->run == 1334894615 && vars->evt == 133569) ||
        (vars->run == 1421791849 && vars->evt ==  67582) ||
        (vars->run == 1439981103 && vars->evt == 409381) ||
        (vars->run == 1456017070 && vars->evt == 574622) ||
        (vars->run == 1489729605 && vars->evt == 808247) ||
        (vars->run == 1562077316 && vars->evt == 527467) ||
        (vars->run == 1571366898 && vars->evt == 391367);
    if (cover) {
        adtreeRH->Fill();
    }

    return true;
}


#endif // __Analyzer_C__
