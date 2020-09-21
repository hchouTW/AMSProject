#ifndef __Analyzer_C__
#define __Analyzer_C__

#include "analyzer.h"

using namespace MGROOT;

void Analyzer::set_environment() {
    std::cout << Format("\n====  Set Environment ====\n");
    LOG(INFO) << Format("\n====  Set Environment ====\n");
    
    //TrSys::PhysEnv::ReadMagAMS("/eos/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    //TrSys::PhysEnv::ReadMatAMS("/eos/user/h/hchou/ExternalLibs/DB/material");
    //
    //if (!TrSys::PhysEnv::IMagStatus()) TrSys::PhysEnv::ReadMagAMS("/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    //if (!TrSys::PhysEnv::IMatStatus()) TrSys::PhysEnv::ReadMatAMS("/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    //
    //if (!TrSys::PhysEnv::IMagStatus()) TrSys::PhysEnv::ReadMagAMS("/afs/cern.ch/work/h/hchou/public/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    //if (!TrSys::PhysEnv::IMatStatus()) TrSys::PhysEnv::ReadMatAMS("/afs/cern.ch/work/h/hchou/public/ExternalLibs/DB/material");
}


void Analyzer::process_events() {
    if (mdst == nullptr || mdst->GetEntries() == 0) return;
    if (file == nullptr || tree == nullptr) return;
    Stopwatch stopwatch;

    if (!build_tree()) return;
    if (!build_hist()) return;

    std::string statement_start;
    statement_start += Format("\n**------------------------**\n");
	statement_start += Format("**  Process Events START  **\n");
	statement_start += Format("**------------------------**\n\n");
    std::cout << statement_start;
    LOG(INFO) << statement_start;

    long num_passed = 0;
	long num_processed = 0;
    long loop_entries = mdst->GetEntries();
	long print_rate =  loop_entries / 50;
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

        if (!process_presel()) continue;
        if (!process_data()) continue;

        if (list != nullptr) {
            if (runev.find(list->run) == runev.end()) runev[list->run] = list->event;
            else runev[list->run] = std::max(runev[list->run], list->event);
        }

        //tree->Fill();
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
}

bool Analyzer::build_tree() {
    return true;
}


static UInt_t utime_cur = 0;
static TGraphAsymmErrors* gTMPpr27 = nullptr;
static TGraphAsymmErrors* gTMPap2pr = nullptr;

bool Analyzer::build_hist() {
    TFile* fTMPpr = TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/antip/others/ams02_pr.root");
    gTMPpr27 = (TGraphAsymmErrors*) (fTMPpr->Get("gr_exp2"))->Clone("gTMPpr27");
    
    TFile* fTMPap2pr = TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/antip/others/ams02_ap2pr.root");
    gTMPap2pr = (TGraphAsymmErrors*) (fTMPap2pr->Get("gr_exp1"))->Clone("gTMPap2pr");
    
    TFile* fTMPcf = TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/antip/others/expt_ISS.root");
    Hist* hTMPcf = Hist::New("hTMPcf", (TH1*) fTMPcf->Get("hPROF_HC_expt_gISS"));

    file->cd();

    std::vector<double> vmom( {
        0.80,
        1.00,    1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
        3.29,    3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
        8.48,    9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
        19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
        41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
        93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 330.00, 525.00,
        800.00, 1200.00, 1800.00, 3600.00 } );
    
    std::vector<double> vmom2( {
        1.00,    2.15,   3.29,   4.43,   5.37,   6.47,   7.76,   9.26,  11.00,  13.00,  
        15.30,  18.00,  21.10,  26.70,  33.50,  41.90,  52.20,  64.80,  80.50, 125.00 } );

    std::vector<double> vmom3( { // e+- binning
        0.50 ,0.65 ,0.82, // 3
        1.01, 2.0, 3.0, 4.12, 5.0, 6.0, 7.1, 8.3, 9.62, 11.04,
        12.59, 14.25, 16.05, 17.98, 20.04, 22.25, 24.62, 27.25, 30.21, 35.36, 40 ,// 23
        41.9 , 45.1 ,48.5 ,52.2 ,56.1 ,60.3 ,64.8 ,69.7 ,74.9 ,80.5  ,86.5, // 11
        93.0 ,100.  ,108. ,116. ,125. ,135. ,147. ,160  ,175. ,192.,//10
        211. ,233.  ,259. ,291. ,330. ,379. ,441. ,525. ,643  ,822. ,1130. ,1800.//12
    } );
    
    std::vector<double> vmom4( {
        1.00, 1.51, 1.92,  2.97,  3.64,  4.43,  5.37,  5.90,  6.47,  7.09, 
        7.76, 8.48, 9.26, 11.00, 13.00, 15.30, 18.00, 21.10, 24.70, 28.80, 
        33.5 } );
    
    std::vector<double> vmom5( {
        1.00,    1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97 } );
    
    Axis AXmom("Momentum [GeV]", vmom);
    Axis AXrig("Rigidity [GV]" , vmom);
    
    Axis AXmom2("Momentum [GeV]", vmom4);
    Axis AXrig2("Rigidity [GV]" , vmom4);
    
    Axis AXrig3("Rigidity [GV]" , 400, 0.6, 800.0, AxisScale::kLog);
    
    Axis AXnchi("nchi", 500, -3.0,  2.0, AxisScale::kLinear);
    Axis AXnext("next", 50,  0, 50.0, AxisScale::kLinear);
    Axis AXchrg("Z", 400, 0.0, 10.0);
    
    Axis AXPROFrhb("Beta", 50, 0.95, 0.98);
   
    Axis AXPROFmass("Mass/Z [GV/c^{2}]", 600, -3.0, 3.0);
    //Axis AXPROFmass("Mass/Z [GV/c^{2}]", 1000, -10.0, 10.0);
    Hist::New("hTF_CK_mass_ISS", AXPROFmass);
    Hist::New("hTF_CK_mass_ISS_CF", AXPROFmass);
    Hist::New("hTF_CK_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hTF_CK_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hTF_CK_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hTF_CK_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hTF_CK_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hTF_CK_mass_cnt_ISS_CF", HistAxis(AXrig3));
    Hist::New("hTF_CK_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hTF_CK_mass_cnt_MC_FLUX10", HistAxis(AXrig3));
    Hist::New("hTF_CK_mass_cnt_MC_FLUX27", HistAxis(AXrig3));
    Hist::New("hTF_CK_mass_cnt_MC_FLUX10_CF", HistAxis(AXrig3));
    Hist::New("hTF_CK_mass_cnt_MC_FLUX27_CF", HistAxis(AXrig3));
    
    Hist::New("hTF_HC_expt_gISS", HistAxis(AXrig3));
    Hist::New("hTF_HC_expt_evt_gISS", HistAxis(AXrig3));
    Hist::New("hTF_HC_expt_evt_ISS", HistAxis(AXrig3));
    Hist::New("hTF_HC_expt_evt_ISS_CF", HistAxis(AXrig3));
    
    Hist::New("hTF_HC_phtrg_ISS", AXPROFmass);
    Hist::New("hTF_HC_altrg_ISS", AXPROFmass);
    Hist::New("hTF_HC_phtrg_ISS_CF", AXPROFmass);
    Hist::New("hTF_HC_altrg_ISS_CF", AXPROFmass);
    
    Hist::New("hTF_HC_rig_phtrg_ISS", AXrig3);
    Hist::New("hTF_HC_rig_altrg_ISS", AXrig3);
    Hist::New("hTF_HC_rig_phtrg_ISS_CF", AXrig3);
    Hist::New("hTF_HC_rig_altrg_ISS_CF", AXrig3);

    Hist::New("hTF_HC_mass_ISS", AXPROFmass);
    Hist::New("hTF_HC_mass_ISS_CF", AXPROFmass);
    Hist::New("hTF_HC_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hTF_HC_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hTF_HC_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hTF_HC_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hTF_HC_mass_MC_PURE_FLUX10", AXPROFmass);
    Hist::New("hTF_HC_mass_MC_PURE_FLUX27", AXPROFmass);
    Hist::New("hTF_HC_mass_MC_PURE_FLUX10_CF", AXPROFmass);
    Hist::New("hTF_HC_mass_MC_PURE_FLUX27_CF", AXPROFmass);
    
    Hist::New("hRH_HC_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hRH_HC_mass_cnt_ISS_CF", HistAxis(AXrig3));
    Hist::New("hRH_HC_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hRH_HC_mass_cnt_MC_FLUX10", HistAxis(AXrig3));
    Hist::New("hRH_HC_mass_cnt_MC_FLUX27", HistAxis(AXrig3));
    Hist::New("hRH_HC_mass_cnt_MC_FLUX10_CF", HistAxis(AXrig3));
    Hist::New("hRH_HC_mass_cnt_MC_FLUX27_CF", HistAxis(AXrig3));
    
    Hist::New("hRH_HC_mass_PI_POS_MC_FLUX10" , AXPROFmass);
    Hist::New("hRH_HC_mass_PI_POS_MC_FLUX27", AXPROFmass);
    Hist::New("hRH_HC_mass_PI_NEG_MC_FLUX10" , AXPROFmass);
    Hist::New("hRH_HC_mass_PI_NEG_MC_FLUX27", AXPROFmass);
    
    Hist::New("hRH_HC_mass_K_POS_MC_FLUX10" , AXPROFmass);
    Hist::New("hRH_HC_mass_K_POS_MC_FLUX27", AXPROFmass);
    Hist::New("hRH_HC_mass_K_NEG_MC_FLUX10" , AXPROFmass);
    Hist::New("hRH_HC_mass_K_NEG_MC_FLUX27", AXPROFmass);
    
    Hist::New("hRH_HC_mass_lchix_ISS", HistAxis(AXPROFmass, AXnchi));
    Hist::New("hRH_HC_mass_lchiy_ISS", HistAxis(AXPROFmass, AXnchi));
    Hist::New("hRH_HC_mass_lchib_ISS", HistAxis(AXPROFmass, AXnchi));
    
    Hist::New("hRH_HC_mass_nseg_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_nhit_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_nvtxx_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_nvtxy_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hRH_HC_mass_tkL1_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_tkL2_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_tkL9_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hRH_HC_mass_inn_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_out_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hRH_HC_mass_ntd1_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_ntd2_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_ntd3_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hRH_CK_mass_ISS", AXPROFmass);
    Hist::New("hRH_CK_mass_ISS_CF", AXPROFmass);
    Hist::New("hRH_CK_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hRH_CK_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hRH_CK_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hRH_CK_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hRH_HC_mass_ISS", AXPROFmass);
    Hist::New("hRH_HC_mass_ISS_CF", AXPROFmass);
    Hist::New("hRH_HC_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hRH_HC_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hRH_HC_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hRH_HC_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hRH_HC_mass_MC_PURE_FLUX10", AXPROFmass);
    Hist::New("hRH_HC_mass_MC_PURE_FLUX27", AXPROFmass);
    Hist::New("hRH_HC_mass_MC_PURE_FLUX10_CF", AXPROFmass);
    Hist::New("hRH_HC_mass_MC_PURE_FLUX27_CF", AXPROFmass);
    
    Hist::New("hRH2_HC_mass_ISS", AXPROFmass);
    Hist::New("hRH2_HC_mass_ISS_CF", AXPROFmass);
    Hist::New("hRH2_HC_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hRH2_HC_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hRH2_HC_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hRH2_HC_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hRH2_HC_mass_MC_PURE_FLUX10", AXPROFmass);
    Hist::New("hRH2_HC_mass_MC_PURE_FLUX27", AXPROFmass);
    Hist::New("hRH2_HC_mass_MC_PURE_FLUX10_CF", AXPROFmass);
    Hist::New("hRH2_HC_mass_MC_PURE_FLUX27_CF", AXPROFmass);
    
    Hist::New("hRH_HC_mass_lchix_ISS", HistAxis(AXPROFmass, AXnchi));
    Hist::New("hRH_HC_mass_lchiy_ISS", HistAxis(AXPROFmass, AXnchi));
    Hist::New("hRH_HC_mass_lchib_ISS", HistAxis(AXPROFmass, AXnchi));
    
    Hist::New("hRH_HC_mass_inn_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_out_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hRH_HC_mass_npmt_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_nhit_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hRH_HC_mass_nvtxx_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hRH_HC_mass_nvtxy_ISS", HistAxis(AXPROFmass, AXnext));
    
    Axis AXPROFxxx("", 1000, 0, 20);
    Hist::New("hRH_HC_mass_tofq_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_tof_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_tof2_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_rh_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_rh2_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    
    Hist::New("hRH_HC_mass_npe_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_chi_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    
    Hist::New("hRH_HC_mass_mij_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_exp_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_chg_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_bdr_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_trc_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hRH_HC_mass_acc_ISS", HistAxis(AXPROFmass, AXPROFxxx));

    return true;
}

bool Analyzer::process_presel() {
    // RTI
    if (CheckType(Type::ISS)) {
        if (rti->is_in_SAA) return false;
        if (rti->live_time < 0.5) return false;
        if (std::abs(rti->tk_align[0][0]) > 35.0) return false;
        if (std::abs(rti->tk_align[0][1]) > 35.0) return false;
        if (std::abs(rti->tk_align[1][0]) > 45.0) return false;
        if (std::abs(rti->tk_align[1][1]) > 45.0) return false;
    }

    // Trigger
    if ((trg->bit&2) != 2 && (trg->bit&8) != 8) return false;

    // Charge
    if (tof->Qall < 0.8 || tof->Qall > 1.3) return false;
    if (trk->QIn < 0.8 || trk->QIn > 1.3) return false;

    // ACC
    if (acc->num_cls != 0) return false;

    // TOF
    //if (tof->num_beta != 1) return false;
    if (!tof->status) return false;

    // TRK
    //if (tkInnext > 6) return false;
    //if (tkL1next > 6) return false;

    // TRD
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    
    return true;
}

bool Analyzer::process_data() {
    process_data_tof();
    process_data_rich();
    return true;
}

bool Analyzer::process_data_tof() {
    double mc_flux10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom, -1.7) : 1.0;
    double mc_flux27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        double crr_flux_pr = gTMPpr27->Eval(g4mc->prm_mom, 0, "");
        if (crr_flux_pr < 0) crr_flux_pr = 0.0;

        double crr_flux_ap = crr_flux_pr * gTMPap2pr->Eval(g4mc->prm_mom, 0, "");
        if (crr_flux_ap < 0) crr_flux_ap = 0.0;

        mc_flux10 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap);
        mc_flux27 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap);
    }
    
    if (tof->extcls_noise != 0) return false;
    if (tof->num_in_time_cls > 4) return false;
    if (trk->lay[1] == 0) return false;

    bool is_tk_ok = 
        !((trk->ext_num_hit[2] > 0 && trk->ext_num_hit[3] > 0) || 
          (trk->ext_num_hit[4] > 0 && trk->ext_num_hit[5] > 0) || 
          (trk->ext_num_hit[6] > 0 && trk->ext_num_hit[7] > 0));
    if (!is_tk_ok) return false;

    bool is_rich_ok = 
        !rich->self_status ||
        (rich->self_is_good_geom &&
        rich->self_num_stone <= 1 &&
        rich->self_num_cloud == 0 &&
        rich->self_num_tumor == 0 &&
        rich->self_num_ghost == 0 &&
        rich->self_nhit_other_inn <= 1 &&
        rich->self_nhit_other_out <= 3 &&
        rich->self_trace > 0.30 &&
        (!rich->self_stn_status || rich->self_stn_dist <= 3.4));
    if (!is_rich_ok) return false;

    bool is_trd_ok = (trd->num_vtx[0] <= 20 && trd->num_vtx[1] <= 20);
    if (!is_rich_ok) return false;
    
    // Trigger
    bool phtrg = (trg->bit&8) == 8;
    bool untrg = (trg->bit&2) == 2;
    int  wgtrg = (phtrg ? 1 : 100);
    
    if (trk->ck_status[0] && 
        std::log(trk->ck_nchi[0][0]) < 2.0 &&
        std::log(trk->ck_nchi[0][1]) < 2.0 &&
        tof->status &&
        std::log(tof->nchi_t) < 2.0 &&
        std::log(tof->nchi_c) < 2.0 &&
        phtrg
        ) {
        
        double sign = (trk->ck_rig[0] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * std::sqrt((trk->ck_rig[0] * (1.0 / tof->beta + 1.0)) * (trk->ck_rig[0] * (1.0 / tof->beta - 1.0)));
        double cfr  = CheckType(Type::ISS) ? std::abs(trk->ck_rig[0]/rti->max_IGRF) : 0.0;

        double wgtcf = (*Hist::Head("hTFcf"))()->GetBinContent( (*Hist::Head("hTFcf"))()->FindBin(std::abs(trk->ck_rig[0])/0.75) );

        if (tof->beta > 0.5 && tof->beta < 0.8) {
            Hist::Head("hTF_CK_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hTF_CK_mass_ISS_CF")->fillH1D(mass, list->weight);
            
            Hist::Head("hTF_CK_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hTF_CK_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hTF_CK_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hTF_CK_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
        }
        if (trk->ck_rig[0] > 0) Hist::Head("hTF_CK_mass_cnt_ISS")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
        if (trk->ck_rig[0] > 0 && cfr > 0.75) Hist::Head("hTF_CK_mass_cnt_ISS_CF")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);

        if (trk->ck_rig[0] > 0) Hist::Head("hTF_CK_mass_cnt_MC")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
        if (trk->ck_rig[0] > 0) Hist::Head("hTF_CK_mass_cnt_MC_FLUX10")->fillH1D(std::abs(trk->ck_rig[0]), mc_flux10 * list->weight);
        if (trk->ck_rig[0] > 0) Hist::Head("hTF_CK_mass_cnt_MC_FLUX27")->fillH1D(std::abs(trk->ck_rig[0]), mc_flux27 * list->weight);
        if (trk->ck_rig[0] > 0) Hist::Head("hTF_CK_mass_cnt_MC_FLUX10_CF")->fillH1D(std::abs(trk->ck_rig[0]), wgtcf * mc_flux10 * list->weight);
        if (trk->ck_rig[0] > 0) Hist::Head("hTF_CK_mass_cnt_MC_FLUX27_CF")->fillH1D(std::abs(trk->ck_rig[0]), wgtcf * mc_flux27 * list->weight);
    }
    
    
    if (hyc->geom_status[0] &&
        std::log(hyc->geom_nchi_x[0]) < 1.75 &&
        std::log(hyc->geom_nchi_y[0]) < 1.75 &&
        hyc->vel_status[1] &&
        std::log(hyc->vel_nchi[1]) < 2.00 &&
        hyc->mutr_status[1] && 
        std::log(hyc->mutr_nchi_x[1]) < 1.50 &&
        std::log(hyc->mutr_nchi_y[1]) < 1.50 &&
        std::log(hyc->mutr_nchi_b[1]) < 1.50
        ) {

        double sign = (hyc->mutr_top_rig[1] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * hyc->mutr_mass[1];
        double cfr  = CheckType(Type::ISS) ? std::abs(hyc->mutr_top_rig[1]/rti->max_IGRF) : 0.0;
        
        double wgtcf = (*Hist::Head("hTFcf"))()->GetBinContent( (*Hist::Head("hTFcf"))()->FindBin(std::abs(hyc->mutr_top_rig[1])/0.75) );
        
        bool pure_mc = CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3]);
            
        bool is_new_utime = false;
        if (CheckType(Type::ISS) && utime_cur != list->utime) { utime_cur = list->utime; is_new_utime = true; }
        TH1D* hexpt = (TH1D*) (*Hist::Head("hTF_HC_expt_gISS"))();
        int   hexpt_rbin = CheckType(Type::ISS) ? hexpt->FindBin(rti->max_IGRF) : 0;

        if (is_new_utime) {
            for (int bin = hexpt_rbin + 1; bin <= hexpt->GetXaxis()->GetNbins(); bin++) {
                hexpt->Fill( hexpt->GetBinCenter(bin), (CheckType(Type::ISS) ? rti->live_time : 1.0));
            }
        }
        if (CheckType(Type::ISS)) {
            for (int bin = hexpt_rbin + 1; bin <= hexpt->GetXaxis()->GetNbins(); bin++) {
                Hist::Head("hTF_HC_expt_evt_gISS")->fillH1D(hexpt->GetBinCenter(bin), rti->live_time * list->weight);
            }
        }
        
        if (hyc->mutr_top_bta[1] > 0.5 && hyc->mutr_top_bta[1] < 0.8) {
            if (phtrg) Hist::Head("hTF_HC_phtrg_ISS")->fillH1D(mass, list->weight);
            Hist::Head("hTF_HC_altrg_ISS")->fillH1D(mass, wgtrg * list->weight);
            if (cfr > 0.75 && phtrg) Hist::Head("hTF_HC_phtrg_ISS_CF")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hTF_HC_altrg_ISS_CF")->fillH1D(mass, wgtrg * list->weight);
            
            if (phtrg) Hist::Head("hTF_HC_rig_phtrg_ISS")->fillH1D(std::abs(hyc->mutr_top_rig[1]), list->weight);
            Hist::Head("hTF_HC_rig_altrg_ISS")->fillH1D(std::abs(hyc->mutr_top_rig[1]), wgtrg * list->weight);
            if (cfr > 0.75 && phtrg) Hist::Head("hTF_HC_rig_phtrg_ISS_CF")->fillH1D(std::abs(hyc->mutr_top_rig[1]), list->weight);
            if (cfr > 0.75) Hist::Head("hTF_HC_rig_altrg_ISS_CF")->fillH1D(std::abs(hyc->mutr_top_rig[1]), wgtrg * list->weight);
        }

        if (phtrg && hyc->mutr_top_bta[1] > 0.5 && hyc->mutr_top_bta[1] < 0.8) {
            if (CheckType(Type::ISS)) {
                for (int bin = hexpt_rbin + 1; bin <= hexpt->GetXaxis()->GetNbins(); bin++) {
                    Hist::Head("hTF_HC_expt_evt_ISS")->fillH1D(hexpt->GetBinCenter(bin), rti->live_time * list->weight);
                    if (cfr > 0.75) Hist::Head("hTF_HC_expt_evt_ISS_CF")->fillH1D(hexpt->GetBinCenter(bin), rti->live_time * list->weight);
                }
            }
            
            Hist::Head("hTF_HC_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hTF_HC_mass_ISS_CF")->fillH1D(mass, list->weight);
            
            Hist::Head("hTF_HC_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hTF_HC_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hTF_HC_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hTF_HC_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
            
            if (pure_mc) Hist::Head("hTF_HC_mass_MC_PURE_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            if (pure_mc) Hist::Head("hTF_HC_mass_MC_PURE_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            if (pure_mc) Hist::Head("hTF_HC_mass_MC_PURE_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            if (pure_mc) Hist::Head("hTF_HC_mass_MC_PURE_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
                
            //if (mass < -1.6) clone_tree->Fill();

            Hist::Head("hTF_HC_mass_lchix_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_x[1]), list->weight);
            Hist::Head("hTF_HC_mass_lchiy_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_y[1]), list->weight);
            Hist::Head("hTF_HC_mass_lchib_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_b[1]), list->weight);
            
            Hist::Head("hTF_HC_mass_nseg_ISS")->fillH2D(mass, trd->num_extra_seg, list->weight);
            Hist::Head("hTF_HC_mass_nhit_ISS")->fillH2D(mass, trd->num_extra_hit, list->weight);
            Hist::Head("hTF_HC_mass_nvtxx_ISS")->fillH2D(mass, trd->num_vtx[0], list->weight);
            Hist::Head("hTF_HC_mass_nvtxy_ISS")->fillH2D(mass, trd->num_vtx[1], list->weight);
            
            Hist::Head("hTF_HC_mass_tkL1_ISS")->fillH2D(mass, trk->ext_num_hit[0], list->weight);
            Hist::Head("hTF_HC_mass_tkL2_ISS")->fillH2D(mass, trk->ext_num_hit[1], list->weight);
            Hist::Head("hTF_HC_mass_tkL9_ISS")->fillH2D(mass, trk->ext_num_hit[8], list->weight);
            
            Hist::Head("hTF_HC_mass_inn_ISS")->fillH2D(mass, rich->self_nhit_other_inn, list->weight);
            Hist::Head("hTF_HC_mass_out_ISS")->fillH2D(mass, rich->self_nhit_other_out, list->weight);
            
            Hist::Head("hTF_HC_mass_ntd1_ISS")->fillH2D(mass, trd->tkLLR_num_hit, list->weight);
            Hist::Head("hTF_HC_mass_ntd2_ISS")->fillH2D(mass, trd->num_tkHit, list->weight);
            Hist::Head("hTF_HC_mass_ntd3_ISS")->fillH2D(mass, trd->tkLLR_num_hit-trd->num_tkHit, list->weight);
                
            // pion 0 ~ 0.3
            if (std::abs(mass) < 0.28) {
                Hist::Head("hTF_HC_mass_PI_POS_MC_FLUX10")->fillH1D( std::abs(mass), mc_flux10 * list->weight);
                Hist::Head("hTF_HC_mass_PI_POS_MC_FLUX27")->fillH1D( std::abs(mass), mc_flux27 * list->weight);
                Hist::Head("hTF_HC_mass_PI_NEG_MC_FLUX10")->fillH1D(-std::abs(mass), mc_flux10 * list->weight);
                Hist::Head("hTF_HC_mass_PI_NEG_MC_FLUX27")->fillH1D(-std::abs(mass), mc_flux27 * list->weight);
            }
            
            // koan -0.35 ~ -0.75
            if (mass > -0.75 && mass < -0.35) {
                Hist::Head("hTF_HC_mass_K_POS_MC_FLUX10")->fillH1D( std::abs(mass), mc_flux10 * list->weight);
                Hist::Head("hTF_HC_mass_K_POS_MC_FLUX27")->fillH1D( std::abs(mass), mc_flux27 * list->weight);
                Hist::Head("hTF_HC_mass_K_NEG_MC_FLUX10")->fillH1D(-std::abs(mass), mc_flux10 * list->weight);
                Hist::Head("hTF_HC_mass_K_NEG_MC_FLUX27")->fillH1D(-std::abs(mass), mc_flux27 * list->weight);
            }
        }
        
        if (phtrg) {
            if (hyc->mutr_top_rig[1] > 0) Hist::Head("hTF_HC_mass_cnt_ISS")->fillH1D(std::abs(hyc->mutr_top_rig[1]), list->weight);
            if (hyc->mutr_top_rig[1] > 0 && cfr > 0.75) Hist::Head("hTF_HC_mass_cnt_ISS_CF")->fillH1D(std::abs(hyc->mutr_top_rig[1]), list->weight);

            if (g4mc != nullptr) Hist::Head("hTF_HC_mass_cnt_MC")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
            if (g4mc != nullptr) Hist::Head("hTF_HC_mass_cnt_MC_FLUX10")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux10 * list->weight);
            if (g4mc != nullptr) Hist::Head("hTF_HC_mass_cnt_MC_FLUX27")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux27 * list->weight);
            if (g4mc != nullptr) Hist::Head("hTF_HC_mass_cnt_MC_FLUX10_CF")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), wgtcf * mc_flux10 * list->weight);
            if (g4mc != nullptr) Hist::Head("hTF_HC_mass_cnt_MC_FLUX27_CF")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), wgtcf * mc_flux27 * list->weight);
        }
    }


    return true;
}

bool Analyzer::process_data_rich() {
    double mc_flux10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom, -1.7) : 1.0;
    double mc_flux27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        double crr_flux_pr = gTMPpr27->Eval(g4mc->prm_mom, 0, "");
        if (crr_flux_pr < 0) crr_flux_pr = 0.0;

        double crr_flux_ap = crr_flux_pr * gTMPap2pr->Eval(g4mc->prm_mom, 0, "");
        if (crr_flux_ap < 0) crr_flux_ap = 0.0;

        mc_flux10 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap);
        mc_flux27 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap);
    }
    
    if (tof->extcls_noise != 0) return false;
    if (tof->num_in_time_cls > 4) return false;
    if (trk->lay[1] == 0) return false;
    
    // Trigger
    bool phtrg = (trg->bit&8) == 8;
    bool untrg = (trg->bit&2) == 2;
    int  wgtrg = (phtrg ? 1 : 100);

    bool is_rich_pr_ok = 
        rich->self_status &&
        rich->self_kind == 1;
        //rich->self_kind == 1 &&
        //rich->self_is_good_geom &&
        //!rich->self_is_bad_tile &&
        //rich->self_num_stone <= 1 &&
        //rich->self_num_cloud == 1 &&
        //rich->self_num_tumor == 0 &&
        //rich->self_num_ghost == 0 &&
        ////rich->self_nhit_other_inn <= 1 &&
        ////rich->self_nhit_other_out <= 3 &&
        //rich->self_trace > 0.30 &&
        //rich->self_cld_status &&
        //(!rich->self_stn_status || rich->self_stn_dist <= 3.4);
    
    bool is_official_rich_pr_ok =
        rich->status &&
        rich->kind == 1;

    if (trk->ck_status[0] && 
        std::log(trk->ck_nchi[0][0]) < 2.0 &&
        std::log(trk->ck_nchi[0][1]) < 2.0 &&
        is_official_rich_pr_ok &&
        phtrg
        ) {
        
        double sign = (trk->ck_rig[0] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * std::sqrt((trk->ck_rig[0] * (1.0 / rich->beta + 1.0)) * (trk->ck_rig[0] * (1.0 / rich->beta - 1.0)));
        double cfr  = CheckType(Type::ISS) ? std::abs(trk->ck_rig[0]/rti->max_IGRF) : 0.0;

        double wgtcf = (*Hist::Head("hRHcf"))()->GetBinContent( (*Hist::Head("hRHcf"))()->FindBin(std::abs(trk->ck_rig[0])/0.75) );

        if (rich->beta > 0.96 && rich->beta < 0.98) {
            Hist::Head("hRH_CK_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hRH_CK_mass_ISS_CF")->fillH1D(mass, list->weight);
            
            Hist::Head("hRH_CK_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hRH_CK_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hRH_CK_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hRH_CK_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
        }
    }
    
    if (hyc->geom_status[0] &&
        std::log(hyc->geom_nchi_x[0]) < 1.75 &&
        std::log(hyc->geom_nchi_y[0]) < 1.75 &&
        hyc->vel_status[2] &&
        std::log(hyc->vel_nchi[2]) < 2.00 &&
        hyc->mutr_status[2] && 
        std::log(hyc->mutr_nchi_x[2]) < 1.50 &&
        std::log(hyc->mutr_nchi_y[2]) < 1.50 &&
        std::log(hyc->mutr_nchi_b[2]) < 1.50 &&
        is_rich_pr_ok
        ) {
        double sign = (hyc->mutr_top_rig[2] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * hyc->mutr_mass[2];
        double cfr  = CheckType(Type::ISS) ? std::abs(hyc->mutr_top_rig[2]/rti->max_IGRF) : 0.0;
        
        double wgtcf = (*Hist::Head("hRHcf"))()->GetBinContent( (*Hist::Head("hRHcf"))()->FindBin(std::abs(hyc->mutr_top_rig[2])/0.75) );
        
        bool pure_mc = CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3]);
            
        bool is_new_utime = false;
        if (CheckType(Type::ISS) && utime_cur != list->utime) { utime_cur = list->utime; is_new_utime = true; }
        TH1D* hexpt = (TH1D*) (*Hist::Head("hRH_HC_expt_gISS"))();
        int   hexpt_rbin = CheckType(Type::ISS) ? hexpt->FindBin(rti->max_IGRF) : 0;
        
        //if (phtrg && hyc->mutr_top_bta[2] > 0.96 && hyc->mutr_top_bta[2] < 0.98) {
        if (phtrg && rich->self_cld_cbta > 0.96 && rich->self_cld_cbta < 0.98) {
            //if (trd->num_vtx[0] <= 4 && 
            //    trd->num_vtx[1] <= 4 &&
            //    rich->self_nhit_other_inn <= 1 &&
            //    rich->self_nhit_other_out <= 3 &&
            //    rich->self_cld_misjudge < 0.1 &&
            //    rich->self_cld_expnpe > 0.1) {
                Hist::Head("hRH_HC_mass_ISS")->fillH1D(mass, list->weight);
                if (cfr > 0.75) Hist::Head("hRH_HC_mass_ISS_CF")->fillH1D(mass, list->weight);

                Hist::Head("hRH_HC_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
                Hist::Head("hRH_HC_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
                Hist::Head("hRH_HC_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
                Hist::Head("hRH_HC_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
                
                if (pure_mc) Hist::Head("hRH_HC_mass_MC_PURE_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
                if (pure_mc) Hist::Head("hRH_HC_mass_MC_PURE_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
                if (pure_mc) Hist::Head("hRH_HC_mass_MC_PURE_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
                if (pure_mc) Hist::Head("hRH_HC_mass_MC_PURE_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
            
            //}
        
            if (is_official_rich_pr_ok && trk->ck_status[0] && std::abs(rich->self_cld_cbta-rich->beta) < 0.004 && rich->cstcb < 0.004 && rich->npmt >= 3) {
                Hist::Head("hRH2_HC_mass_ISS")->fillH1D(mass, list->weight);
                if (cfr > 0.75) Hist::Head("hRH2_HC_mass_ISS_CF")->fillH1D(mass, list->weight);

                Hist::Head("hRH2_HC_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
                Hist::Head("hRH2_HC_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
                Hist::Head("hRH2_HC_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
                Hist::Head("hRH2_HC_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
                
                if (pure_mc) Hist::Head("hRH2_HC_mass_MC_PURE_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
                if (pure_mc) Hist::Head("hRH2_HC_mass_MC_PURE_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
                if (pure_mc) Hist::Head("hRH2_HC_mass_MC_PURE_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
                if (pure_mc) Hist::Head("hRH2_HC_mass_MC_PURE_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
                
                if (mass < 0) clone_tree->Fill();
            }
            
            Hist::Head("hRH_HC_mass_lchix_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_x[2]), list->weight);
            Hist::Head("hRH_HC_mass_lchiy_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_y[2]), list->weight);
            Hist::Head("hRH_HC_mass_lchib_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_b[2]), list->weight);
            
            Hist::Head("hRH_HC_mass_inn_ISS")->fillH2D(mass, rich->self_nhit_other_inn, list->weight);
            Hist::Head("hRH_HC_mass_out_ISS")->fillH2D(mass, rich->self_nhit_other_out, list->weight);
            
            Hist::Head("hRH_HC_mass_nvtxx_ISS")->fillH2D(mass, trd->num_vtx[0], list->weight);
            Hist::Head("hRH_HC_mass_nvtxy_ISS")->fillH2D(mass, trd->num_vtx[1], list->weight);
            
            Hist::Head("hRH_HC_mass_npmt_ISS")->fillH2D(mass, rich->self_cld_nhit, list->weight);
            Hist::Head("hRH_HC_mass_nhit_ISS")->fillH2D(mass, rich->self_cld_npmt, list->weight);
            
            Hist::Head("hRH_HC_mass_tof_ISS")->fillH2D(mass, hyc->vel_top_bta[1], list->weight);
            Hist::Head("hRH_HC_mass_tof2_ISS")->fillH2D(mass, tof->beta, list->weight);
            Hist::Head("hRH_HC_mass_rh_ISS")->fillH2D(mass, hyc->vel_top_bta[2], list->weight);
            if (is_official_rich_pr_ok) Hist::Head("hRH_HC_mass_rh2_ISS")->fillH2D(mass, rich->beta, list->weight);
            
            Hist::Head("hRH_HC_mass_tofq_ISS")->fillH2D(mass, (tof->Q[0]+tof->Q[1])/(tof->Q[2]+tof->Q[3]), list->weight);
            
            Hist::Head("hRH_HC_mass_npe_ISS")->fillH2D(mass, rich->self_cld_npe, list->weight);
            Hist::Head("hRH_HC_mass_chi_ISS")->fillH2D(mass, rich->self_cld_nchi, list->weight);
            
            Hist::Head("hRH_HC_mass_mij_ISS")->fillH2D(mass, rich->self_cld_misjudge, list->weight);
            Hist::Head("hRH_HC_mass_exp_ISS")->fillH2D(mass, rich->self_cld_expnpe, list->weight);
            if (rich->self_cld_expnpe > 0.1) Hist::Head("hRH_HC_mass_chg_ISS")->fillH2D(mass, std::sqrt(rich->self_cld_npe / rich->self_cld_expnpe), list->weight);
            Hist::Head("hRH_HC_mass_bdr_ISS")->fillH2D(mass, rich->self_cld_border, list->weight);
            Hist::Head("hRH_HC_mass_trc_ISS")->fillH2D(mass, rich->self_cld_trace, list->weight);
            Hist::Head("hRH_HC_mass_acc_ISS")->fillH2D(mass, rich->self_cld_accuracy, list->weight);
        }
    }

    return true;
}


#endif // __Analyzer_C__
