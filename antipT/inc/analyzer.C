#ifndef __Analyzer_C__
#define __Analyzer_C__

#include "analyzer.h"

using namespace MGROOT;

void Analyzer::set_environment() {
    std::cout << Format("\n====  Set Environment ====\n");
    LOG(INFO) << Format("\n====  Set Environment ====\n");
    
    TrSys::PhysEnv::ReadMagAMS("/eos/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    TrSys::PhysEnv::ReadMatAMS("/eos/user/h/hchou/ExternalLibs/DB/material");
    
    if (!TrSys::PhysEnv::IMagStatus()) TrSys::PhysEnv::ReadMagAMS("/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    if (!TrSys::PhysEnv::IMatStatus()) TrSys::PhysEnv::ReadMatAMS("/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    
    if (!TrSys::PhysEnv::IMagStatus()) TrSys::PhysEnv::ReadMagAMS("/afs/cern.ch/work/h/hchou/public/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    if (!TrSys::PhysEnv::IMatStatus()) TrSys::PhysEnv::ReadMatAMS("/afs/cern.ch/work/h/hchou/public/ExternalLibs/DB/material");
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

bool Analyzer::build_hist() {
    std::vector<double> vmom( {
        0.50,    0.80,
        1.00,    1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
        3.29,    3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
        8.48,    9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00
    } );
    
    Axis AXmom("Momentum [GeV]", vmom);
    Axis AXrig("Rigidity [GV]" , vmom);
    Axis AXtme("Unix Time", 50, 1305850000, 1546410000);
   
    Axis AXsqrm("Mass^{2}", 1000, -4.0, 16.0);
    Hist::New("hLP_sqrm", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN_sqrm", HistAxis(AXrig, AXsqrm));
    
    Hist::New("hLP_sqrm2", HistAxis(AXrig, AXsqrm));
    Hist::New("hLP_sqrm3", HistAxis(AXrig, AXsqrm));
    Hist::New("hLP_sqrm4", HistAxis(AXrig, AXsqrm));
    Hist::New("hLP_sqrm5", HistAxis(AXrig, AXsqrm));
    
    Hist::New("hLN_sqrm2", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN_sqrm3", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN_sqrm4", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN_sqrm5", HistAxis(AXrig, AXsqrm));
    
    Hist::New("hLP2_sqrm", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN2_sqrm", HistAxis(AXrig, AXsqrm));
    
    Hist::New("hLP2_sqrm2", HistAxis(AXrig, AXsqrm));
    Hist::New("hLP2_sqrm3", HistAxis(AXrig, AXsqrm));
    Hist::New("hLP2_sqrm4", HistAxis(AXrig, AXsqrm));
    Hist::New("hLP2_sqrm5", HistAxis(AXrig, AXsqrm));
    
    Hist::New("hLN2_sqrm2", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN2_sqrm3", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN2_sqrm4", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN2_sqrm5", HistAxis(AXrig, AXsqrm));
    
    Hist::New("hLP_sqrm_pr", HistAxis(AXrig, AXsqrm));
    Hist::New("hLP_sqrm_pr2", HistAxis(AXrig, AXsqrm));
    Hist::New("hLN_sqrm_el", HistAxis(AXrig, AXsqrm));

    //for (int it = 1; it <= AXrig.nbin(); ++it) {
    //    Hist::New(Format("hLP%02d_sqrm", it).c_str(), AXtme, AXsqrm);
    //    Hist::New(Format("hLN%02d_sqrm", it).c_str(), AXtme, AXsqrm);
    //}
    
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

    //// Charge
    if (tof->Qall < 0.8 || tof->Qall > 1.3) return false;
    if (trk->QIn < 0.8 || trk->QIn > 1.3) return false;

    //// ACC
    if (acc->num_cls != 0) return false;

    //// TOF
    if (!tof->status) return false;
    if (tof->nchi_t > 10.) return false;
    if (tof->nchi_c > 10.) return false;

    // TRD
    if (!trd->status || trd->num_cls > 100) return false;
    //if (trd->num_hit > 50 || trd->num_Hhit > 50) return false;
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 10) return false;

    // HYC
    if (!hyc->geom_status[0]) return false;
    if (!hyc->phys_status[0][0]) return false;

    //// cutoff
    if (CheckType(Type::ISS)) {
        //if (std::abs(hyc->geom_top_rig[0]) < 0.8 * rti->max_IGRF) return false;
        if (std::abs(hyc->phys_top_rig[0][0]) < 0.8 * rti->max_IGRF) return false;
    }

    return true;
}

bool Analyzer::process_data() {
    process_data_l();
    //process_data_m();
    return true;
}

bool Analyzer::process_data_l() {
    if (!hyc->geom_status[0]) return false;
    if (!hyc->phys_status[0][0]) return false;
    if (!hyc->mutr_status[0]) return false;

    if (std::log(hyc->geom_nchi_x[0]) > 2.0) return false; 
    if (std::log(hyc->geom_nchi_y[0]) > 2.0) return false; 

    if (std::log(hyc->mutr_nchi_x[0]) > 2.0) return false; 
    if (std::log(hyc->mutr_nchi_y[0]) > 2.0) return false; 
    if (std::log(hyc->mutr_nchi_b[0]) > 2.0) return false; 

    if (std::log(hyc->phys_nchi_x[0][0]) > 2.0) return false;
    if (std::log(hyc->phys_nchi_y[0][0]) > 2.0) return false;
    if (std::log(hyc->phys_nchi_b[0][0]) > 2.0) return false;

    if (!rich->self_status) return false;
    if (rich->self_kind != 1) return false;
    if (!rich->self_is_good_geom) return false;
    //if (rich->self_border < 0.35) return false;
    //if (rich->self_trace  < 0.10) return false;
    if (rich->self_num_stone  > 1) return false;
    if (rich->self_num_tumor != 0) return false;
    if (rich->self_num_ghost != 0) return false;

    if (rich->self_stn_status && rich->self_stn_nchi > 3.5) return false;
    if (rich->self_stn_status && rich->self_stn_chic > 3.0) return false;
    if (rich->self_stn_status && rich->self_stn_dist > 3.4) return false;
    
    if (rich->self_nhit_other_inn != 0) return false;
    if (rich->self_nhit_other_out >  2) return false;

    double rig     = hyc->phys_top_rig[0][0];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cfr     = (cf_rig / abs_rig);
    
    double sqrm  = hyc->mutr_sqrm[0];
    
    //if (CheckType(Type::ISS) && cfr > 1.2) return false;

    //if (signr > 0) Hist::Head("hLP_trd")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
    //if (signr < 0) Hist::Head("hLN_trd")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
 
    if (signr > 0 && trd->tdLLR_ep > 0.75 && !rich->self_cld_status) {
        Hist::Head("hLP_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
        Hist::Head("hLP_sqrm_pr2")->fillH2D(abs_rig, sqrm, list->weight);
    }
    if (signr < 0 && trd->tdLLR_ep < 0.6 && rich->self_cld_status) {
        Hist::Head("hLN_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    }

    if (trd->tdLLR_ep < 0.75) return false;
    if (rich->self_cld_status) return false;

    if (signr > 0) Hist::Head("hLP_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hLN_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
    return true;
}

bool Analyzer::process_data_m() {
    /*
    if (!rich->self_status || rich->self_kind != 1) return false;
    if (!rich->self_is_good_geom) return false;
    if (rich->self_is_bad_tile) return false;

    double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    
    if (!trk->ck_status[0]) return false;
    double rig     = trk->ck_rig[0];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    double cfr     = (cf_rig / abs_rig);

    double nchix = std::log(trk->ck_nchi[0][0]);
    double nchiy = std::log(trk->ck_nchi[0][1]);

    if (!rich->status || rich->kind != 1) return false;
    bool is_good = (rich->npmt >= 3 && rich->prob > 0.01 && rich->cstcq < 5.0 && rich->eft_colPE > 0.4 && rich->num_clsZ1 >= 1 && rich->num_clsZ1 <= 3);
    if (!is_good) return false;

    double beta  = rich->beta;
    double mass2 = rig * rig * (1.0 / beta / beta - 1.0);
    
    if (signr > 0) Hist::Head("hMP_nchix")->fillH2D(abs_rig, nchix, list->weight);
    if (signr < 0) Hist::Head("hMN_nchix")->fillH2D(abs_rig, nchix, list->weight);
    
    if (signr > 0) Hist::Head("hMP_nchiy")->fillH2D(abs_rig, nchiy, list->weight);
    if (signr < 0) Hist::Head("hMN_nchiy")->fillH2D(abs_rig, nchiy, list->weight);
    
    if (nchix > 2.0 || nchiy > 2.0) return false;
    
    if (signr > 0) Hist::Head("hMP_cf")->fillH2D(abs_rig, cf_rig, list->weight);
    if (signr < 0) Hist::Head("hMN_cf")->fillH2D(abs_rig, cf_rig, list->weight);

    if (signr > 0) Hist::Head("hMP_cfr")->fillH2D(cf_rig, cfr, list->weight);
    if (signr < 0) Hist::Head("hMN_cfr")->fillH2D(cf_rig, cfr, list->weight);
    
    if (CheckType(Type::ISS) && cfr > 1.2) return false;
    
    if (signr > 0) Hist::Head("hMP_trd")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
    if (signr < 0) Hist::Head("hMN_trd")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);

    if (signr > 0 && trd->tdLLR_ep > 0.8) {
        Hist::Head("hMP_m2_pr")->fillH2D(abs_rig, mass2, list->weight);
    }
    if (signr < 0 && trd->tdLLR_ep < 0.6) {
        Hist::Head("hMN_m2_el")->fillH2D(abs_rig, mass2, list->weight);
    }

    if (trd->tdLLR_ep < 0.75) return false;

    if (signr > 0) Hist::Head("hMP_cnt")->fillH1D(abs_rig, list->weight);
    if (signr < 0) Hist::Head("hMN_cnt")->fillH1D(abs_rig, list->weight);
    
    if (signr > 0) Hist::Head("hMP_cnt_T")->fillH2D(list->utime, abs_rig, list->weight);
    if (signr < 0) Hist::Head("hMN_cnt_T")->fillH2D(list->utime, abs_rig, list->weight);

    if (signr > 0) Hist::Head("hMP_m2")->fillH2D(abs_rig, mass2, list->weight);
    if (signr < 0) Hist::Head("hMN_m2")->fillH2D(abs_rig, mass2, list->weight);
    
    if (signr > 0 && abs_rig > 50.0) Hist::Head("hMP_beta_T")->fillH2D(list->utime, beta, list->weight);
   */ 
    return true;
}

#endif // __Analyzer_C__
