#ifndef __Selector_C__
#define __Selector_C__

#include "selector.h"

using namespace MGROOT;

void Selector::set_environment() {
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


void Selector::process_events() {
    if (mdst == nullptr || mdst->GetEntries() == 0) return;
    if (file == nullptr || tree == nullptr) return;
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

        if (!process_init()) continue;
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

bool Selector::process_init() {
    varsFLX.init();
    varsLTF.init();
    varsLRH.init();
    varsIIN.init();
    varsIEX.init();
    varsHEX.init();
    varsHFS.init();

    return true;
}

bool Selector::process_presel() {
    // RTI
    if (CheckType(Type::ISS)) {
        if (!rti->good) return false;
        if (rti->is_in_SAA) return false;
        if (rti->zenith > 40.0) return false;
        if (rti->livetime < 0.5) return false;
        if (std::abs(rti->tk_align[0][0]) > 35.0) return false;
        if (std::abs(rti->tk_align[0][1]) > 35.0) return false;
        if (std::abs(rti->tk_align[1][0]) > 45.0) return false;
        if (std::abs(rti->tk_align[1][1]) > 45.0) return false;
    }
    
    // Trigger
    if ((trg->bit&2) != 2 && (trg->bit&8) != 8) return false;
    
    // Charge
    if (tof->Qall < 0.8 || tof->Qall > 1.4) return false;
    if (trk->QIn < 0.8 || trk->QIn > 1.4) return false;
   
    if (hyc->tfQup < 0.75 || hyc->tfQup > 2.0) return false;
    if (hyc->tfQlw < 0.75 || hyc->tfQlw > 2.0) return false;

    // TOF
    int tofu_max_num_extcls = std::max(tof->num_extcls[0], tof->num_extcls[1]);
    int tofl_max_num_extcls = std::max(tof->num_extcls[2], tof->num_extcls[3]);
    if ((tof->num_extcls[0] > 0 && tof->num_extcls[1] > 0) || tofu_max_num_extcls >= 2) return false;
    if ((tof->num_extcls[2] > 0 && tof->num_extcls[3] > 0) || tofl_max_num_extcls >= 2) return false;

    // Tracker
    int tkSL1_max_num_exthit = std::max(trk->ext_num_hit[2], trk->ext_num_hit[3]);
    int tkSL2_max_num_exthit = std::max(trk->ext_num_hit[4], trk->ext_num_hit[5]);
    int tkSL3_max_num_exthit = std::max(trk->ext_num_hit[6], trk->ext_num_hit[7]);
    if (trk->ext_num_hit[2] > 0 && trk->ext_num_hit[3] > 0 && tkSL1_max_num_exthit >= 6) return false;
    if (trk->ext_num_hit[4] > 0 && trk->ext_num_hit[5] > 0 && tkSL2_max_num_exthit >= 6) return false;
    if (trk->ext_num_hit[6] > 0 && trk->ext_num_hit[7] > 0 && tkSL3_max_num_exthit >= 6) return false;
    if (trk->ext_num_hit[1] >= 6) return false;
    if (trk->ext_num_hit[0] >= 6) return false;
    if (trk->ext_num_hit[8] >= 6) return false;

    std::vector<double> sn10;
    for (int il = 1; il <= 7; ++il) {
        if (trk->lay[il] == 0) continue;
        sn10.push_back(trk->sn10[il]);
    }
    double min_sn10 = *std::min_element(sn10.begin(), sn10.end());
    if (min_sn10 < 0.33) return false;

    if (trk->lay[0] != 0 && trk->sn10[0] < 0.65) return false;
    if (trk->lay[8] != 0 && trk->sn10[8] < 0.65) return false;

    // Others
    bool is_part   = ((rich->self_num_stone - rich->self_stn_status) == 0 && rich->self_stn_dist < 3.4);
    bool is_ring   = ((rich->self_num_cloud - rich->self_cld_status) == 0);
    bool is_shower = ((ecal->num_shower - ecal->status) == 0);
   
    if (!is_part  ) return false;
    if (!is_ring  ) return false;
    if (!is_shower) return false;

    return true;
}

bool Selector::process_data() {
    bool pass_flx = process_data_flx();
    bool pass_ltf = process_data_ltf();
    bool pass_lrh = process_data_lrh();
    bool pass_iin = process_data_iin();
    bool pass_iex = process_data_iex();
    bool pass_hex = process_data_hex();
    bool pass_hfs = process_data_hfs();

    if (pass_flx) varsFLX.fill();
    if (pass_ltf) varsLTF.fill();
    if (pass_lrh) varsLRH.fill();
    if (pass_iin) varsIIN.fill();
    if (pass_iex) varsIEX.fill();
    if (pass_hex) varsHEX.fill();
    if (pass_hfs) varsHFS.fill();

    return true;
}

bool Selector::process_data_flx() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // TRD
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    if (!trd->tkLLR_status || trd->tkLLR_num_hit < 8 || trd->num_tkHit < 6) return false;

    // Tracker
    if (!hyc->geom_status[0]) return false;
    //if (!hyc->geom_status[1]) return false;
    //if (!hyc->geom_status[2]) return false;
    //if (!hyc->geom_status[3]) return false;

    if (!trk->ck_status[0]) return false;
    //if (!trk->ck_status[1]) return false;
    //if (!trk->ck_status[2]) return false;
    //if (!trk->ck_status[3]) return false;
    
    if (std::log(hyc->geom_nchi_x[0]) > 1.75) return false;
    if (std::log(hyc->geom_nchi_y[0]) > 1.75) return false;
    //if (std::log(hyc->geom_nchi_x[1]) > 1.75) return false;
    //if (std::log(hyc->geom_nchi_y[1]) > 1.75) return false;
    //if (std::log(hyc->geom_nchi_x[2]) > 1.75) return false;
    //if (std::log(hyc->geom_nchi_y[2]) > 1.75) return false;
        
    if (std::log(trk->ck_nchi[0][0]) > 2.00) return false;
    if (std::log(trk->ck_nchi[0][1]) > 2.00) return false;
    //if (std::log(trk->ck_nchi[1][0]) > 2.00) return false;
    //if (std::log(trk->ck_nchi[1][1]) > 2.00) return false;
    //if (std::log(trk->ck_nchi[2][0]) > 2.00) return false;
    //if (std::log(trk->ck_nchi[2][1]) > 2.00) return false;

    varsFLX.init();
    varsFLX.run = list->run;
    varsFLX.evt = list->event;
    varsFLX.ut  = list->utime;
    varsFLX.wgt = list->weight;

    varsFLX.mc     = (g4mc != nullptr);
    varsFLX.mc_rig = mc_rig;
    varsFLX.mc_w10 = mc_w10;
    varsFLX.mc_w27 = mc_w27;

    varsFLX.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsFLX.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsFLX.lv    = CheckType(Type::ISS) ? rti->livetime : 1;

    varsFLX.trg = ((trg->bit&8) == 8);

    varsFLX.tdllr = trd->tdLLR_ep;
    varsFLX.tkllr = trd->tkLLR_ep;
    
    varsFLX.tkL1 = trk->lay[0];
    varsFLX.tkL9 = trk->lay[8];
    varsFLX.tkL2 = trk->lay[1];

    if (hyc->geom_status[3]) {
        varsFLX.sign  = (hyc->geom_top_rig[3] > 0.0) ? 1 : -1;
        varsFLX.rig   = hyc->geom_top_rig[3];
        varsFLX.rigc  = hyc->geom_top_crr_rig[3];
        varsFLX.lxfs  = std::log(hyc->geom_nchi_x[3]);
        varsFLX.lyfs  = std::log(hyc->geom_nchi_y[3]);
    }
    varsFLX.rigIn = hyc->geom_status[0] ? hyc->geom_top_rig[0] : 0.0;
    varsFLX.rigL1 = hyc->geom_status[1] ? hyc->geom_top_rig[1] : 0.0;
    varsFLX.rigL9 = hyc->geom_status[2] ? hyc->geom_top_rig[2] : 0.0;

    if (trk->ck_status[3]) {
        varsFLX.signC  = (trk->ck_rig[3] > 0.0) ? 1 : -1;
        varsFLX.rigC   = trk->ck_rig[3];
        varsFLX.rigcC  = trk->ck_crr_rig[3];
        varsFLX.lxfsC  = std::log(trk->ck_nchi[3][0]);
        varsFLX.lyfsC  = std::log(trk->ck_nchi[3][1]);
    }
    varsFLX.rigInC = trk->ck_status[0] ? trk->ck_rig[0] : 0.0;
    varsFLX.rigL1C = trk->ck_status[1] ? trk->ck_rig[1] : 0.0;
    varsFLX.rigL9C = trk->ck_status[2] ? trk->ck_rig[2] : 0.0;

    if (hyc->vel_status[1] && hyc->mutr_status[1] && hyc->phys_status[0][1]) {
        varsFLX.tf = true;
        varsFLX.tf_sqrm = hyc->mutr_sqrm[1];
        varsFLX.tf_bta = hyc->phys_top_bta[0][1];
        varsFLX.tf_rig = hyc->phys_top_rig[0][1];
        
        std::array<double, 3> tf_lb({
            std::log(hyc->vel_nchi[1]), 
            std::log(hyc->mutr_nchi_b[1]), 
            std::log(hyc->phys_nchi_b[0][1]) 
        });
        std::array<double, 2> tf_lx({
            std::log(hyc->mutr_nchi_x[1]), 
            std::log(hyc->phys_nchi_x[0][1]) 
        });
        std::array<double, 2> tf_ly({
            std::log(hyc->mutr_nchi_y[1]), 
            std::log(hyc->phys_nchi_y[0][1]) 
        });

        varsFLX.tf_lb = *std::max_element(tf_lb.begin(), tf_lb.end());
        varsFLX.tf_lx = *std::max_element(tf_lx.begin(), tf_lx.end());
        varsFLX.tf_ly = *std::max_element(tf_ly.begin(), tf_ly.end());
    }
    
    if (hyc->vel_status[2] && hyc->mutr_status[2] && hyc->phys_status[0][2]) {
        varsFLX.rh = true;
        varsFLX.rh_sqrm = hyc->mutr_sqrm[2];
        varsFLX.rh_bta = hyc->phys_top_bta[0][2];
        varsFLX.rh_rig = hyc->phys_top_rig[0][2];
        
        std::array<double, 3> rh_lb({
            std::log(hyc->vel_nchi[2]), 
            std::log(hyc->mutr_nchi_b[2]), 
            std::log(hyc->phys_nchi_b[0][2]) 
        });
        std::array<double, 2> rh_lx({
            std::log(hyc->mutr_nchi_x[2]), 
            std::log(hyc->phys_nchi_x[0][2]) 
        });
        std::array<double, 2> rh_ly({
            std::log(hyc->mutr_nchi_y[2]), 
            std::log(hyc->phys_nchi_y[0][2]) 
        });

        varsFLX.rh_lb = *std::max_element(rh_lb.begin(), rh_lb.end());
        varsFLX.rh_lx = *std::max_element(rh_lx.begin(), rh_lx.end());
        varsFLX.rh_ly = *std::max_element(rh_ly.begin(), rh_ly.end());
    }

    if (ecal->status) {
        varsFLX.ecal = ecal->status;
        varsFLX.engE = ecal->engE;
        varsFLX.engD = ecal->engD;
        varsFLX.engH = ecal->hadron_eng;
        varsFLX.apxL = ecal->hadron_apex;
        varsFLX.mvaBDT = ecal->mvaBDT;
    }

    return true;
}


bool Selector::process_data_ltf() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // TRD
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;

    // TOF
    if (tof->extcls_noise != 0) return false;
    if (tof->num_in_time_cls > 4) return false;

    // RICH
    if (!rich->self_status) return false;
    if (rich->self_kind != 1) return false;
    if (!rich->self_is_good_geom) return false;
    if (rich->self_num_stone > 1) return false;
    if (rich->self_num_cloud > 1) return false;
    if (rich->self_num_tumor != 0) return false;
    if (rich->self_num_ghost != 0) return false;
    if (rich->self_trace < 0.30) return false;
   
    if (rich->self_stn_status && rich->self_stn_dist > 3.4) return false;
    if (rich->self_cld_status && rich->self_cld_misjudge > 0.1) return false; 
   
    if (rich->self_expnpe < 0.01) return false;
    if (rich->self_dist < 0.25) return false;

    // TRD extra hit
    int nvtx_xl = (trd->num_vtx[0][0] + trd->num_vtx[1][0]);
    int nvtx_xu = (trd->num_vtx[2][0] + trd->num_vtx[3][0]);
    int nvtx_yl = (trd->num_vtx[0][1] + trd->num_vtx[1][1]);
    int nvtx_yu = (trd->num_vtx[2][1] + trd->num_vtx[3][1]);
    
    // Chi-squared Cut
    const double vel_nchi_cut = 2.00;
    const double trk_nchi_cut = 1.75;
    
    // Velocity
    if (!hyc->vel_status[1]) return false;
    if (std::log(hyc->vel_nchi[1]) > vel_nchi_cut) return false;
   
    // Geometry
    if (!hyc->geom_status[0]) return false;
    if (std::log(hyc->geom_nchi_x[0]) > trk_nchi_cut) return false;
    if (std::log(hyc->geom_nchi_y[0]) > trk_nchi_cut) return false;

    // Mass
    if (!hyc->mutr_status[1]) return false;
    if (std::log(hyc->mutr_nchi_x[1]) > trk_nchi_cut) return false;
    if (std::log(hyc->mutr_nchi_y[1]) > trk_nchi_cut) return false;
    if (std::log(hyc->mutr_nchi_b[1]) > vel_nchi_cut) return false;

    // Physics
    if (!hyc->phys_status[0][1]) return false;
    if (std::log(hyc->phys_nchi_x[0][1]) > trk_nchi_cut) return false;
    if (std::log(hyc->phys_nchi_y[0][1]) > trk_nchi_cut) return false;
    if (std::log(hyc->phys_nchi_b[0][1]) > vel_nchi_cut) return false;
   
    // Choutko track
    if (!trk->ck_status[0]) return false;
    if (std::log(trk->ck_nchi[0][0]) > 3.0) return false;
    if (std::log(trk->ck_nchi[0][1]) > 3.0) return false;
    
    // Variables
    double trdllr  = trd->tdLLR_ep;
    double sqrm    = hyc->mutr_sqrm[1];

    double bta     = hyc->phys_top_bta[0][1];
    double rig     = hyc->phys_top_rig[0][1];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    // Cut on chi-square
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);
   
    bool is_rich_pr = rich->self_cld_status &&
                      hyc->vel_status[2] &&
                      std::log(hyc->vel_nchi[2]) < vel_nchi_cut &&
                      hyc->mutr_status[2] &&
                      std::log(hyc->mutr_nchi_x[2]) < trk_nchi_cut &&
                      std::log(hyc->mutr_nchi_y[2]) < trk_nchi_cut &&
                      std::log(hyc->mutr_nchi_b[2]) < vel_nchi_cut &&
                      hyc->phys_status[0][2] &&
                      std::log(hyc->phys_nchi_x[0][2]) < trk_nchi_cut &&
                      std::log(hyc->phys_nchi_y[0][2]) < trk_nchi_cut &&
                      std::log(hyc->phys_nchi_b[0][2]) < vel_nchi_cut;
    
    bool is_rich_el = rich->self_cld_status &&
                      hyc->mutr_status[2] &&
                      std::log(hyc->mutr_nchi_x[2]) < trk_nchi_cut &&
                      std::log(hyc->mutr_nchi_y[2]) < trk_nchi_cut &&
                      std::log(hyc->mutr_nchi_b[2]) < vel_nchi_cut;

    varsLTF.init();
    varsLTF.run = list->run;
    varsLTF.evt = list->event;
    varsLTF.ut  = list->utime;
    varsLTF.wgt = list->weight;

    varsLTF.mc     = (g4mc != nullptr);
    varsLTF.mc_rig = mc_rig;
    varsLTF.mc_w10 = mc_w10;
    varsLTF.mc_w27 = mc_w27;

    varsLTF.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsLTF.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsLTF.lv    = CheckType(Type::ISS) ? rti->livetime : 1;

    varsLTF.sign    = signr;
    varsLTF.rig     = rig;
    varsLTF.bta     = bta;
    varsLTF.sqrm    = sqrm;
    varsLTF.llr     = trdllr;
    
    varsLTF.lxin    = lchix_in;
    varsLTF.lyin    = lchiy_in;

    varsLTF.rich    = rich->self_cld_status;
    varsLTF.rich_pr = is_rich_pr;
    varsLTF.rich_el = is_rich_el;

    varsLTF.nhinn = rich->self_nhit_other_inn;
    varsLTF.nhout = rich->self_nhit_other_out;
    
    varsLTF.nvtxx[0] = nvtx_xl;
    varsLTF.nvtxx[1] = nvtx_xu;
    varsLTF.nvtxy[0] = nvtx_yl;
    varsLTF.nvtxy[1] = nvtx_yu;

    return true;
}

bool Selector::process_data_lrh() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // TRD
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    // TOF
    if (tof->extcls_noise != 0) return false;
    if (tof->num_in_time_cls > 4) return false;

    // RICH
    if (!rich->self_status) return false;
    if (rich->self_kind != 1) return false;
    if (!rich->self_is_good_geom) return false;
    if (rich->self_is_bad_tile) return false;
    if (rich->self_num_stone > 1) return false;
    if (rich->self_num_cloud > 1) return false;
    if (rich->self_num_tumor != 0) return false;
    if (rich->self_num_ghost != 0) return false;
    if (rich->self_nhit_other_inn > 1) return false;
    if (rich->self_nhit_other_out > 3) return false;
    if (rich->self_trace < 0.30) return false;
    
    if (rich->self_stn_status && rich->self_stn_dist > 3.4) return false;
    if (rich->self_cld_status && rich->self_cld_misjudge > 0.1) return false; 
    
    if (rich->self_expnpe < 0.01) return false;
    if (rich->self_dist < 0.25) return false;
    
    if (!rich->self_cld_status) return false;
    if (!rich->status) return false;

    // Chi-square Cut
    const double vel_nchi_cut = 2.00;
    const double trk_nchi_cut = 1.75;
    
    // Velocity
    if (!hyc->vel_status[2]) return false;
    if (std::log(hyc->vel_nchi[2]) > vel_nchi_cut) return false;
    
    // Geometry
    if (!hyc->geom_status[0]) return false;
    if (std::log(hyc->geom_nchi_x[0]) > trk_nchi_cut) return false;
    if (std::log(hyc->geom_nchi_y[0]) > trk_nchi_cut) return false;
    
    // Mass
    if (!hyc->mutr_status[2]) return false;
    if (std::log(hyc->mutr_nchi_x[2]) > trk_nchi_cut) return false;
    if (std::log(hyc->mutr_nchi_y[2]) > trk_nchi_cut) return false;
    if (std::log(hyc->mutr_nchi_b[2]) > vel_nchi_cut) return false;

    // Physics
    if (!hyc->phys_status[0][2]) return false;
    if (std::log(hyc->phys_nchi_x[0][2]) > trk_nchi_cut) return false;
    if (std::log(hyc->phys_nchi_y[0][2]) > trk_nchi_cut) return false;
    if (std::log(hyc->phys_nchi_b[0][2]) > vel_nchi_cut) return false;
    
    // Choutko track
    if (!trk->ck_status[0]) return false;
    if (std::log(trk->ck_nchi[0][0]) > 3.0) return false;
    if (std::log(trk->ck_nchi[0][1]) > 3.0) return false;

    // Variables
    double trdllr  = trd->tdLLR_ep;
    double sqrm    = hyc->mutr_sqrm[2];
    
    double bta     = hyc->phys_top_bta[0][2];
    double rig     = hyc->phys_top_rig[0][2];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;

    // Cut on in chi-square
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);
   
    varsLRH.init();
    varsLRH.run = list->run;
    varsLRH.evt = list->event;
    varsLRH.ut  = list->utime;
    varsLRH.wgt = list->weight;

    varsLRH.mc     = (g4mc != nullptr);
    varsLRH.mc_rig = mc_rig;
    varsLRH.mc_w10 = mc_w10;
    varsLRH.mc_w27 = mc_w27;

    varsLRH.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsLRH.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsLRH.lv    = CheckType(Type::ISS) ? rti->livetime : 1;

    varsLRH.sign   = signr;
    varsLRH.rig    = rig;
    varsLRH.bta    = bta;
    varsLRH.sqrm   = sqrm;
    varsLRH.llr    = trdllr;
    
    varsLRH.lxin   = lchix_in;
    varsLRH.lyin   = lchiy_in;
    
    return true;
}

bool Selector::process_data_iin() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // TRD
    if (!trd->tkLLR_status || trd->tkLLR_num_hit < 8 || trd->num_tkHit < 6) return false;

    // Geometry
    if (!hyc->geom_status[0]) return false;
    
    // TOF
    if (tof->extcls_noise != 0) return false;
    if (tof->num_in_time_cls > 4) return false;
    
    // Choutko track
    if (!trk->ck_status[0]) return false;
    if (std::log(trk->ck_nchi[0][0]) > 3.0) return false;
    if (std::log(trk->ck_nchi[0][1]) > 3.0) return false;

    // Variables
    double trdllr  = trd->tkLLR_ep;
    
    double rig     = hyc->geom_top_rig[0];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    // Cut on in chi-square
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);
   
    varsIIN.init();
    varsIIN.run = list->run;
    varsIIN.evt = list->event;
    varsIIN.ut  = list->utime;
    varsIIN.wgt = list->weight;

    varsIIN.mc     = (g4mc != nullptr);
    varsIIN.mc_rig = mc_rig;
    varsIIN.mc_w10 = mc_w10;
    varsIIN.mc_w27 = mc_w27;

    varsIIN.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsIIN.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsIIN.lv    = CheckType(Type::ISS) ? rti->livetime : 1;

    varsIIN.sign   = signr;
    varsIIN.rig    = rig;
    varsIIN.llr    = trdllr;
    varsIIN.ecal   = ecal->status;
    varsIIN.engE   = ecal->engE;
    varsIIN.engD   = ecal->engD;
    varsIIN.mvaBDT = ecal->mvaBDT;
    
    varsIIN.lxin = lchix_in;
    varsIIN.lyin = lchiy_in;
    
    return true;
}

bool Selector::process_data_iex() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // TRD
    if (!trd->tkLLR_status || trd->tkLLR_num_hit < 8 || trd->num_tkHit < 6) return false;
    
    int pt = -1;
    if ((pt < 0) && (hyc->geom_status[3] && trk->lay[0] == 3 && trk->lay[8] == 3                    )) pt = 3;
    if ((pt < 0) && (hyc->geom_status[2] &&                     trk->lay[8] == 3 && trk->lay[1] == 3)) pt = 2;
    if ((pt < 0) && (hyc->geom_status[1] && trk->lay[0] == 3                                        )) pt = 1;
    if ((pt < 0)) return false;
    
    if (!hyc->geom_status[0] ) return false;
    if (!hyc->geom_status[pt]) return false;
    
    // Choutko track
    if (!trk->ck_status[0] ) return false;
    if (std::log(trk->ck_nchi[0][0]) > 3.0) return false;
    if (std::log(trk->ck_nchi[0][1]) > 3.0) return false;
    
    if (!trk->ck_status[pt]) return false;
    if (std::log(trk->ck_nchi[pt][0]) > 3.0) return false;
    if (std::log(trk->ck_nchi[pt][1]) > 3.0) return false;

    double trdllr  = trd->tkLLR_ep;
    double rig     = hyc->geom_top_rig[pt];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);
    double lchix_ex = std::log(hyc->geom_nchi_x[pt]);
    double lchiy_ex = std::log(hyc->geom_nchi_y[pt]);
   
    varsIEX.init();
    varsIEX.run = list->run;
    varsIEX.evt = list->event;
    varsIEX.ut  = list->utime;
    varsIEX.wgt = list->weight;

    varsIEX.mc     = (g4mc != nullptr);
    varsIEX.mc_rig = mc_rig;
    varsIEX.mc_w10 = mc_w10;
    varsIEX.mc_w27 = mc_w27;

    varsIEX.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsIEX.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsIEX.lv    = CheckType(Type::ISS) ? rti->livetime : 1;

    varsIEX.sign   = signr;
    varsIEX.rig    = rig;
    varsIEX.llr    = trdllr;
    varsIEX.ecal   = ecal->status;
    varsIEX.engE   = ecal->engE;
    varsIEX.engD   = ecal->engD;
    varsIEX.mvaBDT = ecal->mvaBDT;

    varsIEX.patt = pt;
    varsIEX.lxin = lchix_in;
    varsIEX.lyin = lchiy_in;
    varsIEX.lxex = lchix_ex;
    varsIEX.lyex = lchiy_ex;
    
    return true;
}

bool Selector::process_data_hex() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // TRD
    if (!trd->tkLLR_status || trd->tkLLR_num_hit < 8 || trd->num_tkHit < 6) return false;
    
    // ECAL
    //if (ecal->status && ecal->mvaBDT > -0.6) return false;
   
    int pt = -1;
    if ((pt < 0) && (hyc->geom_status[3] && trk->lay[0] == 3 && trk->lay[8] == 3                    )) pt = 3;
    if ((pt < 0) && (hyc->geom_status[2] &&                     trk->lay[8] == 3 && trk->lay[1] == 3)) pt = 2;
    if ((pt < 0) && (hyc->geom_status[1] && trk->lay[0] == 3                                        )) pt = 1;
    if (pt < 0) return false;
    
    if (!hyc->geom_status[0]) return false;
    if (!hyc->geom_status[pt]) return false;

    // Choutko track
    if (!trk->ck_status[0] ) return false;
    if (!trk->ck_status[pt]) return false;

    double trdllr  = trd->tkLLR_ep;
    double rig     = hyc->geom_top_rig[pt];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);
    double lchix_ex = std::log(hyc->geom_nchi_x[pt]);
    double lchiy_ex = std::log(hyc->geom_nchi_y[pt]);
    
    double lchixC_in = std::log(trk->ck_nchi[0][0]);
    double lchiyC_in = std::log(trk->ck_nchi[0][1]);
    double lchixC_ex = std::log(trk->ck_nchi[pt][0]);
    double lchiyC_ex = std::log(trk->ck_nchi[pt][1]);
    
    double nlxin = hyc->geom_max_norm_lx[0];
    double nlyin = hyc->geom_max_norm_ly[0];
    double nlxex = hyc->geom_max_norm_lx[pt];
    double nlyex = hyc->geom_max_norm_ly[pt];
    
    varsHEX.init();
    varsHEX.run = list->run;
    varsHEX.evt = list->event;
    varsHEX.ut  = list->utime;
    varsHEX.wgt = list->weight;

    varsHEX.mc     = (g4mc != nullptr);
    varsHEX.mc_rig = mc_rig;
    varsHEX.mc_w10 = mc_w10;
    varsHEX.mc_w27 = mc_w27;

    varsHEX.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsHEX.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsHEX.lv    = CheckType(Type::ISS) ? rti->livetime : 1;

    varsHEX.sign   = signr;
    varsHEX.rig    = rig;
    varsHEX.llr    = trdllr;
    varsHEX.ecal   = ecal->status;
    varsHEX.engE   = ecal->engE;
    varsHEX.engD   = ecal->engD;
    varsHEX.mvaBDT = ecal->mvaBDT;
    
    varsHEX.patt = pt;
    varsHEX.lxin = lchix_in;
    varsHEX.lyin = lchiy_in;
    varsHEX.lxex = lchix_ex;
    varsHEX.lyex = lchiy_ex;
    
    varsHEX.nlxin = nlxin;
    varsHEX.nlyin = nlyin;
    varsHEX.nlxex = nlxex;
    varsHEX.nlyex = nlyex;

    varsHEX.lxinC = lchixC_in;
    varsHEX.lyinC = lchiyC_in;
    varsHEX.lxexC = lchixC_ex;
    varsHEX.lyexC = lchiyC_ex;
    
    return true;
}


bool Selector::process_data_hfs() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // TRD
    if (!trd->tkLLR_status || trd->tkLLR_num_hit < 8 || trd->num_tkHit < 6) return false;
    
    // ECAL
    //if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    // Tracker-FS
    if (trk->lay[0] != 3) return false;
    if (trk->lay[8] != 3) return false;
    
    if (!hyc->geom_status[0]) return false;
    if (!hyc->geom_status[1]) return false;
    if (!hyc->geom_status[2]) return false;
    if (!hyc->geom_status[3]) return false;

    // Choutko track
    if (!trk->ck_status[0]) return false;
    if (!trk->ck_status[1]) return false;
    if (!trk->ck_status[2]) return false;
    if (!trk->ck_status[3]) return false;
    
    double trdllr  = trd->tkLLR_ep;
    double rig     = hyc->geom_top_rig[3];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);
    double lchix_l1 = std::log(hyc->geom_nchi_x[1]);
    double lchiy_l1 = std::log(hyc->geom_nchi_y[1]);
    double lchix_l9 = std::log(hyc->geom_nchi_x[2]);
    double lchiy_l9 = std::log(hyc->geom_nchi_y[2]);
    double lchix_fs = std::log(hyc->geom_nchi_x[3]);
    double lchiy_fs = std::log(hyc->geom_nchi_y[3]);
    
    double lchixC_in = std::log(trk->ck_nchi[0][0]);
    double lchiyC_in = std::log(trk->ck_nchi[0][1]);
    double lchixC_l1 = std::log(trk->ck_nchi[1][0]);
    double lchiyC_l1 = std::log(trk->ck_nchi[1][1]);
    double lchixC_l9 = std::log(trk->ck_nchi[2][0]);
    double lchiyC_l9 = std::log(trk->ck_nchi[2][1]);
    double lchixC_fs = std::log(trk->ck_nchi[3][0]);
    double lchiyC_fs = std::log(trk->ck_nchi[3][1]);
    
    double lrin = std::log(1.0 + std::abs(hyc->geom_top_rig[0] / hyc->geom_top_rig[2] - hyc->geom_top_rig[0] / hyc->geom_top_rig[1]));
    double lrfs = std::log(1.0 + std::abs(hyc->geom_top_rig[3] / hyc->geom_top_rig[2] - hyc->geom_top_rig[3] / hyc->geom_top_rig[1]));
    
    double lrinC = std::log(1.0 + std::abs(trk->ck_rig[0] / trk->ck_rig[2] - trk->ck_rig[0] / trk->ck_rig[1]));
    double lrfsC = std::log(1.0 + std::abs(trk->ck_rig[3] / trk->ck_rig[2] - trk->ck_rig[3] / trk->ck_rig[1]));
    
    double nlxin = hyc->geom_max_norm_lx[0];
    double nlyin = hyc->geom_max_norm_ly[0];
    double nlxl1 = hyc->geom_max_norm_lx[1];
    double nlyl1 = hyc->geom_max_norm_ly[1];
    double nlxl9 = hyc->geom_max_norm_lx[2];
    double nlyl9 = hyc->geom_max_norm_ly[2];
    double nlxfs = hyc->geom_max_norm_lx[3];
    double nlyfs = hyc->geom_max_norm_ly[3];
    
    varsHFS.init();
    varsHFS.run = list->run;
    varsHFS.evt = list->event;
    varsHFS.ut  = list->utime;
    varsHFS.wgt = list->weight;

    varsHFS.mc     = (g4mc != nullptr);
    varsHFS.mc_rig = mc_rig;
    varsHFS.mc_w10 = mc_w10;
    varsHFS.mc_w27 = mc_w27;

    varsHFS.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsHFS.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsHFS.lv    = CheckType(Type::ISS) ? rti->livetime : 1;

    varsHFS.sign   = signr;
    varsHFS.rig    = rig;
    varsHFS.llr    = trdllr;
    varsHFS.ecal   = ecal->status;
    varsHFS.engE   = ecal->engE;
    varsHFS.engD   = ecal->engD;
    varsHFS.mvaBDT = ecal->mvaBDT;

    varsHFS.lxin = lchix_in;
    varsHFS.lyin = lchiy_in;
    varsHFS.lxl1 = lchix_l1;
    varsHFS.lyl1 = lchiy_l1;
    varsHFS.lxl9 = lchix_l9;
    varsHFS.lyl9 = lchiy_l9;
    varsHFS.lxfs = lchix_fs;
    varsHFS.lyfs = lchiy_fs;

    varsHFS.nlxin = nlxin;
    varsHFS.nlyin = nlyin;
    varsHFS.nlxl1 = nlxl1;
    varsHFS.nlyl1 = nlyl1;
    varsHFS.nlxl9 = nlxl9;
    varsHFS.nlyl9 = nlyl9;
    varsHFS.nlxfs = nlxfs;
    varsHFS.nlyfs = nlyfs;
    
    varsHFS.lxinC = lchixC_in;
    varsHFS.lyinC = lchiyC_in;
    varsHFS.lxl1C = lchixC_l1;
    varsHFS.lyl1C = lchiyC_l1;
    varsHFS.lxl9C = lchixC_l9;
    varsHFS.lyl9C = lchiyC_l9;
    varsHFS.lxfsC = lchixC_fs;
    varsHFS.lyfsC = lchiyC_fs;
    
    varsHFS.lrin  = lrin;
    varsHFS.lrfs  = lrfs;
    varsHFS.lrinC = lrinC;
    varsHFS.lrfsC = lrfsC;
    
    return true;
}


#endif // __Selector_C__
