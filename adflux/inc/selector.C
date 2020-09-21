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
    varsTF.init();
    varsRH.init();

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
    //if ((trg->bit&2) != 2 && (trg->bit&8) != 8) return false;
    
    // Charge
    short Zall = tof->Zall;
    //if (tof->Qall < 0.8 || tof->Qall > 1.4) return false;
    //if (trk->QIn < 0.8 || trk->QIn > 1.4) return false;
   
    //if (hyc->tfQup > 2.0) return false;
    //if (hyc->tfQlw > 2.0) return false;

    // TOF
    //if (tof->num_extcls[0] > 1 || tof->num_extcls[1] > 1) return false;
    //if (tof->num_extcls[2] > 1 || tof->num_extcls[3] > 1) return false;

    // Tracker
    //if (trk->lay[1] == 0) return false;

    // Tracker SN10
    std::vector<double> sn10;
    for (int il = 2; il <= 7; ++il) {
        if (trk->lay[il] == 0) continue;
        sn10.push_back(trk->sn10[il]);
    }
    double min_sn10 = *std::min_element(sn10.begin(), sn10.end());
    //if (min_sn10 < (1.0/3.0)) return false;

    //if (trk->lay[1] != 0 && trk->sn10[1] < (2.0/3.0)) return false;
    //if (trk->lay[0] != 0 && trk->sn10[0] < (2.0/3.0)) return false;
    //if (trk->lay[8] != 0 && trk->sn10[8] < (2.0/3.0)) return false;
    //
    //if (trk->lay[1] != 0 && trk->ext_num_hit[1] > 3) return false;
    //if (trk->lay[0] != 0 && trk->ext_num_hit[0] > 3) return false;
    //if (trk->lay[8] != 0 && trk->ext_num_hit[8] > 3) return false;
   
    return true;
}

bool Selector::process_data() {
    bool pass_tf = process_data_tf();
    bool pass_rh = process_data_rh();

    if (pass_tf) varsTF.fill();
    if (pass_rh) varsRH.fill();

    return true;
}

bool Selector::process_data_tf() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
 /*   
    // Trigger
    //if ((trg->bit&8) != 8) return false;
    
    // TRD
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;

    // TOF
    if (tof->extcls_noise != 0) return false;
    if (tof->num_in_time_cls > 4) return false;

    //==== RICH Veto ====//
    bool RICHVeto = false;
    while (true) {
        // Self RICH
        if (rich->self_status) {
            if (rich->self_kind == 1 && rich->self_num_cloud != 0) break;
            if (rich->self_kind == 2 && rich->self_num_cloud >  1) break;
            if (rich->self_num_tumor != 0) break;
            if (rich->self_num_ghost != 0) break;
            
            if (rich->self_num_stone > 1) break;
            if (rich->self_stn_status && rich->self_stn_dist > 3.4) break;
        }
        
        // Others RICH
        if (rich->new_status) {
            if (rich->self_kind == 1 && rich->new_num_cloud != 0) break;
            if (rich->self_kind == 2 && rich->new_num_cloud >  1) break;
            if (rich->new_num_ghost != 0) break;
            
            if (rich->new_num_stone > 1) break;
            if (rich->new_stn_status && rich->new_stn_dist > 3.4) break;
        }
        
        // Official RICH
        if (rich->status && rich->kind == 1) break;
        //if (rich->num_clsZ1 > 3) break;

        RICHVeto = true;
        break;
    }
    //if (!RICHVeto) return false;
*/
    // TRD extra hit
    int nvtx_xl = (trd->num_vtx[0][0] + trd->num_vtx[1][0]);
    int nvtx_xu = (trd->num_vtx[2][0] + trd->num_vtx[3][0]);
    int nvtx_yl = (trd->num_vtx[0][1] + trd->num_vtx[1][1]);
    int nvtx_yu = (trd->num_vtx[2][1] + trd->num_vtx[3][1]);

    // Velocity
    if (!hyc->vel_status[0]) return false;
    if (!hyc->vel_status[1]) return false;
   
    // Geometry
    if (!hyc->geom_status[0]) return false;

    // Mass
    if (!hyc->mutr_status[0]) return false;
    if (!hyc->mutr_status[1]) return false;

    // Physics
    if (!hyc->phys_status[0][0]) return false;
    if (!hyc->phys_status[1][0]) return false;
    if (!hyc->phys_status[0][1]) return false;
    if (!hyc->phys_status[1][1]) return false;

    // Choutko track
    if (!trk->ck_status[0]) return false;
   
    // Variables
    short  chrg    = tof->Zall;
    double trdllr  = trd->tdLLR_ep;
    double sqrm    = hyc->mutr_sqrm[1];

    double bta     = hyc->vel_top_bta[1];
    double rig     = hyc->geom_top_rig[0];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;

    //if (hyc->geom_status[1] && (signr * hyc->geom_top_rig[1]) < 0.0) return false;
    //if (hyc->geom_status[2] && (signr * hyc->geom_top_rig[2]) < 0.0) return false;
    //if (hyc->geom_status[3] && (signr * hyc->geom_top_rig[3]) < 0.0) return false;
    //
    //if (trk->ck_status[0] && (signr * trk->ck_rig[0]) < 0.0) return false;
    //if (trk->ck_status[1] && (signr * trk->ck_rig[1]) < 0.0) return false;
    //if (trk->ck_status[2] && (signr * trk->ck_rig[2]) < 0.0) return false;
    //if (trk->ck_status[3] && (signr * trk->ck_rig[3]) < 0.0) return false;
   
    bool   ext   = false;
    double extlx = 0.0;
    double extly = 0.0;
    for (int ip = 1; ip <= 3; ++ip) {
        if (!hyc->geom_status[ip]) continue;
        double iextlx = std::log(hyc->geom_nchi_x[ip]);
        double iextly = std::log(hyc->geom_nchi_y[ip]); 
        extlx = ext ? std::max(extlx, iextlx) : iextlx;
        extly = ext ? std::max(extly, iextly) : iextly;
        ext = true;
    }
    
    // Feet
    double feetraw[7] = {
        (trk->feet[1] > 0) ? trk->feet[1] : 10.0,
        (trk->feet[2] > 0) ? trk->feet[2] : 10.0, (trk->feet[3] > 0) ? trk->feet[3] : 10.0,
        (trk->feet[4] > 0) ? trk->feet[4] : 10.0, (trk->feet[5] > 0) ? trk->feet[5] : 10.0,
        (trk->feet[6] > 0) ? trk->feet[6] : 10.0, (trk->feet[7] > 0) ? trk->feet[7] : 10.0
    };
    double feetL2  = feetraw[0];
    double feetSL1 = std::min(feetraw[1], feetraw[2]);
    double feetSL2 = std::min(feetraw[3], feetraw[4]);
    double feetSL3 = std::min(feetraw[5], feetraw[6]);
    
    varsTF.init();
    varsTF.run = list->run;
    varsTF.evt = list->event;
    varsTF.ut  = list->utime;
    varsTF.wgt = list->weight;

    varsTF.mc     = (g4mc != nullptr);
    varsTF.mc_rig = mc_rig;
    varsTF.mc_w10 = mc_w10;
    varsTF.mc_w27 = mc_w27;

    varsTF.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsTF.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsTF.lv    = CheckType(Type::ISS) ? rti->livetime : 1;
    
    varsTF.trg   = ((trg->bit&8) == 8);

    varsTF.chrg    = chrg;
    varsTF.sign    = signr;
    varsTF.rig     = rig;
    varsTF.bta     = bta;
    varsTF.sqrm    = sqrm;
    varsTF.llr     = trdllr;
    
    varsTF.L2  = trk->lay[1];
    varsTF.nhx = trk->num_inn_x;
    varsTF.nhy = trk->num_inn_y;
   
    varsTF.ftL2  = feetL2;
    varsTF.ftSL1 = feetSL1;
    varsTF.ftSL2 = feetSL2;
    varsTF.ftSL3 = feetSL3;

    varsTF.rich  = rich->self_status;
    //varsTF.veto  = RICHVeto;
    varsTF.ncls  = rich->num_clsZ1;
    varsTF.nhinn = (rich->self_status) ? rich->self_nhit_other_inn : 0;
    varsTF.nhout = (rich->self_status) ? rich->self_nhit_other_out : 0;
    varsTF.nhinn2 = (rich->new_status) ? rich->new_nhit_other_inn : 0;
    varsTF.nhout2 = (rich->new_status) ? rich->new_nhit_other_out : 0;
    
    varsTF.nvtxx[0] = nvtx_xl;
    varsTF.nvtxx[1] = nvtx_xu;
    varsTF.nvtxy[0] = nvtx_yl;
    varsTF.nvtxy[1] = nvtx_yu;
    
    varsTF.ext   = ext;
    varsTF.extlx = extlx;
    varsTF.extly = extly;
    
    varsTF.geom_lx = std::log(hyc->geom_nchi_x[0]);
    varsTF.geom_ly = std::log(hyc->geom_nchi_y[0]);
        
    varsTF.tf_bta[0] = hyc->vel_cen_bta[0];
    varsTF.tf_bta[1] = hyc->vel_cen_bta[1];

    for (int ip = 0; ip <= 1; ++ip) {
        varsTF.vel_lb[ip]  = std::log(hyc->vel_nchi[ip]);
        
        varsTF.mutr_lx[ip] = std::log(hyc->mutr_nchi_x[ip]);
        varsTF.mutr_ly[ip] = std::log(hyc->mutr_nchi_y[ip]);
        varsTF.mutr_lb[ip] = std::log(hyc->mutr_nchi_b[ip]);
        
        varsTF.a_phys_lx[ip] = std::log(hyc->phys_nchi_x[0][ip]);
        varsTF.a_phys_ly[ip] = std::log(hyc->phys_nchi_y[0][ip]);
        varsTF.a_phys_lb[ip] = std::log(hyc->phys_nchi_b[0][ip]);
                      
        varsTF.b_phys_lx[ip] = std::log(hyc->phys_nchi_x[1][ip]);
        varsTF.b_phys_ly[ip] = std::log(hyc->phys_nchi_y[1][ip]);
        varsTF.b_phys_lb[ip] = std::log(hyc->phys_nchi_b[1][ip]);
        
        varsTF.a_top_rig[ip] = hyc->phys_top_rig[0][ip];
        varsTF.a_top_bta[ip] = hyc->phys_top_bta[0][ip];
        varsTF.b_top_rig[ip] = hyc->phys_top_rig[1][ip];
        varsTF.b_top_bta[ip] = hyc->phys_top_bta[1][ip];
        
        varsTF.a_cen_rig[ip] = hyc->phys_cen_rig[0][ip];
        varsTF.a_cen_bta[ip] = hyc->phys_cen_bta[0][ip];
        varsTF.b_cen_rig[ip] = hyc->phys_cen_rig[1][ip];
        varsTF.b_cen_bta[ip] = hyc->phys_cen_bta[1][ip];

        //if ((signr * varsTF.a_cen_rig[ip]) < 0.0) return false;
        //if ((signr * varsTF.b_cen_rig[ip]) < 0.0) return false;
    }

    return true;
}

bool Selector::process_data_rh() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    //if ((trg->bit&8) != 8) return false;
    
    // TRD
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    // TRD extra hit
    int nvtx_xl = (trd->num_vtx[0][0] + trd->num_vtx[1][0]);
    int nvtx_xu = (trd->num_vtx[2][0] + trd->num_vtx[3][0]);
    int nvtx_yl = (trd->num_vtx[0][1] + trd->num_vtx[1][1]);
    int nvtx_yu = (trd->num_vtx[2][1] + trd->num_vtx[3][1]);
    
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
    if (rich->self_trace < 0.30) return false;
    
    if (rich->self_stn_status && rich->self_stn_dist > 3.4) return false;
    if (rich->self_cld_status && rich->self_cld_misjudge > 0.1) return false; 
    
    if (rich->self_expnpe < 2.0) return false;
    if (rich->self_dist < 0.25) return false;
    
    if (!rich->self_cld_status) return false;

    // Velocity
    if (!hyc->vel_status[0]) return false;
    if (!hyc->vel_status[1]) return false;
    if (!hyc->vel_status[2]) return false;
    if (!hyc->vel_status[3]) return false;
    
    // Geometry
    if (!hyc->geom_status[0]) return false;
    
    // Mass
    if (!hyc->mutr_status[0]) return false;
    if (!hyc->mutr_status[1]) return false;
    if (!hyc->mutr_status[2]) return false;
    if (!hyc->mutr_status[3]) return false;

    // Physics
    if (!hyc->phys_status[0][0]) return false;
    if (!hyc->phys_status[1][0]) return false;
    if (!hyc->phys_status[0][1]) return false;
    if (!hyc->phys_status[1][1]) return false;
    if (!hyc->phys_status[0][2]) return false;
    if (!hyc->phys_status[1][2]) return false;
    if (!hyc->phys_status[0][3]) return false;
    if (!hyc->phys_status[1][3]) return false;
    
    // Choutko track
    if (!trk->ck_status[0]) return false;
    
    // Variables
    short  chrg    = tof->Zall;
    double trdllr  = trd->tdLLR_ep;
    double sqrm    = hyc->mutr_sqrm[3];
    
    double bta     = hyc->vel_top_bta[3];
    double rig     = hyc->geom_top_rig[0];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    if (hyc->geom_status[1] && (signr * hyc->geom_top_rig[1]) < 0.0) return false;
    if (hyc->geom_status[2] && (signr * hyc->geom_top_rig[2]) < 0.0) return false;
    if (hyc->geom_status[3] && (signr * hyc->geom_top_rig[3]) < 0.0) return false;
    
    if (trk->ck_status[0] && (signr * trk->ck_rig[0]) < 0.0) return false;
    if (trk->ck_status[1] && (signr * trk->ck_rig[1]) < 0.0) return false;
    if (trk->ck_status[2] && (signr * trk->ck_rig[2]) < 0.0) return false;
    if (trk->ck_status[3] && (signr * trk->ck_rig[3]) < 0.0) return false;
    
    bool   ext   = false;
    double extlx = 0.0;
    double extly = 0.0;
    for (int ip = 1; ip <= 3; ++ip) {
        if (!hyc->geom_status[ip]) continue;
        double iextlx = std::log(hyc->geom_nchi_x[ip]);
        double iextly = std::log(hyc->geom_nchi_y[ip]); 
        extlx = ext ? std::max(extlx, iextlx) : iextlx;
        extly = ext ? std::max(extly, iextly) : iextly;
        ext = true;
    }
    
    // Others RICH
    if (rich->new_status) {
        if (rich->new_num_stone > 1) return false;
        if (rich->new_num_cloud > 1) return false;
        if (rich->new_num_ghost != 0) return false;
        
        if (rich->new_stn_status && rich->new_stn_dist > 3.4) return false;
    }
    else return false;

    // Official RICH
    const short cut_pmt = 3;
    const float cut_prob = 0.01;
    const float cut_betaCstc = 0.005;
    const float cut_chrgCstc = 5.0;
    const float cut_expPhe = 2.0;
    const float cut_collPhe = 0.4;
    if (rich->status && rich->kind == 1) { 
        if (rich->npmt < cut_pmt) return false;
        if (rich->prob < cut_prob) return false;
   
        if (rich->cstcb > cut_betaCstc) return false;
        if (rich->cstcq > cut_chrgCstc) return false;
        if (rich->num_expPE < cut_expPhe) return false;
        if (rich->eft_colPE < cut_collPhe) return false;
    }
    else return false;

    // Compare Official and New
    if (rich->self_cld_status && (rich->status && rich->kind == 1) &&
        std::abs(rich->self_cld_cbta - rich->beta) > cut_betaCstc) return false;

    if (rich->new_cld_status && (rich->status && rich->kind == 1) &&
        std::abs(rich->new_cld_cbta - rich->beta) > cut_betaCstc) return false;
    
    if (rich->self_cld_status && rich->new_cld_status && 
        std::abs(rich->self_cld_cbta - rich->new_cld_cbta) > cut_betaCstc) return false;

    if (std::abs(hyc->vel_cen_bta[2] - hyc->vel_cen_bta[0]) > 0.1) return false;
    if (std::abs(hyc->vel_cen_bta[2] - hyc->vel_cen_bta[1]) > 0.1) return false;
    if (std::abs(hyc->vel_cen_bta[3] - hyc->vel_cen_bta[0]) > 0.1) return false;
    if (std::abs(hyc->vel_cen_bta[3] - hyc->vel_cen_bta[1]) > 0.1) return false;
    
    // Feet
    double feetraw[7] = {
        (trk->feet[1] > 0) ? trk->feet[1] : 10.0,
        (trk->feet[2] > 0) ? trk->feet[2] : 10.0, (trk->feet[3] > 0) ? trk->feet[3] : 10.0,
        (trk->feet[4] > 0) ? trk->feet[4] : 10.0, (trk->feet[5] > 0) ? trk->feet[5] : 10.0,
        (trk->feet[6] > 0) ? trk->feet[6] : 10.0, (trk->feet[7] > 0) ? trk->feet[7] : 10.0
    };
    double feetL2  = feetraw[0];
    double feetSL1 = std::min(feetraw[1], feetraw[2]);
    double feetSL2 = std::min(feetraw[3], feetraw[4]);
    double feetSL3 = std::min(feetraw[5], feetraw[6]);

    varsRH.init();
    varsRH.run = list->run;
    varsRH.evt = list->event;
    varsRH.ut  = list->utime;
    varsRH.wgt = list->weight;

    varsRH.mc     = (g4mc != nullptr);
    varsRH.mc_rig = mc_rig;
    varsRH.mc_w10 = mc_w10;
    varsRH.mc_w27 = mc_w27;

    varsRH.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsRH.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsRH.lv    = CheckType(Type::ISS) ? rti->livetime : 1;
    
    varsRH.trg   = ((trg->bit&8) == 8);

    varsRH.chrg   = chrg;
    varsRH.sign   = signr;
    varsRH.rig    = rig;
    varsRH.bta    = bta;
    varsRH.sqrm   = sqrm;
    varsRH.llr    = trdllr;
    
    varsRH.L2  = trk->lay[1];
    varsRH.nhx = trk->num_inn_x;
    varsRH.nhy = trk->num_inn_y;
    
    varsRH.ftL2  = feetL2;
    varsRH.ftSL1 = feetSL1;
    varsRH.ftSL2 = feetSL2;
    varsRH.ftSL3 = feetSL3;
    
    varsRH.nhinn = rich->self_nhit_other_inn;
    varsRH.nhout = rich->self_nhit_other_out;
    
    varsRH.nhinn2 = rich->new_nhit_other_inn;
    varsRH.nhout2 = rich->new_nhit_other_out;
    
    varsRH.nclsZ1 = rich->num_clsZ1;
    varsRH.selfqq = (rich->self_cld_npe / rich->self_expnpe);
  
    varsRH.rhbta[0] = rich->self_cld_cbta;
    varsRH.rhbta[1] = rich->beta;

    varsRH.nvtxx[0] = nvtx_xl;
    varsRH.nvtxx[1] = nvtx_xu;
    varsRH.nvtxy[0] = nvtx_yl;
    varsRH.nvtxy[1] = nvtx_yu;

    varsRH.ext   = ext;
    varsRH.extlx = extlx;
    varsRH.extly = extly;

    varsRH.geom_lx = std::log(hyc->geom_nchi_x[0]);
    varsRH.geom_ly = std::log(hyc->geom_nchi_y[0]);
    
    varsRH.tf_bta[0] = hyc->vel_cen_bta[0];
    varsRH.tf_bta[1] = hyc->vel_cen_bta[1];
    varsRH.rh_bta[0] = hyc->vel_cen_bta[2];
    varsRH.rh_bta[1] = hyc->vel_cen_bta[3];
    
    for (int ip = 0; ip <= 3; ++ip) {
        varsRH.vel_lb[ip]  = std::log(hyc->vel_nchi[ip]);
        
        varsRH.mutr_lx[ip] = std::log(hyc->mutr_nchi_x[ip]);
        varsRH.mutr_ly[ip] = std::log(hyc->mutr_nchi_y[ip]);
        varsRH.mutr_lb[ip] = std::log(hyc->mutr_nchi_b[ip]);
        
        varsRH.a_phys_lx[ip] = std::log(hyc->phys_nchi_x[0][ip]);
        varsRH.a_phys_ly[ip] = std::log(hyc->phys_nchi_y[0][ip]);
        varsRH.a_phys_lb[ip] = std::log(hyc->phys_nchi_b[0][ip]);
                      
        varsRH.b_phys_lx[ip] = std::log(hyc->phys_nchi_x[1][ip]);
        varsRH.b_phys_ly[ip] = std::log(hyc->phys_nchi_y[1][ip]);
        varsRH.b_phys_lb[ip] = std::log(hyc->phys_nchi_b[1][ip]);
        
        varsRH.a_top_rig[ip] = hyc->phys_top_rig[0][ip];
        varsRH.a_top_bta[ip] = hyc->phys_top_bta[0][ip];
        varsRH.b_top_rig[ip] = hyc->phys_top_rig[1][ip];
        varsRH.b_top_bta[ip] = hyc->phys_top_bta[1][ip];
        
        varsRH.a_cen_rig[ip] = hyc->phys_cen_rig[0][ip];
        varsRH.a_cen_bta[ip] = hyc->phys_cen_bta[0][ip];
        varsRH.b_cen_rig[ip] = hyc->phys_cen_rig[1][ip];
        varsRH.b_cen_bta[ip] = hyc->phys_cen_bta[1][ip];
        
        if ((signr * varsRH.a_cen_rig[ip]) < 0.0) return false;
        if ((signr * varsRH.b_cen_rig[ip]) < 0.0) return false;
    }

    return true;
}

#endif // __Selector_C__
