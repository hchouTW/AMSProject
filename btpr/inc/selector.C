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
    
    TrSys::PhysEnv::ReadMagAMS("/eos/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    TrSys::PhysEnv::ReadMatAMS("/eos/user/h/hchou/ExternalLibs/DB/material");
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
    varsBT.init();

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
    //short Zall = tof->Zall;
    //if (tof->Qall < 0.8 || tof->Qall > 1.4) return false;
    //if (trk->QIn < 0.8 || trk->QIn > 1.4) return false;
   
    //if (hyc->tfQup > 2.0) return false;
    //if (hyc->tfQlw > 2.0) return false;

    //// TOF
    //if (tof->num_extcls[0] > 1 || tof->num_extcls[1] > 1) return false;
    //if (tof->num_extcls[2] > 1 || tof->num_extcls[3] > 1) return false;

    //// Tracker
    //if (trk->lay[1] == 0) return false;

    // Tracker SN10
    std::vector<double> sn10;
    for (int il = 2; il <= 7; ++il) {
        if (trk->lay[il] == 0) continue;
        sn10.push_back(trk->sn10[il]);
    }
    double min_sn10 = *std::min_element(sn10.begin(), sn10.end());
    if (min_sn10 < (1.0/3.0)) return false;

    if (trk->lay[1] != 0 && trk->sn10[1] < (2.0/3.0)) return false;
    if (trk->lay[0] != 0 && trk->sn10[0] < (2.0/3.0)) return false;
    if (trk->lay[8] != 0 && trk->sn10[8] < (2.0/3.0)) return false;
    
    if (trk->lay[1] != 0 && trk->ext_num_hit[1] > 3) return false;
    if (trk->lay[0] != 0 && trk->ext_num_hit[0] > 3) return false;
    if (trk->lay[8] != 0 && trk->ext_num_hit[8] > 3) return false;
   
    return true;
}

bool Selector::process_data() {
    bool pass_bt = process_data_bt();

    if (pass_bt) varsBT.fill();

    return true;
}

bool Selector::process_data_bt() {
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
    
    // TOF
    if (tof->extcls_noise != 0) return false;
    if (tof->num_in_time_cls > 4) return false;

    // Variables
    varsBT.init();
    varsBT.run = list->run;
    varsBT.evt = list->event;
    varsBT.ut  = list->utime;
    varsBT.wgt = list->weight;

    varsBT.mc     = (g4mc != nullptr);
    varsBT.mc_rig = mc_rig;
    varsBT.mc_w10 = mc_w10;
    varsBT.mc_w27 = mc_w27;

    varsBT.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    varsBT.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    varsBT.lv    = CheckType(Type::ISS) ? rti->livetime : 1;
    
    varsBT.trg   = ((trg->bit&8) == 8);

    varsBT.L2  = trk->lay[1];
    varsBT.nhx = trk->num_inn_x;
    varsBT.nhy = trk->num_inn_y;

    varsBT.bt_rig = 400.0;

    varsBT.tdQ = trd->Qall;
    varsBT.tdN = trd->num_tdHit;
    for (int i = 0; i < trd->num_tdHit; ++i) {
        varsBT.tdLQ[i] = std::sqrt(trd->tdHit_amp.at(i) / trd->tdHit_len.at(i));
        varsBT.tdLL[i] = trd->tdHit_len.at(i);
    }

    varsBT.tfN = tof->Qall_nlay;
    varsBT.tfQ = tof->Qall;
    for (int i = 0; i < 4; ++i) {
        if (!tof->lay[i]) continue;
        varsBT.tfL[i]  = tof->lay[i];
        varsBT.tfLQ[i] = tof->Q[i];
    }
   
    varsBT.tkQIn = trk->QIn;
    for (int i = 0; i < 7; ++i) {
        int il = 2 + i;
        if (!trk->lay[il]) continue;
        varsBT.tkL[i]   = trk->lay[il];
        varsBT.tkLQ[i]  = trk->chrg_yj[il][2];
        varsBT.tkLQx[i] = (trk->lay[il]%2==0) ? -1 : trk->chrg_yj[il][0];
        varsBT.tkLQy[i] = (trk->lay[il]/2==0) ? -1 : trk->chrg_yj[il][1];
    }
            
    for (int i = 0; i < 9; ++i) {
        if (!trk->lay[i]) continue;
        varsBT.tkHitL[i] = trk->lay[i];
        varsBT.tkHitX[i] = trk->loc[i][0];
        varsBT.tkHitY[i] = trk->loc[i][1];
        varsBT.tkHitZ[i] = trk->loc[i][2];
    }

    for (int i = 0; i < 4; ++i) {
        if (!trk->ck_status[i]) continue;
        varsBT.ck[i]         = trk->ck_status[i];
        varsBT.ck_rig[i]     = trk->ck_rig[i];
        varsBT.ck_nchi[i][0] = trk->ck_nchi[i][0];
        varsBT.ck_nchi[i][1] = trk->ck_nchi[i][1];
        varsBT.ck_cpu[i]     = trk->ck_cpu_time[i] * 1.0e+3;
    }
    
    for (int i = 0; i < 4; ++i) {
        if (!trk->kf_status[i]) continue;
        varsBT.kf[i]         = trk->kf_status[i];
        varsBT.kf_rig[i]     = trk->kf_cen_rig[i];
        varsBT.kf_nchi[i][0] = trk->kf_nchi[i][0];
        varsBT.kf_nchi[i][1] = trk->kf_nchi[i][1];
        varsBT.kf_cpu[i]     = trk->kf_cpu_time[i] * 1.0e+3;
    }
    
    for (int i = 0; i < 4; ++i) {
        if (!hyc->geom_status[i]) continue;
        varsBT.hc[i]         = hyc->geom_status[i];
        varsBT.hc_rig[i]     = hyc->geom_cen_rig[i];
        varsBT.hc_nchi[i][0] = hyc->geom_nchi_x[i];
        varsBT.hc_nchi[i][1] = hyc->geom_nchi_y[i];
        varsBT.hc_cpu[i]     = hyc->geom_cpu_time[i] * 1.0e+3;
    }
    
    // New Track Fitting
    TrSys::Tracker tracker;
    for (int il = 0; il < 9; ++il) {
        if (!varsBT.tkHitL[il]) continue;
        TrSys::TrackerHit hit(
            il+1, varsBT.tkHitL[il]%2==1, varsBT.tkHitL[il]/2==1, 
            {  varsBT.tkHitX[il], varsBT.tkHitY[il], varsBT.tkHitZ[il] });
        tracker.add_hit(hit);
    }
    
    TrSys::GeomTrFit geom_fit_in(tracker.get_hits_with_l(TrSys::Tracker::Pattern::kInn));
    TrSys::GeomTrFit geom_fit_l1(tracker.get_hits_with_l(TrSys::Tracker::Pattern::kInnL1));
    TrSys::GeomTrFit geom_fit_l9(tracker.get_hits_with_l(TrSys::Tracker::Pattern::kInnL9));
    TrSys::GeomTrFit geom_fit_fs(tracker.get_hits_with_l(TrSys::Tracker::Pattern::kFullSpan));

    if (geom_fit_in.status()) {
        varsBT.new_hc[0]         = geom_fit_in.status();
        varsBT.new_hc_rig[0]     = geom_fit_in.part().rig();
        varsBT.new_hc_nchi[0][0] = geom_fit_in.nchi_x();
        varsBT.new_hc_nchi[0][1] = geom_fit_in.nchi_y();
    }
    if (geom_fit_l1.status()) {
        varsBT.new_hc[1]         = geom_fit_l1.status();
        varsBT.new_hc_rig[1]     = geom_fit_l1.part().rig();
        varsBT.new_hc_nchi[1][0] = geom_fit_l1.nchi_x();
        varsBT.new_hc_nchi[1][1] = geom_fit_l1.nchi_y();
    }
    if (geom_fit_l9.status()) {
        varsBT.new_hc[2]         = geom_fit_l9.status();
        varsBT.new_hc_rig[2]     = geom_fit_l9.part().rig();
        varsBT.new_hc_nchi[2][0] = geom_fit_l9.nchi_x();
        varsBT.new_hc_nchi[2][1] = geom_fit_l9.nchi_y();
    }
    if (geom_fit_fs.status()) {
        varsBT.new_hc[3]         = geom_fit_fs.status();
        varsBT.new_hc_rig[3]     = geom_fit_fs.part().rig();
        varsBT.new_hc_nchi[3][0] = geom_fit_fs.nchi_x();
        varsBT.new_hc_nchi[3][1] = geom_fit_fs.nchi_y();
    }

    // TRD
    std::vector<double> tdHits;
    for (int it = 0; it < varsBT.tdN; ++it) {
        tdHits.push_back(varsBT.tdLQ[it] * varsBT.tdLQ[it]);
    }
    if (tdHits.size() > 1) std::sort(tdHits.begin(), tdHits.end());

    const std::array<double, 3> tdPARS({ 
        0.00099201, 1.51589084, 0.67709517 
    });
    TrSys::LandauGaus tdALG(TrSys::Robust(4.0, 0.5), tdPARS[0], 0.0, tdPARS[2]);
    TrSys::SimpleALGFit tdALGFit(tdHits, tdALG);
    
    if (tdALGFit.status()) varsBT.new_tdQ      = std::sqrt(tdALGFit.param()) * 1.068630;
    if (tdALGFit.status()) varsBT.new_tdQ_nchi = tdALGFit.nchi();
    
    // TOF
    std::vector<double> tfHits;
    for (int it = 0; it < 4; ++it) {
        if (!varsBT.tfL[it]) continue;
        tfHits.push_back(varsBT.tfLQ[it] * varsBT.tfLQ[it]);
    }
    if (tfHits.size() > 1) std::sort(tfHits.begin(), tfHits.end());
    
    const std::array<double, 4> tfPARS({ 
        0.00111238, 1.04406752, 0.10728318, 0.07851389 
    });
    TrSys::LandauGaus tfALG(TrSys::Robust(4.0, 0.5), tfPARS[0], 0.0, tfPARS[2], tfPARS[3]);
    TrSys::SimpleALGFit tfALGFit(tfHits, tfALG);
    
    if (tfALGFit.status()) varsBT.new_tfQ      = std::sqrt(tfALGFit.param()) * 0.987062;
    if (tfALGFit.status()) varsBT.new_tfQ_nchi = tfALGFit.nchi();
   
    // Tracker
    std::vector<double> tkHits;
    for (int it = 0; it < 7; ++it) {
        if (!varsBT.tkL[it]) continue;
        if (varsBT.tkLQx[it] > 0) tkHits.push_back(varsBT.tkLQx[it] * varsBT.tkLQx[it]);
        if (varsBT.tkLQy[it] > 0) tkHits.push_back(varsBT.tkLQy[it] * varsBT.tkLQy[it]);
    }
    if (tkHits.size() > 1) std::sort(tkHits.begin(), tkHits.end());
    
    const std::array<double, 4> tkPARS({ 
        0.00073263, 0.92954739, 0.08338253, 0.15640932 
    });
    TrSys::LandauGaus tkALG(TrSys::Robust(4.0, 0.5), tkPARS[0], 0.0, tkPARS[2], tkPARS[3]);
    TrSys::SimpleALGFit tkALGFit(tkHits, tkALG);
   
    if (tkALGFit.status()) varsBT.new_tkQ      = std::sqrt(tkALGFit.param()) * 1.021950;
    if (tkALGFit.status()) varsBT.new_tkQ_nchi = tkALGFit.nchi();
        
    return true;
}

#endif // __Selector_C__
