#ifndef __Selector_C__
#define __Selector_C__

#include "selector.h"

using namespace MGROOT;

void Selector::set_environment() {
    std::cout << Format("\n====  Set Environment ====\n");
    LOG(INFO) << Format("\n====  Set Environment ====\n");
    
    TrSys::PhysEnv::ReadMagAMS("/eos/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    TrSys::PhysEnv::ReadMatAMS("/eos/user/h/hchou/ExternalLibs/DB/material");
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
    vars.init();

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

    // testcode
    if (tof->num_beta == 0) return false;
    if (trk->num_track == 0) return false;
    
    //if (tof->num_beta != 1) return false;
    //if (trk->num_track != 1) return false;

    if (std::log(trk->ck_nchi[0][0]) > 2.0) return false;
    if (std::log(trk->ck_nchi[0][1]) > 2.0) return false;
   
    //std::cerr << Form("%14.8f %14.8f\n", tof->nchi_c, tof->nchi_t);
    if (std::log(tof->nchi_c) > 2.0) return false;
    if (std::log(tof->nchi_t) > 2.0) return false;

    if (std::abs(tof->beta - 1.0) > 0.03) return false;
    
    //if (tof->Qall < 5) return false;
    //if (trk->QIn  < 5) return false;

    /* 
    // Charge
    //short Zall = tof->Zall;
    //if (tof->Qall < 0.8 || tof->Qall > 1.4) return false;
    //if (trk->QIn < 0.8 || trk->QIn > 1.4) return false;
   
    //if (hyc->tfQup > 2.0) return false;
    //if (hyc->tfQlw > 2.0) return false;
*/
    // TOF
    if (tof->num_extcls[0] > 1 || tof->num_extcls[1] > 1) return false;
    if (tof->num_extcls[2] > 1 || tof->num_extcls[3] > 1) return false;
    

    // Tracker
    if (trk->lay[1] == 0) return false;
   
   
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
    bool pass = process_data_db();

    if (pass) vars.fill();
    return pass;

    return true;
}

bool Selector::process_data_db() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_rig = (g4mc != nullptr) ? g4mc->prm_mom/g4mc->prm_chrg : 0.0;
    double mc_igb = (g4mc != nullptr) ? g4mc->prm_mass/g4mc->prm_mom : 0.0;
    double mc_w10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_w27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_w10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_w27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // TRD
    //if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    //if (trd->tdLLR_ep < 0.75) return false; 

    //// TOF
    //if (tof->extcls_noise != 0) return false;
    //if (tof->num_in_time_cls > 4) return false;
   
    if (!trk->ck_status[0]) return false;
    if (CheckType(Type::ISS)) {
        double cfr = std::abs(trk->ck_rig[0] / rti->max_IGRF);
        if (cfr < 1.2) return false;
    }

    if (std::abs(trk->ck_rig[0]) < 50.0) return false;

    // Variables
    vars.init();
    vars.run = list->run;
    vars.evt = list->event;
    vars.ut  = list->utime;
    vars.wgt = list->weight;

    vars.mc     = (g4mc != nullptr);
    vars.mc_z   = (g4mc != nullptr) ? g4mc->prm_chrg : 0;
    vars.mc_rig = mc_rig;
    vars.mc_igb = mc_igb;
    vars.mc_w10 = mc_w10;
    vars.mc_w27 = mc_w27;

    vars.cfsec = CheckType(Type::ISS) ? rti->max_IGRF : 0;
    vars.cfevt = CheckType(Type::ISS) ? hyc->max_IGRF : 0;
    vars.lv    = CheckType(Type::ISS) ? rti->livetime : 1;
    
    vars.trg   = ((trg->bit&8) == 8);

    vars.L2  = trk->lay[1];
    vars.nhx = trk->num_inn_x;
    vars.nhy = trk->num_inn_y;
    
    for (int i = 0; i < 4; ++i) {
        if (!trk->ck_status[i]) continue;
        vars.ck[i]         = trk->ck_status[i];
        vars.ck_rig[i]     = trk->ck_rig[i];
        vars.ck_nchi[i][0] = trk->ck_nchi[i][0];
        vars.ck_nchi[i][1] = trk->ck_nchi[i][1];
    }

    vars.tfN = tof->Qall_nlay;
    vars.tfQ = tof->Qall;
    for (int i = 0; i < 4; ++i) {
        if (!tof->lay[i]) continue;
        vars.tfL[i]  = tof->lay[i];
        vars.tfLQ[i] = tof->Q[i];
    }
    vars.tfQ_nchi_c = std::log(tof->nchi_c);
    vars.tfQ_nchi_t = std::log(tof->nchi_t);

    vars.tkQhl = trk->QIn_hl;
    vars.tkQyj = trk->QIn_yj;
   
    vars.tkQ = trk->QIn;
    for (int i = 0; i < 9; ++i) {
        if (!trk->lay[i]) continue;
        vars.tkL[i]   = trk->lay[i];
        vars.tkLQ[i]  = trk->chrg_yj[i][2];
        vars.tkLQx[i] = (trk->lay[i]%2==0) ? -1 : trk->chrg_yj[i][0];
        vars.tkLQy[i] = (trk->lay[i]/2==0) ? -1 : trk->chrg_yj[i][1];
    }
            
    for (int i = 0; i < 9; ++i) {
        if (!trk->lay[i]) continue;
        vars.tkHitL[i] = trk->lay[i];
        vars.tkHitX[i] = trk->loc[i][0];
        vars.tkHitY[i] = trk->loc[i][1];
        vars.tkHitZ[i] = trk->loc[i][2];
    }
    
    std::vector<TrSys::TofZHit> tfHits;
    for (int it = 0; it < 4; ++it) {
        if (!vars.tfL[it]) continue;
        TrSys::TofZHit hit(it, vars.tfLQ[it]);
        if (!hit.status()) continue;
        tfHits.push_back(hit);
    }
    
    TrSys::TofZFit tfZFit(tfHits);
    if (tfZFit.status() && tfZFit.clss().at(0).status()) {
        const TrSys::TofZCluster& cls = tfZFit.clss().at(0);
        vars.new_tfQ      = cls.param();
        vars.new_tfQ_nchi = std::log(cls.nchi());
        if (!std::isfinite(vars.new_tfQ_nchi)) vars.new_tfQ_nchi = -5.0;
        
        vars.new_tfQ_nhit = cls.hits().size();
        vars.new_tfQ_ncls = tfZFit.clss().size();

        //std::cerr << Form("%14.8f %14.8f\n", vars.new_tfQ, vars.new_tfQ_nchi);
    }

    
   
    // Tracker
    std::vector<TrSys::TrkZHit> tkHits;
    for (int it = 1; it < 8; ++it) {
        if (!vars.tkL[it]) continue;
        if (vars.tkLQx[it] > 0) {
            TrSys::TrkZHit hit(it, vars.tkLQx[it]);
            if (hit.status()) {
                tkHits.push_back(hit);
            }
        }
        if (vars.tkLQy[it] > 0) {
            TrSys::TrkZHit hit(it, vars.tkLQy[it]);
            if (hit.status()) {
                tkHits.push_back(hit);
            }
        }
    }
    
    TrSys::TrkZFit tkZFit(tkHits);
    if (tkZFit.status() && tkZFit.clss().at(0).status()) {
        const TrSys::TrkZCluster& cls = tkZFit.clss().at(0);
        vars.new_tkQ      = cls.param();
        vars.new_tkQ_nchi = std::log(cls.nchi());
        if (!std::isfinite(vars.new_tkQ_nchi)) vars.new_tkQ_nchi = -5.0;
        
        vars.new_tkQ_nhit = cls.hits().size();
        vars.new_tkQ_ncls = tkZFit.clss().size();
    }
        
    return true;
}


#endif // __Selector_C__
