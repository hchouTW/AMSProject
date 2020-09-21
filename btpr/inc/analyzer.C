#ifndef __Analyzer_C__
#define __Analyzer_C__

#include "analyzer.h"

using namespace MGROOT;
        
static constexpr double Mproton    = 0.938272297;
static constexpr double Mdeuterium = 1.876123915;


void Analyzer::set_environment() {
    std::cout << Format("\n====  Set Environment ====\n");
    LOG(INFO) << Format("\n====  Set Environment ====\n");
    
    //TrSys::PhysEnv::ReadMagAMS("/eos/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    //TrSys::PhysEnv::ReadMatAMS("/eos/user/h/hchou/ExternalLibs/DB/material");
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
    
    std::cout << Format("\n----==== BT ====----\n");
    LOG(INFO) << Format("\n----==== BT ====----\n");
    process_data(varsBT.get_chain(), &Analyzer::process_data_bt); 
    
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

    Axis AXRrso("R_{BT}/R_{Meas}-1", 2000, -10.0, 10.0);
    Axis AXQtd("TRD Q", 1500, 0.0, 3.0);
    Axis AXQtf("TOF Q", 1500, 0.0, 3.0);
    Axis AXQtk("TRK Q", 1500, 0.0, 3.0);
    Axis AXEtd("TRD dE/dx", 2000, 0.0, 60);
    Axis AXEtf("TOF dE/dx", 1000, 0.0, 20);
    Axis AXEtk("TRK dE/dx", 1000, 0.0, 10);
    Axis AXlchi("Log(#chi^{2}/NDF)", 200, -4, 6);

    Hist::New("hRrso_ckIn", HistAxis(AXRrso));
    Hist::New("hRrso_ckL1", HistAxis(AXRrso));
    Hist::New("hRrso_ckL9", HistAxis(AXRrso));
    Hist::New("hRrso_ckFs", HistAxis(AXRrso));
    Hist::New("hRrso_kfIn", HistAxis(AXRrso));
    Hist::New("hRrso_kfL1", HistAxis(AXRrso));
    Hist::New("hRrso_kfL9", HistAxis(AXRrso));
    Hist::New("hRrso_kfFs", HistAxis(AXRrso));
    Hist::New("hRrso_hcIn", HistAxis(AXRrso));
    Hist::New("hRrso_hcL1", HistAxis(AXRrso));
    Hist::New("hRrso_hcL9", HistAxis(AXRrso));
    Hist::New("hRrso_hcFs", HistAxis(AXRrso));
    
    Hist::New("hRrso_hcIn_new", HistAxis(AXRrso));
    Hist::New("hRrso_hcL1_new", HistAxis(AXRrso));
    Hist::New("hRrso_hcL9_new", HistAxis(AXRrso));
    Hist::New("hRrso_hcFs_new", HistAxis(AXRrso));
    
    Hist::New("hEtd", HistAxis(AXEtd));
    Hist::New("hEtf", HistAxis(AXEtf));
    Hist::New("hEtk", HistAxis(AXEtk));
    Hist::New("hEtkx", HistAxis(AXEtk));
    Hist::New("hEtky", HistAxis(AXEtk));
    Hist::New("hEtkxy", HistAxis(AXEtk));

    Hist::New("hQtd_off", HistAxis(AXQtd));
    Hist::New("hQtf_off", HistAxis(AXQtf));
    Hist::New("hQtk_off", HistAxis(AXQtk));
    
    Hist::New("hQtd_alg", HistAxis(AXQtd));
    Hist::New("hQtf_alg", HistAxis(AXQtf));
    Hist::New("hQtk_alg", HistAxis(AXQtk));
    
    Hist::New("hQtd_comp", HistAxis(AXQtd, AXQtd));
    Hist::New("hQtf_comp", HistAxis(AXQtf, AXQtf));
    Hist::New("hQtk_comp", HistAxis(AXQtk, AXQtk));

    Hist::New("hQtd_alg_good", HistAxis(AXQtd));
    Hist::New("hQtf_alg_good", HistAxis(AXQtf));
    Hist::New("hQtk_alg_good", HistAxis(AXQtk));
    
    Hist::New("hQtd_lchi", HistAxis(AXQtd, AXlchi));
    Hist::New("hQtf_lchi", HistAxis(AXQtf, AXlchi));
    Hist::New("hQtk_lchi", HistAxis(AXQtk, AXlchi));
    
    Hist::New("hQtktf_off", HistAxis(AXQtk, AXQtf));
    Hist::New("hQtktf_alg", HistAxis(AXQtk, AXQtf));
    Hist::New("hQtktf_lchi", HistAxis(AXlchi, AXlchi));
    
    return true;
}


bool Analyzer::process_data_bt() {
    VarsBT* vars = &varsBT;
    
    if (vars->ck[0]) Hist::Head("hRrso_ckIn")->fillH1D(vars->bt_rig / vars->ck_rig[0] - 1.0);
    if (vars->ck[1]) Hist::Head("hRrso_ckL1")->fillH1D(vars->bt_rig / vars->ck_rig[1] - 1.0);
    if (vars->ck[2]) Hist::Head("hRrso_ckL9")->fillH1D(vars->bt_rig / vars->ck_rig[2] - 1.0);
    if (vars->ck[3]) Hist::Head("hRrso_ckFs")->fillH1D(vars->bt_rig / vars->ck_rig[3] - 1.0);
    
    if (vars->kf[0]) Hist::Head("hRrso_kfIn")->fillH1D(vars->bt_rig / vars->kf_rig[0] - 1.0);
    if (vars->kf[1]) Hist::Head("hRrso_kfL1")->fillH1D(vars->bt_rig / vars->kf_rig[1] - 1.0);
    if (vars->kf[2]) Hist::Head("hRrso_kfL9")->fillH1D(vars->bt_rig / vars->kf_rig[2] - 1.0);
    if (vars->kf[3]) Hist::Head("hRrso_kfFs")->fillH1D(vars->bt_rig / vars->kf_rig[3] - 1.0);
    
    if (vars->hc[0]) Hist::Head("hRrso_hcIn")->fillH1D(vars->bt_rig / vars->hc_rig[0] - 1.0);
    if (vars->hc[1]) Hist::Head("hRrso_hcL1")->fillH1D(vars->bt_rig / vars->hc_rig[1] - 1.0);
    if (vars->hc[2]) Hist::Head("hRrso_hcL9")->fillH1D(vars->bt_rig / vars->hc_rig[2] - 1.0);
    if (vars->hc[3]) Hist::Head("hRrso_hcFs")->fillH1D(vars->bt_rig / vars->hc_rig[3] - 1.0);
  
    //TrSys::Tracker tracker;
    //for (int il = 0; il < 9; ++il) {
    //    if (!vars->tkHitL[il]) continue;
    //    TrSys::TrackerHit hit(il+1, vars->tkHitL[il]%2==1, vars->tkHitL[il]/2==1, {  vars->tkHitX[il], vars->tkHitY[il], vars->tkHitZ[il] });
    //    tracker.add_hit(hit);
    //}
    //
    //TrSys::GeomTrFit geom_fit(tracker.get_hits_with_l(TrSys::Tracker::Pattern::kInn), TrSys::PartList::kProton, false, false);
    //if (geom_fit.status()) Hist::Head("hRrso_hcTT")->fillH1D(vars->bt_rig / geom_fit.part().rig() - 1.0);
    
    if (vars->new_hc[0]) Hist::Head("hRrso_hcIn_new")->fillH1D(vars->bt_rig / vars->new_hc_rig[0] - 1.0);
    if (vars->new_hc[1]) Hist::Head("hRrso_hcL1_new")->fillH1D(vars->bt_rig / vars->new_hc_rig[1] - 1.0);
    if (vars->new_hc[2]) Hist::Head("hRrso_hcL9_new")->fillH1D(vars->bt_rig / vars->new_hc_rig[2] - 1.0);
    if (vars->new_hc[3]) Hist::Head("hRrso_hcFs_new")->fillH1D(vars->bt_rig / vars->new_hc_rig[3] - 1.0);

    for (int it = 0; it < vars->tdN; ++it) {
        Hist::Head("hEtd")->fillH1D(vars->tdLQ[it] * vars->tdLQ[it]);
    }
    
    for (int it = 0; it < 4; ++it) {
        if (!vars->tfL[it]) continue;
        Hist::Head("hEtf")->fillH1D(vars->tfLQ[it] * vars->tfLQ[it]);
    }
    
    for (int it = 0; it < 7; ++it) {
        if (!vars->tkL[it]) continue;
        if (vars->tkLQ[it]  > 0) Hist::Head("hEtk") ->fillH1D(vars->tkLQ[it]  * vars->tkLQ[it]);
        if (vars->tkLQx[it] > 0) Hist::Head("hEtkx")->fillH1D(vars->tkLQx[it] * vars->tkLQx[it]);
        if (vars->tkLQy[it] > 0) Hist::Head("hEtky")->fillH1D(vars->tkLQy[it] * vars->tkLQy[it]);
        if (vars->tkLQx[it] > 0) Hist::Head("hEtkxy")->fillH1D(vars->tkLQx[it] * vars->tkLQx[it]);
        if (vars->tkLQy[it] > 0) Hist::Head("hEtkxy")->fillH1D(vars->tkLQy[it] * vars->tkLQy[it]);
    }
    
    if (vars->tdQ   > 0) Hist::Head("hQtd_off")->fillH1D(vars->tdQ  );
    if (vars->tfQ   > 0) Hist::Head("hQtf_off")->fillH1D(vars->tfQ  );
    if (vars->tkQIn > 0) Hist::Head("hQtk_off")->fillH1D(vars->tkQIn);
    
    if (vars->tkQIn > 0 && vars->tfQ > 0) Hist::Head("hQtktf_off")->fillH2D(vars->tkQIn, vars->tfQ);
    
    if (vars->new_tdQ > 0) Hist::Head("hQtd_alg")->fillH1D(vars->new_tdQ);
    if (vars->new_tfQ > 0) Hist::Head("hQtf_alg")->fillH1D(vars->new_tfQ);
    if (vars->new_tkQ > 0) Hist::Head("hQtk_alg")->fillH1D(vars->new_tkQ);
    
    if (vars->tdQ   > 0 && vars->new_tdQ > 0) Hist::Head("hQtd_comp")->fillH2D(vars->new_tdQ, vars->tdQ  );
    if (vars->tfQ   > 0 && vars->new_tfQ > 0) Hist::Head("hQtf_comp")->fillH2D(vars->new_tfQ, vars->tfQ  );
    if (vars->tkQIn > 0 && vars->new_tkQ > 0) Hist::Head("hQtk_comp")->fillH2D(vars->new_tkQ, vars->tkQIn);
    
    if (vars->new_tdQ > 0 && std::log(vars->new_tdQ_nchi) < 1.6) Hist::Head("hQtd_alg_good")->fillH1D(vars->new_tdQ);
    if (vars->new_tfQ > 0 && std::log(vars->new_tfQ_nchi) < 1.6) Hist::Head("hQtf_alg_good")->fillH1D(vars->new_tfQ);
    if (vars->new_tkQ > 0 && std::log(vars->new_tkQ_nchi) < 1.6) Hist::Head("hQtk_alg_good")->fillH1D(vars->new_tkQ);
    
    if (vars->new_tdQ > 0) Hist::Head("hQtd_lchi")->fillH2D(vars->new_tdQ, std::log(vars->new_tdQ_nchi));
    if (vars->new_tfQ > 0) Hist::Head("hQtf_lchi")->fillH2D(vars->new_tfQ, std::log(vars->new_tfQ_nchi));
    if (vars->new_tkQ > 0) Hist::Head("hQtk_lchi")->fillH2D(vars->new_tkQ, std::log(vars->new_tkQ_nchi));
    
    if (vars->new_tfQ > 0 && vars->new_tkQ > 0) Hist::Head("hQtktf_alg" )->fillH2D(vars->new_tkQ, vars->new_tfQ);
    if (vars->new_tfQ > 0 && vars->new_tkQ > 0) Hist::Head("hQtktf_lchi")->fillH2D(std::log(vars->new_tkQ_nchi), std::log(vars->new_tfQ_nchi));

    return true;
}


#endif // __Analyzer_C__
