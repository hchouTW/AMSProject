#ifndef __Analyzer_C__
#define __Analyzer_C__

#include "analyzer.h"

#include <CPPLibs.h>
#include <ROOTLibs.h>
//#include <TRACKSys.h>
//#include <TrSys.h>

using namespace MGROOT;

void Analyzer::set_environment() {
    std::cout << Format("\n====  Set Environment ====\n");
    LOG(INFO) << Format("\n====  Set Environment ====\n");
    
    TrSys::PhysEnv::ReadMagAMS("/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    TrSys::PhysEnv::ReadMatAMS("/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    
    if (!TrSys::PhysEnv::IMagStatus()) TrSys::PhysEnv::ReadMagAMS("/eos/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    if (!TrSys::PhysEnv::IMatStatus()) TrSys::PhysEnv::ReadMatAMS("/eos/user/h/hchou/ExternalLibs/DB/material");
    
    if (!TrSys::PhysEnv::IMagStatus()) TrSys::PhysEnv::ReadMagAMS("/afs/cern.ch/work/h/hchou/public/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    if (!TrSys::PhysEnv::IMatStatus()) TrSys::PhysEnv::ReadMatAMS("/afs/cern.ch/work/h/hchou/public/ExternalLibs/DB/material");

    std::cerr << Format("Load Mag %d Mat %d\n", TrSys::PhysEnv::IMagStatus(), TrSys::PhysEnv::IMatStatus());
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
	long print_rate =  loop_entries / 100;
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
    
        //if (ientry%2!=0) continue; // testcode
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
    mdst->GetEvent(0);
    int chrg = std::abs(g4mc->prm_chrg);

    Axis AXrig("Rigidity [GeV/c]", 75, (chrg == 1 ? 0.6 : 1.0), (chrg == 1 ? 3500.0 : 7000.0), AxisScale::kLog);

    Axis AXrigr("", 4500, -3.0, 3.0);
    Hist::New("hCK_inrh_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hCK_in_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hCK_l1_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hCK_l9_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hCK_fs_rig", HistAxis(AXrig, AXrigr));
    
    Hist::New("hCK_inrh_numt", HistAxis(AXrig));
    Hist::New("hCK_in_numt", HistAxis(AXrig));
    Hist::New("hCK_l1_numt", HistAxis(AXrig));
    Hist::New("hCK_l9_numt", HistAxis(AXrig));
    Hist::New("hCK_fs_numt", HistAxis(AXrig));
    
    Hist::New("hCK_inrh_dent", HistAxis(AXrig));
    Hist::New("hCK_in_dent", HistAxis(AXrig));
    Hist::New("hCK_l1_dent", HistAxis(AXrig));
    Hist::New("hCK_l9_dent", HistAxis(AXrig));
    Hist::New("hCK_fs_dent", HistAxis(AXrig));

    Hist::New("hKF_inrh_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hKF_in_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hKF_l1_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hKF_l9_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hKF_fs_rig", HistAxis(AXrig, AXrigr));
    
    Hist::New("hKF_inrh_numt", HistAxis(AXrig));
    Hist::New("hKF_in_numt", HistAxis(AXrig));
    Hist::New("hKF_l1_numt", HistAxis(AXrig));
    Hist::New("hKF_l9_numt", HistAxis(AXrig));
    Hist::New("hKF_fs_numt", HistAxis(AXrig));
    
    Hist::New("hKF_inrh_dent", HistAxis(AXrig));
    Hist::New("hKF_in_dent", HistAxis(AXrig));
    Hist::New("hKF_l1_dent", HistAxis(AXrig));
    Hist::New("hKF_l9_dent", HistAxis(AXrig));
    Hist::New("hKF_fs_dent", HistAxis(AXrig));
    
    Hist::New("hHC_inrh_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hHC_in_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hHC_l1_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hHC_l9_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hHC_fs_rig", HistAxis(AXrig, AXrigr));
    
    Hist::New("hHC_inrh_numt", HistAxis(AXrig));
    Hist::New("hHC_in_numt", HistAxis(AXrig));
    Hist::New("hHC_l1_numt", HistAxis(AXrig));
    Hist::New("hHC_l9_numt", HistAxis(AXrig));
    Hist::New("hHC_fs_numt", HistAxis(AXrig));
    
    Hist::New("hHC_inrh_dent", HistAxis(AXrig));
    Hist::New("hHC_in_dent", HistAxis(AXrig));
    Hist::New("hHC_l1_dent", HistAxis(AXrig));
    Hist::New("hHC_l9_dent", HistAxis(AXrig));
    Hist::New("hHC_fs_dent", HistAxis(AXrig));
    
    Hist::New("hCK_drig", HistAxis(AXrig, AXrigr));
    Hist::New("hKF_drig", HistAxis(AXrig, AXrigr));
    Hist::New("hHC_drig", HistAxis(AXrig, AXrigr));

    Axis AXTbta("", 40, 0.55, 0.99);
    Axis AXTbtar("", 1000, -0.30, 0.30);
    Hist::New("hTF_T_bta" , HistAxis(AXTbta, AXTbtar));
    Hist::New("hHC_T_bta" , HistAxis(AXTbta, AXTbtar));
    Hist::New("hHC_TQ_bta", HistAxis(AXTbta, AXTbtar));
    
    Hist::New("hHC_T_numt", HistAxis(AXTbta));
    Hist::New("hHC_T_dent", HistAxis(AXTbta));
    
    Hist::New("hHC_TQ_numt", HistAxis(AXTbta));
    Hist::New("hHC_TQ_dent", HistAxis(AXTbta));
    
    Axis AXBbta("", 40,  0.96, 0.999);
    Axis AXBbtar("", 1200, -0.02, 0.02);

    Hist::New("hRH_B_bta" , HistAxis(AXBbta, AXBbtar));
    Hist::New("hHC_B_bta" , HistAxis(AXBbta, AXBbtar));
    
    Hist::New("hHC_B_numt", HistAxis(AXBbta));
    Hist::New("hHC_B_dent", HistAxis(AXBbta));

    Axis AXsqrm("Mass^{2}/Z^{2} [(GeV/c^{2})^{2}]", 1600, -12.0 * (chrg == 1 ? 1.0 : 4.0), 12.0 * (chrg == 1 ? 1.0 : 4.0));
    Axis AXmass("Mass/Z [GeV/c^{2}]", 1600, 0.01, 4.0 * (chrg == 1 ? 1.0 : 2.0));
    Hist::New("hTF_T_sqrm", HistAxis(AXTbta, AXsqrm));
    Hist::New("hRH_B_sqrm", HistAxis(AXBbta, AXsqrm));
    
    Hist::New("hTF_T_mass", HistAxis(AXTbta, AXmass));
    Hist::New("hRH_B_mass", HistAxis(AXBbta, AXmass));
    
    Hist::New("hHCmutr_TQ_sqrm", HistAxis(AXTbta, AXsqrm));
    Hist::New("hHCmutr_B_sqrm",  HistAxis(AXBbta, AXsqrm));
    
    Hist::New("hHCmutr_TQ_mass", HistAxis(AXTbta, AXmass));
    Hist::New("hHCmutr_B_mass",  HistAxis(AXBbta, AXmass));
    
    Hist::New("hHCmutr_TQ_numt", HistAxis(AXTbta));
    Hist::New("hHCmutr_TQ_dent", HistAxis(AXTbta));
    
    Hist::New("hHCmutr_B_numt", HistAxis(AXBbta));
    Hist::New("hHCmutr_B_dent", HistAxis(AXBbta));
    
    Hist::New("hHCphys_TQ_rig", HistAxis(AXrig, AXrigr));
    Hist::New("hHCphys_B_rig" , HistAxis(AXrig, AXrigr));
    
    Hist::New("hHCphys_TQ_numt", HistAxis(AXrig));
    Hist::New("hHCphys_TQ_dent", HistAxis(AXrig));
    
    Hist::New("hHCphys_B_numt", HistAxis(AXrig));
    Hist::New("hHCphys_B_dent", HistAxis(AXrig));


    Axis AXigb("1/Gamma/Beta [1]", 75, 0.01, 1.75, AxisScale::kLog);

    Axis AXl("[cm]", 1000, -8., 8.);
    Axis AXu("[1]", 1000, -0.15, 0.15);
    Axis AXe("[MeV]", 4000, 0.1, (chrg == 1 ? 100. : 400.));

    Hist::New("hLxMC", HistAxis(AXigb, AXl));
    Hist::New("hLyMC", HistAxis(AXigb, AXl));
    Hist::New("hUxMC", HistAxis(AXigb, AXu));
    Hist::New("hUyMC", HistAxis(AXigb, AXu));
    Hist::New("hEaMC", HistAxis(AXigb, AXe));
    
    Hist::New("hLxSM", HistAxis(AXigb, AXl));
    Hist::New("hLySM", HistAxis(AXigb, AXl));
    Hist::New("hUxSM", HistAxis(AXigb, AXu));
    Hist::New("hUySM", HistAxis(AXigb, AXu));
    Hist::New("hEaSM", HistAxis(AXigb, AXe));
    
    Hist::New("hEaBB", HistAxis(AXigb, AXe));
    Hist::New("hEaLD", HistAxis(AXigb, AXe));
    
    return true;
}

bool Analyzer::process_presel() {
    // Charge
    //if (tof->Qall < 0.8 || tof->Qall > 1.3) return false;
    //if (trk->QIn < 0.8 || trk->QIn > 1.3) return false;
    //if (tof->Qall < 1.7 || tof->Qall > 2.4) return false;
    //if (trk->QIn < 1.7 || trk->QIn > 2.4) return false;

    // TOF
    //if (tof->nchi_t > 10.) return false;
    //if (tof->nchi_c > 10.) return false;

    // TRK
    if (!trk->ck_status[0]) return false;

    // TRD
    if (trd->tdLLR_num_hit < 8 && trd->tkLLR_num_hit < 8) return false;
    
    // MC
    if (!g4mc->tk[2] && !g4mc->tk[3]) return false;
    if (!g4mc->tk[4] && !g4mc->tk[5]) return false;
    if (!g4mc->tk[6] && !g4mc->tk[7]) return false;

    return true;
}

bool Analyzer::process_data() {
    process_scan();
    process_fit();
    return true;
}

bool Analyzer::process_scan() {
    if (!CheckType(Type::MC)) return false;
    std::array<std::string, 4> trpts({"in", "l1", "l9", "fs"});

    for (int ip = 0; ip < 4; ++ip) {
        if (!trk->ck_status[ip]) continue;
        double sclr = std::sqrt(g4mc->prm_mom/g4mc->prm_chrg);
        double dirg = 1.0/trk->ck_rig[ip] - g4mc->prm_chrg/g4mc->prm_mom;

        Hist::Head(Format("hCK_%s_rig", trpts[ip].c_str()))->fillH2D(g4mc->prm_mom/g4mc->prm_chrg, dirg * sclr, list->weight);
        Hist::Head(Format("hCK_%s_numt", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight * trk->ck_cpu_time[ip] * 1000.0);
        Hist::Head(Format("hCK_%s_dent", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight);

        if (ip == 0 && rich->status && rich->kind == 1 && rich->nhit >= 3 && rich->npmt >= 3) {
            Hist::Head(Format("hCK_%srh_rig", trpts[ip].c_str()))->fillH2D(g4mc->prm_mom/g4mc->prm_chrg, dirg * sclr, list->weight);
            Hist::Head(Format("hCK_%srh_numt", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight * trk->ck_cpu_time[ip] * 1000.0);
            Hist::Head(Format("hCK_%srh_dent", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight);
        }
    }
    
    for (int ip = 0; ip < 4; ++ip) {
        if (!trk->kf_status[ip]) continue;
        double sclr = std::sqrt(g4mc->prm_mom/g4mc->prm_chrg);
        double dirg = 1.0/trk->kf_top_rig[ip] - g4mc->prm_chrg/g4mc->prm_mom;

        Hist::Head(Format("hKF_%s_rig", trpts[ip].c_str()))->fillH2D(g4mc->prm_mom/g4mc->prm_chrg, dirg * sclr, list->weight);
        Hist::Head(Format("hKF_%s_numt", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight * trk->kf_cpu_time[ip] * 1000.0);
        Hist::Head(Format("hKF_%s_dent", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight);

        if (ip == 0 && rich->status && rich->kind == 1 && rich->nhit >= 3 && rich->npmt >= 3) {
            Hist::Head(Format("hKF_%srh_rig", trpts[ip].c_str()))->fillH2D(g4mc->prm_mom/g4mc->prm_chrg, dirg * sclr, list->weight);
            Hist::Head(Format("hKF_%srh_numt", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight * trk->kf_cpu_time[ip] * 1000.0);
            Hist::Head(Format("hKF_%srh_dent", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight);
        }
    }
    
    for (int ip = 0; ip < 4; ++ip) {
        if (!hyc->geom_status[ip]) continue;
        double sclr = std::sqrt(g4mc->prm_mom/g4mc->prm_chrg);
        double dirg = 1.0/hyc->geom_top_rig[ip] - g4mc->prm_chrg/g4mc->prm_mom;

        Hist::Head(Format("hHC_%s_rig", trpts[ip].c_str()))->fillH2D(g4mc->prm_mom/g4mc->prm_chrg, dirg * sclr, list->weight);
        Hist::Head(Format("hHC_%s_numt", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight * hyc->geom_cpu_time[ip] * 1000.0);
        Hist::Head(Format("hHC_%s_dent", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight);
       
        if (ip == 0 && rich->self_status && rich->self_kind == 1 && rich->self_cld_nhit >= 3 && rich->self_cld_npmt >= 3) {
            Hist::Head(Format("hHC_%srh_rig", trpts[ip].c_str()))->fillH2D(g4mc->prm_mom/g4mc->prm_chrg, dirg * sclr, list->weight);
            Hist::Head(Format("hHC_%srh_numt", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight * hyc->geom_cpu_time[ip] * 1000.0);
            Hist::Head(Format("hHC_%srh_dent", trpts[ip].c_str()))->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight);
        }
    }

    if (trk->ck_status[1] && trk->ck_status[2] && trk->ck_status[3] && trk->ck_rig[3] > 0) {
        double sclr = std::sqrt(trk->ck_rig[3]);
        double dirg = 1.0/trk->ck_rig[2] - 1.0/trk->ck_rig[1];
        Hist::Head("hCK_drig")->fillH2D(trk->ck_rig[3], dirg * sclr, list->weight);
    }
    
    if (trk->kf_status[1] && trk->kf_status[2] && trk->kf_status[3] && trk->kf_top_rig[3] > 0) {
        double sclr = std::sqrt(trk->kf_top_rig[3]);
        double dirg = 1.0/trk->kf_top_rig[2] - 1.0/trk->kf_top_rig[1];
        Hist::Head("hKF_drig")->fillH2D(trk->kf_top_rig[3], dirg * sclr, list->weight);
    }
    
    if (hyc->geom_status[1] && hyc->geom_status[2] && hyc->geom_status[3] && hyc->geom_top_rig[3] > 0) {
        double sclr = std::sqrt(hyc->geom_top_rig[3]);
        double dirg = 1.0/hyc->geom_top_rig[2] - 1.0/hyc->geom_top_rig[1];
        Hist::Head("hHC_drig")->fillH2D(hyc->geom_top_rig[3], dirg * sclr, list->weight);
    }
    
    if (tof->status) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double dbta = (tof->beta - mbta) / mbta;
        
        Hist::Head("hTF_T_bta")->fillH2D(mbta, dbta, list->weight);
    }
    
    if (rich->status && rich->kind == 1 && rich->nhit >= 3 && rich->npmt >= 3) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double dbta = (rich->beta - mbta) / mbta;
        
        Hist::Head("hRH_B_bta") ->fillH2D(mbta, dbta, list->weight);
    }

    if (hyc->vel_status[0]) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double dbta = (hyc->vel_top_bta[0] - mbta) / mbta;

        Hist::Head("hHC_T_bta")->fillH2D(mbta, dbta, list->weight);
        Hist::Head("hHC_T_numt")->fillH1D(mbta, list->weight * hyc->vel_cpu_time[0] * 1000.0);
        Hist::Head("hHC_T_dent")->fillH1D(mbta, list->weight);
    }
    
    if (hyc->vel_status[1]) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double dbta = (hyc->vel_top_bta[1] - mbta) / mbta;

        Hist::Head("hHC_TQ_bta")->fillH2D(mbta, dbta, list->weight);
        Hist::Head("hHC_TQ_numt")->fillH1D(mbta, list->weight * hyc->vel_cpu_time[1] * 1000.0);
        Hist::Head("hHC_TQ_dent")->fillH1D(mbta, list->weight);
    }
    
    if (rich->self_status && rich->self_kind == 1 && rich->self_cld_nhit >= 3 && rich->self_cld_npmt >= 3 && hyc->vel_status[2]) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double dbta = (hyc->vel_top_bta[2] - mbta) / mbta;

        Hist::Head("hHC_B_bta") ->fillH2D(mbta, dbta, list->weight);
        Hist::Head("hHC_B_numt")->fillH1D(mbta, list->weight * hyc->vel_cpu_time[2] * 1000.0);
        Hist::Head("hHC_B_dent")->fillH1D(mbta, list->weight);
    }

    if (tof->status && trk->ck_status[0]) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double sqrm = (trk->ck_rig[0] * (1.0/tof->beta + 1.0)) * (trk->ck_rig[0] * (1.0/tof->beta - 1.0));
        double mass = (tof->beta < 1.0) ? std::sqrt(sqrm) : 0.0;
        
        Hist::Head("hTF_T_sqrm")->fillH2D(mbta, sqrm, list->weight);
        Hist::Head("hTF_T_mass")->fillH2D(mbta, mass, list->weight);
    }
    
    if (rich->status && rich->kind == 1 && rich->nhit >= 3 && rich->npmt >= 3 && trk->ck_status[0]) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double sqrm = (trk->ck_rig[0] * (1.0/rich->beta + 1.0)) * (trk->ck_rig[0] * (1.0/rich->beta - 1.0));
        double mass = (rich->beta < 1.0) ? std::sqrt(sqrm) : 0.0;
        
        Hist::Head("hRH_B_sqrm")->fillH2D(mbta, sqrm, list->weight);
        Hist::Head("hRH_B_mass")->fillH2D(mbta, mass, list->weight);
    }
    
    if (hyc->mutr_status[0]) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double sqrm = hyc->mutr_sqrm[0] / g4mc->prm_chrg / g4mc->prm_chrg;
        double mass = (hyc->mutr_sqrm[0] > 0) ? std::sqrt(sqrm) : 0.0;

        Hist::Head("hHCmutr_TQ_sqrm")->fillH2D(mbta, sqrm, list->weight);
        Hist::Head("hHCmutr_TQ_mass")->fillH2D(mbta, mass, list->weight);
        Hist::Head("hHCmutr_TQ_numt")->fillH1D(mbta, list->weight * hyc->mutr_cpu_time[0] * 1000.0);
        Hist::Head("hHCmutr_TQ_dent")->fillH1D(mbta, list->weight);
    }
    
    if (rich->self_status && rich->self_kind == 1 && rich->self_cld_nhit >= 3 && rich->self_cld_npmt >= 3 && hyc->mutr_status[1]) {
        double mbta = 1.0/std::sqrt((g4mc->prm_mass/g4mc->prm_mom)*(g4mc->prm_mass/g4mc->prm_mom)+1.0);
        double sqrm = hyc->mutr_sqrm[1] / g4mc->prm_chrg / g4mc->prm_chrg;
        double mass = (hyc->mutr_sqrm[1] > 0) ? std::sqrt(sqrm) : 0.0;

        Hist::Head("hHCmutr_B_sqrm")->fillH2D(mbta, sqrm, list->weight);
        Hist::Head("hHCmutr_B_mass")->fillH2D(mbta, mass, list->weight);
        Hist::Head("hHCmutr_B_numt")->fillH1D(mbta, list->weight * hyc->mutr_cpu_time[1] * 1000.0);
        Hist::Head("hHCmutr_B_dent")->fillH1D(mbta, list->weight);
    }
    
    if (hyc->phys_status[0][0]) {
        double sclr = std::sqrt(g4mc->prm_mom/g4mc->prm_chrg);
        double dirg = 1.0/hyc->phys_top_rig[0][0] - g4mc->prm_chrg/g4mc->prm_mom;
        
        Hist::Head("hHCphys_TQ_rig")->fillH2D(g4mc->prm_mom/g4mc->prm_chrg, dirg * sclr, list->weight);
        Hist::Head("hHCphys_TQ_numt")->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight * hyc->phys_cpu_time[0][0] * 1000.0);
        Hist::Head("hHCphys_TQ_dent")->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight);
    }
    
    if (rich->self_status && rich->self_kind == 1 && rich->self_cld_nhit >= 3 && rich->self_cld_npmt >= 3 && hyc->phys_status[0][1]) {
        double sclr = std::sqrt(g4mc->prm_mom/g4mc->prm_chrg);
        double dirg = 1.0/hyc->phys_top_rig[0][1] - g4mc->prm_chrg/g4mc->prm_mom;
        
        Hist::Head("hHCphys_B_rig") ->fillH2D(g4mc->prm_mom/g4mc->prm_chrg, dirg * sclr, list->weight);
        Hist::Head("hHCphys_B_numt")->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight * hyc->phys_cpu_time[0][1] * 1000.0);
        Hist::Head("hHCphys_B_dent")->fillH1D(g4mc->prm_mom/g4mc->prm_chrg, list->weight);
    }

    return false;
}

bool Analyzer::process_fit() {
    if (!CheckType(Type::MC)) return false;
    
    int lsat = 0;
    int lend = 1;
    if (g4mc->tk[lsat] && g4mc->tk[lend]) {
        bool mscat = true;
        bool eloss = true;

        TrSys::Part org_part(TrSys::PartType(g4mc->prm_chrg, g4mc->prm_mass));
        org_part.set_location (g4mc->tk_loc[lsat][0], g4mc->tk_loc[lsat][1], g4mc->tk_loc[lsat][2]);
        org_part.set_direction(g4mc->tk_dir[lsat][0], g4mc->tk_dir[lsat][1], g4mc->tk_dir[lsat][2]);
        org_part.set_mom(g4mc->tk_mom[lsat]);

        TrSys::Part mc_part(org_part);
        auto&& mc_rlt = TrSys::Prop::DoPropToZ(g4mc->tk_loc[lend][2], mc_part, TrSys::PropArgs(mscat, eloss));

        double mc_lx = g4mc->tk_loc[lend][0] - mc_part.lx();
        double mc_ly = g4mc->tk_loc[lend][1] - mc_part.ly();
        double mc_ux = g4mc->tk_dir[lend][0] - mc_part.ux();
        double mc_uy = g4mc->tk_dir[lend][1] - mc_part.uy();
        double mc_ea = (std::hypot(g4mc->prm_mass, g4mc->tk_mom[lsat]) - std::hypot(g4mc->prm_mass, g4mc->tk_mom[lend])) * 1000.;

        TrSys::Part sm_part(org_part);
        auto sm_rlt = TrSys::Prop::DoPropToZ(g4mc->tk_loc[lend][2], sm_part, TrSys::PropArgs(mscat, eloss), TrSys::InteractionArgs(true));
        
        double sm_lx = sm_part.lx() - mc_part.lx();
        double sm_ly = sm_part.ly() - mc_part.ly();
        double sm_ux = sm_part.ux() - mc_part.ux();
        double sm_uy = sm_part.uy() - mc_part.uy();
        double sm_ea = (std::hypot(g4mc->prm_mass, g4mc->tk_mom[lsat]) - sm_part.eng()) * 1000.;
        
        TrSys::Part bb_part(org_part);
        TrSys::Ionization::kDefaultFormula = TrSys::Ionization::Formula::kBetheBloch;
        auto bb_rlt = TrSys::Prop::DoPropToZ(g4mc->tk_loc[lend][2], bb_part, TrSys::PropArgs(mscat, eloss), TrSys::InteractionArgs(true));
        double bb_ea = (std::hypot(g4mc->prm_mass, g4mc->tk_mom[lsat]) - bb_part.eng()) * 1000.;
   
        TrSys::Part ld_part(org_part);
        TrSys::Ionization::kDefaultFormula = TrSys::Ionization::Formula::kLandau;
        auto ld_rlt = TrSys::Prop::DoPropToZ(g4mc->tk_loc[lend][2], ld_part, TrSys::PropArgs(mscat, eloss), TrSys::InteractionArgs(true));
        double ld_ea = (std::hypot(g4mc->prm_mass, g4mc->tk_mom[lsat]) - ld_part.eng()) * 1000.;
        
        TrSys::Ionization::kDefaultFormula = TrSys::Ionization::Formula::kMixture;
 
        double igb = g4mc->prm_mass/g4mc->tk_mom[lsat];
        double scale_lu = (g4mc->tk_mom[lsat] / g4mc->prm_mass);
        double scale_ee = 1.0 / (1.0 + (g4mc->prm_mass / g4mc->tk_mom[lsat]) * (g4mc->prm_mass / g4mc->tk_mom[lsat]));
        
        Hist::Head("hLxMC")->fillH2D(igb, scale_lu * mc_lx, list->weight);
        Hist::Head("hLyMC")->fillH2D(igb, scale_lu * mc_ly, list->weight);
        Hist::Head("hUxMC")->fillH2D(igb, scale_lu * mc_ux, list->weight);
        Hist::Head("hUyMC")->fillH2D(igb, scale_lu * mc_uy, list->weight);
        Hist::Head("hEaMC")->fillH2D(igb, scale_ee * mc_ea, list->weight);
        
        Hist::Head("hLxSM")->fillH2D(igb, scale_lu * sm_lx, list->weight);
        Hist::Head("hLySM")->fillH2D(igb, scale_lu * sm_ly, list->weight);
        Hist::Head("hUxSM")->fillH2D(igb, scale_lu * sm_ux, list->weight);
        Hist::Head("hUySM")->fillH2D(igb, scale_lu * sm_uy, list->weight);
        Hist::Head("hEaSM")->fillH2D(igb, scale_ee * sm_ea, list->weight);
        
        Hist::Head("hEaBB")->fillH2D(igb, scale_ee * bb_ea, list->weight);
        Hist::Head("hEaLD")->fillH2D(igb, scale_ee * ld_ea, list->weight);
    }



    /* 
    int lcen = -1;
    if (lcen < 0 && g4mc->tk[4]) lcen = 4;
    if (lcen < 0 && g4mc->tk[5]) lcen = 5;
    if (lcen < 0) return false;

    int linn = -1;
    if (linn < 0 && g4mc->tk[1]) linn = 1;
    if (linn < 0 && g4mc->tk[2]) linn = 2;
    if (linn < 0 && g4mc->tk[3]) linn = 3;
    if (linn < 0) return false;

    for (int il = 0; il < 20; ++il) {
        if (trd->num_tkHit < 8) continue;
        if (!g4mc->tk[il]) continue;

        double igbta = g4mc->prm_mass / g4mc->td_mom[il];
        //double ibta  = std::sqrt(igbta * igbta + 1.0);
        for (int jl = 0; jl < trd->num_tkHit; ++jl) {
            if (trd->tkHit_lay[jl] != il) continue;
            if (trd->tkHit_len[jl] < 0.3) continue;
            double dEdx = trd->tkHit_amp[jl] / trd->tkHit_len[jl];
            Hist::Head("hTd")->fillH2D(igbta, dEdx, list->weight);
            Hist::Head("hTdL")->fillH2D(igbta, std::log(dEdx), list->weight);
        }
    }

    int lsat = 0;
    int lend = 1;
    if (g4mc->tk[lsat] && g4mc->tk[lend]) {
        bool mscat = true;
        bool eloss = true;

        TrSys::Part org_part(TrSys::PartList::kProton);
        org_part.set_location (g4mc->tk_loc[lsat][0], g4mc->tk_loc[lsat][1], g4mc->tk_loc[lsat][2]);
        org_part.set_direction(g4mc->tk_dir[lsat][0], g4mc->tk_dir[lsat][1], g4mc->tk_dir[lsat][2]);
        org_part.set_mom(g4mc->tk_mom[lsat]);

        TrSys::Part mc_part(org_part);
        auto&& mc_rlt = TrSys::Prop::DoPropToZ(g4mc->tk_loc[lend][2], mc_part, TrSys::PropArgs(mscat, eloss));

        double mc_lx = g4mc->tk_loc[lend][0] - mc_part.lx();
        double mc_ly = g4mc->tk_loc[lend][1] - mc_part.ly();
        double mc_ux = g4mc->tk_dir[lend][0] - mc_part.ux();
        double mc_uy = g4mc->tk_dir[lend][1] - mc_part.uy();
        double mc_ea = (std::hypot(g4mc->prm_mass, g4mc->tk_mom[lsat]) - std::hypot(g4mc->prm_mass, g4mc->tk_mom[lend])) * 1000.;

        TrSys::Part sm_part(org_part);
        auto sm_rlt = TrSys::Prop::DoPropToZ(g4mc->tk_loc[lend][2], sm_part, TrSys::PropArgs(mscat, eloss), TrSys::InteractionArgs(true));
        
        double sm_lx = sm_part.lx() - mc_part.lx();
        double sm_ly = sm_part.ly() - mc_part.ly();
        double sm_ux = sm_part.ux() - mc_part.ux();
        double sm_uy = sm_part.uy() - mc_part.uy();
        double sm_ea = (std::hypot(g4mc->prm_mass, g4mc->tk_mom[lsat]) - sm_part.eng()) * 1000.;
        
        TrSys::Part bb_part(org_part);
        TrSys::Ionization::kDefaultFormula = TrSys::Ionization::Formula::kBetheBloch;
        auto bb_rlt = TrSys::Prop::DoPropToZ(g4mc->tk_loc[lend][2], bb_part, TrSys::PropArgs(mscat, eloss), TrSys::InteractionArgs(true));
        double bb_ea = (std::hypot(g4mc->prm_mass, g4mc->tk_mom[lsat]) - bb_part.eng()) * 1000.;
   
        TrSys::Part ld_part(org_part);
        TrSys::Ionization::kDefaultFormula = TrSys::Ionization::Formula::kLandau;
        auto ld_rlt = TrSys::Prop::DoPropToZ(g4mc->tk_loc[lend][2], ld_part, TrSys::PropArgs(mscat, eloss), TrSys::InteractionArgs(true));
        double ld_ea = (std::hypot(g4mc->prm_mass, g4mc->tk_mom[lsat]) - ld_part.eng()) * 1000.;
        
        TrSys::Ionization::kDefaultFormula = TrSys::Ionization::Formula::kMixture;
 
        double mom = g4mc->tk_mom[lsat];
        double scale_lu = g4mc->tk_mom[lsat];
        double scale_ee = 1.0 / (1.0 + (g4mc->prm_mass / g4mc->tk_mom[lsat]) * (g4mc->prm_mass / g4mc->tk_mom[lsat]));
        
        Hist::Head("hLxMC")->fillH2D(mom, scale_lu * mc_lx, list->weight);
        Hist::Head("hLyMC")->fillH2D(mom, scale_lu * mc_ly, list->weight);
        Hist::Head("hUxMC")->fillH2D(mom, scale_lu * mc_ux, list->weight);
        Hist::Head("hUyMC")->fillH2D(mom, scale_lu * mc_uy, list->weight);
        Hist::Head("hEaMC")->fillH2D(mom, scale_ee * mc_ea, list->weight);
        
        Hist::Head("hLxSM")->fillH2D(mom, scale_lu * sm_lx, list->weight);
        Hist::Head("hLySM")->fillH2D(mom, scale_lu * sm_ly, list->weight);
        Hist::Head("hUxSM")->fillH2D(mom, scale_lu * sm_ux, list->weight);
        Hist::Head("hUySM")->fillH2D(mom, scale_lu * sm_uy, list->weight);
        Hist::Head("hEaSM")->fillH2D(mom, scale_ee * sm_ea, list->weight);
        
        Hist::Head("hEaBB")->fillH2D(mom, scale_ee * bb_ea, list->weight);
        Hist::Head("hEaLD")->fillH2D(mom, scale_ee * ld_ea, list->weight);
    }
    */

    /*
    TrSys::Tracker data_tracker;
    for (int il = 0; il < 9; ++il) {
        if (trk->lay[il] == 0) continue;
        TrSys::TrackerHit hit(il+1, trk->lay[il]%2==1, trk->lay[il]/2==1, { trk->loc[il][0], trk->loc[il][1], trk->loc[il][2] }, { trk->chrg[il][0], trk->chrg[il][1], trk->chrg[il][2] });
        data_tracker.add_hit(hit);
    }
    
    TrSys::Tof data_tof;
    for (int il = 0; il < 4; ++il) {
        if (!tof->lay[il]) continue;
        TrSys::TofHit hit(il+1, { tof->loc[il][0], tof->loc[il][1], tof->loc[il][2] }, tof->T[il] * TrSys::TofHit::kTransNsToCm, tof->Q[il] );
        data_tof.add_hit(hit);
    }
    
    TrSys::Trd data_trd;
    for (int ih = 0; ih < trd->num_tdHit; ++ih) {
        if (trd->tdHit_len.at(ih) < 0.3) continue;
        double dEdx = trd->tdHit_amp.at(ih) / trd->tdHit_len.at(ih);
        double chrg = std::sqrt(dEdx);
        TrSys::TrdRawHit hit(trd->tdHit_lay.at(ih), { 0.0, 0.0, trd->tdHit_lz.at(ih) }, chrg);
        data_trd.add_rawhit(hit);
    }
    if (data_trd.rawhits().size() < 8) data_trd = TrSys::Trd();
   
    TrSys::Rich data_rich(rich->self_kind);
    for (int ih = 0; ih < rich->self_cld_nhit; ++ih) {
        TrSys::RichHit hit(rich->self_kind, { rich->self_rad_loc[0], rich->self_rad_loc[1], rich->self_rad_loc[2] }, 1.0 / rich->self_cldhit_beta.at(ih) / rich->self_beta_crr);
        data_rich.add_hit(hit);
    }

    TrSys::Part cand_part;
    if (data_tracker.get_num_hit_with_l(TrSys::Tracker::Pattern::kInn) > 0) {
        TrSys::Sys::Stopwatch sw_geom;
        TrSys::GeomTrFit geom_fit(data_tracker.get_hits_with_l(TrSys::Tracker::Pattern::kInn), TrSys::PartList::kProton);
        TrSys::Part part = geom_fit.interpolate_to_z(195.);
        sw_geom.stop();
        double sclr = std::sqrt(g4mc->prm_mom);
        if (geom_fit.status()) {
            double drig = part.irig() - 1.0/g4mc->prm_mom;
            Hist::Head("hNEWr_in_n")->fillH1D(g4mc->prm_mom, list->weight);
            Hist::Head("hNEWr_in_t")->fillH1D(g4mc->prm_mom, sw_geom.time() * list->weight);
            Hist::Head("hNEWr_in_rig")->fillH2D(g4mc->prm_mom, sclr * drig, list->weight);
            Hist::Head("hNEWr_in_nchix")->fillH2D(g4mc->prm_mom, std::log(geom_fit.nchi_x()), list->weight);
            Hist::Head("hNEWr_in_nchiy")->fillH2D(g4mc->prm_mom, std::log(geom_fit.nchi_y()), list->weight);
        }
        //TrSys::Part part2 = geom_fit.interpolate_to_z(g4mc->tk_loc[linn][2]);
        //if (geom_fit.status()) {
        //    double drig2 = part2.irig() - 1.0/g4mc->tk_mom[linn];
        //    Hist::Head("hNEWr2_in_rig")->fillH2D(g4mc->prm_mom, sclr * drig2, list->weight);
        //}

        if (geom_fit.status()) cand_part = geom_fit.part();
    }
    if (trk->ck_status[0]) {
        double sclr = std::sqrt(g4mc->prm_mom);
        double drig = 1.0/trk->ck_rig[0] - 1.0/g4mc->prm_mom;
        Hist::Head("hOFFr_in_rig")->fillH2D(g4mc->prm_mom, sclr * drig, list->weight);
        //double drig2 = 1.0/trk->ck_rig[0] - 1.0/g4mc->tk_mom[linn];
        //Hist::Head("hOFFr2_in_rig")->fillH2D(g4mc->prm_mom, sclr * drig2, list->weight);
    }
    */
    
    return true;
}


#endif // __Analyzer_C__
