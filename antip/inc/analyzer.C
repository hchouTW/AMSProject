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
static TGraphAsymmErrors* gPROFpr27 = nullptr;
static TGraphAsymmErrors* gPROFap2pr = nullptr;

bool Analyzer::build_hist() {
    TFile* fPROPpr = TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/antip/others/ams02_pr.root");
    gPROFpr27 = (TGraphAsymmErrors*) (fPROPpr->Get("gr_exp2"))->Clone("gPROFpr27");
    
    TFile* fPROPap2pr = TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/antip/others/ams02_ap2pr.root");
    gPROFap2pr = (TGraphAsymmErrors*) (fPROPap2pr->Get("gr_exp1"))->Clone("gPROFap2pr");
    
    TFile* fPROPcf = TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/antip/others/expt_ISS.root");
    Hist* hPROFcf = Hist::New("hPROFcf", (TH1*) fPROPcf->Get("hPROF_HC_expt_gISS"));

    file->cd();

    std::vector<double> vtme( {
        1305417600, 1312416000, 1319414400, 1326412800, 1333411200,
        1340409600, 1347408000, 1354406400, 1361404800, 1368403200,
        1375401600, 1382400000, 1389398400, 1396396800, 1403395200,
        1410393600, 1417392000, 1424390400, 1431388800, 1438387200,
        1445385600, 1452384000, 1459382400, 1466380800, 1473379200, 
        1480377600, 1487376000, 1494374400, 1501372800, 1508371200,
        1515369600, 1522368000, 1529366400, 1536364800, 1543363200,
        1550361600 } );

    Axis AXtme("Time", vtme);

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
    Hist::New("hPROF_CK_mass_ISS", AXPROFmass);
    Hist::New("hPROF_CK_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF_CK_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_ISS_CF", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC_FLUX10", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC_FLUX27", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC_FLUX10_CF", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC_FLUX27_CF", HistAxis(AXrig3));
    
    Hist::New("hPROF_KF_mass_ISS", AXPROFmass);
    Hist::New("hPROF_KF_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF_KF_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_ISS_CF", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC_FLUX10", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC_FLUX27", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC_FLUX10_CF", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC_FLUX27_CF", HistAxis(AXrig3));
   
    Hist::New("hPROF_HC_expt_gISS", HistAxis(AXrig3));
    Hist::New("hPROF_HC_expt_evt_gISS", HistAxis(AXrig3));
    Hist::New("hPROF_HC_expt_evt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF_HC_expt_evt_ISS_CF", HistAxis(AXrig3));
    
    Hist::New("hPROF_HC_phtrg_ISS", AXPROFmass);
    Hist::New("hPROF_HC_altrg_ISS", AXPROFmass);
    Hist::New("hPROF_HC_phtrg_ISS_CF", AXPROFmass);
    Hist::New("hPROF_HC_altrg_ISS_CF", AXPROFmass);
    
    Hist::New("hPROF_HC_rig_phtrg_ISS", AXrig3);
    Hist::New("hPROF_HC_rig_altrg_ISS", AXrig3);
    Hist::New("hPROF_HC_rig_phtrg_ISS_CF", AXrig3);
    Hist::New("hPROF_HC_rig_altrg_ISS_CF", AXrig3);

    Hist::New("hPROF_HC_mass_ISS", AXPROFmass);
    Hist::New("hPROF_HC_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF_HC_mass_MC_PURE_FLUX10", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_PURE_FLUX27", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_PURE_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_PURE_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF2_HC_mass_ISS", AXPROFmass);
    Hist::New("hPROF2_HC_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF3_HC_mass_ISS", AXPROFmass);
    Hist::New("hPROF3_HC_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF3_HC_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF3_HC_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF3_HC_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF3_HC_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF_HC_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_ISS_CF", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC_FLUX10", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC_FLUX27", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC_FLUX10_CF", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC_FLUX27_CF", HistAxis(AXrig3));
    
    Hist::New("hPROF_HC_mass_PI_POS_MC_FLUX10" , AXPROFmass);
    Hist::New("hPROF_HC_mass_PI_POS_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF_HC_mass_PI_NEG_MC_FLUX10" , AXPROFmass);
    Hist::New("hPROF_HC_mass_PI_NEG_MC_FLUX27", AXPROFmass);
                                       
    Hist::New("hPROF_HC_mass_K_POS_MC_FLUX10" , AXPROFmass);
    Hist::New("hPROF_HC_mass_K_POS_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF_HC_mass_K_NEG_MC_FLUX10" , AXPROFmass);
    Hist::New("hPROF_HC_mass_K_NEG_MC_FLUX27", AXPROFmass);
    
    Hist::New("hPROF_HC_mass_lchix_ISS", HistAxis(AXPROFmass, AXnchi));
    Hist::New("hPROF_HC_mass_lchiy_ISS", HistAxis(AXPROFmass, AXnchi));
    Hist::New("hPROF_HC_mass_lchib_ISS", HistAxis(AXPROFmass, AXnchi));
    
    Hist::New("hPROF_HC_mass_nseg_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_mass_nhit_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_mass_nvtxx_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_mass_nvtxy_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hPROF_HC_mass_tkL1_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_mass_tkL2_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_mass_tkL9_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hPROF_HC_mass_inn_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_mass_out_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hPROF_HC_mass_ntd1_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_mass_ntd2_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_mass_ntd3_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hPROF_CK_RH_mass_ISS", AXPROFmass);
    Hist::New("hPROF_CK_RH_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF_CK_RH_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF_CK_RH_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF_CK_RH_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF_CK_RH_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF_KF_RH_mass_ISS", AXPROFmass);
    Hist::New("hPROF_KF_RH_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF_KF_RH_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF_KF_RH_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF_KF_RH_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF_KF_RH_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF_HC_RH_mass_ISS", AXPROFmass);
    Hist::New("hPROF_HC_RH_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF_HC_RH_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF_HC_RH_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF_HC_RH_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF_HC_RH_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF_HC_RH_mass_MC_PURE_FLUX10", AXPROFmass);
    Hist::New("hPROF_HC_RH_mass_MC_PURE_FLUX27", AXPROFmass);
    Hist::New("hPROF_HC_RH_mass_MC_PURE_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF_HC_RH_mass_MC_PURE_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF2_HC_RH_mass_ISS", AXPROFmass);
    Hist::New("hPROF2_HC_RH_mass_ISS_CF", AXPROFmass);
    Hist::New("hPROF2_HC_RH_mass_MC_FLUX10", AXPROFmass);
    Hist::New("hPROF2_HC_RH_mass_MC_FLUX27", AXPROFmass);
    Hist::New("hPROF2_HC_RH_mass_MC_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF2_HC_RH_mass_MC_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF2_HC_RH_mass_MC_PURE_FLUX10", AXPROFmass);
    Hist::New("hPROF2_HC_RH_mass_MC_PURE_FLUX27", AXPROFmass);
    Hist::New("hPROF2_HC_RH_mass_MC_PURE_FLUX10_CF", AXPROFmass);
    Hist::New("hPROF2_HC_RH_mass_MC_PURE_FLUX27_CF", AXPROFmass);
    
    Hist::New("hPROF_HC_RH_mass_lchix_ISS", HistAxis(AXPROFmass, AXnchi));
    Hist::New("hPROF_HC_RH_mass_lchiy_ISS", HistAxis(AXPROFmass, AXnchi));
    Hist::New("hPROF_HC_RH_mass_lchib_ISS", HistAxis(AXPROFmass, AXnchi));
    
    Hist::New("hPROF_HC_RH_mass_inn_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_RH_mass_out_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hPROF_HC_RH_mass_npmt_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_RH_mass_nhit_ISS", HistAxis(AXPROFmass, AXnext));
    
    Hist::New("hPROF_HC_RH_mass_nvtxx_ISS", HistAxis(AXPROFmass, AXnext));
    Hist::New("hPROF_HC_RH_mass_nvtxy_ISS", HistAxis(AXPROFmass, AXnext));
    
    Axis AXPROFxxx("", 1000, 0, 20);
    Hist::New("hPROF_HC_RH_mass_tofq_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_tof_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_tof2_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_rh_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_rh2_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    
    Hist::New("hPROF_HC_RH_mass_npe_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_chi_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    
    Hist::New("hPROF_HC_RH_mass_mij_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_exp_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_chg_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_bdr_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_trc_ISS", HistAxis(AXPROFmass, AXPROFxxx));
    Hist::New("hPROF_HC_RH_mass_acc_ISS", HistAxis(AXPROFmass, AXPROFxxx));

    Axis AXPROFrvar("RVAR", 1000, -2, 2);
    Hist::New("hPROF_CK_rvar", HistAxis(AXPROFrvar));
    Hist::New("hPROF_CK_rvar_MC_FLUX10", HistAxis(AXPROFrvar));
    Hist::New("hPROF_CK_rvar_MC_FLUX27", HistAxis(AXPROFrvar));
    Hist::New("hPROF_CK_rvar_cut", HistAxis(AXPROFrvar));
    Hist::New("hPROF_CK_rvar_cut_MC_FLUX10", HistAxis(AXPROFrvar));
    Hist::New("hPROF_CK_rvar_cut_MC_FLUX27", HistAxis(AXPROFrvar));
    
    Hist::New("hPROF_KF_rvar", HistAxis(AXPROFrvar));
    Hist::New("hPROF_KF_rvar_MC_FLUX10", HistAxis(AXPROFrvar));
    Hist::New("hPROF_KF_rvar_MC_FLUX27", HistAxis(AXPROFrvar));
    Hist::New("hPROF_KF_rvar_cut", HistAxis(AXPROFrvar));
    Hist::New("hPROF_KF_rvar_cut_MC_FLUX10", HistAxis(AXPROFrvar));
    Hist::New("hPROF_KF_rvar_cut_MC_FLUX27", HistAxis(AXPROFrvar));
    
    Hist::New("hPROF_HC_rvar", HistAxis(AXPROFrvar));
    Hist::New("hPROF_HC_rvar_MC_FLUX10", HistAxis(AXPROFrvar));
    Hist::New("hPROF_HC_rvar_MC_FLUX27", HistAxis(AXPROFrvar));
    Hist::New("hPROF_HC_rvar_cut", HistAxis(AXPROFrvar));
    Hist::New("hPROF_HC_rvar_cut_MC_FLUX10", HistAxis(AXPROFrvar));
    Hist::New("hPROF_HC_rvar_cut_MC_FLUX27", HistAxis(AXPROFrvar));
    
    Hist::New("h_MC_cnt", HistAxis(AXrig));
    Hist::New("hT_MC_cnt", HistAxis(AXrig2));
    
    Axis AXLllr("TRD estimator", 200, 0.0, 1.6);
    Hist::New("hLP_llr", HistAxis(AXrig, AXLllr));
    Hist::New("hLN_llr", HistAxis(AXrig, AXLllr));

    Axis AXLsqrm("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 75, -2.5, 3.5);
    Axis AXLsqrm2("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 50, -2.5, 3.5);
    Hist::New("hLP_sqrm_pr",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_sqrm_el",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_sqrm_pi",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLP_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_sqrm",      HistAxis(AXrig, AXLsqrm));
    
    Hist::New("hTLP_sqrm_pr",   HistAxis(AXtme, AXrig2, AXLsqrm2));
    Hist::New("hTLN_sqrm_el",   HistAxis(AXtme, AXrig2, AXLsqrm2));
    Hist::New("hTLN_sqrm_pi",   HistAxis(AXtme, AXrig2, AXLsqrm2));
    Hist::New("hTLP_sqrm",      HistAxis(AXtme, AXrig2, AXLsqrm2));
    Hist::New("hTLN_sqrm",      HistAxis(AXtme, AXrig2, AXLsqrm2));
    
    Hist::New("hLP_cf80_sqrm_pr",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf80_sqrm_el",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf80_sqrm_pi",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLP_cf80_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf80_sqrm",      HistAxis(AXrig, AXLsqrm));
    
    Hist::New("hLP_cf85_sqrm_pr",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf85_sqrm_el",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf85_sqrm_pi",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLP_cf85_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf85_sqrm",      HistAxis(AXrig, AXLsqrm));
    
    Hist::New("hLP_cf90_sqrm_pr",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf90_sqrm_el",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf90_sqrm_pi",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLP_cf90_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf90_sqrm",      HistAxis(AXrig, AXLsqrm));
    
    Hist::New("hLP_cnt", HistAxis(AXrig));
    Hist::New("hLN_cnt", HistAxis(AXrig));
    
    Hist::New("hL_MC_cnt", HistAxis(AXrig));
    Hist::New("hTL_MC_cnt", HistAxis(AXrig2));
    
    Axis AXMllr("TRD estimator", 100, 0.0, 1.6);
    Hist::New("hMP_llr", HistAxis(AXrig, AXMllr));
    Hist::New("hMN_llr", HistAxis(AXrig, AXMllr));
    
    Axis AXMsqrm("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 100, -1.75, 3.25);
    Axis AXMsqrm2("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 75, -1.75, 3.25);
    Hist::New("hMP_sqrm_pr", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_sqrm_el", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMP_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_sqrm",    HistAxis(AXrig, AXMsqrm));
    
    Hist::New("hTMP_sqrm_pr", HistAxis(AXtme, AXrig2, AXMsqrm2));
    Hist::New("hTMN_sqrm_el", HistAxis(AXtme, AXrig2, AXMsqrm2));
    Hist::New("hTMP_sqrm",    HistAxis(AXtme, AXrig2, AXMsqrm2));
    Hist::New("hTMN_sqrm",    HistAxis(AXtme, AXrig2, AXMsqrm2));
    
    Hist::New("hMP_cf80_sqrm_pr", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf80_sqrm_el", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMP_cf80_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf80_sqrm",    HistAxis(AXrig, AXMsqrm));
    
    Hist::New("hMP_cf85_sqrm_pr", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf85_sqrm_el", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMP_cf85_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf85_sqrm",    HistAxis(AXrig, AXMsqrm));
    
    Hist::New("hMP_cf90_sqrm_pr", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf90_sqrm_el", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMP_cf90_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf90_sqrm",    HistAxis(AXrig, AXMsqrm));
    
    Hist::New("hMP_cnt", HistAxis(AXrig));
    Hist::New("hMN_cnt", HistAxis(AXrig));
    
    Hist::New("hM_MC_cnt", HistAxis(AXrig));
    Hist::New("hTM_MC_cnt", HistAxis(AXrig2));
    
    Axis AXIlchi("Log(#chi^{2}/DOF)", 100, -4, 8);
    Hist::New("hIP_lchiy", HistAxis(AXrig, AXIlchi));
    Hist::New("hIN_lchiy", HistAxis(AXrig, AXIlchi));
    
    Axis AXIllr("TRD estimator",  150, 0.0, 1.6);
    Axis AXIllr2("TRD estimator", 100, 0.0, 1.6);
    Hist::New("hIP_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIP_llr_pr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_llr_el", HistAxis(AXrig, AXIllr));
    
    Hist::New("hTIP_llr", HistAxis(AXtme, AXrig2, AXIllr2));
    Hist::New("hTIN_llr", HistAxis(AXtme, AXrig2, AXIllr2));
    Hist::New("hTIP_llr_pr", HistAxis(AXtme, AXrig2, AXIllr2));
    Hist::New("hTIN_llr_el", HistAxis(AXtme, AXrig2, AXIllr2));
    
    Hist::New("hIP_cf80_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf80_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIP_cf80_llr_pr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf80_llr_el", HistAxis(AXrig, AXIllr));
    
    Hist::New("hIP_cf85_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf85_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIP_cf85_llr_pr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf85_llr_el", HistAxis(AXrig, AXIllr));
    
    Hist::New("hIP_cf90_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf90_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIP_cf90_llr_pr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf90_llr_el", HistAxis(AXrig, AXIllr));
    
    Hist::New("hIP_cnt", HistAxis(AXrig));
    Hist::New("hIN_cnt", HistAxis(AXrig));
    
    Hist::New("hI_MC_cnt", HistAxis(AXrig));
    Hist::New("hTI_MC_cnt", HistAxis(AXrig2));
    
    Axis AXHllr("TRD estimator", 100, 0.0, 1.6);
    Axis AXHlchi2("Log(#chi^{2}/DOF)", 200, -4, 8);
    Axis AXHdchi2("#Delta Log(#chi^{2}/DOF)", 200, -8, 8);
    Axis AXHlrvar("|RVAR|", 200, 0, 5);
    
    Axis AXHlchi2_cut("Log(#chi^{2}/DOF)", 100, -3, 3);
    Axis AXHlrvar_cut("|RVAR|", 100, 0, 3);

    Axis AXHnumtk("Ext-SuperL", 4, 0.0, 4.0);

    Hist::New("hHPl1_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNl1_llr", HistAxis(AXrig, AXHllr));

    Hist::New("hHPl1_lchix_in", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl1_lchix_in", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl1_lchix_l1", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl1_lchix_l1", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl1_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl1_lchix", HistAxis(AXrig, AXHlchi2));
    
    Hist::New("hHPl1_dchix_l1", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNl1_dchix_l1", HistAxis(AXrig, AXHdchi2));

    Hist::New("hHPl1_lchiy_in", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl1_lchiy_in", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl1_lchiy_l1", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl1_lchiy_l1", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl1_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl1_lchiy", HistAxis(AXrig, AXHlchi2));
    
    Hist::New("hHPl1_dchiy_l1", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNl1_dchiy_l1", HistAxis(AXrig, AXHdchi2));
    
    Hist::New("hHPl1_lrvar", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNl1_lrvar", HistAxis(AXrig, AXHlrvar));

    Hist::New("hHPl1_num_tk", HistAxis(AXrig, AXHnumtk));
    Hist::New("hHNl1_num_tk", HistAxis(AXrig, AXHnumtk));
    
    Hist::New("hHPl1_lchiy_lrvar", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNl1_lchiy_lrvar", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));

    Hist::New("hHPl1_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNl1_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHPl1_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNl1_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));

    Hist::New("hHPl1_lchiy_cut", HistAxis(AXrig, AXHlchi2_cut));
    Hist::New("hHNl1_lchiy_cut", HistAxis(AXrig, AXHlchi2_cut));
    Hist::New("hHPl1_lrvar_cut", HistAxis(AXrig, AXHlrvar_cut));
    Hist::New("hHNl1_lrvar_cut", HistAxis(AXrig, AXHlrvar_cut));

    Hist::New("hHPl1_cnt", HistAxis(AXrig));
    Hist::New("hHNl1_cnt", HistAxis(AXrig));
    
    Hist::New("hHl1_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hHl1_cnt_MC_FLUX27", HistAxis(AXrig));
    
    
    Hist::New("hHPl9_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNl9_llr", HistAxis(AXrig, AXHllr));

    Hist::New("hHPl9_lchix_in", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl9_lchix_in", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl9_lchix_l9", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl9_lchix_l9", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl9_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl9_lchix", HistAxis(AXrig, AXHlchi2));
    
    Hist::New("hHPl9_dchix_l9", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNl9_dchix_l9", HistAxis(AXrig, AXHdchi2));

    Hist::New("hHPl9_lchiy_in", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl9_lchiy_in", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl9_lchiy_l9", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl9_lchiy_l9", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl9_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl9_lchiy", HistAxis(AXrig, AXHlchi2));
    
    Hist::New("hHPl9_dchiy_l9", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNl9_dchiy_l9", HistAxis(AXrig, AXHdchi2));

    Hist::New("hHPl9_lrvar", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNl9_lrvar", HistAxis(AXrig, AXHlrvar));
    
    Hist::New("hHPl9_num_tk", HistAxis(AXrig, AXHnumtk));
    Hist::New("hHNl9_num_tk", HistAxis(AXrig, AXHnumtk));
    
    Hist::New("hHPl9_lchiy_lrvar", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNl9_lchiy_lrvar", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    
    Hist::New("hHPl9_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNl9_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHPl9_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNl9_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    
    Hist::New("hHPl9_lchiy_cut", HistAxis(AXrig, AXHlchi2_cut));
    Hist::New("hHNl9_lchiy_cut", HistAxis(AXrig, AXHlchi2_cut));
    Hist::New("hHPl9_lrvar_cut", HistAxis(AXrig, AXHlrvar_cut));
    Hist::New("hHNl9_lrvar_cut", HistAxis(AXrig, AXHlrvar_cut));

    Hist::New("hHPl9_cnt", HistAxis(AXrig));
    Hist::New("hHNl9_cnt", HistAxis(AXrig));
    
    Hist::New("hHl9_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hHl9_cnt_MC_FLUX27", HistAxis(AXrig));
    
    
    Hist::New("hHPfs_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNfs_llr", HistAxis(AXrig, AXHllr));

    Hist::New("hHPfs_lchix_l1", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchix_l1", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPfs_lchix_l9", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchix_l9", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPfs_lchix_fs", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchix_fs", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPfs_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchix", HistAxis(AXrig, AXHlchi2));
    
    Hist::New("hHPfs_dchix_l1", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNfs_dchix_l1", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHPfs_dchix_l9", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNfs_dchix_l9", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHPfs_dchix_fs", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNfs_dchix_fs", HistAxis(AXrig, AXHdchi2));
    
    Hist::New("hHPfs_lchiy_l1", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchiy_l1", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPfs_lchiy_l9", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchiy_l9", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPfs_lchiy_fs", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchiy_fs", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPfs_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchiy", HistAxis(AXrig, AXHlchi2));
    
    Hist::New("hHPfs_dchiy_l1", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNfs_dchiy_l1", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHPfs_dchiy_l9", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNfs_dchiy_l9", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHPfs_dchiy_fs", HistAxis(AXrig, AXHdchi2));
    Hist::New("hHNfs_dchiy_fs", HistAxis(AXrig, AXHdchi2));
    
    Hist::New("hHPfs_lrvar_l1", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNfs_lrvar_l1", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHPfs_lrvar_l9", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNfs_lrvar_l9", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHPfs_lrvar_fs", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNfs_lrvar_fs", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHPfs_lrvar", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNfs_lrvar", HistAxis(AXrig, AXHlrvar));
    
    Hist::New("hHPfs_num_tk", HistAxis(AXrig, AXHnumtk));
    Hist::New("hHNfs_num_tk", HistAxis(AXrig, AXHnumtk));
    
    Hist::New("hHPfs_lchiy_lrvar", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNfs_lchiy_lrvar", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    
    Hist::New("hHPfs_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNfs_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHPfs_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    Hist::New("hHNfs_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2_cut, AXHlrvar_cut));
    
    Hist::New("hHPfs_lchiy_cut", HistAxis(AXrig, AXHlchi2_cut));
    Hist::New("hHNfs_lchiy_cut", HistAxis(AXrig, AXHlchi2_cut));
    Hist::New("hHPfs_lrvar_cut", HistAxis(AXrig, AXHlrvar_cut));
    Hist::New("hHNfs_lrvar_cut", HistAxis(AXrig, AXHlrvar_cut));

    Hist::New("hHPfs_cnt", HistAxis(AXrig));
    Hist::New("hHNfs_cnt", HistAxis(AXrig));
    
    Hist::New("hHfs_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hHfs_cnt_MC_FLUX27", HistAxis(AXrig));

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
    //if (tof->nchi_t > 10.0) return false;
    //if (tof->nchi_c > 10.0) return false;

    // TRK
    //int tkL1next = trk->ext_num_hit[0];
    //int tkL9next = trk->ext_num_hit[8];
    //int tkL2next = trk->ext_num_hit[1];
    //int tkInnext = trk->ext_num_hit[2] + trk->ext_num_hit[3] + trk->ext_num_hit[4] + trk->ext_num_hit[5] + trk->ext_num_hit[6] + trk->ext_num_hit[7];
    //if (tkInnext > 6) return false;
    //if (tkL1next > 6) return false;
    //if (tkL9next > 2) return false;
    //if (tkL2next > 4) return false;

    //bool is_all_inn = 
    //    (trk->lay[2] == 3 &&
    //     trk->lay[3] == 3 &&
    //     trk->lay[4] == 3 &&
    //     trk->lay[5] == 3 &&
    //     trk->lay[6] == 3 &&
    //     trk->lay[7] == 3);
    //if (!is_all_inn) return false;

    // TRD
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    
    // ECAL
    //if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    if (g4mc != nullptr) Hist::Head("h_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
    if (g4mc != nullptr) Hist::Head("hT_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);

    return true;
}

bool Analyzer::process_data() {
    process_data_prof();
    process_data_l();
    process_data_m();
    process_data_i();
    process_data_h();
    return true;
}

bool Analyzer::process_data_prof() {
    double mc_flux10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom, -1.7) : 1.0;
    double mc_flux27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        double crr_flux_pr = gPROFpr27->Eval(g4mc->prm_mom, 0, "");
        if (crr_flux_pr < 0) crr_flux_pr = 0.0;

        double crr_flux_ap = crr_flux_pr * gPROFap2pr->Eval(g4mc->prm_mom, 0, "");
        if (crr_flux_ap < 0) crr_flux_ap = 0.0;

        mc_flux10 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap);
        mc_flux27 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap);
    }
    
    bool is_clean_tk = 
        !((trk->ext_num_hit[2] > 0 && trk->ext_num_hit[3] > 0) || 
          (trk->ext_num_hit[4] > 0 && trk->ext_num_hit[5] > 0) || 
          (trk->ext_num_hit[6] > 0 && trk->ext_num_hit[7] > 0));

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

    bool is_trd_ok = (trd->num_vtx[0] <= 20 && trd->num_vtx[1] <= 20);
    
    // Trigger
    bool phtrg = (trg->bit&8) == 8;
    bool untrg = (trg->bit&2) == 2;
    int  wgtrg = (phtrg ? 1 : 100);
    
    if (trk->lay[1] != 0 &&
        trk->ck_status[0] && 
        std::log(trk->ck_nchi[0][0]) < 2.0 &&
        std::log(trk->ck_nchi[0][1]) < 2.0 &&
        tof->status &&
        std::log(tof->nchi_t) < 2.0 &&
        std::log(tof->nchi_c) < 2.0 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        is_rich_ok &&
        is_clean_tk &&
        is_trd_ok &&
        phtrg
        ) {
        
        double sign = (trk->ck_rig[0] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * std::sqrt((trk->ck_rig[0] * (1.0 / tof->beta + 1.0)) * (trk->ck_rig[0] * (1.0 / tof->beta - 1.0)));
        double cfr  = CheckType(Type::ISS) ? std::abs(trk->ck_rig[0]/rti->max_IGRF) : 0.0;

        double wgtcf = (*Hist::Head("hPROFcf"))()->GetBinContent( (*Hist::Head("hPROFcf"))()->FindBin(std::abs(trk->ck_rig[0])/0.75) );

        if (tof->beta > 0.5 && tof->beta < 0.8) {
            Hist::Head("hPROF_CK_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hPROF_CK_mass_ISS_CF")->fillH1D(mass, list->weight);
            
            Hist::Head("hPROF_CK_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hPROF_CK_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hPROF_CK_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hPROF_CK_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
        }
        if (trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_ISS")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
        if (trk->ck_rig[0] > 0 && cfr > 0.75) Hist::Head("hPROF_CK_mass_cnt_ISS_CF")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);

        if (trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
        if (trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC_FLUX10")->fillH1D(std::abs(trk->ck_rig[0]), mc_flux10 * list->weight);
        if (trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC_FLUX27")->fillH1D(std::abs(trk->ck_rig[0]), mc_flux27 * list->weight);
        if (trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC_FLUX10_CF")->fillH1D(std::abs(trk->ck_rig[0]), wgtcf * mc_flux10 * list->weight);
        if (trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC_FLUX27_CF")->fillH1D(std::abs(trk->ck_rig[0]), wgtcf * mc_flux27 * list->weight);
    }
    
    if (trk->lay[1] != 0 &&
        trk->kf_status[0] && 
        std::log(trk->kf_nchi[0][0]) < 2.0 &&
        std::log(trk->kf_nchi[0][1]) < 2.0 &&
        tof->status &&
        std::log(tof->nchi_t) < 2.0 &&
        std::log(tof->nchi_c) < 2.0 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        is_rich_ok &&
        is_clean_tk &&
        is_trd_ok &&
        phtrg
        ) {
        
        double sign = (trk->kf_cen_rig[0] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * std::sqrt((trk->kf_cen_rig[0] * (1.0 / tof->beta + 1.0)) * (trk->kf_cen_rig[0] * (1.0 / tof->beta - 1.0)));
        double cfr  = CheckType(Type::ISS) ? std::abs(trk->kf_top_rig[0]/rti->max_IGRF) : 0.0;
        
        double wgtcf = (*Hist::Head("hPROFcf"))()->GetBinContent( (*Hist::Head("hPROFcf"))()->FindBin(std::abs(trk->kf_top_rig[0])/0.75) );

        if (tof->beta > 0.5 && tof->beta < 0.8) {
            Hist::Head("hPROF_KF_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hPROF_KF_mass_ISS_CF")->fillH1D(mass, list->weight);

            Hist::Head("hPROF_KF_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hPROF_KF_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hPROF_KF_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hPROF_KF_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
        }
        if (trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_ISS")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
        if (trk->kf_top_rig[0] > 0 && cfr > 0.75) Hist::Head("hPROF_KF_mass_cnt_ISS_CF")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
        
        if (trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
        if (trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC_FLUX10")->fillH1D(std::abs(trk->kf_top_rig[0]), mc_flux10 * list->weight);
        if (trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC_FLUX27")->fillH1D(std::abs(trk->kf_top_rig[0]), mc_flux27 * list->weight);
        if (trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC_FLUX10_CF")->fillH1D(std::abs(trk->kf_top_rig[0]), wgtcf * mc_flux10 * list->weight);
        if (trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC_FLUX27_CF")->fillH1D(std::abs(trk->kf_top_rig[0]), wgtcf * mc_flux27 * list->weight);
    }
    
    if (trk->lay[1] != 0 &&
        hyc->geom_status[0] &&
        std::log(hyc->geom_nchi_x[0]) < 1.75 &&
        std::log(hyc->geom_nchi_y[0]) < 1.75 &&
        hyc->vel_status[1] &&
        std::log(hyc->vel_nchi[1]) < 2.00 &&
        hyc->mutr_status[1] && 
        //std::log(hyc->mutr_nchi_x[1]) < 1.25 &&
        //std::log(hyc->mutr_nchi_y[1]) < 1.25 &&
        //std::log(hyc->mutr_nchi_b[1]) < 1.50 &&
        std::log(hyc->mutr_nchi_x[1]) < 1.50 &&
        std::log(hyc->mutr_nchi_y[1]) < 1.50 &&
        std::log(hyc->mutr_nchi_b[1]) < 1.50 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        is_rich_ok &&
        is_clean_tk &&
        is_trd_ok
        ) {

        double sign = (hyc->mutr_top_rig[1] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * hyc->mutr_mass[1];
        double cfr  = CheckType(Type::ISS) ? std::abs(hyc->mutr_top_rig[1]/rti->max_IGRF) : 0.0;
        
        double wgtcf = (*Hist::Head("hPROFcf"))()->GetBinContent( (*Hist::Head("hPROFcf"))()->FindBin(std::abs(hyc->mutr_top_rig[1])/0.75) );
        
        bool pure_mc = CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3]);
            
        bool is_new_utime = false;
        if (CheckType(Type::ISS) && utime_cur != list->utime) { utime_cur = list->utime; is_new_utime = true; }
        TH1D* hexpt = (TH1D*) (*Hist::Head("hPROF_HC_expt_gISS"))();
        int   hexpt_rbin = CheckType(Type::ISS) ? hexpt->FindBin(rti->max_IGRF) : 0;

        if (is_new_utime) {
            for (int bin = hexpt_rbin + 1; bin <= hexpt->GetXaxis()->GetNbins(); bin++) {
                hexpt->Fill( hexpt->GetBinCenter(bin), (CheckType(Type::ISS) ? rti->live_time : 1.0));
            }
        }
        if (CheckType(Type::ISS)) {
            for (int bin = hexpt_rbin + 1; bin <= hexpt->GetXaxis()->GetNbins(); bin++) {
                Hist::Head("hPROF_HC_expt_evt_gISS")->fillH1D(hexpt->GetBinCenter(bin), rti->live_time * list->weight);
            }
        }
        
        if (hyc->mutr_top_bta[1] > 0.5 && hyc->mutr_top_bta[1] < 0.8) {
            if (phtrg) Hist::Head("hPROF_HC_phtrg_ISS")->fillH1D(mass, list->weight);
            Hist::Head("hPROF_HC_altrg_ISS")->fillH1D(mass, wgtrg * list->weight);
            if (cfr > 0.75 && phtrg) Hist::Head("hPROF_HC_phtrg_ISS_CF")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hPROF_HC_altrg_ISS_CF")->fillH1D(mass, wgtrg * list->weight);
            
            if (phtrg) Hist::Head("hPROF_HC_rig_phtrg_ISS")->fillH1D(std::abs(hyc->mutr_top_rig[1]), list->weight);
            Hist::Head("hPROF_HC_rig_altrg_ISS")->fillH1D(std::abs(hyc->mutr_top_rig[1]), wgtrg * list->weight);
            if (cfr > 0.75 && phtrg) Hist::Head("hPROF_HC_rig_phtrg_ISS_CF")->fillH1D(std::abs(hyc->mutr_top_rig[1]), list->weight);
            if (cfr > 0.75) Hist::Head("hPROF_HC_rig_altrg_ISS_CF")->fillH1D(std::abs(hyc->mutr_top_rig[1]), wgtrg * list->weight);
        }

        if (phtrg && hyc->mutr_top_bta[1] > 0.5 && hyc->mutr_top_bta[1] < 0.8) {
            if (CheckType(Type::ISS)) {
                for (int bin = hexpt_rbin + 1; bin <= hexpt->GetXaxis()->GetNbins(); bin++) {
                    Hist::Head("hPROF_HC_expt_evt_ISS")->fillH1D(hexpt->GetBinCenter(bin), rti->live_time * list->weight);
                    if (cfr > 0.75) Hist::Head("hPROF_HC_expt_evt_ISS_CF")->fillH1D(hexpt->GetBinCenter(bin), rti->live_time * list->weight);
                }
            }
            
            Hist::Head("hPROF_HC_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hPROF_HC_mass_ISS_CF")->fillH1D(mass, list->weight);
            
            Hist::Head("hPROF_HC_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hPROF_HC_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hPROF_HC_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hPROF_HC_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
            
            if (pure_mc) Hist::Head("hPROF_HC_mass_MC_PURE_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            if (pure_mc) Hist::Head("hPROF_HC_mass_MC_PURE_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            if (pure_mc) Hist::Head("hPROF_HC_mass_MC_PURE_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            if (pure_mc) Hist::Head("hPROF_HC_mass_MC_PURE_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
                
            //if (mass < -1.6) clone_tree->Fill();

            Hist::Head("hPROF_HC_mass_lchix_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_x[1]), list->weight);
            Hist::Head("hPROF_HC_mass_lchiy_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_y[1]), list->weight);
            Hist::Head("hPROF_HC_mass_lchib_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_b[1]), list->weight);
            
            Hist::Head("hPROF_HC_mass_nseg_ISS")->fillH2D(mass, trd->num_extra_seg, list->weight);
            Hist::Head("hPROF_HC_mass_nhit_ISS")->fillH2D(mass, trd->num_extra_hit, list->weight);
            Hist::Head("hPROF_HC_mass_nvtxx_ISS")->fillH2D(mass, trd->num_vtx[0], list->weight);
            Hist::Head("hPROF_HC_mass_nvtxy_ISS")->fillH2D(mass, trd->num_vtx[1], list->weight);
            
            Hist::Head("hPROF_HC_mass_tkL1_ISS")->fillH2D(mass, trk->ext_num_hit[0], list->weight);
            Hist::Head("hPROF_HC_mass_tkL2_ISS")->fillH2D(mass, trk->ext_num_hit[1], list->weight);
            Hist::Head("hPROF_HC_mass_tkL9_ISS")->fillH2D(mass, trk->ext_num_hit[8], list->weight);
            
            Hist::Head("hPROF_HC_mass_inn_ISS")->fillH2D(mass, rich->self_nhit_other_inn, list->weight);
            Hist::Head("hPROF_HC_mass_out_ISS")->fillH2D(mass, rich->self_nhit_other_out, list->weight);
            
            Hist::Head("hPROF_HC_mass_ntd1_ISS")->fillH2D(mass, trd->tkLLR_num_hit, list->weight);
            Hist::Head("hPROF_HC_mass_ntd2_ISS")->fillH2D(mass, trd->num_tkHit, list->weight);
            Hist::Head("hPROF_HC_mass_ntd3_ISS")->fillH2D(mass, trd->tkLLR_num_hit-trd->num_tkHit, list->weight);
                
            // pion 0 ~ 0.3
            if (std::abs(mass) < 0.28) {
                Hist::Head("hPROF_HC_mass_PI_POS_MC_FLUX10")->fillH1D( std::abs(mass), mc_flux10 * list->weight);
                Hist::Head("hPROF_HC_mass_PI_POS_MC_FLUX27")->fillH1D( std::abs(mass), mc_flux27 * list->weight);
                Hist::Head("hPROF_HC_mass_PI_NEG_MC_FLUX10")->fillH1D(-std::abs(mass), mc_flux10 * list->weight);
                Hist::Head("hPROF_HC_mass_PI_NEG_MC_FLUX27")->fillH1D(-std::abs(mass), mc_flux27 * list->weight);
            }
            
            // koan -0.35 ~ -0.75
            if (mass > -0.75 && mass < -0.35) {
                Hist::Head("hPROF_HC_mass_K_POS_MC_FLUX10")->fillH1D( std::abs(mass), mc_flux10 * list->weight);
                Hist::Head("hPROF_HC_mass_K_POS_MC_FLUX27")->fillH1D( std::abs(mass), mc_flux27 * list->weight);
                Hist::Head("hPROF_HC_mass_K_NEG_MC_FLUX10")->fillH1D(-std::abs(mass), mc_flux10 * list->weight);
                Hist::Head("hPROF_HC_mass_K_NEG_MC_FLUX27")->fillH1D(-std::abs(mass), mc_flux27 * list->weight);
            }
        }
        
        if (phtrg && hyc->mutr_top_bta[1] > 0.5 && hyc->mutr_top_bta[1] < 0.8 && trd->num_vtx[0] <= 12 && trd->num_vtx[1] <= 12) {
            Hist::Head("hPROF2_HC_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hPROF2_HC_mass_ISS_CF")->fillH1D(mass, list->weight);
            
            Hist::Head("hPROF2_HC_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hPROF2_HC_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hPROF2_HC_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hPROF2_HC_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
        }
        if (phtrg && hyc->mutr_top_bta[1] > 0.5 && hyc->mutr_top_bta[1] < 0.8 && trd->num_vtx[0] <= 4 && trd->num_vtx[1] <= 8) {
            Hist::Head("hPROF3_HC_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hPROF3_HC_mass_ISS_CF")->fillH1D(mass, list->weight);
            
            Hist::Head("hPROF3_HC_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hPROF3_HC_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hPROF3_HC_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hPROF3_HC_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
        }

        if (phtrg) {
            if (hyc->mutr_top_rig[1] > 0) Hist::Head("hPROF_HC_mass_cnt_ISS")->fillH1D(std::abs(hyc->mutr_top_rig[1]), list->weight);
            if (hyc->mutr_top_rig[1] > 0 && cfr > 0.75) Hist::Head("hPROF_HC_mass_cnt_ISS_CF")->fillH1D(std::abs(hyc->mutr_top_rig[1]), list->weight);

            if (g4mc != nullptr) Hist::Head("hPROF_HC_mass_cnt_MC")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
            if (g4mc != nullptr) Hist::Head("hPROF_HC_mass_cnt_MC_FLUX10")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux10 * list->weight);
            if (g4mc != nullptr) Hist::Head("hPROF_HC_mass_cnt_MC_FLUX27")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux27 * list->weight);
            if (g4mc != nullptr) Hist::Head("hPROF_HC_mass_cnt_MC_FLUX10_CF")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), wgtcf * mc_flux10 * list->weight);
            if (g4mc != nullptr) Hist::Head("hPROF_HC_mass_cnt_MC_FLUX27_CF")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), wgtcf * mc_flux27 * list->weight);
        }
    }


    // RICH
    bool is_all_inn = 
        (trk->lay[1] > 0 && 
         trk->lay[2] == 3 &&
         trk->lay[3] == 3 &&
         trk->lay[4] == 3 &&
         trk->lay[5] == 3 &&
         trk->lay[6] == 3 &&
         trk->lay[7] == 3);
    is_all_inn = true;

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

    if (is_all_inn &&
        trk->ck_status[0] && 
        std::log(trk->ck_nchi[0][0]) < 2.0 &&
        std::log(trk->ck_nchi[0][1]) < 2.0 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        is_official_rich_pr_ok &&
        //is_clean_tk &&
        phtrg
        ) {
        
        double sign = (trk->ck_rig[0] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * std::sqrt((trk->ck_rig[0] * (1.0 / rich->beta + 1.0)) * (trk->ck_rig[0] * (1.0 / rich->beta - 1.0)));
        double cfr  = CheckType(Type::ISS) ? std::abs(trk->ck_rig[0]/rti->max_IGRF) : 0.0;

        double wgtcf = (*Hist::Head("hPROFcf"))()->GetBinContent( (*Hist::Head("hPROFcf"))()->FindBin(std::abs(trk->ck_rig[0])/0.75) );

        if (rich->beta > 0.96 && rich->beta < 0.98) {
            Hist::Head("hPROF_CK_RH_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hPROF_CK_RH_mass_ISS_CF")->fillH1D(mass, list->weight);
            
            Hist::Head("hPROF_CK_RH_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hPROF_CK_RH_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hPROF_CK_RH_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hPROF_CK_RH_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
        }
    }
    
    if (is_all_inn &&
        trk->kf_status[0] && 
        std::log(trk->kf_nchi[0][0]) < 2.0 &&
        std::log(trk->kf_nchi[0][1]) < 2.0 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        is_official_rich_pr_ok &&
        //is_clean_tk &&
        phtrg
        ) {
        
        double sign = (trk->kf_cen_rig[0] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * std::sqrt((trk->kf_cen_rig[0] * (1.0 / rich->beta + 1.0)) * (trk->kf_cen_rig[0] * (1.0 / rich->beta - 1.0)));
        double cfr  = CheckType(Type::ISS) ? std::abs(trk->kf_top_rig[0]/rti->max_IGRF) : 0.0;
        
        double wgtcf = (*Hist::Head("hPROFcf"))()->GetBinContent( (*Hist::Head("hPROFcf"))()->FindBin(std::abs(trk->kf_top_rig[0])/0.75) );

        if (rich->beta > 0.96 && rich->beta < 0.98) {
            Hist::Head("hPROF_KF_RH_mass_ISS")->fillH1D(mass, list->weight);
            if (cfr > 0.75) Hist::Head("hPROF_KF_RH_mass_ISS_CF")->fillH1D(mass, list->weight);

            Hist::Head("hPROF_KF_RH_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
            Hist::Head("hPROF_KF_RH_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
            Hist::Head("hPROF_KF_RH_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
            Hist::Head("hPROF_KF_RH_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
        }
    }
    
    if (is_all_inn &&
        hyc->geom_status[0] &&
        std::log(hyc->geom_nchi_x[0]) < 1.75 &&
        std::log(hyc->geom_nchi_y[0]) < 1.75 &&
        hyc->vel_status[2] &&
        std::log(hyc->vel_nchi[2]) < 2.00 &&
        hyc->mutr_status[2] && 
        std::log(hyc->mutr_nchi_x[2]) < 1.50 &&
        std::log(hyc->mutr_nchi_y[2]) < 1.50 &&
        std::log(hyc->mutr_nchi_b[2]) < 1.50 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        is_rich_pr_ok
        //is_rich_pr_ok &&
        //is_clean_tk
        ) {
        double sign = (hyc->mutr_top_rig[2] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * hyc->mutr_mass[2];
        double cfr  = CheckType(Type::ISS) ? std::abs(hyc->mutr_top_rig[2]/rti->max_IGRF) : 0.0;
        
        double wgtcf = (*Hist::Head("hPROFcf"))()->GetBinContent( (*Hist::Head("hPROFcf"))()->FindBin(std::abs(hyc->mutr_top_rig[2])/0.75) );
        
        bool pure_mc = CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3]);
            
        bool is_new_utime = false;
        if (CheckType(Type::ISS) && utime_cur != list->utime) { utime_cur = list->utime; is_new_utime = true; }
        TH1D* hexpt = (TH1D*) (*Hist::Head("hPROF_HC_expt_gISS"))();
        int   hexpt_rbin = CheckType(Type::ISS) ? hexpt->FindBin(rti->max_IGRF) : 0;
        
        //if (phtrg && hyc->mutr_top_bta[2] > 0.96 && hyc->mutr_top_bta[2] < 0.98) {
        if (phtrg && rich->self_cld_cbta > 0.96 && rich->self_cld_cbta < 0.98) {
            //if (trd->num_vtx[0] <= 4 && 
            //    trd->num_vtx[1] <= 4 &&
            //    rich->self_nhit_other_inn <= 1 &&
            //    rich->self_nhit_other_out <= 3 &&
            //    rich->self_cld_misjudge < 0.1 &&
            //    rich->self_cld_expnpe > 0.1) {
                Hist::Head("hPROF_HC_RH_mass_ISS")->fillH1D(mass, list->weight);
                if (cfr > 0.75) Hist::Head("hPROF_HC_RH_mass_ISS_CF")->fillH1D(mass, list->weight);

                Hist::Head("hPROF_HC_RH_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
                Hist::Head("hPROF_HC_RH_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
                Hist::Head("hPROF_HC_RH_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
                Hist::Head("hPROF_HC_RH_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
                
                if (pure_mc) Hist::Head("hPROF_HC_RH_mass_MC_PURE_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
                if (pure_mc) Hist::Head("hPROF_HC_RH_mass_MC_PURE_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
                if (pure_mc) Hist::Head("hPROF_HC_RH_mass_MC_PURE_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
                if (pure_mc) Hist::Head("hPROF_HC_RH_mass_MC_PURE_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
            
            //}
        
            if (is_official_rich_pr_ok && trk->ck_status[0] && std::abs(rich->self_cld_cbta-rich->beta) < 0.004 && rich->cstcb < 0.004 && rich->npmt >= 3) {
                Hist::Head("hPROF2_HC_RH_mass_ISS")->fillH1D(mass, list->weight);
                if (cfr > 0.75) Hist::Head("hPROF2_HC_RH_mass_ISS_CF")->fillH1D(mass, list->weight);

                Hist::Head("hPROF2_HC_RH_mass_MC_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
                Hist::Head("hPROF2_HC_RH_mass_MC_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
                Hist::Head("hPROF2_HC_RH_mass_MC_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
                Hist::Head("hPROF2_HC_RH_mass_MC_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
                
                if (pure_mc) Hist::Head("hPROF2_HC_RH_mass_MC_PURE_FLUX10")->fillH1D(mass, mc_flux10 * list->weight);
                if (pure_mc) Hist::Head("hPROF2_HC_RH_mass_MC_PURE_FLUX27")->fillH1D(mass, mc_flux27 * list->weight);
                if (pure_mc) Hist::Head("hPROF2_HC_RH_mass_MC_PURE_FLUX10_CF")->fillH1D(mass, wgtcf * mc_flux10 * list->weight);
                if (pure_mc) Hist::Head("hPROF2_HC_RH_mass_MC_PURE_FLUX27_CF")->fillH1D(mass, wgtcf * mc_flux27 * list->weight);
                
                if (mass < 0) clone_tree->Fill();
            }
            
            Hist::Head("hPROF_HC_RH_mass_lchix_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_x[2]), list->weight);
            Hist::Head("hPROF_HC_RH_mass_lchiy_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_y[2]), list->weight);
            Hist::Head("hPROF_HC_RH_mass_lchib_ISS")->fillH2D(mass, std::log(hyc->mutr_nchi_b[2]), list->weight);
            
            Hist::Head("hPROF_HC_RH_mass_inn_ISS")->fillH2D(mass, rich->self_nhit_other_inn, list->weight);
            Hist::Head("hPROF_HC_RH_mass_out_ISS")->fillH2D(mass, rich->self_nhit_other_out, list->weight);
            
            Hist::Head("hPROF_HC_RH_mass_nvtxx_ISS")->fillH2D(mass, trd->num_vtx[0], list->weight);
            Hist::Head("hPROF_HC_RH_mass_nvtxy_ISS")->fillH2D(mass, trd->num_vtx[1], list->weight);
            
            Hist::Head("hPROF_HC_RH_mass_npmt_ISS")->fillH2D(mass, rich->self_cld_nhit, list->weight);
            Hist::Head("hPROF_HC_RH_mass_nhit_ISS")->fillH2D(mass, rich->self_cld_npmt, list->weight);
            
            Hist::Head("hPROF_HC_RH_mass_tof_ISS")->fillH2D(mass, hyc->vel_top_bta[1], list->weight);
            Hist::Head("hPROF_HC_RH_mass_tof2_ISS")->fillH2D(mass, tof->beta, list->weight);
            Hist::Head("hPROF_HC_RH_mass_rh_ISS")->fillH2D(mass, hyc->vel_top_bta[2], list->weight);
            if (is_official_rich_pr_ok) Hist::Head("hPROF_HC_RH_mass_rh2_ISS")->fillH2D(mass, rich->beta, list->weight);
            
            Hist::Head("hPROF_HC_RH_mass_tofq_ISS")->fillH2D(mass, (tof->Q[0]+tof->Q[1])/(tof->Q[2]+tof->Q[3]), list->weight);
            
            Hist::Head("hPROF_HC_RH_mass_npe_ISS")->fillH2D(mass, rich->self_cld_npe, list->weight);
            Hist::Head("hPROF_HC_RH_mass_chi_ISS")->fillH2D(mass, rich->self_cld_nchi, list->weight);
            
            Hist::Head("hPROF_HC_RH_mass_mij_ISS")->fillH2D(mass, rich->self_cld_misjudge, list->weight);
            Hist::Head("hPROF_HC_RH_mass_exp_ISS")->fillH2D(mass, rich->self_cld_expnpe, list->weight);
            if (rich->self_cld_expnpe > 0.1) Hist::Head("hPROF_HC_RH_mass_chg_ISS")->fillH2D(mass, std::sqrt(rich->self_cld_npe / rich->self_cld_expnpe), list->weight);
            Hist::Head("hPROF_HC_RH_mass_bdr_ISS")->fillH2D(mass, rich->self_cld_border, list->weight);
            Hist::Head("hPROF_HC_RH_mass_trc_ISS")->fillH2D(mass, rich->self_cld_trace, list->weight);
            Hist::Head("hPROF_HC_RH_mass_acc_ISS")->fillH2D(mass, rich->self_cld_accuracy, list->weight);
        }
    }

    if (!phtrg) return false;
    const double thres_rig = 100.0;

    if (trk->ck_status[0] && 
        trk->ck_status[1] && 
        trk->ck_status[2] && 
        trk->ck_status[3]) {
        double rvar = thres_rig * (1.0/trk->ck_rig[2] - 1.0/trk->ck_rig[1]);
        if (std::fabs(trk->ck_rig[3]) > thres_rig) Hist::Head("hPROF_CK_rvar")->fillH1D(rvar, list->weight);
        if (g4mc != nullptr && std::fabs(trk->ck_rig[3]) > thres_rig) Hist::Head("hPROF_CK_rvar_MC_FLUX10")->fillH1D(rvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && std::fabs(trk->ck_rig[3]) > thres_rig) Hist::Head("hPROF_CK_rvar_MC_FLUX27")->fillH1D(rvar, mc_flux27 * list->weight);
        if (std::log(trk->ck_nchi[1][0]) < 1.75 &&
            std::log(trk->ck_nchi[1][1]) < 1.75 &&
            std::log(trk->ck_nchi[2][0]) < 1.75 &&
            std::log(trk->ck_nchi[2][1]) < 1.75) {
            if (std::fabs(trk->ck_rig[3]) > thres_rig) Hist::Head("hPROF_CK_rvar_cut")->fillH1D(rvar, list->weight);
            if (g4mc != nullptr && std::fabs(trk->ck_rig[3]) > thres_rig) Hist::Head("hPROF_CK_rvar_cut_MC_FLUX10")->fillH1D(rvar, mc_flux10 * list->weight);
            if (g4mc != nullptr && std::fabs(trk->ck_rig[3]) > thres_rig) Hist::Head("hPROF_CK_rvar_cut_MC_FLUX27")->fillH1D(rvar, mc_flux27 * list->weight);
        }
    }
    if (trk->kf_status[0] && 
        trk->kf_status[1] && 
        trk->kf_status[2] && 
        trk->kf_status[3]) {
        double rvar = thres_rig * (1.0/trk->kf_top_rig[2] - 1.0/trk->kf_top_rig[1]);
        if (std::fabs(trk->kf_top_rig[3]) > thres_rig) Hist::Head("hPROF_KF_rvar")->fillH1D(rvar, list->weight);
        if (g4mc != nullptr && std::fabs(trk->kf_top_rig[3]) > thres_rig) Hist::Head("hPROF_KF_rvar_MC_FLUX10")->fillH1D(rvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && std::fabs(trk->kf_top_rig[3]) > thres_rig) Hist::Head("hPROF_KF_rvar_MC_FLUX27")->fillH1D(rvar, mc_flux27 * list->weight);
        if (std::log(trk->kf_nchi[1][0]) < 1.75 &&
            std::log(trk->kf_nchi[1][1]) < 1.75 &&
            std::log(trk->kf_nchi[2][0]) < 1.75 &&
            std::log(trk->kf_nchi[2][1]) < 1.75) {
            if (std::fabs(trk->kf_top_rig[3]) > thres_rig) Hist::Head("hPROF_KF_rvar")->fillH1D(rvar, list->weight);
            if (g4mc != nullptr && std::fabs(trk->kf_top_rig[3]) > thres_rig) Hist::Head("hPROF_KF_rvar_cut_MC_FLUX10")->fillH1D(rvar, mc_flux10 * list->weight);
            if (g4mc != nullptr && std::fabs(trk->kf_top_rig[3]) > thres_rig) Hist::Head("hPROF_KF_rvar_cut_MC_FLUX27")->fillH1D(rvar, mc_flux27 * list->weight);
        }
    }
    if (hyc->geom_status[0] &&
        hyc->geom_status[1] && 
        hyc->geom_status[2] && 
        hyc->geom_status[3]) {
        double rvar = thres_rig * (1.0/hyc->geom_top_rig[2] - 1.0/hyc->geom_top_rig[1]);
        if (std::fabs(hyc->geom_top_rig[3]) > thres_rig) Hist::Head("hPROF_HC_rvar")->fillH1D(rvar, list->weight);
        if (g4mc != nullptr && std::fabs(hyc->geom_top_rig[3]) > thres_rig) Hist::Head("hPROF_HC_rvar_MC_FLUX10")->fillH1D(rvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && std::fabs(hyc->geom_top_rig[3]) > thres_rig) Hist::Head("hPROF_HC_rvar_MC_FLUX27")->fillH1D(rvar, mc_flux27 * list->weight);
        if (std::log(hyc->geom_nchi_x[1]) < 1.75 &&
            std::log(hyc->geom_nchi_y[1]) < 1.75 &&
            std::log(hyc->geom_nchi_x[2]) < 1.75 &&
            std::log(hyc->geom_nchi_y[2]) < 1.75) {
            if (std::fabs(hyc->geom_top_rig[3]) > thres_rig) Hist::Head("hPROF_HC_rvar_cut")->fillH1D(rvar, list->weight);
            if (g4mc != nullptr && std::fabs(hyc->geom_top_rig[3]) > thres_rig) Hist::Head("hPROF_HC_rvar_cut_MC_FLUX10")->fillH1D(rvar, mc_flux10 * list->weight);
            if (g4mc != nullptr && std::fabs(hyc->geom_top_rig[3]) > thres_rig) Hist::Head("hPROF_HC_rvar_cut_MC_FLUX27")->fillH1D(rvar, mc_flux27 * list->weight);
        }
    }
    
    return true;
}

bool Analyzer::process_data_l() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
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
    if (rich->self_nhit_other_inn > 1) return false;
    if (rich->self_nhit_other_out > 3) return false;
    if (rich->self_trace < 0.30) return false;

    if (rich->self_stn_status && rich->self_stn_dist > 3.4) return false;
    if (rich->self_cld_status && rich->self_cld_misjudge > 0.1) return false; 
    
    // Chi-square Cut
    const double vel_nchi_cut = 2.00;
    const double trk_nchi_cut = 1.75;
    
    // Velocity
    if (!hyc->vel_status[1]) return false;
    if (std::log(hyc->vel_nchi[1]) > vel_nchi_cut) return false;
   
    // Geometry
    if (trk->lay[1] == 0) return false;
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
    
    // Variables
    double trdllr  = trd->tdLLR_ep;
    double sqrm    = hyc->mutr_sqrm[1];

    double bta     = hyc->phys_top_bta[0][1];
    double rig     = hyc->phys_top_rig[0][1];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
    
    bool is_like_pr = (trd->tdLLR_ep > 0.75 && !rich->self_cld_status);
    bool is_like_el = (trd->tdLLR_ep < 0.65 &&  rich->self_cld_status);
    bool is_like_pi = (rich->self_cld_status && 
                       hyc->mutr_status[2] &&
                       std::log(hyc->mutr_nchi_x[2]) < trk_nchi_cut &&
                       std::log(hyc->mutr_nchi_y[2]) < trk_nchi_cut &&
                       std::log(hyc->mutr_nchi_b[2]) < vel_nchi_cut &&
                       hyc->mutr_sqrm[2] < 0.4);
    
    // Cut on in chi-square
    const double nchi_geom_cut = 0.75 + 0.5 * std::erfc((abs_rig - 25.0) / 20.0);
    if (std::log(hyc->geom_nchi_y[0]) > nchi_geom_cut) return false;
   
    // TRD extra hit
    int nvtx_x = 1 + static_cast<int>(3.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig));
    int nvtx_y = 3 + 2 * static_cast<int>(3.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig));
    if (trd->num_vtx[0] > nvtx_x) return false;
    if (trd->num_vtx[1] > nvtx_y) return false;
    
    // Cutoff
    if (CheckType(Type::ISS) && cfr < 0.75) return false;
    
    if (signr > 0 && !rich->self_cld_status) Hist::Head("hLP_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && !rich->self_cld_status) Hist::Head("hLN_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hLP_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hLN_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_pi) Hist::Head("hLN_sqrm_pi")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && is_like_pr && cfr > 0.80) Hist::Head("hLP_cf80_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.80) Hist::Head("hLN_cf80_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_pi && cfr > 0.80) Hist::Head("hLN_cf80_sqrm_pi")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.85) Hist::Head("hLP_cf85_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.85) Hist::Head("hLN_cf85_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_pi && cfr > 0.85) Hist::Head("hLN_cf85_sqrm_pi")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.90) Hist::Head("hLP_cf90_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.90) Hist::Head("hLN_cf90_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_pi && cfr > 0.90) Hist::Head("hLN_cf90_sqrm_pi")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTLP_sqrm_pr")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTLN_sqrm_el")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_pi) Hist::Head("hTLN_sqrm_pi")->fillH3D(list->utime, abs_rig, sqrm, list->weight);

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

    bool is_sel_pr = trd->tdLLR_ep > 0.75 && (!rich->self_cld_status || is_rich_pr);

    if (!is_sel_pr) return false;

    if (signr > 0) Hist::Head("hLP_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hLN_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && cfr > 0.80) Hist::Head("hLP_cf80_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.80) Hist::Head("hLN_cf80_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && cfr > 0.85) Hist::Head("hLP_cf85_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.85) Hist::Head("hLN_cf85_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && cfr > 0.90) Hist::Head("hLP_cf90_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.90) Hist::Head("hLN_cf90_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0) Hist::Head("hTLP_sqrm")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hTLN_sqrm")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    
    if (signr > 0) Hist::Head("hLP_cnt")->fillH1D(abs_rig, list->weight);
    if (signr < 0) Hist::Head("hLN_cnt")->fillH1D(abs_rig, list->weight);
    
    if (g4mc != nullptr) Hist::Head("hL_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
    if (g4mc != nullptr) Hist::Head("hTL_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
  
    return true;
}

bool Analyzer::process_data_m() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
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
    if (rich->self_trace < 0.30) return false;
    
    if (rich->self_stn_status && rich->self_stn_dist > 3.4) return false;
    if (rich->self_cld_status && rich->self_cld_misjudge > 0.1) return false; 
   
    if (!rich->self_cld_status) return false;
    
    // Chi-square Cut
    const double vel_nchi_cut = 2.00;
    const double trk_nchi_cut = 1.75;
    
    // Velocity
    if (!hyc->vel_status[2]) return false;
    if (std::log(hyc->vel_nchi[2]) > vel_nchi_cut) return false;
    
    // Geometry
    if (trk->lay[1] == 0) return false;
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

    // Variables
    double trdllr  = trd->tdLLR_ep;
    double sqrm    = hyc->mutr_sqrm[2];
    
    double rig     = hyc->phys_top_rig[0][2];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;

    bool is_like_pr = (trdllr > 0.75);
    bool is_like_el = (trdllr < 0.65);
    
    // Cut on in chi-square
    const double nchi_geom_cut = 0.75 + 0.5 * std::erfc((abs_rig - 25.0) / 20.0);
    if (std::log(hyc->geom_nchi_y[0]) > nchi_geom_cut) return false;

    // Cutoff
    if (CheckType(Type::ISS) && cfr < 0.75) return false;
    
    if (signr > 0) Hist::Head("hMP_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0) Hist::Head("hMN_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hMP_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hMN_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && is_like_pr && cfr > 0.80) Hist::Head("hMP_cf80_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.80) Hist::Head("hMN_cf80_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.85) Hist::Head("hMP_cf85_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.85) Hist::Head("hMN_cf85_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.90) Hist::Head("hMP_cf90_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.90) Hist::Head("hMN_cf90_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTMP_sqrm_pr")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTMN_sqrm_el")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
  
    if (trdllr < 0.75) return false;

    if (signr > 0) Hist::Head("hMP_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hMN_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && cfr > 0.80) Hist::Head("hMP_cf80_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.80) Hist::Head("hMN_cf80_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && cfr > 0.85) Hist::Head("hMP_cf85_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.85) Hist::Head("hMN_cf85_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && cfr > 0.90) Hist::Head("hMP_cf90_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.90) Hist::Head("hMN_cf90_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0) Hist::Head("hTMP_sqrm")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hTMN_sqrm")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
  
    if (signr > 0) Hist::Head("hMP_cnt")->fillH1D(abs_rig, list->weight);
    if (signr < 0) Hist::Head("hMN_cnt")->fillH1D(abs_rig, list->weight);
    
    if (g4mc != nullptr) Hist::Head("hM_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
    if (g4mc != nullptr) Hist::Head("hTM_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
    
    return true;
}

bool Analyzer::process_data_i() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    // Trigger
    if ((trg->bit&8) != 8) return false;

    // Chi-square Cut
    const double vel_nchi_cut = 2.00;
    const double trk_nchi_cut = 1.75;

    // Geometry
    if (trk->lay[1] == 0) return false;
    if (!hyc->geom_status[0]) return false;
    
    // TOF
    if (tof->extcls_noise != 0) return false;
    if (tof->num_in_time_cls > 4) return false;

    // Variables
    double trdllr  = trd->tdLLR_ep;
    
    double rig     = hyc->geom_top_rig[0];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;

    // Cutoff
    if (CheckType(Type::ISS) && cfr < 0.75) return false;
    
    // Cut on in chi-square
    double lchix = std::log(hyc->geom_nchi_x[0]);
    double lchiy = std::log(hyc->geom_nchi_y[0]);

    if (lchix > trk_nchi_cut) return false;
    
    if (signr > 0) Hist::Head("hIP_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr < 0) Hist::Head("hIN_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    
    const double nchi_geom_cut = 0.75 + 0.5 * std::erfc((abs_rig - 25.0) / 20.0);
    if (lchiy > nchi_geom_cut) return false;
  
    bool is_like_pr = (ecal->status && ecal->mvaBDT < -0.6);
    bool is_like_el = (ecal->status && ecal->mvaBDT >  0.2);
    
    if (signr > 0 && is_like_pr) Hist::Head("hIP_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hIN_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr && cfr > 0.80) Hist::Head("hIP_cf80_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.80) Hist::Head("hIN_cf80_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.85) Hist::Head("hIP_cf85_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.85) Hist::Head("hIN_cf85_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.90) Hist::Head("hIP_cf90_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.90) Hist::Head("hIN_cf90_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTIP_llr_pr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTIN_llr_el")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
   
    // RICH 
    if (rich->self_status && rich->self_kind == 1 && hyc->mutr_status[2]) {
        if (std::log(hyc->mutr_nchi_x[2]) > trk_nchi_cut) return false;
        if (std::log(hyc->mutr_nchi_y[2]) > trk_nchi_cut) return false;
        if (std::log(hyc->mutr_nchi_b[2]) > vel_nchi_cut) return false;
    }
    if (rich->self_status && rich->self_kind == 1 && hyc->phys_status[0][2]) {
        if (std::log(hyc->phys_nchi_x[0][2]) > trk_nchi_cut) return false;
        if (std::log(hyc->phys_nchi_y[0][2]) > trk_nchi_cut) return false;
        if (std::log(hyc->phys_nchi_b[0][2]) > vel_nchi_cut) return false;
    }
    
    if (signr > 0) Hist::Head("hIP_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0) Hist::Head("hIN_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && cfr > 0.80) Hist::Head("hIP_cf80_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && cfr > 0.80) Hist::Head("hIN_cf80_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr > 0 && cfr > 0.85) Hist::Head("hIP_cf85_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && cfr > 0.85) Hist::Head("hIN_cf85_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr > 0 && cfr > 0.90) Hist::Head("hIP_cf90_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && cfr > 0.90) Hist::Head("hIN_cf90_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0) Hist::Head("hTIP_llr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr < 0) Hist::Head("hTIN_llr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    
    if (signr > 0) Hist::Head("hIP_cnt")->fillH1D(abs_rig, list->weight);
    if (signr < 0) Hist::Head("hIN_cnt")->fillH1D(abs_rig, list->weight);

    if (g4mc != nullptr) Hist::Head("hI_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
    if (g4mc != nullptr) Hist::Head("hTI_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);

    return true;
}

bool Analyzer::process_data_h() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_flux10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_flux27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        double crr_flux_pr = gPROFpr27->Eval(g4mc->prm_mom, 0, "");
        if (crr_flux_pr < 0) crr_flux_pr = 0.0;

        double crr_flux_ap = crr_flux_pr * gPROFap2pr->Eval(g4mc->prm_mom, 0, "");
        if (crr_flux_ap < 0) crr_flux_ap = 0.0;

        mc_flux10 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap);
        mc_flux27 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    // Variables
    double trdllr = trd->tdLLR_ep;
  
    // MDR: Proton
    const std::array<double, 2> frso_pr_in_ck({ 1.03308e-01, 2.74374e+02 });
    const std::array<double, 2> frso_pr_l1_ck({ 2.16622e-01, 6.64065e+02 });
    const std::array<double, 2> frso_pr_l9_ck({ 1.46456e-01, 8.72824e+02 });
    const std::array<double, 2> frso_pr_fs_ck({ 1.47755e-01, 1.82999e+03 });
    
    const std::array<double, 2> frso_pr_in_hc({ 1.03081e-01, 2.96552e+02 });
    const std::array<double, 2> frso_pr_l1_hc({ 2.17407e-01, 7.29381e+02 });
    const std::array<double, 2> frso_pr_l9_hc({ 1.49239e-01, 9.52511e+02 });
    const std::array<double, 2> frso_pr_fs_hc({ 1.47495e-01, 1.92807e+03 });

    // MDR: Helium
    const std::array<double, 2> frso_he_in_ck({ 1.09510e-01, 6.00266e+02 });
    const std::array<double, 2> frso_he_l1_ck({ 2.24697e-01, 1.23937e+03 });
    const std::array<double, 2> frso_he_l9_ck({ 1.44779e-01, 1.57583e+03 });
    const std::array<double, 2> frso_he_fs_ck({ 1.71372e-01, 2.70841e+03 });

    const std::array<double, 2> frso_he_in_hc({ 1.07304e-01, 6.54203e+02 });
    const std::array<double, 2> frso_he_l1_hc({ 2.29644e-01, 1.37133e+03 });
    const std::array<double, 2> frso_he_l9_hc({ 1.44710e-01, 1.69553e+03 });
    const std::array<double, 2> frso_he_fs_hc({ 1.71859e-01, 2.80689e+03 });

    bool ptfs = hyc->geom_status[0] && hyc->geom_status[1] && hyc->geom_status[2] && hyc->geom_status[3] && (trk->lay[0] == 3) && (trk->lay[8] == 3);
    bool ptl9 = hyc->geom_status[0] && hyc->geom_status[2] && (trk->lay[8] == 3);
    bool ptl1 = hyc->geom_status[0] && hyc->geom_status[1] && (trk->lay[0] == 3 && trk->lay[1] > 0);

    bool swfs = ptfs;
    bool swl9 = ptl9 && !ptfs;
    bool swl1 = ptl1 && !ptfs;
   
    if (!swfs && !swl9 && !swl1) return false;
    int num_tk = (trk->ext_num_hit[2] > 0 && trk->ext_num_hit[3] > 0) + 
                 (trk->ext_num_hit[4] > 0 && trk->ext_num_hit[5] > 0) + 
                 (trk->ext_num_hit[6] > 0 && trk->ext_num_hit[7] > 0);

    // Chi-square Cut
    const double nchix_cut = 3.00;
    //const double nchiy_cut = 2.75;
    const double nchiy_cut = 2.50;

    while (swl1) {
        double rig     = hyc->geom_top_rig[1];
        double abs_rig = std::abs(rig);
        short  signr   = (rig > 0) ? 1 : -1;
    
        double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
        double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
        
        double lchix_in = std::log(hyc->geom_nchi_x[0]);
        double lchiy_in = std::log(hyc->geom_nchi_y[0]);

        double lchix_l1 = std::log(hyc->geom_nchi_x[1]);
        double lchiy_l1 = std::log(hyc->geom_nchi_y[1]);

        double lchix = (9.0/25.0) * lchix_in + (16.0/25.0) * lchix_l1;
        double lchiy = (9.0/25.0) * lchiy_in + (16.0/25.0) * lchiy_l1;
        
        double dchix_l1 = lchix_l1 - lchix_in;
        double dchiy_l1 = lchiy_l1 - lchiy_in;
        
        double lrvar = std::log(1.0 + std::abs(hyc->geom_top_rig[1] / hyc->geom_top_rig[0] - 1.0));

        // Cutoff
        if (CheckType(Type::ISS) && cfr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPl1_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (signr < 0) Hist::Head("hHNl1_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (trdllr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPl1_lchix_in")->fillH2D(abs_rig, lchix_in, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchix_in")->fillH2D(abs_rig, lchix_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchix_l1")->fillH2D(abs_rig, lchix_l1, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchix_l1")->fillH2D(abs_rig, lchix_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchix")->fillH2D(abs_rig, lchix, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchix")->fillH2D(abs_rig, lchix, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_dchix_l1")->fillH2D(abs_rig, dchix_l1, list->weight);
        if (signr < 0) Hist::Head("hHNl1_dchix_l1")->fillH2D(abs_rig, dchix_l1, list->weight);

        if (lchix_in > nchix_cut) break;
        if (lchix_l1 > nchix_cut) break;

        if (signr > 0) Hist::Head("hHPl1_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy")->fillH2D(abs_rig, lchiy, list->weight);

        if (signr > 0) Hist::Head("hHPl1_dchiy_l1")->fillH2D(abs_rig, dchiy_l1, list->weight);
        if (signr < 0) Hist::Head("hHNl1_dchiy_l1")->fillH2D(abs_rig, dchiy_l1, list->weight);
        
        if (lchiy_in > nchiy_cut) break;
        if (lchiy_l1 > nchiy_cut) break;

        if (signr > 0) Hist::Head("hHPl1_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        
        double lrvar_cut = std::sqrt(0.05*0.05 + 1.0 * std::log(1.0+abs_rig/100.0) * std::log(1.0+abs_rig/100.0));
        if (lrvar > 1.0) break;
        
        if (signr > 0) Hist::Head("hHPl1_num_tk")->fillH2D(abs_rig, num_tk, list->weight);
        if (signr < 0) Hist::Head("hHNl1_num_tk")->fillH2D(abs_rig, num_tk, list->weight);

        if (signr > 0) Hist::Head("hHPl1_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        
        if (g4mc != nullptr && signr > 0) Hist::Head("hHPl1_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && signr < 0) Hist::Head("hHNl1_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && signr > 0) Hist::Head("hHPl1_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
        if (g4mc != nullptr && signr < 0) Hist::Head("hHNl1_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_cnt")->fillH1D(abs_rig, list->weight);
        if (signr < 0) Hist::Head("hHNl1_cnt")->fillH1D(abs_rig, list->weight);
        
        if (g4mc != nullptr) Hist::Head("hHl1_cnt_MC_FLUX10")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux10 * list->weight);
        if (g4mc != nullptr) Hist::Head("hHl1_cnt_MC_FLUX27")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux27 * list->weight);

        break;
    }

    while (swl9) {
        double rig     = hyc->geom_top_rig[2];
        double abs_rig = std::abs(rig);
        short  signr   = (rig > 0) ? 1 : -1;
    
        double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
        double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
        
        double lchix_in = std::log(hyc->geom_nchi_x[0]);
        double lchiy_in = std::log(hyc->geom_nchi_y[0]);
        
        double lchix_l9 = std::log(hyc->geom_nchi_x[2]);
        double lchiy_l9 = std::log(hyc->geom_nchi_y[2]);
        
        double lchix = (1.0/5.0) * lchix_in + (4.0/5.0) * lchix_l9;
        double lchiy = (1.0/5.0) * lchiy_in + (4.0/5.0) * lchiy_l9;
        
        double dchix_l9 = lchix_l9 - lchix_in;
        double dchiy_l9 = lchiy_l9 - lchiy_in;
        
        double lrvar = std::log(1.0 + std::abs(hyc->geom_top_rig[2] / hyc->geom_top_rig[0] - 1.0));

        // Cutoff
        if (CheckType(Type::ISS) && cfr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPl9_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (signr < 0) Hist::Head("hHNl9_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (trdllr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPl9_lchix_in")->fillH2D(abs_rig, lchix_in, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchix_in")->fillH2D(abs_rig, lchix_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchix_l9")->fillH2D(abs_rig, lchix_l9, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchix_l9")->fillH2D(abs_rig, lchix_l9, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchix")->fillH2D(abs_rig, lchix, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchix")->fillH2D(abs_rig, lchix, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_dchix_l9")->fillH2D(abs_rig, dchix_l9, list->weight);
        if (signr < 0) Hist::Head("hHNl9_dchix_l9")->fillH2D(abs_rig, dchix_l9, list->weight);

        if (lchix_in > nchix_cut) break;
        if (lchix_l9 > nchix_cut) break;

        if (signr > 0) Hist::Head("hHPl9_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_dchiy_l9")->fillH2D(abs_rig, dchiy_l9, list->weight);
        if (signr < 0) Hist::Head("hHNl9_dchiy_l9")->fillH2D(abs_rig, dchiy_l9, list->weight);
        
        if (lchiy_in > nchiy_cut) break;
        if (lchiy_l9 > nchiy_cut) break;
        
        if (signr > 0) Hist::Head("hHPl9_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        
        double lrvar_cut = std::sqrt(0.06*0.06 + 1.0 * std::log(1.0+abs_rig/80.0) * std::log(1.0+abs_rig/80.0));
        if (lrvar > 1.0) break;
        
        if (signr > 0) Hist::Head("hHPl9_num_tk")->fillH2D(abs_rig, num_tk, list->weight);
        if (signr < 0) Hist::Head("hHNl9_num_tk")->fillH2D(abs_rig, num_tk, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        
        if (g4mc != nullptr && signr > 0) Hist::Head("hHPl9_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && signr < 0) Hist::Head("hHNl9_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && signr > 0) Hist::Head("hHPl9_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
        if (g4mc != nullptr && signr < 0) Hist::Head("hHNl9_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_cnt")->fillH1D(abs_rig, list->weight);
        if (signr < 0) Hist::Head("hHNl9_cnt")->fillH1D(abs_rig, list->weight);
        
        if (g4mc != nullptr) Hist::Head("hHl9_cnt_MC_FLUX10")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux10 * list->weight);
        if (g4mc != nullptr) Hist::Head("hHl9_cnt_MC_FLUX27")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux27 * list->weight);

        break;
    }

    while (swfs) {
        double rig     = hyc->geom_top_rig[3];
        double abs_rig = std::abs(rig);
        short  signr   = (rig > 0) ? 1 : -1;
    
        double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
        double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
        
        double lchix_in = std::log(hyc->geom_nchi_x[0]);
        double lchiy_in = std::log(hyc->geom_nchi_y[0]);
        
        double lchix_l1 = std::log(hyc->geom_nchi_x[1]);
        double lchiy_l1 = std::log(hyc->geom_nchi_y[1]);
        
        double lchix_l9 = std::log(hyc->geom_nchi_x[2]);
        double lchiy_l9 = std::log(hyc->geom_nchi_y[2]);

        double lchix_fs = std::log(hyc->geom_nchi_x[3]);
        double lchiy_fs = std::log(hyc->geom_nchi_y[3]);
        
        double lchix = (9.0/61.0) * lchix_l1 + (16.0/61.0) * lchix_l9 + (36.0/61.0) * lchix_fs;
        double lchiy = (9.0/61.0) * lchiy_l1 + (16.0/61.0) * lchiy_l9 + (36.0/61.0) * lchiy_fs;

        double dchix_l1 = lchix_l1 - lchix_in;
        double dchiy_l1 = lchiy_l1 - lchiy_in;

        double dchix_l9 = lchix_l9 - lchix_in;
        double dchiy_l9 = lchiy_l9 - lchiy_in;

        double dchix_fs = lchix_fs - lchix_in;
        double dchiy_fs = lchiy_fs - lchiy_in;
        
        double lrvar_l1 = std::log(1.0 + std::abs(hyc->geom_top_rig[3] / hyc->geom_top_rig[1] - 1.0));
        double lrvar_l9 = std::log(1.0 + std::abs(hyc->geom_top_rig[3] / hyc->geom_top_rig[2] - 1.0));
        double lrvar_fs = std::log(1.0 + std::abs(hyc->geom_top_rig[3] / hyc->geom_top_rig[2] - hyc->geom_top_rig[3] / hyc->geom_top_rig[1]));

        double lrvar = (9.0/61.0) * lrvar_l1 + (16.0/61.0) * lrvar_l9 + (36.0/61.0) * lrvar_fs;
        
        // Cutoff
        if (CheckType(Type::ISS) && cfr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPfs_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (signr < 0) Hist::Head("hHNfs_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (trdllr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPfs_lchix_l1")->fillH2D(abs_rig, lchix_l1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchix_l1")->fillH2D(abs_rig, lchix_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchix_l9")->fillH2D(abs_rig, lchix_l9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchix_l9")->fillH2D(abs_rig, lchix_l9, list->weight);
       
        if (signr > 0) Hist::Head("hHPfs_lchix_fs")->fillH2D(abs_rig, lchix_fs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchix_fs")->fillH2D(abs_rig, lchix_fs, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchix")->fillH2D(abs_rig, lchix, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchix")->fillH2D(abs_rig, lchix, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_dchix_l1")->fillH2D(abs_rig, dchix_l1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_dchix_l1")->fillH2D(abs_rig, dchix_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_dchix_l9")->fillH2D(abs_rig, dchix_l9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_dchix_l9")->fillH2D(abs_rig, dchix_l9, list->weight);
       
        if (signr > 0) Hist::Head("hHPfs_dchix_fs")->fillH2D(abs_rig, dchix_fs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_dchix_fs")->fillH2D(abs_rig, dchix_fs, list->weight);
        
        if (lchix_l1 > nchix_cut) break;
        if (lchix_l9 > nchix_cut) break;
        if (lchix_fs > nchix_cut) break;
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_fs")->fillH2D(abs_rig, lchiy_fs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_fs")->fillH2D(abs_rig, lchiy_fs, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy")->fillH2D(abs_rig, lchiy, list->weight);

        if (lchiy_l1 > nchiy_cut) break;
        if (lchiy_l9 > nchiy_cut) break;
        if (lchiy_fs > nchiy_cut) break;
        
        if (signr > 0) Hist::Head("hHPfs_lrvar_l1")->fillH2D(abs_rig, lrvar_l1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lrvar_l1")->fillH2D(abs_rig, lrvar_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lrvar_l9")->fillH2D(abs_rig, lrvar_l9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lrvar_l9")->fillH2D(abs_rig, lrvar_l9, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lrvar_fs")->fillH2D(abs_rig, lrvar_fs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lrvar_fs")->fillH2D(abs_rig, lrvar_fs, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_dchiy_l1")->fillH2D(abs_rig, dchiy_l1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_dchiy_l1")->fillH2D(abs_rig, dchiy_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_dchiy_l9")->fillH2D(abs_rig, dchiy_l9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_dchiy_l9")->fillH2D(abs_rig, dchiy_l9, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_dchiy_fs")->fillH2D(abs_rig, dchiy_fs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_dchiy_fs")->fillH2D(abs_rig, dchiy_fs, list->weight);

        //double lrvar_cut = std::sqrt(0.06*0.06 + 1.0 * std::log(1.0+abs_rig/80.0) * std::log(1.0+abs_rig/80.0));
        //if (lrvar_l1 > 1.0) break;
        //if (lrvar_l9 > 1.0) break;
        //if (lrvar_fs > 1.0) break;
        //if (lrvar > 1.0) break;
        
        if (signr > 0) Hist::Head("hHPfs_num_tk")->fillH2D(abs_rig, num_tk, list->weight);
        if (signr < 0) Hist::Head("hHNfs_num_tk")->fillH2D(abs_rig, num_tk, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        
        if (g4mc != nullptr && signr > 0) Hist::Head("hHPfs_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && signr < 0) Hist::Head("hHNfs_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
        if (g4mc != nullptr && signr > 0) Hist::Head("hHPfs_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
        if (g4mc != nullptr && signr < 0) Hist::Head("hHNfs_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);

        if (signr > 0) Hist::Head("hHPfs_cnt")->fillH1D(abs_rig, list->weight);
        if (signr < 0) Hist::Head("hHNfs_cnt")->fillH1D(abs_rig, list->weight);
        
        if (g4mc != nullptr) Hist::Head("hHfs_cnt_MC_FLUX10")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux10 * list->weight);
        if (g4mc != nullptr) Hist::Head("hHfs_cnt_MC_FLUX27")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux27 * list->weight);

        break;
    }

    return true;
}


#endif // __Analyzer_C__
