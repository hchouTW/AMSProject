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


bool Analyzer::build_hist() {
    // test
    TFile* hPROPfile = TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/antip/others/mcflux.root");
    Hist* hPROFwgt = Hist::New("hPROFwgt", (TH1D*) hPROPfile->Get("hPROF_HC_mass_cnt_iss"));
    
    TFile* hPROPapp = TFile::Open("/afs/cern.ch/user/h/hchou/AMSProject/antip/others/ams02_acceptance_corrected_20150330.root");
    Hist* hPROFapp = Hist::New("hPROFapp", (TH1D*) hPROPapp->Get("havg_ratio_vs_rig"));

    //std::cerr << Form("INFO %14.8f %14.8f\n", (*hPROFwgt)()->GetEntries(), (*hPROFwgt)()->Integral());
    //(*hPROFwgt)()->Scale(1.0/(*hPROFwgt)()->GetMaximum());
    //hPROPfile->Close();
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

    Axis AXmom("Momentum [GeV]", vmom);
    Axis AXrig("Rigidity [GV]" , vmom);
    
    Axis AXmom2("Momentum [GeV]", vmom4);
    Axis AXrig2("Rigidity [GV]" , vmom4);
    
    Axis AXrig3("Rigidity [GV]" , 400, 0.6, 800.0, AxisScale::kLog);
    
    Axis AXnchi("nchi", 50, -3,  2.0, AxisScale::kLinear);
    Axis AXnext("next", 50,  0, 50.0, AxisScale::kLinear);
    Axis AXchrg("Z", 400, 0.0, 10.0);
   
    Axis AXPROFmass("Mass/Z [GV/c^{2}]", 350, -3.5, 3.5);
    //Axis AXPROFmass("Mass/Z [GV/c^{2}]", 1000, -10.0, 10.0);
    Hist::New("hPROF_CK_mass_ISS" , AXPROFmass);
    Hist::New("hPROF_CK_mass_ISS2", AXPROFmass);
    Hist::New("hPROF_CK_mass_ISS3", AXPROFmass);
    Hist::New("hPROF_CK_mass_ISS4", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC_FLUX", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC_FLUX2", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC_FLUX3", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC2_FLUX", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC2_FLUX2", AXPROFmass);
    Hist::New("hPROF_CK_mass_MC2_FLUX3", AXPROFmass);
    Hist::New("hPROF_CK_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_ISS2", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_ISS3", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_ISS4", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC_FLUX", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC_FLUX2", HistAxis(AXrig3));
    Hist::New("hPROF_CK_mass_cnt_MC_FLUX3", HistAxis(AXrig3));
    
    Hist::New("hPROF2_CK_mass_ISS" , AXPROFmass);
    Hist::New("hPROF2_CK_mass_ISS2", AXPROFmass);
    Hist::New("hPROF2_CK_mass_ISS3", AXPROFmass);
    Hist::New("hPROF2_CK_mass_ISS4", AXPROFmass);
    Hist::New("hPROF2_CK_mass_MC_FLUX", AXPROFmass);
    Hist::New("hPROF2_CK_mass_MC_FLUX2", AXPROFmass);
    Hist::New("hPROF2_CK_mass_MC_FLUX3", AXPROFmass);
    Hist::New("hPROF2_CK_mass_MC2_FLUX", AXPROFmass);
    Hist::New("hPROF2_CK_mass_MC2_FLUX2", AXPROFmass);
    Hist::New("hPROF2_CK_mass_MC2_FLUX3", AXPROFmass);
    Hist::New("hPROF2_CK_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF2_CK_mass_cnt_ISS2", HistAxis(AXrig3));
    Hist::New("hPROF2_CK_mass_cnt_ISS3", HistAxis(AXrig3));
    Hist::New("hPROF2_CK_mass_cnt_ISS4", HistAxis(AXrig3));
    Hist::New("hPROF2_CK_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF2_CK_mass_cnt_MC_FLUX", HistAxis(AXrig3));
    Hist::New("hPROF2_CK_mass_cnt_MC_FLUX2", HistAxis(AXrig3));
    Hist::New("hPROF2_CK_mass_cnt_MC_FLUX3", HistAxis(AXrig3));
    
    Hist::New("hPROF_KF_mass_ISS" , AXPROFmass);
    Hist::New("hPROF_KF_mass_ISS2", AXPROFmass);
    Hist::New("hPROF_KF_mass_ISS3", AXPROFmass);
    Hist::New("hPROF_KF_mass_ISS4", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC_FLUX", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC_FLUX2", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC_FLUX3", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC2_FLUX", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC2_FLUX2", AXPROFmass);
    Hist::New("hPROF_KF_mass_MC2_FLUX3", AXPROFmass);
    Hist::New("hPROF_KF_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_ISS2", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_ISS3", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_ISS4", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC_FLUX", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC_FLUX2", HistAxis(AXrig3));
    Hist::New("hPROF_KF_mass_cnt_MC_FLUX3", HistAxis(AXrig3));
    
    Hist::New("hPROF2_KF_mass_ISS" , AXPROFmass);
    Hist::New("hPROF2_KF_mass_ISS2", AXPROFmass);
    Hist::New("hPROF2_KF_mass_ISS3", AXPROFmass);
    Hist::New("hPROF2_KF_mass_ISS4", AXPROFmass);
    Hist::New("hPROF2_KF_mass_MC_FLUX", AXPROFmass);
    Hist::New("hPROF2_KF_mass_MC_FLUX2", AXPROFmass);
    Hist::New("hPROF2_KF_mass_MC_FLUX3", AXPROFmass);
    Hist::New("hPROF2_KF_mass_MC2_FLUX", AXPROFmass);
    Hist::New("hPROF2_KF_mass_MC2_FLUX2", AXPROFmass);
    Hist::New("hPROF2_KF_mass_MC2_FLUX3", AXPROFmass);
    Hist::New("hPROF2_KF_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF2_KF_mass_cnt_ISS2", HistAxis(AXrig3));
    Hist::New("hPROF2_KF_mass_cnt_ISS3", HistAxis(AXrig3));
    Hist::New("hPROF2_KF_mass_cnt_ISS4", HistAxis(AXrig3));
    Hist::New("hPROF2_KF_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF2_KF_mass_cnt_MC_FLUX", HistAxis(AXrig3));
    Hist::New("hPROF2_KF_mass_cnt_MC_FLUX2", HistAxis(AXrig3));
    Hist::New("hPROF2_KF_mass_cnt_MC_FLUX3", HistAxis(AXrig3));
    
    Hist::New("hPROF_HC_mass_ISS" , AXPROFmass);
    Hist::New("hPROF_HC_mass_ISS2", AXPROFmass);
    Hist::New("hPROF_HC_mass_ISS3", AXPROFmass);
    Hist::New("hPROF_HC_mass_ISS4", AXPROFmass);
    Hist::New("hPROF_HC_mass_ISS5", AXPROFmass);
    Hist::New("hPROF_HC_mass_ISS6", AXPROFmass);
    Hist::New("hPROF_HC_mass_map_ISS" , HistAxis(AXPROFmass, AXPROFmass));
    Hist::New("hPROF_HC_mass_map_ISS2", HistAxis(AXPROFmass, AXPROFmass));
    Hist::New("hPROF_HC_mass_map_ISS3", HistAxis(AXPROFmass, AXPROFmass));
    Hist::New("hPROF_HC_mass_map_ISS4", HistAxis(AXPROFmass, AXPROFmass));
    Hist::New("hPROF_HC_mass_MC_FLUX", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_FLUX2", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC_FLUX3", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC2_FLUX", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC2_FLUX2", AXPROFmass);
    Hist::New("hPROF_HC_mass_MC2_FLUX3", AXPROFmass);
    Hist::New("hPROF_HC_mass_cnt_ISS" , HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_ISS2", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_ISS3", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_ISS4", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC_FLUX", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC_FLUX2", HistAxis(AXrig3));
    Hist::New("hPROF_HC_mass_cnt_MC_FLUX3", HistAxis(AXrig3));
    
    Hist::New("hPROF2_HC_mass_ISS" , AXPROFmass);
    Hist::New("hPROF2_HC_mass_ISS2", AXPROFmass);
    Hist::New("hPROF2_HC_mass_ISS3", AXPROFmass);
    Hist::New("hPROF2_HC_mass_ISS4", AXPROFmass);
    Hist::New("hPROF2_HC_mass_ISS5", AXPROFmass);
    Hist::New("hPROF2_HC_mass_ISS6", AXPROFmass);
    Hist::New("hPROF2_HC_mass_map_ISS" , HistAxis(AXPROFmass, AXPROFmass));
    Hist::New("hPROF2_HC_mass_map_ISS2", HistAxis(AXPROFmass, AXPROFmass));
    Hist::New("hPROF2_HC_mass_map_ISS3", HistAxis(AXPROFmass, AXPROFmass));
    Hist::New("hPROF2_HC_mass_map_ISS4", HistAxis(AXPROFmass, AXPROFmass));
    Hist::New("hPROF2_HC_mass_MC_FLUX", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_FLUX2", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_FLUX3", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC2_FLUX", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC2_FLUX2", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC2_FLUX3", AXPROFmass);
    Hist::New("hPROF2_HC_mass_cnt_ISS", HistAxis(AXrig3));
    Hist::New("hPROF2_HC_mass_cnt_ISS2", HistAxis(AXrig3));
    Hist::New("hPROF2_HC_mass_cnt_ISS3", HistAxis(AXrig3));
    Hist::New("hPROF2_HC_mass_cnt_ISS4", HistAxis(AXrig3));
    Hist::New("hPROF2_HC_mass_cnt_MC", HistAxis(AXrig3));
    Hist::New("hPROF2_HC_mass_cnt_MC_FLUX", HistAxis(AXrig3));
    Hist::New("hPROF2_HC_mass_cnt_MC_FLUX2", HistAxis(AXrig3));
    Hist::New("hPROF2_HC_mass_cnt_MC_FLUX3", HistAxis(AXrig3));
    
    Hist::New("hPROF2_HC_mass_MC_PI_POS_FLUX" , AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_PI_POS_FLUX2", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_PI_POS_FLUX3", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_PI_NEG_FLUX" , AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_PI_NEG_FLUX2", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_PI_NEG_FLUX3", AXPROFmass);
    
    Hist::New("hPROF2_HC_mass_MC_K_POS_FLUX" , AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_K_POS_FLUX2", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_K_POS_FLUX3", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_K_NEG_FLUX" , AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_K_NEG_FLUX2", AXPROFmass);
    Hist::New("hPROF2_HC_mass_MC_K_NEG_FLUX3", AXPROFmass);
    
    Axis AXPROFrvar1("RVAR1", 1000, -3, 3);
    Axis AXPROFrvar("RVAR", 2000, -20, 20);
    Hist::New("hPROF_CK_rvar1", HistAxis(AXPROFrvar1));
    Hist::New("hPROF_CK_rvar2", HistAxis(AXrig, AXPROFrvar));
    
    Hist::New("hPROF_KF_rvar1", HistAxis(AXPROFrvar1));
    Hist::New("hPROF_KF_rvar2", HistAxis(AXrig, AXPROFrvar));
    
    Hist::New("hPROF_HC_rvar1", HistAxis(AXPROFrvar1));
    Hist::New("hPROF_HC_rvar2", HistAxis(AXrig, AXPROFrvar));
    
    Hist::New("hPROF2_CK_rvar1", HistAxis(AXPROFrvar1));
    Hist::New("hPROF2_CK_rvar2", HistAxis(AXrig, AXPROFrvar));
    
    Hist::New("hPROF2_KF_rvar1", HistAxis(AXPROFrvar1));
    Hist::New("hPROF2_KF_rvar2", HistAxis(AXrig, AXPROFrvar));
    
    Hist::New("hPROF2_HC_rvar1", HistAxis(AXPROFrvar1));
    Hist::New("hPROF2_HC_rvar2", HistAxis(AXrig, AXPROFrvar));
    
    Axis AXPROFbta("beta", 800, 0.92, 1.08);
    Hist::New("hPROF_bta",  HistAxis(AXtme, AXPROFbta));
    Hist::New("hPROF_bta1", HistAxis(AXtme, AXPROFbta));
    Hist::New("hPROF_bta2", HistAxis(AXtme, AXPROFbta));
    Hist::New("hPROF_bta3", HistAxis(AXtme, AXPROFbta));
    Hist::New("hPROF_bta4", HistAxis(AXtme, AXPROFbta));
    
    Axis AXPROFtile("tile", 122, 0, 122);
    Axis AXPROFrbta("beta", 800, 0.985, 1.015);
    Hist::New("hPROF_tile_rbta", HistAxis(AXPROFtile, AXPROFrbta));
    Hist::New("hPROF_rbta", HistAxis(AXtme, AXPROFrbta));
    Hist::New("hPROF_rbta25", HistAxis(AXtme, AXPROFrbta));
    Hist::New("hPROF_rbta26", HistAxis(AXtme, AXPROFrbta));
    Hist::New("hPROF_rbta27", HistAxis(AXtme, AXPROFrbta));
    Hist::New("hPROF_rbta28", HistAxis(AXtme, AXPROFrbta));
    
    Hist::New("h_MC_cnt", HistAxis(AXrig));
    Hist::New("hT_MC_cnt", HistAxis(AXrig2));
    
    Axis AXLextseg("TRD extseg", 20, 0.0, 20);
    Hist::New("hLP_extseg", HistAxis(AXrig, AXLextseg));
    Hist::New("hLN_extseg", HistAxis(AXrig, AXLextseg));
    
    Axis AXLexthit("TRD exthit", 200, 0.0, 200);
    Hist::New("hLP_exthit", HistAxis(AXrig, AXLexthit));
    Hist::New("hLN_exthit", HistAxis(AXrig, AXLexthit));
    
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
    
    Hist::New("hIP_cnt", HistAxis(AXrig));
    Hist::New("hIN_cnt", HistAxis(AXrig));
    
    Hist::New("hI_MC_cnt", HistAxis(AXrig));
    Hist::New("hTI_MC_cnt", HistAxis(AXrig2));
    
    Axis AXHllr("TRD estimator", 100, 0.0, 1.6);
    Axis AXHrvar("RVAR", 400, -10, 10);
    Axis AXHlchi("Log(#chi^{2}/DOF)", 100, -4, 8);
    
    Axis AXHrvar2("|RVAR|", 100, 0, 10);
    Axis AXHlchi2("Log(#chi^{2}/DOF)", 100, -3, 3);

    Hist::New("hHPl1_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNl1_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHPl1_rvar", HistAxis(AXrig, AXHrvar));
    Hist::New("hHNl1_rvar", HistAxis(AXrig, AXHrvar));
    Hist::New("hHPl1_lchix_in", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl1_lchix_in", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPl1_lchiy_in", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl1_lchiy_in", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPl1_lchix", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl1_lchix", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPl1_lchiy", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl1_lchiy", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPl1_cnt", HistAxis(AXrig));
    Hist::New("hHNl1_cnt", HistAxis(AXrig));
    
    Hist::New("hHl1_MC_cnt", HistAxis(AXrig));
    
    Hist::New("hHPl1_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNl1_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Hist::New("hHPl1_MC_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNl1_MC_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Hist::New("hHPl1_lrvar_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHNl1_lrvar_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHPl1_lchiy_cut", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl1_lchiy_cut", HistAxis(AXrig, AXHlchi));
    
    Hist::New("hHPl9_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNl9_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHPl9_rvar", HistAxis(AXrig, AXHrvar));
    Hist::New("hHNl9_rvar", HistAxis(AXrig, AXHrvar));
    Hist::New("hHPl9_lchix_in", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl9_lchix_in", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPl9_lchiy_in", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl9_lchiy_in", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPl9_lchix", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl9_lchix", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPl9_lchiy", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl9_lchiy", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPl9_cnt", HistAxis(AXrig));
    Hist::New("hHNl9_cnt", HistAxis(AXrig));
    
    Hist::New("hHl9_MC_cnt", HistAxis(AXrig));
    
    Hist::New("hHPl9_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNl9_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Hist::New("hHPl9_MC_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNl9_MC_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Hist::New("hHPl9_lrvar_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHNl9_lrvar_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHPl9_lchiy_cut", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNl9_lchiy_cut", HistAxis(AXrig, AXHlchi));
    
    Hist::New("hHPfs_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNfs_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHPfs_rvar_l1", HistAxis(AXrig, AXHrvar));
    Hist::New("hHNfs_rvar_l1", HistAxis(AXrig, AXHrvar));
    Hist::New("hHPfs_rvar_l9", HistAxis(AXrig, AXHrvar));
    Hist::New("hHNfs_rvar_l9", HistAxis(AXrig, AXHrvar));
    Hist::New("hHPfs_rvar", HistAxis(AXrig, AXHrvar));
    Hist::New("hHNfs_rvar", HistAxis(AXrig, AXHrvar));
    Hist::New("hHPfs_lchix_l1", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchix_l1", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPfs_lchiy_l1", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchiy_l1", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPfs_lchix_l9", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchix_l9", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPfs_lchiy_l9", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchiy_l9", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPfs_lchix", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchix", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPfs_lchiy", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchiy", HistAxis(AXrig, AXHlchi));
    Hist::New("hHPfs_cnt", HistAxis(AXrig));
    Hist::New("hHNfs_cnt", HistAxis(AXrig));
    
    Hist::New("hHfs_MC_cnt", HistAxis(AXrig));
    
    Hist::New("hHPfs_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNfs_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Hist::New("hHPfs_lchiy_lrvar_v1", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNfs_lchiy_lrvar_v1", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Hist::New("hHPfs_lchiy_lrvar_v2", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNfs_lchiy_lrvar_v2", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Hist::New("hHPfs_lchiy_lrvar_v3", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNfs_lchiy_lrvar_v3", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Axis AXHrvar3("|RVAR|", 75, 0, 8);
    Axis AXHlchi3("Log(#chi^{2}/DOF)", 75, -3, 3);
    Hist::New("hHPfs_lchiy_lrvar2", HistAxis(AXrig, AXHlchi3, AXHrvar3));
    Hist::New("hHNfs_lchiy_lrvar2", HistAxis(AXrig, AXHlchi3, AXHrvar3));

    Hist::New("hHPfs_lchiy_lrvar2_v1", HistAxis(AXrig, AXHlchi3, AXHrvar3));
    Hist::New("hHNfs_lchiy_lrvar2_v1", HistAxis(AXrig, AXHlchi3, AXHrvar3));
    
    Hist::New("hHPfs_lchiy_lrvar2_v2", HistAxis(AXrig, AXHlchi3, AXHrvar3));
    Hist::New("hHNfs_lchiy_lrvar2_v2", HistAxis(AXrig, AXHlchi3, AXHrvar3));
    
    Hist::New("hHPfs_lchiy_lrvar2_v3", HistAxis(AXrig, AXHlchi3, AXHrvar3));
    Hist::New("hHNfs_lchiy_lrvar2_v3", HistAxis(AXrig, AXHlchi3, AXHrvar3));
    
    Axis AXHrvar4("|RVAR|", 50, 0, 8);
    Axis AXHlchi4("Log(#chi^{2}/DOF)", 50, -3, 3);
    Hist::New("hHPfs_lchiy_lrvar3", HistAxis(AXrig, AXHlchi4, AXHrvar4));
    Hist::New("hHNfs_lchiy_lrvar3", HistAxis(AXrig, AXHlchi4, AXHrvar4));
    
    Hist::New("hHPfs_lchiy_lrvar3_v1", HistAxis(AXrig, AXHlchi4, AXHrvar4));
    Hist::New("hHNfs_lchiy_lrvar3_v1", HistAxis(AXrig, AXHlchi4, AXHrvar4));
    
    Hist::New("hHPfs_lchiy_lrvar3_v2", HistAxis(AXrig, AXHlchi4, AXHrvar4));
    Hist::New("hHNfs_lchiy_lrvar3_v2", HistAxis(AXrig, AXHlchi4, AXHrvar4));
    
    Hist::New("hHPfs_lchiy_lrvar3_v3", HistAxis(AXrig, AXHlchi4, AXHrvar4));
    Hist::New("hHNfs_lchiy_lrvar3_v3", HistAxis(AXrig, AXHlchi4, AXHrvar4));
    
    Hist::New("hHPfs_MC_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    Hist::New("hHNfs_MC_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHrvar2));
    
    Hist::New("hHPfs_lrvarl1_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHNfs_lrvarl1_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHPfs_lrvarl9_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHNfs_lrvarl9_cut", HistAxis(AXrig, AXHrvar2));
    
    Hist::New("hHPfs_lrvar_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHNfs_lrvar_cut", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHPfs_lchiy_cut", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchiy_cut", HistAxis(AXrig, AXHlchi));

    Hist::New("hHPfs_lrvar_cut1", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHNfs_lrvar_cut1", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHPfs_lchiy_cut1", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchiy_cut1", HistAxis(AXrig, AXHlchi));

    Hist::New("hHPfs_lrvar_cut2", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHNfs_lrvar_cut2", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHPfs_lchiy_cut2", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchiy_cut2", HistAxis(AXrig, AXHlchi));

    Hist::New("hHPfs_lrvar_cut3", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHNfs_lrvar_cut3", HistAxis(AXrig, AXHrvar2));
    Hist::New("hHPfs_lchiy_cut3", HistAxis(AXrig, AXHlchi));
    Hist::New("hHNfs_lchiy_cut3", HistAxis(AXrig, AXHlchi));

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
    if ((trg->bit&8) != 8) return false;

    // Charge
    if (tof->Qall < 0.8 || tof->Qall > 1.3) return false;
    if (trk->QIn < 0.8 || trk->QIn > 1.3) return false;

    // ACC
    if (acc->num_cls != 0) return false;

    // TOF
    if (tof->num_beta != 1) return false;
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
            
    // TRD
    //if (trd->num_track == 0) return false;
    //if (trd->num_track != 1) return false;
    //if (!trd->status) return false;
    //if (trd->status && trd->num_extra_hit > 12) return false;
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8) return false;
    
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
    double mc_wgt_flux = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    
    double mc_mass_wgt_flux  = (g4mc != nullptr) ? std::pow(g4mc->prm_mom, -1.7) : 1.0;
    double mc_mass_wgt_flux2 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom, -1.7) : 1.0;
    double mc_mass_wgt_flux3 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom, -1.7) : 1.0;
    if (g4mc != nullptr) {
        double crr = (*Hist::Head("hPROFwgt"))()->Interpolate(g4mc->prm_mom);
        if (crr < 0) crr = 0.0;
        mc_mass_wgt_flux2 *= crr;
        
        double apcrr = (*Hist::Head("hPROFapp"))()->Interpolate(g4mc->prm_mom);
        if (apcrr < 0) apcrr = 0.0;
        mc_mass_wgt_flux3 *= (crr * apcrr);
    }
    
    int ntkL1 = trk->ext_num_hit[0];
    int ntkL9 = trk->ext_num_hit[8];
    int ntkL2 = trk->ext_num_hit[1];
    int ntkIn = (trk->ext_num_hit[2] + trk->ext_num_hit[3] + trk->ext_num_hit[4] + trk->ext_num_hit[5] + trk->ext_num_hit[6] + trk->ext_num_hit[7]);
    bool is_clean_tk = (ntkL1 == 0 && ntkL2 == 0 && ntkIn == 0);
    
    if (trk->lay[1] != 0 &&
        trk->ck_status[0] && 
        std::log(trk->ck_nchi[0][0]) < 1.6 &&
        std::log(trk->ck_nchi[0][1]) < 1.6 &&
        tof->status &&
        std::log(tof->nchi_t) < 1.75 &&
        std::log(tof->nchi_c) < 1.75 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        trd->status
        ) {
        double sign = (trk->ck_rig[0] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * std::sqrt((trk->ck_rig[0] * (1.0 / tof->beta + 1.0)) * (trk->ck_rig[0] * (1.0 / tof->beta - 1.0)));
        if (tof->beta > 0.5 && tof->beta < 0.8) {
            Hist::Head("hPROF_CK_mass_ISS")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.70 : true) Hist::Head("hPROF_CK_mass_ISS2")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.75 : true) Hist::Head("hPROF_CK_mass_ISS3")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.80 : true) Hist::Head("hPROF_CK_mass_ISS4")->fillH1D(mass, list->weight);
            Hist::Head("hPROF_CK_mass_MC_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
            Hist::Head("hPROF_CK_mass_MC_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
            Hist::Head("hPROF_CK_mass_MC_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_CK_mass_MC2_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_CK_mass_MC2_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_CK_mass_MC2_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
        }
        if (trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_ISS")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
        if (trk->ck_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.70) Hist::Head("hPROF_CK_mass_cnt_ISS2")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
        if (trk->ck_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.75) Hist::Head("hPROF_CK_mass_cnt_ISS3")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
        if (trk->ck_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.80) Hist::Head("hPROF_CK_mass_cnt_ISS4")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
        if (g4mc != nullptr && trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
        if (g4mc != nullptr && trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC_FLUX")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux * list->weight);
        if (g4mc != nullptr && trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC_FLUX2")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux2 * list->weight);
        if (g4mc != nullptr && trk->ck_rig[0] > 0) Hist::Head("hPROF_CK_mass_cnt_MC_FLUX3")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux3 * list->weight);
        if (trd->num_vtx[0] == 0 && trd->num_vtx[1] == 0 && is_clean_tk) {
            if (tof->beta > 0.5 && tof->beta < 0.8) {
                Hist::Head("hPROF2_CK_mass_ISS")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.70 : true) Hist::Head("hPROF2_CK_mass_ISS2")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.75 : true) Hist::Head("hPROF2_CK_mass_ISS3")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.80 : true) Hist::Head("hPROF2_CK_mass_ISS4")->fillH1D(mass, list->weight);
                Hist::Head("hPROF2_CK_mass_MC_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
                Hist::Head("hPROF2_CK_mass_MC_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
                Hist::Head("hPROF2_CK_mass_MC_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_CK_mass_MC2_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_CK_mass_MC2_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_CK_mass_MC2_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
            }
            if (trk->ck_rig[0] > 0) Hist::Head("hPROF2_CK_mass_cnt_ISS")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
            if (trk->ck_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.70) Hist::Head("hPROF2_CK_mass_cnt_ISS2")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
            if (trk->ck_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.75) Hist::Head("hPROF2_CK_mass_cnt_ISS3")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
            if (trk->ck_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->ck_rig[0])/rti->max_IGRF)>0.80) Hist::Head("hPROF2_CK_mass_cnt_ISS4")->fillH1D(std::abs(trk->ck_rig[0]), list->weight);
            if (g4mc != nullptr && trk->ck_rig[0] > 0) Hist::Head("hPROF2_CK_mass_cnt_MC")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
            if (g4mc != nullptr && trk->ck_rig[0] > 0) Hist::Head("hPROF2_CK_mass_cnt_MC_FLUX")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux * list->weight);
            if (g4mc != nullptr && trk->ck_rig[0] > 0) Hist::Head("hPROF2_CK_mass_cnt_MC_FLUX2")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux2 * list->weight);
            if (g4mc != nullptr && trk->ck_rig[0] > 0) Hist::Head("hPROF2_CK_mass_cnt_MC_FLUX3")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux3 * list->weight);
        }
    }
    
    if (trk->lay[1] != 0 &&
        trk->kf_status[0] && 
        std::log(trk->kf_nchi[0][0]) < 1.6 &&
        std::log(trk->kf_nchi[0][1]) < 1.6 &&
        tof->status &&
        std::log(tof->nchi_t) < 1.75 &&
        std::log(tof->nchi_c) < 1.75 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        trd->status
        ) {
        double sign = (trk->kf_cen_rig[0] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * std::sqrt((trk->kf_cen_rig[0] * (1.0 / tof->beta + 1.0)) * (trk->kf_cen_rig[0] * (1.0 / tof->beta - 1.0)));
        if (tof->beta > 0.5 && tof->beta < 0.8) {
            Hist::Head("hPROF_KF_mass_ISS")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.70 : true) Hist::Head("hPROF_KF_mass_ISS2")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.75 : true) Hist::Head("hPROF_KF_mass_ISS3")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.80 : true) Hist::Head("hPROF_KF_mass_ISS4")->fillH1D(mass, list->weight);
            Hist::Head("hPROF_KF_mass_MC_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
            Hist::Head("hPROF_KF_mass_MC_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
            Hist::Head("hPROF_KF_mass_MC_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_KF_mass_MC2_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_KF_mass_MC2_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_KF_mass_MC2_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
        }
        if (trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_ISS")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
        if (trk->kf_top_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.70) Hist::Head("hPROF_KF_mass_cnt_ISS2")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
        if (trk->kf_top_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.75) Hist::Head("hPROF_KF_mass_cnt_ISS3")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
        if (trk->kf_top_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.80) Hist::Head("hPROF_KF_mass_cnt_ISS4")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
        if (g4mc != nullptr && trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
        if (g4mc != nullptr && trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC_FLUX")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux * list->weight);
        if (g4mc != nullptr && trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC_FLUX2")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux2 * list->weight);
        if (g4mc != nullptr && trk->kf_top_rig[0] > 0) Hist::Head("hPROF_KF_mass_cnt_MC_FLUX3")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux3 * list->weight);
        if (trd->num_vtx[0] == 0 && trd->num_vtx[1] == 0 && is_clean_tk) {
            if (tof->beta > 0.5 && tof->beta < 0.8) {
                Hist::Head("hPROF2_KF_mass_ISS")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.70 : true) Hist::Head("hPROF2_KF_mass_ISS2")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.75 : true) Hist::Head("hPROF2_KF_mass_ISS3")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.80 : true) Hist::Head("hPROF2_KF_mass_ISS4")->fillH1D(mass, list->weight);
                Hist::Head("hPROF2_KF_mass_MC_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
                Hist::Head("hPROF2_KF_mass_MC_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
                Hist::Head("hPROF2_KF_mass_MC_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_KF_mass_MC2_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_KF_mass_MC2_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_KF_mass_MC2_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
            }
            if (trk->kf_top_rig[0] > 0) Hist::Head("hPROF2_KF_mass_cnt_ISS")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
            if (trk->kf_top_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.70) Hist::Head("hPROF2_KF_mass_cnt_ISS2")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
            if (trk->kf_top_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.75) Hist::Head("hPROF2_KF_mass_cnt_ISS3")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
            if (trk->kf_top_rig[0] > 0 && CheckType(Type::ISS) && (std::abs(trk->kf_top_rig[0])/rti->max_IGRF)>0.80) Hist::Head("hPROF2_KF_mass_cnt_ISS4")->fillH1D(std::abs(trk->kf_top_rig[0]), list->weight);
            if (g4mc != nullptr && trk->kf_top_rig[0] > 0) Hist::Head("hPROF2_KF_mass_cnt_MC")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
            if (g4mc != nullptr && trk->kf_top_rig[0] > 0) Hist::Head("hPROF2_KF_mass_cnt_MC_FLUX")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux * list->weight);
            if (g4mc != nullptr && trk->kf_top_rig[0] > 0) Hist::Head("hPROF2_KF_mass_cnt_MC_FLUX2")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux2 * list->weight);
            if (g4mc != nullptr && trk->kf_top_rig[0] > 0) Hist::Head("hPROF2_KF_mass_cnt_MC_FLUX3")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux3 * list->weight);
        }
    }
    
    if (trk->lay[1] != 0 &&
        hyc->geom_status[0] &&
        std::log(hyc->geom_nchi_x[0]) < 1.75 &&
        std::log(hyc->geom_nchi_y[0]) < 1.75 &&
        hyc->vel_status[1] &&
        std::log(hyc->vel_nchi[1]) < 2.00 &&
        hyc->mutr_status[1] && 
        std::log(hyc->mutr_nchi_x[1]) < 1.00 &&
        std::log(hyc->mutr_nchi_y[1]) < 1.00 &&
        std::log(hyc->mutr_nchi_b[1]) < 1.25 &&
        tof->extcls_noise == 0 &&
        tof->num_in_time_cls <= 4 &&
        trd->status
        ) {
        double sign = (hyc->mutr_top_rig[1] >= 0.0) ? 1.0 : -1.0;
        double mass = sign * hyc->mutr_mass[1];
        
        double ck_sign = (trk->ck_rig[0] >= 0.0) ? 1.0 : -1.0;
        double ck_mass = ck_sign * std::sqrt((trk->ck_rig[0] * (1.0 / tof->beta + 1.0)) * (trk->ck_rig[0] * (1.0 / tof->beta - 1.0)));
        
        if (hyc->mutr_top_bta[1] > 0.5 && hyc->mutr_top_bta[1] < 0.8) {
            //int ntkL1 = trk->ext_num_hit[0];
            //int ntkL9 = trk->ext_num_hit[8];
            //int ntkL2 = trk->ext_num_hit[1];
            //int ntkIn = (trk->ext_num_hit[2] + trk->ext_num_hit[3] + trk->ext_num_hit[4] + trk->ext_num_hit[5] + trk->ext_num_hit[6] + trk->ext_num_hit[7]);
            //Hist::Head("hPROF_HC_mass_ISS_ntkL1")->fillH2D(mass, ntkL1, list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_ntkL9")->fillH2D(mass, ntkL9, list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_ntkL2")->fillH2D(mass, ntkL2, list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_ntkIn")->fillH2D(mass, ntkIn, list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_extvtxx")->fillH2D(mass, trd->num_vtx[0], list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_extvtxy")->fillH2D(mass, trd->num_vtx[1], list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_extseg")->fillH2D(mass, trd->num_extra_seg, list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_exthit")->fillH2D(mass, trd->num_extra_hit, list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_nchix")->fillH2D(mass, std::log(hyc->geom_nchi_x[0]), list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_nchiy")->fillH2D(mass, std::log(hyc->geom_nchi_y[0]), list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_nchib")->fillH2D(mass, std::log(hyc->vel_nchi[1]), list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_nchix2")->fillH2D(mass, std::log(hyc->mutr_nchi_x[1]), list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_nchiy2")->fillH2D(mass, std::log(hyc->mutr_nchi_y[1]), list->weight);
            //Hist::Head("hPROF_HC_mass_ISS_nchib2")->fillH2D(mass, std::log(hyc->mutr_nchi_b[1]), list->weight);

            Hist::Head("hPROF_HC_mass_ISS")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.70 : true) Hist::Head("hPROF_HC_mass_ISS2")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.75 : true) Hist::Head("hPROF_HC_mass_ISS3")->fillH1D(mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.80 : true) Hist::Head("hPROF_HC_mass_ISS4")->fillH1D(mass, list->weight);
            
            Hist::Head("hPROF_HC_mass_map_ISS")->fillH2D(mass, ck_mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.70 : true) Hist::Head("hPROF_HC_mass_map_ISS2")->fillH2D(mass, ck_mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.75 : true) Hist::Head("hPROF_HC_mass_map_ISS3")->fillH2D(mass, ck_mass, list->weight);
            if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.80 : true) Hist::Head("hPROF_HC_mass_map_ISS4")->fillH2D(mass, ck_mass, list->weight);

            Hist::Head("hPROF_HC_mass_MC_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
            Hist::Head("hPROF_HC_mass_MC_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
            Hist::Head("hPROF_HC_mass_MC_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_HC_mass_MC2_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_HC_mass_MC2_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
            if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF_HC_mass_MC2_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
        }

        if (hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF_HC_mass_cnt_ISS")->fillH1D(std::abs(hyc->phys_top_rig[0][0]), list->weight);
        if (hyc->phys_top_rig[0][0] > 0 && CheckType(Type::ISS) && (std::abs(hyc->phys_top_rig[0][0])/rti->max_IGRF)>0.70) Hist::Head("hPROF_HC_mass_cnt_ISS2")->fillH1D(std::abs(hyc->phys_top_rig[0][0]), list->weight);
        if (hyc->phys_top_rig[0][0] > 0 && CheckType(Type::ISS) && (std::abs(hyc->phys_top_rig[0][0])/rti->max_IGRF)>0.75) Hist::Head("hPROF_HC_mass_cnt_ISS3")->fillH1D(std::abs(hyc->phys_top_rig[0][0]), list->weight);
        if (hyc->phys_top_rig[0][0] > 0 && CheckType(Type::ISS) && (std::abs(hyc->phys_top_rig[0][0])/rti->max_IGRF)>0.80) Hist::Head("hPROF_HC_mass_cnt_ISS4")->fillH1D(std::abs(hyc->phys_top_rig[0][0]), list->weight);
        if (g4mc != nullptr && hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF_HC_mass_cnt_MC")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
        if (g4mc != nullptr && hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF_HC_mass_cnt_MC_FLUX")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux * list->weight);
        if (g4mc != nullptr && hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF_HC_mass_cnt_MC_FLUX2")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux2 * list->weight);
        if (g4mc != nullptr && hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF_HC_mass_cnt_MC_FLUX3")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux3 * list->weight);
        if (trd->num_vtx[0] == 0 && trd->num_vtx[1] == 0 && is_clean_tk) {
            if (hyc->mutr_top_bta[1] > 0.5 && hyc->mutr_top_bta[1] < 0.8) {
                Hist::Head("hPROF2_HC_mass_ISS")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.70 : true) Hist::Head("hPROF2_HC_mass_ISS2")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.75 : true) Hist::Head("hPROF2_HC_mass_ISS3")->fillH1D(mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.80 : true) Hist::Head("hPROF2_HC_mass_ISS4")->fillH1D(mass, list->weight);
                
                if (mass < -2.2) clone_tree->Fill();
                
                Hist::Head("hPROF2_HC_mass_map_ISS")->fillH2D(mass, ck_mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.70 : true) Hist::Head("hPROF2_HC_mass_map_ISS2")->fillH2D(mass, ck_mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.75 : true) Hist::Head("hPROF2_HC_mass_map_ISS3")->fillH2D(mass, ck_mass, list->weight);
                if (CheckType(Type::ISS) ? (std::abs(hyc->mutr_top_rig[1])/rti->max_IGRF)>0.80 : true) Hist::Head("hPROF2_HC_mass_map_ISS4")->fillH2D(mass, ck_mass, list->weight);
                
                Hist::Head("hPROF2_HC_mass_MC_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
                Hist::Head("hPROF2_HC_mass_MC_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
                Hist::Head("hPROF2_HC_mass_MC_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_HC_mass_MC2_FLUX")->fillH1D(mass, mc_mass_wgt_flux * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_HC_mass_MC2_FLUX2")->fillH1D(mass, mc_mass_wgt_flux2 * list->weight);
                if (CheckType(Type::MC) && (g4mc->tk[6] && g4mc->tk[7]) && (g4mc->tf[2] || g4mc->tf[3])) Hist::Head("hPROF2_HC_mass_MC2_FLUX3")->fillH1D(mass, mc_mass_wgt_flux3 * list->weight);

                // pion 0 ~ 0.3
                if (std::abs(mass) < 0.28) {
                    Hist::Head("hPROF2_HC_mass_MC_PI_POS_FLUX" )->fillH1D(std::abs(mass), mc_mass_wgt_flux * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_PI_POS_FLUX2")->fillH1D(std::abs(mass), mc_mass_wgt_flux2 * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_PI_POS_FLUX3")->fillH1D(std::abs(mass), mc_mass_wgt_flux3 * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_PI_NEG_FLUX" )->fillH1D(-std::abs(mass), mc_mass_wgt_flux * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_PI_NEG_FLUX2")->fillH1D(-std::abs(mass), mc_mass_wgt_flux2 * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_PI_NEG_FLUX3")->fillH1D(-std::abs(mass), mc_mass_wgt_flux3 * list->weight);
                }
                
                // koan -0.35 ~ -0.75
                if (mass > -0.75 && mass < -0.35) {
                    Hist::Head("hPROF2_HC_mass_MC_K_POS_FLUX" )->fillH1D(std::abs(mass), mc_mass_wgt_flux * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_K_POS_FLUX2")->fillH1D(std::abs(mass), mc_mass_wgt_flux2 * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_K_POS_FLUX3")->fillH1D(std::abs(mass), mc_mass_wgt_flux3 * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_K_NEG_FLUX" )->fillH1D(-std::abs(mass), mc_mass_wgt_flux * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_K_NEG_FLUX2")->fillH1D(-std::abs(mass), mc_mass_wgt_flux2 * list->weight);
                    Hist::Head("hPROF2_HC_mass_MC_K_NEG_FLUX3")->fillH1D(-std::abs(mass), mc_mass_wgt_flux3 * list->weight);
                }
            }
            if (hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF2_HC_mass_cnt_ISS")->fillH1D(std::abs(hyc->phys_top_rig[0][0]), list->weight);
            if (hyc->phys_top_rig[0][0] > 0 && CheckType(Type::ISS) && (std::abs(hyc->phys_top_rig[0][0])/rti->max_IGRF)>0.70) Hist::Head("hPROF2_HC_mass_cnt_ISS2")->fillH1D(std::abs(hyc->phys_top_rig[0][0]), list->weight);
            if (hyc->phys_top_rig[0][0] > 0 && CheckType(Type::ISS) && (std::abs(hyc->phys_top_rig[0][0])/rti->max_IGRF)>0.75) Hist::Head("hPROF2_HC_mass_cnt_ISS3")->fillH1D(std::abs(hyc->phys_top_rig[0][0]), list->weight);
            if (hyc->phys_top_rig[0][0] > 0 && CheckType(Type::ISS) && (std::abs(hyc->phys_top_rig[0][0])/rti->max_IGRF)>0.80) Hist::Head("hPROF2_HC_mass_cnt_ISS4")->fillH1D(std::abs(hyc->phys_top_rig[0][0]), list->weight);
            if (g4mc != nullptr && hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF2_HC_mass_cnt_MC")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
            if (g4mc != nullptr && hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF2_HC_mass_cnt_MC_FLUX")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux * list->weight);
            if (g4mc != nullptr && hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF2_HC_mass_cnt_MC_FLUX2")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux2 * list->weight);
            if (g4mc != nullptr && hyc->phys_top_rig[0][0] > 0) Hist::Head("hPROF2_HC_mass_cnt_MC_FLUX3")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_mass_wgt_flux3 * list->weight);
        }
    }

    if (ecal->status && ecal->mvaBDT > -0.6) return false;
    const std::array<double, 2> frso_pr_l1({ 2.17407e-01, 7.29381e+02 });
    const std::array<double, 2> frso_pr_l9({ 1.49239e-01, 9.52511e+02 });
    const std::array<double, 2> frso_pr_fs({ 1.47495e-01, 1.92807e+03 });

    if (trk->ck_status[0] && 
        trk->ck_status[1] && 
        trk->ck_status[2] && 
        trk->ck_status[3]) {
        //double mrig = std::sqrt(std::abs(trk->ck_rig[2] * trk->ck_rig[1]));
        double mrig = 100.;
        double drig = mrig * (1.0/trk->ck_rig[2] - 1.0/trk->ck_rig[1]);
        //double rsl1 = std::sqrt(frso_pr_l1[0]*frso_pr_l1[0] + (mrig/frso_pr_l1[1])*(mrig/frso_pr_l1[1]));
        //double rsl9 = std::sqrt(frso_pr_l9[0]*frso_pr_l9[0] + (mrig/frso_pr_l9[1])*(mrig/frso_pr_l9[1]));
        //double rs19 = std::sqrt(rsl1 * rsl1 + rsl9 * rsl9);
        //double rvar = drig / rs19;
        double rvar = (trk->ck_rig[1]/trk->ck_rig[2] - 1.0);
        if (std::fabs(trk->ck_rig[3]) > mrig) Hist::Head("hPROF_CK_rvar1")->fillH1D(drig, list->weight);
        if (std::fabs(trk->ck_rig[3]) > 0) Hist::Head("hPROF_CK_rvar2")->fillH2D(trk->ck_rig[3], rvar, list->weight);
        if (std::log(trk->ck_nchi[1][0]) < 1.75 &&
            std::log(trk->ck_nchi[1][1]) < 1.75 &&
            std::log(trk->ck_nchi[2][0]) < 1.75 &&
            std::log(trk->ck_nchi[2][1]) < 1.75) {
            if (std::fabs(trk->ck_rig[3]) > mrig) Hist::Head("hPROF2_CK_rvar1")->fillH1D(drig, list->weight);
            if (std::fabs(trk->ck_rig[3]) > 0) Hist::Head("hPROF2_CK_rvar2")->fillH2D(trk->ck_rig[3], rvar, list->weight);
        }
    }
    if (trk->kf_status[0] && 
        trk->kf_status[1] && 
        trk->kf_status[2] && 
        trk->kf_status[3]) {
        //double mrig = std::sqrt(std::abs(trk->kf_top_rig[2] * trk->kf_top_rig[1]));
        double mrig = 100.;
        double drig = mrig * (1.0/trk->kf_top_rig[2] - 1.0/trk->kf_top_rig[1]);
        //double rsl1 = std::sqrt(frso_pr_l1[0]*frso_pr_l1[0] + (mrig/frso_pr_l1[1])*(mrig/frso_pr_l1[1]));
        //double rsl9 = std::sqrt(frso_pr_l9[0]*frso_pr_l9[0] + (mrig/frso_pr_l9[1])*(mrig/frso_pr_l9[1]));
        //double rs19 = std::sqrt(rsl1 * rsl1 + rsl9 * rsl9);
        //double rvar = drig / rs19;
        double rvar = (trk->kf_top_rig[1]/trk->kf_top_rig[2] - 1.0);
        if (std::fabs(trk->kf_top_rig[3]) > mrig) Hist::Head("hPROF_KF_rvar1")->fillH1D(drig, list->weight);
        if (std::fabs(trk->kf_top_rig[3]) > 0) Hist::Head("hPROF_KF_rvar2")->fillH2D(trk->kf_top_rig[3], rvar, list->weight);
        
        if (std::log(trk->kf_nchi[1][0]) < 1.75 &&
            std::log(trk->kf_nchi[1][1]) < 1.75 &&
            std::log(trk->kf_nchi[2][0]) < 1.75 &&
            std::log(trk->kf_nchi[2][1]) < 1.75) {
            if (std::fabs(trk->kf_top_rig[3]) > mrig) Hist::Head("hPROF2_KF_rvar1")->fillH1D(drig, list->weight);
            if (std::fabs(trk->kf_top_rig[3]) > 0) Hist::Head("hPROF2_KF_rvar2")->fillH2D(trk->kf_top_rig[3], rvar, list->weight);
        }
    }
    if (hyc->geom_status[0] &&
        hyc->geom_status[1] && 
        hyc->geom_status[2] && 
        hyc->geom_status[3]) {
        //double mrig = std::sqrt(std::abs(hyc->geom_top_rig[2] * hyc->geom_top_rig[1]));
        double mrig = 100.;
        double drig = mrig * (1.0/hyc->geom_top_rig[2] - 1.0/hyc->geom_top_rig[1]);
        //double rsl1 = std::sqrt(frso_pr_l1[0]*frso_pr_l1[0] + (mrig/frso_pr_l1[1])*(mrig/frso_pr_l1[1]));
        //double rsl9 = std::sqrt(frso_pr_l9[0]*frso_pr_l9[0] + (mrig/frso_pr_l9[1])*(mrig/frso_pr_l9[1]));
        //double rs19 = std::sqrt(rsl1 * rsl1 + rsl9 * rsl9);
        //double rvar = drig / rs19;
        double rvar = (hyc->geom_top_rig[1]/hyc->geom_top_rig[2] - 1.0);
        if (std::fabs(hyc->geom_top_rig[3]) > mrig) Hist::Head("hPROF_HC_rvar1")->fillH1D(drig, list->weight);
        if (hyc->geom_top_rig[3] > 0) Hist::Head("hPROF_HC_rvar2")->fillH2D(hyc->geom_top_rig[3], rvar, list->weight);
        if (std::log(hyc->geom_nchi_x[1]) < 1.7 &&
            std::log(hyc->geom_nchi_y[1]) < 1.7 &&
            std::log(hyc->geom_nchi_x[2]) < 1.7 &&
            std::log(hyc->geom_nchi_y[2]) < 1.7) {
            if (std::fabs(hyc->geom_top_rig[3]) > mrig) Hist::Head("hPROF2_HC_rvar1")->fillH1D(drig, list->weight);
            if (hyc->geom_top_rig[3] > 0) Hist::Head("hPROF2_HC_rvar2")->fillH2D(hyc->geom_top_rig[3], rvar, list->weight);
        }
    }
    
    if (hyc->geom_status[0] && hyc->geom_top_rig[0] > 30.0 && hyc->vel_status[0]) {
        double bta = hyc->vel_top_bta[0];
        Hist::Head("hPROF_bta")->fillH2D(list->utime, bta, list->weight);
        if (std::abs(tof->loc[0][0]) > 40) Hist::Head("hPROF_bta1")->fillH2D(list->utime, bta, list->weight);
        if (std::abs(tof->loc[1][0]) > 40) Hist::Head("hPROF_bta2")->fillH2D(list->utime, bta, list->weight);
        if (std::abs(tof->loc[2][0]) > 40) Hist::Head("hPROF_bta3")->fillH2D(list->utime, bta, list->weight);
        if (std::abs(tof->loc[3][0]) > 40) Hist::Head("hPROF_bta4")->fillH2D(list->utime, bta, list->weight);
    }
    
    if (hyc->geom_status[0] && hyc->geom_top_rig[0] > 30.0 && hyc->vel_status[2] && rich->self_kind == 1) {
        double bta = hyc->vel_top_bta[2];
        Hist::Head("hPROF_tile_rbta")->fillH2D(rich->self_tile, bta, list->weight);
        Hist::Head("hPROF_rbta")->fillH2D(list->utime, bta, list->weight);
        if (rich->self_tile == 25) Hist::Head("hPROF_rbta25")->fillH2D(list->utime, bta, list->weight);
        if (rich->self_tile == 26) Hist::Head("hPROF_rbta26")->fillH2D(list->utime, bta, list->weight);
        if (rich->self_tile == 27) Hist::Head("hPROF_rbta27")->fillH2D(list->utime, bta, list->weight);
        if (rich->self_tile == 28) Hist::Head("hPROF_rbta28")->fillH2D(list->utime, bta, list->weight);
    }
    
    return true;
}

bool Analyzer::process_data_l() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    // TRD
    if (!trd->status) return false;
    //if (trd->status && trd->num_extra_hit > 8) return false;

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
    if (rich->self_trace < 0.30) return false;
    
    if (rich->self_stn_status && rich->self_stn_dist > 3.4) return false;
    
    // Chi-square Cut
    const double nchi_cut = 1.75;
    
    // Velocity
    //if (!hyc->vel_status[0]) return false;
    //if (std::log(hyc->vel_nchi[0]) > nchi_cut) return false;
    
    if (!hyc->vel_status[1]) return false;
    if (std::log(hyc->vel_nchi[1]) > 2.00) return false;
   
    // Geometry
    if (trk->lay[1] == 0) return false;
    if (!hyc->geom_status[0]) return false;
    if (std::log(hyc->geom_nchi_x[0]) > nchi_cut) return false;
    if (std::log(hyc->geom_nchi_y[0]) > nchi_cut) return false;

    // Mass
    //if (!hyc->mutr_status[0]) return false;
    //if (std::log(hyc->mutr_nchi_x[0]) > nchi_cut) return false;
    //if (std::log(hyc->mutr_nchi_y[0]) > nchi_cut) return false;
    //if (std::log(hyc->mutr_nchi_b[0]) > nchi_cut) return false;
    
    if (!hyc->mutr_status[1]) return false;
    if (std::log(hyc->mutr_nchi_x[1]) > nchi_cut) return false;
    if (std::log(hyc->mutr_nchi_y[1]) > nchi_cut) return false;
    if (std::log(hyc->mutr_nchi_b[1]) > 2.00) return false;

    // Physics
    //if (!hyc->phys_status[0][0]) return false;
    //if (std::log(hyc->phys_nchi_x[0][0]) > nchi_cut) return false;
    //if (std::log(hyc->phys_nchi_y[0][0]) > nchi_cut) return false;
    //if (std::log(hyc->phys_nchi_b[0][0]) > nchi_cut) return false;
    
    if (!hyc->phys_status[0][1]) return false;
    if (std::log(hyc->phys_nchi_x[0][1]) > nchi_cut) return false;
    if (std::log(hyc->phys_nchi_y[0][1]) > nchi_cut) return false;
    if (std::log(hyc->phys_nchi_b[0][1]) > 2.00) return false;
    
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
                       std::log(hyc->mutr_nchi_x[2]) < nchi_cut &&
                       std::log(hyc->mutr_nchi_y[2]) < nchi_cut &&
                       std::log(hyc->mutr_nchi_b[2]) < nchi_cut &&
                       hyc->mutr_sqrm[2] < 0.4);
    
    // Cut on in chi-square
    const double nchi_geom_cut = 0.75 + 0.5 * std::erfc((abs_rig - 25.0) / 20.0);
    if (std::log(hyc->geom_nchi_y[0]) > nchi_geom_cut) return false;
   
    // TRD extra hit
    if (signr > 0) Hist::Head("hLP_extseg")->fillH2D(abs_rig, trd->num_extra_seg, list->weight);
    if (signr < 0) Hist::Head("hLN_extseg")->fillH2D(abs_rig, trd->num_extra_seg, list->weight);
    if (signr > 0) Hist::Head("hLP_exthit")->fillH2D(abs_rig, trd->num_extra_hit, list->weight);
    if (signr < 0) Hist::Head("hLN_exthit")->fillH2D(abs_rig, trd->num_extra_hit, list->weight);
    //if (trd->num_extra_hit >= static_cast<int>(8.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig * abs_rig))) return false;
    
    // testcode
    //if (trd->num_vtx[0] > 2) return false;
    //if (trd->num_vtx[1] > 7) return false;
    //int nvtx_x = static_cast<int>(3.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig * abs_rig));
    //int nvtx_y = 2 + 2 * static_cast<int>(3.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig * abs_rig));
    int nvtx_x = static_cast<int>(3.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig));
    int nvtx_y = 2 + 2 * static_cast<int>(3.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig));
    if (trd->num_vtx[0] >= nvtx_x) return false;
    if (trd->num_vtx[1] >= nvtx_y) return false;
    
    // Cutoff
    //if (CheckType(Type::ISS) && cfr < 0.80) return false;
    //if (CheckType(Type::ISS) && cfr < 0.75) return false; // testcode
    if (CheckType(Type::ISS) && cfr < 0.70) return false; // testcode
    
    if (signr > 0 && !rich->self_cld_status) Hist::Head("hLP_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && !rich->self_cld_status) Hist::Head("hLN_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hLP_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hLN_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_pi) Hist::Head("hLN_sqrm_pi")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTLP_sqrm_pr")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTLN_sqrm_el")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_pi) Hist::Head("hTLN_sqrm_pi")->fillH3D(list->utime, abs_rig, sqrm, list->weight);

    bool is_rich_pr = rich->self_cld_status &&
                      hyc->vel_status[2] &&
                      std::log(hyc->vel_nchi[2]) < 2.00 &&
                      hyc->mutr_status[2] &&
                      std::log(hyc->mutr_nchi_x[2]) < nchi_cut &&
                      std::log(hyc->mutr_nchi_y[2]) < nchi_cut &&
                      std::log(hyc->mutr_nchi_b[2]) < 2.00 &&
                      hyc->phys_status[0][2] &&
                      std::log(hyc->phys_nchi_x[0][2]) < nchi_cut &&
                      std::log(hyc->phys_nchi_y[0][2]) < nchi_cut &&
                      std::log(hyc->phys_nchi_b[0][2]) < 2.00;

    bool is_sel_pr = trd->tdLLR_ep > 0.75 && (!rich->self_cld_status || is_rich_pr);

    if (!is_sel_pr) return false;

    if (signr > 0) Hist::Head("hLP_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hLN_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
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
   
    if (!rich->self_cld_status) return false;
    
    // Chi-square Cut
    const double nchi_cut = 1.75;
    
    // Velocity
    //if (!hyc->vel_status[1]) return false;
    //if (std::log(hyc->vel_nchi[1]) > nchi_cut) return false;
    
    if (!hyc->vel_status[2]) return false;
    if (std::log(hyc->vel_nchi[2]) > 2.00) return false;
    
    // Geometry
    if (trk->lay[1] == 0) return false;
    if (!hyc->geom_status[0]) return false;
    if (std::log(hyc->geom_nchi_x[0]) > nchi_cut) return false;
    if (std::log(hyc->geom_nchi_y[0]) > nchi_cut) return false;
    
    // Mass
    //if (!hyc->mutr_status[1]) return false;
    //if (std::log(hyc->mutr_nchi_x[1]) > nchi_cut) return false;
    //if (std::log(hyc->mutr_nchi_y[1]) > nchi_cut) return false;
    //if (std::log(hyc->mutr_nchi_b[1]) > nchi_cut) return false;

    if (!hyc->mutr_status[2]) return false;
    if (std::log(hyc->mutr_nchi_x[2]) > nchi_cut) return false;
    if (std::log(hyc->mutr_nchi_y[2]) > nchi_cut) return false;
    if (std::log(hyc->mutr_nchi_b[2]) > 2.00) return false;

    // Physics
    //if (!hyc->phys_status[0][1]) return false;
    //if (std::log(hyc->phys_nchi_x[0][1]) > nchi_cut) return false;
    //if (std::log(hyc->phys_nchi_y[0][1]) > nchi_cut) return false;
    //if (std::log(hyc->phys_nchi_b[0][1]) > nchi_cut) return false;

    if (!hyc->phys_status[0][2]) return false;
    if (std::log(hyc->phys_nchi_x[0][2]) > nchi_cut) return false;
    if (std::log(hyc->phys_nchi_y[0][2]) > nchi_cut) return false;
    if (std::log(hyc->phys_nchi_b[0][2]) > 2.00) return false;

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
    //if (CheckType(Type::ISS) && cfr < 0.80) return false;
    //if (CheckType(Type::ISS) && cfr < 0.75) return false; // testcode
    if (CheckType(Type::ISS) && cfr < 0.70) return false; // testcode
    
    if (signr > 0) Hist::Head("hMP_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0) Hist::Head("hMN_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hMP_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hMN_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTMP_sqrm_pr")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTMN_sqrm_el")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
  
    if (trdllr < 0.75) return false;

    if (signr > 0) Hist::Head("hMP_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hMN_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
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

    // Chi-square Cut
    const double nchi_cut = 1.75;
    
    // Geometry
    if (trk->lay[1] == 0) return false;
    if (!hyc->geom_status[0]) return false;

    // Variables
    double trdllr  = trd->tdLLR_ep;
    
    double rig     = hyc->geom_top_rig[0];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;

    // Cutoff
    //if (CheckType(Type::ISS) && cfr < 0.80) return false;
    //if (CheckType(Type::ISS) && cfr < 0.75) return false; // testcode
    if (CheckType(Type::ISS) && cfr < 0.70) return false; // testcode
    
    // Cut on in chi-square
    double lchix = std::log(hyc->geom_nchi_x[0]);
    double lchiy = std::log(hyc->geom_nchi_y[0]);

    if (lchix > nchi_cut) return false;
    
    if (signr > 0) Hist::Head("hIP_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr < 0) Hist::Head("hIN_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    
    const double nchi_geom_cut = 0.75 + 0.5 * std::erfc((abs_rig - 25.0) / 20.0);
    if (lchiy > nchi_geom_cut) return false;
  
    /*
    if (rich->self_status && rich->self_kind == 1 && rich->self_cld_status) {
        // Velocity
        if (hyc->vel_status[2]) {
            if (std::log(hyc->vel_nchi[2]) > nchi_cut) return false;
        }
        
        // Mass
        if (hyc->mutr_status[2]) {
            if (std::log(hyc->mutr_nchi_x[2]) > nchi_cut) return false;
            if (std::log(hyc->mutr_nchi_y[2]) > nchi_cut) return false;
            if (std::log(hyc->mutr_nchi_b[2]) > nchi_cut) return false;
        }

        // Physics
        if (hyc->phys_status[0][2]) {
            if (std::log(hyc->phys_nchi_x[0][2]) > nchi_cut) return false;
            if (std::log(hyc->phys_nchi_y[0][2]) > nchi_cut) return false;
            if (std::log(hyc->phys_nchi_b[0][2]) > nchi_cut) return false;
        }
    }
    */
    
    bool is_like_pr = (ecal->status && ecal->mvaBDT < -0.6);
    bool is_like_el = (ecal->status && ecal->mvaBDT >  0.2);
    
    if (signr > 0 && is_like_pr) Hist::Head("hIP_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hIN_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTIP_llr_pr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTIN_llr_el")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    if (signr > 0) Hist::Head("hIP_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0) Hist::Head("hIN_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
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
    
    double mc_wgt_flux = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    
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

    bool ptfs = hyc->geom_status[0] && hyc->geom_status[1] && hyc->geom_status[2] && hyc->geom_status[3];
    bool ptl9 = hyc->geom_status[0] && hyc->geom_status[2];
    bool ptl1 = hyc->geom_status[0] && hyc->geom_status[1];

    bool swfs = ptfs;
    bool swl9 = ptl9 && !ptfs;
    bool swl1 = ptl1 && !ptfs;
   
    if (!swfs && !swl9 && !swl1) return false;

    // Chi-square Cut
    const double nchix_cut = 2.25;
    const double nchiy_cut = 2.00;

    while (swl1) {
        double rig     = hyc->geom_top_rig[1];
        double abs_rig = std::abs(rig);
        short  signr   = (rig > 0) ? 1 : -1;
    
        double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
        double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
        
        double lchix_in = std::log(hyc->geom_nchi_x[0]);
        double lchiy_in = std::log(hyc->geom_nchi_y[0]);

        double lchix = std::log(hyc->geom_nchi_x[1]);
        double lchiy = std::log(hyc->geom_nchi_y[1]);
    
        double drig = hyc->geom_top_rig[1] / hyc->geom_top_rig[0] - 1.0;
        double rsin = std::sqrt(frso_pr_in_hc[0]*frso_pr_in_hc[0] + (hyc->geom_top_rig[0]/frso_pr_in_hc[1])*(hyc->geom_top_rig[0]/frso_pr_in_hc[1]));
        double rsl1 = std::sqrt(frso_pr_l1_hc[0]*frso_pr_l1_hc[0] + (hyc->geom_top_rig[1]/frso_pr_l1_hc[1])*(hyc->geom_top_rig[1]/frso_pr_l1_hc[1]));
        double reso = std::sqrt(rsl1 * rsl1 + rsin * rsin);
        double rvar = drig / reso;
        
        // Cutoff
        //if (CheckType(Type::ISS) && cfr < 0.80) break;
        //if (CheckType(Type::ISS) && cfr < 0.75) return false; // testcode
        if (CheckType(Type::ISS) && cfr < 0.70) return false; // testcode
        
        if (signr > 0) Hist::Head("hHPl1_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (signr < 0) Hist::Head("hHNl1_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (trdllr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPl1_lchix_in")->fillH2D(abs_rig, lchix_in, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchix_in")->fillH2D(abs_rig, lchix_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchix")->fillH2D(abs_rig, lchix, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchix")->fillH2D(abs_rig, lchix, list->weight);

        if (lchix_in > nchix_cut) break;
        if (lchix    > nchix_cut) break;

        if (signr > 0) Hist::Head("hHPl1_rvar")->fillH2D(abs_rig, rvar, list->weight);
        if (signr < 0) Hist::Head("hHNl1_rvar")->fillH2D(abs_rig, rvar, list->weight);

        if (signr > 0) Hist::Head("hHPl1_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy")->fillH2D(abs_rig, lchiy, list->weight);

        if (lchiy_in > nchiy_cut) break;
       
        double lrvar = std::abs(rvar);
        if (signr > 0) Hist::Head("hHPl1_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_MC_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight * mc_wgt_flux);
        if (signr < 0) Hist::Head("hHNl1_MC_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight * mc_wgt_flux);
        
        if (signr > 0) Hist::Head("hHPl1_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_cnt")->fillH1D(abs_rig, list->weight);
        if (signr < 0) Hist::Head("hHNl1_cnt")->fillH1D(abs_rig, list->weight);
        
        if (g4mc != nullptr) Hist::Head("hHl1_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);

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
        
        double lchix = std::log(hyc->geom_nchi_x[2]);
        double lchiy = std::log(hyc->geom_nchi_y[2]);
    
        double drig = hyc->geom_top_rig[2] / hyc->geom_top_rig[0] - 1.0;
        double rsin = std::sqrt(frso_pr_in_hc[0]*frso_pr_in_hc[0] + (hyc->geom_top_rig[0]/frso_pr_in_hc[1])*(hyc->geom_top_rig[0]/frso_pr_in_hc[1]));
        double rsl9 = std::sqrt(frso_pr_l9_hc[0]*frso_pr_l9_hc[0] + (hyc->geom_top_rig[2]/frso_pr_l9_hc[1])*(hyc->geom_top_rig[2]/frso_pr_l9_hc[1]));
        double reso = std::sqrt(rsl9 * rsl9 + rsin * rsin);
        double rvar = drig / reso;

        // Cutoff
        //if (CheckType(Type::ISS) && cfr < 0.80) break;
        //if (CheckType(Type::ISS) && cfr < 0.75) return false; // testcode
        if (CheckType(Type::ISS) && cfr < 0.70) return false; // testcode
        
        if (signr > 0) Hist::Head("hHPl9_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (signr < 0) Hist::Head("hHNl9_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (trdllr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPl9_lchix_in")->fillH2D(abs_rig, lchix_in, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchix_in")->fillH2D(abs_rig, lchix_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchix")->fillH2D(abs_rig, lchix, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchix")->fillH2D(abs_rig, lchix, list->weight);

        if (lchix_in > nchix_cut) break;
        if (lchix    > nchix_cut) break;

        if (signr > 0) Hist::Head("hHPl9_rvar")->fillH2D(abs_rig, rvar, list->weight);
        if (signr < 0) Hist::Head("hHNl9_rvar")->fillH2D(abs_rig, rvar, list->weight);

        if (signr > 0) Hist::Head("hHPl9_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (lchiy_in > nchiy_cut) break;

        double lrvar = std::abs(rvar);
        if (signr > 0) Hist::Head("hHPl9_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_MC_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight * mc_wgt_flux);
        if (signr < 0) Hist::Head("hHNl9_MC_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight * mc_wgt_flux);
        
        if (signr > 0) Hist::Head("hHPl9_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lrvar_cut")->fillH2D(abs_rig, lrvar, list->weight);

        if (signr > 0) Hist::Head("hHPl9_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_cnt")->fillH1D(abs_rig, list->weight);
        if (signr < 0) Hist::Head("hHNl9_cnt")->fillH1D(abs_rig, list->weight);
        
        if (g4mc != nullptr) Hist::Head("hHl9_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);

        break;
    }

    while (swfs) {
        double rig     = hyc->geom_top_rig[3];
        double abs_rig = std::abs(rig);
        short  signr   = (rig > 0) ? 1 : -1;
    
        double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
        double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
        
        double lchix_l1 = std::log(hyc->geom_nchi_x[1]);
        double lchiy_l1 = std::log(hyc->geom_nchi_y[1]);
        
        double lchix_l9 = std::log(hyc->geom_nchi_x[2]);
        double lchiy_l9 = std::log(hyc->geom_nchi_y[2]);

        double lchix = std::log(hyc->geom_nchi_x[3]);
        double lchiy = std::log(hyc->geom_nchi_y[3]);
        
        double rvarl1 = 0;
        {
            double drig = hyc->geom_top_rig[3] / hyc->geom_top_rig[1] - 1.0;
            double rsl1 = std::sqrt(frso_pr_l1_hc[0]*frso_pr_l1_hc[0] + (hyc->geom_top_rig[1]/frso_pr_l1_hc[1])*(hyc->geom_top_rig[1]/frso_pr_l1_hc[1]));
            double rsfs = std::sqrt(frso_pr_fs_hc[0]*frso_pr_fs_hc[0] + (hyc->geom_top_rig[3]/frso_pr_fs_hc[1])*(hyc->geom_top_rig[3]/frso_pr_fs_hc[1]));
            double reso = std::sqrt(rsl1 * rsl1 + rsfs * rsfs);
            rvarl1 = drig / reso;
        }

        double rvarl9 = 0;
        {
            double drig = hyc->geom_top_rig[3] / hyc->geom_top_rig[2] - 1.0;
            double rsl9 = std::sqrt(frso_pr_l9_hc[0]*frso_pr_l9_hc[0] + (hyc->geom_top_rig[2]/frso_pr_l9_hc[1])*(hyc->geom_top_rig[2]/frso_pr_l9_hc[1]));
            double rsfs = std::sqrt(frso_pr_fs_hc[0]*frso_pr_fs_hc[0] + (hyc->geom_top_rig[3]/frso_pr_fs_hc[1])*(hyc->geom_top_rig[3]/frso_pr_fs_hc[1]));
            double reso = std::sqrt(rsl9 * rsl9 + rsfs * rsfs);
            rvarl9 = drig / reso;

        }

        double rvarfs = 0;
        {
            double drig = rvarfs = hyc->geom_top_rig[1] / hyc->geom_top_rig[2] - 1.0;
            double rsl1 = std::sqrt(frso_pr_l1_hc[0]*frso_pr_l1_hc[0] + (hyc->geom_top_rig[1]/frso_pr_l1_hc[1])*(hyc->geom_top_rig[1]/frso_pr_l1_hc[1]));
            double rsl9 = std::sqrt(frso_pr_l9_hc[0]*frso_pr_l9_hc[0] + (hyc->geom_top_rig[2]/frso_pr_l9_hc[1])*(hyc->geom_top_rig[2]/frso_pr_l9_hc[1]));
            double rsfs = std::sqrt(frso_pr_fs_hc[0]*frso_pr_fs_hc[0] + (hyc->geom_top_rig[3]/frso_pr_fs_hc[1])*(hyc->geom_top_rig[3]/frso_pr_fs_hc[1]));
            double rs19 = std::sqrt(rsl1 * rsl1 + rsl9 * rsl9);
            double reso = std::sqrt(rs19 * rs19 + rsfs * rsfs);
            rvarfs = drig / reso;
        }

        // Cutoff
        //if (CheckType(Type::ISS) && cfr < 0.80) break;
        //if (CheckType(Type::ISS) && cfr < 0.75) return false; // testcode
        if (CheckType(Type::ISS) && cfr < 0.70) return false; // testcode
        
        if (signr > 0) Hist::Head("hHPfs_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (signr < 0) Hist::Head("hHNfs_llr")->fillH2D(abs_rig, trdllr, list->weight);
        if (trdllr < 0.75) break;
        
        if (signr > 0) Hist::Head("hHPfs_lchix_l1")->fillH2D(abs_rig, lchix_l1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchix_l1")->fillH2D(abs_rig, lchix_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchix_l9")->fillH2D(abs_rig, lchix_l9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchix_l9")->fillH2D(abs_rig, lchix_l9, list->weight);
       
        if (signr > 0) Hist::Head("hHPfs_lchix")->fillH2D(abs_rig, lchix, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchix")->fillH2D(abs_rig, lchix, list->weight);
        
        if (lchix_l1 > nchix_cut) break;
        if (lchix_l9 > nchix_cut) break;
        if (lchix    > nchix_cut) break;
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (lchiy_l1 > nchiy_cut) break;
        if (lchiy_l9 > nchiy_cut) break;
        
        if (signr > 0) Hist::Head("hHPfs_rvar_l1")->fillH2D(abs_rig, rvarl1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_rvar_l1")->fillH2D(abs_rig, rvarl1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_rvar_l9")->fillH2D(abs_rig, rvarl9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_rvar_l9")->fillH2D(abs_rig, rvarl9, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_rvar")->fillH2D(abs_rig, rvarfs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_rvar")->fillH2D(abs_rig, rvarfs, list->weight);
    
        //if (std::abs(rvarl1) > rvar_cut) break;
        //if (std::abs(rvarl9) > rvar_cut) break;
        
        double lrvarl1 = std::abs(rvarl1);
        double lrvarl9 = std::abs(rvarl9);
        double lrvarfs = std::abs(rvarfs);
        if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        
        if (std::abs(rvarl1) < 4.0 && std::abs(rvarl9) < 4.0) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar_v1")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar_v1")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }
        if (std::abs(rvarl1) < 3.5 && std::abs(rvarl9) < 3.5) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar_v2")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar_v2")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }
        if (std::abs(rvarl1) < 3.0 && std::abs(rvarl9) < 3.0) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar_v3")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar_v3")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }

        //////////////////////////////////////////////////////////////////////////////////
        if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar2")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar2")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        
        if (std::abs(rvarl1) < 4.0 && std::abs(rvarl9) < 4.0) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar2_v1")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar2_v1")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }
        if (std::abs(rvarl1) < 3.5 && std::abs(rvarl9) < 3.5) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar2_v2")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar2_v2")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }
        if (std::abs(rvarl1) < 3.0 && std::abs(rvarl9) < 3.0) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar2_v3")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar2_v3")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }
        //////////////////////////////////////////////////////////////////////////////////
        
        //////////////////////////////////////////////////////////////////////////////////
        if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar3")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar3")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        
        if (std::abs(rvarl1) < 4.0 && std::abs(rvarl9) < 4.0) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar3_v1")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar3_v1")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }
        if (std::abs(rvarl1) < 3.5 && std::abs(rvarl9) < 3.5) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar3_v2")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar3_v2")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }
        if (std::abs(rvarl1) < 3.0 && std::abs(rvarl9) < 3.0) {
            if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar3_v3")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar3_v3")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight);
        }
        //////////////////////////////////////////////////////////////////////////////////

        if (signr > 0) Hist::Head("hHPfs_MC_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight * mc_wgt_flux);
        if (signr < 0) Hist::Head("hHNfs_MC_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvarfs, list->weight * mc_wgt_flux);
        
        if (signr > 0) Hist::Head("hHPfs_lrvarl1_cut")->fillH2D(abs_rig, lrvarl1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lrvarl1_cut")->fillH2D(abs_rig, lrvarl1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lrvarl9_cut")->fillH2D(abs_rig, lrvarl9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lrvarl9_cut")->fillH2D(abs_rig, lrvarl9, list->weight);

        if (signr > 0) Hist::Head("hHPfs_lrvar_cut")->fillH2D(abs_rig, lrvarfs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lrvar_cut")->fillH2D(abs_rig, lrvarfs, list->weight);

        if (signr > 0) Hist::Head("hHPfs_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_cut")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (std::abs(rvarl1) < 5.0 && std::abs(rvarl9) < 5.0) {
            if (signr > 0) Hist::Head("hHPfs_lrvar_cut1")->fillH2D(abs_rig, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lrvar_cut1")->fillH2D(abs_rig, lrvarfs, list->weight);

            if (signr > 0) Hist::Head("hHPfs_lchiy_cut1")->fillH2D(abs_rig, lchiy, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_cut1")->fillH2D(abs_rig, lchiy, list->weight);
        }
        if (std::abs(rvarl1) < 4.0 && std::abs(rvarl9) < 4.0) {
            if (signr > 0) Hist::Head("hHPfs_lrvar_cut2")->fillH2D(abs_rig, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lrvar_cut2")->fillH2D(abs_rig, lrvarfs, list->weight);

            if (signr > 0) Hist::Head("hHPfs_lchiy_cut2")->fillH2D(abs_rig, lchiy, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_cut2")->fillH2D(abs_rig, lchiy, list->weight);
        }
        if (std::abs(rvarl1) < 3.0 && std::abs(rvarl9) < 3.0) {
            if (signr > 0) Hist::Head("hHPfs_lrvar_cut3")->fillH2D(abs_rig, lrvarfs, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lrvar_cut3")->fillH2D(abs_rig, lrvarfs, list->weight);

            if (signr > 0) Hist::Head("hHPfs_lchiy_cut3")->fillH2D(abs_rig, lchiy, list->weight);
            if (signr < 0) Hist::Head("hHNfs_lchiy_cut3")->fillH2D(abs_rig, lchiy, list->weight);
        }

        if (signr > 0) Hist::Head("hHPfs_cnt")->fillH1D(abs_rig, list->weight);
        if (signr < 0) Hist::Head("hHNfs_cnt")->fillH1D(abs_rig, list->weight);
        
        if (g4mc != nullptr) Hist::Head("hHfs_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);

        break;
    }

    return true;
}


#endif // __Analyzer_C__
