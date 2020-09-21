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
    
    std::vector<double> vtme2( {
        1305417600, 1307750400, 1310083200, 1312416000, 1314748800, 
        1317081600, 1319414400, 1321747200, 1324080000, 1326412800, 
        1328745600, 1331078400, 1333411200, 1335744000, 1338076800, 
        1340409600, 1342742400, 1345075200, 1347408000, 1349740800, 
        1352073600, 1354406400, 1356739200, 1359072000, 1361404800, 
        1363737600, 1366070400, 1368403200, 1370736000, 1373068800, 
        1375401600, 1377734400, 1380067200, 1382400000, 1384732800, 
        1387065600, 1389398400, 1391731200, 1394064000, 1396396800, 
        1398729600, 1401062400, 1403395200, 1405728000, 1408060800, 
        1410393600, 1412726400, 1415059200, 1417392000, 1419724800, 
        1422057600, 1424390400, 1426723200, 1429056000, 1431388800, 
        1433721600, 1436054400, 1438387200, 1440720000, 1443052800, 
        1445385600, 1447718400, 1450051200, 1452384000, 1454716800, 
        1457049600, 1459382400, 1461715200, 1464048000, 1466380800, 
        1468713600, 1471046400, 1473379200, 1475712000, 1478044800, 
        1480377600, 1482710400, 1485043200, 1487376000, 1489708800, 
        1492041600, 1494374400, 1496707200, 1499040000, 1501372800, 
        1503705600, 1506038400, 1508371200, 1510704000, 1513036800, 
        1515369600, 1517702400, 1520035200, 1522368000, 1524700800, 
        1527033600, 1529366400, 1531699200, 1534032000, 1536364800, 
        1538697600, 1541030400, 1543363200, 1545696000, 1548028800, 
        1550361600 } );
    vtme = vtme2;

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
    
    std::vector<double> vmom_el( { // e+- binning
        0.50 ,0.65 ,0.82, // 3
        1.01, 2.0, 3.0, 4.12, 5.0, 6.0, 7.1, 8.3, 9.62, 11.04,
        12.59, 14.25, 16.05, 17.98, 20.04, 22.25, 24.62, 27.25, 30.21, 35.36, 40 ,// 23
        41.9 , 45.1 ,48.5 ,52.2 ,56.1 ,60.3 ,64.8 ,69.7 ,74.9 ,80.5  ,86.5, // 11
        93.0 ,100.  ,108. ,116. ,125. ,135. ,147. ,160  ,175. ,192.,//10
        211. ,233.  ,259. ,291. ,330. ,379. ,441. ,525. ,643  ,822. ,1130. ,1800.//12
    } );
    
    std::vector<double> vmom_ap( {
        1.00, 1.51, 1.92,  2.97,  3.64,  4.43,  5.37,  5.90,  6.47,  7.09, 
        7.76, 8.48, 9.26, 11.00, 13.00, 15.30, 18.00, 21.10, 24.70, 28.80, 
        33.5 } );
    
    Axis AXrig ("Rigidity [GV]", vmom);
    Axis AXTrig("Rigidity [GV]", vmom_ap);
    
    Axis AXLllr("TRD estimator", 200, 0.0, 1.6);
    Hist::New("hLP_llr", HistAxis(AXrig, AXLllr));
    Hist::New("hLN_llr", HistAxis(AXrig, AXLllr));

    Axis AXLsqrm ("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 75, -2.5, 3.5);
    Axis AXTLsqrm("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 50, -2.5, 3.5);
    
    Hist::New("hLP_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLP_sqrm_pr",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_sqrm_el",   HistAxis(AXrig, AXLsqrm));
    
    Hist::New("hTLP_sqrm",      HistAxis(AXtme, AXTrig, AXTLsqrm));
    Hist::New("hTLN_sqrm",      HistAxis(AXtme, AXTrig, AXTLsqrm));
    Hist::New("hTLP_sqrm_pr",   HistAxis(AXtme, AXTrig, AXTLsqrm));
    Hist::New("hTLN_sqrm_el",   HistAxis(AXtme, AXTrig, AXTLsqrm));
    Hist::New("hTRLP_sqrm_pr",   HistAxis(AXTrig, AXTLsqrm));
    Hist::New("hTRLN_sqrm_el",   HistAxis(AXTrig, AXTLsqrm));

    Hist::New("hLP_cf85_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf85_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLP_cf85_sqrm_pr",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf85_sqrm_el",   HistAxis(AXrig, AXLsqrm));
    
    Hist::New("hLP_cf95_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf95_sqrm",      HistAxis(AXrig, AXLsqrm));
    Hist::New("hLP_cf95_sqrm_pr",   HistAxis(AXrig, AXLsqrm));
    Hist::New("hLN_cf95_sqrm_el",   HistAxis(AXrig, AXLsqrm));
    
    Hist::New("hLP_cnt", HistAxis(AXrig));
    Hist::New("hLN_cnt", HistAxis(AXrig));
    
    Hist::New("hL_MC_cnt", HistAxis(AXrig));
    Hist::New("hTL_MC_cnt", HistAxis(AXTrig));
    
    Axis AXMllr("TRD estimator", 100, 0.0, 1.6);
    Hist::New("hMP_llr", HistAxis(AXrig, AXMllr));
    Hist::New("hMN_llr", HistAxis(AXrig, AXMllr));
    
    Axis AXMsqrm ("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 100, -1.75, 3.25);
    Axis AXTMsqrm("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 75, -1.75, 3.25);
    Hist::New("hMP_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMP_sqrm_pr", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_sqrm_el", HistAxis(AXrig, AXMsqrm));

    Hist::New("hTMP_sqrm",    HistAxis(AXtme, AXTrig, AXTMsqrm));
    Hist::New("hTMN_sqrm",    HistAxis(AXtme, AXTrig, AXTMsqrm));
    Hist::New("hTMP_sqrm_pr", HistAxis(AXtme, AXTrig, AXTMsqrm));
    Hist::New("hTMN_sqrm_el", HistAxis(AXtme, AXTrig, AXTMsqrm));
    Hist::New("hTRMP_sqrm_pr", HistAxis(AXTrig, AXTMsqrm));
    Hist::New("hTRMN_sqrm_el", HistAxis(AXTrig, AXTMsqrm));

    Hist::New("hMP_cf85_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf85_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMP_cf85_sqrm_pr", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf85_sqrm_el", HistAxis(AXrig, AXMsqrm));

    Hist::New("hMP_cf95_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf95_sqrm",    HistAxis(AXrig, AXMsqrm));
    Hist::New("hMP_cf95_sqrm_pr", HistAxis(AXrig, AXMsqrm));
    Hist::New("hMN_cf95_sqrm_el", HistAxis(AXrig, AXMsqrm));
    
    Hist::New("hMP_cnt", HistAxis(AXrig));
    Hist::New("hMN_cnt", HistAxis(AXrig));
    
    Hist::New("hM_MC_cnt", HistAxis(AXrig));
    Hist::New("hTM_MC_cnt", HistAxis(AXTrig));
    
    Axis AXIlchi("Log(#chi^{2}/DOF)", 100, -4, 8);
    Hist::New("hIP_lchiy", HistAxis(AXrig, AXIlchi));
    Hist::New("hIN_lchiy", HistAxis(AXrig, AXIlchi));
    
    Axis AXIllr ("TRD estimator", 150, 0.0, 1.6);
    Axis AXTIllr("TRD estimator", 100, 0.0, 1.6);
    Hist::New("hIP_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIP_llr_pr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_llr_el", HistAxis(AXrig, AXIllr));
    
    Hist::New("hTIP_llr", HistAxis(AXtme, AXTrig, AXTIllr));
    Hist::New("hTIN_llr", HistAxis(AXtme, AXTrig, AXTIllr));
    Hist::New("hTIP_llr_pr", HistAxis(AXtme, AXTrig, AXTIllr));
    Hist::New("hTIN_llr_el", HistAxis(AXtme, AXTrig, AXTIllr));
    Hist::New("hTRIP_llr_pr", HistAxis(AXTrig, AXTIllr));
    Hist::New("hTRIN_llr_el", HistAxis(AXTrig, AXTIllr));
    
    Hist::New("hIP_cf85_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf85_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIP_cf85_llr_pr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf85_llr_el", HistAxis(AXrig, AXIllr));
    
    Hist::New("hIP_cf95_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf95_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIP_cf95_llr_pr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_cf95_llr_el", HistAxis(AXrig, AXIllr));
    
    Hist::New("hIP_RH_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_RH_llr", HistAxis(AXrig, AXIllr));
    Hist::New("hIP_RH_llr_pr", HistAxis(AXrig, AXIllr));
    Hist::New("hIN_RH_llr_el", HistAxis(AXrig, AXIllr));
    
    Hist::New("hTIP_RH_llr", HistAxis(AXtme, AXTrig, AXTIllr));
    Hist::New("hTIN_RH_llr", HistAxis(AXtme, AXTrig, AXTIllr));
    Hist::New("hTIP_RH_llr_pr", HistAxis(AXtme, AXTrig, AXTIllr));
    Hist::New("hTIN_RH_llr_el", HistAxis(AXtme, AXTrig, AXTIllr));
    Hist::New("hTRIP_RH_llr_pr", HistAxis(AXTrig, AXTIllr));
    Hist::New("hTRIN_RH_llr_el", HistAxis(AXTrig, AXTIllr));
    
    Hist::New("hIP_cnt", HistAxis(AXrig));
    Hist::New("hIN_cnt", HistAxis(AXrig));
    
    Hist::New("hI_MC_cnt", HistAxis(AXrig));
    Hist::New("hTI_MC_cnt", HistAxis(AXTrig));
    
    Axis AXHllr("TRD estimator", 100, 0.0, 1.6);
    Axis AXHlchi2("Log(#chi^{2}/DOF)", 200, -4, 4);
    Axis AXHlrvar("|RVAR|", 200, 0, 8);
    Axis AXHnorm("Norm", 200, 0, 8);
    
    Hist::New("hHPl1_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNl1_llr", HistAxis(AXrig, AXHllr));
    
    Hist::New("hHPl1_lx",  HistAxis(AXrig, AXHnorm));
    Hist::New("hHNl1_ly",  HistAxis(AXrig, AXHnorm));
    Hist::New("hHPl1_tau", HistAxis(AXrig, AXHnorm));
    Hist::New("hHNl1_rho", HistAxis(AXrig, AXHnorm));
    
    Hist::New("hHPl1_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl1_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl1_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl1_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl1_lrvar", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNl1_lrvar", HistAxis(AXrig, AXHlrvar));
    
    Hist::New("hHPl1_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNl1_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    
    Hist::New("hHPl1_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNl1_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHPl1_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNl1_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2, AXHlrvar));

    Hist::New("hHPl1_cnt", HistAxis(AXrig));
    Hist::New("hHNl1_cnt", HistAxis(AXrig));
    
    Hist::New("hHl1_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hHl1_cnt_MC_FLUX27", HistAxis(AXrig));
    
    Hist::New("hHPl9_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNl9_llr", HistAxis(AXrig, AXHllr));
    
    Hist::New("hHPl9_lx",  HistAxis(AXrig, AXHnorm));
    Hist::New("hHNl9_ly",  HistAxis(AXrig, AXHnorm));
    Hist::New("hHPl9_tau", HistAxis(AXrig, AXHnorm));
    Hist::New("hHNl9_rho", HistAxis(AXrig, AXHnorm));
    
    Hist::New("hHPl9_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl9_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl9_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNl9_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPl9_lrvar", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNl9_lrvar", HistAxis(AXrig, AXHlrvar));
    
    Hist::New("hHPl9_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNl9_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    
    Hist::New("hHPl9_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNl9_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHPl9_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNl9_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2, AXHlrvar));

    Hist::New("hHPl9_cnt", HistAxis(AXrig));
    Hist::New("hHNl9_cnt", HistAxis(AXrig));
    
    Hist::New("hHl9_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hHl9_cnt_MC_FLUX27", HistAxis(AXrig));

    Hist::New("hHPfs_llr", HistAxis(AXrig, AXHllr));
    Hist::New("hHNfs_llr", HistAxis(AXrig, AXHllr));
    
    Hist::New("hHPfs_lx",  HistAxis(AXrig, AXHnorm));
    Hist::New("hHNfs_ly",  HistAxis(AXrig, AXHnorm));
    Hist::New("hHPfs_tau", HistAxis(AXrig, AXHnorm));
    Hist::New("hHNfs_rho", HistAxis(AXrig, AXHnorm));
    
    Hist::New("hHPfs_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchix", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPfs_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHNfs_lchiy", HistAxis(AXrig, AXHlchi2));
    Hist::New("hHPfs_lrvar", HistAxis(AXrig, AXHlrvar));
    Hist::New("hHNfs_lrvar", HistAxis(AXrig, AXHlrvar));
    
    Hist::New("hHPfs_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNfs_lchiy_lrvar", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    
    Hist::New("hHPfs_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNfs_lchiy_lrvar_MC_FLUX10", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHPfs_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2, AXHlrvar));
    Hist::New("hHNfs_lchiy_lrvar_MC_FLUX27", HistAxis(AXrig, AXHlchi2, AXHlrvar));

    Hist::New("hHPfs_cnt", HistAxis(AXrig));
    Hist::New("hHNfs_cnt", HistAxis(AXrig));
    
    Hist::New("hHfs_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hHfs_cnt_MC_FLUX27", HistAxis(AXrig));

    return true;
}

bool Analyzer::process_init() {
    varsL1.init();
    varsL9.init();
    varsFS.init();
    return true;
}

bool Analyzer::process_presel() {
    // RTI
    if (CheckType(Type::ISS)) {
        if (rti->is_in_SAA) return false;
        if (rti->livetime < 0.5) return false;
        if (std::abs(rti->tk_align[0][0]) > 35.0) return false;
        if (std::abs(rti->tk_align[0][1]) > 35.0) return false;
        if (std::abs(rti->tk_align[1][0]) > 45.0) return false;
        if (std::abs(rti->tk_align[1][1]) > 45.0) return false;
    }
    
    // Charge
    if (tof->Qall < 0.8 || tof->Qall > 1.4) return false;
    if (trk->QIn < 0.8 || trk->QIn > 1.4) return false;

    // TRD
    if (!trd->tdLLR_status || trd->tdLLR_num_hit < 8 || trd->num_tdHit < 6) return false;
    
    // Trigger
    if ((trg->bit&2) != 2 && (trg->bit&8) != 8) return false;
    
    return true;
}

bool Analyzer::process_data() {
    process_data_l();
    process_data_m();
    process_data_i();
    process_data_hl1();
    process_data_hl9();
    process_data_hfs();
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
    
    // Chi-squared Cut
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
    
    //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
    
    // Cut on chi-square
    const double nchi_geom_cut = 0.75 + 0.5 * std::erfc((abs_rig - 25.0) / 20.0);
    if (std::log(hyc->geom_nchi_y[0]) > nchi_geom_cut) return false;
   
    // TRD extra hit
    int td_nvtx_x = (trd->num_vtx[2][0] + trd->num_vtx[3][0]);
    int td_nvtx_y = (trd->num_vtx[2][1] + trd->num_vtx[3][1]);
    int tdcut_nvtx_x = 2 + static_cast<int>(2.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig));
    int tdcut_nvtx_y = 2 + static_cast<int>(2.0 * std::sqrt(abs_rig) * std::log(1.0 + abs_rig * abs_rig));
    if (td_nvtx_x > tdcut_nvtx_x) return false;
    if (td_nvtx_y > tdcut_nvtx_y) return false;
    
    bool is_rich_el = rich->self_cld_status &&
                      hyc->mutr_status[2] &&
                      std::log(hyc->mutr_nchi_x[2]) < trk_nchi_cut &&
                      std::log(hyc->mutr_nchi_y[2]) < trk_nchi_cut &&
                      std::log(hyc->mutr_nchi_b[2]) < vel_nchi_cut;
    
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
    
    double trdllr_crr = hyc->phys_top_bta[0][1];
    bool is_like_pr = (trd->tdLLR_ep > (0.75 * trdllr_crr) && (!rich->self_cld_status || is_rich_pr));
    bool is_like_el = (trd->tdLLR_ep < (0.75 * trdllr_crr) && is_rich_el);
     
    // Cutoff
    if (CheckType(Type::ISS) && cfr < 0.75) return false;

    if (signr > 0 && !rich->self_cld_status) Hist::Head("hLP_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && !rich->self_cld_status) Hist::Head("hLN_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hLP_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hLN_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);

    if (signr > 0 && is_like_pr && cfr > 0.85) Hist::Head("hLP_cf85_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.85) Hist::Head("hLN_cf85_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.95) Hist::Head("hLP_cf95_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.95) Hist::Head("hLN_cf95_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTLP_sqrm_pr")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTLN_sqrm_el")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr > 0 && is_like_pr) Hist::Head("hTRLP_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTRLN_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);

    if (!is_like_pr) return false;

    if (signr > 0) Hist::Head("hLP_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hLN_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && cfr > 0.85) Hist::Head("hLP_cf85_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.85) Hist::Head("hLN_cf85_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && cfr > 0.95) Hist::Head("hLP_cf95_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.95) Hist::Head("hLN_cf95_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
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
    
    //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
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
    
    if (signr > 0 && is_like_pr && cfr > 0.85) Hist::Head("hMP_cf85_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.85) Hist::Head("hMN_cf85_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.95) Hist::Head("hMP_cf95_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.95) Hist::Head("hMN_cf95_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTMP_sqrm_pr")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTMN_sqrm_el")->fillH3D(list->utime, abs_rig, sqrm, list->weight);
    if (signr > 0 && is_like_pr) Hist::Head("hTRMP_sqrm_pr")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTRMN_sqrm_el")->fillH2D(abs_rig, sqrm, list->weight);
  
    if (trdllr < 0.75) return false;

    if (signr > 0) Hist::Head("hMP_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0) Hist::Head("hMN_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
    if (signr > 0 && cfr > 0.85) Hist::Head("hMP_cf85_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.85) Hist::Head("hMN_cf85_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr > 0 && cfr > 0.95) Hist::Head("hMP_cf95_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    if (signr < 0 && cfr > 0.95) Hist::Head("hMN_cf95_sqrm")->fillH2D(abs_rig, sqrm, list->weight);
    
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
    
    //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
    
    // Cut on in chi-square
    double lchix = std::log(hyc->geom_nchi_x[0]);
    double lchiy = std::log(hyc->geom_nchi_y[0]);

    if (lchix > trk_nchi_cut) return false;
    
    if (signr > 0) Hist::Head("hIP_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr < 0) Hist::Head("hIN_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    
    const double nchi_geom_cut = 0.75 + 0.5 * std::erfc((abs_rig - 25.0) / 20.0);
    if (lchiy > nchi_geom_cut) return false;
    
    // RICH
    bool has_rich = rich->self_status && rich->self_kind == 1;
    if (has_rich && hyc->mutr_status[2]) {
        if (std::log(hyc->mutr_nchi_x[2]) > trk_nchi_cut) return false;
        if (std::log(hyc->mutr_nchi_y[2]) > trk_nchi_cut) return false;
        if (std::log(hyc->mutr_nchi_b[2]) > vel_nchi_cut) return false;
    }
    if (has_rich && hyc->phys_status[0][2]) {
        if (std::log(hyc->phys_nchi_x[0][2]) > trk_nchi_cut) return false;
        if (std::log(hyc->phys_nchi_y[0][2]) > trk_nchi_cut) return false;
        if (std::log(hyc->phys_nchi_b[0][2]) > vel_nchi_cut) return false;
    }
  
    bool is_like_pr = (ecal->status && ecal->mvaBDT < -0.6);
    bool is_like_el = (ecal->status && ecal->mvaBDT >  0.2);
    
    // Cutoff
    if (CheckType(Type::ISS) && cfr < 0.75) return false;
    
    if (signr > 0 && is_like_pr) Hist::Head("hIP_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hIN_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr && has_rich) Hist::Head("hIP_RH_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el && has_rich) Hist::Head("hIN_RH_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr && cfr > 0.85) Hist::Head("hIP_cf85_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.85) Hist::Head("hIN_cf85_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr > 0 && is_like_pr && cfr > 0.95) Hist::Head("hIP_cf95_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el && cfr > 0.95) Hist::Head("hIN_cf95_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr) Hist::Head("hTIP_llr_pr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTIN_llr_el")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr > 0 && is_like_pr) Hist::Head("hTRIP_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el) Hist::Head("hTRIN_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && is_like_pr && has_rich) Hist::Head("hTIP_RH_llr_pr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el && has_rich) Hist::Head("hTIN_RH_llr_el")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr > 0 && is_like_pr && has_rich) Hist::Head("hTRIP_RH_llr_pr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && is_like_el && has_rich) Hist::Head("hTRIN_RH_llr_el")->fillH2D(abs_rig, trdllr, list->weight);
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
   
    if (signr > 0) Hist::Head("hIP_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0) Hist::Head("hIN_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && has_rich) Hist::Head("hIP_RH_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && has_rich) Hist::Head("hIN_RH_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0 && cfr > 0.85) Hist::Head("hIP_cf85_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && cfr > 0.85) Hist::Head("hIN_cf85_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr > 0 && cfr > 0.95) Hist::Head("hIP_cf95_llr")->fillH2D(abs_rig, trdllr, list->weight);
    if (signr < 0 && cfr > 0.95) Hist::Head("hIN_cf95_llr")->fillH2D(abs_rig, trdllr, list->weight);
    
    if (signr > 0) Hist::Head("hTIP_llr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr < 0) Hist::Head("hTIN_llr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    
    if (signr > 0 && has_rich) Hist::Head("hTIP_RH_llr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    if (signr < 0 && has_rich) Hist::Head("hTIN_RH_llr")->fillH3D(list->utime, abs_rig, trdllr, list->weight);
    
    if (signr > 0) Hist::Head("hIP_cnt")->fillH1D(abs_rig, list->weight);
    if (signr < 0) Hist::Head("hIN_cnt")->fillH1D(abs_rig, list->weight);

    if (g4mc != nullptr) Hist::Head("hI_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);
    if (g4mc != nullptr) Hist::Head("hTI_MC_cnt")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), list->weight);

    return true;
}

bool Analyzer::process_data_hl1() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_flux10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_flux27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_flux10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_flux27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    // Tracker-L1
    if (!hyc->geom_status[0]) return false;
    if (!hyc->geom_status[1]) return false;
    if (trk->lay[0] != 3 || trk->lay[1] == 0) return false;
    if (trk->num_inn_x < 4) return false;
    if (trk->num_inn_y < 5) return false;

    double rig     = hyc->geom_top_rig[1];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
        
    // Cutoff
    if (CheckType(Type::ISS) && cfr < 0.75) return false;
    
    if (signr > 0) Hist::Head("hHPl1_llr")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
    if (signr < 0) Hist::Head("hHNl1_llr")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
    
    // TRD
    if (trd->tdLLR_ep < 0.75) return false;
    
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);

    double lchix_l1 = std::log(hyc->geom_nchi_x[1]);
    double lchiy_l1 = std::log(hyc->geom_nchi_y[1]);
    
    const double lchi_cut = 4.00;
    if (lchix_in > lchi_cut || lchix_l1 > lchi_cut) return false;
    if (lchiy_in > lchi_cut || lchiy_l1 > lchi_cut) return false;

    double lchix = (9.0/25.0) * lchix_in + (16.0/25.0) * lchix_l1;
    double lchiy = (9.0/25.0) * lchiy_in + (16.0/25.0) * lchiy_l1;
    
    double lrvar = std::log(1.0 + std::abs(hyc->geom_top_rig[1] / hyc->geom_top_rig[0] - 1.0));
    
    double nrm_lx  = hyc->geom_max_norm_lx[1];
    double nrm_ly  = hyc->geom_max_norm_ly[1];
    double nrm_tau = hyc->geom_max_norm_tau[1];
    double nrm_rho = hyc->geom_max_norm_rho[1];
    
    if (signr > 0) Hist::Head("hHPl1_lx")->fillH2D(abs_rig, nrm_lx, list->weight);
    if (signr < 0) Hist::Head("hHNl1_ly")->fillH2D(abs_rig, nrm_ly, list->weight);
    if (signr > 0) Hist::Head("hHPl1_tau")->fillH2D(abs_rig, nrm_tau, list->weight);
    if (signr < 0) Hist::Head("hHNl1_rho")->fillH2D(abs_rig, nrm_rho, list->weight);
    
    if (signr > 0) Hist::Head("hHPl1_lchix")->fillH2D(abs_rig, lchix, list->weight);
    if (signr < 0) Hist::Head("hHNl1_lchix")->fillH2D(abs_rig, lchix, list->weight);
    if (signr > 0) Hist::Head("hHPl1_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr < 0) Hist::Head("hHNl1_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr > 0) Hist::Head("hHPl1_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
    if (signr < 0) Hist::Head("hHNl1_lrvar")->fillH2D(abs_rig, lrvar, list->weight);

    if (signr > 0) Hist::Head("hHPl1_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
    if (signr < 0) Hist::Head("hHNl1_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
    
    if (g4mc != nullptr && signr > 0) Hist::Head("hHPl1_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
    if (g4mc != nullptr && signr < 0) Hist::Head("hHNl1_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
    if (g4mc != nullptr && signr > 0) Hist::Head("hHPl1_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
    if (g4mc != nullptr && signr < 0) Hist::Head("hHNl1_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
    
    if (signr > 0) Hist::Head("hHPl1_cnt")->fillH1D(abs_rig, list->weight);
    if (signr < 0) Hist::Head("hHNl1_cnt")->fillH1D(abs_rig, list->weight);
    
    if (g4mc != nullptr) Hist::Head("hHl1_cnt_MC_FLUX10")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux10 * list->weight);
    if (g4mc != nullptr) Hist::Head("hHl1_cnt_MC_FLUX27")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux27 * list->weight);
    
    varsL1.run = list->run;
    varsL1.ev  = list->event;
    varsL1.pt[0] = hyc->geom_status[0]; 
    varsL1.pt[1] = hyc->geom_status[1]; 
    varsL1.pt[2] = hyc->geom_status[2]; 
    varsL1.pt[3] = hyc->geom_status[3]; 
    varsL1.w10 = mc_flux10; 
    varsL1.w27 = mc_flux27; 
    varsL1.rig = rig; 
    varsL1.chix[0] = lchix_in;
    varsL1.chix[1] = lchix_l1;
    varsL1.chiy[0] = lchiy_in;
    varsL1.chiy[1] = lchiy_l1;
    varsL1.rvar[0] = lrvar;
    varsL1.norm[0] = nrm_lx;
    varsL1.norm[1] = nrm_ly;
    varsL1.norm[2] = std::log(nrm_rho / nrm_tau);
    varsL1.fill();
    
    return true;
}

bool Analyzer::process_data_hl9() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_flux10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_flux27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_flux10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_flux27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    // Tracker-L9
    if (!hyc->geom_status[0]) return false;
    if (!hyc->geom_status[2]) return false;
    if (trk->lay[8] != 3 || trk->lay[1] == 0) return false;
    if (trk->num_inn_x < 4) return false;
    if (trk->num_inn_y < 5) return false;

    double rig     = hyc->geom_top_rig[2];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
        
    // Cutoff
    if (CheckType(Type::ISS) && cfr < 0.75) return false;
    
    if (signr > 0) Hist::Head("hHPl9_llr")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
    if (signr < 0) Hist::Head("hHNl9_llr")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
    
    // TRD
    if (trd->tdLLR_ep < 0.75) return false;
    
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);

    double lchix_l9 = std::log(hyc->geom_nchi_x[2]);
    double lchiy_l9 = std::log(hyc->geom_nchi_y[2]);
    
    const double lchi_cut = 4.00;
    if (lchix_in > lchi_cut || lchix_l9 > lchi_cut) return false;
    if (lchiy_in > lchi_cut || lchiy_l9 > lchi_cut) return false;

    double lchix = (1.0/5.0) * lchix_in + (4.0/5.0) * lchix_l9;
    double lchiy = (1.0/5.0) * lchiy_in + (4.0/5.0) * lchiy_l9;
    
    double lrvar = std::log(1.0 + std::abs(hyc->geom_top_rig[2] / hyc->geom_top_rig[0] - 1.0));
    
    double nrm_lx  = hyc->geom_max_norm_lx[2];
    double nrm_ly  = hyc->geom_max_norm_ly[2];
    double nrm_tau = hyc->geom_max_norm_tau[2];
    double nrm_rho = hyc->geom_max_norm_rho[2];
    
    if (signr > 0) Hist::Head("hHPl9_lx")->fillH2D(abs_rig, nrm_lx, list->weight);
    if (signr < 0) Hist::Head("hHNl9_ly")->fillH2D(abs_rig, nrm_ly, list->weight);
    if (signr > 0) Hist::Head("hHPl9_tau")->fillH2D(abs_rig, nrm_tau, list->weight);
    if (signr < 0) Hist::Head("hHNl9_rho")->fillH2D(abs_rig, nrm_rho, list->weight);
    
    if (signr > 0) Hist::Head("hHPl9_lchix")->fillH2D(abs_rig, lchix, list->weight);
    if (signr < 0) Hist::Head("hHNl9_lchix")->fillH2D(abs_rig, lchix, list->weight);
    if (signr > 0) Hist::Head("hHPl9_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr < 0) Hist::Head("hHNl9_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr > 0) Hist::Head("hHPl9_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
    if (signr < 0) Hist::Head("hHNl9_lrvar")->fillH2D(abs_rig, lrvar, list->weight);

    if (signr > 0) Hist::Head("hHPl9_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
    if (signr < 0) Hist::Head("hHNl9_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
    
    if (g4mc != nullptr && signr > 0) Hist::Head("hHPl9_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
    if (g4mc != nullptr && signr < 0) Hist::Head("hHNl9_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
    if (g4mc != nullptr && signr > 0) Hist::Head("hHPl9_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
    if (g4mc != nullptr && signr < 0) Hist::Head("hHNl9_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
    
    if (signr > 0) Hist::Head("hHPl9_cnt")->fillH1D(abs_rig, list->weight);
    if (signr < 0) Hist::Head("hHNl9_cnt")->fillH1D(abs_rig, list->weight);
    
    if (g4mc != nullptr) Hist::Head("hHl9_cnt_MC_FLUX10")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux10 * list->weight);
    if (g4mc != nullptr) Hist::Head("hHl9_cnt_MC_FLUX27")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux27 * list->weight);
    
    varsL9.run = list->run;
    varsL9.ev  = list->event;
    varsL9.pt[0] = hyc->geom_status[0]; 
    varsL9.pt[1] = hyc->geom_status[1]; 
    varsL9.pt[2] = hyc->geom_status[2]; 
    varsL9.pt[3] = hyc->geom_status[3]; 
    varsL9.w10 = mc_flux10; 
    varsL9.w27 = mc_flux27; 
    varsL9.rig = rig; 
    varsL9.chix[0] = lchix_in;
    varsL9.chix[1] = lchix_l9;
    varsL9.chiy[0] = lchiy_in;
    varsL9.chiy[1] = lchiy_l9;
    varsL9.rvar[0] = lrvar;
    varsL9.norm[0] = nrm_lx;
    varsL9.norm[1] = nrm_ly;
    varsL9.norm[2] = std::log(nrm_rho / nrm_tau);
    varsL9.fill();

    return true;
}

bool Analyzer::process_data_hfs() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_flux10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_flux27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        mc_flux10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_flux27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
    }
    
    // Trigger
    if ((trg->bit&8) != 8) return false;
    
    // ECAL
    if (ecal->status && ecal->mvaBDT > -0.6) return false;
    
    // Tracker-FS
    if (!hyc->geom_status[0]) return false;
    if (!hyc->geom_status[1]) return false;
    if (!hyc->geom_status[2]) return false;
    if (!hyc->geom_status[3]) return false;
    if (trk->lay[0] != 3 || trk->lay[8] != 3) return false;
    if (trk->num_inn_x < 4) return false;
    if (trk->num_inn_y < 5) return false;

    double rig     = hyc->geom_top_rig[3];
    double abs_rig = std::abs(rig);
    short  signr   = (rig > 0) ? 1 : -1;
    
    //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
    double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
    double cfr     = CheckType(Type::ISS) ? (abs_rig / cf_rig) : 0.0;
        
    // Cutoff
    if (CheckType(Type::ISS) && cfr < 0.75) return false;
    
    if (signr > 0) Hist::Head("hHPfs_llr")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
    if (signr < 0) Hist::Head("hHNfs_llr")->fillH2D(abs_rig, trd->tdLLR_ep, list->weight);
    
    // TRD
    if (trd->tdLLR_ep < 0.75) return false;
    
    double lchix_in = std::log(hyc->geom_nchi_x[0]);
    double lchiy_in = std::log(hyc->geom_nchi_y[0]);

    double lchix_l1 = std::log(hyc->geom_nchi_x[1]);
    double lchiy_l1 = std::log(hyc->geom_nchi_y[1]);

    double lchix_l9 = std::log(hyc->geom_nchi_x[2]);
    double lchiy_l9 = std::log(hyc->geom_nchi_y[2]);

    double lchix_fs = std::log(hyc->geom_nchi_x[3]);
    double lchiy_fs = std::log(hyc->geom_nchi_y[3]);
    
    const double lchi_cut = 4.00;
    if (lchix_in > lchi_cut || lchix_l1 > lchi_cut || lchix_l9 > lchi_cut || lchix_fs > lchi_cut) return false;
    if (lchiy_in > lchi_cut || lchiy_l1 > lchi_cut || lchiy_l9 > lchi_cut || lchiy_fs > lchi_cut) return false;

    double lchix = (9.0/61.0) * lchix_l1 + (16.0/61.0) * lchix_l9 + (36.0/61.0) * lchix_fs;
    double lchiy = (9.0/61.0) * lchiy_l1 + (16.0/61.0) * lchiy_l9 + (36.0/61.0) * lchiy_fs;
        
    double lrvar_l1 = std::log(1.0 + std::abs(hyc->geom_top_rig[3] / hyc->geom_top_rig[1] - 1.0));
    double lrvar_l9 = std::log(1.0 + std::abs(hyc->geom_top_rig[3] / hyc->geom_top_rig[2] - 1.0));
    double lrvar_fs = std::log(1.0 + std::abs(hyc->geom_top_rig[3] / hyc->geom_top_rig[2] - hyc->geom_top_rig[3] / hyc->geom_top_rig[1]));

    double lrvar = (9.0/61.0) * lrvar_l1 + (16.0/61.0) * lrvar_l9 + (36.0/61.0) * lrvar_fs;
    
    double nrm_lx  = hyc->geom_max_norm_lx[3];
    double nrm_ly  = hyc->geom_max_norm_ly[3];
    double nrm_tau = hyc->geom_max_norm_tau[3];
    double nrm_rho = hyc->geom_max_norm_rho[3];
    
    if (signr > 0) Hist::Head("hHPfs_lx")->fillH2D(abs_rig, nrm_lx, list->weight);
    if (signr < 0) Hist::Head("hHNfs_ly")->fillH2D(abs_rig, nrm_ly, list->weight);
    if (signr > 0) Hist::Head("hHPfs_tau")->fillH2D(abs_rig, nrm_tau, list->weight);
    if (signr < 0) Hist::Head("hHNfs_rho")->fillH2D(abs_rig, nrm_rho, list->weight);
    
    if (signr > 0) Hist::Head("hHPfs_lchix")->fillH2D(abs_rig, lchix, list->weight);
    if (signr < 0) Hist::Head("hHNfs_lchix")->fillH2D(abs_rig, lchix, list->weight);
    if (signr > 0) Hist::Head("hHPfs_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr < 0) Hist::Head("hHNfs_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
    if (signr > 0) Hist::Head("hHPfs_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
    if (signr < 0) Hist::Head("hHNfs_lrvar")->fillH2D(abs_rig, lrvar, list->weight);

    if (signr > 0) Hist::Head("hHPfs_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
    if (signr < 0) Hist::Head("hHNfs_lchiy_lrvar")->fillH3D(abs_rig, lchiy, lrvar, list->weight);
    
    if (g4mc != nullptr && signr > 0) Hist::Head("hHPfs_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
    if (g4mc != nullptr && signr < 0) Hist::Head("hHNfs_lchiy_lrvar_MC_FLUX10")->fillH3D(abs_rig, lchiy, lrvar, mc_flux10 * list->weight);
    if (g4mc != nullptr && signr > 0) Hist::Head("hHPfs_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
    if (g4mc != nullptr && signr < 0) Hist::Head("hHNfs_lchiy_lrvar_MC_FLUX27")->fillH3D(abs_rig, lchiy, lrvar, mc_flux27 * list->weight);
    
    if (signr > 0) Hist::Head("hHPfs_cnt")->fillH1D(abs_rig, list->weight);
    if (signr < 0) Hist::Head("hHNfs_cnt")->fillH1D(abs_rig, list->weight);
    
    if (g4mc != nullptr) Hist::Head("hHfs_cnt_MC_FLUX10")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux10 * list->weight);
    if (g4mc != nullptr) Hist::Head("hHfs_cnt_MC_FLUX27")->fillH1D(std::abs(g4mc->prm_mom/g4mc->prm_chrg), mc_flux27 * list->weight);
    
    varsFS.run = list->run;
    varsFS.ev  = list->event;
    varsFS.pt[0] = hyc->geom_status[0]; 
    varsFS.pt[1] = hyc->geom_status[1]; 
    varsFS.pt[2] = hyc->geom_status[2]; 
    varsFS.pt[3] = hyc->geom_status[3]; 
    varsFS.w10 = mc_flux10; 
    varsFS.w27 = mc_flux27; 
    varsFS.rig = rig; 
    varsFS.chix[0] = lchix_in;
    varsFS.chix[1] = lchix_l1;
    varsFS.chix[2] = lchix_l9;
    varsFS.chix[3] = lchix_fs;
    varsFS.chiy[0] = lchiy_in;
    varsFS.chiy[1] = lchiy_l1;
    varsFS.chiy[2] = lchiy_l9;
    varsFS.chiy[3] = lchiy_fs;
    varsFS.rvar[0] = lrvar_l1;
    varsFS.rvar[1] = lrvar_l9;
    varsFS.rvar[2] = lrvar_fs;
    varsFS.norm[0] = nrm_lx;
    varsFS.norm[1] = nrm_ly;
    varsFS.norm[2] = std::log(nrm_rho / nrm_tau);
    varsFS.fill();

    //double BDTD = modelFS_BDTD->EvaluateMVA("BDTD");
    //double BDTG = modelFS_BDTG->EvaluateMVA("BDTG");

    return true;
}




























/*
bool Analyzer::process_data_h() {
    TrSys::PartType type(TrSys::PartList::kProton);
    
    double mc_flux10 = (g4mc != nullptr) ? std::pow(g4mc->prm_mom/100.0, -1.7) : 1.0;
    double mc_flux27 = (g4mc != nullptr) ? 1.0 : 1.0;
    if (g4mc != nullptr) {
        //double crr_flux_pr = gPROFpr27->Eval(g4mc->prm_mom, 0, "");
        //if (crr_flux_pr < 0) crr_flux_pr = 0.0;

        //double crr_flux_ap = crr_flux_pr * gPROFap2pr->Eval(g4mc->prm_mom, 0, "");
        //if (crr_flux_ap < 0) crr_flux_ap = 0.0;

        //mc_flux10 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap) * (list->antimatter_sw_trigger ? 5.0 : 1.0);
        //mc_flux27 *= (g4mc->prm_chrg > 0 ? crr_flux_pr : crr_flux_ap) * (list->antimatter_sw_trigger ? 5.0 : 1.0);
        
        mc_flux10 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
        mc_flux27 *= (list->antimatter_sw_trigger ? 5.0 : 1.0);
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
    const double nchiy_cut = 2.50;
    //const double nchiy_cut = 2.0;

    while (swl1) {
        double rig     = hyc->geom_top_rig[1];
        double abs_rig = std::abs(rig);
        short  signr   = (rig > 0) ? 1 : -1;
    
        //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
        double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
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

        //if (lchix_in > nchix_cut) break;
        //if (lchix_l1 > nchix_cut) break;

        if (signr > 0) Hist::Head("hHPl1_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPl1_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lchiy")->fillH2D(abs_rig, lchiy, list->weight);

        if (signr > 0) Hist::Head("hHPl1_dchiy_l1")->fillH2D(abs_rig, dchiy_l1, list->weight);
        if (signr < 0) Hist::Head("hHNl1_dchiy_l1")->fillH2D(abs_rig, dchiy_l1, list->weight);
        
        //if (lchiy_in > nchiy_cut) break;
        //if (lchiy_l1 > nchiy_cut) break;

        if (signr > 0) Hist::Head("hHPl1_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl1_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        
        double lrvar_cut = std::sqrt(0.05*0.05 + 1.0 * std::log(1.0+abs_rig/100.0) * std::log(1.0+abs_rig/100.0));
        //if (lrvar > 1.0) break;
        
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
    
        //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
        double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
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

        //if (lchix_in > nchix_cut) break;
        //if (lchix_l9 > nchix_cut) break;

        if (signr > 0) Hist::Head("hHPl9_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_in")->fillH2D(abs_rig, lchiy_in, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        
        if (signr > 0) Hist::Head("hHPl9_dchiy_l9")->fillH2D(abs_rig, dchiy_l9, list->weight);
        if (signr < 0) Hist::Head("hHNl9_dchiy_l9")->fillH2D(abs_rig, dchiy_l9, list->weight);
        
        //if (lchiy_in > nchiy_cut) break;
        //if (lchiy_l9 > nchiy_cut) break;
        
        if (signr > 0) Hist::Head("hHPl9_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        if (signr < 0) Hist::Head("hHNl9_lrvar")->fillH2D(abs_rig, lrvar, list->weight);
        
        double lrvar_cut = std::sqrt(0.06*0.06 + 1.0 * std::log(1.0+abs_rig/80.0) * std::log(1.0+abs_rig/80.0));
        //if (lrvar > 1.0) break;
        
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
    
        //double cf_rig  = CheckType(Type::ISS) ? rti->max_IGRF : 0.0;
        double cf_rig  = CheckType(Type::ISS) ? hyc->max_IGRF : 0.0;
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
        
        //if (lchix_l1 > nchix_cut) break;
        //if (lchix_l9 > nchix_cut) break;
        //if (lchix_fs > nchix_cut) break;
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_l1")->fillH2D(abs_rig, lchiy_l1, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_l9")->fillH2D(abs_rig, lchiy_l9, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy_fs")->fillH2D(abs_rig, lchiy_fs, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy_fs")->fillH2D(abs_rig, lchiy_fs, list->weight);
        
        if (signr > 0) Hist::Head("hHPfs_lchiy")->fillH2D(abs_rig, lchiy, list->weight);
        if (signr < 0) Hist::Head("hHNfs_lchiy")->fillH2D(abs_rig, lchiy, list->weight);

        //if (lchiy_l1 > nchiy_cut) break;
        //if (lchiy_l9 > nchiy_cut) break;
        //if (lchiy_fs > nchiy_cut) break;
        
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

        ana_tree->Fill();

        break;
    }

    return true;
}
*/

#endif // __Analyzer_C__
