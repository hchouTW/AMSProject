#ifndef __Analyzer_C__
#define __Analyzer_C__

#include "analyzer.h"

using namespace MGROOT;

void Analyzer::set_environment() {
    std::cout << Format("\n====  Set Environment ====\n");
    LOG(INFO) << Format("\n====  Set Environment ====\n");
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
    
    std::cout << Format("\n----==== FLX ====----\n");
    LOG(INFO) << Format("\n----==== FLX ====----\n");
    process_data(varsFLX.get_chain(), &Analyzer::process_data_flx); 

    std::cout << Format("\n----==== LTF ====----\n");
    LOG(INFO) << Format("\n----==== LTF ====----\n");
    process_data(varsLTF.get_chain(), &Analyzer::process_data_ltf); 
    
    std::cout << Format("\n----==== LRH ====----\n");
    LOG(INFO) << Format("\n----==== LRH ====----\n");
    process_data(varsLRH.get_chain(), &Analyzer::process_data_lrh); 
    
    std::cout << Format("\n----==== IIN ====----\n");
    LOG(INFO) << Format("\n----==== IIN ====----\n");
    process_data(varsIIN.get_chain(), &Analyzer::process_data_iin); 
    
    std::cout << Format("\n----==== IEX ====----\n");
    LOG(INFO) << Format("\n----==== IEX ====----\n");
    process_data(varsIEX.get_chain(), &Analyzer::process_data_iex); 
    
    std::cout << Format("\n----==== HEX ====----\n");
    LOG(INFO) << Format("\n----==== HEX ====----\n");
    process_data(varsHEX.get_chain(), &Analyzer::process_data_hex); 
    
    std::cout << Format("\n----==== HFS ====----\n");
    LOG(INFO) << Format("\n----==== HFS ====----\n");
    process_data(varsHFS.get_chain(), &Analyzer::process_data_hfs); 
    
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
    file->cd();

    std::vector<double> vtme1M( { // 1-month
        1307750400, 1310083200, 1312416000, 1314748800, 1317081600, 1319414400, 1321747200, 1324080000, 1326412800, 1328745600, 
        1331078400, 1333411200, 1335744000, 1338076800, 1340409600, 1342742400, 1345075200, 1347408000, 1349740800, 1352073600, 
        1354406400, 1356739200, 1359072000, 1361404800, 1363737600, 1366070400, 1368403200, 1370736000, 1373068800, 1375401600, 
        1377734400, 1380067200, 1382400000, 1384732800, 1387065600, 1389398400, 1391731200, 1394064000, 1396396800, 1398729600, 
        1401062400, 1403395200, 1405728000, 1408060800, 1410393600, 1412726400, 1415059200, 1417392000, 1419724800, 1422057600, 
        1424390400, 1426723200, 1429056000, 1431388800, 1433721600, 1436054400, 1438387200, 1440720000, 1443052800, 1445385600, 
        1447718400, 1450051200, 1452384000, 1454716800, 1457049600, 1459382400, 1461715200, 1464048000, 1466380800, 1468713600, 
        1471046400, 1473379200, 1475712000, 1478044800, 1480377600, 1482710400, 1485043200, 1487376000, 1489708800, 1492041600, 
        1494374400, 1496707200, 1499040000, 1501372800, 1503705600, 1506038400, 1508371200, 1510704000, 1513036800, 1515369600, 
        1517702400, 1520035200, 1522368000, 1524700800, 1527033600, 1529366400, 1531699200, 1534032000, 1536364800, 1538697600, 
        1541030400, 1543363200, 1545696000 } );
    std::vector<double> vtme2M( { // 2-month
        1307750400, 1312416000, 1317081600, 1321747200, 1326412800, 
        1331078400, 1335744000, 1340409600, 1345075200, 1349740800,
        1354406400, 1359072000, 1363737600, 1368403200, 1373068800,
        1377734400, 1382400000, 1387065600, 1391731200, 1396396800,
        1401062400, 1405728000, 1410393600, 1415059200, 1419724800,
        1424390400, 1429056000, 1433721600, 1438387200, 1443052800,
        1447718400, 1452384000, 1457049600, 1461715200, 1466380800,
        1471046400, 1475712000, 1480377600, 1485043200, 1489708800,
        1494374400, 1499040000, 1503705600, 1508371200, 1513036800,
        1517702400, 1522368000, 1527033600, 1531699200, 1536364800,
        1541030400, 1545696000 } );
    std::vector<double> vtme3M( { // 3-month
        1307750400, 1314748800, 1321747200, 1328745600, 1335744000, 
        1342742400, 1349740800, 1356739200, 1363737600, 1370736000, 
        1377734400, 1384732800, 1391731200, 1398729600, 1405728000, 
        1412726400, 1419724800, 1426723200, 1433721600, 1440720000, 
        1447718400, 1454716800, 1461715200, 1468713600, 1475712000, 
        1482710400, 1489708800, 1496707200, 1503705600, 1510704000, 
        1517702400, 1524700800, 1531699200, 1538697600, 1545696000 } );
    std::vector<double> vtme6M( { // 6-month
        1307750400, 1321747200, 1335744000, 1349740800, 1363737600, 
        1377734400, 1391731200, 1405728000, 1419724800, 1433721600,
        1447718400, 1461715200, 1475712000, 1489708800, 1503705600, 
        1517702400, 1531699200, 1545696000 } );
    Axis AXtme1M("Time", vtme1M);
    Axis AXtme2M("Time", vtme2M);
    Axis AXtme3M("Time", vtme3M);
    Axis AXtme6M("Time", vtme6M);

    std::vector<double> vrig( {
        1.00,    1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
        3.29,    3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
        8.48,    9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
        19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
        41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
        93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 330.00, 525.00 } );
    
    std::vector<double> vrig_tme_v1( {
        1.00,    1.51,   2.15,   2.67,   3.29,    4.02,   4.88,   5.90,   7.09,   8.48,    
        10.10,  12.00,  14.10,  16.60,  19.50,   22.80,  26.70,  31.10,  36.10,  41.90,  
        48.50,  56.10 } );
    std::vector<double> vrig_tme_v2( {
        1.00,   1.51,   1.92,   2.40,   2.97,   3.64,   4.43,   5.37,   6.47,   7.76,      
        9.26,  11.00,  13.00,  15.30,  18.00,  21.10,  24.70,  28.80,  33.50,  38.90,  
       45.10,  52.20,  60.30  } );
    
    Axis AXrig ("|Rigidity| [GV]", vrig);
    //Axis AXTrig("|Rigidity| [GV]", vrig_tme_v1);
    Axis AXTrig("|Rigidity| [GV]", vrig_tme_v2);
    
    Hist::New("hFlx_lv", HistAxis(AXrig));
    Hist::New("hT1MFlx_lv", HistAxis(AXtme1M, AXTrig));
    Hist::New("hT2MFlx_lv", HistAxis(AXtme2M, AXTrig));
    Hist::New("hT3MFlx_lv", HistAxis(AXtme3M, AXTrig));
    Hist::New("hT6MFlx_lv", HistAxis(AXtme6M, AXTrig));
    
    Hist::New("hFlxP_cnt_alltrg", HistAxis(AXrig));
    Hist::New("hFlxN_cnt_alltrg", HistAxis(AXrig));
    Hist::New("hT1MFlxP_cnt_alltrg", HistAxis(AXtme1M, AXTrig));
    Hist::New("hT1MFlxN_cnt_alltrg", HistAxis(AXtme1M, AXTrig));
    Hist::New("hT2MFlxP_cnt_alltrg", HistAxis(AXtme2M, AXTrig));
    Hist::New("hT2MFlxN_cnt_alltrg", HistAxis(AXtme2M, AXTrig));
    Hist::New("hT3MFlxP_cnt_alltrg", HistAxis(AXtme3M, AXTrig));
    Hist::New("hT3MFlxN_cnt_alltrg", HistAxis(AXtme3M, AXTrig));
    Hist::New("hT6MFlxP_cnt_alltrg", HistAxis(AXtme6M, AXTrig));
    Hist::New("hT6MFlxN_cnt_alltrg", HistAxis(AXtme6M, AXTrig));
    
    Hist::New("hFlxP_cnt", HistAxis(AXrig));
    Hist::New("hFlxN_cnt", HistAxis(AXrig));
    Hist::New("hT1MFlxP_cnt", HistAxis(AXtme1M, AXTrig));
    Hist::New("hT1MFlxN_cnt", HistAxis(AXtme1M, AXTrig));
    Hist::New("hT2MFlxP_cnt", HistAxis(AXtme2M, AXTrig));
    Hist::New("hT2MFlxN_cnt", HistAxis(AXtme2M, AXTrig));
    Hist::New("hT3MFlxP_cnt", HistAxis(AXtme3M, AXTrig));
    Hist::New("hT3MFlxN_cnt", HistAxis(AXtme3M, AXTrig));
    Hist::New("hT6MFlxP_cnt", HistAxis(AXtme6M, AXTrig));
    Hist::New("hT6MFlxN_cnt", HistAxis(AXtme6M, AXTrig));
    
    Hist::New("hFlx_cnt_MC"       , HistAxis(AXrig));
    Hist::New("hFlx_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hFlx_cnt_MC_FLUX27", HistAxis(AXrig));
    
    Hist::New("hTFlx_cnt_MC"       , HistAxis(AXTrig));
    Hist::New("hTFlx_cnt_MC_FLUX10", HistAxis(AXTrig));
    Hist::New("hTFlx_cnt_MC_FLUX27", HistAxis(AXTrig));

    Axis AXLTFsqrm ("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 100, -3.0, 3.0);
    Axis AXTLTFsqrm("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",  75, -3.0, 3.0);
    
    Hist::New("hLtfP_sqrm",    HistAxis(AXrig, AXLTFsqrm));
    Hist::New("hLtfN_sqrm",    HistAxis(AXrig, AXLTFsqrm));
    Hist::New("hLtfP_sqrm_pr", HistAxis(AXrig, AXLTFsqrm));
    Hist::New("hLtfN_sqrm_el", HistAxis(AXrig, AXLTFsqrm));
    
    Hist::New("hLtfP_CF_sqrm",    HistAxis(AXrig, AXLTFsqrm));
    Hist::New("hLtfN_CF_sqrm",    HistAxis(AXrig, AXLTFsqrm));
    Hist::New("hLtfP_CF_sqrm_pr", HistAxis(AXrig, AXLTFsqrm));
    Hist::New("hLtfN_CF_sqrm_el", HistAxis(AXrig, AXLTFsqrm));
    
    Hist::New("hT1MLtfP_sqrm",    HistAxis(AXtme1M, AXTrig, AXTLTFsqrm));
    Hist::New("hT1MLtfN_sqrm",    HistAxis(AXtme1M, AXTrig, AXTLTFsqrm));
    Hist::New("hT1MLtfP_sqrm_pr", HistAxis(AXtme1M, AXTrig, AXTLTFsqrm));
    Hist::New("hT1MLtfN_sqrm_el", HistAxis(AXtme1M, AXTrig, AXTLTFsqrm));
    Hist::New("hT2MLtfP_sqrm",    HistAxis(AXtme2M, AXTrig, AXTLTFsqrm));
    Hist::New("hT2MLtfN_sqrm",    HistAxis(AXtme2M, AXTrig, AXTLTFsqrm));
    Hist::New("hT2MLtfP_sqrm_pr", HistAxis(AXtme2M, AXTrig, AXTLTFsqrm));
    Hist::New("hT2MLtfN_sqrm_el", HistAxis(AXtme2M, AXTrig, AXTLTFsqrm));
    Hist::New("hT3MLtfP_sqrm",    HistAxis(AXtme3M, AXTrig, AXTLTFsqrm));
    Hist::New("hT3MLtfN_sqrm",    HistAxis(AXtme3M, AXTrig, AXTLTFsqrm));
    Hist::New("hT3MLtfP_sqrm_pr", HistAxis(AXtme3M, AXTrig, AXTLTFsqrm));
    Hist::New("hT3MLtfN_sqrm_el", HistAxis(AXtme3M, AXTrig, AXTLTFsqrm));
    Hist::New("hT6MLtfP_sqrm",    HistAxis(AXtme6M, AXTrig, AXTLTFsqrm));
    Hist::New("hT6MLtfN_sqrm",    HistAxis(AXtme6M, AXTrig, AXTLTFsqrm));
    Hist::New("hT6MLtfP_sqrm_pr", HistAxis(AXtme6M, AXTrig, AXTLTFsqrm));
    Hist::New("hT6MLtfN_sqrm_el", HistAxis(AXtme6M, AXTrig, AXTLTFsqrm));
    Hist::New("hRLtfP_sqrm_pr", HistAxis(AXTrig, AXTLTFsqrm));
    Hist::New("hRLtfN_sqrm_el", HistAxis(AXTrig, AXTLTFsqrm));
    
    Hist::New("hLtfP_cnt", HistAxis(AXrig));
    Hist::New("hLtfN_cnt", HistAxis(AXrig));

    Hist::New("hLtf_cnt_MC", HistAxis(AXrig));
    Hist::New("hLtf_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hLtf_cnt_MC_FLUX27", HistAxis(AXrig));
    
    Axis AXLRHsqrm ("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", 125, -1.5, 3.0);
    Axis AXTLRHsqrm("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",  75, -1.5, 3.0);
    
    Hist::New("hLrhP_sqrm",    HistAxis(AXrig, AXLRHsqrm));
    Hist::New("hLrhN_sqrm",    HistAxis(AXrig, AXLRHsqrm));
    Hist::New("hLrhP_sqrm_pr", HistAxis(AXrig, AXLRHsqrm));
    Hist::New("hLrhN_sqrm_el", HistAxis(AXrig, AXLRHsqrm));
    
    Hist::New("hLrhP_CF_sqrm",    HistAxis(AXrig, AXLRHsqrm));
    Hist::New("hLrhN_CF_sqrm",    HistAxis(AXrig, AXLRHsqrm));
    Hist::New("hLrhP_CF_sqrm_pr", HistAxis(AXrig, AXLRHsqrm));
    Hist::New("hLrhN_CF_sqrm_el", HistAxis(AXrig, AXLRHsqrm));

    Hist::New("hT1MLrhP_sqrm",    HistAxis(AXtme1M, AXTrig, AXTLRHsqrm));
    Hist::New("hT1MLrhN_sqrm",    HistAxis(AXtme1M, AXTrig, AXTLRHsqrm));
    Hist::New("hT1MLrhP_sqrm_pr", HistAxis(AXtme1M, AXTrig, AXTLRHsqrm));
    Hist::New("hT1MLrhN_sqrm_el", HistAxis(AXtme1M, AXTrig, AXTLRHsqrm));
    Hist::New("hT2MLrhP_sqrm",    HistAxis(AXtme2M, AXTrig, AXTLRHsqrm));
    Hist::New("hT2MLrhN_sqrm",    HistAxis(AXtme2M, AXTrig, AXTLRHsqrm));
    Hist::New("hT2MLrhP_sqrm_pr", HistAxis(AXtme2M, AXTrig, AXTLRHsqrm));
    Hist::New("hT2MLrhN_sqrm_el", HistAxis(AXtme2M, AXTrig, AXTLRHsqrm));
    Hist::New("hT3MLrhP_sqrm",    HistAxis(AXtme3M, AXTrig, AXTLRHsqrm));
    Hist::New("hT3MLrhN_sqrm",    HistAxis(AXtme3M, AXTrig, AXTLRHsqrm));
    Hist::New("hT3MLrhP_sqrm_pr", HistAxis(AXtme3M, AXTrig, AXTLRHsqrm));
    Hist::New("hT3MLrhN_sqrm_el", HistAxis(AXtme3M, AXTrig, AXTLRHsqrm));
    Hist::New("hT6MLrhP_sqrm",    HistAxis(AXtme6M, AXTrig, AXTLRHsqrm));
    Hist::New("hT6MLrhN_sqrm",    HistAxis(AXtme6M, AXTrig, AXTLRHsqrm));
    Hist::New("hT6MLrhP_sqrm_pr", HistAxis(AXtme6M, AXTrig, AXTLRHsqrm));
    Hist::New("hT6MLrhN_sqrm_el", HistAxis(AXtme6M, AXTrig, AXTLRHsqrm));
    Hist::New("hRLrhP_sqrm_pr", HistAxis(AXTrig, AXTLRHsqrm));
    Hist::New("hRLrhN_sqrm_el", HistAxis(AXTrig, AXTLRHsqrm));

    Hist::New("hLrhP_cnt", HistAxis(AXrig));
    Hist::New("hLrhN_cnt", HistAxis(AXrig));
    
    Hist::New("hLrh_cnt_MC", HistAxis(AXrig));
    Hist::New("hLrh_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hLrh_cnt_MC_FLUX27", HistAxis(AXrig));
    
    Axis AXIINllr ("TRD Estimator", 150, 0.0, 1.6);
    Axis AXTIINllr("TRD Estimator", 100, 0.0, 1.6);
    
    Hist::New("hIinP_llr",    HistAxis(AXrig, AXIINllr));
    Hist::New("hIinN_llr",    HistAxis(AXrig, AXIINllr));
    Hist::New("hIinP_llr_pr", HistAxis(AXrig, AXIINllr));
    Hist::New("hIinN_llr_el", HistAxis(AXrig, AXIINllr));
    
    Hist::New("hIinP_CF_llr",    HistAxis(AXrig, AXIINllr));
    Hist::New("hIinN_CF_llr",    HistAxis(AXrig, AXIINllr));
    Hist::New("hIinP_CF_llr_pr", HistAxis(AXrig, AXIINllr));
    Hist::New("hIinN_CF_llr_el", HistAxis(AXrig, AXIINllr));
    
    Hist::New("hT1MIinP_llr",    HistAxis(AXtme1M, AXTrig, AXTIINllr));
    Hist::New("hT1MIinN_llr",    HistAxis(AXtme1M, AXTrig, AXTIINllr));
    Hist::New("hT1MIinP_llr_pr", HistAxis(AXtme1M, AXTrig, AXTIINllr));
    Hist::New("hT1MIinN_llr_el", HistAxis(AXtme1M, AXTrig, AXTIINllr));
    Hist::New("hT2MIinP_llr",    HistAxis(AXtme2M, AXTrig, AXTIINllr));
    Hist::New("hT2MIinN_llr",    HistAxis(AXtme2M, AXTrig, AXTIINllr));
    Hist::New("hT2MIinP_llr_pr", HistAxis(AXtme2M, AXTrig, AXTIINllr));
    Hist::New("hT2MIinN_llr_el", HistAxis(AXtme2M, AXTrig, AXTIINllr));
    Hist::New("hT3MIinP_llr",    HistAxis(AXtme3M, AXTrig, AXTIINllr));
    Hist::New("hT3MIinN_llr",    HistAxis(AXtme3M, AXTrig, AXTIINllr));
    Hist::New("hT3MIinP_llr_pr", HistAxis(AXtme3M, AXTrig, AXTIINllr));
    Hist::New("hT3MIinN_llr_el", HistAxis(AXtme3M, AXTrig, AXTIINllr));
    Hist::New("hT6MIinP_llr",    HistAxis(AXtme6M, AXTrig, AXTIINllr));
    Hist::New("hT6MIinN_llr",    HistAxis(AXtme6M, AXTrig, AXTIINllr));
    Hist::New("hT6MIinP_llr_pr", HistAxis(AXtme6M, AXTrig, AXTIINllr));
    Hist::New("hT6MIinN_llr_el", HistAxis(AXtme6M, AXTrig, AXTIINllr));
    Hist::New("hRIinP_llr_pr", HistAxis(AXTrig, AXTIINllr));
    Hist::New("hRIinN_llr_el", HistAxis(AXTrig, AXTIINllr));
    
    Hist::New("hIinP_cnt", HistAxis(AXrig));
    Hist::New("hIinN_cnt", HistAxis(AXrig));
    
    Hist::New("hIin_cnt_MC", HistAxis(AXrig));
    Hist::New("hIin_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hIin_cnt_MC_FLUX27", HistAxis(AXrig));
   
    Axis AXIEXllr ("TRD Estimator", 100, 0.0, 1.6);
    Axis AXTIEXllr("TRD Estimator", 100, 0.0, 1.6);
    
    Hist::New("hIexP_llr", HistAxis(AXrig, AXIEXllr));
    Hist::New("hIexN_llr", HistAxis(AXrig, AXIEXllr));
    Hist::New("hIexP_llr_pr", HistAxis(AXrig, AXIEXllr));
    Hist::New("hIexN_llr_el", HistAxis(AXrig, AXIEXllr));
    
    Hist::New("hIexP_CF_llr", HistAxis(AXrig, AXIEXllr));
    Hist::New("hIexN_CF_llr", HistAxis(AXrig, AXIEXllr));
    Hist::New("hIexP_CF_llr_pr", HistAxis(AXrig, AXIEXllr));
    Hist::New("hIexN_CF_llr_el", HistAxis(AXrig, AXIEXllr));
    
    Hist::New("hT1MIexP_llr",    HistAxis(AXtme1M, AXTrig, AXTIEXllr));
    Hist::New("hT1MIexN_llr",    HistAxis(AXtme1M, AXTrig, AXTIEXllr));
    Hist::New("hT1MIexP_llr_pr", HistAxis(AXtme1M, AXTrig, AXTIEXllr));
    Hist::New("hT1MIexN_llr_el", HistAxis(AXtme1M, AXTrig, AXTIEXllr));
    Hist::New("hT2MIexP_llr",    HistAxis(AXtme2M, AXTrig, AXTIEXllr));
    Hist::New("hT2MIexN_llr",    HistAxis(AXtme2M, AXTrig, AXTIEXllr));
    Hist::New("hT2MIexP_llr_pr", HistAxis(AXtme2M, AXTrig, AXTIEXllr));
    Hist::New("hT2MIexN_llr_el", HistAxis(AXtme2M, AXTrig, AXTIEXllr));
    Hist::New("hT3MIexP_llr",    HistAxis(AXtme3M, AXTrig, AXTIEXllr));
    Hist::New("hT3MIexN_llr",    HistAxis(AXtme3M, AXTrig, AXTIEXllr));
    Hist::New("hT3MIexP_llr_pr", HistAxis(AXtme3M, AXTrig, AXTIEXllr));
    Hist::New("hT3MIexN_llr_el", HistAxis(AXtme3M, AXTrig, AXTIEXllr));
    Hist::New("hT6MIexP_llr",    HistAxis(AXtme6M, AXTrig, AXTIEXllr));
    Hist::New("hT6MIexN_llr",    HistAxis(AXtme6M, AXTrig, AXTIEXllr));
    Hist::New("hT6MIexP_llr_pr", HistAxis(AXtme6M, AXTrig, AXTIEXllr));
    Hist::New("hT6MIexN_llr_el", HistAxis(AXtme6M, AXTrig, AXTIEXllr));
    Hist::New("hRIexP_llr_pr", HistAxis(AXTrig, AXTIEXllr));
    Hist::New("hRIexN_llr_el", HistAxis(AXTrig, AXTIEXllr));
    
    Hist::New("hIexP_cnt", HistAxis(AXrig));
    Hist::New("hIexN_cnt", HistAxis(AXrig));
    
    Hist::New("hIex_cnt_MC", HistAxis(AXrig));
    Hist::New("hIex_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hIex_cnt_MC_FLUX27", HistAxis(AXrig));
    
    Axis AXHEXllr("TRD Estimator", 200, 0.0, 1.6);
    Axis AXHEXcc ("CC Estimator", 50, -1, 1);
    Axis AXHEXchi("#chi^{2}", 50, -4, 4);
    
    Hist::New("hHexP_llr", HistAxis(AXrig, AXHEXllr));
    Hist::New("hHexN_llr", HistAxis(AXrig, AXHEXllr));
    
    Hist::New("hHexP_cc", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexN_cc", HistAxis(AXrig, AXHEXcc));
    
    Hist::New("hHexP_cc_MC", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexN_cc_MC", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexP_cc_MC_FLUX10", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexN_cc_MC_FLUX10", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexP_cc_MC_FLUX27", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexN_cc_MC_FLUX27", HistAxis(AXrig, AXHEXcc));
    
    Hist::New("hHexP_CF_cc", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexN_CF_cc", HistAxis(AXrig, AXHEXcc));
    
    Hist::New("hHexP_CF_cc_MC", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexN_CF_cc_MC", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexP_CF_cc_MC_FLUX10", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexN_CF_cc_MC_FLUX10", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexP_CF_cc_MC_FLUX27", HistAxis(AXrig, AXHEXcc));
    Hist::New("hHexN_CF_cc_MC_FLUX27", HistAxis(AXrig, AXHEXcc));
    
    Hist::New("hHexP_cnt", HistAxis(AXrig));
    Hist::New("hHexN_cnt", HistAxis(AXrig));
    
    Hist::New("hHex_cnt_MC", HistAxis(AXrig));
    Hist::New("hHex_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hHex_cnt_MC_FLUX27", HistAxis(AXrig));
    
    Axis AXHFSllr("TRD Estimator", 200, 0.0, 1.6);
    Axis AXHFScc ("CC Estimator", 50, -1, 1);
    Axis AXHFSchi("#chi^{2}", 50, -4, 4);
    
    Hist::New("hHfsP_llr", HistAxis(AXrig, AXHFSllr));
    Hist::New("hHfsN_llr", HistAxis(AXrig, AXHFSllr));

    Hist::New("hHfsP_cc", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsN_cc", HistAxis(AXrig, AXHFScc));
    
    Hist::New("hHfsP_cc_MC", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsN_cc_MC", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsP_cc_MC_FLUX10", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsN_cc_MC_FLUX10", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsP_cc_MC_FLUX27", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsN_cc_MC_FLUX27", HistAxis(AXrig, AXHFScc));
    
    Hist::New("hHfsP_CF_cc", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsN_CF_cc", HistAxis(AXrig, AXHFScc));
    
    Hist::New("hHfsP_CF_cc_MC", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsN_CF_cc_MC", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsP_CF_cc_MC_FLUX10", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsN_CF_cc_MC_FLUX10", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsP_CF_cc_MC_FLUX27", HistAxis(AXrig, AXHFScc));
    Hist::New("hHfsN_CF_cc_MC_FLUX27", HistAxis(AXrig, AXHFScc));
    
    Hist::New("hHfsP_cnt", HistAxis(AXrig));
    Hist::New("hHfsN_cnt", HistAxis(AXrig));
    
    Hist::New("hHfs_cnt_MC", HistAxis(AXrig));
    Hist::New("hHfs_cnt_MC_FLUX10", HistAxis(AXrig));
    Hist::New("hHfs_cnt_MC_FLUX27", HistAxis(AXrig));
    
    return true;
}

bool Analyzer::process_data_flx() {
    VarsFLX* vars = &varsFLX;

    const Axis& AXrig  = Hist::Head("hFlxP_cnt")->xaxis();
    const Axis& AXTrig = Hist::Head("hT1MFlxP_cnt")->xaxis();

    if (vars->tkL1 != 3) return false;
    if (vars->tkL9 != 3) return false;
    if (vars->tkL2 == 0) return false;

    short  sign = vars->sign;
    double arig = std::abs(vars->rig);

    bool overcf = false;
    int icfsec = AXrig.find(1.20 * vars->cfsec);
    if (arig >= AXrig.max() || arig > AXrig(icfsec)) overcf = true;
    
    bool overcfT = false;
    int icfsecT = AXTrig.find(1.20 * vars->cfsec);
    if (arig >= AXTrig.max() || arig > AXTrig(icfsecT)) overcfT = true;

    bool update_ut = (FlxUTime != vars->ut);
    FlxUTime = vars->ut;

    if (update_ut) {
        for (int ir = icfsec + 1; ir <= AXrig.nbin(); ++ir) {
            double cenr = AXrig.center(ir, AxisScale::kLog);
            Hist::Head("hFlx_lv")->fillH1D(cenr, vars->lv);
        }
        for (int ir = icfsecT + 1; ir <= AXTrig.nbin(); ++ir) {
            double cenr = AXTrig.center(ir, AxisScale::kLog);
            Hist::Head("hT1MFlx_lv")->fillH2D(vars->ut, cenr, vars->lv);
            Hist::Head("hT2MFlx_lv")->fillH2D(vars->ut, cenr, vars->lv);
            Hist::Head("hT3MFlx_lv")->fillH2D(vars->ut, cenr, vars->lv);
            Hist::Head("hT6MFlx_lv")->fillH2D(vars->ut, cenr, vars->lv);
        }
    }
   
    if (sign == 0) return false;

    if (vars->mc || overcf) {
    if (sign > 0) Hist::Head("hFlxP_cnt_alltrg")->fillH1D(arig, vars->wgt * (vars->trg>0?1.0:100.0));
    if (sign < 0) Hist::Head("hFlxN_cnt_alltrg")->fillH1D(arig, vars->wgt * (vars->trg>0?1.0:100.0));
    } 
    
    if (vars->mc || overcfT) {
    if (sign > 0) Hist::Head("hT1MFlxP_cnt_alltrg")->fillH2D(vars->ut, arig, vars->wgt * (vars->trg>0?1.0:100.0));
    if (sign < 0) Hist::Head("hT1MFlxN_cnt_alltrg")->fillH2D(vars->ut, arig, vars->wgt * (vars->trg>0?1.0:100.0));
    if (sign > 0) Hist::Head("hT2MFlxP_cnt_alltrg")->fillH2D(vars->ut, arig, vars->wgt * (vars->trg>0?1.0:100.0));
    if (sign < 0) Hist::Head("hT2MFlxN_cnt_alltrg")->fillH2D(vars->ut, arig, vars->wgt * (vars->trg>0?1.0:100.0));
    if (sign > 0) Hist::Head("hT3MFlxP_cnt_alltrg")->fillH2D(vars->ut, arig, vars->wgt * (vars->trg>0?1.0:100.0));
    if (sign < 0) Hist::Head("hT3MFlxN_cnt_alltrg")->fillH2D(vars->ut, arig, vars->wgt * (vars->trg>0?1.0:100.0));
    if (sign > 0) Hist::Head("hT6MFlxP_cnt_alltrg")->fillH2D(vars->ut, arig, vars->wgt * (vars->trg>0?1.0:100.0));
    if (sign < 0) Hist::Head("hT6MFlxN_cnt_alltrg")->fillH2D(vars->ut, arig, vars->wgt * (vars->trg>0?1.0:100.0));
    }

    if (vars->trg == 0) return false;

    if (vars->mc || overcf) {
    if (sign > 0) Hist::Head("hFlxP_cnt")->fillH1D(arig, vars->wgt);
    if (sign < 0) Hist::Head("hFlxN_cnt")->fillH1D(arig, vars->wgt);
    }

    if (vars->mc || overcfT) {
    if (sign > 0) Hist::Head("hT1MFlxP_cnt")->fillH2D(vars->ut, arig, vars->wgt);
    if (sign < 0) Hist::Head("hT1MFlxN_cnt")->fillH2D(vars->ut, arig, vars->wgt);
    if (sign > 0) Hist::Head("hT2MFlxP_cnt")->fillH2D(vars->ut, arig, vars->wgt);
    if (sign < 0) Hist::Head("hT2MFlxN_cnt")->fillH2D(vars->ut, arig, vars->wgt);
    if (sign > 0) Hist::Head("hT3MFlxP_cnt")->fillH2D(vars->ut, arig, vars->wgt);
    if (sign < 0) Hist::Head("hT3MFlxN_cnt")->fillH2D(vars->ut, arig, vars->wgt);
    if (sign > 0) Hist::Head("hT6MFlxP_cnt")->fillH2D(vars->ut, arig, vars->wgt);
    if (sign < 0) Hist::Head("hT6MFlxN_cnt")->fillH2D(vars->ut, arig, vars->wgt);
    }
    
    if (vars->mc) Hist::Head("hFlx_cnt_MC")->fillH1D(std::abs(vars->mc_rig), vars->wgt);
    if (vars->mc) Hist::Head("hFlx_cnt_MC_FLUX10")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w10);
    if (vars->mc) Hist::Head("hFlx_cnt_MC_FLUX27")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w27);
    
    if (vars->mc) Hist::Head("hTFlx_cnt_MC")->fillH1D(std::abs(vars->mc_rig), vars->wgt);
    if (vars->mc) Hist::Head("hTFlx_cnt_MC_FLUX10")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w10);
    if (vars->mc) Hist::Head("hTFlx_cnt_MC_FLUX27")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w27);

    return true;
}

bool Analyzer::process_data_ltf() {
    VarsLTF* vars = &varsLTF;
    
    short  sign = vars->sign;
    double rig  = vars->rig;
    double bta  = vars->bta;
    double arig = std::abs(vars->rig);
    double cfr  = vars->mc ? 0.0 : (arig / vars->cfevt);
    //double cfr  = vars->mc ? 0.0 : (arig / vars->cfsec);

    double mass  = 0.938272297;
    double beta  = 1.0 / std::hypot(1.0, mass / arig);
    double gamma = 1.0 / sqrt((1.0 - beta) * (1.0 + beta));
    double sqrm_cut = 0.5 - 2.0 * std::pow(gamma - 1.0, 6.0);
    if (vars->sqrm < sqrm_cut) return false;
    
    const double lchix_cut = 1.75;
    const double lchiy_cut = 0.75 + 0.5 * std::erfc((arig - 25.0) / 20.0);
    if (vars->lxin > lchix_cut) return false;
    if (vars->lyin > lchiy_cut) return false;

    bool is_like_pr_not_ring = (!vars->rich    && vars->nhinn == 0 && vars->nhout <= 3);
    bool is_like_pr_has_ring = ( vars->rich_pr && vars->nhinn <= 1 && vars->nhout <= 3);
    bool is_like_pr_rich = (is_like_pr_not_ring || is_like_pr_has_ring);
    bool is_like_el_rich = (vars->rich_el && vars->nhinn <= 1 && vars->nhout <= 3);
    bool is_like_pr = (vars->llr > (0.75 * vars->bta) && is_like_pr_rich);
    bool is_like_el = (vars->llr < (0.75 * vars->bta) && is_like_el_rich);
     
    // Cutoff
    double lmtcfr = 1.3;
    double limit_cfr = 0.80 * (1.0 + 1.0 * std::pow((1.0/bta/bta - 1.0), 3.0));
    if (!vars->mc && cfr < limit_cfr) return false;
    
    if (vars->nvtxx[0] >= 3) return false;
    if (vars->nvtxx[1] >= 3) return false;
    if (vars->nvtxy[0] >= 7) return false;
    if (vars->nvtxy[1] >= 7) return false;

    if (sign > 0 && is_like_pr) Hist::Head("hLtfP_sqrm_pr")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hLtfN_sqrm_el")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (sign > 0 && is_like_pr && cfr > lmtcfr) Hist::Head("hLtfP_CF_sqrm_pr")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el && cfr > lmtcfr) Hist::Head("hLtfN_CF_sqrm_el")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (sign > 0 && is_like_pr) Hist::Head("hT1MLtfP_sqrm_pr")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT1MLtfN_sqrm_el")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT2MLtfP_sqrm_pr")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT2MLtfN_sqrm_el")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT3MLtfP_sqrm_pr")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT3MLtfN_sqrm_el")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT6MLtfP_sqrm_pr")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT6MLtfN_sqrm_el")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hRLtfP_sqrm_pr")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hRLtfN_sqrm_el")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (!is_like_pr) return false;
    
    if (sign > 0) Hist::Head("hLtfP_sqrm")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hLtfN_sqrm")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (sign > 0 && cfr > lmtcfr) Hist::Head("hLtfP_CF_sqrm")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0 && cfr > lmtcfr) Hist::Head("hLtfN_CF_sqrm")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (sign > 0) Hist::Head("hT1MLtfP_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hT1MLtfN_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0) Hist::Head("hT2MLtfP_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hT2MLtfN_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0) Hist::Head("hT3MLtfP_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hT3MLtfN_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0) Hist::Head("hT6MLtfP_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hT6MLtfN_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    
    if (sign > 0) Hist::Head("hLtfP_cnt")->fillH1D(arig, vars->wgt);
    if (sign < 0) Hist::Head("hLtfN_cnt")->fillH1D(arig, vars->wgt);
    
    if (vars->mc) Hist::Head("hLtf_cnt_MC")->fillH1D(std::abs(vars->mc_rig), vars->wgt);
    if (vars->mc) Hist::Head("hLtf_cnt_MC_FLUX10")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w10);
    if (vars->mc) Hist::Head("hLtf_cnt_MC_FLUX27")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w27);

    return true;
}


bool Analyzer::process_data_lrh() {
    VarsLRH* vars = &varsLRH;

    short  sign = vars->sign;
    double rig  = vars->rig;
    double bta  = vars->bta;
    double arig = std::abs(vars->rig);
    double cfr  = vars->mc ? 0.0 : (arig / vars->cfevt);
    //double cfr  = vars->mc ? 0.0 : (arig / vars->cfsec);
    
    const double lchix_cut = 1.75;
    const double lchiy_cut = 0.75 + 0.5 * std::erfc((arig - 25.0) / 20.0);
    if (vars->lxin > lchix_cut) return false;
    if (vars->lyin > lchiy_cut) return false;
    
    bool is_like_pr = (vars->llr > 0.75);
    bool is_like_el = (vars->llr < 0.65);

    // Cutoff
    double lmtcfr = 1.3;
    double limit_cfr = 0.80 * (1.0 + 0.8 * std::pow((1.0/bta/bta - 1.0), 3.0));
    if (!vars->mc && cfr < limit_cfr) return false;
    
    if (sign > 0 && is_like_pr) Hist::Head("hLrhP_sqrm_pr")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hLrhN_sqrm_el")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (sign > 0 && is_like_pr && cfr > lmtcfr) Hist::Head("hLrhP_CF_sqrm_pr")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el && cfr > lmtcfr) Hist::Head("hLrhN_CF_sqrm_el")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (sign > 0 && is_like_pr) Hist::Head("hT1MLrhP_sqrm_pr")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT1MLrhN_sqrm_el")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT2MLrhP_sqrm_pr")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT2MLrhN_sqrm_el")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT3MLrhP_sqrm_pr")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT3MLrhN_sqrm_el")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT6MLrhP_sqrm_pr")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT6MLrhN_sqrm_el")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hRLrhP_sqrm_pr")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hRLrhN_sqrm_el")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (!is_like_pr) return false;
    
    if (sign > 0) Hist::Head("hLrhP_sqrm")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hLrhN_sqrm")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (sign > 0 && cfr > lmtcfr) Hist::Head("hLrhP_CF_sqrm")->fillH2D(arig, vars->sqrm, vars->wgt);
    if (sign < 0 && cfr > lmtcfr) Hist::Head("hLrhN_CF_sqrm")->fillH2D(arig, vars->sqrm, vars->wgt);
    
    if (sign > 0) Hist::Head("hT1MLrhP_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hT1MLrhN_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0) Hist::Head("hT2MLrhP_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hT2MLrhN_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0) Hist::Head("hT3MLrhP_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hT3MLrhN_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign > 0) Hist::Head("hT6MLrhP_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
    if (sign < 0) Hist::Head("hT6MLrhN_sqrm")->fillH3D(vars->ut, arig, vars->sqrm, vars->wgt);
  
    if (sign > 0) Hist::Head("hLrhP_cnt")->fillH1D(arig, vars->wgt);
    if (sign < 0) Hist::Head("hLrhN_cnt")->fillH1D(arig, vars->wgt);
    
    if (vars->mc) Hist::Head("hLrh_cnt_MC")->fillH1D(std::abs(vars->mc_rig), vars->wgt);
    if (vars->mc) Hist::Head("hLrh_cnt_MC_FLUX10")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w10);
    if (vars->mc) Hist::Head("hLrh_cnt_MC_FLUX27")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w27);
    
    return true;
}

bool Analyzer::process_data_iin() {
    VarsIIN* vars = &varsIIN;

    short  sign = vars->sign;
    double rig  = vars->rig;
    double arig = std::abs(vars->rig);
    double cfr  = vars->mc ? 0.0 : (arig / vars->cfevt);
    //double cfr  = vars->mc ? 0.0 : (arig / vars->cfsec);
    
    const double lchix_cut = 1.75;
    const double lchiy_cut = 0.75 + 0.5 * std::erfc((arig - 25.0) / 20.0);
    if (vars->lxin > lchix_cut) return false;
    if (vars->lyin > lchiy_cut) return false;

    bool is_like_pr = (vars->ecal && vars->mvaBDT < -0.6 && std::abs(vars->engD/vars->rig) < 0.4);
    bool is_like_el = (vars->ecal && vars->mvaBDT >  0.2);
    
    // Cutoff
    double lmtcfr = 1.3;
    double limit_cfr = 0.80;
    if (!vars->mc && cfr < limit_cfr) return false;
    
    if (sign > 0 && is_like_pr) Hist::Head("hIinP_llr_pr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hIinN_llr_el")->fillH2D(arig, vars->llr, vars->wgt);
    
    if (sign > 0 && is_like_pr && cfr > lmtcfr) Hist::Head("hIinP_CF_llr_pr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el && cfr > lmtcfr) Hist::Head("hIinN_CF_llr_el")->fillH2D(arig, vars->llr, vars->wgt);
    
    if (sign > 0 && is_like_pr) Hist::Head("hT1MIinP_llr_pr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT1MIinN_llr_el")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT2MIinP_llr_pr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT2MIinN_llr_el")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT3MIinP_llr_pr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT3MIinN_llr_el")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT6MIinP_llr_pr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT6MIinN_llr_el")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hRIinP_llr_pr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hRIinN_llr_el")->fillH2D(arig, vars->llr, vars->wgt);
    
    // ECAL
    if (vars->ecal && (vars->mvaBDT > -0.6 || std::abs(vars->engD/vars->rig) > 0.4)) return false;
   
    if (sign > 0) Hist::Head("hIinP_llr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hIinN_llr")->fillH2D(arig, vars->llr, vars->wgt);
    
    if (sign > 0 && cfr > lmtcfr) Hist::Head("hIinP_CF_llr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0 && cfr > lmtcfr) Hist::Head("hIinN_CF_llr")->fillH2D(arig, vars->llr, vars->wgt);
    
    if (sign > 0) Hist::Head("hT1MIinP_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hT1MIinN_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0) Hist::Head("hT2MIinP_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hT2MIinN_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0) Hist::Head("hT3MIinP_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hT3MIinN_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0) Hist::Head("hT6MIinP_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hT6MIinN_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    
    if (sign > 0) Hist::Head("hIinP_cnt")->fillH1D(arig, vars->wgt);
    if (sign < 0) Hist::Head("hIinN_cnt")->fillH1D(arig, vars->wgt);

    if (vars->mc) Hist::Head("hIin_cnt_MC")->fillH1D(std::abs(vars->mc_rig), vars->wgt);
    if (vars->mc) Hist::Head("hIin_cnt_MC_FLUX10")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w10);
    if (vars->mc) Hist::Head("hIin_cnt_MC_FLUX27")->fillH1D(std::abs(vars->mc_rig), vars->wgt * vars->mc_w27);
    
    return true;
}


bool Analyzer::process_data_iex() {
    VarsIEX* vars = &varsIEX;

    short  sign = vars->sign;
    double rig  = vars->rig;
    double arig = std::abs(vars->rig);
    double cfr  = vars->mc ? 0.0 : (arig / vars->cfevt);
    //double cfr  = vars->mc ? 0.0 : (arig / vars->cfsec);
    
    // Tracker
    const double lchi_cut  = 1.75;
    const double lchix_cut = 1.75;
    const double lchiy_cut = 0.50 + 0.625 * std::erfc((arig - 45.0) / 20.0);
    if (vars->lxin > lchi_cut) return false;
    if (vars->lyin > lchi_cut) return false;
    if (vars->lxex > lchix_cut) return false;
    if (vars->lyex > lchiy_cut) return false;

    // Cutoff
    double lmtcfr = 1.3;
    double limit_cfr = 0.80;
    if (!vars->mc && cfr < limit_cfr) return false;
    
    bool is_like_pr = (vars->ecal && vars->mvaBDT < -0.6 && std::abs(vars->engD/vars->rig) < 0.4);
    bool is_like_el = (vars->ecal && vars->mvaBDT >  0.2);

    if (sign > 0 && is_like_pr) Hist::Head("hIexP_llr_pr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hIexN_llr_el")->fillH2D(arig, vars->llr, vars->wgt);
    
    if (sign > 0 && is_like_pr && cfr > lmtcfr) Hist::Head("hIexP_CF_llr_pr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el && cfr > lmtcfr) Hist::Head("hIexN_CF_llr_el")->fillH2D(arig, vars->llr, vars->wgt);
    
    if (sign > 0 && is_like_pr) Hist::Head("hT1MIexP_llr_pr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT1MIexN_llr_el")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT2MIexP_llr_pr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT2MIexN_llr_el")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT3MIexP_llr_pr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT3MIexN_llr_el")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hT6MIexP_llr_pr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hT6MIexN_llr_el")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0 && is_like_pr) Hist::Head("hRIexP_llr_pr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0 && is_like_el) Hist::Head("hRIexN_llr_el")->fillH2D(arig, vars->llr, vars->wgt);
    
    // ECAL
    if (vars->ecal && (vars->mvaBDT > -0.6 || std::abs(vars->engD/vars->rig) > 0.4)) return false;
    
    if (sign > 0) Hist::Head("hIexP_llr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hIexN_llr")->fillH2D(arig, vars->llr, vars->wgt);
    
    if (sign > 0 && cfr > lmtcfr) Hist::Head("hIexP_CF_llr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0 && cfr > lmtcfr) Hist::Head("hIexN_CF_llr")->fillH2D(arig, vars->llr, vars->wgt);
    
    if (sign > 0) Hist::Head("hT1MIexP_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hT1MIexN_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0) Hist::Head("hT2MIexP_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hT2MIexN_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0) Hist::Head("hT3MIexP_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hT3MIexN_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign > 0) Hist::Head("hT6MIexP_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hT6MIexN_llr")->fillH3D(vars->ut, arig, vars->llr, vars->wgt);
    
    if (sign > 0) Hist::Head("hIexP_cnt")->fillH1D(arig, vars->wgt);
    if (sign < 0) Hist::Head("hIexN_cnt")->fillH1D(arig, vars->wgt);
    
    if (vars->mc) Hist::Head("hIex_cnt_MC")->fillH1D(std::abs(vars->mc_rig), vars->wgt);
    if (vars->mc) Hist::Head("hIex_cnt_MC_FLUX10")->fillH1D(std::abs(vars->mc_rig), vars->mc_w10 * vars->wgt);
    if (vars->mc) Hist::Head("hIex_cnt_MC_FLUX27")->fillH1D(std::abs(vars->mc_rig), vars->mc_w27 * vars->wgt);
    
    return true;
}
    

bool Analyzer::process_data_hex() {
    VarsHEX* vars = &varsHEX;

    short  sign = vars->sign;
    double rig  = vars->rig;
    double arig = std::abs(vars->rig);
    double cfr  = vars->mc ? 0.0 : (arig / vars->cfevt);
    //double cfr  = vars->mc ? 0.0 : (arig / vars->cfsec);

    double cc = vars->get_tmva_value();
    
    // ECAL
    if (vars->ecal && (vars->mvaBDT > -0.6 || std::abs(vars->engD/vars->rig) > 0.4)) return false;
    
    // Tracker
    const double lchi_cut  = 2.25;
    if (vars->lxin > lchi_cut) return false;
    if (vars->lyin > lchi_cut) return false;
    if (vars->lxex > lchi_cut) return false;
    if (vars->lyex > lchi_cut) return false;

    const double lchiC_cut  = 3.0;
    if (vars->lxinC > lchiC_cut) return false;
    if (vars->lyinC > lchiC_cut) return false;
    if (vars->lxexC > lchiC_cut) return false;
    if (vars->lyexC > lchiC_cut) return false;
    
    // Cutoff
    double lmtcfr = 1.3;
    double limit_cfr = 0.80;
    if (!vars->mc && cfr < limit_cfr) return false;
   
    if (sign > 0) Hist::Head("hHexP_llr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hHexN_llr")->fillH2D(arig, vars->llr, vars->wgt);
    
    // TRD
    double llr_cut = 0.75 - 0.1 * std::erfc(3.0 * std::log(350.0 / arig));
    if (vars->llr < llr_cut) return false;
    
    if (sign > 0) Hist::Head("hHexP_cc")->fillH2D(arig, cc, vars->wgt);
    if (sign < 0) Hist::Head("hHexN_cc")->fillH2D(arig, cc, vars->wgt);
    
    if (vars->mc && sign > 0) Hist::Head("hHexP_cc_MC")->fillH2D(arig, cc, vars->wgt);
    if (vars->mc && sign < 0) Hist::Head("hHexN_cc_MC")->fillH2D(arig, cc, vars->wgt);
    if (vars->mc && sign > 0) Hist::Head("hHexP_cc_MC_FLUX10")->fillH2D(arig, cc, vars->mc_w10 * vars->wgt);
    if (vars->mc && sign < 0) Hist::Head("hHexN_cc_MC_FLUX10")->fillH2D(arig, cc, vars->mc_w10 * vars->wgt);
    if (vars->mc && sign > 0) Hist::Head("hHexP_cc_MC_FLUX27")->fillH2D(arig, cc, vars->mc_w27 * vars->wgt);
    if (vars->mc && sign < 0) Hist::Head("hHexN_cc_MC_FLUX27")->fillH2D(arig, cc, vars->mc_w27 * vars->wgt);
    
    if (sign > 0 && cfr > lmtcfr) Hist::Head("hHexP_CF_cc")->fillH2D(arig, cc, vars->wgt);
    if (sign < 0 && cfr > lmtcfr) Hist::Head("hHexN_CF_cc")->fillH2D(arig, cc, vars->wgt);
   
    if (vars->mc && sign > 0 && cfr > lmtcfr) Hist::Head("hHexP_CF_cc_MC")->fillH2D(arig, cc, vars->wgt);
    if (vars->mc && sign < 0 && cfr > lmtcfr) Hist::Head("hHexN_CF_cc_MC")->fillH2D(arig, cc, vars->wgt);
    if (vars->mc && sign > 0 && cfr > lmtcfr) Hist::Head("hHexP_CF_cc_MC_FLUX10")->fillH2D(arig, cc, vars->mc_w10 * vars->wgt);
    if (vars->mc && sign < 0 && cfr > lmtcfr) Hist::Head("hHexN_CF_cc_MC_FLUX10")->fillH2D(arig, cc, vars->mc_w10 * vars->wgt);
    if (vars->mc && sign > 0 && cfr > lmtcfr) Hist::Head("hHexP_CF_cc_MC_FLUX27")->fillH2D(arig, cc, vars->mc_w27 * vars->wgt);
    if (vars->mc && sign < 0 && cfr > lmtcfr) Hist::Head("hHexN_CF_cc_MC_FLUX27")->fillH2D(arig, cc, vars->mc_w27 * vars->wgt);
    
    if (sign > 0) Hist::Head("hHexP_cnt")->fillH1D(arig, vars->wgt);
    if (sign < 0) Hist::Head("hHexN_cnt")->fillH1D(arig, vars->wgt);
    
    if (vars->mc) Hist::Head("hHex_cnt_MC")->fillH1D(std::abs(vars->mc_rig), vars->wgt);
    if (vars->mc) Hist::Head("hHex_cnt_MC_FLUX10")->fillH1D(std::abs(vars->mc_rig), vars->mc_w10 * vars->wgt);
    if (vars->mc) Hist::Head("hHex_cnt_MC_FLUX27")->fillH1D(std::abs(vars->mc_rig), vars->mc_w27 * vars->wgt);
    
    return true;
}


bool Analyzer::process_data_hfs() {
    VarsHFS* vars = &varsHFS;

    short  sign = vars->sign;
    double rig  = vars->rig;
    double arig = std::abs(vars->rig);
    double cfr  = vars->mc ? 0.0 : (arig / vars->cfevt);
    //double cfr  = vars->mc ? 0.0 : (arig / vars->cfsec);

    double cc = vars->get_tmva_value();
   
    // ECAL
    if (vars->ecal && (vars->mvaBDT > -0.6 || std::abs(vars->engD/vars->rig) > 0.4)) return false;
    
    // Tracker
    const double lchi_cut  = 2.25;
    if (vars->lxin > lchi_cut) return false;
    if (vars->lyin > lchi_cut) return false;
    if (vars->lxl1 > lchi_cut) return false;
    if (vars->lyl1 > lchi_cut) return false;
    if (vars->lxl9 > lchi_cut) return false;
    if (vars->lyl9 > lchi_cut) return false;
    if (vars->lxfs > lchi_cut) return false;
    if (vars->lyfs > lchi_cut) return false;
    
    const double lchiC_cut  = 3.0;
    if (vars->lxinC > lchiC_cut) return false;
    if (vars->lyinC > lchiC_cut) return false;
    if (vars->lxl1C > lchiC_cut) return false;
    if (vars->lyl1C > lchiC_cut) return false;
    if (vars->lxl9C > lchiC_cut) return false;
    if (vars->lyl9C > lchiC_cut) return false;
    if (vars->lxfsC > lchiC_cut) return false;
    if (vars->lyfsC > lchiC_cut) return false;
    
    // Cutoff
    double lmtcfr = 1.3;
    double limit_cfr = 0.80;
    if (!vars->mc && cfr < limit_cfr) return false;
    
    if (sign > 0) Hist::Head("hHfsP_llr")->fillH2D(arig, vars->llr, vars->wgt);
    if (sign < 0) Hist::Head("hHfsN_llr")->fillH2D(arig, vars->llr, vars->wgt);
    
    // TRD
    double llr_cut = 0.75 - 0.1 * std::erfc(3.0 * std::log(350.0 / arig));
    if (vars->llr < llr_cut) return false;
    //if (vars->llr < 0.75) return false;
     
    if (sign > 0) Hist::Head("hHfsP_cc")->fillH2D(arig, cc, vars->wgt);
    if (sign < 0) Hist::Head("hHfsN_cc")->fillH2D(arig, cc, vars->wgt);
    
    if (vars->mc && sign > 0) Hist::Head("hHfsP_cc_MC")->fillH2D(arig, cc, vars->wgt);
    if (vars->mc && sign < 0) Hist::Head("hHfsN_cc_MC")->fillH2D(arig, cc, vars->wgt);
    if (vars->mc && sign > 0) Hist::Head("hHfsP_cc_MC_FLUX10")->fillH2D(arig, cc, vars->mc_w10 * vars->wgt);
    if (vars->mc && sign < 0) Hist::Head("hHfsN_cc_MC_FLUX10")->fillH2D(arig, cc, vars->mc_w10 * vars->wgt);
    if (vars->mc && sign > 0) Hist::Head("hHfsP_cc_MC_FLUX27")->fillH2D(arig, cc, vars->mc_w27 * vars->wgt);
    if (vars->mc && sign < 0) Hist::Head("hHfsN_cc_MC_FLUX27")->fillH2D(arig, cc, vars->mc_w27 * vars->wgt);
    
    if (sign > 0 && cfr > lmtcfr) Hist::Head("hHfsP_CF_cc")->fillH2D(arig, cc, vars->wgt);
    if (sign < 0 && cfr > lmtcfr) Hist::Head("hHfsN_CF_cc")->fillH2D(arig, cc, vars->wgt);
    
    if (vars->mc && sign > 0 && cfr > lmtcfr) Hist::Head("hHfsP_CF_cc_MC")->fillH2D(arig, cc, vars->wgt);
    if (vars->mc && sign < 0 && cfr > lmtcfr) Hist::Head("hHfsN_CF_cc_MC")->fillH2D(arig, cc, vars->wgt);
    if (vars->mc && sign > 0 && cfr > lmtcfr) Hist::Head("hHfsP_CF_cc_MC_FLUX10")->fillH2D(arig, cc, vars->mc_w10 * vars->wgt);
    if (vars->mc && sign < 0 && cfr > lmtcfr) Hist::Head("hHfsN_CF_cc_MC_FLUX10")->fillH2D(arig, cc, vars->mc_w10 * vars->wgt);
    if (vars->mc && sign > 0 && cfr > lmtcfr) Hist::Head("hHfsP_CF_cc_MC_FLUX27")->fillH2D(arig, cc, vars->mc_w27 * vars->wgt);
    if (vars->mc && sign < 0 && cfr > lmtcfr) Hist::Head("hHfsN_CF_cc_MC_FLUX27")->fillH2D(arig, cc, vars->mc_w27 * vars->wgt);
    
    if (sign > 0) Hist::Head("hHfsP_cnt")->fillH1D(arig, vars->wgt);
    if (sign < 0) Hist::Head("hHfsN_cnt")->fillH1D(arig, vars->wgt);
    
    if (vars->mc) Hist::Head("hHfs_cnt_MC")->fillH1D(std::abs(vars->mc_rig), vars->wgt);
    if (vars->mc) Hist::Head("hHfs_cnt_MC_FLUX10")->fillH1D(std::abs(vars->mc_rig), vars->mc_w10 * vars->wgt);
    if (vars->mc) Hist::Head("hHfs_cnt_MC_FLUX27")->fillH1D(std::abs(vars->mc_rig), vars->mc_w27 * vars->wgt);
   
    return true;
}


#endif // __Analyzer_C__
