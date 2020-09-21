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
    
    std::cout << Format("\n----==== DB ====----\n");
    LOG(INFO) << Format("\n----==== DB ====----\n");
    process_data(varsDB.get_chain(), &Analyzer::process_data_db); 
    
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
  
    Axis AXz("Z",   30, 0.5, 30.5);
    Axis AXq("Q", 3100, 0.0, 31.0);
    Axis AXex("(dE/dx)/Z^2", 2000, 0.0, 8.0);
    Axis AXlchi("Log(#chi^{2}/NDF)", 400, -4.0, 8.0);

    
    //Axis AXz("Z",   10, 0.5, 10.5);
    //Axis AXq("Q", 1000, 0.0, 5.0);
    //Axis AXex("(dE/dx)/Z^2", 1000, 0.0, 4.0);
    
    
    Hist::New("hQtf", HistAxis(AXq));
    Hist::New("hQtk", HistAxis(AXq));
    
    Hist::New("hQtf_lchi_c", HistAxis(AXq, AXlchi));
    Hist::New("hQtf_lchi_t", HistAxis(AXq, AXlchi));
    
    Hist::New("hQtkQtf", HistAxis(AXq, AXq));
    
    Hist::New("hZtfQtk", HistAxis(AXz, AXq));
    Hist::New("hZtkQtf", HistAxis(AXz, AXq));
    
    Hist::New("hZtfEXtk", HistAxis(AXz, AXex));
    Hist::New("hZtkEXtf", HistAxis(AXz, AXex));
    
    Hist::New("hZtfEXtk_x", HistAxis(AXz, AXex));
    Hist::New("hZtfEXtk_y", HistAxis(AXz, AXex));
    
    Hist::New("hQtf_new", HistAxis(AXq));
    Hist::New("hQtk_new", HistAxis(AXq));
    
    Hist::New("hQtf_n2", HistAxis(AXq));
    Hist::New("hQtf_n3", HistAxis(AXq));
    
    Hist::New("hQtf_new_n2", HistAxis(AXq));
    Hist::New("hQtf_new_n3", HistAxis(AXq));
    Hist::New("hQtf_new_n4", HistAxis(AXq));
    
    Hist::New("hQtk_new_n1", HistAxis(AXq));
    Hist::New("hQtk_new_n2", HistAxis(AXq));
    Hist::New("hQtk_new_n3", HistAxis(AXq));

    Hist::New("hQtf_lchi_new", HistAxis(AXq, AXlchi));
    Hist::New("hQtk_lchi_new", HistAxis(AXq, AXlchi));
    
    Hist::New("hQtf_lchi_new_n2", HistAxis(AXq, AXlchi));
    Hist::New("hQtf_lchi_new_n3", HistAxis(AXq, AXlchi));
    Hist::New("hQtf_lchi_new_n4", HistAxis(AXq, AXlchi));
    
    Hist::New("hQtk_lchi_new_n1", HistAxis(AXq, AXlchi));
    Hist::New("hQtk_lchi_new_n2", HistAxis(AXq, AXlchi));
    Hist::New("hQtk_lchi_new_n3", HistAxis(AXq, AXlchi));
    
    Hist::New("hQtkQtf_new", HistAxis(AXq, AXq));
    
    Hist::New("hQtf_good", HistAxis(AXq));
    Hist::New("hQtk_good", HistAxis(AXq));
    Hist::New("hQtf_good_new", HistAxis(AXq));
    Hist::New("hQtk_good_new", HistAxis(AXq));
    
    Hist::New("hQtf_good_n2", HistAxis(AXq));
    Hist::New("hQtf_good_n3", HistAxis(AXq));
    
    Hist::New("hQtf_good_new_n2", HistAxis(AXq));
    Hist::New("hQtf_good_new_n3", HistAxis(AXq));
    Hist::New("hQtf_good_new_n4", HistAxis(AXq));
    
    Hist::New("hQtk_good_new_n1", HistAxis(AXq));
    Hist::New("hQtk_good_new_n2", HistAxis(AXq));
    Hist::New("hQtk_good_new_n3", HistAxis(AXq));
    
    Hist::New("hQtfQtf", HistAxis(AXq, AXq));
    Hist::New("hQtkQtk", HistAxis(AXq, AXq));
    
    Hist::New("hZ1Qtf", HistAxis(AXq));
    Hist::New("hZ1Qtf_new", HistAxis(AXq));
    
    Hist::New("hZ1Qtf_good", HistAxis(AXq));
    Hist::New("hZ1Qtf_good_new", HistAxis(AXq));
    
    Hist::New("hZ1Qtk", HistAxis(AXq));
    Hist::New("hZ1Qtk_new", HistAxis(AXq));
    
    Hist::New("hZ1Qtk_good", HistAxis(AXq));
    Hist::New("hZ1Qtk_good_new", HistAxis(AXq));

    return true;
}


bool Analyzer::process_data_db() {
    VarsDB* vars = &varsDB;
 
    float     Qtf = vars->tfQ;
    float     Ztf = std::lrint(Qtf);
    bool  goodZtf = std::abs(Qtf - Ztf) < 0.2;
    std::vector<float> tfEX;
    for (int it = 0; it < 4; ++it) {
        if (!vars->tfL[it]) continue;
        tfEX.push_back(vars->tfLQ[it] * vars->tfLQ[it]);
    }
    
    float    Qtk = vars->tkQ;
    float    Ztk = std::lrint(Qtk);
    bool goodZtk = std::abs(Qtk - Ztk) < 0.2;
    std::vector<float> tkEXx;
    std::vector<float> tkEXy;
    std::vector<float> tkEXxy;
    for (int it = 1; it < 8; ++it) {
        if (!vars->tkL[it]) continue;
        if (vars->tkLQx[it] > 0) tkEXx .push_back(vars->tkLQx[it] * vars->tkLQx[it]);
        if (vars->tkLQy[it] > 0) tkEXy .push_back(vars->tkLQy[it] * vars->tkLQy[it]);
        if (vars->tkLQx[it] > 0) tkEXxy.push_back(vars->tkLQx[it] * vars->tkLQx[it]);
        if (vars->tkLQy[it] > 0) tkEXxy.push_back(vars->tkLQy[it] * vars->tkLQy[it]);
    }
        
    Hist::Head("hQtkQtf")->fillH2D(Qtk, Qtf, vars->wgt);

    Ztf     = std::lrint(vars->new_tfQ);
    goodZtf = (vars->new_tfQ_nhit == 4) && (vars->new_tfQ_nchi < 2) && (std::abs(vars->new_tfQ - Ztf) < 0.3);

    Hist::Head("hQtf")->fillH1D(Qtf, vars->wgt);
    Hist::Head("hQtf_lchi_c")->fillH2D(Qtf, vars->tfQ_nchi_c, vars->wgt);
    Hist::Head("hQtf_lchi_t")->fillH2D(Qtf, vars->tfQ_nchi_t, vars->wgt);
    if (goodZtf) {
        Hist::Head("hZtfQtk")->fillH2D(Ztf, Qtk, vars->wgt);
        for (auto&& dedx : tkEXxy) {
            Hist::Head("hZtfEXtk")->fillH2D(Ztf, dedx / (Ztf * Ztf), vars->wgt);
        }
        for (auto&& dedx : tkEXx) {
            Hist::Head("hZtfEXtk_x")->fillH2D(Ztf, dedx / (Ztf * Ztf), vars->wgt);
        }
        for (auto&& dedx : tkEXy) {
            Hist::Head("hZtfEXtk_y")->fillH2D(Ztf, dedx / (Ztf * Ztf), vars->wgt);
        }
    }

    Hist::Head("hQtk")->fillH1D(Qtk, vars->wgt);
    if (goodZtk) {
        Hist::Head("hZtkQtf")->fillH2D(Ztk, Qtf, vars->wgt);
        for (auto&& dedx : tfEX) {
            Hist::Head("hZtkEXtf")->fillH2D(Ztk, dedx / (Ztk * Ztk), vars->wgt);
        }
    }
    
    Hist::Head("hQtf_new")->fillH1D(vars->new_tfQ, vars->wgt);
    Hist::Head("hQtk_new")->fillH1D(vars->new_tkQ, vars->wgt);
    
    if (vars->tfN == 2) Hist::Head("hQtf_n2")->fillH1D(Qtf, vars->wgt);
    if (vars->tfN == 3) Hist::Head("hQtf_n3")->fillH1D(Qtf, vars->wgt);
    
    if (vars->new_tfQ_nhit == 2) Hist::Head("hQtf_new_n2")->fillH1D(vars->new_tfQ, vars->wgt);
    if (vars->new_tfQ_nhit == 3) Hist::Head("hQtf_new_n3")->fillH1D(vars->new_tfQ, vars->wgt);
    if (vars->new_tfQ_nhit == 4) Hist::Head("hQtf_new_n4")->fillH1D(vars->new_tfQ, vars->wgt);
    
    if (vars->new_tkQ_ncls == 1) Hist::Head("hQtk_new_n1")->fillH1D(vars->new_tkQ, vars->wgt);
    if (vars->new_tkQ_ncls == 2) Hist::Head("hQtk_new_n2")->fillH1D(vars->new_tkQ, vars->wgt);
    if (vars->new_tkQ_ncls >= 3) Hist::Head("hQtk_new_n3")->fillH1D(vars->new_tkQ, vars->wgt);
    
    Hist::Head("hQtkQtf_new")->fillH2D(vars->new_tkQ, vars->new_tfQ, vars->wgt);
    
    if (vars->new_tfQ > 0) Hist::Head("hQtf_lchi_new")->fillH2D(vars->new_tfQ, vars->new_tfQ_nchi, vars->wgt);
    if (vars->new_tkQ > 0) Hist::Head("hQtk_lchi_new")->fillH2D(vars->new_tkQ, vars->new_tkQ_nchi, vars->wgt);
    
    if (vars->new_tfQ_nhit == 2) Hist::Head("hQtf_lchi_new_n2")->fillH2D(vars->new_tfQ, vars->new_tfQ_nchi, vars->wgt);
    if (vars->new_tfQ_nhit == 3) Hist::Head("hQtf_lchi_new_n3")->fillH2D(vars->new_tfQ, vars->new_tfQ_nchi, vars->wgt);
    if (vars->new_tfQ_nhit == 4) Hist::Head("hQtf_lchi_new_n4")->fillH2D(vars->new_tfQ, vars->new_tfQ_nchi, vars->wgt);
    
    if (vars->new_tkQ_ncls == 1) Hist::Head("hQtk_lchi_new_n1")->fillH2D(vars->new_tkQ, vars->new_tkQ_nchi, vars->wgt);
    if (vars->new_tkQ_ncls == 2) Hist::Head("hQtk_lchi_new_n2")->fillH2D(vars->new_tkQ, vars->new_tkQ_nchi, vars->wgt);
    if (vars->new_tkQ_ncls >= 3) Hist::Head("hQtk_lchi_new_n3")->fillH2D(vars->new_tkQ, vars->new_tkQ_nchi, vars->wgt);
    
    if (vars->new_tfQ_nchi < 2) Hist::Head("hQtf_good")->fillH1D(Qtf, vars->wgt);
    if (vars->new_tkQ_nchi < 2) Hist::Head("hQtk_good")->fillH1D(Qtk, vars->wgt);
    if (vars->new_tfQ_nchi < 2) Hist::Head("hQtf_good_new")->fillH1D(vars->new_tfQ, vars->wgt);
    if (vars->new_tkQ_nchi < 2) Hist::Head("hQtk_good_new")->fillH1D(vars->new_tkQ, vars->wgt);
    
    if (vars->tfN == 2 && vars->new_tfQ_nchi < 2) Hist::Head("hQtf_good_n2")->fillH1D(Qtf, vars->wgt);
    if (vars->tfN == 3 && vars->new_tfQ_nchi < 2) Hist::Head("hQtf_good_n3")->fillH1D(Qtf, vars->wgt);
    
    if (vars->new_tfQ_nhit == 2 && vars->new_tfQ_nchi < 2) Hist::Head("hQtf_good_new_n2")->fillH1D(vars->new_tfQ, vars->wgt);
    if (vars->new_tfQ_nhit == 3 && vars->new_tfQ_nchi < 2) Hist::Head("hQtf_good_new_n3")->fillH1D(vars->new_tfQ, vars->wgt);
    if (vars->new_tfQ_nhit == 4 && vars->new_tfQ_nchi < 2) Hist::Head("hQtf_good_new_n4")->fillH1D(vars->new_tfQ, vars->wgt);
    
    if (vars->new_tkQ_ncls == 1 && vars->new_tkQ_nchi < 2) Hist::Head("hQtk_good_new_n1")->fillH1D(vars->new_tkQ, vars->wgt);
    if (vars->new_tkQ_ncls == 2 && vars->new_tkQ_nchi < 2) Hist::Head("hQtk_good_new_n2")->fillH1D(vars->new_tkQ, vars->wgt);
    if (vars->new_tkQ_ncls >= 3 && vars->new_tkQ_nchi < 2) Hist::Head("hQtk_good_new_n3")->fillH1D(vars->new_tkQ, vars->wgt);
    
    Hist::Head("hQtfQtf")->fillH2D(Qtf, vars->new_tfQ, vars->wgt);
    Hist::Head("hQtkQtk")->fillH2D(Qtk, vars->new_tkQ, vars->wgt);
   

    bool tkZ1 = (Qtk > 0.8 && Qtk < 1.4);

    if (tkZ1 && vars->tfN == 3 && vars->new_tfQ_nhit == 4) Hist::Head("hZ1Qtf")->fillH1D(Qtf);
    if (tkZ1 && vars->tfN == 3 && vars->new_tfQ_nhit == 4) Hist::Head("hZ1Qtf_new")->fillH1D(vars->new_tfQ);
    
    if (tkZ1 && vars->tfN == 3 && vars->new_tfQ_nhit == 4 && vars->new_tfQ_nchi < 2) Hist::Head("hZ1Qtf_good")->fillH1D(Qtf);
    if (tkZ1 && vars->tfN == 3 && vars->new_tfQ_nhit == 4 && vars->new_tfQ_nchi < 2) Hist::Head("hZ1Qtf_good_new")->fillH1D(vars->new_tfQ);
    
    bool tfZ1 = (Qtf > 0.8 && Qtf < 1.4);
    
    if (tfZ1) Hist::Head("hZ1Qtk")->fillH1D(Qtk, vars->wgt);
    if (tfZ1) Hist::Head("hZ1Qtk_new")->fillH1D(vars->new_tkQ, vars->wgt);
    
    if (tfZ1 && vars->new_tkQ_nchi < 2) Hist::Head("hZ1Qtk_good")->fillH1D(Qtk, vars->wgt);
    if (tfZ1 && vars->new_tkQ_nchi < 2) Hist::Head("hZ1Qtk_good_new")->fillH1D(vars->new_tkQ, vars->wgt);
    
    return true;
}

#endif // __Analyzer_C__
