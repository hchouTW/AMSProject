#ifndef __Selector_C__
#define __Selector_C__

#include "selector.h"

void Selector::set_environment() {
    std::cout << Format("\n====  Set Environment ====\n");
    LOG(INFO) << Format("\n====  Set Environment ====\n");

	if (!CheckType(Type::ISS)) AMSSetupR::SlowControlR::ReadFromExternalFile = false;
	if (CheckType(Type::MC))   AMSSetupR::LoadISSMC = false;
	if (CheckType(Type::ISS))  AMSSetupR::RTI::UseLatest(7); // pass7

    // Enable latest alignment
    if (!CheckType(Type::BT)) TkDBc::UseFinal();

	// Disable overwriting of datacards from file
    TRMCFFKEY_DEF::ReadFromFile = 0;
    TRFITFFKEY_DEF::ReadFromFile = 0;
    TRMCFFKEY.ReadFromFile = 0;
    TRFITFFKEY.ReadFromFile = 0;
	
    if (CheckType(Type::BT)) {
		TrdKCluster::IsReadGlobalAlignment = false;
    }
	if (CheckType(Type::MC)) {
		TrdKCluster::IsReadGlobalAlignment = false;
		TrdKCluster::ForceReadAlignment    = false;
		TrdKCluster::ForceReadCalibration  = false;
		TrdKCluster::ForceReadXePressure   = false;
        TrdKCluster::SetDefaultMCXePressure(900);
	}
	
    if (CheckType(Type::ISS)) {
        RichRingR::reloadRunTag = true;
        RichRingR::useTemperatureCorrections = true;
        RichRingR::loadChargeUniformityCorrection = true;
		RichRingR::setBetaCorrection(RichRingR::fullUniformityCorrection);
	}
	if (CheckType(Type::BT) || CheckType(Type::MC)) {
		RichRingR::setBetaCorrection(RichRingR::noCorrection);
	}
	
    EcalShowerR::enableAutomaticRecoveryOfBrokenBuilds = false;
}


void Selector::process_events() {
    if (amsch == nullptr || amsch->GetEntries() == 0) return;
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
    long loop_entries = amsch->GetEntries();
	long print_rate =  loop_entries / 50;
	for (long ientry = 0; ientry < loop_entries; ++ientry) {
        bool last_one_of_files = ((ientry + 1) == loop_entries);
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
        
        AMSEventR* evt_nxt = (last_one_of_files ? amsch->GetEvent(ientry) : amsch->GetEvent(ientry+1));
        if (evt_nxt == nullptr) continue;
        UInt_t utime_nxt = evt_nxt->UTime();

        AMSEventR* evt_cur = amsch->GetEvent(ientry);
        if (evt_cur == nullptr) continue;
        UInt_t utime_cur = evt_cur->UTime();

        bool last_one_of_sec = ((utime_cur != utime_nxt) || last_one_of_files);
        if (!last_one_of_sec) continue;
        process_init();
        event = evt_cur;

        if (!process_prefix()) continue;
        if (!process_data()) continue;
    
        tree->Fill();
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


bool Selector::process_prefix() {
    // Resolution tuning (by Qi Yan)
	if (CheckType(Type::ISS) || CheckType(Type::BT)) {
        TrLinearEtaDB::SetLinearCluster(); // Enable new Eta uniformity(Z=1-26 and above)
        TRFITFFKEY.Zshift = 2; // Enable the dZ correction
    }
    else if (CheckType(Type::MC)) {
        event->SetDefaultMCTuningParameters();
        TrExtAlignDB::SmearExtAlign(); // MC Smear Ext-Layer
        TRCLFFKEY.UseSensorAlign = 0;
        TRCLFFKEY.ClusterCofGOpt = 1;
        TRFITFFKEY.Zshift = -1; // Disable the dZ correction
	}
    TRFITFFKEY.ErcHeY = 0;
    
    if (CheckType(Type::MC)) {
        MCEventgR* primaryMC = event->GetPrimaryMC();
	    if (primaryMC == nullptr) return false;
    }
    
	//-----------------------------//
	//----  Fast Preselection  ----//
	//-----------------------------//
	
    // ~1~ (Based on BetaH(Beta))
    TofRecH::BuildOpt = 0; // normal

    return true;
}

bool Selector::process_data() {
    if (!process_list()) return false;
    if (!process_rti() ) return false;
    
    return true;
}

bool Selector::process_list() {
    data_list.file   = amsch->GetCurrentFile()->GetName();
    data_list.run    = event->Run();
    data_list.event  = event->Event();
    data_list.entry  = amsch->get_tree_entry();
    data_list.utime  = event->UTime();

    return true;
}

bool Selector::process_rti() {
    if (!CheckType(Type::ISS)) return true;
	
    AMSSetupR * setup = AMSSetupR::gethead();
	
    AMSSetupR::RTI rti;
	event->GetRTI(rti);
	
    data_rti.flag     = rti.good;
	data_rti.zenith   = rti.zenith;
	data_rti.livetime = rti.lf * rti.nev / (rti.nev + rti.nerr);

    // good second
	bool good_second = true;
	if ((rti.ntrig/rti.nev) < 0.98 ||
			rti.nerr < 0 || (rti.nerr/rti.nev) > 0.1 ||
			(rti.npart/rti.ntrig) < (0.07/1600*rti.ntrig) || (rti.npart/rti.ntrig) > 0.25 ||
			rti.npart <= 0 || rti.nev > 1800)
		good_second = false;
	else
		good_second = true;
	data_rti.good = good_second;
	
    data_rti.GTOD[0] = rti.r;
	data_rti.GTOD[1] = rti.theta;
	data_rti.GTOD[2] = rti.phi;
	data_rti.GM[0]   = rti.getthetam();
	data_rti.GM[1]   = rti.getphim();
	data_rti.GAT[0]  = rti.glat * TMath::DegToRad();
	data_rti.GAT[1]  = rti.glong * TMath::DegToRad();

    // geomagnetic cuttof
    Float_t mincf_Stoermer[4] = { 0. };
    Float_t maxcf_Stoermer[4] = { 0. };
    Float_t mincf_IGRF[4] = { 0. };
    Float_t maxcf_IGRF[4] = { 0. };
    for (int i = 0; i < 4; i++) {
		mincf_Stoermer[i] = (std::fabs(rti.cf[i][0]) < std::fabs(rti.cf[i][1])) ?
			std::fabs(rti.cf[i][0]) : std::fabs(rti.cf[i][1]);
		maxcf_Stoermer[i] = (std::fabs(rti.cf[i][0]) > std::fabs(rti.cf[i][1])) ?
			std::fabs(rti.cf[i][0]) : std::fabs(rti.cf[i][1]);
		mincf_IGRF[i] = (std::fabs(rti.cfi[i][0]) < std::fabs(rti.cfi[i][1])) ?
			std::fabs(rti.cfi[i][0]) : std::fabs(rti.cfi[i][1]);
		maxcf_IGRF[i] = (std::fabs(rti.cfi[i][0]) > std::fabs(rti.cfi[i][1])) ?
			std::fabs(rti.cfi[i][0]) : std::fabs(rti.cfi[i][1]);

        data_rti.Stoermer[i][0] = std::fabs(rti.cf[i][0]);
        data_rti.Stoermer[i][1] = std::fabs(rti.cf[i][1]);
        data_rti.IGRF[i][0] = std::fabs(rti.cfi[i][0]);
        data_rti.IGRF[i][1] = std::fabs(rti.cfi[i][1]);
	}
    data_rti.min_Stoermer = (*std::min_element(mincf_Stoermer, mincf_Stoermer+4));
    data_rti.max_Stoermer = (*std::max_element(maxcf_Stoermer, maxcf_Stoermer+4));
    data_rti.min_IGRF = (*std::min_element(mincf_IGRF, mincf_IGRF+4));
    data_rti.max_IGRF = (*std::max_element(maxcf_IGRF, maxcf_IGRF+4));

    // |PG-MD| < 35e-4 (L1), 45e-4 (L9) [cm]
	AMSPoint pn1, pn9, pd1, pd9;
	event->GetRTIdL1L9(0, pn1, pd1, event->UTime(), 60);
	event->GetRTIdL1L9(1, pn9, pd9, event->UTime(), 60);
	data_rti.tk_align[0][0] = pd1.x();
	data_rti.tk_align[0][1] = pd1.y();
	data_rti.tk_align[1][0] = pd9.x();
	data_rti.tk_align[1][1] = pd9.y();
	
    // Inner Tracker Temperature (Sensor A)
	std::vector<float> tempSensA;
    if (!setup->fSlowControl.GetData("Sensor A", event->Run(), 0, tempSensA) && 
	    tempSensA.size() > 0) {
		data_rti.tk_temp = tempSensA.at(0);
	}

    data_rti.is_in_SAA = rti.IsInSAA();

    return true;
}

#endif // __Selector_C__
