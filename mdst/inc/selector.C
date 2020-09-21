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
    
    TrSys::PhysEnv::ReadMagAMS("/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    TrSys::PhysEnv::ReadMatAMS("/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    
    if (!TrSys::PhysEnv::IMagStatus()) TrSys::PhysEnv::ReadMagAMS("/eos/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    if (!TrSys::PhysEnv::IMatStatus()) TrSys::PhysEnv::ReadMatAMS("/eos/user/h/hchou/ExternalLibs/DB/material");
    
    if (!TrSys::PhysEnv::IMagStatus()) TrSys::PhysEnv::ReadMagAMS("/afs/cern.ch/work/h/hchou/public/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    if (!TrSys::PhysEnv::IMatStatus()) TrSys::PhysEnv::ReadMatAMS("/afs/cern.ch/work/h/hchou/public/ExternalLibs/DB/material");

    LOG_IF(ERROR, !TrSys::PhysEnv::IMagStatus()) << Format("NO MAG DATABASE.");
    LOG_IF(ERROR, !TrSys::PhysEnv::IMatStatus()) << Format("NO MAT DATABASE.");
    if (!TrSys::PhysEnv::IMagStatus() || !TrSys::PhysEnv::IMatStatus()) exit(-1);
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
        
        AMSEventR* evt = amsch->GetEvent(ientry);
        if (evt == nullptr) continue;
        //if (ientry>500) break; // testcode

        process_init();
        event = evt;

        utime_pre = utime_cur;
        utime_cur = event->UTime();

        if (!process_prefix()) continue;

        if (!process_data()) continue;
        //if (!process_sel()) continue;
        //if (!process_prd()) continue;
    
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

	// ~2~ (Based on TrTrack)
    //if (event->NTrTrack() != 1) return false;
    if (event->NTrTrack() == 0) return false;
    
	// ~3~ (Based on TrdTrack)
    //if (event->NTrdTrack() != 1) return false;
    //if (event->NTrdTrack() == 0) return false;

	// ~4~ (Based on Particle)
	ParticleR   * partSIG = (event->NParticle() > 0) ? event->pParticle(0) : nullptr;
	TrTrackR    * trtkSIG = (partSIG != nullptr) ? partSIG->pTrTrack() : nullptr;
	BetaHR      * btahSIG = (partSIG != nullptr) ? partSIG->pBetaH() : nullptr;
	//EcalShowerR * ecalSIG = (partSIG != nullptr) ? partSIG->pEcalShower() : nullptr;
	if (partSIG == nullptr) return false;
	if (trtkSIG == nullptr || btahSIG == nullptr) return false;
    //if (ecalSIG == nullptr) return false;
    
	// ~5~ (Based on BetaH)
	if (event->NBetaH() == 0) return false;
	if (btahSIG->GetBetaPattern() != 4444) return false;

    double betah = btahSIG->GetBeta();
    if (betah <= 0.30) return false; // keep down-going
    if (betah <= 0.95) return false; // keep down-going with high velocity

	// ~6~ (Based on Track Hits)
	const unsigned short TrPtL2  =   2; //  2
	const unsigned short TrPtL1  =   1; //  1
	const unsigned short TrPtL9  = 256; //  256
	const unsigned short TrPtL34 =  12; //  4 +   8
	const unsigned short TrPtL56 =  48; // 16 +  32
	const unsigned short TrPtL78 = 192; // 64 + 128
	const unsigned short TrPtL19 = 257; //  1 + 256
	unsigned short trBitPattJ   = trtkSIG->GetBitPatternJ();
	unsigned short trBitPattXYJ = trtkSIG->GetBitPatternXYJ();
	
	bool isTrInner   = ((trBitPattJ&TrPtL34) > 0 && 
	                    (trBitPattJ&TrPtL56) > 0 && 
					    (trBitPattJ&TrPtL78) > 0);
	bool isTrInnerXY = ((trBitPattXYJ&TrPtL34) > 0 && 
	                    (trBitPattXYJ&TrPtL56) > 0 && 
					    (trBitPattXYJ&TrPtL78) > 0);
	if (!isTrInner) return false;
	if (!isTrInnerXY) return false;
	
    bool hasTrL2   = ((trBitPattJ&TrPtL2)  > 0); 
    bool hasTrL1o9 = ((trBitPattJ&TrPtL19) > 0); 
    if (!hasTrL2 && !hasTrL1o9) return false;
	
	Int_t numOfTrInX = 0;
	Int_t numOfTrInY = 0;
	for (Int_t ilay = 2; ilay <= 7; ilay++) {
		bool hasY  = ((trBitPattJ&(1<<ilay)) > 0);
		bool hasXY = ((trBitPattXYJ&(1<<ilay)) > 0);
		if (hasY)  numOfTrInY++;
		if (hasXY) numOfTrInX++;
	}
	if (numOfTrInX <= 2 || numOfTrInY <= 3) return false;

    // Set Reference Pointer
    pParticle   = event->pParticle(0);
    pBetaH      = pParticle->pBetaH()     ;
    pTrTrack    = pParticle->pTrTrack()   ;
    pTrdTrack   = pParticle->pTrdTrack()  ;
    pTrdHTrack  = pParticle->pTrdHTrack() ;
    pEcalShower = pParticle->pEcalShower();
    pRichRing   = pParticle->pRichRing()  ;
	
    // Beta Information
	float beta = pBetaH->GetBeta();
    
    // TRD Track Information
    //if (pTrdTrack == nullptr) return false;

    // Tracker Track Information
    if (pTrTrack == nullptr) return false;
    int idMax = pTrTrack->iTrTrackPar(1, 0, 23); // Rebuild coordinate align
    if (idMax < 0) pTrTrack->iTrTrackPar(1, 3, 23); // Rebuild coordinate align
    
    const int refit = 22;
    int   tkID = pTrTrack->iTrTrackPar(1, 3, refit);
    //float qin  = pTrTrack->GetInnerQH(2, beta, tkID);
    float QYJrms; int QYJpatt;
    float qin = pTrTrack->GetInnerQYJ(QYJrms, QYJpatt, 2, beta, tkID);
    if (tkID < 0 || qin <= 0) return false;

    int   zin  = (qin <= 1.0) ? 1 : std::lrint(qin);
    float mass = (zin <    2) ? TrFit::Mproton : (0.5 * (TrFit::Mhelium) * zin);
    tkID = pTrTrack->iTrTrackPar(1, 3, refit, mass, zin);
    if (tkID < 0) return false;
		
    float ckin_rig = pTrTrack->GetRigidity(tkID, 1); // z = 0
    //if (std::abs(ckin_rig) < 10) return false; // testcode

    data_zin  = zin;
    data_mass = mass;
    data_beta = beta;

    // preselection in particle charge
    //if (data_zin >= 3) return false;
    
    // Keep Q=1 particles
    //if (data_zin != 1) return false;
    //if (CheckType(Type::ISS) && data_zin != 1) return false;

	// ECAL Information
	// pre-selection (ECAL)
	const float lmtrECAL = 7; // 7 cm
	if (pEcalShower != nullptr) {
		AMSPoint EcalPnt(pEcalShower->CofG);
		AMSPoint pnt; AMSDir dir;
		pTrTrack->Interpolate(EcalPnt.z(), pnt, dir, tkID);
		float drPnt = std::abs(pnt.x() - EcalPnt.x());
        dist_tk_ecal = drPnt;
		if (drPnt > lmtrECAL) { pEcalShower = nullptr; }
        //if (dist_tk_ecal > lmtrECAL) return false;
	}

	// TRD Information
	// pre-selection (TRD)
	const float lmtrTRD = 7; // 7 cm
	if (pTrdTrack != nullptr) {
		AMSPoint TrdPnt(pTrdTrack->Coo);
		AMSPoint pnt; AMSDir dir;
		pTrTrack->Interpolate(TrdPnt.z(), pnt, dir, tkID);
		float drPnt = std::abs(pnt.x() - TrdPnt.x());
        dist_tk_td = drPnt;
		if (drPnt > lmtrTRD) { pTrdTrack = nullptr; }
        //if (dist_tk_td > lmtrTRD) return false;
	}

	if (pTrdHTrack != nullptr) {
		AMSPoint TrdHPnt(pTrdHTrack->Coo);
		AMSPoint pnt; AMSDir dir;
		pTrTrack->Interpolate(TrdHPnt.z(), pnt, dir, tkID);
		float drPnt = std::abs(pnt.x() - TrdHPnt.x());
        dist_tk_tdh = drPnt;
		if (drPnt > lmtrTRD) { pTrdHTrack == nullptr; }
        //if (dist_tk_tdh > lmtrTRD) return false;
	}

    return true;
}

bool Selector::process_data() {
    if (!process_list()) return false;
    if (!process_g4mc()) return false;
    if (!process_rti() ) return false;
    if (!process_trg() ) return false;
    if (!process_acc() ) return false;
    if (!process_tof() ) return false;
    if (!process_trk() ) return false;
    if (!process_trd() ) return false;
    if (!process_ecal()) return false;
   
    /*
    //if (!process_rich()) return false;
    //if (!process_trk_kf() ) return false;
    //if (!process_hyc() ) return false;
    
    TrSys::Part org_part(TrSys::PartList::kElectron);
    //TrSys::Part org_part(TrSys::PartList::kDeuterium);
    org_part.set_location (0, 0, 158);
    org_part.set_direction(0, 0, -1);
    
    TrSys::Part st1(org_part);
    st1.set_rig(5.00);
    std::cerr << Form("BEFORE MOM %8.4f BETA %8.4f\n", st1.mom(), st1.bta());
    auto&& rlt1 = TrSys::Prop::DoPropToZ(55, st1, TrSys::PropArgs(1, 1, false, true));
    //auto&& rlt1 = TrSys::Prop::DoProp(300, st1, TrSys::PropArgs(1, 1));
    //std::cerr << Form("AFTER  MOM %8.4f BETA %8.4f   LOC %8.2f %8.2f %8.2f\n", st1.mom(), st1.bta(), st1.lx(), st1.ly(), st1.lz());
    std::cerr << Form("STATUS %d\n\n", rlt1.status);
    std::cerr << Form("RIG %14.8f\n\n", st1.rig());

    //TrSys::Part st2(org_part);
    //st2.set_mom(0.45);
    //std::cerr << Form("BEFORE MOM %8.4f BETA %8.4f\n", st2.mom(), st2.bta());
    //auto&& rlt2 = TrSys::Prop::DoPropToZ(-60, st2, TrSys::PropArgs(1, 1));
    ////auto&& rlt2 = TrSys::Prop::DoProp(300, st2, TrSys::PropArgs(1, 1));
    //std::cerr << Form("AFTER  MOM %8.4f BETA %8.4f   LOC %8.2f %8.2f %8.2f\n", st2.mom(), st2.bta(), st2.lx(), st2.ly(), st2.lz());
    //std::cerr << Form("STATUS %d\n\n", rlt2.status);
    //
    //TrSys::Part st3(org_part);
    //st3.set_mom(0.425);
    //std::cerr << Form("BEFORE MOM %8.4f BETA %8.4f\n", st3.mom(), st3.bta());
    //auto&& rlt3 = TrSys::Prop::DoPropToZ(-60, st3, TrSys::PropArgs(1, 1));
    ////auto&& rlt3 = TrSys::Prop::DoProp(300, st3, TrSys::PropArgs(1, 1));
    //std::cerr << Form("AFTER  MOM %8.4f BETA %8.4f   LOC %8.2f %8.2f %8.2f\n", st3.mom(), st3.bta(), st3.lx(), st3.ly(), st3.lz());
    //std::cerr << Form("STATUS %d\n\n", rlt3.status);
    //
    //TrSys::Part st4(org_part);
    //st4.set_mom(0.40);
    //std::cerr << Form("BEFORE MOM %8.4f BETA %8.4f\n", st4.mom(), st4.bta());
    //auto&& rlt4 = TrSys::Prop::DoPropToZ(-60, st4, TrSys::PropArgs(1, 1));
    ////auto&& rlt4 = TrSys::Prop::DoProp(300, st4, TrSys::PropArgs(1, 1));
    //std::cerr << Form("AFTER  MOM %8.4f BETA %8.4f   LOC %8.2f %8.2f %8.2f\n", st4.mom(), st4.bta(), st4.lx(), st4.ly(), st4.lz());
    //std::cerr << Form("STATUS %d\n\n", rlt4.status);

    //exit(-1);
    */

    return true;
}

bool Selector::process_prd() {
    if (!process_rich()) return false;
    
    // kalman
    //if (!process_trk_kf()) return false;
    
    // HYC
    //if (!process_hyc()) return false;

    return true;
}

bool Selector::process_sel() {
    return true;

    // Scale events by geomagnetic cutoff (only positive rigidity)
    if (CheckType(Type::ISS)) {
        const float maximum_cutoff = 50.0;
        float max_cf = data_rti.max_IGRF;
        float eft_cf = std::sqrt(maximum_cutoff * max_cf); 
        
        bool signr = false;
        float max_rig = 0;
        for (int it = 0; it < 4; ++it) {
            if (!data_trk.ck_status[it]) continue;
            if (data_trk.ck_rig[it] <= 0.0) { signr = false; break; }
            max_rig = std::max(max_rig, std::abs(data_trk.ck_rig[it]));
            signr = true;
        }

        if (max_cf > 0.0 && signr) {
            float rndm = RndmGenerator.Rndm();
            float wpar[2] = { 0.01, 0.99 };
            if (data_zin > 1) wpar[0] = 0.05;
            if (data_zin > 1) wpar[1] = 0.95;

            float norm  = ((max_rig / eft_cf) - 2.5) / 0.20;
            float thres = wpar[0] + wpar[1] * 0.5 * (1.0 + std::erf(norm));

            if (rndm > thres) return false;
            data_list.weight /= thres;
        }
    }
    
    // Electron study (require EM shower)
    //if (!data_ecal.status) return false;
    
    return true;
}

bool Selector::process_list() {
    data_list.file   = amsch->GetCurrentFile()->GetName();
    data_list.run    = event->Run();
    data_list.event  = event->Event();
    data_list.entry  = amsch->get_tree_entry();
    data_list.utime  = event->UTime();
    data_list.weight = 1.0;

    data_list.header_error = event->fHeader.Error;
    data_list.antimatter_sw_trigger = event->AntiMatteriSWTrigger();

    // testcode
    //if (CheckType(Type::ISS) && data_zin <= 2) {
    //    float rndm = RndmGenerator.Rndm();
    //    float thres = (data_zin == 1) ? 0.01 : 0.1;

    //    if (rndm > thres) return false;
    //    data_list.weight /= thres;
    //}

    return true;
}

bool Selector::process_g4mc() {
    if (!CheckType(Type::MC)) return true;
	
    MCEventgR * prm = event->GetPrimaryMC();
	if (prm == nullptr) return false;
		
    data_g4mc.prm_chrg   = prm->Charge;
	data_g4mc.prm_mass   = prm->Mass;
	data_g4mc.prm_mom    = prm->Momentum;
	data_g4mc.prm_loc[0] = prm->Coo[0];
	data_g4mc.prm_loc[1] = prm->Coo[1];
	data_g4mc.prm_loc[2] = prm->Coo[2];
	data_g4mc.prm_dir[0] = prm->Dir[0];
	data_g4mc.prm_dir[1] = prm->Dir[1];
	data_g4mc.prm_dir[2] = prm->Dir[2];
        
    // Only For Primary Particle
    constexpr Int_t Range[2] = { -1000, -1020 };
    constexpr Int_t TkLay[9] = { -1000, -1007, -1008, -1009, -1010, -1011, -1012, -1013, -1018 };
    constexpr Int_t TfLay[4] = { -1004, -1005, -1015, -1016 };
    constexpr Int_t TdLay[2] = { -1001, -1002 };
    constexpr Int_t EcLay[2] = { -1019, -1020 };
    constexpr Int_t RhLay    =   -1017;

    constexpr Float_t kEngTh = 5.0e-2; 
    std::vector<MCEventgR*> mcev_tk(9, nullptr);
    std::vector<MCEventgR*> mcev_tf(4, nullptr);
    std::vector<MCEventgR*> mcev_td(2, nullptr);
    std::vector<MCEventgR*> mcev_ec(2, nullptr);
    std::vector<MCEventgR*> mcev_rh(1, nullptr);

    for (UInt_t it = 0; it < event->NMCEventg(); it++) {
	    MCEventgR* mcev = event->pMCEventg(it);
        Int_t id = (mcev->trkID == prm->trkID) + (mcev->parentID == prm->trkID) * 2;
        if (mcev == nullptr || id != 1) continue;
	  
        Short_t dec = -1, lay = -1;
        if (mcev->Nskip > Range[0] || mcev->Nskip < Range[1]) continue;
        for (Short_t il = 0; il < 9 && dec < 0; ++il) { if (mcev->Nskip == TkLay[il]) { dec = 0; lay = il; break; } } // Silicon
        for (Short_t il = 0; il < 4 && dec < 0; ++il) { if (mcev->Nskip == TfLay[il]) { dec = 1; lay = il; break; } } // TOF
        for (Short_t il = 0; il < 2 && dec < 0; ++il) { if (mcev->Nskip == TdLay[il]) { dec = 2; lay = il; break; } } // TRD
        for (Short_t il = 0; il < 2 && dec < 0; ++il) { if (mcev->Nskip == EcLay[il]) { dec = 3; lay = il; break; } } // ECAL
        if (mcev->Nskip == RhLay) { dec = 4; lay = 0; } // RICH
        if (dec < 0 || lay < 0) continue;

        if (dec == 0) mcev_tk.at(lay) = mcev;
        if (dec == 1) mcev_tf.at(lay) = mcev;
        if (dec == 2) mcev_td.at(lay) = mcev;
        if (dec == 3) mcev_ec.at(lay) = mcev;
        if (dec == 4) mcev_rh.at(lay) = mcev;
    }

	for (UInt_t icls = 0; icls < event->NTrMCCluster(); icls++) {
		TrMCClusterR* cluster = event->pTrMCCluster(icls);
		if (cluster->GetGtrkID() != prm->trkID) continue;
        
        double keng = std::hypot(cluster->GetMomentum(), prm->Mass) - prm->Mass;
        if (keng < kEngTh) continue;

		int   layJ   = TkDBc::Head->GetJFromLayer(std::fabs(cluster->GetTkId()/100));
        float mes[3] = { cluster->GetXgl()[0], cluster->GetXgl()[1], cluster->GetXgl()[2] };
        float smr[2] = { 0.0, 0.0 };
        float edep   = cluster->GetEdep();
        float mom    = cluster->GetMomentum();

        MCEventgR* mcev = mcev_tk.at(layJ-1);
        if (mcev == nullptr) continue;
        float loc[3] = { mcev->Coo[0], mcev->Coo[1], mcev->Coo[2] };
        float dir[3] = { mcev->Dir[0], mcev->Dir[1], mcev->Dir[2] };
        float tx     = (mcev->Dir[0] / mcev->Dir[2]);
        float ty     = (mcev->Dir[1] / mcev->Dir[2]);
        
        float dz = (mes[2] - loc[2]);
        mes[0] = loc[0] + dz * tx;
        mes[1] = loc[1] + dz * ty;

        int ilay = layJ - 1;
        data_g4mc.tk[ilay]        = true;
        data_g4mc.tk_mom[ilay]    = mom;
        data_g4mc.tk_mes[ilay][0] = mes[0];
        data_g4mc.tk_mes[ilay][1] = mes[1];
        data_g4mc.tk_mes[ilay][2] = mes[2];
        data_g4mc.tk_loc[ilay][0] = loc[0];
        data_g4mc.tk_loc[ilay][1] = loc[1];
        data_g4mc.tk_loc[ilay][2] = loc[2];
        data_g4mc.tk_dir[ilay][0] = dir[0];
        data_g4mc.tk_dir[ilay][1] = dir[1];
        data_g4mc.tk_dir[ilay][2] = dir[2];
        data_g4mc.tk_edep[ilay]   = edep;
	}
    
    for (int it = 0; it < mcev_tk.size(); ++it) {
        MCEventgR* mcev = mcev_tk.at(it);
        if (mcev == nullptr) continue;

        data_g4mc.tkL[it] = true;
        data_g4mc.tkL_mom[it] = mcev->Momentum;
        data_g4mc.tkL_beta[it] = 1.0 / std::sqrt((mcev->Mass/mcev->Momentum) * (mcev->Mass/mcev->Momentum) + 1.0);
    }
		
    if (pBetaH != nullptr) {
        Bool_t  has[4]  = { false, false, false, false };
        Float_t cooz[4] = { 0, 0, 0, 0 };
        for (int il = 0; il < 4; ++il) {
	    	if (!pBetaH->TestExistHL(il)) continue;
	    	TofClusterHR* cls = pBetaH->GetClusterHL(il);
	    	if (!cls->IsGoodTime()) continue;
            has[il]  = true;
            cooz[il] = cls->Coo[2];
        }
        
        std::vector<Short_t> hitsTfI;
        std::vector<Short_t> hitsTfL;
        std::vector<Float_t> hitsTfD;
        for (UInt_t it = 0; it < event->NTofMCCluster(); ++it) {
            TofMCClusterR* cls = event->pTofMCCluster(it);
            if (cls == nullptr) continue;
	    	if (cls->GtrkID != prm->trkID) continue;
            Short_t lay = cls->GetLayer();
            if (!has[lay]) continue;
            hitsTfI.push_back(it);
            hitsTfL.push_back(lay);
            hitsTfD.push_back(std::fabs(cls->Coo[2] - cooz[lay]));
        }

        Short_t hitid[4] = { -1, -1, -1, -1 };
        Float_t distz[4] = { 99., 99., 99., 99. };
        for (UInt_t ih = 0; ih < hitsTfI.size(); ++ih) {
            Short_t lay = hitsTfL.at(ih);
            if (hitsTfD.at(ih) > distz[lay]) continue;
            hitid[lay] = hitsTfI.at(ih);
            distz[lay] = hitsTfD.at(ih);
        }
        
        for (int il = 0; il < 4; ++il) {
            if (!has[il] || hitid[il] < 0) continue;
            TofMCClusterR* cls = event->pTofMCCluster( hitid[il] );
            
            data_g4mc.tf[il]        = true;
            data_g4mc.tf_beta[il]   = cls->Beta;
            data_g4mc.tf_time[il]   = cls->TOF * 1.0e-9;
            data_g4mc.tf_loc[il][0] = cls->Coo[0];
            data_g4mc.tf_loc[il][1] = cls->Coo[1];
            data_g4mc.tf_loc[il][2] = cls->Coo[2];
        }
    }
    
    for (UInt_t icls = 0; icls < event->NTrdMCCluster(); ++icls) {
	    TrdMCClusterR* cluster = event->pTrdMCCluster(icls);
	    if (cluster->GtrkID != prm->trkID) continue;
        int il = cluster->Layer;
        data_g4mc.td[il]        = true;
        data_g4mc.td_mom[il]    = std::sqrt(cluster->Ekin * (cluster->Ekin + 2.0 * prm->Mass));
        //data_g4mc.td_loc[il][0] = cluster->Xgl[0];
        //data_g4mc.td_loc[il][1] = cluster->Xgl[1];
        //data_g4mc.td_loc[il][2] = cluster->Xgl[2];
    }
    
    for (int it = 0; it < mcev_td.size(); ++it) {
        MCEventgR* mcev = mcev_td.at(it);
        if (mcev == nullptr) continue;

        data_g4mc.tdL[it] = true;
        data_g4mc.tdL_mom[it] = mcev->Momentum;
        data_g4mc.tdL_beta[it] = 1.0 / std::sqrt((mcev->Mass/mcev->Momentum) * (mcev->Mass/mcev->Momentum) + 1.0);
        data_g4mc.tdL_loc[it][0] = mcev->Coo[0];
        data_g4mc.tdL_loc[it][1] = mcev->Coo[1];
        data_g4mc.tdL_loc[it][2] = mcev->Coo[2];
        data_g4mc.tdL_dir[it][0] = mcev->Dir[0];
        data_g4mc.tdL_dir[it][1] = mcev->Dir[1];
        data_g4mc.tdL_dir[it][2] = mcev->Dir[2];
    }

    if (mcev_rh.at(0) != nullptr) {
        MCEventgR* mcev = mcev_rh.at(0);

        data_g4mc.rh = true;
        data_g4mc.rh_mom = mcev->Momentum;
        data_g4mc.rh_beta = 1.0 / std::sqrt((mcev->Mass/mcev->Momentum) * (mcev->Mass/mcev->Momentum) + 1.0);
        data_g4mc.rh_loc[0] = mcev->Coo[0];
        data_g4mc.rh_loc[1] = mcev->Coo[1];
        data_g4mc.rh_loc[2] = mcev->Coo[2];
        data_g4mc.rh_dir[0] = mcev->Dir[0];
        data_g4mc.rh_dir[1] = mcev->Dir[1];
        data_g4mc.rh_dir[2] = mcev->Dir[2];
    }
    for (UInt_t icls = 0; icls < event->NRichMCCluster(); ++icls) {
        RichMCClusterR* cls = event->pRichMCCluster(icls);
        if (cls == nullptr) continue;

        RichHitR* hit = (cls->fRichHit >= 0) ? event->pRichHit(cls->fRichHit) : nullptr;
        if (hit == nullptr) continue;
        
        if (cls->Origin[2] < -80.0) continue;
        bool is_prim  = (cls->GtrkID == prm->trkID);
        bool is_noise = (cls->Id == -666);
        bool is_photo = (cls->Id == 50);
        
        int type = -1;
        if (is_prim)  type = 0;
        if (is_photo) type = 1;
        if (is_noise) type = 2;
        if (type < 0) continue;
        data_g4mc.rh_hit_type[type]++;
    }
   
    for (int it = 0; it < mcev_ec.size(); ++it) {
        MCEventgR* mcev = mcev_ec.at(it);
        if (mcev == nullptr) continue;

        data_g4mc.ec[it] = true;
        data_g4mc.ec_mom[it] = mcev->Momentum;
        data_g4mc.ec_loc[it][0] = mcev->Coo[0];
        data_g4mc.ec_loc[it][1] = mcev->Coo[1];
        data_g4mc.ec_loc[it][2] = mcev->Coo[2];
        data_g4mc.ec_dir[it][0] = mcev->Dir[0];
        data_g4mc.ec_dir[it][1] = mcev->Dir[1];
        data_g4mc.ec_dir[it][2] = mcev->Dir[2];
    }

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

	// ISS solar array && backtracing (based on particle)
	bool is_in_shadow = false;
	if (event->NParticle() > 0) {
		AMSPoint ic;
		int idx = event->isInShadow(ic, 0);
		if      (idx == 0) is_in_shadow = false;
		else if (idx == 1) is_in_shadow = true;
	}
	data_rti.is_in_shadow = is_in_shadow;
    if (data_rti.is_in_shadow) return false;

    return true;
}

bool Selector::process_trg() {
	Level1R * lvl1 = event->pLevel1(0);
	if (lvl1 == nullptr) return false;

    int logic = 0;
    int physics = 0;
	if (CheckType(Type::MC)) {
		// Rebuild according to Flight tr.setup
		lvl1->RebuildTrigPatt(logic, physics);
	}
	else {
		lvl1->RestorePhysBPat();
		physics = lvl1->PhysBPatt;
		logic   = lvl1->JMembPatt;
	}

	// trigger info
	bool extTrg        = ((physics&0x80) > 0);
	bool unBiasTrgTOF  = ((physics&0x01) > 0);
	bool unBiasTrgECAL = ((physics&0x40) > 0);
	bool physTrg       = ((physics&0x3e) > 0);

	data_trg.bit = extTrg * 1 + unBiasTrgTOF * 2 + unBiasTrgECAL * 4 + physTrg * 8;
    data_trg.logic = logic;
    data_trg.physics = physics;

    bool ubsbit = ((data_trg.bit&2) == 2);
    bool phybit = ((data_trg.bit&8) == 8);
    //if (!ubsbit && !phybit) return false;

    return true;
}

bool Selector::process_acc() {
	const double TimeOfOneM = 3.335640e+00; // (speed of light)
	if (pBetaH != nullptr) {
		std::vector<float> time;
		for (int il = 0; il < 4; il++) {
			if (!pBetaH->TestExistHL(il)) continue;
			TofClusterHR * cls = pBetaH->GetClusterHL(il);
			if (cls == nullptr) continue;
			time.push_back(cls->Time);
		}
		if (time.size() != 0) {
		    std::sort(time.begin(), time.end());
		    double timeRange[2] = { (time.front() - 5. * TimeOfOneM), (time.back()  + 10. * TimeOfOneM) };
		    double minTimeOfTOF = time.front();

            int count = 0;
		    for (UInt_t icls = 0; icls < event->NAntiCluster(); ++icls) {
		    	AntiClusterR * cls = event->pAntiCluster(icls);
		    	if (cls == nullptr) continue;
		    	if (cls->time < timeRange[0] || cls->time > timeRange[1]) continue;
                count++;
            }

            data_acc.num_cls = count;
        }
	}
    //if (data_acc.num_cls != 0) return false;

    return true;
}

bool Selector::process_tof() {
	data_tof.num_cls  = event->NTofClusterH();
	data_tof.num_beta = event->NBetaH();

    if (pBetaH == nullptr) return true;
    
    const int qopt = TofRecH::kThetaCor|TofRecH::kBirkCor|TofRecH::kReAttCor|TofRecH::kDAWeight|TofRecH::kQ2Q;
	const short pattIdx[4] = { 1, 2, 4, 8 };
	
	data_tof.status = true;
	data_tof.bit = ((pBetaH->pTrTrack   ()) ? 1 : 0) +
                   ((pBetaH->pTrdTrack  ()) ? 2 : 0) +
                   ((pBetaH->pEcalShower()) ? 4 : 0);
	data_tof.patt = 0;
	data_tof.beta = pBetaH->GetBeta();
    data_tof.mass = pBetaH->GetMass();

    if (CheckType(Type::MC)) data_tof.mc_beta = pBetaH->GetMCBeta();
    data_tof.is_tk_match = pBetaH->IsTkTofMatch();

    data_tof.nchi_t = pBetaH->GetNormChi2T();
	data_tof.nchi_c = pBetaH->GetNormChi2C();
	
    int Qall_nlay = 0;
	float Qall_RMS = 0;
	data_tof.Qall = pBetaH->GetQ(Qall_nlay, Qall_RMS);
	data_tof.Qall_nlay = Qall_nlay;
   
    int Zall_nlay = 0;
    float Zall_prob = 0;
    data_tof.Zall = pBetaH->GetZ(Zall_nlay, Zall_prob, 0);
    
    int Qup_nlay = 0;
	float Qup_RMS = 0;
	data_tof.Qup = pBetaH->GetQ(Qup_nlay, Qup_RMS, 2, TofClusterHR::DefaultQOptIonW, 1100);
    
    int Qlw_nlay = 0;
	float Qlw_RMS = 0;
	data_tof.Qlw = pBetaH->GetQ(Qlw_nlay, Qlw_RMS, 2, TofClusterHR::DefaultQOptIonW, 11);
    
    Float_t minT = 0.;
    Float_t btaT[4] = { 0. };
	for (int il = 0; il < 4; ++il) {
		if (!pBetaH->TestExistHL(il)) continue;
		TofClusterHR* cls = pBetaH->GetClusterHL(il);
		if (cls == nullptr) continue;
		if (!cls->IsGoodTime()) continue;
        data_tof.lay[il]    = true;
        data_tof.loc[il][0] = cls->Coo[0];
        data_tof.loc[il][1] = cls->Coo[1];
        data_tof.loc[il][2] = cls->Coo[2];
        data_tof.T[il]      = pBetaH->GetTime(il);
		data_tof.Q[il]      = pBetaH->GetQL(il, 2, qopt);
        data_tof.QL[il]     = pBetaH->GetQL(il);
		data_tof.patt += pattIdx[il];
        minT = std::min(minT, data_tof.T[il]);
        btaT[il] = data_tof.T[il];
	}
	for (int il = 0; il < 4; il++) {
        if (data_tof.T[il] >= 0) { data_tof.T[il] = -1; continue; }
        data_tof.T[il] -= minT;
    }

    // Find Hits in the TOF supper layer (Time)
	std::vector<int> betaHClsId(4, -1);
	if (data_tof.status) {
		for (int it = 0; it < pBetaH->NTofClusterH(); ++it)
			betaHClsId.at(pBetaH->pTofClusterH(it)->Layer) = pBetaH->iTofClusterH(it);
	}

    const double ChrgLimit = 0.875;
	short nearHitNm[4] = { 0, 0, 0, 0 };
	const double TimeOfOneM = 3.335640e+00; // (speed of light)
	for (UInt_t it = 0; it < event->NTofClusterH(); ++it) {
		TofClusterHR * cls = event->pTofClusterH(it);
		if (cls == nullptr) continue;
		if (betaHClsId.at(cls->Layer) == Int_t(it)) continue;
		int    lay  = cls->Layer;
		double chrg = cls->GetQSignal();
        double dtme = std::fabs(cls->Time - btaT[lay]);
		if (chrg > ChrgLimit && dtme < TimeOfOneM) nearHitNm[lay]++;
	}
	for (int il = 0; il < 4; ++il) {
		data_tof.num_extcls[il] = nearHitNm[il];
	}

    bool noiseALL  = (data_tof.num_extcls[0] > 0 && data_tof.num_extcls[1] > 0 && data_tof.num_extcls[2] > 0 && data_tof.num_extcls[3] > 0);
    bool noiseUTOF = (data_tof.num_extcls[0] > 0 && data_tof.num_extcls[1] > 0) && (data_tof.num_extcls[0] + data_tof.num_extcls[1] >= 4);
    bool noiseLTOF = (data_tof.num_extcls[2] > 0 && data_tof.num_extcls[3] > 0) && (data_tof.num_extcls[2] + data_tof.num_extcls[3] >= 4);
    data_tof.extcls_noise = noiseALL * 1 + noiseUTOF * 2 + noiseLTOF * 4;
    
    int ncls[4] = {0};
	data_tof.num_in_time_cls = event->GetNTofClustersInTime(pBetaH, ncls);

    if (!data_tof.status) return false;
    //if (data_tof.patt != 15) return false;
    //if (data_tof.extcls_noise > 0) return false;
    //if (data_tof.num_in_time_cls > 4) return false;

    return true;
}


bool Selector::process_trk() {
    data_trk.num_track = event->NTrTrack();
    
    int fitidInn = (pTrTrack != nullptr) ? pTrTrack->iTrTrackPar(1, 3, 21, data_mass, data_zin) : -1;
    if (pTrTrack == nullptr || fitidInn < 0) return true;

	const unsigned short _hasL1  =   1;
	const unsigned short _hasL2  =   2;
	const unsigned short _hasL34 =  12;
	const unsigned short _hasL56 =  48;
	const unsigned short _hasL78 = 192;
	const unsigned short _hasL9  = 256;

	unsigned short bitPattJ   = pTrTrack->GetBitPatternJ();
	unsigned short bitPattXYJ = pTrTrack->GetBitPatternXYJ();
    
    short isInner   = ((bitPattJ&_hasL34) > 0 &&
	                   (bitPattJ&_hasL56) > 0 &&
					   (bitPattJ&_hasL78) > 0) ? 1 : 0;
	short isL2      = ((bitPattJ&_hasL2) > 0) ?  2 : 0;
	short isL1      = ((bitPattJ&_hasL1) > 0) ?  4 : 0;
	short isL9      = ((bitPattJ&_hasL9) > 0) ?  8 : 0;
	short bitPatt   = isInner + isL2 + isL1 + isL9;

    short isInnerXY = ((bitPattXYJ&_hasL34) > 0 &&
	                   (bitPattXYJ&_hasL56) > 0 &&
					   (bitPattXYJ&_hasL78) > 0) ? 1 : 0;
	short isL2XY    = ((bitPattXYJ&_hasL2) > 0) ?  2 : 0;
	short isL1XY    = ((bitPattXYJ&_hasL1) > 0) ?  4 : 0;
	short isL9XY    = ((bitPattXYJ&_hasL9) > 0) ?  8 : 0;
	short bitPattXY = isInnerXY + isL2XY + isL1XY + isL9XY;

	data_trk.patt   = bitPatt; 
	data_trk.pattXY = bitPattXY;

    // Qrecon: Hu Liu
    //data_trk.QIn = pTrTrack->GetInnerQH(2, data_beta, fitidInn);
	//data_trk.QL2 = (isL2>0) ? pTrTrack->GetLayerJQH(2, 2, data_beta, fitidInn) : 0.0;
	//data_trk.QL1 = (isL1>0) ? pTrTrack->GetLayerJQH(1, 2, data_beta, fitidInn) : 0.0;
	//data_trk.QL9 = (isL9>0) ? pTrTrack->GetLayerJQH(9, 2, data_beta, fitidInn) : 0.0;
    
    float xxxrms;
    int xxxpatt;
    data_trk.QIn = pTrTrack->GetInnerQYJ(xxxrms, xxxpatt, 2, data_beta, fitidInn);
	data_trk.QL2 = (isL2>0) ? pTrTrack->GetLayerQYJ(2, 2, data_beta, fitidInn) : 0.0;
	data_trk.QL1 = (isL1>0) ? pTrTrack->GetLayerQYJ(1, 2, data_beta, fitidInn) : 0.0;
	data_trk.QL9 = (isL9>0) ? pTrTrack->GetLayerQYJ(9, 2, data_beta, fitidInn) : 0.0;
    
    // Qrecon: Hu Liu
    data_trk.QIn_yj = data_trk.QIn;
    data_trk.QIn_hl = pTrTrack->GetInnerQH(2, data_beta, fitidInn);

    double qinmin = 0.0;
    for (int il = 2; il <= 7; ++il) {
        if ((bitPattJ&(1<<il)) == 0) continue;
        double ql = pTrTrack->GetLayerJQH(il+1, 2, data_beta, fitidInn);
        if (ql <= 0.) continue;
        if (qinmin <= 0.) qinmin = ql;
        else              qinmin = std::min(qinmin, ql);
    }
    data_trk.QInMin = qinmin;
   
    const float NoiseQIn = 0.70;
    const float NoiseQL1 = 0.65;
    const float NoiseQL9 = 0.70;
    if (data_trk.QInMin > 0.0 && data_trk.QInMin < NoiseQIn) return false;
    if (data_trk.QL2 > 0.0 && data_trk.QL2 < NoiseQIn) return false;
    if (data_trk.QL1 > 0.0 && data_trk.QL1 < NoiseQL1) return false;
    if (data_trk.QL9 > 0.0 && data_trk.QL9 < NoiseQL9) return false;

    float sendx[9] = { 0.0 };
    float sendy[9] = { 0.0 };
    if (event->GetTrSensorDistance(pTrTrack, sendx, sendy) == 0) {
        for (int il = 0; il < 9; ++il) {
            if (sendx[il] >= 0.0) data_trk.sen[il][0] = sendx[il];
            if (sendy[il] >= 0.0) data_trk.sen[il][1] = sendy[il];
        }
    }

    std::array<std::pair<TrClusterR*, TrClusterR*>, 9> trhits;
    trhits.fill({nullptr, nullptr});
	for (int ilay = 0; ilay < 9; ++ilay) {
		if (!pTrTrack->TestHitLayerJ(ilay+1)) continue;
		TrRecHitR* recHit = pTrTrack->GetHitLJ(ilay+1);
		if (recHit == nullptr) continue;

		int tkid = recHit->GetTkId();
		int mult = recHit->GetResolvedMultiplicity(); //  -1 resolved multiplicty coordinates
		                                              // > 0 requested multiplicty coordinates
        //AMSPoint coo = (ilay==0 || ilay==8) ? (pTrTrack->GetHitCooLJ(ilay+1, 0) + pTrTrack->GetHitCooLJ(ilay+1, 1))*0.5 : pTrTrack->GetHitCooLJ(ilay+1); // (CIEMAT+PG)/2
        AMSPoint coo = TrTrackR::FitCoo[ilay]; // (CIEMAT+PG)/2 after TrTrackR maxspan refit 23
	
        TrClusterR* xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : nullptr;
		TrClusterR* ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : nullptr;
        bool xside = (xcls != nullptr);
        bool yside = (ycls != nullptr);
	
        trhits[ilay].first  = xcls;
        trhits[ilay].second = ycls;

		TkSens tksens(coo, CheckType(Type::MC));
		int sens = (tksens.LadFound()) ? tksens.GetSensor() : -1;

        int xstrip = (xcls == nullptr || !tksens.LadFound() || tksens.GetStripX() < 0) ? -1.0 : tksens.GetStripX();
        int ystrip = (ycls == nullptr || !tksens.LadFound() || tksens.GetStripY() < 0) ? -1.0 : tksens.GetStripY();
        double xeta = (xcls == nullptr || !tksens.LadFound() || std::fabs(tksens.GetImpactPointX()) > 0.5) ? -1.0 : tksens.GetImpactPointX();
        double yeta = (ycls == nullptr || !tksens.LadFound() || std::fabs(tksens.GetImpactPointY()) > 0.5) ? -1.0 : tksens.GetImpactPointY();

        // Qrecon: Hu Liu
        float hl_xchrg  = (!xside) ? -1.0 : pTrTrack->GetLayerJQH(ilay+1, 0, 1, fitidInn);
		float hl_ychrg  = (!yside) ? -1.0 : pTrTrack->GetLayerJQH(ilay+1, 1, 1, fitidInn);
		float hl_xychrg = (!xside || !yside) ? -1.0 : pTrTrack->GetLayerJQH(ilay+1, 2, 1, fitidInn);
        if (hl_xchrg  <= 0.0) hl_xchrg  = 0.0;
        if (hl_ychrg  <= 0.0) hl_ychrg  = 0.0;
        if (hl_xychrg <= 0.0) hl_xychrg = 0.0;

        float yj_xchrg  = (!xside) ? -1.0 : pTrTrack->GetLayerQYJ(ilay+1, 0, 1, fitidInn);
		float yj_ychrg  = (!yside) ? -1.0 : pTrTrack->GetLayerQYJ(ilay+1, 1, 1, fitidInn);
		float yj_xychrg = (!xside || !yside) ? -1.0 : pTrTrack->GetLayerQYJ(ilay+1, 2, 1, fitidInn);
        if (yj_xchrg  <= 0.0) yj_xchrg  = 0.0;
        if (yj_ychrg  <= 0.0) yj_ychrg  = 0.0;
        if (yj_xychrg <= 0.0) yj_xychrg = 0.0;
            
        float sn10 = event->GetTrackerRawSignalRatio(ilay+1, 10, pTrTrack);
        float feet = event->GetTkFeetDist(ilay+1, pTrTrack);

        data_trk.lay[ilay] = (xcls != nullptr) + (ycls != nullptr) * 2;

        data_trk.strip[ilay][0] = xstrip;
        data_trk.strip[ilay][1] = ystrip;
        data_trk.eta[ilay][0] = xeta;
        data_trk.eta[ilay][1] = yeta;

        data_trk.loc[ilay][0] = coo[0];
        data_trk.loc[ilay][1] = coo[1];
        data_trk.loc[ilay][2] = coo[2];

        data_trk.chrg_hl[ilay][0] = hl_xchrg;
        data_trk.chrg_hl[ilay][1] = hl_ychrg;
        data_trk.chrg_hl[ilay][2] = hl_xychrg;
        
        data_trk.chrg_yj[ilay][0] = yj_xchrg;
        data_trk.chrg_yj[ilay][1] = yj_ychrg;
        data_trk.chrg_yj[ilay][2] = yj_xychrg;
        
        data_trk.sn10[ilay] = sn10;
        data_trk.feet[ilay] = feet;
	} // for loop - layer

    int num_inn_x = 0;
    int num_inn_y = 0;
    for (int il = 2; il <= 7; ++il) {
        num_inn_x += ((data_trk.lay[il] % 2) == 1);
        num_inn_y += ((data_trk.lay[il] / 2) == 1);
    }
    data_trk.num_inn_x = num_inn_x;
    data_trk.num_inn_y = num_inn_y;

    // External hits
    std::array<int, 9> ext_min_nhit; ext_min_nhit.fill(0);
    std::array<int, 9> ext_max_nhit; ext_max_nhit.fill(0);
    std::array<double, 9> ext_maxq_hl; ext_maxq_hl.fill(0);
    std::array<double, 9> ext_maxq_yj; ext_maxq_yj.fill(0);
    std::map<int, std::vector<std::array<int, 3>>> ext_hit_idx;
    for (int ih = 0; ih < event->NTrRecHit(); ++ih) {
        TrRecHitR* hit = event->pTrRecHit(ih);
        if (hit == nullptr) continue;
        if (hit->GetXCluster() == nullptr) continue;
        if (hit->GetYCluster() == nullptr) continue;
        
        int xidx = hit->GetXClusterIndex();
        int yidx = hit->GetYClusterIndex();
        
        AMSPoint coo = hit->GetCoord(-1, 2);
        
        double chrgJ = hit->GetQYJ(2);
        double chrgH = hit->GetQH(2);
        if (chrgJ < NoiseQIn) continue;
        if (chrgH < NoiseQIn) continue;
        
        int ilay = hit->GetLayerJ() - 1;
		int tkid = hit->GetTkId();
        if (data_trk.lay[ilay]%2==1 && hit->GetXCluster() == trhits[ilay].first ) continue;
        if (data_trk.lay[ilay]/2==1 && hit->GetYCluster() == trhits[ilay].second) continue;

        ext_maxq_yj[ilay] = std::max(ext_maxq_yj[ilay], chrgJ);
        ext_maxq_hl[ilay] = std::max(ext_maxq_hl[ilay], chrgH);

        ext_hit_idx[tkid].push_back({ilay, xidx, yidx});
    }

    for (auto&& hits : ext_hit_idx) {
        int ilay = -1;
        std::set<int> xcls;
        std::set<int> ycls;
        for (auto&& hit : hits.second) {
            ilay = hit[0];
            xcls.insert(hit[1]);
            ycls.insert(hit[2]);
        }
        int min_nhit = std::min(xcls.size(), ycls.size());
        int max_nhit = std::max(xcls.size(), ycls.size());
        ext_min_nhit[ilay] += min_nhit;
        ext_max_nhit[ilay] += max_nhit;
    }

    for (int il = 0; il < 9; ++il) {
        data_trk.ext_num_hit[il] = ext_min_nhit[il];
        data_trk.ext_chrg_hl[il] = ext_maxq_hl[il];
        data_trk.ext_chrg_yj[il] = ext_maxq_yj[il];
    }

    // Track Pattern
	const short _npatt = 4;
	const short _patt[_npatt] = { 3, 5, 6, 7 };

    // Choutko
    Bool_t ckSwOpt = true;
    Int_t ckRefit = 22; // check fit in recEv
	for (int patt = 0; patt < _npatt && ckSwOpt; ++patt) {
        Stopwatch sw; sw.start();

		int fitid = pTrTrack->iTrTrackPar(1, _patt[patt], ckRefit, data_mass, data_zin);
		if (fitid < 0) continue;

        data_trk.ck_status[patt] = true;
        data_trk.ck_ndof[patt][0] = pTrTrack->GetNdofX(fitid);
        data_trk.ck_ndof[patt][1] = pTrTrack->GetNdofY(fitid);
		data_trk.ck_nchi[patt][0] = pTrTrack->GetNormChisqX(fitid);
		data_trk.ck_nchi[patt][1] = pTrTrack->GetNormChisqY(fitid);
		data_trk.ck_rig[patt]     = pTrTrack->GetRigidity(fitid, 1); // z = 0
		data_trk.ck_crr_rig[patt] = pTrTrack->GetCorrectedRigidity(fitid, 3, 1); // z = 0, 7years new rigidity-scale
   
        if (_patt[patt] == 3) {
            AMSPoint pnt; AMSDir dir;
            pTrTrack->Interpolate(0.0, pnt, dir, fitid);
            data_trk.ck_cen_state[0] =  pnt[0];
		    data_trk.ck_cen_state[1] =  pnt[1];
		    data_trk.ck_cen_state[2] =  pnt[2];
		    data_trk.ck_cen_state[3] = -dir[0];
		    data_trk.ck_cen_state[4] = -dir[1];
		    data_trk.ck_cen_state[5] = -dir[2];
        }
        if (_patt[patt] == 3) {
            AMSPoint pnt; AMSDir dir;
            pTrTrack->Interpolate(195.0, pnt, dir, fitid);
            data_trk.ck_top_state[0] =  pnt[0];
		    data_trk.ck_top_state[1] =  pnt[1];
		    data_trk.ck_top_state[2] =  pnt[2];
		    data_trk.ck_top_state[3] = -dir[0];
		    data_trk.ck_top_state[4] = -dir[1];
		    data_trk.ck_top_state[5] = -dir[2];
        }

        sw.stop();
        data_trk.ck_cpu_time[patt] = sw.time();
    }
    
    return true;
}

bool Selector::process_trk_kf() {
    int fitidInn = (pTrTrack != nullptr) ? pTrTrack->iTrTrackPar(1, 3, 21, data_mass, data_zin) : -1;
    if (pTrTrack == nullptr || fitidInn < 0) return true;
	
    const short _npatt = 4;
	const short _patt[_npatt] = { 3, 5, 6, 7 };
    
    Bool_t kfSwOpt = true;
    Int_t kfRefit = 22; // check fit in recEv
	for (int patt = 0; patt < _npatt && kfSwOpt; ++patt) {
        Stopwatch sw; sw.start();
        
        TrFit trFit;
		int fitid = pTrTrack->iTrTrackPar(trFit, 6, _patt[patt], kfRefit, data_mass, data_zin);
		if (fitid < 0) continue;

        data_trk.kf_status[patt] = true;
        data_trk.kf_ndof[patt][0] = pTrTrack->GetNdofX(fitid);
        data_trk.kf_ndof[patt][1] = pTrTrack->GetNdofY(fitid);
		data_trk.kf_nchi[patt][0] = pTrTrack->GetNormChisqX(fitid);
		data_trk.kf_nchi[patt][1] = pTrTrack->GetNormChisqY(fitid);
		data_trk.kf_cen_rig[patt] = pTrTrack->GetRigidity(fitid, 1); // z = 0
		data_trk.kf_top_rig[patt] = pTrTrack->GetRigidity(fitid, 0); // z = 195
		data_trk.kf_cen_crr_rig[patt] = pTrTrack->GetCorrectedRigidity(fitid, 3, 1); // z = 0, 7years new rigidity-scale
		data_trk.kf_top_crr_rig[patt] = pTrTrack->GetCorrectedRigidity(fitid, 3, 0); // z = 0, 7years new rigidity-scale

        const int ustate = 0; // KALMAN
        if (_patt[patt] == 3) {
            AMSPoint pnt; AMSDir dir; double rig = 0;
            trFit.InterpolateKalman(0.0, pnt, dir, rig, ustate);
            if (rig != 0.0) {
                data_trk.kf_cen_state[0] =  pnt[0];
		        data_trk.kf_cen_state[1] =  pnt[1];
		        data_trk.kf_cen_state[2] =  pnt[2];
		        data_trk.kf_cen_state[3] = -dir[0];
		        data_trk.kf_cen_state[4] = -dir[1];
		        data_trk.kf_cen_state[5] = -dir[2];
            }
        }
        if (_patt[patt] == 3) {
            AMSPoint pnt; AMSDir dir; double rig = 0;
            trFit.InterpolateKalman(195.0, pnt, dir, rig, ustate);
            if (rig != 0.0) {
                data_trk.kf_top_state[0] =  pnt[0];
		        data_trk.kf_top_state[1] =  pnt[1];
		        data_trk.kf_top_state[2] =  pnt[2];
		        data_trk.kf_top_state[3] = -dir[0];
		        data_trk.kf_top_state[4] = -dir[1];
		        data_trk.kf_top_state[5] = -dir[2];
            }
        }
        
        sw.stop();
        data_trk.kf_cpu_time[patt] = sw.time();
    }

    return true;
}

bool Selector::process_trd() {
    data_trd.num_cls = event->NTrdCluster();

    data_trd.num_track = event->NTrdTrack();
    for (int is = 0; is < event->NTrdSegment(); ++is) {
        TrdSegmentR* seg = event->pTrdSegment(is);
        if (seg == nullptr) continue;
        data_trd.num_hit += seg->NTrdCluster();
    }
    
    data_trd.num_Htrack = event->NTrdHTrack();
    for (int is = 0; is < event->NTrdHSegment(); ++is) {
        TrdHSegmentR* seg = event->pTrdHSegment(is);
        if (seg == nullptr) continue;
        data_trd.num_Hhit += seg->NTrdRawHit();
    }

    // Track
    if (pTrdTrack != nullptr) {
		AMSPoint trd_coo;
		AMSDir trd_dir;
		
		// Base on TrdTrack
		trd_coo.setp(pTrdTrack->Coo[0], pTrdTrack->Coo[1], pTrdTrack->Coo[2]);
		trd_dir.SetTheta(pTrdTrack->Theta);
		trd_dir.SetPhi(pTrdTrack->Phi);
		
        data_trd.status   = true;
		data_trd.state[0] = trd_coo[0];
		data_trd.state[1] = trd_coo[1];
		data_trd.state[2] = trd_coo[2];
		data_trd.state[3] = trd_dir[0];
		data_trd.state[4] = trd_dir[1];
		data_trd.state[5] = trd_dir[2];
        data_trd.Qall = pTrdTrack->Q;
        data_trd.Qall_crr = pTrdTrack->Q / 1.2;
    }

    int num_seg_vtx[4][2] = { {0,0}, {0,0}, {0,0}, {0,0} };
    for (int is = 0; is < event->NTrdSegment(); ++is) {
    for (int js = is+1; js < event->NTrdSegment(); ++js) {
        TrdSegmentR* iseg = event->pTrdSegment(is);
        TrdSegmentR* jseg = event->pTrdSegment(js);
        
        if (iseg->NTrdCluster() < 4) continue;
        if (jseg->NTrdCluster() < 4) continue;
        if (iseg->Orientation != jseg->Orientation) continue;
            
        double iseg_tr = (iseg->Orientation == 0) ? -iseg->FitPar[0] : iseg->FitPar[0];
        double iseg_r0 = (iseg->Orientation == 0) ? -iseg->FitPar[1] : iseg->FitPar[1];
        double jseg_tr = (jseg->Orientation == 0) ? -jseg->FitPar[0] : jseg->FitPar[0];
        double jseg_r0 = (jseg->Orientation == 0) ? -jseg->FitPar[1] : jseg->FitPar[1];
        
        double iseg_agl = std::atan(iseg_tr);
        double jseg_agl = std::atan(jseg_tr);

        double vtx_r = (iseg_tr * jseg_r0 - jseg_tr * iseg_r0) / (iseg_tr - jseg_tr);
        double vtx_z = (jseg_r0 - iseg_r0) / (iseg_tr - jseg_tr);
        if (!std::isfinite(vtx_r)) continue;
        if (!std::isfinite(vtx_z)) continue;

        double dr = 0.0;
        if      (vtx_z > 200.0) dr = std::abs((iseg_r0 + iseg_tr * 200.0) - (jseg_r0 + jseg_tr * 200.0));
        else if (vtx_z <  40.0) dr = std::abs((iseg_r0 + iseg_tr *  40.0) - (jseg_r0 + jseg_tr *  40.0));
        if (!std::isfinite(dr)) continue;
        if (dr > 3.0) continue;

        bool within_ortt = std::abs(iseg_tr - jseg_tr) < 0.07;
        if (within_ortt) continue;

        float imax = 0, imin = 200;
        for (int ic = 0; ic < iseg->NTrdCluster(); ++ic) {
            imax = std::max(imax, iseg->pTrdCluster(ic)->Coo[2]);
            imin = std::min(imin, iseg->pTrdCluster(ic)->Coo[2]);
        }
        double iseg_rr = iseg_r0 + 0.5 * iseg_tr * (imin + imax);

        float jmax = 0, jmin = 200;
        for (int ic = 0; ic < jseg->NTrdCluster(); ++ic) {
            jmax = std::max(jmax, jseg->pTrdCluster(ic)->Coo[2]);
            jmin = std::min(jmin, jseg->pTrdCluster(ic)->Coo[2]);
        }
        double jseg_rr = jseg_r0 + 0.5 * jseg_tr * (jmin + jmax);

        bool env = (vtx_z > imax || vtx_z < imin) && (vtx_z > jmax || vtx_z < jmin);

        int zseg = static_cast<int>((vtx_z - 40.0) / 40.0);
        if (zseg < 0) zseg = 0;
        if (zseg > 3) zseg = 3;
        num_seg_vtx[zseg][jseg->Orientation]++;
    }}
    for (int iz = 0; iz < 4; ++iz) {
        data_trd.num_vtx[iz][0] = num_seg_vtx[iz][0];
        data_trd.num_vtx[iz][1] = num_seg_vtx[iz][1];
    }

	const float threshold = 15; //ADC above which will be taken into account in Likelihood Calculation,  15 ADC is the recommended value for the moment.

	if (pTrdTrack != nullptr) {
        TrdKCluster* trdkcls = TrdKCluster::gethead();
        trdkcls->Build(pTrdTrack);
	
        // parameters from TrdKCluster
        float ECAL_Energy_Hypothesis = 0;
        int fitmethod = 1;
        int particle_hypothesis = 1;

        int nhit = 0; //To be filled with number of hits taken into account in Likelihood Calculation
		double llr[3] = {-1, -1, -1}; //To be filled with 3 LikelihoodRatio :  e/P, e/H, P/H
        double llv[3] = {-1, -1, -1}; // To be filled with 3 Likelihood : e, p, H
		trdkcls->GetLikelihoodRatio_TRDRefit(threshold, llr, nhit, ECAL_Energy_Hypothesis, llv, fitmethod, particle_hypothesis);
		if (llr[0] > 0 && llr[1] > 0 && llr[2] > 0) {
            data_trd.tdLLR_status = true;
            data_trd.tdLLR_num_hit = nhit;
            data_trd.tdLLR_ep = llr[0];
            data_trd.tdLLR_eh = llr[1];
            data_trd.tdLLR_ph = llr[2];
            data_trd.tdLL_el = llv[0];
            data_trd.tdLL_pr = llv[1];
            data_trd.tdLL_he = llv[2];
        }
        
        if (data_trd.tdLLR_status) {
            AMSPoint trdKP0 (data_trd.state[0], data_trd.state[1], data_trd.state[2]);
            AMSDir   trdKDir(data_trd.state[3], data_trd.state[4], data_trd.state[5]);

            bool   hlay[20] = { false };
            double hlen[20] = { 0.0 };
            double hamp[20] = { 0.0 };
            double hlz[20]  = { 0.0 };
            for (int ih = 0; ih < nhit; ih++) {
			    TrdKHit* hit = trdkcls->GetHit(ih);
                
                short  lay = hit->TRDHit_Layer;
                double amp = 0.01 * hit->TRDHit_Amp;
                double len = hit->Tube_Track_3DLength_New(&trdKP0, &trdKDir);
                double lz  = hit->TRDHit_z;
                if (len <= 0.05) continue;
                if (amp <= 0.01) continue;

                hlay[lay] = true;
                hlen[lay] = len;
                hamp[lay] = amp;
                hlz [lay] = lz;
            }

            for (int il = 0; il < 20; ++il) {
                if (!hlay[il]) continue;
                data_trd.tdHit_lay.push_back(il);
                data_trd.tdHit_len.push_back(hlen[il]);
                data_trd.tdHit_amp.push_back(hamp[il]);
                data_trd.tdHit_lz .push_back(hlz[il]);
                data_trd.num_tdHit++;
            }

            double tdHitQ_sum = 0.0;
            for (int ih = 0; ih < data_trd.tdHit_amp.size(); ++ih) {
                tdHitQ_sum += data_trd.tdHit_amp.at(ih) / data_trd.tdHit_len.at(ih);
            }
            double tdHitQ = (data_trd.tdHit_amp.size() == 0) ? 0.0 : std::sqrt(tdHitQ_sum / static_cast<double>(data_trd.tdHit_amp.size()));
            data_trd.tdHitQ = tdHitQ;

            const double qsgm_std = (data_zin * data_zin) * ((data_beta > 0) ? (0.60 * std::pow(data_beta, -0.75)) : 0.60);
            std::vector<std::tuple<double, double, int>> hzq;
            for (int ih = 0; ih < data_trd.num_tdHit; ++ih) {
                hzq.push_back(std::make_tuple(data_trd.tdHit_lz.at(ih), data_trd.tdHit_amp.at(ih) / data_trd.tdHit_len.at(ih), data_trd.tdHit_lay.at(ih)));
            }

            int hzq_cntl = 0;
            int hzq_cntu = 0;
            bool is_reach_condition_1st = (hzq.size() <= 2);
            while (!is_reach_condition_1st) {
                if (hzq.size() <= 2) { is_reach_condition_1st = true; break; }
                double nseed = static_cast<double>(hzq.size() - 1);
                double width = std::sqrt((2.75 * 2.75) * (1.0 + 1.0 / nseed) + 2.0 * std::log(nseed)); // 3-sigma

                std::vector<std::pair<double, int>> cand_noise;
                for (int it = 0; it < hzq.size(); ++it) {
                    double qavg = 0.0;
                    double qsgm = qsgm_std * width;
                    for (int ih = 0; ih < hzq.size(); ++ih) {
                        if (it == ih) continue;
                        qavg += std::get<1>(hzq.at(ih));
                    }
                    qavg /= static_cast<double>(hzq.size() - 1);
                    double res = (std::get<1>(hzq.at(it)) - qavg) / qsgm;
                    cand_noise.push_back(std::make_pair(res, it));
                }
                std::sort(cand_noise.begin(), cand_noise.end());
                std::pair<double, int>& cand = (std::abs(std::get<0>(cand_noise.front())) > std::abs(std::get<0>(cand_noise.back()))) ? cand_noise.front() : cand_noise.back();
                
                if (std::abs(cand.first) < 1.0) is_reach_condition_1st = true;
                else {
                    if (cand.first < 0) hzq_cntl++;
                    if (cand.first > 0) hzq_cntu++;
                    hzq.erase(hzq.begin() + cand.second);
                }
            }

            if (hzq.size() > 0) std::sort(hzq.begin(), hzq.end());
            if (hzq.size() > 0) std::reverse(hzq.begin(), hzq.end());
           
            double tdQv = 0.0;
            data_trd.num_tdQ = hzq.size();
            data_trd.num_tdQl = hzq_cntl;
            data_trd.num_tdQu = hzq_cntu;
            for (auto&& hit : hzq) {
                data_trd.tdQl.push_back(std::get<2>(hit));
                data_trd.tdQz.push_back(std::get<0>(hit));
                data_trd.tdQq.push_back(std::sqrt(std::get<1>(hit)));
                tdQv += data_trd.tdQq.back() * data_trd.tdQq.back();
                
                double gbta = (CheckType(Type::MC) && data_g4mc.td[data_trd.tdQl.back()]) ? data_g4mc.td_mom[data_trd.tdQl.back()]/data_g4mc.prm_mass : 0.0;
                data_trd.tdQgb.push_back(gbta);
            }
            data_trd.tdQv = (hzq.size() == 0) ? 0.0 : std::sqrt(tdQv / static_cast<double>(hzq.size()));
            data_trd.tdQv_crr = data_trd.tdQv * ((data_beta > 0) ? 0.943 / (9.63454e-01 * std::erfc(3.03807e+00 * data_beta * data_beta * data_beta) + 1.28178e+00) : 1.0);
        }
	}
	
    if (pTrTrack != nullptr) {
        TrdKCluster* trdkcls = TrdKCluster::gethead();
		int fitid_max = pTrTrack->iTrTrackPar(1, 0, 21, data_mass, data_zin);
		if (fitid_max >= 0) trdkcls->SetTrTrack(pTrTrack, fitid_max);
        
        // parameters from TrdKCluster
        float ECAL_Energy_Hypothesis = 0;
		
        int nhit = 0; //To be filled with number of hits taken into account in Likelihood Calculation
		double llr[3] = {-1, -1, -1}; // To be filled with 3 LikelihoodRatio :  e/P, e/H, P/H
        double llv[3] = {-1, -1, -1}; // To be filled with 3 Likelihood : e, p, H
		trdkcls->GetLikelihoodRatio_TrTrack(threshold, llr, nhit, ECAL_Energy_Hypothesis, llv);
		if (llr[0] > 0 && llr[1] > 0 && llr[2] > 0) {
            data_trd.tkLLR_status = true;
            data_trd.tkLLR_num_hit = nhit;
            data_trd.tkLLR_ep = llr[0];
            data_trd.tkLLR_eh = llr[1];
            data_trd.tkLLR_ph = llr[2];
            data_trd.tkLL_el = llv[0];
            data_trd.tkLL_pr = llv[1];
            data_trd.tkLL_he = llv[2];
        }

        if (data_trd.tkLLR_status) {
            AMSPoint trdKP0  = trdkcls->GetPropogated_TrTrack_P0();
            AMSDir   trdKDir = trdkcls->GetPropogated_TrTrack_Dir();

            bool   hlay[20] = { false };
            double hlen[20] = { 0.0 };
            double hamp[20] = { 0.0 };
            double hlz[20]  = { 0.0 };
            for (int ih = 0; ih < nhit; ih++) {
                TrdKHit* hit = trdkcls->GetHit(ih);

                short  lay = hit->TRDHit_Layer;
                double amp = 0.01 * hit->TRDHit_Amp;
                double len = hit->Tube_Track_3DLength_New(&trdKP0, &trdKDir);
                double lz  = hit->TRDHit_z;
                if (len <= 0.05) continue;
                if (amp <= 0.01) continue;

                hlay[lay] = true;
                hlen[lay] = len;
                hamp[lay] = amp;
                hlz [lay] = lz;
            }

            for (int il = 0; il < 20; ++il) {
                if (!hlay[il]) continue;
                data_trd.tkHit_lay.push_back(il);
                data_trd.tkHit_len.push_back(hlen[il]);
                data_trd.tkHit_amp.push_back(hamp[il]);
                data_trd.tkHit_lz .push_back(hlz[il]);
                data_trd.num_tkHit++;
            }
            
            double tkHitQ_sum = 0.0;
            for (int ih = 0; ih < data_trd.tkHit_amp.size(); ++ih) {
                tkHitQ_sum += data_trd.tkHit_amp.at(ih) / data_trd.tkHit_len.at(ih);
            }
            double tkHitQ = (data_trd.tkHit_amp.size() == 0) ? 0.0 : std::sqrt(tkHitQ_sum / static_cast<double>(data_trd.tkHit_amp.size()));
            data_trd.tkHitQ = tkHitQ;

            const double qsgm_std = (data_zin * data_zin) * ((data_beta > 0) ? (0.60 * std::pow(data_beta, -0.75)) : 0.60);
            std::vector<std::tuple<double, double, int>> hzq;
            for (int ih = 0; ih < data_trd.num_tkHit; ++ih) {
                hzq.push_back(std::make_tuple(data_trd.tkHit_lz.at(ih), data_trd.tkHit_amp.at(ih) / data_trd.tkHit_len.at(ih), data_trd.tkHit_lay.at(ih)));
            }

            int hzq_cntl = 0;
            int hzq_cntu = 0;
            bool is_reach_condition_1st = (hzq.size() <= 2);
            while (!is_reach_condition_1st) {
                if (hzq.size() <= 2) { is_reach_condition_1st = true; break; }
                double nseed = static_cast<double>(hzq.size() - 1);
                double width = std::sqrt((2.75 * 2.75) * (1.0 + 1.0 / nseed) + 2.0 * std::log(nseed)); // 3-sigma

                std::vector<std::pair<double, int>> cand_noise;
                for (int it = 0; it < hzq.size(); ++it) {
                    double qavg = 0.0;
                    double qsgm = qsgm_std * width;
                    for (int jh = 0; jh < hzq.size(); ++jh) {
                        if (it == jh) continue;
                        qavg += std::get<1>(hzq.at(jh));
                    }
                    qavg /= static_cast<double>(hzq.size() - 1);
                    double res = std::abs(std::get<1>(hzq.at(it)) - qavg) / qsgm;
                    cand_noise.push_back(std::make_pair(res, it));
                }
                std::sort(cand_noise.begin(), cand_noise.end());
                std::pair<double, int>& cand = (std::abs(std::get<0>(cand_noise.front())) > std::abs(std::get<0>(cand_noise.back()))) ? cand_noise.front() : cand_noise.back();
                
                if (std::abs(cand.first) < 1.0) is_reach_condition_1st = true;
                else {
                    if (cand.first < 0) hzq_cntl++;
                    if (cand.first > 0) hzq_cntu++;
                    hzq.erase(hzq.begin() + cand.second);
                }
            }

            if (hzq.size() > 0) std::sort(hzq.begin(), hzq.end());
            if (hzq.size() > 0) std::reverse(hzq.begin(), hzq.end());
            
            double tkQv = 0.0;
            data_trd.num_tkQ = hzq.size();
            data_trd.num_tkQl = hzq_cntl;
            data_trd.num_tkQu = hzq_cntu;
            for (auto&& hit : hzq) {
                data_trd.tkQl.push_back(std::get<2>(hit));
                data_trd.tkQz.push_back(std::get<0>(hit));
                data_trd.tkQq.push_back(std::sqrt(std::get<1>(hit)));
                tkQv += data_trd.tkQq.back() * data_trd.tkQq.back();

                double gbta = (CheckType(Type::MC) && data_g4mc.td[data_trd.tkQl.back()]) ? data_g4mc.td_mom[data_trd.tkQl.back()]/data_g4mc.prm_mass : 0.0;
                data_trd.tkQgb.push_back(gbta);
            } 
            data_trd.tkQv = (hzq.size() == 0) ? 0.0 : std::sqrt(tkQv / static_cast<double>(hzq.size()));
            data_trd.tkQv_crr = data_trd.tkQv * ((data_beta > 0) ? 0.943 / (9.63454e-01 * std::erfc(3.03807e+00 * data_beta * data_beta * data_beta) + 1.28178e+00) : 1.0);
        }
	}

    //if (!data_trd.tdLLR_status && !data_trd.tkLLR_status) return false;

    return true;
}

bool Selector::process_ecal() {
	data_ecal.num_shower = event->NEcalShower();
   
    const int NCell = 72;
    for (int ih = 0; ih < event->NEcalHit(); ++ih) {
        EcalHitR* hit = event->pEcalHit(ih);
        int idx = hit->Plane * NCell + hit->Cell;
        if (hit->Edep < 0.01) continue;
        data_ecal.hit_idx.push_back(idx);
        data_ecal.hit_edep.push_back(hit->Edep * 0.001);
        data_ecal.num_hit++;
    }

    if (pEcalShower == nullptr) return true;

    data_ecal.status = true;
	data_ecal.engD   = 1.e-3 * pEcalShower->EnergyD;
	data_ecal.engE   = pEcalShower->EnergyE;
	data_ecal.mvaBDT = pEcalShower->GetEcalBDT();
	
	// Charge estimator based on ECAl-only;
	// it is advised to discard low-rigidity events
	// and to use events having only one EcalShower.
	data_ecal.chrg = pEcalShower->EcalChargeEstimator();
	if (data_ecal.chrg < 1e-3) data_ecal.chrg = 0.0;

	// hadron shower
	EcalHadron::Build(pEcalShower, data_zin);
	data_ecal.hadron_apex = EcalHadron::EcalApex;
	data_ecal.hadron_eng  = EcalHadron::EcalRigidity;

    return true;
}


bool Selector::process_rich() {
    data_rich.num_ring = event->NRichRing();
    
	// official RichRingR - start
	if (pRichRing != nullptr) {
		int kindOfRad = pRichRing->IsNaF() ? 2 : 1;
		int tileOfRad = pRichRing->getTileIndex();

		data_rich.status = true;
		data_rich.beta = pRichRing->getBeta();
        data_rich.kind = pRichRing->IsNaF() ? 2 : 1;
        data_rich.tile = pRichRing->getTileIndex();
        data_rich.refz = pRichRing->AMSTrPars[2];
        data_rich.dist = pRichRing->DistanceTileBorder();
		data_rich.chrg = pRichRing->getCharge2Estimate(true);
		data_rich.chrg = (data_rich.chrg > 1.0e-3) ? std::sqrt(data_rich.chrg) : -1;

        data_rich.isGood = pRichRing->IsGood();
        data_rich.nhit = pRichRing->getHits();
        data_rich.npmt = pRichRing->getPMTs();
        data_rich.prob = pRichRing->getProb();
        data_rich.cstcb = pRichRing->getBetaConsistency();
        data_rich.cstcq = pRichRing->getPMTChargeConsistency();
        
        data_rich.num_clsPE = pRichRing->getPhotoElectrons();
        data_rich.num_expPE = pRichRing->getExpectedPhotoelectrons();
        data_rich.num_colPE = RichHitR::getCollectedPhotoElectrons();
        data_rich.eft_colPE = pRichRing->getPhotoElectrons()/RichHitR::getCollectedPhotoElectrons();
            
        data_rich.num_clsZ1 = pRichRing->ClusterizeZ1();

	    // Rich cuts from Javie/Jorge
        //const short cut_pmt = 3;                     // number of PMTs
        //const float cut_prob = 0.01;                 // Kolmogorov test prob
	    //const float cut_betaCstc[2] = {0.01, 0.005}; // beta_lip vs beta_ciemat consistency (NaF, Aerogel)
	    //const float cut_chrgCstc = 5;                // hit/PMT charge consistency test
	    //const float cut_expPhe[2] = {1, 2};          // expected number of phe (NaF, Aerogel)
	    //const float cut_collPhe[2] = {0.6, 0.4};     // ring phe / total phe (NaF, Aerogel)
		
        //data_rich.isGood = (data_rich.npmt >= cut_pmt) && 
        //                   (data_rich.prob > cut_prob) && 
        //                   (data_rich.cstcb < cut_betaCstc[kindOfRad-1]) && 
        //                   (data_rich.cstcq < cut_chrgCstc) && 
        //                   (data_rich.num_expPE > cut_expPhe[kindOfRad-1]) &&
        //                   (data_rich.eft_colPE > cut_collPhe[kindOfRad-1]);
	}
	// official RichRingR - end
   
    AmsRich rich(event);
    if (rich.status()) {
        data_rich.self_status       = rich.status();
        data_rich.self_kind         = rich.kind();
        data_rich.self_tile         = rich.tile();
        data_rich.self_index        = rich.index();
        data_rich.self_dist         = rich.dist();
        data_rich.self_loc_lx       = rich.locx();
        data_rich.self_loc_ly       = rich.locy();
        data_rich.self_loc_tha      = rich.loctha();
        data_rich.self_loc_phi      = rich.locphi();
        data_rich.self_beta_crr     = rich.beta_crr();
        data_rich.self_is_good_geom = rich.is_good_geom();
        data_rich.self_is_bad_tile  = rich.is_bad_tile();
        data_rich.self_rad_loc[0]   = rich.radp()[0];
        data_rich.self_rad_loc[1]   = rich.radp()[1];
        data_rich.self_rad_loc[2]   = rich.radp()[2];
        data_rich.self_rad_dir[0]   = rich.radd()[0];
        data_rich.self_rad_dir[1]   = rich.radd()[1];
        data_rich.self_rad_dir[2]   = rich.radd()[2];
        data_rich.self_pmt_loc[0]   = rich.pmtp()[0];
        data_rich.self_pmt_loc[1]   = rich.pmtp()[1];
        data_rich.self_pmt_loc[2]   = rich.pmtp()[2];
       
        CherenkovRayTrace ray_trace(
            { rich.radp()[0], rich.radp()[1], rich.radp()[2], rich.radd()[0], rich.radd()[1], rich.radd()[2] }, 
            1.0, rich.index(), rich.kind(), rich.tile()
        ); 
        data_rich.self_expnpe = (pTrTrack != nullptr) ? RichRingR::ComputeNpExp(pTrTrack, 1.0, data_zin) : 0.0;
        data_rich.self_border = ray_trace.border();
        data_rich.self_trace  = ray_trace.trace();
    }

    //for (auto&& hit : rich.hits()) {
    //    if (!hit.status()) continue;
    //    data_rich.hit_chann .push_back(hit.chann());
    //    data_rich.hit_type  .push_back(hit.type());
    //    data_rich.hit_dbeta .push_back(hit.dbeta());
    //    data_rich.hit_rbetaA.push_back(hit.rbetaA());
    //    data_rich.hit_rbetaB.push_back(hit.rbetaB());
    //    data_rich.hit_npe   .push_back(hit.npe());
    //    data_rich.hit_lx    .push_back(hit.cx());
    //    data_rich.hit_ly    .push_back(hit.cy());
    //    data_rich.num_hit++;
    //}
 
    if (rich.status() && rich.kind() != 0) {
        std::vector<CherenkovHit> hits;
        for (auto&& hit : rich.hits()) {
            CherenkovHit chhit(hit.chann(), hit.chann()/16, hit.dbeta(), hit.rbetaA(), hit.rbetaB(), hit.npe(), hit.cx(), hit.cy());
            if (!chhit.status()) continue;
            hits.push_back(chhit);
        }
    
        CherenkovFit chfit(hits, { rich.pmtp()[0], rich.pmtp()[1] }, rich.index(), (rich.kind() == 1) ? CherenkovFit::AGL_BETA_WIDTH : CherenkovFit::NAF_BETA_WIDTH, rich.beta_crr());
        if (chfit.status()) {
            data_rich.self_num_stone = chfit.stns().size();
            data_rich.self_num_cloud = chfit.clds().size();
            data_rich.self_num_tumor = chfit.tmrs().size();
            data_rich.self_num_ghost = chfit.gsts().size();

            data_rich.self_nhit_total = chfit.nhit_total();
            data_rich.self_nhit_stone = chfit.nhit_stone();
            data_rich.self_nhit_cloud = chfit.nhit_cloud();
            data_rich.self_nhit_tumor = chfit.nhit_tumor();
            data_rich.self_nhit_ghost = chfit.nhit_ghost();
            data_rich.self_nhit_other = chfit.nhit_other();
            data_rich.self_nhit_other_inn = chfit.nhit_other_inn();
            data_rich.self_nhit_other_out = chfit.nhit_other_out();
            
            data_rich.self_npe_total = chfit.npe_total();
            data_rich.self_npe_stone = chfit.npe_stone();
            data_rich.self_npe_cloud = chfit.npe_cloud();
            data_rich.self_npe_tumor = chfit.npe_tumor();
            data_rich.self_npe_ghost = chfit.npe_ghost();
            data_rich.self_npe_other = chfit.npe_other();
            data_rich.self_npe_other_inn = chfit.npe_other_inn();
            data_rich.self_npe_other_out = chfit.npe_other_out();

            if (chfit.stns().size() != 0) {
                const CherenkovStone& stn = chfit.stns().at(0);
                data_rich.self_stn_status = stn.status();
                data_rich.self_stn_nhit   = stn.nhit();
                data_rich.self_stn_npmt   = stn.npmt();
                data_rich.self_stn_lx     = stn.cx();
                data_rich.self_stn_ly     = stn.cy();
                data_rich.self_stn_npe    = stn.npe();
                data_rich.self_stn_dist   = stn.dist();
                data_rich.self_stn_nchi   = stn.nchi();
                data_rich.self_stn_chic   = stn.chic();
            }

            if (chfit.clds().size() != 0) {
                const CherenkovCloud& cld = chfit.clds().at(0);
                
                data_rich.self_cld_status   = cld.status();
                data_rich.self_cld_nhit     = cld.nhit();
                data_rich.self_cld_npmt     = cld.npmt();
                data_rich.self_cld_beta     = cld.beta();
                data_rich.self_cld_cbta     = cld.cbta();
                data_rich.self_cld_npe      = cld.npe();
                data_rich.self_cld_nchi     = cld.nchi();
                data_rich.self_cld_misjudge = cld.misjudge();
                data_rich.self_cld_expnpe   = RichRingR::ComputeNpExp(pTrTrack, cld.cbta(), data_zin);
           
                CherenkovRayTrace ray_trace(
                    { rich.radp()[0], rich.radp()[1], rich.radp()[2], rich.radd()[0], rich.radd()[1], rich.radd()[2] }, 
                    cld.cbta(), rich.index(), rich.kind(), rich.tile()
                ); 
                ray_trace.cal(&cld);

                data_rich.self_cld_border   = ray_trace.border();
                data_rich.self_cld_trace    = ray_trace.trace();
                data_rich.self_cld_accuracy = ray_trace.accuracy();

                for (auto&& hit : cld.hits()) {
                    data_rich.self_cldhit_chann.push_back(hit.chann());
                    data_rich.self_cldhit_beta .push_back(hit.beta() );
                    data_rich.self_cldhit_npe  .push_back(hit.npe()  );
                    data_rich.self_cldhit_lx   .push_back(hit.cx()   );
                    data_rich.self_cldhit_ly   .push_back(hit.cy()   );
                }

                if (ray_trace.phis().size() == cld.hits().size()) {
                    for (auto&& phi : ray_trace.phis()) {
                        data_rich.self_cldhit_phi.push_back(phi);
                    }
                }
                else {
                    data_rich.self_cldhit_phi = std::vector<Float_t>(cld.hits().size(), -1.0);
                }
            }
        }
    }
   
    /*
    if (rich.status() && rich.kind() != 0) {
        std::vector<ChHit> hits;
        for (auto&& hit : rich.hits()) {
            ChHit chhit(hit.chann(), hit.chann()/16, hit.dbeta(), hit.rbetaA(), hit.rbetaB(), hit.npe(), hit.cx(), hit.cy());
            if (!chhit.status()) continue;
            hits.push_back(chhit);
        }

        Stopwatch sw; sw.start();
        ChFit chfit(hits, { rich.pmtp()[0], rich.pmtp()[1] }, rich.index(), (rich.kind() == 1) ? ChFit::AGL_BETA_WIDTH : ChFit::NAF_BETA_WIDTH, rich.beta_crr());
        sw.stop();
        
        if (chfit.status()) {
            data_rich.new_status = true;
            data_rich.new_num_stone = chfit.stns().size();
            data_rich.new_num_cloud = chfit.clds().size();
            data_rich.new_num_ghost = chfit.ghts().size();
            
            data_rich.new_nhit_total = chfit.nhit_total();
            data_rich.new_nhit_stone = chfit.nhit_stone();
            data_rich.new_nhit_cloud = chfit.nhit_cloud();
            data_rich.new_nhit_ghost = chfit.nhit_ghost();
            data_rich.new_nhit_other_inn = chfit.nhit_other_inn();
            data_rich.new_nhit_other_out = chfit.nhit_other_out();
            
            data_rich.new_npe_total = chfit.npe_total();
            data_rich.new_npe_stone = chfit.npe_stone();
            data_rich.new_npe_cloud = chfit.npe_cloud();
            data_rich.new_npe_ghost = chfit.npe_ghost();
            data_rich.new_npe_other_inn = chfit.npe_other_inn();
            data_rich.new_npe_other_out = chfit.npe_other_out();

            if (chfit.stns().size() > 0) {
                const ChStone& stn = chfit.stns().at(0);
                data_rich.new_stn_status = stn.status();
                data_rich.new_stn_nhit   = stn.nhit();
                data_rich.new_stn_npmt   = stn.npmt();
                data_rich.new_stn_lx     = stn.lx();
                data_rich.new_stn_ly     = stn.ly();
                data_rich.new_stn_npe    = stn.npe();
                data_rich.new_stn_dist   = stn.dist();
            }
            if (chfit.clds().size() > 0) {
                const ChCloud& cld = chfit.clds().at(0);
                data_rich.new_cld_status = cld.status();
                data_rich.new_cld_nhit   = cld.nhit();
                data_rich.new_cld_npmt   = cld.npmt();
                data_rich.new_cld_nhit_dir = cld.nhit_dir();
                data_rich.new_cld_nhit_rfl = cld.nhit_rfl();
                data_rich.new_cld_nhit_ght = cld.nhit_ght();
                data_rich.new_cld_beta   = cld.beta();
                data_rich.new_cld_cbta   = cld.cbta();
                data_rich.new_cld_nchi   = cld.nchi();
                data_rich.new_cld_npe    = cld.npe();

                for (auto&& cldhit : cld.hits()) {
                    data_rich.new_cldhit_chann.push_back(cldhit.chann());
                    data_rich.new_cldhit_beta .push_back(cldhit.beta() );
                    data_rich.new_cldhit_npe  .push_back(cldhit.npe()  );
                    data_rich.new_cldhit_lx   .push_back(cldhit.lx()   );
                    data_rich.new_cldhit_ly   .push_back(cldhit.ly()   );
                }
            }
        }

        //testcode
        //if (chfit.status() && chfit.clds().size() > 1) {
        //    std::cerr << Form("\n============ RICH  ================= N(%d %d) S(%d %d) B(%14.8f %14.8f)   TIME %14.8f\n", 
        //            data_rich.num_ring, 
        //            data_rich.self_num_cloud, 
        //            data_rich.status, 
        //            data_rich.self_cld_status, 
        //            data_rich.beta, 
        //            data_rich.self_cld_cbta,
        //            sw.time() * 1.0e+3);
        //    std::cerr << Form("NEW RICH STONE %d CLOUD %d GHOST %d\n", chfit.stns().size(), chfit.clds().size(), chfit.ghts().size());

        //    for (auto&& stn : chfit.stns()) {
        //        std::cerr << Form("STONE LOC %14.8f %14.8f (%14.8f) NCHI %14.8f NPE %14.8f\n", stn.lx(), stn.ly(), stn.dist(), stn.nchi(), stn.npe());
        //    }
        //    for (auto&& cld : chfit.clds()) {
        //        std::cerr << Form("CLOUD BTA %14.8f %14.8f NCHI %14.8f NPE %14.8f NHIT %d %d %d\n", cld.beta(), cld.cbta(), cld.nchi(), cld.npe(), cld.nhit_dir(), cld.nhit_rfl(), cld.nhit_ght());
        //        std::cerr << Form("HIT ");
        //        for (auto&& hit : cld.hits()) {
        //            std::cerr << Form("%d ", hit.chann());
        //        }
        //        std::cerr << std::endl;
        //    }
        //}
    }
    */

    return true;
}


bool Selector::process_hyc() {
    TrSys::PartType              part_type  = (data_zin == 1) ? TrSys::PartList::kProton    : TrSys::PartList::kHelium4;
    TrSys::PartType              secp_type  = (data_zin == 1) ? TrSys::PartList::kDeuterium : TrSys::PartList::kHelium3;
    std::vector<TrSys::PartType> part_types = (data_zin == 1) ? TrSys::PartList::kQ1        : TrSys::PartList::kQ2;

    data_hyc.chrg = part_type.chrg();
    data_hyc.mass = part_type.mass();

    TrSys::Tracker ams_trk;
    for (int il = 0; il < 9; ++il) {
        if (data_trk.lay[il] == 0) continue;
        TrSys::TrackerHit hit(il+1, data_trk.lay[il]%2==1, data_trk.lay[il]/2==1, 
            { data_trk.loc[il][0], data_trk.loc[il][1], data_trk.loc[il][2] }, 
            { data_trk.chrg_yj[il][0], data_trk.chrg_yj[il][1], data_trk.chrg_yj[il][2] }
        );
        ams_trk.add_hit(hit);
    }
    
    TrSys::Tof ams_tof;
    for (int il = 0; il < 4; ++il) {
        if (!data_tof.lay[il]) continue;
        TrSys::TofHit hit(il+1,
            { data_tof.loc[il][0], data_tof.loc[il][1], data_tof.loc[il][2] },
            data_tof.T[il] * TrSys::TofHit::kTransNsToCm,
            data_tof.Q[il]
        );
        ams_tof.add_hit(hit);
    }
    
    //TrSys::Trd ams_trd;
    //for (int ih = 0; ih < data_trd.num_tdHit; ++ih) {
    //    if (data_trd.tdHit_len.at(ih) < 0.3) continue;
    //    double dEdx = data_trd.tdHit_amp.at(ih) / data_trd.tdHit_len.at(ih);
    //    double chrg = std::sqrt(dEdx);
    //    TrSys::TrdHit hit(data_trd.tdHit_lay.at(ih), { 0.0, 0.0, data_trd.tdHit_lz.at(ih) }, chrg);
    //    ams_trd.add_hit(hit);
    //}
    //if (ams_trd.hits().size() < 3) ams_trd = TrSys::Trd();
    
    TrSys::Rich ams_rich(data_rich.self_kind);
    for (int ih = 0; ih < data_rich.self_cld_nhit; ++ih) {
        TrSys::RichHit hit(data_rich.self_kind, { 0.0, 0.0, data_rich.self_rad_loc[2] }, 1.0 / data_rich.self_cldhit_beta.at(ih) / data_rich.self_beta_crr);
        ams_rich.add_hit(hit);
    }
    
    TrSys::Part cand_vel_part;

    const int geom_npatt = 4;
    TrSys::Tracker::Pattern geom_patt[geom_npatt] = { TrSys::Tracker::Pattern::kInn, TrSys::Tracker::Pattern::kInnL1, TrSys::Tracker::Pattern::kInnL9, TrSys::Tracker::Pattern::kFullSpan };
    for (int ip = 0; ip < geom_npatt; ++ip) {
        if (ams_trk.get_num_hit_with_l(geom_patt[ip]) == 0) continue;
        Stopwatch sw; sw.start();

        TrSys::GeomTrFit geom_fit(ams_trk.get_hits_with_l(geom_patt[ip]), part_type);
        if (!geom_fit.status()) continue;
        if (ip == 0 && geom_fit.part().is_dynamic()) cand_vel_part = geom_fit.part();
        
        TrSys::Part part_inn = geom_fit.interpolate_to_z(0.);
        TrSys::Part part_top = geom_fit.interpolate_to_z(195.);

        if (!part_inn.is_dynamic()) continue;
        if (!part_top.is_dynamic()) continue;
        
        data_hyc.geom_status[ip] = true;
        data_hyc.geom_ndof_x[ip] = geom_fit.ndof_x();
        data_hyc.geom_ndof_y[ip] = geom_fit.ndof_y();
        data_hyc.geom_nchi_x[ip] = geom_fit.nchi_x();
        data_hyc.geom_nchi_y[ip] = geom_fit.nchi_y();
        data_hyc.geom_nchi_lx[ip]  = geom_fit.nchi_lx();
        data_hyc.geom_nchi_ly[ip]  = geom_fit.nchi_ly();
        data_hyc.geom_nchi_tau[ip] = geom_fit.nchi_tau();
        data_hyc.geom_nchi_rho[ip] = geom_fit.nchi_rho();
        
        double max_norm_lx = 0.0;
        double max_norm_ly = 0.0;
        for (auto&& hit : geom_fit.hits()) {
            if (hit.slx()) max_norm_lx = std::max(max_norm_lx, std::abs(hit.nlx()));
            if (hit.sly()) max_norm_ly = std::max(max_norm_ly, std::abs(hit.nly()));
        }
        data_hyc.geom_max_norm_lx[ip] = max_norm_lx;
        data_hyc.geom_max_norm_ly[ip] = max_norm_ly;

        double max_norm_tau = 0;
        double max_norm_rho = 0;
        for (auto&& arg : geom_fit.args()) {
            max_norm_tau = std::max(max_norm_tau, std::sqrt(arg.tauu * arg.tauu + arg.taul * arg.taul));
            max_norm_rho = std::max(max_norm_rho, std::sqrt(arg.rhou * arg.rhou + arg.rhol * arg.rhol));
        }
        data_hyc.geom_max_norm_tau[ip] = max_norm_tau;
        data_hyc.geom_max_norm_rho[ip] = max_norm_rho;
        
        data_hyc.geom_cen_loc[ip][0] = part_inn.lx();
        data_hyc.geom_cen_loc[ip][1] = part_inn.ly();
        data_hyc.geom_cen_loc[ip][2] = part_inn.lz();
        data_hyc.geom_cen_dir[ip][0] = part_inn.ux();
        data_hyc.geom_cen_dir[ip][1] = part_inn.uy();
        data_hyc.geom_cen_dir[ip][2] = part_inn.uz();
        data_hyc.geom_cen_rig[ip]    = part_inn.rig();
        data_hyc.geom_cen_crr_rig[ip] = event->GetCorrectedRigidity(data_hyc.geom_cen_rig[ip], 0, 3, 2); // PG+CIEMAT, 7years new rigidity-scale(QY) 
        
        data_hyc.geom_top_loc[ip][0] = part_top.lx();
        data_hyc.geom_top_loc[ip][1] = part_top.ly();
        data_hyc.geom_top_loc[ip][2] = part_top.lz();
        data_hyc.geom_top_dir[ip][0] = part_top.ux();
        data_hyc.geom_top_dir[ip][1] = part_top.uy();
        data_hyc.geom_top_dir[ip][2] = part_top.uz();
        data_hyc.geom_top_rig[ip]    = part_top.rig();
        data_hyc.geom_top_crr_rig[ip] = event->GetCorrectedRigidity(data_hyc.geom_top_rig[ip], 0, 3, 2); // PG+CIEMAT, 7years new rigidity-scale(QY) 

        sw.stop();
        data_hyc.geom_cpu_time[ip] = sw.time();
    } 
    //if (!data_hyc.geom_status[0]) return false;

    const int vel_npatt = 4;
    for (int ip = 0; ip < vel_npatt; ++ip) {
        if (!cand_vel_part.is_dynamic()) continue;
        if (ip == 2 && ams_rich.get_num_hit_with_b() == 0) continue;
        if (ip == 3 && ams_rich.get_num_hit_with_b() == 0) continue;
        Stopwatch sw; sw.start();

        std::vector<TrSys::TofHit>     tof_hits  = std::vector<TrSys::TofHit>();
        std::vector<TrSys::RichHit>    rich_hits = std::vector<TrSys::RichHit>();
        std::vector<TrSys::TrackerHit> trk_hits  = std::vector<TrSys::TrackerHit>();
        
        switch (ip) {
            case 0 :
                tof_hits = ams_tof.get_hits_with_t();
                break;
            case 1 :
                tof_hits = ams_tof.get_hits_with_tq();
                trk_hits = ams_trk.get_hits_with_q(TrSys::Tracker::Pattern::kInn);
                break;
            case 2 :
                rich_hits = ams_rich.get_hits_with_b();
                break;
            case 3 :
                tof_hits  = ams_tof.get_hits_with_t();
                rich_hits = ams_rich.get_hits_with_b();
                break;
            default :
                break;
        }

        TrSys::VelocityTrFit vel_fit(tof_hits, rich_hits, trk_hits, cand_vel_part);
        if (!vel_fit.status()) continue;
        
        std::tuple<TrSys::Part, double>&& part_inn = vel_fit.interpolate_to_z(0.);
        std::tuple<TrSys::Part, double>&& part_top = vel_fit.interpolate_to_z(195.);
        
        if (!std::get<0>(part_inn).is_dynamic()) continue;
        if (!std::get<0>(part_top).is_dynamic()) continue;

        data_hyc.vel_status[ip] = true;
        data_hyc.vel_ndof[ip] = vel_fit.ndof();
        data_hyc.vel_nchi[ip] = vel_fit.nchi();
            
        data_hyc.vel_cen_loc[ip][0] = std::get<0>(part_inn).lx();
        data_hyc.vel_cen_loc[ip][1] = std::get<0>(part_inn).ly();
        data_hyc.vel_cen_loc[ip][2] = std::get<0>(part_inn).lz();
        data_hyc.vel_cen_dir[ip][0] = std::get<0>(part_inn).ux();
        data_hyc.vel_cen_dir[ip][1] = std::get<0>(part_inn).uy();
        data_hyc.vel_cen_dir[ip][2] = std::get<0>(part_inn).uz();
        data_hyc.vel_cen_bta[ip]    = 1.0 / std::get<1>(part_inn);
            
        data_hyc.vel_top_loc[ip][0] = std::get<0>(part_top).lx();
        data_hyc.vel_top_loc[ip][1] = std::get<0>(part_top).ly();
        data_hyc.vel_top_loc[ip][2] = std::get<0>(part_top).lz();
        data_hyc.vel_top_dir[ip][0] = std::get<0>(part_top).ux();
        data_hyc.vel_top_dir[ip][1] = std::get<0>(part_top).uy();
        data_hyc.vel_top_dir[ip][2] = std::get<0>(part_top).uz();
        data_hyc.vel_top_bta[ip]    = 1.0 / std::get<1>(part_top);
      
        sw.stop();
        data_hyc.vel_cpu_time[ip] = sw.time();
    }
    
    const int phys_npatt = 4;
    for (int ipart = 0; ipart < 2; ++ipart) {
    for (int ip = 0; ip < phys_npatt; ++ip) {
        if (ip == 2 && ams_rich.get_num_hit_with_b() == 0) continue;
        if (ip == 3 && ams_rich.get_num_hit_with_b() == 0) continue;
        Stopwatch sw; sw.start();
        
        std::vector<TrSys::TrackerHit> trk_hits  = std::vector<TrSys::TrackerHit>();
        std::vector<TrSys::TofHit>     tof_hits  = std::vector<TrSys::TofHit>();
        std::vector<TrSys::RichHit>    rich_hits = std::vector<TrSys::RichHit>();
        
        switch (ip) {
            case 0 :
                trk_hits = ams_trk.get_hits_with_l(TrSys::Tracker::Pattern::kInn);
                tof_hits = ams_tof.get_hits_with_t();
                break;
            case 1 :
                trk_hits = ams_trk.get_hits_with_lq(TrSys::Tracker::Pattern::kInn);
                tof_hits = ams_tof.get_hits_with_tq();
                break;
            case 2 :
                trk_hits  = ams_trk.get_hits_with_l(TrSys::Tracker::Pattern::kInn);
                rich_hits = ams_rich.get_hits_with_b();
                break;
            case 3 :
                trk_hits  = ams_trk.get_hits_with_l(TrSys::Tracker::Pattern::kInn);
                tof_hits  = ams_tof.get_hits_with_t();
                rich_hits = ams_rich.get_hits_with_b();
                break;
            default :
                break;
        }

        
        TrSys::PhysTrFit phys_fit(trk_hits, tof_hits, rich_hits, ((ipart == 0) ? part_type : secp_type));
        if (!phys_fit.status()) continue;
        
        TrSys::Part&& part_inn = phys_fit.interpolate_to_z(0.);
        TrSys::Part&& part_top = phys_fit.interpolate_to_z(195.);

        if (!part_inn.is_dynamic()) continue;
        if (!part_top.is_dynamic()) continue;

        data_hyc.phys_status[ipart][ip] = true;
        data_hyc.phys_ndof_x[ipart][ip] = phys_fit.ndof_x();
        data_hyc.phys_ndof_y[ipart][ip] = phys_fit.ndof_y();
        data_hyc.phys_ndof_b[ipart][ip] = phys_fit.ndof_b();
        data_hyc.phys_nchi_x[ipart][ip] = phys_fit.nchi_x();
        data_hyc.phys_nchi_y[ipart][ip] = phys_fit.nchi_y();
        data_hyc.phys_nchi_b[ipart][ip] = phys_fit.nchi_b();
            
        data_hyc.phys_cen_loc[ipart][ip][0] = part_inn.lx();
        data_hyc.phys_cen_loc[ipart][ip][1] = part_inn.ly();
        data_hyc.phys_cen_loc[ipart][ip][2] = part_inn.lz();
        data_hyc.phys_cen_dir[ipart][ip][0] = part_inn.ux();
        data_hyc.phys_cen_dir[ipart][ip][1] = part_inn.uy();
        data_hyc.phys_cen_dir[ipart][ip][2] = part_inn.uz();
        data_hyc.phys_cen_rig[ipart][ip]    = part_inn.rig();
        data_hyc.phys_cen_bta[ipart][ip]    = part_inn.bta();
        
        data_hyc.phys_top_loc[ipart][ip][0] = part_top.lx();
        data_hyc.phys_top_loc[ipart][ip][1] = part_top.ly();
        data_hyc.phys_top_loc[ipart][ip][2] = part_top.lz();
        data_hyc.phys_top_dir[ipart][ip][0] = part_top.ux();
        data_hyc.phys_top_dir[ipart][ip][1] = part_top.uy();
        data_hyc.phys_top_dir[ipart][ip][2] = part_top.uz();
        data_hyc.phys_top_rig[ipart][ip]    = part_top.rig();
        data_hyc.phys_top_bta[ipart][ip]    = part_top.bta();
        
        sw.stop();
        data_hyc.phys_cpu_time[ipart][ip] = sw.time();
    }}
   
    const int mutr_npatt = 4;
    for (int ip = 0; ip < mutr_npatt; ++ip) {
        if (ip == 2 && ams_rich.get_num_hit_with_b() == 0) continue;
        if (ip == 3 && ams_rich.get_num_hit_with_b() == 0) continue;
        Stopwatch sw; sw.start();

        std::vector<TrSys::TrackerHit> trk_hits  = std::vector<TrSys::TrackerHit>();
        std::vector<TrSys::TofHit>     tof_hits  = std::vector<TrSys::TofHit>();
        std::vector<TrSys::RichHit>    rich_hits = std::vector<TrSys::RichHit>();
        
        switch (ip) {
            case 0 :
                trk_hits = ams_trk.get_hits_with_l(TrSys::Tracker::Pattern::kInn);
                tof_hits = ams_tof.get_hits_with_t();
                break;
            case 1 :
                trk_hits = ams_trk.get_hits_with_lq(TrSys::Tracker::Pattern::kInn);
                tof_hits = ams_tof.get_hits_with_tq();
                break;
            case 2 :
                trk_hits  = ams_trk.get_hits_with_l(TrSys::Tracker::Pattern::kInn);
                rich_hits = ams_rich.get_hits_with_b();
                break;
            case 3 :
                trk_hits  = ams_trk.get_hits_with_l(TrSys::Tracker::Pattern::kInn);
                tof_hits  = ams_tof.get_hits_with_t();
                rich_hits = ams_rich.get_hits_with_b();
                break;
            default :
                break;
        }
        
        TrSys::MassTrFit mutr_fit(trk_hits, tof_hits, rich_hits, part_types);
        
        if (!mutr_fit.status()) continue;
        if (!mutr_fit.phys() || !mutr_fit.phys()->status()) continue;
        
        TrSys::Part&& part_inn = mutr_fit.phys()->interpolate_to_z(0.);
        TrSys::Part&& part_top = mutr_fit.phys()->interpolate_to_z(195.);

        if (!part_inn.is_dynamic()) continue;
        if (!part_top.is_dynamic()) continue;

        data_hyc.mutr_status[ip] = true;
        data_hyc.mutr_ndof_x[ip] = mutr_fit.phys()->ndof_x();
        data_hyc.mutr_ndof_y[ip] = mutr_fit.phys()->ndof_y();
        data_hyc.mutr_ndof_b[ip] = mutr_fit.phys()->ndof_b();
        data_hyc.mutr_nchi_x[ip] = mutr_fit.phys()->nchi_x();
        data_hyc.mutr_nchi_y[ip] = mutr_fit.phys()->nchi_y();
        data_hyc.mutr_nchi_b[ip] = mutr_fit.phys()->nchi_b();
        
        data_hyc.mutr_mass[ip] = mutr_fit.part().mass();
        data_hyc.mutr_sqrm[ip] = mutr_fit.square_mass();
        
        data_hyc.mutr_cen_loc[ip][0] = part_inn.lx();
        data_hyc.mutr_cen_loc[ip][1] = part_inn.ly();
        data_hyc.mutr_cen_loc[ip][2] = part_inn.lz();
        data_hyc.mutr_cen_dir[ip][0] = part_inn.ux();
        data_hyc.mutr_cen_dir[ip][1] = part_inn.uy();
        data_hyc.mutr_cen_dir[ip][2] = part_inn.uz();
        data_hyc.mutr_cen_rig[ip]    = part_inn.rig();
        data_hyc.mutr_cen_bta[ip]    = part_inn.bta();
        
        data_hyc.mutr_top_loc[ip][0] = part_top.lx();
        data_hyc.mutr_top_loc[ip][1] = part_top.ly();
        data_hyc.mutr_top_loc[ip][2] = part_top.lz();
        data_hyc.mutr_top_dir[ip][0] = part_top.ux();
        data_hyc.mutr_top_dir[ip][1] = part_top.uy();
        data_hyc.mutr_top_dir[ip][2] = part_top.uz();
        data_hyc.mutr_top_rig[ip]    = part_top.rig();
        data_hyc.mutr_top_bta[ip]    = part_top.bta();
        
        sw.stop();
        data_hyc.mutr_cpu_time[ip] = sw.time();
    }


    // Cutoff
    if (CheckType(Type::ISS) && data_hyc.geom_status[0]) {
        AMSDir dir(double(data_hyc.geom_top_dir[0][0]), double(data_hyc.geom_top_dir[0][1]), double(data_hyc.geom_top_dir[0][2]));
        data_hyc.TOI_theta = dir.gettheta();
        data_hyc.TOI_phi   = dir.getphi();

        bool opt_Stoermer = true;
        if (opt_Stoermer) {
            double rcutp = 0;
            double rcutn = 0;
            int rltp = event->GetStoermerCutoff(rcutp,  1, dir);
            int rltn = event->GetStoermerCutoff(rcutn, -1, dir);
            if (rltp >= 0 && rltn >= 0) {
                double min_rcut_ev = std::min(std::abs(rcutp), std::abs(rcutn));
                double max_rcut_ev = std::max(std::abs(rcutp), std::abs(rcutn));
                data_hyc.min_Stoermer = min_rcut_ev;
                data_hyc.max_Stoermer = max_rcut_ev;
            }
        }
       
        bool opt_IGRF = true;
        if (opt_IGRF) {
            double rcutp = 0;
            double rcutn = 0;
            int rltp = event->GetIGRFCutoff(rcutp,  1, dir);
            int rltn = event->GetIGRFCutoff(rcutn, -1, dir);
            if (rltp >= 0 && rltn >= 0) {
                double min_rcut_ev = std::min(std::abs(rcutp), std::abs(rcutn));
                double max_rcut_ev = std::max(std::abs(rcutp), std::abs(rcutn));
                data_hyc.min_IGRF = min_rcut_ev;
                data_hyc.max_IGRF = max_rcut_ev;
            }
        }
    }
    
    if (data_hyc.geom_status[0]) {
        int Qup_nlay = 0;
	    float Qup_RMS = 0;
	    float Qup = pBetaH->GetQ(Qup_nlay, Qup_RMS, 2, TofClusterHR::DefaultQOptIonW, 1100, 0, data_hyc.geom_cen_rig[0]);
        
        int Qlw_nlay = 0;
	    float Qlw_RMS = 0;
	    float Qlw = pBetaH->GetQ(Qlw_nlay, Qlw_RMS, 2, TofClusterHR::DefaultQOptIonW, 11, 0, data_hyc.geom_cen_rig[0]);

        data_hyc.tfQup = Qup;
        data_hyc.tfQlw = Qlw;
    }

    // TrMass
    if (true) {
        int TrM_ntktd = TrMass::GetNtrdSegTrk(event);
        
        float TrM_min_prob[9] = {0};
        float TrM_max_prob[9] = {0};
        float TrM_hit_prob[9] = {0};
        float TrM_npick = TrMass::GetNpick(event, pTrTrack, 1, TrM_min_prob, TrM_max_prob, TrM_hit_prob);

        double TrM_beta_tf     = TrMass::GetBeta(event,   1, pTrTrack);
        double TrM_beta_tftk   = TrMass::GetBeta(event,  11, pTrTrack);
        double TrM_beta_tftktd = TrMass::GetBeta(event, 111, pTrTrack);

        float TrM_beta = 0;
        float TrM_mqlv[TrMass::MQL_Nv] = {0};
        double TrM_mql  = TrMass::GetMQL(event, pTrTrack);
        int    TrM_nvar = TrMass::GetMQLv(event, TrM_mqlv, TrM_beta, pTrTrack, 0);
        double TrM_prob = TrMass::GetProb(TrM_mql, TrM_beta);
        
        data_hyc.trM_ntdseg = TrM_ntktd;
        data_hyc.trM_npick  = TrM_npick;
        data_hyc.trM_bta[0] = TrM_beta_tf;
        data_hyc.trM_bta[1] = TrM_beta_tftk;
        data_hyc.trM_bta[2] = TrM_beta_tftktd;
        data_hyc.trM_mql    = TrM_mql;
        data_hyc.trM_prb    = TrM_prob;
    }

    return true;
}

#endif // __Selector_C__
