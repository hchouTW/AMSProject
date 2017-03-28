/**********************************************************
 * Author        : Hsin-Yi Chou
 * Email         : hchou@cern.ch
 * Last modified : 2015-07-22 16:29
 * Filename      : YiProdNtuple.tcc
 * Description   : xxx
 * *******************************************************/
#ifndef __YiProdNtuple_TCC__
#define __YiProdNtuple_TCC__
#include "YiProdNtuple.h"


//---- RecEvent ----//
void RecEvent::init() {
	iBeta       = -1;
	iBetaH      = -1;
	iTrTrack    = -1;
	iEcalShower = -1;
	iTrdTrack   = -1;
	iTrdHTrack  = -1;
	iRichRing   = -1;

	std::fill_n(trackerZJ, 9, 0);
}

bool RecEvent::rebuild(AMSEventR * event) {
	if (event == nullptr) return false;
	fStopwatch.start();
	init();

	int npar    = event->NParticle();
	int nbeta   = event->NBeta();
	int nbetaH  = event->NBetaH();
	int ntrtk   = event->NTrTrack();
	int nshower = event->NEcalShower();
	int ntrdtk  = event->NTrdTrack();
	int ntrdhtk = event->NTrdHTrack();
	int nring   = event->NRichRing();


	/** Particle **/
	bool isParticle = false;
	if (event->NParticle() > 0 && event->pParticle(0) != 0) {
		ParticleR * par = event->pParticle(0);

		iBeta = par->iBeta();
		iBetaH = par->iBetaH();
		iTrTrack = par->iTrTrack();
		iEcalShower = par->iEcalShower();
		iTrdTrack = par->iTrdTrack();
		iTrdHTrack = par->iTrdHTrack();
		iRichRing = par->iRichRing();

		if (!EventBase::checkEventMode(EventBase::MC))
			for (int it = 0; it < event->NBetaH() && iTrTrack>=0; ++it)
				if (event->pBetaH(it)->iTrTrack() == iTrTrack) iBetaH = it;

		isParticle = true;
	}
	if (!isParticle) { init(); fStopwatch.stop(); return false; }

	// Beta Information
	float Beta = 0;
	if      (iBetaH >= 0) Beta = event->pBetaH(iBetaH)->GetBeta();
	else if (iBeta  >= 0) Beta = event->pBeta(iBeta)->Beta;
	else                  Beta = 1;

	// Tracker Information
	if (EventBase::checkEventMode(EventBase::ISS))
		for (int layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerAJ(layJ);
	else
		for (int layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerJ(layJ);

	TrTrackR * TkStPar = 0;
	int        TkStID = -1;
	AMSPoint   TkStCoo;
	AMSDir     TkStDir;
	float      TkStRig = 0;
	if (iTrTrack >= 0) {
		TkStPar = event->pTrTrack(iTrTrack);
		TkStID = TkStPar->iTrTrackPar(1, 3, 21);
		if (TkStID >= 0) {
			TkStCoo = TkStPar->GetP0(TkStID);
			TkStDir = TkStPar->GetDir(TkStID);
			TkStRig = TkStPar->GetRigidity(TkStID);
		}
		else { fStopwatch.stop(); return false; }
	}
		
	// ECAL Information
	// pre-selection (ECAL)
	if (TkStID >= 0 && iEcalShower >= 0) {
		EcalShowerR * ecal = event->pEcalShower(iEcalShower);
		AMSPoint EPnt(ecal->CofG);
		AMSPoint EPntP;
		AMSDir   EPntD;
		TkStPar->Interpolate(EPnt.z(), EPntP, EPntD, TkStID);
		
		float dxPnt = std::fabs(EPntP.x() - EPnt.x());
		float dyPnt = std::fabs(EPntP.y() - EPnt.y());

		float lmtx = 5, lmty = 5; // 5 cm, 5 cm
		if (dxPnt > lmtx || dyPnt > lmty) iEcalShower = -1;
	}

	// TRD Information
	// pre-selection (TRD)
	if (TkStID >= 0 && iTrdTrack >= 0) {
		AMSPoint TrdP(event->pTrdTrack(iTrdTrack)->Coo);
		AMSDir   TrdD(event->pTrdTrack(iTrdTrack)->Theta, event->pTrdTrack(iTrdTrack)->Phi);

		AMSPoint TrdEP;
		AMSDir   TrdED;
		TkStPar->Interpolate(TrdP.z(), TrdEP, TrdED, TkStID);
		float dxP = std::fabs(TrdEP.x() - TrdP.x());
		float dyP = std::fabs(TrdEP.y() - TrdP.y());
		
		float lmtx = 5, lmty = 5; // 5 cm, 5 cm
		if (dxP > lmtx || dyP > lmty) iTrdTrack = -1;
	}

	if (TkStID >= 0 && iTrdHTrack >= 0) {
		AMSPoint TrdHP(event->pTrdHTrack(iTrdHTrack)->Coo);
		AMSDir   TrdHD(event->pTrdHTrack(iTrdHTrack)->Dir);

		AMSPoint TrdHEP;
		AMSDir   TrdHED;
		TkStPar->Interpolate(TrdHP.z(), TrdHEP, TrdHED, TkStID);
		float dxHP = std::fabs(TrdHEP.x() - TrdHP.x());
		float dyHP = std::fabs(TrdHEP.y() - TrdHP.y());
		
		float lmtx = 5, lmty = 5; // 5 cm, 5 cm
		if (dxHP > lmtx || dyHP > lmty) iTrdHTrack = -1;
	}

	if (isParticle) {
		fStopwatch.stop();
		return true;
	}
	else init();
	
	fStopwatch.stop();
	return false;
}


//---- EventBase ----//
EventBase::MODE EventBase::eventMode = EventBase::ISS;
EventBase::VERSION EventBase::eventVersion = EventBase::B950;

EventBase::EventBase() {
	isTreeSelf = false;
	tree = 0;
}

EventBase::~EventBase() {
	if (isTreeSelf) deleteTree();
}

inline void EventBase::fill() {
	if (isTreeSelf == false || tree == nullptr) return;
	tree->Fill();
}

inline void EventBase::setTree(const std::string& name, const std::string& title) {
	isTreeSelf = true;
	tree = new TTree(name.c_str(), title.c_str());
}

inline void EventBase::deleteTree() {
	if (tree != 0) delete tree;
	tree = 0;
}


//---- EventList ----//
EventList::EventList() : EventBase() {
}

EventList::~EventList() {
}

void EventList::initEvent() {
	fList.init();
	fG4mc.init();
}

void EventList::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventList::setEventTree()\n";
#endif

	if (evTree == nullptr) setTree("EventList", "AMS-02 Event List");
	else tree = evTree;

	tree->Branch("list", &fList);
	if (checkEventMode(EventBase::MC))
		tree->Branch("g4mc", &fG4mc);
}

void EventList::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventList::setEnvironment()\n";
#endif

	if (!checkEventMode(EventBase::ISS)) {
		AMSSetupR::SlowControlR::ReadFromExternalFile = false;
	}
	if (checkEventMode(EventBase::MC)) {
		AMSSetupR::LoadISSMC = false;
	}
}

bool EventList::processEvent(AMSEventR * event, AMSChain * chain) {
	initEvent();
	if (event == nullptr)	return false;
	fStopwatch.start();

	// LIST
	fList.run    = event->Run();
	fList.event  = event->Event();
	fList.entry  = chain->Entry();
	fList.weight = EventList::Weight;

	// G4MC
	if (checkEventMode(EventBase::MC)) {
		MCEventgR * primary = event->GetPrimaryMC();
		if (primary == nullptr) return false;

		fG4mc.beamID = primary->TBl + 1;
		fG4mc.primPart.partID   = primary->Particle;
		fG4mc.primPart.chrg     = primary->Charge;
		fG4mc.primPart.mass     = primary->Mass;
		fG4mc.primPart.mom      = primary->Momentum;
		fG4mc.primPart.kEng     = std::sqrt(primary->Momentum * primary->Momentum + primary->Mass * primary->Mass) - primary->Mass;
		fG4mc.primPart.coo[0]   = primary->Coo[0];
		fG4mc.primPart.coo[1]   = primary->Coo[1];
		fG4mc.primPart.coo[2]   = primary->Coo[2];
		fG4mc.primPart.dir[0]   = primary->Dir[0];
		fG4mc.primPart.dir[1]   = primary->Dir[1];
		fG4mc.primPart.dir[2]   = primary->Dir[2];

		std::set<int> PrimSET;
		for (int it = 0; it < event->NMCEventg(); it++) {
			MCEventgR * mcev = event->pMCEventg(it);
			if (mcev == nullptr) continue;
			if (mcev->parentID != 1) continue; // ONLY SAVE PRIM PART
			if (mcev->Momentum < 5.0e-2) continue;
			double momRat = (mcev->Momentum / fG4mc.primPart.mom);
			if (mcev->parentID == primary->trkID && momRat > 5.0e-2) { PrimSET.insert(mcev->trkID); continue; }
			if (momRat < 0.368) continue;
			PrimSET.insert(mcev->trkID);
		}

		std::set<int> TrSET;
		TrSET.insert(primary->trkID);
		std::vector<short> SecTrkID;
		for (int iev = 0; iev < event->NMCEventg(); iev++) {
			MCEventgR * mcev = event->pMCEventg(iev);
			if (mcev == nullptr) continue;
			if (mcev->Momentum < 5.0e-2) continue;
			if (PrimSET.find(mcev->parentID) == PrimSET.end()) continue;
			TrSET.insert(mcev->trkID);
			SecTrkID.push_back(mcev->trkID);
			
			PartMCInfo part;
			part.partID   = mcev->Particle;
			part.chrg     = mcev->Charge;
			part.mass     = mcev->Mass;
			part.mom      = mcev->Momentum;
			part.kEng     = std::sqrt(mcev->Momentum * mcev->Momentum + mcev->Mass * mcev->Mass) - mcev->Mass;
			part.coo[0]   = mcev->Coo[0];
			part.coo[1]   = mcev->Coo[1];
			part.coo[2]   = mcev->Coo[2];
			part.dir[0]   = mcev->Dir[0];
			part.dir[1]   = mcev->Dir[1];
			part.dir[2]   = mcev->Dir[2];
			
			fG4mc.secParts.push_back(part);
		}
		
		for (int icls = 0; icls < event->NTrMCCluster(); icls++) {
			TrMCClusterR * cluster = event->pTrMCCluster(icls);
			if (cluster == nullptr) continue;
			if (cluster->GetMomentum() < 5.0e-2) continue;
			if (TrSET.find(cluster->GetGtrkID()) == TrSET.end()) continue;
			
			int tkid = cluster->GetTkId();
			int layJ = TkDBc::Head->GetJFromLayer(std::fabs(cluster->GetTkId()/100));
			int trkid = cluster->GetGtrkID();
			AMSPoint coo = cluster->GetXgl();
			AMSDir dir = cluster->GetDir();

			HitTRKMCInfo hit;
			hit.layJ    = layJ;
			hit.tkid    = tkid;
			hit.edep    = cluster->GetEdep();
			hit.mom     = cluster->GetMomentum();
			hit.coo[0]  = coo[0];
			hit.coo[1]  = coo[1];
			hit.coo[2]  = coo[2];
			hit.dir[0]  = dir[0];
			hit.dir[1]  = dir[1];
			hit.dir[2]  = dir[2];
			
			if (cluster->GetGtrkID() == primary->trkID) fG4mc.primPart.hits.push_back(hit);
			else {
				for (int ipart = 0; ipart < fG4mc.secParts.size(); ++ipart) {
					if (SecTrkID.at(ipart) != cluster->GetGtrkID()) continue;
					fG4mc.secParts.at(ipart).hits.push_back(hit);
					break;
				}
			}
		}

		std::sort(fG4mc.primPart.hits.begin(), fG4mc.primPart.hits.end(), HitTRKMCInfo_sort());
		
		std::sort(fG4mc.secParts.begin(), fG4mc.secParts.end(), PartMCInfo_sort());
		for (auto && secPart : fG4mc.secParts)
			std::sort(secPart.hits.begin(), secPart.hits.end(), HitTRKMCInfo_sort());

		if (fG4mc.secParts.size() >= 2) {
			VertexMCInfo vertex;
			vertex.status = true;
			vertex.coo[0] = fG4mc.secParts.at(0).coo[0];
			vertex.coo[1] = fG4mc.secParts.at(0).coo[1];
			vertex.coo[2] = fG4mc.secParts.at(0).coo[2];
			fG4mc.primVtx = vertex;
		}
	}

	fStopwatch.stop();
	return selectEvent(event);
}

bool EventList::selectEvent(AMSEventR * event) {
	return true;
}


//---- EventRti ----//
EventRti::EventRti() : EventBase() {
}

EventRti::~EventRti() {
}

void EventRti::initEvent() {
	fRti.init();
}

void EventRti::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventRti::setEventTree()\n";
#endif

	if (!checkEventMode(EventBase::ISS)) return;

	if (evTree == nullptr) setTree("EventRti", "Run time information");
	else tree = evTree;

	tree->Branch("rti", &fRti);
}

void EventRti::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventRti::setEnvironment()\n";
#endif
	if (!checkEventMode(EventBase::ISS)) return;

	AMSSetupR::RTI::UseLatest(6); // pass6
}

bool EventRti::processEvent(AMSEventR * event, AMSChain * chain) {
	if (!checkEventMode(EventBase::ISS)) return true;
	if (event == nullptr)	return false;

	AMSSetupR::RTI rti;
	event->GetRTI(rti);

	bool rebuildRTI = true;
	if (rti.utime == CurrUTime) rebuildRTI = false;
	else CurrUTime = rti.utime;
	
	fRti.uTime = event->UTime();
	
	MGClock::TTime * ttime = MGClock::ConvertFromUTimeToTTime(event->UTime(), MGClock::ClockType::UTC);
	fRti.dateUTC = (ttime->tm_year + 1900) * 10000 + (ttime->tm_mon+1) * 100 + (ttime->tm_mday);
	fRti.timeUTC = (ttime->tm_hour) * 10000 + (ttime->tm_min) * 100 + (ttime->tm_sec);

	// ISS information
	AMSSetupR * setup = AMSSetupR::gethead();
	double RPT[3] = {0};   // ISS coordinates (R, Phi, Theta) (GTOD)
	double VelPT[3] = {0}; // ISS velocity (Vel rad/sec, VelPhi rad, VelTheta rad)
	double YPR[3] = {0};   // ISS attitude (Yaw, Pitch, Roll)
	if (setup != 0) {
		float rpt[3], velpt[3], yaw = 0, pitch = 0, roll = 0;
		setup->getISSTLE(rpt, velpt, fRti.uTime);
		setup->getISSAtt(roll, pitch, yaw, fRti.uTime);
		RPT[0] = rpt[0]; RPT[1] = rpt[1]; RPT[2] = rpt[2];
		VelPT[0] = velpt[0]; VelPT[1] = velpt[1]; VelPT[2] = velpt[2];
		YPR[0] = yaw; YPR[1] = pitch; YPR[2] = roll;
	}
	for (int dir = 0; dir < 3; ++dir) {
		fRti.rptISS[dir] = RPT[dir];
		fRti.velISS[dir] = VelPT[dir];
		fRti.yprISS[dir] = YPR[dir];
	}

	// ISS solar array && backtracing (based on particle)
	int isInShadow = -1;
	if (event->NParticle() > 0) {
		AMSPoint ic;
		int idx = event->isInShadow(ic, 0);
		if (idx == 1) isInShadow = 1;
		else          isInShadow = 0;
	}
	fRti.isInShadow = isInShadow;

	/* Backtracing (based on particle)
	int isFromSpace = -1;
	int backtrace[2][3] = { {-1, -1, -1}, {-1, -1, -1} };
	bool isOpenBacktrace = false;
	while (setup != 0) {
		if (event->NParticle() == 0 || event->pParticle(0) == nullptr) break;
		ParticleR * part = event->pParticle(0);
		TrTrackR * trtk = part->pTrTrack();
		BetaHR * betaH = part->pBetaH();
		if (trtk == nullptr || betaH == nullptr) break;
		int fitid = trtk->iTrTrackPar(1, 3, 21);
		if (fitid < 0) break;
	
		double AMSTheta = part->Theta;
		double AMSPhi   = part->Phi;
		int Chrg[2]     = { 1, -1 };
		double Mom      = std::fabs(trtk->GetRigidity(fitid));
		double Beta     = std::fabs(betaH->GetBeta());

		HeaderR header;
		double RPTO[3] = {0}; // particle final GTOD coo (R, Phi, Theta)
		double GalLong = 0, GalLat = 0; // Galctic coo
		double TraceTime = 0; // time of flight
		double GPT[2] = {0}; // particle final GTOD direction (Phi, Theta)
		int result = -1; // 0 unercutoff (i.e. atmospheric origin), 1 over cutoff (i.e. coming from space), 2 trapped, -1 error
		
		const int nStable = 3;
		double stableFT[nStable] = { 1.00, 1.15, 1.30 };
		int    isOver[nStable] = { -1, -1, -1 }; // -1 error, 0 other, 1 weak, 2 strong
		if (isOpenBacktrace) {
			for (int ift = 0; ift < nStable; ++ift) {
				result = header.do_backtracing(GalLong, GalLat, TraceTime, RPTO, GPT, // output
				                               AMSTheta, AMSPhi, Mom/stableFT[ift], Beta, Chrg[0], // input particle info
																			 RPT, VelPT, YPR, fRti.uTime // input ISS info
																			);
				backtrace[0][ift] = result;

				result = header.do_backtracing(GalLong, GalLat, TraceTime, RPTO, GPT, // output
				                               AMSTheta, AMSPhi, Mom/stableFT[ift], Beta, Chrg[1], // input particle info
																			 RPT, VelPT, YPR, fRti.uTime // input ISS info
																			);
				backtrace[1][ift] = result;

				if      (backtrace[0][ift] == -1 || backtrace[1][ift] == -1) isOver[ift] = -1;
				else if (backtrace[0][ift] ==  1 && backtrace[1][ift] ==  1) isOver[ift] =  2;
				else if (backtrace[0][ift] ==  1 || backtrace[1][ift] ==  1) isOver[ift] =  1;
				else                                                         isOver[ift] =  0;
			}
		}
	
		int successTrails = 0;
		int countFromSpace = 0;
		for (int ift = 0; ift < nStable; ++ift) {
			if (isOver[ift] != -1) successTrails++;
			if (isOver[ift] ==  2) countFromSpace++;
		}

		if (successTrails == 3) {
			if      (countFromSpace == 0) isFromSpace = 0;
			else if (countFromSpace == 1) isFromSpace = 1;
			else if (countFromSpace == 2) isFromSpace = 2;
			else if (countFromSpace == 3) isFromSpace = 3;
		}

		break;
	}

	fRti.isFromSpace     = isFromSpace;
	fRti.backtrace[0][0] = backtrace[0][0];
	fRti.backtrace[0][1] = backtrace[0][1];
	fRti.backtrace[0][2] = backtrace[0][2];
	fRti.backtrace[1][0] = backtrace[1][0];
	fRti.backtrace[1][1] = backtrace[1][1];
	fRti.backtrace[1][2] = backtrace[1][2];
	*/

	if (!rebuildRTI) return selectEvent(event);

	initEvent();
	fStopwatch.start();

	fRti.isInShadow      = isInShadow;
	//fRti.isFromSpace     = isFromSpace;
	//fRti.backtrace[0][0] = backtrace[0][0];
	//fRti.backtrace[0][1] = backtrace[0][1];
	//fRti.backtrace[0][2] = backtrace[0][2];
	//fRti.backtrace[1][0] = backtrace[1][0];
	//fRti.backtrace[1][1] = backtrace[1][1];
	//fRti.backtrace[1][2] = backtrace[1][2];

	fRti.flagRun = true;
	if (event->GetRTIStat() != 0) fRti.flagRun = false;
	if (event->IsBadRun("") != 0) fRti.flagRun = false;
	if (rti.good != 0) fRti.flagRun = false;

	// good second
	bool isGoodSecond = true;
	if ((rti.ntrig/rti.nev) < 0.98 ||
			rti.nerr < 0 || (rti.nerr/rti.nev) > 0.1 ||
			(rti.npart/rti.ntrig) < (0.07/1600*rti.ntrig) || (rti.npart/rti.ntrig) > 0.25 ||
			rti.npart <= 0 || rti.nev > 1800)
		isGoodSecond = false;
	else
		isGoodSecond = true;
	fRti.isGoodSecond = isGoodSecond;

	fRti.zenith = rti.zenith;
	for (int i = 0; i < 4; i++) {
		fRti.cutoffStormer[i] = (std::fabs(rti.cf[i][0]) > std::fabs(rti.cf[i][1])) ?
			std::fabs(rti.cf[i][0]) : std::fabs(rti.cf[i][1]);
		fRti.cutoffIGRF[i] = (std::fabs(rti.cfi[i][0]) > std::fabs(rti.cfi[i][1])) ?
			std::fabs(rti.cfi[i][0]) : std::fabs(rti.cfi[i][1]);
	}
	fRti.radiusGTOD = rti.r;
	fRti.thetaGTOD  = rti.theta;
	fRti.phiGTOD    = rti.phi;
	fRti.latGXY     = rti.glat * TMath::DegToRad();
	fRti.longGXY    = rti.glong * TMath::DegToRad();
	fRti.isInSAA    = rti.IsInSAA();
	fRti.uTime      = event->UTime();
	fRti.liveTime   = rti.lf * rti.nev / (rti.nev + rti.nerr);
	fRti.thetaMAG   = rti.getthetam();
	fRti.phiMAG     = rti.getphim();
	
	for (int dir = 0; dir < 3; ++dir) {
		fRti.rptISS[dir] = RPT[dir];
		fRti.velISS[dir] = VelPT[dir];
		fRti.yprISS[dir] = YPR[dir];
	}

	// |PG-MD| < 35e-4 (L1), 45e-4 (L9) [cm]
	AMSPoint pn1, pn9, pd1, pd9;
	event->GetRTIdL1L9(0, pn1, pd1, event->UTime(), 60);
	event->GetRTIdL1L9(1, pn9, pd9, event->UTime(), 60);
	fRti.trackerAlign[0][0] = pd1.x();
	fRti.trackerAlign[0][1] = pd1.y();
	fRti.trackerAlign[1][0] = pd9.x();
	fRti.trackerAlign[1][1] = pd9.y();

	// Inner Tracker Temperature (Sensor A)
	std::vector<float> tempSensA;
	if (!setup->fSlowControl.GetData("Sensor A", event->Run(), 0, tempSensA) && 
	    tempSensA.size() > 0) {
		fRti.trackerTemp = tempSensA.at(0);
	}
	
	fStopwatch.stop();
	return selectEvent(event);
}

bool EventRti::selectEvent(AMSEventR * event) {
	// RTI cut
	if (!EventBase::checkEventMode(EventBase::ISS)) return true;

	if (!fRti.flagRun) return false;
	if (!fRti.isGoodSecond) return false;
	if (fRti.zenith > 40) return false;
	if (fRti.isInSAA) return false;
	if (fRti.liveTime < 0.5) return false;
	if (fRti.trackerAlign[0][1] > 35. || fRti.trackerAlign[1][1] > 45.) return false;
	
	return true;
}


//---- EventTrg ----//
EventTrg::EventTrg() : EventBase() {
}

EventTrg::~EventTrg() {
}

void EventTrg::initEvent() {
	fTrg.init();
}

void EventTrg::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventTrg::setEventTree()\n";
#endif

	if (evTree == nullptr) setTree("EventTrg", "Anti-Coincidence Counter");
	else tree = evTree;

	tree->Branch("trg", &fTrg);
}

void EventTrg::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventTrg::setEnvironment()\n";
#endif
}

bool EventTrg::processEvent(AMSEventR * event, AMSChain * chain) {
	initEvent();
	if (event == nullptr)	return false;
	fStopwatch.start();

	Level1R * lvl1 = event->pLevel1(0);
	if (lvl1 == nullptr) return false;
	if (checkEventMode(EventBase::MC)) {
		// Rebuild according to Flight tr.setup
		lvl1->RebuildTrigPatt(fTrg.logicPatt, fTrg.physicalPatt);
	}
	else {
		lvl1->RestorePhysBPat();
		fTrg.physicalPatt = lvl1->PhysBPatt;
		fTrg.logicPatt = lvl1->JMembPatt;
	}

	// trigger info
	bool extTrg        = ((fTrg.physicalPatt&0x80) > 0);
	bool unBiasTrgTOF  = ((fTrg.physicalPatt&0x01) > 0);
	bool unBiasTrgECAL = ((fTrg.physicalPatt&0x40) > 0);
	bool physTrg       = ((fTrg.physicalPatt&0x3e) > 0);
	fTrg.bit = extTrg * 1 + unBiasTrgTOF * 2 + unBiasTrgECAL * 4 + physTrg * 8;

	fStopwatch.stop();
	return selectEvent(event);
}

bool EventTrg::selectEvent(AMSEventR * event) {
	if ((fTrg.bit&1) > 0) return false; // external

	return true;
}


//---- EventTof ----//
EventTof::EventTof() : EventBase() {
}

EventTof::~EventTof() {
}

void EventTof::initEvent() {
	fTof.init();
}

void EventTof::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventTof::setEventTree()\n";
#endif

	if (evTree == nullptr) setTree("EventTof", "Time of Flight");
	else tree = evTree;

	tree->Branch("tof", &fTof);
}

void EventTof::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventTof::setEnvironment()\n";
#endif
}

bool EventTof::processEvent(AMSEventR * event, AMSChain * chain) {
	initEvent();
	if (event == nullptr)	return false;
	fStopwatch.start();

	fTof.numOfCluster = event->NTofCluster();
	fTof.numOfClusterH = event->NTofClusterH();
	fTof.numOfBeta = event->NBeta();
	fTof.numOfBetaH = event->NBetaH();

	int tofLayerProj[4] = {0, 1, 1, 0}; // 0,1 := x,y
		
	while (recEv.iBeta >= 0) {
		BetaR * beta = event->pBeta(recEv.iBeta);
		if (beta == nullptr) break;
		fTof.statusBeta = true;
		fTof.beta = beta->Beta;
		short betaPatt = beta->Pattern;
		break;
	}


	TofRecH::BuildOpt = 0; // normal
	const short pattIdx[4] = { 1, 2, 4, 8 };
	while (recEv.iBetaH >= 0 && event->pBetaH(recEv.iBetaH) != nullptr) {
		BetaHR * betaH = event->pBetaH(recEv.iBetaH);

		int ncls[4] = {0};
		fTof.numOfInTimeCluster = event->GetNTofClustersInTime(betaH, ncls);

		fTof.statusBetaH = true;
		fTof.betaH = betaH->GetBeta();
		fTof.betaHBit = ((betaH->pTrTrack   ()) ? 1 : 0) +
                    ((betaH->pTrdTrack  ()) ? 2 : 0) +
                    ((betaH->pEcalShower()) ? 4 : 0);
		fTof.normChisqT  = betaH->GetNormChi2T();
		fTof.normChisqC  = betaH->GetNormChi2C();

		fTof.betaHPatt = 0;
		for (int il = 0; il < 4; il++) {
			if (!betaH->TestExistHL(il)) continue;
			TofClusterHR * cls = betaH->GetClusterHL(il);
			if (cls == nullptr) continue;
			fTof.T[il] = cls->Time;
			fTof.Q[il] = betaH->GetQL(il);
			fTof.betaHPatt += pattIdx[il];
			if (cls->IsGoodTime())
				fTof.betaHGoodTime += pattIdx[il];
		}
		float minTime = (*std::min_element(fTof.T, fTof.T+4));
		for (int il = 0; il < 4; il++) {
			if (fTof.Q[il] < 0) continue;
			fTof.T[il] -= minTime; 
		}

		int nlay = 0;
		float Qall_RMS = 0;
		fTof.Qall = betaH->GetQ(nlay, Qall_RMS);
		
		//float Qupper, Qupper_RMS, Qlower, Qlower_RMS;
		//Qupper = betaH->GetQ(nlay, Qupper_RMS, 2, TofClusterHR::DefaultQOpt, 1100);
		//Qlower = betaH->GetQ(nlay, Qlower_RMS, 2, TofClusterHR::DefaultQOpt, 11);

		//TofChargeHR tofQH = betaH->gTofCharge();
		//float Z[2], Zupper, Zlower, Z_prob[2], Zupper_prob, Zlower_prob;
	  //Z[0] = tofQH.GetZ(nlay, Z_prob[0], 0); // max prob charge (int)
		//Z[1] = tofQH.GetZ(nlay, Z_prob[1], 1); // next to max prob charge (int)
		//Zupper = tofQH.GetZ(nlay, Zupper_prob, 0, 1100);
		//Zlower = tofQH.GetZ(nlay, Zlower_prob, 0, 11);

		break;
	} // while loop - ibetaH > 0


	// Find Hits in the TOF supper layer (Time)
	std::vector<int> betaHClsId(4, -1);
	if (fTof.statusBetaH) {
		BetaHR * betaH = event->pBetaH(recEv.iBetaH);
		for (int it = 0; it < betaH->NTofClusterH(); ++it)
			betaHClsId.at(betaH->pTofClusterH(it)->Layer) = betaH->iTofClusterH(it);
	}

	bool   isHasTime[2] = { false, false };
	double avgTime[2] = {0, 0};
	double avgChrg[2] = {0, 0};
	for (int it = 0; it < 2; ++it) {
		isHasTime[it] = (betaHClsId.at(2*it+0) >= 0 || betaHClsId.at(2*it+1) >= 0);
		if (isHasTime[it]) {
			TofClusterHR * ucls = (betaHClsId.at(2*it+0) >= 0) ? event->pTofClusterH(betaHClsId.at(2*it+0)) : nullptr;
			TofClusterHR * lcls = (betaHClsId.at(2*it+1) >= 0) ? event->pTofClusterH(betaHClsId.at(2*it+1)) : nullptr;
			double chrg  = (((ucls) ? ucls->GetQSignal() : 0.) + ((lcls) ? lcls->GetQSignal() : 0.)) / ((ucls!=nullptr) + (lcls!=nullptr));
			double wgval = ((ucls) ? (ucls->Time/ucls->ETime/ucls->ETime) : 0.) + ((lcls) ? (lcls->Time/lcls->ETime/lcls->ETime) : 0.);
			double sumwg = ((ucls) ? (1./ucls->ETime/ucls->ETime) : 0.) + ((lcls) ? (1./lcls->ETime/lcls->ETime) : 0.);
			double value = wgval / sumwg;
			avgChrg[it] = chrg;
			avgTime[it] = value;
		}
	}
	
	int    nearHitId[2] = { -1, -1 };
	double nearHitDt[2] = { 0, 0 };
	for (int it = 0; it < event->NTofClusterH(); ++it) {
		TofClusterHR * cls = event->pTofClusterH(it);
		if (cls == nullptr) continue;
		if (betaHClsId.at(cls->Layer) < 0) continue;
		if (betaHClsId.at(cls->Layer) == it) continue;
		int    slay = (cls->Layer / 2);
		double dltT = std::fabs(cls->Time - avgTime[slay]) / cls->ETime;
		if (nearHitId[slay] < 0 || dltT < nearHitDt[slay]) {
			nearHitId[slay] = it;
			nearHitDt[slay] = dltT;
		} 
	}

	for (int it = 0; it < 2; ++it) {
		if (nearHitId[it] < 0) continue;
		TofClusterHR * cls = event->pTofClusterH(nearHitId[it]);
		int    lay  = cls->Layer;
		double chrg = cls->GetQSignal();
		double time = (cls->Time - avgTime[it]);
		
		fTof.extClsL[it] = lay;
		fTof.extClsQ[it] = chrg;
		fTof.extClsT[it] = time;
	}


	/*
	TofRecH::BuildOpt = 1; // TrTrack independent
	TofRecH::ReBuild(0);
	while (event->NBetaH() > 0) {
		BetaHR * betaHs = event->pBetaH(0);
		if (betaHs == nullptr) break;

		fTof.statusBetaHs = true;
		fTof.betaHs = betaHs->GetBeta();
		fTof.betaHBits = ((betaHs->pTrTrack   ()) ? 1 : 0) +
                     ((betaHs->pTrdTrack  ()) ? 2 : 0) +
                     ((betaHs->pEcalShower()) ? 4 : 0);
		fTof.normChisqTs  = betaHs->GetNormChi2T();
		fTof.normChisqCs  = betaHs->GetNormChi2C();

		fTof.betaHPatts = 0;
		for (int il = 0; il < 4; il++) {
			if (!betaHs->TestExistHL(il)) continue;
			TofClusterHR * cls = betaHs->GetClusterHL(il);
			if (cls == nullptr) continue;
			fTof.Ts[il] = cls->Time;
			fTof.Qs[il] = betaHs->GetQL(il);
			fTof.betaHPatts += pattIdx[il];
			if (cls->IsGoodTime())
				fTof.betaHGoodTimes += pattIdx[il];
		}
		float minTime = (*std::min_element(fTof.Ts, fTof.Ts+4));
		for (int il = 0; il < 4; il++) {
			if (fTof.Qs[il] < 0) continue;
			fTof.Ts[il] -= minTime; 
		}

		int nlays = 0;
		float Qalls_RMS = 0;
		fTof.Qalls = betaHs->GetQ(nlays, Qalls_RMS);

		if ((fTof.betaHBits&2) == 2) {
			TrdTrackR * trd = betaHs->pTrdTrack();
			AMSDir trd_dir(trd->Theta, trd->Phi); trd_dir = trd_dir * -1;
			float dir[3] = { float(trd_dir[0]), float(trd_dir[1]), float(trd_dir[2]) };
			float stt[6] = { trd->Coo[0], trd->Coo[1], trd->Coo[2], -dir[0], -dir[1], -dir[2] };
			std::copy(stt, stt+6, fTof.betaHStates);
		}

		break;
	}
	TofRecH::BuildOpt = 0; // normal
	*/

	fStopwatch.stop();
	return selectEvent(event);
}

bool EventTof::selectEvent(AMSEventR * event) {
	return true;
}


//---- EventAcc ----//
EventAcc::EventAcc() : EventBase() {
}

EventAcc::~EventAcc() {
}

void EventAcc::initEvent() {
	fAcc.init();
}

void EventAcc::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventAcc::setEventTree()\n";
#endif

	if (evTree == nullptr) setTree("EventAcc", "Anti-Coincidence Counter");
	else tree = evTree;

	tree->Branch("acc", &fAcc);
}

void EventAcc::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventAcc::setEnvironment()\n";
#endif
}

bool EventAcc::processEvent(AMSEventR * event, AMSChain * chain) {
	initEvent();
	if (event == nullptr)	return false;
	fStopwatch.start();
	
	event->RebuildAntiClusters();
	fAcc.numOfCluster = event->NAntiCluster();

	const double TimeOfOneM = 3.335640e+00; // (speed of light)
	while (recEv.iBetaH >= 0 && event->pBetaH(recEv.iBetaH) != nullptr) {
		BetaHR * betaH = event->pBetaH(recEv.iBetaH);
		std::vector<float> time;
		for (int il = 0; il < 4; il++) {
			if (!betaH->TestExistHL(il)) continue;
			TofClusterHR * cls = betaH->GetClusterHL(il);
			if (cls == nullptr) continue;
			time.push_back(cls->Time);
		}
		if (time.size() == 0) break;
		std::sort(time.begin(), time.end());
		float timeRange[2] = { (time.front() - 5. * TimeOfOneM), (time.back()  + 10. * TimeOfOneM) };
		float minTimeOfTOF = time.front();

		for (int icls = 0; icls < event->NAntiCluster(); ++icls) {
			AntiClusterR * cls = event->pAntiCluster(icls);
			if (cls == nullptr) continue;
			if (cls->time < timeRange[0] || cls->time > timeRange[1]) continue;
			ClsACCInfo clsACC;
			clsACC.sector = cls->Sector;
			clsACC.time   = cls->time - minTimeOfTOF;
			clsACC.rawQ   = (MGNumc::EqualToZero(cls->rawq) ? -1 : cls->rawq);
			clsACC.coo[0] = cls->AntiCoo[0];
			clsACC.coo[1] = cls->AntiCoo[1];
			clsACC.coo[2] = cls->AntiCoo[2];
			fAcc.clusters.push_back(clsACC);
		}
		if (fAcc.clusters.size() >= 2)
			std::sort(fAcc.clusters.begin(), fAcc.clusters.end(), ClsACCInfo_sort());

		break;
	}


	fStopwatch.stop();
	return selectEvent(event);
}

bool EventAcc::selectEvent(AMSEventR * event) {
	return true;
}


//---- EventTrk ----//
EventTrk::EventTrk() : EventBase() {
}

EventTrk::~EventTrk() {
}

void EventTrk::initEvent() {
	fTrk.init();
}

void EventTrk::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventTrk::setEventTree()\n";
#endif

	if (evTree == nullptr) setTree("EventTrk", "Silicon Tracker");
	else tree = evTree;

	tree->Branch("trk", &fTrk);
}

void EventTrk::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventTrk::setEnvironment()\n";
#endif

	// Disable overwriting of datacards from file
	TRMCFFKEY_DEF::ReadFromFile = 0;
	TRFITFFKEY_DEF::ReadFromFile = 0;
	TRFITFFKEY.magtemp = 0;

	// Enable latest alignment
	if (checkEventMode(EventBase::ISS)) {
		TkDBc::UseFinal();
	}
}

bool EventTrk::processEvent(AMSEventR * event, AMSChain * chain) {
	initEvent();
	if (event == nullptr)	return false;
	fStopwatch.start();

	// Beta Info
	float Beta = 0;
	if      (recEv.iBetaH >= 0) Beta = event->pBetaH(recEv.iBetaH)->GetBeta();
	else if (recEv.iBeta  >= 0) Beta = event->pBeta(recEv.iBeta)->Beta;
	else                        Beta = 1;

	const int qopt = TrClusterR::kAsym | TrClusterR::kGain | TrClusterR::kLoss | TrClusterR::kMIP;
	std::map<int, int> trackIDMap;
	for (int itr = 0; (itr <= event->NTrTrack() && recEv.iTrTrack >= 0); ++itr) {
		int jtr = itr - 1;
		if (jtr == recEv.iTrTrack) continue;
		TrTrackR * trtk = event->pTrTrack( ((jtr==-1)?recEv.iTrTrack:jtr) );
		if (trtk == nullptr) continue;
		
		if (jtr == -1) {
			double BTdist = 0;
			fTrk.beamID = event->GetBeamPos(BTdist, trtk);
			fTrk.beamDist = BTdist;
		}
		
		TrackInfo track;

		const unsigned short _hasL1  =   1;
		const unsigned short _hasL2  =   2;
		const unsigned short _hasL34 =  12;
		const unsigned short _hasL56 =  48;
		const unsigned short _hasL78 = 192;
		const unsigned short _hasL9  = 256;
		track.bitPattJ = trtk->GetBitPatternJ();
		track.bitPattXYJ = trtk->GetBitPatternXYJ();
	
		short isInner   = ((track.bitPattJ&_hasL34) > 0 &&
		                   (track.bitPattJ&_hasL56) > 0 &&
											 (track.bitPattJ&_hasL78) > 0) ? 1 : 0;
		short isInnerXY = ((track.bitPattXYJ&_hasL34) > 0 &&
		                   (track.bitPattXYJ&_hasL56) > 0 &&
											 (track.bitPattXYJ&_hasL78) > 0) ? 2 : 0;
		short isL2    = ((track.bitPattJ  &_hasL2) > 0) ?   4 : 0;
		short isL2XY  = ((track.bitPattXYJ&_hasL2) > 0) ?   8 : 0;
		short isL1    = ((track.bitPattJ  &_hasL1) > 0) ?  16 : 0;
		short isL1XY  = ((track.bitPattXYJ&_hasL1) > 0) ?  32 : 0;
		short isL9    = ((track.bitPattJ  &_hasL9) > 0) ?  64 : 0;
		short isL9XY  = ((track.bitPattXYJ&_hasL9) > 0) ? 128 : 0;
		short bitPatt = isInner + isInnerXY + isL2 + isL2XY + isL1 + isL1XY + isL9 + isL9XY;
	
		short fitidInn = trtk->iTrTrackPar(1, 3, 21);
		if (fitidInn < 0) continue;
		track.bitPatt = bitPatt; 
		track.QIn = trtk->GetInnerQ_all(Beta, fitidInn).Mean; 
		track.QL2 = (isL2>0) ? trtk->GetLayerJQ(2, Beta) : -1;
		track.QL1 = (isL1>0) ? trtk->GetLayerJQ(1, Beta) : -1;
		track.QL9 = (isL9>0) ? trtk->GetLayerJQ(9, Beta) : -1;

		const short _nalgo = 2;
		const short _algo[_nalgo] = { 1, 3 };
		const short _npatt = 4;
		const short _patt[_npatt] = { 3, 5, 6, 7 };
		for (int algo = 0; algo < _nalgo; algo++) {
			for (int patt = 0; patt < _npatt; patt++) {
				int fitid = trtk->iTrTrackPar(_algo[algo], _patt[patt], 21);
				if (fitid < 0) continue;
				AMSPoint coo = trtk->GetP0(fitid);
				AMSDir   dir = trtk->GetDir(fitid);
	
				track.status[algo][patt] = true;
				track.rigidity[algo][patt] = trtk->GetRigidity(fitid);
				track.chisq[algo][patt][0] = trtk->GetNormChisqX(fitid);
				track.chisq[algo][patt][1] = trtk->GetNormChisqY(fitid);
				track.state[algo][patt][0] = coo[0];
				track.state[algo][patt][1] = coo[1];
				track.state[algo][patt][2] = coo[2];
				track.state[algo][patt][3] = -dir[0];
				track.state[algo][patt][4] = -dir[1];
				track.state[algo][patt][5] = -dir[2];
			
				for (int il = 0; il < 9; ++il) {
					AMSPoint pntLJ;
					AMSDir   dirLJ;
					trtk->InterpolateLayerJ(il+1, pntLJ, dirLJ, fitid);
					
					track.stateLJ[algo][patt][il][0] = pntLJ[0];
					track.stateLJ[algo][patt][il][1] = pntLJ[1];
					track.stateLJ[algo][patt][il][2] = pntLJ[2];
					track.stateLJ[algo][patt][il][3] = -dirLJ[0];
					track.stateLJ[algo][patt][il][4] = -dirLJ[1];
					track.stateLJ[algo][patt][il][5] = -dirLJ[2];

					TkSens tksens(pntLJ, EventBase::checkEventMode(EventBase::MC));
					if (tksens.LadFound()) {
						short tkid = tksens.GetLadTkID();
						short sens = tksens.GetSensor();
						short mult = tksens.GetMultIndex();
						float xloc = (tksens.GetStripX() + tksens.GetImpactPointX() + 640);
						float yloc = (tksens.GetStripY() + tksens.GetImpactPointY());
						
						track.localID[algo][patt][il][0] = tkid;
						track.localID[algo][patt][il][1] = sens;
						track.localID[algo][patt][il][2] = mult;
						track.localLJ[algo][patt][il][0] = xloc;
						track.localLJ[algo][patt][il][1] = yloc;
					}
				}
			} // for loop - pattern
		}

		for (int ilay = 0; ilay < 9; ++ilay) {
			if (!trtk->TestHitLayerJ(ilay+1)) continue;
			TrRecHitR * recHit = trtk->GetHitLJ(ilay+1);
			if (recHit == nullptr) continue;

			int tkid = recHit->GetTkId();
			int mult = recHit->GetResolvedMultiplicity(); //  -1 resolved multiplicty coordinates
			                                             // > 0 requested multiplicty coordinates
			AMSPoint coo = recHit->GetCoord(mult, 3); // (CIEMAT+PG)/2
  		
			TkSens tksens(coo, EventBase::checkEventMode(EventBase::MC));
			int sens = (tksens.LadFound()) ? tksens.GetSensor() : -1;
		
			TrClusterR * xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : nullptr;
			TrClusterR * ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : nullptr;

			short clsIdX = (xcls) ? recHit->GetXClusterIndex() : -1;
			short clsIdY = (ycls) ? recHit->GetYClusterIndex() : -1;

			short side = (xcls ? 1 : 0) + (ycls ? 2 : 0);

			float xchrg = (xcls) ? std::sqrt(recHit->GetSignalCombination(0, qopt, 1)) : -1.;
			float ychrg = std::sqrt(recHit->GetSignalCombination(1, qopt, 1));

			float xloc = (xcls) ? (xcls->GetCofG() + xcls->GetSeedAddress()) : (recHit->GetDummyX() + 640);
			float yloc = (ycls->GetCofG() + ycls->GetSeedAddress());
  		
			// origin info
			short xseedAddr = (xcls) ? xcls->GetSeedAddress() : -1;
			short yseedAddr = (ycls) ? ycls->GetSeedAddress() : -1;

			short xseedIndx = (xcls) ? xcls->GetSeedIndex() : -1;
			short yseedIndx = (ycls) ? ycls->GetSeedIndex() : -1;

			std::vector<float> xstripSig;
			std::vector<float> xstripSgm;
			for (int it = 0; (xcls!=nullptr) && (it < xcls->GetLength()); ++it) {
				xstripSig.push_back(xcls->GetSignal(it));
				xstripSgm.push_back(xcls->GetNoise(it));
			}
			
			std::vector<float> ystripSig;
			std::vector<float> ystripSgm;
			for (int it = 0; (ycls!=nullptr) && (it < ycls->GetLength()); ++it) {
				ystripSig.push_back(ycls->GetSignal(it));
				ystripSgm.push_back(ycls->GetNoise(it));
			}
			
			HitTRKInfo hit;
			hit.clsId[0] = clsIdX;
			hit.clsId[1] = clsIdY;
			hit.layJ     = ilay+1;
			hit.tkid     = tkid;
			hit.sens     = sens;
			hit.mult     = mult; 
			hit.side     = side;
			hit.coo[0]   = coo[0];
			hit.coo[1]   = coo[1];
			hit.coo[2]   = coo[2];
			hit.chrg[0]  = xchrg;
			hit.chrg[1]  = ychrg;
			
			hit.cofgX     = xloc;
			hit.seedAddrX = xseedAddr;
			hit.seedIndxX = xseedIndx;
			hit.stripSigX = xstripSig;
			hit.stripSgmX = xstripSgm;
			
			hit.cofgY     = yloc;
			hit.seedAddrY = yseedAddr;
			hit.seedIndxY = yseedIndx;
			hit.stripSigY = ystripSig;
			hit.stripSgmY = ystripSgm;
	
			track.hits.push_back(hit);
		} // for loop - layer
		if (track.hits.size() > 1) std::sort(track.hits.begin(), track.hits.end(), HitTRKInfo_sort());
	
		int offTrID = ((jtr==-1)?recEv.iTrTrack:jtr);
		int selTrID = fTrk.tracks.size();
		trackIDMap[offTrID] = selTrID;
		
		fTrk.tracks.push_back(track);
	}

	// It isnot match with official one.
	if (fTrk.tracks.size() != event->NTrTrack()) return false;

	// Vertex
	TrReconQ reconQ(event);
	reconQ.BuildVertex();
	AMSPoint recQ_CVtx, recQ_DVtx;
	int recQ_TrID[2] = { -1, -1 };
	if (reconQ.GetVertex(recQ_CVtx, recQ_DVtx, recQ_TrID) != 0) {
		std::map<int, int>::iterator ita = trackIDMap.find(recQ_TrID[0]);
		std::map<int, int>::iterator itb = trackIDMap.find(recQ_TrID[1]);
		if (ita != trackIDMap.end() && itb != trackIDMap.end()) {
			std::pair<int, int> trId = std::minmax(ita->second, itb->second);
			VertexInfo vertex;
			
			vertex.vtxZqu[0] = recQ_CVtx[0];
			vertex.vtxZqu[1] = recQ_CVtx[1];
			vertex.vtxZqu[2] = recQ_CVtx[2];
			vertex.vtxZqu[3] = recQ_DVtx[0];
			vertex.vtxZqu[4] = recQ_DVtx[1];

			VertexR vtxR(event->pTrTrack(recQ_TrID[0]), event->pTrTrack(recQ_TrID[1]));
			if (vtxR.IsFilled()) {
				AMSPoint cooR = vtxR.getvert();
				AMSDir   dirR(vtxR.gettheta(), vtxR.getphi());
				vertex.vtxR[0] = cooR[0];
				vertex.vtxR[1] = cooR[1];
				vertex.vtxR[2] = cooR[2];
				vertex.vtxR[3] = -dirR[0];
				vertex.vtxR[4] = -dirR[1];
				vertex.vtxR[5] = -dirR[2];
				vertex.vtxR[6] = vtxR.getmom();
				vertex.vtxR[7] = (vtxR.getchi2() / vtxR.getndof());
			}

			vertex.trackID.push_back(trId.first);
			vertex.trackID.push_back(trId.second);
			fTrk.vertices.push_back(vertex);
		}
	}

	//---- Find External Hits ----//
	struct CandHit {
		int   hitId;
		short xcls, ycls;
		float xchg, ychg;
		float mtch;
	};
	struct CandHit_sort {
		bool operator() (const CandHit & hit1, const CandHit & hit2) {
			return (hit1.mtch < hit2.mtch);
		}
	};

	std::map<short, std::vector<CandHit> > candHits;

	for (int ih = 0; ih < event->NTrRecHit(); ++ih) {
		TrRecHitR * recHit = event->pTrRecHit(ih);
		if (recHit == nullptr) continue;
		
		TrClusterR * xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : nullptr;
		TrClusterR * ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : nullptr;
		if (xcls == nullptr || ycls == nullptr) continue;
		if (xcls->Used() || ycls->Used()) continue;
		
		bool isNoise = (xcls->GetSeedSN() < 4.5 || ycls->GetSeedSN() < 4.5); // now, 4.5 is best
		if (isNoise) continue;
		
		float xchrg = std::sqrt(recHit->GetSignalCombination(0, qopt, 1));
  	float ychrg = std::sqrt(recHit->GetSignalCombination(1, qopt, 1));
		float dchrg = std::fabs(xchrg - ychrg) / (xchrg + ychrg);
		if (xchrg < 0.7 || ychrg < 0.7) continue;
		
		const float NSfluc = 0.08;
		float chrgNS = std::sqrt(2.0 * (xcls->GetSeedNoise()/xcls->GetSeedSignal()/xcls->GetSeedSignal() + 
		                                  ycls->GetSeedNoise()/ycls->GetSeedSignal()/ycls->GetSeedSignal()) +
															 NSfluc * NSfluc);
		if (dchrg > chrgNS) continue;

		CandHit candh;
		candh.hitId = ih;
		candh.xcls = recHit->GetXClusterIndex();
		candh.ycls = recHit->GetYClusterIndex();
		candh.xchg = xchrg;
		candh.ychg = ychrg;
		candh.mtch = dchrg;

		candHits[recHit->GetTkId()].push_back(candh);
	}

	for(auto &elmap : candHits) {
		std::sort(elmap.second.begin(), elmap.second.end(), CandHit_sort());
		std::vector<CandHit> & hits = elmap.second;
		if (hits.size() < 2) continue;
		float              fineScore = 0;
		std::vector<short> fineCand;
		for (int ih = 0; ih < hits.size(); ++ih) {
			float sumscore = 0;
			std::vector<short> cand;
			std::set<short> xset;
			std::set<short> yset;
			for (int jh = ih; jh < hits.size(); ++jh) {
				CandHit & hit = hits.at(jh);
				if (xset.find(hit.xcls) != xset.end()) continue;
				else xset.insert(hit.xcls);
				if (yset.find(hit.ycls) != yset.end()) continue;
				else yset.insert(hit.ycls);
				sumscore += hit.mtch;
				cand.push_back(jh);
			}
			float score = sumscore / cand.size();
			bool change = false;
			if (cand.size() < fineCand.size()) continue;
			else if (cand.size() > fineCand.size()) change = true;
			else if (score < fineScore) change = true;
			if (!change) continue;
			fineScore = score;
			fineCand = cand;
		}
		std::vector<CandHit> fineHits;
		for (auto && it : fineCand) {
			fineHits.push_back( hits.at(it) );
		}
		elmap.second = fineHits;
	}

	for(auto &elmap : candHits) {
		std::sort(elmap.second.begin(), elmap.second.end(), CandHit_sort());
		for(auto &candHit : elmap.second) {
			TrRecHitR * recHit = event->pTrRecHit(candHit.hitId);
			if (recHit == nullptr) continue;
			TrClusterR * xcls = recHit->GetXCluster();
			TrClusterR * ycls = recHit->GetYCluster();

			short clsIdX = candHit.xcls;
			short clsIdY = candHit.ycls;
			short layJ   = recHit->GetLayerJ();
			short tkid   = recHit->GetTkId();  	
			short side = 3;
		
			float xchrg = candHit.xchg;
			float ychrg = candHit.ychg;

			float xloc = (xcls->GetCofG() + xcls->GetSeedAddress());
			float yloc = (ycls->GetCofG() + ycls->GetSeedAddress());
			
			short xseedAddr = xcls->GetSeedAddress();
			short yseedAddr = ycls->GetSeedAddress();

			short xseedIndx = xcls->GetSeedIndex();
			short yseedIndx = ycls->GetSeedIndex();

			std::vector<float> xstripSig;
			std::vector<float> xstripSgm;
			for (int it = 0; (xcls!=nullptr) && (it < xcls->GetLength()); ++it) {
				xstripSig.push_back(xcls->GetSignal(it));
				xstripSgm.push_back(xcls->GetNoise(it));
			}
			
			std::vector<float> ystripSig;
			std::vector<float> ystripSgm;
			for (int it = 0; (ycls!=nullptr) && (it < ycls->GetLength()); ++it) {
				ystripSig.push_back(ycls->GetSignal(it));
				ystripSgm.push_back(ycls->GetNoise(it));
			}
			
			AMSPoint coo;
			int mult = -1;
			float distLJ = 500;
			if (fTrk.tracks.size() > 0 && fTrk.tracks.at(0).status[0][0]) {
				TrackInfo & track = fTrk.tracks.at(0);
				for (int im = 0; im < recHit->GetMultiplicity(); im++) {
					AMSPoint refcoo = recHit->GetCoord(im, 3);
					float dx = (refcoo[0] - track.stateLJ[0][0][layJ-1][0]);
					float dy = (refcoo[1] - track.stateLJ[0][0][layJ-1][1]);
					float dist = std::sqrt(dx * dx + dy * dy);
					if (dist > distLJ) continue;
					distLJ = dist;
					coo = refcoo;
					mult = im;
				}
			}
			else {
				mult = 0;
				coo = recHit->GetCoord(0, 3);
			}
  		
			TkSens tksens(coo, EventBase::checkEventMode(EventBase::MC));
			int sens = (tksens.LadFound()) ? tksens.GetSensor() : -1;
			
			HitTRKInfo hit;
			hit.clsId[0] = clsIdX;
			hit.clsId[1] = clsIdY;
			hit.layJ     = layJ;
			hit.tkid     = tkid;
			hit.sens     = sens;
			hit.mult     = mult; 
			hit.side     = side;
			hit.coo[0]   = coo[0];
			hit.coo[1]   = coo[1];
			hit.coo[2]   = coo[2];
			hit.chrg[0]  = xchrg;
			hit.chrg[1]  = ychrg;
			
			hit.cofgX     = xloc;
			hit.seedAddrX = xseedAddr;
			hit.seedIndxX = xseedIndx;
			hit.stripSigX = xstripSig;
			hit.stripSgmX = xstripSgm;
			
			hit.cofgY     = yloc;
			hit.seedAddrY = yseedAddr;
			hit.seedIndxY = yseedIndx;
			hit.stripSigY = ystripSig;
			hit.stripSgmY = ystripSgm;
	
			fTrk.otherHits.push_back(hit);
		}
	}
	if (fTrk.otherHits.size() > 1) 
		std::sort(fTrk.otherHits.begin(), fTrk.otherHits.end(), HitTRKInfo_sort());

	fStopwatch.stop();
	return selectEvent(event);
}

bool EventTrk::selectEvent(AMSEventR * event) {
	return true;
}


//---- EventTrd ----//
EventTrd::EventTrd() : EventBase() {
}

EventTrd::~EventTrd() {
}

void EventTrd::initEvent() {
	fTrd.init();
}

void EventTrd::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventTrd::setEventTree()\n";
#endif

	if (evTree == nullptr) setTree("EventTrd", "Transition Radiation");
	else tree = evTree;

	tree->Branch("trd", &fTrd);
}

void EventTrd::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventTrd::setEnvironment()\n";
#endif

	if (checkEventMode(EventBase::BT)) {
		TrdKCluster::IsReadGlobalAlignment = false;
	}
	if (checkEventMode(EventBase::MC)) {
		TrdKCluster::ForceReadAlignment = false;
		TrdKCluster::ForceReadCalibration = false;
	}
}

bool EventTrd::processEvent(AMSEventR * event, AMSChain * chain) {
	initEvent();
	if (event == nullptr)	return false;
	fStopwatch.start();

	fTrd.numOfTrack = event->NTrdTrack();
	fTrd.numOfHTrack = event->NTrdHTrack();
		
	// numOfHSegVtx (by HY.Chou)
	// number of TRDH segments that make a vertex with TrTrack
	// between Z = 80 cm and 200 cm
	int nseg = event->nTrdHSegment();
	TrTrackR * trtk = (recEv.iTrTrack >= 0) ? event->pTrTrack(recEv.iTrTrack) : nullptr;
	int        trId = (trtk) ? trtk->iTrTrackPar(1, 3, 21) : -1;
	if (trtk && trId >= 0) {
		const double trCooZ = 120.;
		AMSPoint trPnt; AMSDir trDir;
		trtk->Interpolate(trCooZ, trPnt, trDir);
		double trX = trPnt[0];
		double trY = trPnt[1];
		double trTanX = trDir[0] / trDir[2];
		double trTanY = trDir[1] / trDir[2];
		double trX0   = trX - trTanX * trCooZ;
		double trY0   = trY - trTanY * trCooZ;
		double trNRig = std::fabs(trtk->GetRigidity(trId));
		if (trNRig < 0.8) trNRig = 0.8;

		const double MsCooSgmFact = 2.843291e-01;
		const double MsDirSgmFact = 7.108227e-03;
		double msSgmC = MsCooSgmFact / trNRig;
		double msSgmM = MsDirSgmFact / trNRig;
		double msSgmR = std::sqrt(msSgmC*msSgmC + msSgmM*msSgmM*trCooZ*trCooZ);

		short numOfVtx[2] = { 0, 0 };
		std::vector<std::pair<double, double> > vtxCooXZ;
		std::vector<std::pair<double, double> > vtxCooYZ;

		const double SgmLimit = 7.;
		const double ZLimit[2] = { 200., 60. };
		for (int iseg = 0; iseg < event->NTrdHSegment(); ++iseg) {
			TrdHSegmentR * seg = event->pTrdHSegment(iseg);
			if (seg == nullptr) continue;
			double segR0 = (seg->r - seg->m * seg->z);
			double trR   = ((seg->d==0) ? (trX + trTanX * (seg->z - trCooZ)) : (trY + trTanY * (seg->z - trCooZ)));
			double trM   = ((seg->d==0) ? trTanX : trTanY);
			double dr    = std::fabs((trR - seg->r) / std::sqrt(seg->er*seg->er + msSgmC*msSgmC));
			double dm    = std::fabs((trM - seg->m) / std::sqrt(seg->em*seg->em + msSgmM*msSgmM));
			bool isTrSeg = (dr < SgmLimit && dm < SgmLimit);
			if (isTrSeg) continue;

			double sgmSegR = std::sqrt(seg->er*seg->er+seg->em*seg->em*seg->z*seg->z);
			double sgmSegM = seg->em;
			double sgmVtxR = std::sqrt(msSgmR*msSgmR+sgmSegR*sgmSegR);
			double sgmVtxM = std::sqrt(msSgmM*msSgmM+sgmSegM*sgmSegM);

			double vtxDR = (segR0  - ((seg->d==0) ? trX0   : trY0));
			double vtxDM = (seg->m - ((seg->d==0) ? trTanX : trTanY));

			double vtxZ    = -(vtxDR / vtxDM);
			double sgmVtxZ = std::fabs(vtxZ) * std::sqrt((sgmVtxR*sgmVtxR)/(vtxDR*vtxDR) + (sgmVtxM*sgmVtxM)/(vtxDM*vtxDM+1e-4));
			if (sgmVtxZ > 10.) sgmVtxZ = 10.;

			double lmtZu = (ZLimit[0] + sgmVtxZ * SgmLimit);
			double lmtZl = (ZLimit[1] - sgmVtxZ * SgmLimit);

			if (vtxZ > lmtZu || vtxZ < lmtZl) continue;

			numOfVtx[seg->d]++;
			if (seg->d==0) vtxCooXZ.push_back(std::make_pair(vtxZ, sgmVtxZ));
			else           vtxCooYZ.push_back(std::make_pair(vtxZ, sgmVtxZ));
		}

		fTrd.numOfHSegVtx[0] = numOfVtx[0];
		fTrd.numOfHSegVtx[1] = numOfVtx[1];
	}


	// TrdKCluster
	float TOF_Beta = 1;
	if      (recEv.iBetaH >= 0) TOF_Beta = std::fabs(event->pBetaH(recEv.iBetaH)->GetBeta());
	else if (recEv.iBeta  >= 0) TOF_Beta = std::fabs(event->pBeta(recEv.iBeta)->Beta);
	
	TrdKCluster * trdkcls = TrdKCluster::gethead();
	for (int kindOfFit = 0; kindOfFit <= 1; ++kindOfFit) {
		if (!(checkEventMode(EventBase::BT) ||
					(trdkcls->IsReadAlignmentOK == 2 && trdkcls->IsReadCalibOK == 1))) break;
		bool isOK = false;
		switch(kindOfFit) {
			case 0 :
				{
				  if      (recEv.iTrdHTrack < 0 && recEv.iTrdTrack  < 0) break;
				  else if (recEv.iTrdHTrack >= 0) {
				  	TrdHTrackR * trdh = event->pTrdHTrack(recEv.iTrdHTrack);
				  	trdkcls->Build(trdh);
				  }
				  else if (recEv.iTrdTrack  >= 0) {
				  	TrdTrackR * trd = event->pTrdTrack(recEv.iTrdTrack);
				  	trdkcls->Build(trd);
				  }
					else break;
				  isOK = true;
				}
				break;
			case 1 :
				{
				  if (recEv.iTrTrack < 0) break;
				  TrTrackR * trtk = event->pTrTrack(recEv.iTrTrack);
				  int fitid_max = trtk->iTrTrackPar(1, 0, 21);
				  if (fitid_max < 0) break;
				  trdkcls->SetTrTrack(trtk, fitid_max);
				  isOK = true;
				}
				break;
			default :
				break;
		}
		if (!isOK) continue;

		// It speeds lots of time 0.02 sec
		float Q = -1;
		float Qerror = -1;
		int   QnumberOfHit = -1;
		int   Qstatus = trdkcls->CalculateTRDCharge(0, TOF_Beta);
		if (Qstatus >= 0) {
			Q = trdkcls->GetTRDCharge();
			Qerror = trdkcls->GetTRDChargeError();
			QnumberOfHit = trdkcls->GetQNHit();
		}
		else continue;

		int nhits = 0; //To be filled with number of hits taken into account in Likelihood Calculation
		float threshold = 15; //ADC above which will be taken into account in Likelihood Calculation,  15 ADC is the recommended value for the moment.
		double llr[3] = {-1, -1, -1}; //To be filled with 3 LikelihoodRatio :  e/P, e/H, P/H
		if      (kindOfFit == 0) trdkcls->GetLikelihoodRatio_TRDRefit(threshold, llr, nhits);
		else if (kindOfFit == 1) trdkcls->GetLikelihoodRatio_TrTrack(threshold, llr, nhits);
		if (llr[0] < 0 || llr[1] < 0 || llr[2] < 0) continue;

		int numberOfHit = 0;
		for (int ih = 0; ih < nhits; ih++) {
			TrdKHit * hit = trdkcls->GetHit(ih);
			if (hit == nullptr) continue;
			if (checkEventMode(EventBase::ISS) || checkEventMode(EventBase::BT))
				if (!hit->IsCalibrated) continue;
			if (checkEventMode(EventBase::ISS))
				if (!hit->IsAligned) continue;
			//int lay = hit->TRDHit_Layer;
			//float amp = hit->TRDHit_Amp;
			numberOfHit++;
		}
		if (numberOfHit <= 0) continue;

		fTrd.statusKCls[kindOfFit] = true;
		fTrd.Q[kindOfFit]          = Q;
		fTrd.LLR[kindOfFit][0]     = llr[0];
		fTrd.LLR[kindOfFit][1]     = llr[1];
		fTrd.LLR[kindOfFit][2]     = llr[2];
		fTrd.LLR_nhit[kindOfFit]   = numberOfHit;
	}

	short trdIdx = -1;
	if      (recEv.iTrdHTrack >= 0) trdIdx = 1;
	else if (recEv.iTrdTrack  >= 0) trdIdx = 0;

	while (trdIdx >= 0) {
		bool isOK = false;
		AMSPoint trd_coo;
		AMSDir trd_dir;
		// Base on TrdTrack
		if (trdIdx == 0) {
			TrdTrackR * trd = event->pTrdTrack(recEv.iTrdTrack);
			trd_coo.setp(trd->Coo[0], trd->Coo[1], trd->Coo[2]);
			trd_dir.SetTheta(trd->Theta);
			trd_dir.SetPhi(trd->Phi);
			trd_dir[0] *= -1; trd_dir[1] *= -1; trd_dir[2] *= -1;
			isOK = true;
		}
		// Base on TrdHTrack
		if (trdIdx == 1) {
			TrdHTrackR * trdh = event->pTrdHTrack(recEv.iTrdHTrack);
			trd_coo.setp(trdh->Coo[0], trdh->Coo[1], trdh->Coo[2]);
			trd_dir.setd(trdh->Dir[0], trdh->Dir[1], trdh->Dir[2]);
			trd_dir[0] *= -1; trd_dir[1] *= -1; trd_dir[2] *= -1;
			isOK = true;
		}
		if (!isOK) break;
		else isOK = true;

		fTrd.trackStatus = true;
		fTrd.trackState[0] = trd_coo[0];
		fTrd.trackState[1] = trd_coo[1];
		fTrd.trackState[2] = trd_coo[2];
		fTrd.trackState[3] = trd_dir[0];
		fTrd.trackState[4] = trd_dir[1];
		fTrd.trackState[5] = trd_dir[2];

		break;
	} // while loop --- trd track

	fStopwatch.stop();
	return selectEvent(event);
}

bool EventTrd::selectEvent(AMSEventR * event) {
	return true;
}


//---- EventRich ----//
EventRich::EventRich() : EventBase() {
}

EventRich::~EventRich() {
}

void EventRich::initEvent() {
	fRich.init();
}

void EventRich::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventRich::setEventTree()\n";
#endif

	if (evTree == nullptr) setTree("EventRich", "Ring Imaging Cherenkov");
	else tree = evTree;

	tree->Branch("rich", &fRich);
}

void EventRich::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventRich::setEnvironment()\n";
#endif

	if (checkEventMode(EventBase::ISS) || checkEventMode(EventBase::MC)) {
		RichRingR::setBetaCorrection(RichRingR::fullUniformityCorrection);
	}
	if (checkEventMode(EventBase::BT)) {
		RichRingR::setBetaCorrection(RichRingR::noCorrection);
	}
		
	TString RichDefaultAGLTables = Form("%s/v5.00/RichDefaultAGLTables.04.dat",getenv("AMSDataDir"));
	RichOffline::RichRadiatorTileManager::Init((char*)RichDefaultAGLTables.Data());
}

bool EventRich::processEvent(AMSEventR * event, AMSChain * chain) {
	initEvent();
	if (event == nullptr)	return false;
	fStopwatch.start();
		
	// Rich cuts from Javie/Jorge
  const float richRadZ[2] = {-73.65, -74.65};          // aerogel / NaF
	const float cut_prob = 0.01;                         // Kolmogorov test prob
	const float cut_pmt = 3;                             // number of PMTs
	const float cut_collPhe[2] = {0.6, 0.4};             // ring phe / total phe (NaF, Aerogel)
	const float cut_chargeConsistency = 5;               // hit / PMT charge consistency test
	const float cut_betaConsistency[2] = {0.01, 0.005};  // beta_lip vs beta_ciemat consistency (NaF, Aerogel)
	const float cut_expPhe[2] = {1, 2};                  // expected number of phe (NaF, Aerogel)
	const float cut_aerogelExternalBorder = 3350;        // aerogel external border (r**2)
	// modify 3500 -> 3360 (by S.H.)
	// modify 3500 -> 3350 (by HY.Chou)
	const float cut_PMTExternalBorder = 4050;             // PMT external border (r**2) (by HY.Chou)
	const float cut_aerogelNafBorder[2] = {17.25, 17.75}; // aerogel/NaF border (NaF, Aerogel)
	// modify (17, 18) -> (17.25, 17.75) (by HY.Chou)
	const float cut_distToTileBorder[2] = {0.08, 0.05};   // distance to tile border 
	
	const int nBadTile = 5;
	const int kBadTile_Offical[nBadTile]  = { 3, 7, 87, 100, 108 }; // tiles with bad beta reconstruction
	const int kBadTile_MgntTile[nBadTile] = { 13, 23, 58, 86, 91 }; // tiles with bad beta reconstruction
	
	//fRich.numOfRing = event->NRichRing();
	//fRich.numOfHit = event->NRichHit();

	// RichVeto - start
	while (recEv.iTrTrack >= 0) {
		TrTrackR * trtk = event->pTrTrack(recEv.iTrTrack);
		if (trtk == nullptr) break;
		int fitid = trtk->iTrTrackPar(1, 3, 21);
		if (fitid < 0) break;

		AMSPoint ems_coo;
		AMSDir   ems_dir;
		int    tmp_kindOfRad = 0;
		for (int it = 0; it < 2; ++it) {
			trtk->Interpolate(richRadZ[tmp_kindOfRad], ems_coo, ems_dir, fitid);
			RichOffline::TrTrack tmp_track(ems_coo, ems_dir);
			RichOffline::RichRadiatorTileManager tmp_mgntTile(&tmp_track);
			tmp_kindOfRad = tmp_mgntTile.getkind() - 1;
			if (tmp_kindOfRad < 0 || tmp_kindOfRad > 1) break;
		}
		if (tmp_kindOfRad < 0 || tmp_kindOfRad > 1) break;
		
		float ems[6] = {0};
		ems[0] =  ems_coo[0]; ems[1] =  ems_coo[1]; ems[2] =  ems_coo[2];
		ems[3] = -ems_dir[0]; ems[4] = -ems_dir[1]; ems[5] = -ems_dir[2];

		// Radiator tile manager
		RichOffline::TrTrack track(ems_coo, ems_dir);
		RichOffline::RichRadiatorTileManager mgntTile(&track);
	
		if (mgntTile.getkind() < 1 || mgntTile.getkind() > 2) break;

		short kindOfRad    = mgntTile.getkind() - 1;
		short tileOfRad    = mgntTile.getcurrenttile();
		float rfrIndex     = mgntTile.getindex();
		AMSPoint richems   = mgntTile.getemissionpoint();
		float distToBorder = mgntTile.getdistance();
			
		int countKB = 0;
		for (int kb = 0; kb < nBadTile; kb++) {
			if (tileOfRad == kBadTile_MgntTile[kb]) break;
			countKB++;
		}
		bool isGoodTile = (countKB == nBadTile);

		// geometry cut
		bool isStruct = false;
		if ((kindOfRad == 1 && (std::max(std::fabs(richems[0]), std::fabs(richems[1])) > cut_aerogelNafBorder[0])) ||
				(kindOfRad == 0 && (std::max(std::fabs(richems[0]), std::fabs(richems[1])) < cut_aerogelNafBorder[1] ||
					(richems[0]*richems[0]+richems[1]*richems[1]) > cut_aerogelExternalBorder))) isStruct = true;
		//if (distToBorder < cut_distToTileBorder[kindOfRad] * std::fabs(1./ems[5])) isStruct = true;
		bool isInFiducialVolume = !isStruct;
	
		// Number of photoelectrons expected for a given track, beta and charge.
		const int    npart = 5;
		const bool   openCal[npart] = { 1, 0, 0, 1, 0 };
		const double chrg[npart] = { 1., 1., 1., 1., 1. };
		const double mass[npart] = { 0.000510999,   // electron
		                               0.1395701835,  // pion
																	 0.493667,      // kaon
																	 0.938272297,   // proton
																	 1.876123915 }; // deuterium
		float numOfExpPE[npart] = {0};
		std::fill_n(numOfExpPE, npart, -1);
		double rigAbs = std::fabs(trtk->GetRigidity(fitid));
		for (int i = 0; i < npart; ++i) {
			if (!openCal[i]) continue;
			double massChrg = mass[i] / chrg[i];
			double beta = 1. / std::sqrt((massChrg * massChrg / rigAbs / rigAbs) + 1); 
			double exppe = RichRingR::ComputeNpExp(trtk, beta, chrg[i]);
			numOfExpPE[i] = exppe;
		}
	
		fRich.kindOfRad = kindOfRad;
		fRich.tileOfRad = tileOfRad;
		fRich.rfrIndex  = rfrIndex;
		std::copy(ems, ems+6, fRich.emission);
		fRich.distToBorder = distToBorder;
		std::copy(numOfExpPE, numOfExpPE+npart, fRich.numOfExpPE);
		fRich.isGoodTile = isGoodTile; 
		fRich.isInFiducialVolume = isInFiducialVolume;

		break;
	}
	// RichVeto - end
	
	// official RichRingR - start
	bool isSuccRing = false;
	while (recEv.iRichRing >= 0 && fRich.kindOfRad >= 0) {
		// RichRingR
		RichRingR * rich = event->pRichRing(recEv.iRichRing);
		int kindOfRad = rich->IsNaF() ? 1 : 0;
		int tileOfRad = rich->getTileIndex();
		//float diffDist = std::fabs(rich->DistanceTileBorder() - fRich.distToBorder);
		if (fRich.kindOfRad != kindOfRad) break;
		//if (diffDist > cut_distToTileBorder[kindOfRad]) break;

		fRich.status = true;
		fRich.beta = rich->getBeta();
		fRich.Q = rich->getCharge2Estimate(true);
		fRich.Q = (fRich.Q > 1.0e-3) ? std::sqrt(fRich.Q) : -1;

		fRich.isGoodRecon = true;
		if (!rich->IsGood() || !rich->IsClean() ||
				rich->getProb() < cut_prob ||
				rich->getPMTs() < cut_pmt ||
				rich->getPMTChargeConsistency() > cut_chargeConsistency)
			fRich.isGoodRecon = false;

	 if (rich->IsNaF()) {
			if (rich->getExpectedPhotoelectrons() < cut_expPhe[0] ||
					rich->getBetaConsistency() > cut_betaConsistency[0] ||
					rich->getPhotoElectrons()/RichHitR::getCollectedPhotoElectrons() < cut_collPhe[0])
				fRich.isGoodRecon = false;
		}
		else {
			if (rich->getExpectedPhotoelectrons() < cut_expPhe[1] ||
					rich->getBetaConsistency() > cut_betaConsistency[1] ||
					rich->getPhotoElectrons()/RichHitR::getCollectedPhotoElectrons() < cut_collPhe[1])
				fRich.isGoodRecon = false;
		}

		isSuccRing = true;
		break;
	}
	// official RichRingR - end


	// Rich Hits
	short numOfPrimHit[2] = {0, 0};
	float numOfPrimPE[2]  = {0, 0}; 
	short numOfOthHit[2] = {0, 0};
	float numOfOthPE[2]  = {0, 0};
	for (int it = 0; it < event->NRichHit(); ++it) {
		RichHitR * hit = event->pRichHit(it);
		if (hit == nullptr) continue;
		bool used[2] = { false, false };
		for (int iring = 0; iring < event->NRichRing(); ++iring) {
			bool isUsed = hit->UsedInRingNumber(iring);
			if (!isUsed) continue;
			if (iring == recEv.iRichRing && isSuccRing) used[0] = true;
			else                                        used[1] = true;
		}
		bool isOthers = (!used[0]);
		
		short cross = (hit->IsCrossed() ? 0 : 1);
		float npe   = hit->Npe;

		if (used[0])  { numOfPrimHit[cross]++;  numOfPrimPE[cross] += npe;  }
		if (isOthers) { numOfOthHit[cross]++; numOfOthPE[cross] += npe; }
	}

	fRich.numOfPrimHit[0] = numOfPrimHit[0];
	fRich.numOfPrimHit[1] = numOfPrimHit[1];
	fRich.numOfPrimPE[0]  = numOfPrimPE[0];
	fRich.numOfPrimPE[1]  = numOfPrimPE[1];
	
	fRich.numOfOthHit[0] = numOfOthHit[0];
	fRich.numOfOthHit[1] = numOfOthHit[1];
	fRich.numOfOthPE[0]  = numOfOthPE[0];
	fRich.numOfOthPE[1]  = numOfOthPE[1];

	fStopwatch.stop();
	return selectEvent(event);
}

bool EventRich::selectEvent(AMSEventR * event) {
	return true;
}

//---- EventEcal ----//
EventEcal::EventEcal() : EventBase() {
}

EventEcal::~EventEcal() {
}

void EventEcal::initEvent() {
	fEcal.init();
}

void EventEcal::setEventTree(TTree * evTree) {
#if Debug == true
	std::cerr << "Debug : Now, EventEcal::setEventTree()\n";
#endif

	if (evTree == nullptr) setTree("EventEcal", "Electromagnetic Calorimeter");
	else tree = evTree;

	tree->Branch("ecal", &fEcal);
}

void EventEcal::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, EventEcal::setEnvironment()\n";
#endif

	//EcalShowerR::enableAutomaticRecoveryOfBrokenBuilds = true;
	EcalShowerR::enableAutomaticRecoveryOfBrokenBuilds = false;
	EcalHadron::InitParameters();
}

bool EventEcal::processEvent(AMSEventR * event, AMSChain * chain) {
	initEvent();
	if (event == nullptr)	return false;
	fStopwatch.start();

	// threshold : remove MIPs and low energy particles (lepton study 0.250)
	// threshold : remove very low energy particles (proton study 0.050)
	float threshold = 0.050;

	fEcal.numOfShower = event->NEcalShower();
	
	while (recEv.iEcalShower >= 0) {
		
		for (int ish = 0; ish < event->NEcalShower()+1; ++ish) {
			int jsh = ish - 1;
			if (jsh == recEv.iEcalShower) continue;
			EcalShowerR * ecal = event->pEcalShower( ((jsh==-1)?recEv.iEcalShower:jsh) );
			if (ecal == nullptr) continue;
		
			ShowerInfo shower;
		
			shower.energyD = 1.e-3 * ecal->EnergyD;
			float energyA = ecal->EnergyA;
			shower.energyE = ecal->EnergyE;
			float energyC = ecal->EnergyC;
			shower.energyP = ecal->EnergyP(2);
			if (energyC > 1.0e10 || shower.energyP < 0) {
				shower.energyP = -1;
			}
			
			shower.PisaBDT = ecal->GetEcalBDT();
			
			// Charge estimator based on ECAl-only;
			// it is advised to discard low-rigidity events
			// and to use events having only one EcalShower.
			shower.Q = ecal->EcalChargeEstimator();
			if (shower.Q < 1e-3) shower.Q = -1;

			shower.showerAxis[0] = ecal->CofG[0];
			shower.showerAxis[1] = ecal->CofG[1];
			shower.showerAxis[2] = ecal->CofG[2];
			shower.showerAxis[3] = ecal->Dir[0];
			shower.showerAxis[4] = ecal->Dir[1];
			shower.showerAxis[5] = ecal->Dir[2];

			// hadron shower
			int zeh = 1;
			if (jsh == -1) {
				if (recEv.iTrTrack >= 0) {
					TrTrackR * trtk = event->pTrTrack(recEv.iTrTrack);
					std::vector<like_t> probLvl;
					trtk->GetInnerZ(probLvl);
    		  if (probLvl.size() != 0)
    		    zeh = probLvl.at(0).Z;
				}
				else if (recEv.iBetaH >= 0 && event->pBetaH(recEv.iBetaH) != 0) {
					BetaHR * betaH = event->pBetaH(recEv.iBetaH);
					int nlay = 0; float prob = 0;
					TofChargeHR tofQH = betaH->gTofCharge();
					zeh = tofQH.GetZ(nlay, prob, 0);
				}
			}

			EcalHadron::Build(ecal, zeh);
			shower.hadronApex = EcalHadron::EcalApex;
			shower.hadronEnergy = EcalHadron::EcalRigidity;

			fEcal.showers.push_back(shower);
		}

		break;
	} // while loop --- iEcalShower

	/* Raw Hits
	for (int ih = 0; ih < event->NEcalHit(); ++ih) {
		EcalHitR * ecalHit = event->pEcalHit(ih);
		if (ecalHit == nullptr) continue;

		HitECALInfo hit;
		hit.id = ecalHit->Plane * 100 + ecalHit->Cell;
		hit.side = ecalHit->Proj;
		hit.edep = ecalHit->Edep;
		hit.coo[0] = ecalHit->Coo[0];
		hit.coo[1] = ecalHit->Coo[1];
		hit.coo[2] = ecalHit->Coo[2];
		fEcal.rawHits.push_back(hit);
	}
	std::sort(fEcal.rawHits.begin(), fEcal.rawHits.end(), HitECALInfo_sort());
	*/
	
	fStopwatch.stop();
	return selectEvent(event);
}

bool EventEcal::selectEvent(AMSEventR * event) {
	return true;
}


//---- DataSelection ----//
DataSelection::SWITCH DataSelection::option[DataSelection::NUMBER] = {
	DataSelection::ON
};

DataSelection::DataSelection() {
#if Debug == true
	std::cerr << "Debug : Now, DataSelection::DataSelection()\n";
#endif

	isMultiTree = false;
	evTree = 0;
}

DataSelection::~DataSelection() {
#if Debug == true
	std::cerr << "Debug : Now, DataSelection::~DataSelection()\n";
#endif

	if (evTree != 0) delete evTree;
	evTree = 0;
}

void DataSelection::setMultiTree(bool isMTree) {
#if Debug == true
	std::cerr << "Debug : Now, DataSelection::setMultiTree()\n";
#endif

	isMultiTree = isMTree;
}

void DataSelection::setEventTree() {
#if Debug == true
	std::cerr << "Debug : Now, DataSelection::setEventTree()\n";
#endif

	if (isMultiTree == false) evTree = new TTree("data", "AMS data");
	else evTree = 0;

	if (checkOption(DataSelection::LIST)) list.setEventTree(evTree);
	if (checkOption(DataSelection::RTI)) rti.setEventTree(evTree);
	if (checkOption(DataSelection::TRG)) trg.setEventTree(evTree);
	if (checkOption(DataSelection::TOF)) tof.setEventTree(evTree);
	if (checkOption(DataSelection::ACC)) acc.setEventTree(evTree);
	if (checkOption(DataSelection::TRK)) trk.setEventTree(evTree);
	if (checkOption(DataSelection::TRD)) trd.setEventTree(evTree);
	if (checkOption(DataSelection::RICH)) rich.setEventTree(evTree);
	if (checkOption(DataSelection::ECAL)) ecal.setEventTree(evTree);
}

void DataSelection::setEnvironment() {
#if Debug == true
	std::cerr << "Debug : Now, DataSelection::setEnvironment()\n";
#endif

	list.setEnvironment();
	rti.setEnvironment();
	trg.setEnvironment();
	tof.setEnvironment();
	acc.setEnvironment();
	trk.setEnvironment();
	trd.setEnvironment();
	rich.setEnvironment();
	ecal.setEnvironment();

  // set scale function
  DataSelection::gScaleFact = 0.01;
  DataSelection::gScaleFunc1D.SetParameter(0, DataSelection::gScaleFact);
  DataSelection::gScaleFunc2D.SetParameter(0, DataSelection::gScaleFact);
}

int DataSelection::processEvent(AMSEventR * event, AMSChain * chain) {
#if Debug == true
	//std::cerr << "Debug : Now, DataSelection::processEvent()\n";
#endif

	if (event == nullptr) return -1;

	bool statusList = true;
	//bool statusRti = true;
	bool statusTrg = true;
	bool statusTof = true;
	bool statusAcc = true;
	bool statusTrk = true;
	bool statusTrd = true;
	bool statusRich = true;
	bool statusEcal = true;

	MGClock::HrsStopwatch fStopwatch;
	fStopwatch.start();

	if (checkOption(DataSelection::LIST)) statusList = list.processEvent(event, chain);
	if (!statusList) return -101;

	//if (checkOption(DataSelection::RTI)) statusRti = rti.processEvent(event);
	//if (!statusRti) return -102;

	if (checkOption(DataSelection::TRG)) statusTrg = trg.processEvent(event);
	if (!statusTrg) return -103;

	if (checkOption(DataSelection::ACC)) statusAcc = acc.processEvent(event);
	if (!statusAcc) return -105;

	if (checkOption(DataSelection::TRK)) statusTrk = trk.processEvent(event);
	if (!statusTrk) return -106;

	if (checkOption(DataSelection::TRD)) statusTrd = trd.processEvent(event);
	if (!statusTrd) return -107;

	if (checkOption(DataSelection::RICH)) statusRich = rich.processEvent(event);
	if (!statusRich) return -108;

	if (checkOption(DataSelection::ECAL)) statusEcal = ecal.processEvent(event);
	if (!statusEcal) return -109;

	// due to rebuild BetaH (put on last one)
	if (checkOption(DataSelection::TOF)) statusTof = tof.processEvent(event);
	if (!statusTof) return -104;

	fStopwatch.stop();

#if Debug == true
	const float limitfStopwatch = 5.0;
	float totlTime = fStopwatch.time() + recEv.time() + rti.time();
	if (totlTime > limitfStopwatch) {
		COUT("\nRUN %u  EVENT %u\n", event->Run(), event->Event());
		COUT("REAL TIME : %14.8f (SEC)   100.00%\n", totlTime);
		COUT("    RECON   %14.8f (SEC)   %6.2f%\n", recEv.time(), recEv.time() / totlTime * 100);
		COUT("     LIST   %14.8f (SEC)   %6.2f%\n", list.time(),  list.time() / totlTime * 100);
		COUT("      RTI   %14.8f (SEC)   %6.2f%\n", rti.time(),   rti.time() / totlTime * 100);
		COUT("      TRG   %14.8f (SEC)   %6.2f%\n", trg.time(),   trg.time() / totlTime * 100);
		COUT("      TOF   %14.8f (SEC)   %6.2f%\n", tof.time(),   tof.time() / totlTime * 100);
		COUT("      ACC   %14.8f (SEC)   %6.2f%\n", acc.time(),   acc.time() / totlTime * 100);
		COUT("      TRK   %14.8f (SEC)   %6.2f%\n", trk.time(),   trk.time() / totlTime * 100);
		COUT("      TRD   %14.8f (SEC)   %6.2f%\n", trd.time(),   trd.time() / totlTime * 100);
		COUT("     RICH   %14.8f (SEC)   %6.2f%\n", rich.time(),  rich.time() / totlTime * 100);
		COUT("     ECAL   %14.8f (SEC)   %6.2f%\n", ecal.time(),  ecal.time() / totlTime * 100);
	}
#endif

	return 0;
}

void DataSelection::fill() {
	if (isMultiTree == false) evTree->Fill();
	else {
		if (checkOption(DataSelection::LIST)) list.fill();
		if (checkOption(DataSelection::RTI)) rti.fill();
		if (checkOption(DataSelection::TRG)) trg.fill();
		if (checkOption(DataSelection::TOF)) tof.fill();
		if (checkOption(DataSelection::ACC)) acc.fill();
		if (checkOption(DataSelection::TRK)) trk.fill();
		if (checkOption(DataSelection::TRD)) trd.fill();
		if (checkOption(DataSelection::RICH)) rich.fill();
		if (checkOption(DataSelection::ECAL)) ecal.fill();
	}
}

int DataSelection::preselectEvent(AMSEventR * event, const std::string& officialDir) {
	if (event == nullptr)	return -1;
	EventList::Weight = 1.;

	// Resolution tuning
	if (EventBase::checkEventMode(EventBase::MC)) {
		event->SetDefaultMCTuningParameters();
		TrExtAlignDB::SmearExtAlign();
		TRCLFFKEY.UseSensorAlign = 0;
	}

	//-----------------------------//
	//----  Fast Preselection  ----//
	//-----------------------------//
	//COUT("============= <Run %u Event %u> =============\n", event->Run(), event->Event());

	// ~1~ (Based on BetaH(Beta))
	bool isDownBeta = false;
	for (int ibta = 0; ibta < event->NBeta(); ++ibta)
		if (event->pBeta(ibta)->Beta > 1.0e-2) { isDownBeta = true; break; }

	TofRecH::BuildOpt = 0; // normal
	bool isDownBetaH = false;
	for (int ibta = 0; ibta < event->NBetaH(); ++ibta)
		if (event->pBetaH(ibta)->GetBeta() > 1.0e-2) { isDownBetaH = true; break; }
	
	if (!isDownBeta && !isDownBetaH) return -1001;
	
	// ~2~ (Based on TrTrack)
	if (event->NTrTrack() != 1) return -2001;
	
	// ~3~ (Based on TrdTrack)
	if (event->NTrdTrack() == 0 && event->NTrdHTrack() == 0) return -3001;

	// ~4~ (Based on Particle)
	ParticleR  * partSIG = (event->NParticle() > 0) ? event->pParticle(0) : nullptr;
	TrTrackR   * trtkSIG = (partSIG != nullptr) ? partSIG->pTrTrack() : nullptr;
	BetaHR     * btahSIG = (partSIG != nullptr) ? partSIG->pBetaH()   : nullptr;
	TrdTrackR  * trdSIG  = (partSIG != nullptr) ? partSIG->pTrdTrack()  : nullptr;
	TrdHTrackR * trdhSIG = (partSIG != nullptr) ? partSIG->pTrdHTrack() : nullptr;
	if (partSIG == nullptr) return -4001;
	if (trtkSIG == nullptr || btahSIG == nullptr) return -4002;
	if ( trdSIG == nullptr && trdhSIG == nullptr) return -4003;

	// ~5~ (Based on BetaH)
	if (btahSIG->GetBetaPattern() != 4444) return -5001;
	if (btahSIG->GetBeta() < 0) return -5002;

	// ~6~ (Based on Track Hits)
	const unsigned short _hasTrL34 =  12;
	const unsigned short _hasTrL56 =  48;
	const unsigned short _hasTrL78 = 192;
	unsigned short trBitPattJ = trtkSIG->GetBitPatternJ();
	bool     isTrInner  = ((trBitPattJ&_hasTrL34) > 0 && 
	                       (trBitPattJ&_hasTrL56) > 0 && 
												 (trBitPattJ&_hasTrL78) > 0);
	if (!isTrInner) return -6001;
	//unsigned short trBitPattXYJ = trtkSIG->GetBitPatternXYJ();
	//bool     isTrInnerXY  = ((trBitPattXYJ&_hasTrL34) > 0 && 
	//                         (trBitPattXYJ&_hasTrL56) > 0 && 
	//											   (trBitPattXYJ&_hasTrL78) > 0);

	//if (!isTrInnerXY) return -6001;

	// ~7~ (Only for Antiproton to Proton Flux Ratio Study)
	bool isAppStudy = true;
	if (isAppStudy) {
		//isAppStudy = EventBase::checkEventMode(EventBase::ISS);
		isAppStudy = ( EventBase::checkEventMode(EventBase::ISS) || (officialDir.find("pr.pl1.flux.l1a9.2016000") != std::string::npos) || (officialDir.find("pr.pl1.flux.l1o9.2016000") != std::string::npos) );
	}

	if (isAppStudy) {
		bool   hasTr    = false;
		bool   isScale  = true;
		double sclRig = 0.0;
		const int npatt = 4;
		const int spatt[npatt] = { 3, 5, 6, 7 };
		for (int ip = 0; ip < npatt; ++ip) {
			int fitid = trtkSIG->iTrTrackPar(1, spatt[ip], 21);
			if (fitid < 0) continue;
			hasTr = true;
			sclRig = trtkSIG->GetRigidity(fitid);
			if (sclRig < 0.0) { isScale = false; break; }
		}
		
		int fitidInn = trtkSIG->iTrTrackPar(1, 3, 21);
		double trQIn = (fitidInn>=0) ? trtkSIG->GetInnerQ_all(1.0, fitidInn).Mean : 1.0;

		if (hasTr && isScale) {
			//double scaleProb = gScaleFunc1D.Eval(sclRig); // (by Rigidity)
			double scaleProb = gScaleFunc2D.Eval(sclRig, trQIn); // (by Rigidity & Charge)
    	if (MGNumc::Compare(MGRndm::DecimalUniform(), scaleProb) > 0) return -7001;
    	else EventList::Weight *= (1. / scaleProb);
		}
	}
	
	// ~8~ (Based on RTI)
	if (EventBase::checkEventMode(EventBase::ISS) && checkOption(DataSelection::RTI)) {
		if (!rti.processEvent(event)) return -8001;

		double minStormer = *std::min_element(rti.fRti.cutoffStormer, rti.fRti.cutoffStormer+4);
		double minIGRF    = *std::min_element(rti.fRti.cutoffIGRF, rti.fRti.cutoffIGRF+4);
		double minCf      =  std::min(minStormer, minIGRF);
		double maxRig     = 0.0;

		bool hasTr = false;
		const int npatt = 4;
		const int spatt[npatt] = { 3, 5, 6, 7 };
		for (int ip = 0; ip < npatt; ++ip) {
			int fitid = trtkSIG->iTrTrackPar(1, spatt[ip], 21);
			if (fitid < 0) continue;
			hasTr = true;
			maxRig = std::max(maxRig, std::fabs(trtkSIG->GetRigidity(fitid)));
		}

		const double minFact = 1.2;
		if ( hasTr && (maxRig < (minFact * minCf)) ) return -8002;
	}

	//--------------------------//
	//----  Reconstruction  ----//
	//--------------------------//
	if (!recEv.rebuild(event)) return -9999;

	return 0;
}

int DataSelection::selectEvent(AMSEventR * event) {
	if (event == nullptr) return -1;

	// User Define

	return 0;
}

int DataSelection::analysisEvent(AMSEventR * event) {
	if (event == nullptr) return -1;

	// User Define

	return 0;
}


//---- RunTagOperator ----//
RunTagOperator::RunTagOperator() { init(); }

RunTagOperator::~RunTagOperator() { init(); }

void RunTagOperator::init() {
	fRunTag.clear();
}

bool RunTagOperator::processEvent(AMSEventR * event, AMSChain * chain) {
	if (chain == nullptr) return false;
	TFile * file = chain->GetFile();
	std::string filePath = file->GetName();
	if (event == nullptr) return false;
	UInt_t runID   = event->Run();
	UInt_t eventID = event->Event();
		
	std::map<UInt_t, RunTagInfo>::iterator it = fRunTag.find(runID);
	if (it == fRunTag.end()) {
		RunTagInfo info;
		info.run = runID;
		info.eventFT = eventID;
		info.eventLT = eventID;
		info.numOfSelEvent = 1;
		info.file.push_back(filePath);
		
		if (EventBase::checkEventMode(EventBase::ISS)) {
			MGClock::TTime * ttime = MGClock::ConvertFromUTimeToTTime(runID, MGClock::ClockType::UTC);
			info.dateUTC = (ttime->tm_year + 1900) * 10000 + (ttime->tm_mon+1) * 100 + (ttime->tm_mday);
			info.timeUTC = (ttime->tm_hour) * 10000 + (ttime->tm_min) * 100 + (ttime->tm_sec);
		}
		
		fRunTag[runID] = info;
	}
	else {
		RunTagInfo & info = it->second;
		if (info.eventFT > eventID) info.eventFT = eventID;
		if (info.eventLT < eventID) info.eventLT = eventID;
		if (std::find(info.file.begin(), info.file.end(), filePath) == info.file.end()) info.file.push_back(filePath);
		info.numOfSelEvent++;
	}
	return true;
}

void RunTagOperator::save(TFile * file) {
	if (file == nullptr) return;
	file->cd();
	TTree * tree = new TTree("runTag", "RunTag information");
	RunTagInfo info;
	tree->Branch("runTag", &info);
	for (std::map<UInt_t, RunTagInfo>::iterator it = fRunTag.begin(); it != fRunTag.end(); ++it) {
		info = it->second;
		UInt_t nTrgEv = UInt_t( (info.numOfSelEvent == 1) ? 2 * info.eventFT : ((double(info.numOfSelEvent) / double(info.numOfSelEvent - 1)) * double(info.eventLT - info.eventFT)) );
		info.numOfTrgEvent = nTrgEv;
		tree->Fill();
	}
	file->cd();
}


//---- YiNtuple ----//
YiNtuple::MODE YiNtuple::selectionMode = YiNtuple::NORM;

YiNtuple::YiNtuple() {
#if Debug == true
	std::cerr << "Debug : Now, YiNtuple::YiNtuple()\n";
#endif
	fStopwatch.start();
}

YiNtuple::~YiNtuple() {
#if Debug == true
	std::cerr << "Debug : Now, YiNtuple::~YiNtuple()\n";
#endif
	init();

	fStopwatch.stop();
	fStopwatch.print();
}

inline void YiNtuple::init() {
#if Debug == true
	std::cerr << "Debug : Now, YiNtuple::init()\n";
#endif

	fGroup.first = 0;
	fGroup.second = -1;
	fFileList.clear();
	fFileDir = "";
	if (fChain != 0) delete fChain;
	fChain = 0;
	fFileName = "";
}

void YiNtuple::setOutputFile(const std::string& file_name, const std::string& path, bool isMultiTree) {
#if Debug == true
	std::cerr << "Debug : Now, YiNtuple::setOutputFile()\n";
#endif

	fFileName = std::string(path) + "/" + std::string(file_name);
	fData = new DataSelection();
	fData->setMultiTree(isMultiTree);
	fRunTagOp = new RunTagOperator;
}

void YiNtuple::readDataFrom(const std::string& file_list, Long64_t group_id, Long64_t group_size) {
#if Debug == true
	std::cerr << "Debug : Now, YiNtuple::readDataFrom()\n";
#endif

	COUT("\n**--------------------------------------------**\n");
	COUT("\n**    Read Data Form Source File List Info    **\n");
	COUT("\n**--------------------------------------------**\n");

	// start check sourceFileList.txt
	std::vector<std::string>&& flist = MGIO::ReadFileContent(file_list);
	if (flist.size() == 0)
		MGSys::ShowErrorAndExit(LocAddr(), "ROOT file list cannot be opend! Exiting ...");
	// end check sourceFileList.txt

	// start load data with group
	if (group_id == 0 && group_size == -1) group_size = flist.size();
	if (group_size <= 0 || group_size > flist.size() || group_id < 0 || group_id >= flist.size())
		MGSys::ShowErrorAndExit(LocAddr(), "Group format has error(1)! Exiting ...");
	
	Long64_t begin = group_id * group_size;
	Long64_t end   = (group_id + 1) * group_size;
	if (begin >= 0 && begin < flist.size() && end > flist.size()) {
		end = flist.size();
	}
	else if (begin < 0 || begin >= flist.size() || end < 1 || end > flist.size())
		MGSys::ShowErrorAndExit(LocAddr(), "ERROR : Group format has error(2)! Exiting ...");

	fGroup = std::make_pair(group_id, group_size);
	for (int it = begin; it < end; it++) {
		fFileList.push_back(flist.at(it));
	}

	COUT("\n---- Loading Root Files ----\n");
	COUT("Group : %ld th   [%ld files/group],    Total of Load Files : %ld \n", fGroup.first, fGroup.second, fFileList.size());
	for (Long64_t  it = 0; it < fFileList.size(); it++) {
		COUT("    Number : %ld,   %s\n", it, fFileList.at(it).c_str());
	}

	if (fFileList.size() != 0) {
		std::vector<std::string> && strs = MGRegex::Split(fFileList.at(0), MGRegex::Formula::Slash);
		if (strs.size() > 2) fFileDir = strs.at(strs.size()-2);
	}
	// end load data with group

	// start read source file list
	bool stagedonly = true;
	unsigned int timeout = 10;
	fChain = new AMSChain("AMSRoot");
	int fileStatus = fChain->AddFromFile(file_list.c_str(), begin, end, stagedonly, timeout);
	if (fileStatus == -1)
		MGSys::ShowErrorAndExit(LocAddr(), "ROOT file list cannot be opend! Exiting ...");

	COUT("FileStatus : %d\n", fileStatus);
	COUT("Totally : %ld data events.\n", fChain->GetEntries());
	// end read source file list

	COUT("\n**-------------------------------------------**\n");
	COUT("\n**    Read Data Form Source File List End    **\n");
	COUT("\n**-------------------------------------------**\n");
}

void YiNtuple::saveInputFileList(TFile * file) {
	if (file == nullptr || fFileList.size() == 0) return;
	file->cd();
	TTree * tree = new TTree("fileList", "List Of Input File Info");
	std::string filePath;
	tree->Branch("file", &filePath);
	for (int i = 0; i < fFileList.size(); ++i) {
		filePath = fFileList.at(i);
		tree->Fill();
	}
	file->cd();
}

void YiNtuple::loopEventChain() {
	COUT("\n**-----------------------------**\n");
	COUT("\n**    Loop Event Chain Info    **\n");
	COUT("\n**-----------------------------**\n");

	TFile * file = 0;
	if (YiNtuple::checkSelectionMode(YiNtuple::NORM)) {
		file = new TFile(fFileName.c_str(), "RECREATE");
		fData->setEventTree();
	}
	else if (YiNtuple::checkSelectionMode(YiNtuple::COPY)) {
		fChain->OpenOutputFile(fFileName.c_str());
	}

	fData->setEnvironment(); // it must be before event loop. (before get event !)

	// check event type
	if (fChain->GetEntries() <= 0)
		MGSys::ShowErrorAndExit(LocAddr(), "ERROR : Don't have event! Exiting ...");

	AMSEventR * ev = fChain->GetEvent(0);

	if (ev != 0) {
		bool isMC  = ev->nMCEventg() > 0 &&
			EventBase::checkEventMode(EventBase::MC);
		bool isBT  = (ev->nMCEventg() <= 0) &&
			(ev->Run() < 1305795600) &&
			EventBase::checkEventMode(EventBase::BT);
		bool isISS = (ev->nMCEventg() <= 0) &&
			(ev->Run() > 1305795600) &&
			EventBase::checkEventMode(EventBase::ISS);
		if (!isMC && !isBT && !isISS)
			MGSys::ShowErrorAndExit(LocAddr(), "Event type (ISS, BT, MC) is failed! Exiting ...");
	}

	Long64_t loop_entries = fChain->GetEntries();
	Long64_t first_entry = 0;
	Long64_t last_entry = loop_entries;

	Long64_t npassed = 0;
	Long64_t nprocessed = 0;
	const Long64_t printLimit = 25000;
	Long64_t printRate = loop_entries / 100;
	if (printRate < printLimit) printRate = printLimit;
	if (printRate > printLimit * 5) printRate = printLimit * 5;
	for (Long64_t ientry = first_entry; ientry < last_entry; ++ientry){
		if (nprocessed%printRate == 0) {
			fStopwatch.stop();

			const unsigned int MemSize = 1024;
			ProcInfo_t procinfo;
			gSystem->GetProcInfo(&procinfo);
			Long64_t memRes = procinfo.fMemResident / MemSize;
			Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
			
			COUT("Info :: %lf %\n", 100. * float(nprocessed)/float(loop_entries));
			COUT("        Processed       : %ld / %ld\n", nprocessed, loop_entries);
			COUT("        Passed          : %ld / %ld\n", npassed, nprocessed);
			COUT("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
			COUT("        Real Time       : %9.2f (second)\n", fStopwatch.time());
			COUT("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fStopwatch.time());
			COUT("        Cpu    System   : %4.1f %\n", procinfo.fCpuSys);
			COUT("               User     : %4.1f %\n", procinfo.fCpuUser);
			COUT("        Memory Resident : %2ld GB %4ld MB\n", memRes / MemSize, memRes % MemSize);
			COUT("               Virtual  : %2ld GB %4ld MB\n", memVrl / MemSize, memVrl % MemSize);
		}
		nprocessed++;

		AMSEventR * event = fChain->GetEvent(ientry);

		//if (nprocessed > 10000) break; // testcode
		
		fRunTagOp->processEvent(event, fChain);

		int preselectEventStatus = fData->preselectEvent(event, fFileDir);
		if (preselectEventStatus < 0) continue;

		if (YiNtuple::checkSelectionMode(YiNtuple::NORM)) {
			int processEventStatus = fData->processEvent(event, fChain);
			if (processEventStatus < 0) continue;

			int selectEventStatus = fData->selectEvent(event);
			if (selectEventStatus < 0) continue;

			int analysisEventStatus = fData->analysisEvent(event);
			if (analysisEventStatus < 0) continue;

			fData->fill();
		}
		else if (YiNtuple::checkSelectionMode(YiNtuple::COPY)) {
			fChain->SaveCurrentEvent();
		}
		
		npassed++;
	}

	if (nprocessed == loop_entries) {
		fStopwatch.stop();
		const unsigned int MemSize = 1024;
		ProcInfo_t procinfo;
		gSystem->GetProcInfo(&procinfo);
		Long64_t memRes = procinfo.fMemResident / MemSize;
		Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
		
		COUT("Info :: %lf %\n", 100. * float(nprocessed)/float(loop_entries));
		COUT("        Processed       : %ld / %ld\n", nprocessed, loop_entries);
		COUT("        Passed          : %ld / %ld\n", npassed, nprocessed);
		COUT("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
		COUT("        Real Time       : %9.2f (second)\n", fStopwatch.time());
		COUT("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fStopwatch.time());
		COUT("        Cpu    System   : %4.1f %\n", procinfo.fCpuSys);
		COUT("               User     : %4.1f %\n", procinfo.fCpuUser);
		COUT("        Memory Resident : %2ld GB %4ld MB\n", memRes / MemSize, memRes % MemSize);
		COUT("               Virtual  : %2ld GB %4ld MB\n", memVrl / MemSize, memVrl % MemSize);
		COUT("Info :: AMSRoot Files Processed Successfully Finished.\n");
	}
	else {
		COUT("Info :: AMSRoot Files Processed Seems Failed.\n");
		COUT("        Processed %ld in %ld\n", nprocessed, loop_entries);
	}

	if (YiNtuple::checkSelectionMode(YiNtuple::NORM)) {
		fRunTagOp->save(file);
		saveInputFileList(file);
		
		file->cd();
		file->Write();
		file->Close();
		if (file) delete file;
	}
	else if (YiNtuple::checkSelectionMode(YiNtuple::COPY)) {
		fChain->CloseOutputFile();
	}

	COUT("\n**----------------------------**\n");
	COUT("\n**    Loop Event Chain End    **\n");
	COUT("\n**----------------------------**\n");
}
#endif // __YiProdNtuple_TCC__
