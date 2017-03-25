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
}

bool RecEvent::rebuild(AMSEventR * event) {
	if (event == 0) return false;
	timer.start();
	init();

	Int_t npar    = event->NParticle();
	Int_t nbeta   = event->NBeta();
	Int_t nbetaH  = event->NBetaH();
	Int_t ntrtk   = event->NTrTrack();
	Int_t nshower = event->NEcalShower();
	Int_t ntrdtk  = event->NTrdTrack();
	Int_t ntrdhtk = event->NTrdHTrack();
	Int_t nring   = event->NRichRing();


	// Known feature of the current track finding and new developments (works by Z.Qu)
	// By default all the TrTrack parameters are refitted
	//for (UInt_t i = 0; i < event->NTrTrack(); ++i) {
	//	TrTrackR * trtk = event->pTrTrack(i);
	//	if (trtk == 0) continue;
	//	trtk->AddLostHits();
	//}

	/** Particle **/
	Bool_t isParticle = false;
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
			for (Int_t it = 0; it < event->NBetaH() && iTrTrack>=0; ++it)
				if (event->pBetaH(it)->iTrTrack() == iTrTrack) iBetaH = it;

		isParticle = true;
	}
	if (!isParticle) { init(); timer.stop(); return false; }

	// Beta Information
	Float_t Beta = 0;
	if      (iBetaH >= 0) Beta = event->pBetaH(iBetaH)->GetBeta();
	else if (iBeta  >= 0) Beta = event->pBeta(iBeta)->Beta;
	else                  Beta = 1;

	// Tracker Information
	Float_t trackerZJ[9];
	if (EventBase::checkEventMode(EventBase::ISS))
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerAJ(layJ);
	else
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerJ(layJ);

	TrTrackR * TkStPar = 0;
	Int_t    TkStID = -1;
	AMSPoint TkStCoo;
	AMSDir   TkStDir;
	Float_t  TkStRig = 0;
	Float_t  TkStQ = 0;
	if (iTrTrack >= 0) {
		TkStPar = event->pTrTrack(iTrTrack);
		TkStID = TkStPar->iTrTrackPar(1, 3, 21);
		if (TkStID >= 0) {
			TkStCoo = TkStPar->GetP0(TkStID);
			TkStDir = TkStPar->GetDir(TkStID);
			TkStRig = TkStPar->GetRigidity(TkStID);
			TkStQ   = TkStPar->GetInnerQ_all(Beta, TkStID).Mean;
		}
		else { timer.stop(); return false; }
	}
		
	// ECAL Information
	// pre-selection (ECAL)
	if (TkStID >= 0 && iEcalShower >= 0) {
		EcalShowerR * ecal = event->pEcalShower(iEcalShower);
		AMSPoint EPnt(ecal->CofG);
		AMSPoint EPntP;
		AMSDir   EPntD;
		TkStPar->Interpolate(EPnt.z(), EPntP, EPntD, TkStID);
		
		Float_t dxPnt = std::fabs(EPntP.x() - EPnt.x());
		Float_t dyPnt = std::fabs(EPntP.y() - EPnt.y());

		Float_t lmtx = 5, lmty = 5; // 5 cm, 5 cm
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
		Float_t dxP = std::fabs(TrdEP.x() - TrdP.x());
		Float_t dyP = std::fabs(TrdEP.y() - TrdP.y());
		
		Float_t lmtx = 5, lmty = 5; // 5 cm, 5 cm
		if (dxP > lmtx || dyP > lmty) iTrdTrack = -1;
	}

	if (TkStID >= 0 && iTrdHTrack >= 0) {
		AMSPoint TrdHP(event->pTrdHTrack(iTrdHTrack)->Coo);
		AMSDir   TrdHD(event->pTrdHTrack(iTrdHTrack)->Dir);

		AMSPoint TrdHEP;
		AMSDir   TrdHED;
		TkStPar->Interpolate(TrdHP.z(), TrdHEP, TrdHED, TkStID);
		Float_t dxHP = std::fabs(TrdHEP.x() - TrdHP.x());
		Float_t dyHP = std::fabs(TrdHEP.y() - TrdHP.y());
		
		Float_t lmtx = 5, lmty = 5; // 5 cm, 5 cm
		if (dxHP > lmtx || dyHP > lmty) iTrdHTrack = -1;
	}

	if (isParticle) {
    MgntTrHit::Load(event);
		timer.stop();
		return true;
	}
	else init();
	
	timer.stop();
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
	if (isTreeSelf == false || tree == 0) return;
	tree->Fill();
}

inline void EventBase::setTree(const char * name, const char * title) {
	isTreeSelf = true;
	tree = new TTree(name, title);
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

	if (evTree == 0) setTree("EventList", "AMS-02 Event List");
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
	if (event == 0)	return false;
	timer.start();

	// LIST
	fList.runID   = event->Run();
	fList.eventID = event->Event();
	fList.entryID = chain->Entry();
	fList.weight  = EventList::gWeight;

	// G4MC
	if (checkEventMode(EventBase::MC)) {
		MCEventgR * primary = event->GetPrimaryMC();
		if (primary == 0) return false;

		// Tracker Information
		Float_t trackerZJ[9] = {0};
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerJ(layJ);

		fG4mc.beamID = primary->TBl + 1;
		fG4mc.primPart.trackID  = primary->trkID;
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

		std::set<Int_t> PrimSET;
		for (Int_t it = 0; it < event->NMCEventg(); it++) {
			MCEventgR * mcev = event->pMCEventg(it);
			if (mcev == 0) continue;
			if (mcev->parentID != 1) continue; // ONLY SAVE PRIM PART
			if (mcev->Momentum < 5.0e-2) continue;
			Double_t momRat = (mcev->Momentum / fG4mc.primPart.mom);
			if (mcev->parentID == primary->trkID && momRat > 5.0e-2) { PrimSET.insert(mcev->trkID); continue; }
			if (momRat < 0.368) continue;
			PrimSET.insert(mcev->trkID);
		}

		std::set<Int_t> TrSET;
		TrSET.insert(primary->trkID);
		for (Int_t iev = 0; iev < event->NMCEventg(); iev++) {
			MCEventgR * mcev = event->pMCEventg(iev);
			if (mcev == 0) continue;
			if (mcev->Momentum < 5.0e-2) continue;
			if (PrimSET.find(mcev->parentID) == PrimSET.end()) continue;
			TrSET.insert(mcev->trkID);
			
			PartMCInfo part;
			part.trackID  = mcev->trkID;
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
		std::sort(fG4mc.secParts.begin(), fG4mc.secParts.end(), PartMCInfo_sort());
		
		for (Int_t icls = 0; icls < event->NTrMCCluster(); icls++) {
			TrMCClusterR * cluster = event->pTrMCCluster(icls);
			if (cluster == 0) continue;
			if (cluster->GetMomentum() < 5.0e-2) continue;
			if (TrSET.find(cluster->GetGtrkID()) == TrSET.end()) continue;
			
			Int_t tkid = cluster->GetTkId();
			Int_t layJ = TkDBc::Head->GetJFromLayer(std::fabs(cluster->GetTkId()/100));
			Int_t trkid = cluster->GetGtrkID();
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
				for (Int_t ipart = 0; ipart < fG4mc.secParts.size(); ++ipart) {
					if (fG4mc.secParts.at(ipart).trackID != cluster->GetGtrkID()) continue;
					fG4mc.secParts.at(ipart).hits.push_back(hit);
					break;
				}
			}
		}

		std::sort(fG4mc.primPart.hits.begin(), fG4mc.primPart.hits.end(), HitTRKMCInfo_sort());
		for (Int_t ipart = 0; ipart < fG4mc.secParts.size(); ++ipart)
			std::sort(fG4mc.secParts.at(ipart).hits.begin(), fG4mc.secParts.at(ipart).hits.end(), HitTRKMCInfo_sort());

		if (fG4mc.secParts.size() >= 2) {
			VertexMCInfo vertex;
			vertex.status = true;
			vertex.coo[0] = fG4mc.secParts.at(0).coo[0];
			vertex.coo[1] = fG4mc.secParts.at(0).coo[1];
			vertex.coo[2] = fG4mc.secParts.at(0).coo[2];
			fG4mc.primVtx = vertex;
		}
	}

	timer.stop();
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

	if (evTree == 0) setTree("EventRti", "Run time information");
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
	if (event == 0)	return false;

	AMSSetupR::RTI rti;
	event->GetRTI(rti);

	Bool_t rebuildRTI = true;
	if (rti.utime == CurrUTime) rebuildRTI = false;
	else CurrUTime = rti.utime;
	
	fRti.uTime = event->UTime();
	
	MgntClock::TTime * ttime = MgntClock::ConvertFromUTimeToTTime(event->UTime(), MgntClock::ClockType::UTC);
	fRti.dateUTC = (ttime->tm_year + 1900) * 10000 + (ttime->tm_mon+1) * 100 + (ttime->tm_mday);
	fRti.timeUTC = (ttime->tm_hour) * 10000 + (ttime->tm_min) * 100 + (ttime->tm_sec);

	// ISS information
	AMSSetupR * setup = AMSSetupR::gethead();
	Double_t RPT[3] = {0};   // ISS coordinates (R, Phi, Theta) (GTOD)
	Double_t VelPT[3] = {0}; // ISS velocity (Vel rad/sec, VelPhi rad, VelTheta rad)
	Double_t YPR[3] = {0};   // ISS attitude (Yaw, Pitch, Roll)
	if (setup != 0) {
		Float_t rpt[3], velpt[3], yaw = 0, pitch = 0, roll = 0;
		setup->getISSTLE(rpt, velpt, fRti.uTime);
		setup->getISSAtt(roll, pitch, yaw, fRti.uTime);
		RPT[0] = rpt[0]; RPT[1] = rpt[1]; RPT[2] = rpt[2];
		VelPT[0] = velpt[0]; VelPT[1] = velpt[1]; VelPT[2] = velpt[2];
		YPR[0] = yaw; YPR[1] = pitch; YPR[2] = roll;
	}
	for (Int_t dir = 0; dir < 3; ++dir) {
		fRti.rptISS[dir] = RPT[dir];
		fRti.velISS[dir] = VelPT[dir];
		fRti.yprISS[dir] = YPR[dir];
	}

	// ISS solar array && backtracing (based on particle)
	Int_t isInShadow = -1;
	if (event->NParticle() > 0) {
		AMSPoint ic;
		Int_t idx = event->isInShadow(ic, 0);
		if (idx == 1) isInShadow = 1;
		else          isInShadow = 0;
	}
	fRti.isInShadow = isInShadow;

	// Backtracing (based on particle)
	Int_t isFromSpace = -1;
	Int_t backtrace[2][3] = { {-1, -1, -1}, {-1, -1, -1} };
	Bool_t isOpenBacktrace = false;
	while (setup != 0) {
		if (event->NParticle() == 0 || event->pParticle(0) == nullptr) break;
		ParticleR * part = event->pParticle(0);
		TrTrackR * trtk = part->pTrTrack();
		BetaHR * betaH = part->pBetaH();
		if (trtk == nullptr || betaH == nullptr) break;
		Int_t fitid = trtk->iTrTrackPar(1, 3, 21);
		if (fitid < 0) break;
	
		Double_t AMSTheta = part->Theta;
		Double_t AMSPhi   = part->Phi;
		Int_t Chrg[2]     = { 1, -1 };
		Double_t Mom      = std::fabs(trtk->GetRigidity(fitid));
		Double_t Beta     = std::fabs(betaH->GetBeta());

		HeaderR header;
		Double_t RPTO[3] = {0}; // particle final GTOD coo (R, Phi, Theta)
		Double_t GalLong = 0, GalLat = 0; // Galctic coo
		Double_t TraceTime = 0; // time of flight
		Double_t GPT[2] = {0}; // particle final GTOD direction (Phi, Theta)
		Int_t result = -1; // 0 unercutoff (i.e. atmospheric origin), 1 over cutoff (i.e. coming from space), 2 trapped, -1 error
		
		const Int_t nStable = 3;
		Double_t stableFT[nStable] = { 1.00, 1.15, 1.30 };
		Int_t    isOver[nStable] = { -1, -1, -1 }; // -1 error, 0 other, 1 weak, 2 strong
		if (isOpenBacktrace) {
			for (Int_t ift = 0; ift < nStable; ++ift) {
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
	
		Int_t successTrails = 0;
		Int_t countFromSpace = 0;
		for (Int_t ift = 0; ift < nStable; ++ift) {
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

	if (!rebuildRTI) return selectEvent(event);

	initEvent();
	timer.start();

	fRti.isInShadow      = isInShadow;
	fRti.isFromSpace     = isFromSpace;
	fRti.backtrace[0][0] = backtrace[0][0];
	fRti.backtrace[0][1] = backtrace[0][1];
	fRti.backtrace[0][2] = backtrace[0][2];
	fRti.backtrace[1][0] = backtrace[1][0];
	fRti.backtrace[1][1] = backtrace[1][1];
	fRti.backtrace[1][2] = backtrace[1][2];

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
	for (Int_t i = 0; i < 4; i++) {
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
	
	for (Int_t dir = 0; dir < 3; ++dir) {
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

	timer.stop();
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

	if (evTree == 0) setTree("EventTrg", "Anti-Coincidence Counter");
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
	if (event == 0)	return false;
	timer.start();

	Level1R * lvl1 = event->pLevel1(0);
	if (lvl1 == 0) return false;
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
	Bool_t extTrg        = ((fTrg.physicalPatt&0x80) > 0);
	Bool_t unBiasTrgTOF  = ((fTrg.physicalPatt&0x01) > 0);
	Bool_t unBiasTrgECAL = ((fTrg.physicalPatt&0x40) > 0);
	Bool_t physTrg       = ((fTrg.physicalPatt&0x3e) > 0);
	fTrg.bit = extTrg * 1 + unBiasTrgTOF * 2 + unBiasTrgECAL * 4 + physTrg * 8;

	timer.stop();
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

	if (evTree == 0) setTree("EventTof", "Time of Flight");
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
	if (event == 0)	return false;
	timer.start();

	fTof.numOfCluster = event->NTofCluster();
	fTof.numOfClusterH = event->NTofClusterH();
	fTof.numOfBeta = event->NBeta();
	fTof.numOfBetaH = event->NBetaH();

	Int_t tofLayerProj[4] = {0, 1, 1, 0}; // 0,1 := x,y
		
	// Tracker Information
	Float_t trackerZJ[9] = {0};
	if (checkEventMode(EventBase::ISS))
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerAJ(layJ);
	else
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerJ(layJ);

	while (recEv.iBeta >= 0) {
		BetaR * beta = event->pBeta(recEv.iBeta);
		if (beta == 0) break;
		fTof.statusBeta = true;
		fTof.beta = beta->Beta;
		Short_t betaPatt = beta->Pattern;
		break;
	}


	TofRecH::BuildOpt = 0; // normal
	const Short_t pattIdx[4] = { 1, 2, 4, 8 };
	while (recEv.iBetaH >= 0 && event->pBetaH(recEv.iBetaH) != 0) {
		BetaHR * betaH = event->pBetaH(recEv.iBetaH);
		if (betaH == 0) break;

		Int_t ncls[4] = {0};
		fTof.numOfInTimeCluster = event->GetNTofClustersInTime(betaH, ncls);

		fTof.statusBetaH = true;
		fTof.betaH = betaH->GetBeta();
		fTof.betaHBit = ((betaH->pTrTrack   ()) ? 1 : 0) +
                    ((betaH->pTrdTrack  ()) ? 2 : 0) +
                    ((betaH->pEcalShower()) ? 4 : 0);
		fTof.normChisqT  = betaH->GetNormChi2T();
		fTof.normChisqC  = betaH->GetNormChi2C();

		fTof.betaHPatt = 0;
		for (Int_t il = 0; il < 4; il++) {
			if (!betaH->TestExistHL(il)) continue;
			TofClusterHR * cls = betaH->GetClusterHL(il);
			if (cls == 0) continue;
			fTof.Q[il] = betaH->GetQL(il);
			fTof.betaHPatt += pattIdx[il];
			if (cls->IsGoodTime())
				fTof.betaHGoodTime += pattIdx[il];
		}

		Int_t nlay = 0;
		Float_t Qall_RMS = 0;
		fTof.Qall = betaH->GetQ(nlay, Qall_RMS);
		
		//Float_t Qupper, Qupper_RMS, Qlower, Qlower_RMS;
		//Qupper = betaH->GetQ(nlay, Qupper_RMS, 2, TofClusterHR::DefaultQOpt, 1100);
		//Qlower = betaH->GetQ(nlay, Qlower_RMS, 2, TofClusterHR::DefaultQOpt, 11);

		//TofChargeHR tofQH = betaH->gTofCharge();
		//Float_t Z[2], Zupper, Zlower, Z_prob[2], Zupper_prob, Zlower_prob;
	  //Z[0] = tofQH.GetZ(nlay, Z_prob[0], 0); // max prob charge (int)
		//Z[1] = tofQH.GetZ(nlay, Z_prob[1], 1); // next to max prob charge (int)
		//Zupper = tofQH.GetZ(nlay, Zupper_prob, 0, 1100);
		//Zlower = tofQH.GetZ(nlay, Zlower_prob, 0, 11);

		break;
	} // while loop - ibetaH > 0


	// Find Hits in the TOF supper layer (Time)
	std::vector<Int_t> betaHClsId(4, -1);
	if (fTof.statusBetaH) {
		BetaHR * betaH = event->pBetaH(recEv.iBetaH);
		for (Int_t it = 0; it < betaH->NTofClusterH(); ++it)
			betaHClsId.at(betaH->pTofClusterH(it)->Layer) = betaH->iTofClusterH(it);
	}

	Bool_t   isHasTime[2] = { false, false };
	Double_t avgTime[2][2] = { {0, 0}, {0, 0} };
	Double_t avgChrg[2] = {0, 0};
	for (Int_t it = 0; it < 2; ++it) {
		isHasTime[it] = (betaHClsId.at(2*it+0) >= 0 || betaHClsId.at(2*it+1) >= 0);
		if (isHasTime[it]) {
			TofClusterHR * ucls = (betaHClsId.at(2*it+0) >= 0) ? event->pTofClusterH(betaHClsId.at(2*it+0)) : nullptr;
			TofClusterHR * lcls = (betaHClsId.at(2*it+1) >= 0) ? event->pTofClusterH(betaHClsId.at(2*it+1)) : nullptr;
			Double_t chrg  = (((ucls) ? ucls->GetQSignal() : 0.) + ((lcls) ? lcls->GetQSignal() : 0.)) / ((ucls!=nullptr) + (lcls!=nullptr));
			Double_t wgval = ((ucls) ? (ucls->Time/ucls->ETime/ucls->ETime) : 0.) + ((lcls) ? (lcls->Time/lcls->ETime/lcls->ETime) : 0.);
			Double_t sumwg = ((ucls) ? (1./ucls->ETime/ucls->ETime) : 0.) + ((lcls) ? (1./lcls->ETime/lcls->ETime) : 0.);
			Double_t value = wgval / sumwg;
			Double_t sigma = 1./std::sqrt(sumwg);
			avgChrg[it]    = chrg;
			avgTime[it][0] = value;
			avgTime[it][1] = sigma;
		}
	}
	
	Int_t    nearHitId[2] = { -1, -1 };
	Double_t nearHitDt[2] = { 0, 0 };
	for (Int_t it = 0; it < event->NTofClusterH(); ++it) {
		TofClusterHR * cls = event->pTofClusterH(it);
		if (cls == nullptr) continue;
		if (betaHClsId.at(cls->Layer) < 0) continue;
		if (betaHClsId.at(cls->Layer) == it) continue;
		Int_t    slay = (cls->Layer / 2);
		Double_t dltT = std::fabs(cls->Time - avgTime[slay][0]) / cls->ETime;
		if (nearHitId[slay] < 0 || dltT < nearHitDt[slay]) {
			nearHitId[slay] = it;
			nearHitDt[slay] = dltT;
		} 
	}

	const Double_t TimeOneM = 3.335640e+00; // (speed of light)
	for (Int_t it = 0; it < 2; ++it) {
		if (nearHitId[it] < 0) continue;
		TofClusterHR * cls = event->pTofClusterH(nearHitId[it]);
		Int_t    lay  = cls->Layer;
		Double_t chrg = cls->GetQSignal();
		Double_t time = (cls->Time - avgTime[it][0]) / TimeOneM;
		
		fTof.statusExtCls[it] = true;
		fTof.extClsL[it] = lay;
		fTof.extClsQ[it] = chrg;
		fTof.extClsT[it] = time;
	}


	/*
	TofRecH::BuildOpt = 1; // TrTrack independent
	TofRecH::ReBuild(0);
	while (event->NBetaH() > 0) {
		BetaHR * betaHs = event->pBetaH(0);
		if (betaHs == 0) break;

		fTof.statusBetaHs = true;
		fTof.betaHs = betaHs->GetBeta();
		fTof.betaHBits = ((betaHs->pTrTrack   ()) ? 1 : 0) +
                     ((betaHs->pTrdTrack  ()) ? 2 : 0) +
                     ((betaHs->pEcalShower()) ? 4 : 0);
		fTof.normChisqTs  = betaHs->GetNormChi2T();
		fTof.normChisqCs  = betaHs->GetNormChi2C();

		fTof.betaHPatts = 0;
		for (Int_t il = 0; il < 4; il++) {
			if (!betaHs->TestExistHL(il)) continue;
			TofClusterHR * cls = betaHs->GetClusterHL(il);
			if (cls == 0) continue;
			fTof.Qs[il] = betaHs->GetQL(il);
			fTof.betaHPatts += pattIdx[il];
			if (cls->IsGoodTime())
				fTof.betaHGoodTimes += pattIdx[il];
		}

		Int_t nlays = 0;
		Float_t Qalls_RMS = 0;
		fTof.Qalls = betaHs->GetQ(nlays, Qalls_RMS);

		if ((fTof.betaHBits&2) == 2) {
			TrdTrackR * trd = betaHs->pTrdTrack();
			AMSDir trd_dir(trd->Theta, trd->Phi); trd_dir = trd_dir * -1;
			Float_t dir[3] = { Float_t(trd_dir[0]), Float_t(trd_dir[1]), Float_t(trd_dir[2]) };
			Float_t stt[6] = { trd->Coo[0], trd->Coo[1], trd->Coo[2], -dir[0], -dir[1], -dir[2] };
			std::copy(stt, stt+6, fTof.betaHStates);
		}

		break;
	}
	TofRecH::BuildOpt = 0; // normal
	*/

	timer.stop();
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

	if (evTree == 0) setTree("EventAcc", "Anti-Coincidence Counter");
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
	if (event == 0)	return false;
	timer.start();

	fAcc.numOfCluster = event->NAntiCluster();

	timer.stop();
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

	if (evTree == 0) setTree("EventTrk", "Silicon Tracker");
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
	if (event == 0)	return false;
	timer.start();

	// Beta Info
	Float_t Beta = 0;
	if      (recEv.iBetaH >= 0) Beta = event->pBetaH(recEv.iBetaH)->GetBeta();
	else if (recEv.iBeta  >= 0) Beta = event->pBeta(recEv.iBeta)->Beta;
	else                        Beta = 1;

	fTrk.numOfTrack = event->NTrTrack();

	Float_t expLJ[9][6] = {0};
	const Int_t qopt = TrClusterR::kAsym | TrClusterR::kGain | TrClusterR::kLoss | TrClusterR::kMIP;
	std::map<Int_t, Int_t> trackIDMap;
	while (recEv.iTrTrack >= 0) {
		for (Int_t itr = 0; itr <= event->NTrTrack(); ++itr) {
			Int_t jtr = itr - 1;
			if (jtr == recEv.iTrTrack) continue;
			TrTrackR * trtk = event->pTrTrack( ((jtr==-1)?recEv.iTrTrack:jtr) );
			if (trtk == 0) continue;
			
			if (jtr == -1) {
				Double_t BTdist = 0;
				fTrk.beamID = event->GetBeamPos(BTdist, trtk);
				fTrk.beamDist = BTdist;
			}
			
			TrackInfo track;

			const UShort_t _hasL1  =   1;
			const UShort_t _hasL2  =   2;
			const UShort_t _hasL34 =  12;
			const UShort_t _hasL56 =  48;
			const UShort_t _hasL78 = 192;
			const UShort_t _hasL9  = 256;
			track.bitPattJ = trtk->GetBitPatternJ();
			track.bitPattXYJ = trtk->GetBitPatternXYJ();
	
			Short_t isInner   = ((track.bitPattJ&_hasL34) > 0 &&
			                     (track.bitPattJ&_hasL56) > 0 &&
													 (track.bitPattJ&_hasL78) > 0) ? 1 : 0;
			Short_t isInnerXY = ((track.bitPattXYJ&_hasL34) > 0 &&
			                     (track.bitPattXYJ&_hasL56) > 0 &&
													 (track.bitPattXYJ&_hasL78) > 0) ? 2 : 0;
			Short_t isL2    = ((track.bitPattJ  &_hasL2) > 0) ?   4 : 0;
			Short_t isL2XY  = ((track.bitPattXYJ&_hasL2) > 0) ?   8 : 0;
			Short_t isL1    = ((track.bitPattJ  &_hasL1) > 0) ?  16 : 0;
			Short_t isL1XY  = ((track.bitPattXYJ&_hasL1) > 0) ?  32 : 0;
			Short_t isL9    = ((track.bitPattJ  &_hasL9) > 0) ?  64 : 0;
			Short_t isL9XY  = ((track.bitPattXYJ&_hasL9) > 0) ? 128 : 0;
			Short_t bitPatt = isInner + isInnerXY + isL2 + isL2XY + isL1 + isL1XY + isL9 + isL9XY;
	
			Short_t fitidInn = trtk->iTrTrackPar(1, 3, 21);
			if (fitidInn < 0) continue;
			track.bitPatt = bitPatt; 
			track.Qinner  = trtk->GetInnerQ_all(Beta, fitidInn).Mean; 

			const Short_t _nalgo = 2;
			const Short_t _algo[_nalgo] = { 1, 3 };
			const Short_t _npatt = 4;
			const Short_t _patt[_npatt] = { 3, 5, 6, 7 };
			for (Int_t algo = 0; algo < _nalgo; algo++) {
				for (Int_t patt = 0; patt < _npatt; patt++) {
					Int_t fitid = trtk->iTrTrackPar(_algo[algo], _patt[patt], 21);
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
				
					for (Int_t il = 0; il < 9; ++il) {
						AMSPoint pntLJ;
						AMSDir   dirLJ;
						trtk->InterpolateLayerJ(il+1, pntLJ, dirLJ, fitid);
						track.stateLJ[algo][patt][il][0] = pntLJ[0];
						track.stateLJ[algo][patt][il][1] = pntLJ[1];
						track.stateLJ[algo][patt][il][2] = pntLJ[2];
						track.stateLJ[algo][patt][il][3] = -dirLJ[0];
						track.stateLJ[algo][patt][il][4] = -dirLJ[1];
						track.stateLJ[algo][patt][il][5] = -dirLJ[2];
					}
				} // for loop - pattern
			}

			for (Int_t layJ = 1; layJ <= 9; layJ++) {
				AMSPoint pnt;
				AMSDir   dir;
				trtk->InterpolateLayerJ(layJ, pnt, dir, fitidInn);
				expLJ[layJ-1][0] = pnt[0];
				expLJ[layJ-1][1] = pnt[1];
				expLJ[layJ-1][2] = pnt[2];
				expLJ[layJ-1][3] = -dir[0];
				expLJ[layJ-1][4] = -dir[1];
				expLJ[layJ-1][5] = -dir[2];
			}

			Bool_t hasHitLJ[9] = {false};
			Int_t hitClsIDx[9]; std::fill_n(hitClsIDx, 9, -1);
			Int_t hitClsIDy[9]; std::fill_n(hitClsIDy, 9, -1);
			HitTRKInfo * trackHit[9] = {0};
			Bool_t isDeadChanX[9] = {false};
			for (UInt_t ilay = 0; ilay < 9; ++ilay) {
				if (!trtk->TestHitLayerJ(ilay+1)) continue;
				TrRecHitR * recHit = trtk->GetHitLJ(ilay+1);
				if (recHit == 0) continue;

				Int_t tkid = recHit->GetTkId();
				
				TrClusterR * xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : 0;
				TrClusterR * ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : 0;

				Short_t clsIdX = (xcls) ? recHit->GetXClusterIndex() : -1;
				Short_t clsIdY = (ycls) ? recHit->GetYClusterIndex() : -1;

				Int_t   side = (xcls ? 1 : 0) + (ycls ? 2 : 0);
				Int_t   mult = recHit->GetResolvedMultiplicity(); // -1  resolved multiplicty coordinates
				                                                  // > 0 requested multiplicty coordinates
				Float_t xloc = (xcls) ? (xcls->GetXCofG() + xcls->GetSeedAddress()) : (recHit->GetDummyX() + 640);
				Float_t yloc = (ycls->GetXCofG() + ycls->GetSeedAddress());
  			
				Float_t xchrg = (xcls) ? TMath::Sqrt(recHit->GetSignalCombination(0, qopt, 1)) : -1.;
  			Float_t ychrg = TMath::Sqrt(recHit->GetSignalCombination(1, qopt, 1));
				
				AMSPoint coo = recHit->GetCoord(mult, 3); // (CIEMAT+PG)/2

    		TkSens tksens(coo, EventBase::checkEventMode(EventBase::MC));
				Int_t sens = (tksens.GetLadTkID() != 0) ? tksens.GetSensor() : -1;

				hasHitLJ[ilay] = true;
				hitClsIDx[ilay] = (xcls) ? recHit->GetXClusterIndex() : -1;
				hitClsIDy[ilay] = (ycls) ? recHit->GetYClusterIndex() : -1;
				isDeadChanX[ilay] = (tksens.GetLadTkID() != 0 && tksens.ReadChanX == -1); 

				HitTRKInfo hit;
				hit.clsId[0] = clsIdX;
				hit.clsId[1] = clsIdY;
				hit.layJ     = ilay+1;
				hit.tkid     = tkid;
				hit.sens     = sens;
				hit.mult     = mult; 
				hit.side     = side;
				hit.cofg[0]  = xloc;
				hit.cofg[1]  = yloc;
				hit.chrg[0]  = xchrg;
				hit.chrg[1]  = ychrg;
				hit.coo[0]   = coo[0];
				hit.coo[1]   = coo[1];
				hit.coo[2]   = coo[2];
				
				trackHit[ilay] = &hit;

				track.hits.push_back(hit);
			} // for loop - layer
			if (track.hits.size() > 1) std::sort(track.hits.begin(), track.hits.end(), HitTRKInfo_sort());
		
			Int_t offTrID = ((jtr==-1)?recEv.iTrTrack:jtr);
			Int_t selTrID = fTrk.tracks.size();
			trackIDMap[offTrID] = selTrID;
			
			fTrk.tracks.push_back(track);
		}

		break;
	}

	// Vertex
	TrReconQ reconQ(event);
	reconQ.BuildVertex();
	AMSPoint recQ_CVtx, recQ_DVtx;
	Int_t recQ_TrID[2] = { -1, -1 };
	if (reconQ.GetVertex(recQ_CVtx, recQ_DVtx, recQ_TrID) != 0) {
		std::map<Int_t, Int_t>::iterator ita = trackIDMap.find(recQ_TrID[0]);
		std::map<Int_t, Int_t>::iterator itb = trackIDMap.find(recQ_TrID[1]);
		if (ita != trackIDMap.end() && itb != trackIDMap.end()) {
			std::pair<Int_t, Int_t> trId = std::minmax(ita->second, itb->second);
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

	Int_t   mapL[9] = {0, 1, 2, 2, 2, 2, 2, 2, 3};
	Float_t qMax[4] = { -1, -1, -1, -1 };
	Bool_t  onlyY[4] = { false, false, false, false };
	TrRecHitR * mapH[8] = {0};

	for (Int_t il = 0; il < 9; il++) {
		if (MgntTrHit::MapHits.find(il+1) == MgntTrHit::MapHits.end()) continue;
		std::vector<TrHit *> & itHits = MgntTrHit::MapHits[il+1];
		for (std::vector<TrHit *>::iterator ht = itHits.begin(); ht != itHits.end(); ++ht) {
			TrRecHitR * recHit = (*ht)->Hit();
			Float_t chrg = (*ht)->Chrg();
 			if (chrg <= qMax[mapL[il]]) continue;
			qMax[mapL[il]] = chrg;
			mapH[mapL[il]] = recHit;
			onlyY[mapL[il]] = recHit->OnlyY();
		} // for loop --- hit
	}	// for loop --- layerJ

	for (Int_t it = 0; it < 4; it++) {
		Float_t chrg = qMax[it];
		if (chrg < 0) continue;

		TrRecHitR * recHit = mapH[it];
		Short_t layJ = recHit->GetLayerJ();

		Int_t mult = -1;
		AMSPoint coo;
		Float_t distLJ = 500;
		for (Int_t im = 0; im < recHit->GetMultiplicity(); im++) {
			AMSPoint refcoo = recHit->GetCoord(im, 3);
			Float_t dx = (refcoo[0] - expLJ[layJ-1][0]);
			Float_t dy = (refcoo[1] - expLJ[layJ-1][1]);
			Float_t dist = (onlyY[it]) ? std::fabs(dy) : std::sqrt(dx * dx + dy * dy);
			if (dist > distLJ) continue;
			mult = im;
			distLJ = dist;
			coo = refcoo;
		}
			
		TrClusterR * xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : 0;
		TrClusterR * ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : 0;
		Short_t clsIdX = (xcls) ? recHit->GetXClusterIndex() : -1;
		Short_t clsIdY = (ycls) ? recHit->GetYClusterIndex() : -1;

		Int_t   side = (xcls ? 1 : 0) + (ycls ? 2 : 0);
		Float_t xloc = (xcls) ? (xcls->GetXCofG() + xcls->GetSeedAddress()) : (recHit->GetDummyX() + 640);
		Float_t yloc = (ycls->GetXCofG() + ycls->GetSeedAddress());
				
		Float_t xchrg = (xcls) ? TMath::Sqrt(recHit->GetSignalCombination(0, qopt, 1)) : -1.;
  	Float_t ychrg = TMath::Sqrt(recHit->GetSignalCombination(1, qopt, 1));
					
		TkSens tksens(coo, EventBase::checkEventMode(EventBase::MC));
		Int_t sens = (tksens.GetLadTkID() != 0 && xcls) ? tksens.GetSensor() : -1;

		HitTRKInfo hit;
		hit.clsId[0] = clsIdX;
		hit.clsId[1] = clsIdY;
		hit.layJ     = layJ;
		hit.tkid     = recHit->GetTkId();
		hit.sens     = sens;
		hit.mult     = mult; 
		hit.side     = side;
		hit.cofg[0]  = xloc;
		hit.cofg[1]  = yloc;
		hit.chrg[0]  = xchrg;
		hit.chrg[1]  = ychrg;
		hit.coo[0]   = ((hit.side&1)==1) ? coo[0] : 0;
		hit.coo[1]   = coo[1];
		hit.coo[2]   = coo[2];

		fTrk.maxQHit[it] = hit; 
	}

	timer.stop();
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

	if (evTree == 0) setTree("EventTrd", "Transition Radiation");
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
	if (event == 0)	return false;
	timer.start();

	fTrd.numOfTrack = event->NTrdTrack();
	fTrd.numOfHTrack = event->NTrdHTrack();
		
	// numOfHSegVtx (by HY.Chou)
	// number of TRDH segments that make a vertex with TrTrack
	// between Z = 80 cm and 200 cm
	Int_t nseg = event->nTrdHSegment();
	TrTrackR * trtk = (recEv.iTrTrack >= 0) ? event->pTrTrack(recEv.iTrTrack) : nullptr;
	Int_t      trId = (trtk) ? trtk->iTrTrackPar(1, 3, 21) : -1;
	if (trtk && trId >= 0) {
		const Double_t trCooZ = 120.;
		AMSPoint trPnt; AMSDir trDir;
		trtk->Interpolate(trCooZ, trPnt, trDir);
		Double_t trX = trPnt[0];
		Double_t trY = trPnt[1];
		Double_t trTanX = trDir[0] / trDir[2];
		Double_t trTanY = trDir[1] / trDir[2];
		Double_t trX0   = trX - trTanX * trCooZ;
		Double_t trY0   = trY - trTanY * trCooZ;
		Double_t trNRig = std::fabs(trtk->GetRigidity(trId));
		if (trNRig < 0.8) trNRig = 0.8;

		const Double_t MsCooSgmFact = 2.843291e-01;
		const Double_t MsDirSgmFact = 7.108227e-03;
		Double_t msSgmC = MsCooSgmFact / trNRig;
		Double_t msSgmM = MsDirSgmFact / trNRig;
		Double_t msSgmR = std::sqrt(msSgmC*msSgmC + msSgmM*msSgmM*trCooZ*trCooZ);

		Short_t numOfVtx[2] = { 0, 0 };
		std::vector<std::pair<Double_t, Double_t> > vtxCooXZ;
		std::vector<std::pair<Double_t, Double_t> > vtxCooYZ;

		const Double_t SgmLimit = 7.;
		const Double_t ZLimit[2] = { 200., 60. };
		for (Int_t iseg = 0; iseg < event->NTrdHSegment(); ++iseg) {
			TrdHSegmentR * seg = event->pTrdHSegment(iseg);
			if (seg == nullptr) continue;
			Double_t segR0 = (seg->r - seg->m * seg->z);
			Double_t trR   = ((seg->d==0) ? (trX + trTanX * (seg->z - trCooZ)) : (trY + trTanY * (seg->z - trCooZ)));
			Double_t trM   = ((seg->d==0) ? trTanX : trTanY);
			Double_t dr    = std::fabs((trR - seg->r) / std::sqrt(seg->er*seg->er + msSgmC*msSgmC));
			Double_t dm    = std::fabs((trM - seg->m) / std::sqrt(seg->em*seg->em + msSgmM*msSgmM));
			Bool_t isTrSeg = (dr < SgmLimit && dm < SgmLimit);
			if (isTrSeg) continue;

			Double_t sgmSegR = std::sqrt(seg->er*seg->er+seg->em*seg->em*seg->z*seg->z);
			Double_t sgmSegM = seg->em;
			Double_t sgmVtxR = std::sqrt(msSgmR*msSgmR+sgmSegR*sgmSegR);
			Double_t sgmVtxM = std::sqrt(msSgmM*msSgmM+sgmSegM*sgmSegM);

			Double_t vtxDR = (segR0  - ((seg->d==0) ? trX0   : trY0));
			Double_t vtxDM = (seg->m - ((seg->d==0) ? trTanX : trTanY));

			Double_t vtxZ    = -(vtxDR / vtxDM);
			Double_t sgmVtxZ = std::fabs(vtxZ) * std::sqrt((sgmVtxR*sgmVtxR)/(vtxDR*vtxDR) + (sgmVtxM*sgmVtxM)/(vtxDM*vtxDM+1e-4));
			if (sgmVtxZ > 10.) sgmVtxZ = 10.;

			Double_t lmtZu = (ZLimit[0] + sgmVtxZ * SgmLimit);
			Double_t lmtZl = (ZLimit[1] - sgmVtxZ * SgmLimit);

			if (vtxZ > lmtZu || vtxZ < lmtZl) continue;

			numOfVtx[seg->d]++;
			if (seg->d==0) vtxCooXZ.push_back(std::make_pair(vtxZ, sgmVtxZ));
			else           vtxCooYZ.push_back(std::make_pair(vtxZ, sgmVtxZ));
		}

		fTrd.numOfHSegVtx[0] = numOfVtx[0];
		fTrd.numOfHSegVtx[1] = numOfVtx[1];
	}


	// TrdKCluster
	Float_t TOF_Beta = 1;
	if      (recEv.iBetaH >= 0) TOF_Beta = std::fabs(event->pBetaH(recEv.iBetaH)->GetBeta());
	else if (recEv.iBeta  >= 0) TOF_Beta = std::fabs(event->pBeta(recEv.iBeta)->Beta);
	
	TrdKCluster * trdkcls = TrdKCluster::gethead();
	for (Int_t kindOfFit = 0; kindOfFit <= 1; ++kindOfFit) {
		if (!(checkEventMode(EventBase::BT) ||
					(trdkcls->IsReadAlignmentOK == 2 && trdkcls->IsReadCalibOK == 1))) break;
		Bool_t isOK = false;
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
				  Int_t fitid_max = trtk->iTrTrackPar(1, 0, 21);
				  if (fitid_max < 0) break;
				  trdkcls->SetTrTrack(trtk, fitid_max);
				  isOK = true;
        }
				break;
			default :
				MgntSys::Error(LocAddr(), MgntSys::MESSAGE("TrdKCluster failed Building! Exiting ..."));
				MgntSys::Exit(EXIT_FAILURE);
		}
		if (!isOK) continue;

		// It speeds lots of time 0.02 sec
		Float_t Q = -1;
		Float_t Qerror = -1;
		Int_t QnumberOfHit = -1;
		Int_t Qstatus = trdkcls->CalculateTRDCharge(0, TOF_Beta);
		if (Qstatus >= 0) {
			Q = trdkcls->GetTRDCharge();
			Qerror = trdkcls->GetTRDChargeError();
			QnumberOfHit = trdkcls->GetQNHit();
		}
		else continue;

		Int_t nhits = 0; //To be filled with number of hits taken into account in Likelihood Calculation
		Float_t threshold = 15; //ADC above which will be taken into account in Likelihood Calculation,  15 ADC is the recommended value for the moment.
		Double_t llr[3] = {-1, -1, -1}; //To be filled with 3 LikelihoodRatio :  e/P, e/H, P/H
		if      (kindOfFit == 0) trdkcls->GetLikelihoodRatio_TRDRefit(threshold, llr, nhits);
		else if (kindOfFit == 1) trdkcls->GetLikelihoodRatio_TrTrack(threshold, llr, nhits);
		if (llr[0] < 0 || llr[1] < 0 || llr[2] < 0) continue;

		Int_t numberOfHit = 0;
		for (Int_t ih = 0; ih < nhits; ih++) {
			TrdKHit * hit = trdkcls->GetHit(ih);
			if (hit == nullptr) continue;
			if (checkEventMode(EventBase::ISS) || checkEventMode(EventBase::BT))
				if (!hit->IsCalibrated) continue;
			if (checkEventMode(EventBase::ISS))
				if (!hit->IsAligned) continue;
			//Int_t lay = hit->TRDHit_Layer;
			//Float_t amp = hit->TRDHit_Amp;
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

	// Tracker Information
	Float_t trackerZJ[9] = {0};
	if (checkEventMode(EventBase::ISS))
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerAJ(layJ);
	else
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerJ(layJ);

	Short_t trdIdx = -1;
	if      (recEv.iTrdHTrack >= 0) trdIdx = 1;
	else if (recEv.iTrdTrack  >= 0) trdIdx = 0;

	while (trdIdx >= 0) {
		Bool_t isOK = false;
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

	timer.stop();
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

	if (evTree == 0) setTree("EventRich", "Ring Imaging Cherenkov");
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
	if (event == 0)	return false;
	timer.start();
		
	// Rich cuts from Javie/Jorge
  const Float_t richRadZ[2] = {-73.65, -74.65};          // aerogel / NaF
	const Float_t cut_prob = 0.01;                         // Kolmogorov test prob
	const Float_t cut_pmt = 3;                             // number of PMTs
	const Float_t cut_collPhe[2] = {0.6, 0.4};             // ring phe / total phe (NaF, Aerogel)
	const Float_t cut_chargeConsistency = 5;               // hit / PMT charge consistency test
	const Float_t cut_betaConsistency[2] = {0.01, 0.005};  // beta_lip vs beta_ciemat consistency (NaF, Aerogel)
	const Float_t cut_expPhe[2] = {1, 2};                  // expected number of phe (NaF, Aerogel)
	const Float_t cut_aerogelExternalBorder = 3350;        // aerogel external border (r**2)
	// modify 3500 -> 3360 (by S.H.)
	// modify 3500 -> 3350 (by HY.Chou)
	const Float_t cut_PMTExternalBorder = 4050;             // PMT external border (r**2) (by HY.Chou)
	const Float_t cut_aerogelNafBorder[2] = {17.25, 17.75}; // aerogel/NaF border (NaF, Aerogel)
	// modify (17, 18) -> (17.25, 17.75) (by HY.Chou)
	const Float_t cut_distToTileBorder[2] = {0.08, 0.05};   // distance to tile border 
	
	const Int_t nBadTile = 5;
	const Int_t kBadTile_Offical[nBadTile]  = { 3, 7, 87, 100, 108 }; // tiles with bad beta reconstruction
	const Int_t kBadTile_MgntTile[nBadTile] = { 13, 23, 58, 86, 91 }; // tiles with bad beta reconstruction
	
	//fRich.numOfRing = event->NRichRing();
	//fRich.numOfHit = event->NRichHit();

	// RichVeto - start
	while (recEv.iTrTrack >= 0) {
		TrTrackR * trtk = event->pTrTrack(recEv.iTrTrack);
		if (trtk == 0) break;
		Int_t fitid = trtk->iTrTrackPar(1, 3, 21);
		if (fitid < 0) break;

		AMSPoint ems_coo;
		AMSDir   ems_dir;
		Int_t    tmp_kindOfRad = 0;
		for (Int_t it = 0; it < 4; ++it) {
			trtk->Interpolate(richRadZ[tmp_kindOfRad], ems_coo, ems_dir, fitid);
			RichOffline::TrTrack tmp_track(ems_coo, ems_dir);
			RichOffline::RichRadiatorTileManager tmp_mgntTile(&tmp_track);
			tmp_kindOfRad = tmp_mgntTile.getkind() - 1;
			if (tmp_kindOfRad < 0 || tmp_kindOfRad > 1) break;
		}
		if (tmp_kindOfRad < 0 || tmp_kindOfRad > 1) break;
		
		Float_t ems[6] = {0};
		ems[0] =  ems_coo[0]; ems[1] =  ems_coo[1]; ems[2] =  ems_coo[2];
		ems[3] = -ems_dir[0]; ems[4] = -ems_dir[1]; ems[5] = -ems_dir[2];

		// Radiator tile manager
		RichOffline::TrTrack track(ems_coo, ems_dir);
		RichOffline::RichRadiatorTileManager mgntTile(&track);
	
		if (mgntTile.getkind() < 1 || mgntTile.getkind() > 2) break;

		Short_t kindOfRad      = mgntTile.getkind() - 1;
		Short_t tileOfRad      = mgntTile.getcurrenttile();
		Float_t rfrIndex     = mgntTile.getindex();
		AMSPoint richems     = mgntTile.getemissionpoint();
		Float_t distToBorder = mgntTile.getdistance();
			
		Int_t countKB = 0;
		for (Int_t kb = 0; kb < nBadTile; kb++) {
			if (tileOfRad == kBadTile_MgntTile[kb]) break;
			countKB++;
		}
		Bool_t isGoodTile = (countKB == nBadTile);

		// geometry cut
		Bool_t isStruct = false;
		if ((kindOfRad == 1 && (std::max(std::fabs(richems[0]), std::fabs(richems[1])) > cut_aerogelNafBorder[0])) ||
				(kindOfRad == 0 && (std::max(std::fabs(richems[0]), std::fabs(richems[1])) < cut_aerogelNafBorder[1] ||
					(richems[0]*richems[0]+richems[1]*richems[1]) > cut_aerogelExternalBorder))) isStruct = true;
		//if (distToBorder < cut_distToTileBorder[kindOfRad] * std::fabs(1./ems[5])) isStruct = true;
		Bool_t isInFiducialVolume = !isStruct;
	
		// Number of photoelectrons expected for a given track, beta and charge.
		const Int_t npart = 5;
		const Bool_t openCal[npart] = { 1, 0, 0, 1, 0 };
		const Double_t chrg[npart] = { 1., 1., 1., 1., 1. };
		const Double_t mass[npart] = { 0.000510999,   // electron
		                               0.1395701835,  // pion
																	 0.493667,      // kaon
																	 0.938272297,   // proton
																	 1.876123915 }; // deuterium
		Float_t numOfExpPE[npart] = {0};
		std::fill_n(numOfExpPE, npart, -1);
		Double_t rigAbs = std::fabs(trtk->GetRigidity(fitid));
		for (UInt_t i = 0; i < npart; ++i) {
			if (!openCal[i]) continue;
			Double_t massChrg = mass[i] / chrg[i];
			Double_t beta = 1. / std::sqrt((massChrg * massChrg / rigAbs / rigAbs) + 1); 
			Double_t exppe = RichRingR::ComputeNpExp(trtk, beta, chrg[i]);
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
	Bool_t isSuccRing = false;
	while (recEv.iRichRing >= 0 && fRich.kindOfRad >= 0) {
		// RichRingR
		RichRingR * rich = event->pRichRing(recEv.iRichRing);
		Int_t kindOfRad = rich->IsNaF() ? 1 : 0;
		Int_t tileOfRad = rich->getTileIndex();
		//Float_t diffDist = std::fabs(rich->DistanceTileBorder() - fRich.distToBorder);
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
	Short_t numOfPrimHit[2] = {0, 0};
	Float_t numOfPrimPE[2]  = {0, 0}; 
	Short_t numOfOthHit[2] = {0, 0};
	Float_t numOfOthPE[2]  = {0, 0};
	for (Int_t it = 0; it < event->NRichHit(); ++it) {
		RichHitR * hit = event->pRichHit(it);
		if (hit == nullptr) continue;
		Bool_t  used[2] = { false, false };
		for (Int_t iring = 0; iring < event->NRichRing(); ++iring) {
			Bool_t isUsed = hit->UsedInRingNumber(iring);
			if (!isUsed) continue;
			if (iring == recEv.iRichRing && isSuccRing) used[0] = true;
			else                                        used[1] = true;
		}
		Bool_t isOthers = (!used[0]);
		
		Short_t cross = (hit->IsCrossed() ? 0 : 1);
		Float_t npe   = hit->Npe;

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

	timer.stop();
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

	if (evTree == 0) setTree("EventEcal", "Electromagnetic Calorimeter");
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
	if (event == 0)	return false;
	timer.start();

	// Tracker Information
	Float_t trackerZJ[9] = {0};
	if (checkEventMode(EventBase::ISS))
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerAJ(layJ);
	else
		for (Int_t layJ = 1; layJ <= 9; layJ++)
			trackerZJ[layJ-1] = TkDBc::Head->GetZlayerJ(layJ);
	
	// threshold : remove MIPs and low energy particles (lepton study 0.250)
	// threshold : remove very low energy particles (proton study 0.050)
	Float_t threshold = 0.050;

	//fEcal.numOfShower = event->NEcalShower();
	
	while (recEv.iEcalShower >= 0) {
		
		for (Int_t ish = 0; ish < event->NEcalShower()+1; ++ish) {
			Int_t jsh = ish - 1;
			if (jsh == recEv.iEcalShower) continue;
			EcalShowerR * ecal = event->pEcalShower( ((jsh==-1)?recEv.iEcalShower:jsh) );
			if (ecal == 0) continue;
		
			ShowerInfo shower;
		
			shower.energyD = 1.e-3 * ecal->EnergyD;
			Float_t energyA = ecal->EnergyA;
			shower.energyE = ecal->EnergyE;
			Float_t energyC = ecal->EnergyC;
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
			Int_t zeh = 1;
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
					Int_t nlay = 0; Float_t prob = 0;
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
	for (UInt_t ih = 0; ih < event->NEcalHit(); ++ih) {
		EcalHitR * ecalHit = event->pEcalHit(ih);
		if (ecalHit == 0) continue;

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
	
	timer.stop();
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
  DataSelection::gScaleFunc.SetParameter(0, DataSelection::gScaleFact);
}

int DataSelection::processEvent(AMSEventR * event, AMSChain * chain) {
#if Debug == true
	//std::cerr << "Debug : Now, DataSelection::processEvent()\n";
#endif

	if (event == 0) return -1;

	bool statusList = true;
	//bool statusRti = true;
	bool statusTrg = true;
	bool statusTof = true;
	bool statusAcc = true;
	bool statusTrk = true;
	bool statusTrd = true;
	bool statusRich = true;
	bool statusEcal = true;

	MgntClock::HrsTimer timer;
	timer.start();

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

	timer.stop();

#if Debug == true
	const Float_t limitTimer = 5.0;
	Float_t totTimer = timer.time() + recEv.time() + rti.time();
	if (totTimer > limitTimer) {
		std::cout << Form("\nRUN %u  EVENT %u\n", event->Run(), event->Event());
		std::cout << Form("REAL TIME : %14.8f (SEC)   100.00%\n", totTimer);
		std::cout << Form("    RECON   %14.8f (SEC)   %6.2f%\n",
		                  recEv.time(), recEv.time() / totTimer * 100);
		std::cout << Form("     LIST   %14.8f (SEC)   %6.2f%\n",
		                   list.time(),  list.time() / totTimer * 100);
		std::cout << Form("      RTI   %14.8f (SEC)   %6.2f%\n",
		                    rti.time(),   rti.time() / totTimer * 100);
		std::cout << Form("      TRG   %14.8f (SEC)   %6.2f%\n",
		                    trg.time(),   trg.time() / totTimer * 100);
		std::cout << Form("      TOF   %14.8f (SEC)   %6.2f%\n",
		                    tof.time(),   tof.time() / totTimer * 100);
		std::cout << Form("      ACC   %14.8f (SEC)   %6.2f%\n",
		                    acc.time(),   acc.time() / totTimer * 100);
		std::cout << Form("      TRK   %14.8f (SEC)   %6.2f%\n",
		                    trk.time(),   trk.time() / totTimer * 100);
		std::cout << Form("      TRD   %14.8f (SEC)   %6.2f%\n",
		                    trd.time(),   trd.time() / totTimer * 100);
		std::cout << Form("     RICH   %14.8f (SEC)   %6.2f%\n",
		                   rich.time(),  rich.time() / totTimer * 100);
		std::cout << Form("     ECAL   %14.8f (SEC)   %6.2f%\n",
		                   ecal.time(),  ecal.time() / totTimer * 100);
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

int DataSelection::preselectEvent(AMSEventR * event) {
	if (event == 0)	return -1;
	EventList::gWeight = 1.;

	// Resolution tuning
	if (EventBase::checkEventMode(EventBase::MC)) {
		event->SetDefaultMCTuningParameters();
		TrExtAlignDB::SmearExtAlign();
		TRCLFFKEY.UseSensorAlign = 0;
	}

	//-----------------------------//
	//----  Fast Preselection  ----//
	//-----------------------------//
  //std::cout << Form("============= <Run %u Event %u> =============\n", event->Run(), event->Event());

	// ~1~ (Based on BetaH(Beta))
	Bool_t isDownBeta = false;
	for (Int_t ibta = 0; ibta < event->NBeta(); ++ibta)
		if (event->pBeta(ibta)->Beta > 1.0e-2) { isDownBeta = true; break; }

	TofRecH::BuildOpt = 0; // normal
	Bool_t isDownBetaH = false;
	for (Int_t ibta = 0; ibta < event->NBetaH(); ++ibta)
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
	const UShort_t _hasTrL34 =  12;
	const UShort_t _hasTrL56 =  48;
	const UShort_t _hasTrL78 = 192;
	UShort_t trBitPattJ = trtkSIG->GetBitPatternJ();
	Bool_t   isTrInner  = ((trBitPattJ&_hasTrL34) > 0 && 
	                       (trBitPattJ&_hasTrL56) > 0 && 
												 (trBitPattJ&_hasTrL78) > 0);
	if (!isTrInner) return -6001;

	// ~7~ (Only for Antiproton to Proton Flux Ratio Study)
	Bool_t isAppStudy = true;
	if (EventBase::checkEventMode(EventBase::ISS) && isAppStudy) {
		Bool_t hasTr    = false;
		Bool_t isScale  = true;
		Double_t sclRig = 0.0;
		const Int_t npatt = 4;
		const Int_t spatt[npatt] = { 3, 5, 6, 7 };
		for (Int_t ip = 0; ip < npatt; ++ip) {
			Int_t fitid = trtkSIG->iTrTrackPar(1, spatt[ip], 21);
			if (fitid < 0) continue;
			hasTr = true;
			sclRig = trtkSIG->GetRigidity(fitid);
			if (sclRig < 0.0) { isScale = false; break; }
		}
		if (hasTr && isScale) {
			Double_t scaleProb = gScaleFunc.Eval(sclRig);
    	if (MgntNum::Compare(gRandom.Uniform(0, 1), scaleProb) > 0) return -7001;
    	else EventList::gWeight *= (1. / scaleProb);
		}
	}
	
	// ~8~ (Based on RTI)
	if (EventBase::checkEventMode(EventBase::ISS) && checkOption(DataSelection::RTI)) {
		if (!rti.processEvent(event)) return -8001;

		Double_t minStormer = *std::min_element(rti.fRti.cutoffStormer, rti.fRti.cutoffStormer+4);
		Double_t minIGRF    = *std::min_element(rti.fRti.cutoffIGRF, rti.fRti.cutoffIGRF+4);
		Double_t minCf      =  std::min(minStormer, minIGRF);
		Double_t maxRig     = 0.0;

		Bool_t hasTr = false;
		const Int_t npatt = 4;
		const Int_t spatt[npatt] = { 3, 5, 6, 7 };
		for (Int_t ip = 0; ip < npatt; ++ip) {
			Int_t fitid = trtkSIG->iTrTrackPar(1, spatt[ip], 21);
			if (fitid < 0) continue;
			hasTr = true;
			maxRig = std::max(maxRig, std::fabs(trtkSIG->GetRigidity(fitid)));
		}

		const Double_t minFact = 1.2;
		if ( hasTr && (maxRig < (minFact * minCf)) ) return -8002;
	}

	//--------------------------//
	//----  Reconstruction  ----//
	//--------------------------//
	if (!recEv.rebuild(event)) return -9999;

	return 0;
}

int DataSelection::selectEvent(AMSEventR * event) {
	if (event == 0) return -1;

	// User Define

	return 0;
}

int DataSelection::analysisEvent(AMSEventR * event) {
	if (event == 0) return -1;

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
	if (chain == 0) return false;
	TFile * file = chain->GetFile();
	std::string filePath = file->GetName();
	if (event == 0) return false;
	UInt_t runID   = event->Run();
	UInt_t eventID = event->Event();
		
	std::map<UInt_t, RunTagInfo>::iterator it = fRunTag.find(runID);
	if (it == fRunTag.end()) {
		RunTagInfo info;
		info.runID = runID;
		info.eventFT = eventID;
		info.eventLT = eventID;
		info.numOfSelEvent = 1;
		info.file.push_back(filePath);
		
		if (EventBase::checkEventMode(EventBase::ISS)) {
			MgntClock::TTime * ttime = MgntClock::ConvertFromUTimeToTTime(runID, MgntClock::ClockType::UTC);
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
	if (file == 0) return;
	file->cd();
	TTree * tree = new TTree("runTag", "RunTag information");
	RunTagInfo info;
	tree->Branch("runTag", &info);
	for (std::map<UInt_t, RunTagInfo>::iterator it = fRunTag.begin(); it != fRunTag.end(); ++it) {
		info = it->second;
		UInt_t nTrgEv = UInt_t(
		                (info.numOfSelEvent == 1) ? 
										2 * info.eventFT : 
										((Double_t(info.numOfSelEvent) / Double_t(info.numOfSelEvent - 1)) * 
										  Double_t(info.eventLT - info.eventFT)) 
										);
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
	fTimer.start();
}

YiNtuple::~YiNtuple() {
#if Debug == true
	std::cerr << "Debug : Now, YiNtuple::~YiNtuple()\n";
#endif
	init();

	fTimer.stop();
	fTimer.print();
}

inline void YiNtuple::init() {
#if Debug == true
	std::cerr << "Debug : Now, YiNtuple::init()\n";
#endif

	fGroup.first = 0;
	fGroup.second = -1;
	fFileList.clear();
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

void YiNtuple::readDataFrom(const std::string& file_list, Long64_t group_th, Long64_t group_size) {
#if Debug == true
	std::cerr << "Debug : Now, YiNtuple::readDataFrom()\n";
#endif

	std::cout << "\n**--------------------------------------------**\n";
	std::cout << "\n**    Read Data Form Source File List Info    **\n";
	std::cout << "\n**--------------------------------------------**\n";

	// start check sourceFileList.txt
	std::vector<std::string>&& flist = MgntIO::ReadFileContent(file_list);
	if (flist.size() == 0) {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("ROOT file list cannot be opend! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}
	// end check sourceFileList.txt

	// start load data with group
	if (group_th == 0 && group_size == -1) group_size = flist.size();
	if (group_size <= 0 || group_size > flist.size() || group_th < 0 || group_th >= flist.size()) {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("Group format has error(1)! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}
	Long64_t begin = group_th * group_size;
	Long64_t end   = (group_th + 1) * group_size;
	if (begin >= 0 && begin < flist.size() && end > flist.size()) {
		end = flist.size();
	}
	else if (begin < 0 || begin >= flist.size() || end < 1 || end > flist.size()) {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("ERROR : Group format has error(2)! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}

	fGroup = std::make_pair(group_th, group_size);
	for (int it = begin; it < end; it++) {
		fFileList.push_back(flist.at(it));
	}

	std::cout << "\n---- Loading Root Files ----\n";
	std::cout << Form("Group : %ld th   [%ld files/group],    Total of Load Files : %ld \n", fGroup.first, fGroup.second, fFileList.size());
	for (Long64_t  it = 0; it < fFileList.size(); it++) {
		std::cout << Form("    Number : %ld,   %s\n", it, fFileList.at(it).c_str());
	}
	// end load data with group

	// start read source file list
	bool stagedonly = true;
	unsigned int timeout = 10;
	fChain = new AMSChain("AMSRoot");
	Int_t fileStatus = fChain->AddFromFile(file_list.c_str(), begin, end, stagedonly, timeout);
	if (fileStatus == -1) {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("ROOT file list cannot be opend! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}

	std::cout << "FileStatus : "<< fileStatus << std::endl ;
	std::cout << "Totally : " << fChain->GetEntries() << " data events.\n";
	// end read source file list

	std::cout << "\n**-------------------------------------------**\n";
	std::cout << "\n**    Read Data Form Source File List End    **\n";
	std::cout << "\n**-------------------------------------------**\n";
}

void YiNtuple::saveInputFileList(TFile * file) {
	if (file == 0 || fFileList.size() == 0) return;
	file->cd();
	TTree * tree = new TTree("fileList", "List Of Input File Info");
	std::string filePath;
	tree->Branch("file", &filePath);
	for (Int_t i = 0; i < fFileList.size(); ++i) {
		filePath = fFileList.at(i);
		tree->Fill();
	}
	file->cd();
}

void YiNtuple::loopEventChain() {
	std::cout << "\n**-----------------------------**\n";
	std::cout << "\n**    Loop Event Chain Info    **\n";
	std::cout << "\n**-----------------------------**\n";

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
	if (fChain->GetEntries() <= 0) {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("ERROR : Don't have event! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}

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
		if (!isMC && !isBT && !isISS) {
			MgntSys::Error(LocAddr(), MgntSys::MESSAGE("Event type (ISS, BT, MC) is failed! Exiting ..."));
			MgntSys::Exit(EXIT_FAILURE);
		}
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
			fTimer.stop();

			const UInt_t MemSize = 1024;
			ProcInfo_t procinfo;
			gSystem->GetProcInfo(&procinfo);
			Long64_t memRes = procinfo.fMemResident / MemSize;
			Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
			
			std::cout << Form("Info :: %lf %\n", 100. * float(nprocessed)/float(loop_entries));
			std::cout << Form("        Processed       : %ld / %ld\n", nprocessed, loop_entries);
			std::cout << Form("        Passed          : %ld / %ld\n", npassed, nprocessed);
			std::cout << Form("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
			std::cout << Form("        Real Time       : %9.2f (second)\n", fTimer.time());
			std::cout << Form("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fTimer.time());
			std::cout << Form("        Cpu    System   : %4.1f %\n", procinfo.fCpuSys);
			std::cout << Form("               User     : %4.1f %\n", procinfo.fCpuUser);
			std::cout << Form("        Memory Resident : %2ld GB %4ld MB\n", memRes / MemSize, memRes % MemSize);
			std::cout << Form("               Virtual  : %2ld GB %4ld MB\n", memVrl / MemSize, memVrl % MemSize);
		}
		nprocessed++;

		AMSEventR * event = fChain->GetEvent(ientry);

		//if (nprocessed > 10000) break; // testcode
		
		fRunTagOp->processEvent(event, fChain);

		int preselectEventStatus = fData->preselectEvent(event);
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
		fTimer.stop();
		const UInt_t MemSize = 1024;
		ProcInfo_t procinfo;
		gSystem->GetProcInfo(&procinfo);
		Long64_t memRes = procinfo.fMemResident / MemSize;
		Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
		
		std::cout << Form("Info :: %lf %\n", 100. * float(nprocessed)/float(loop_entries));
		std::cout << Form("        Processed       : %ld / %ld\n", nprocessed, loop_entries);
		std::cout << Form("        Passed          : %ld / %ld\n", npassed, nprocessed);
		std::cout << Form("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
		std::cout << Form("        Real Time       : %9.2f (second)\n", fTimer.time());
		std::cout << Form("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fTimer.time());
		std::cout << Form("        Cpu    System   : %4.1f %\n", procinfo.fCpuSys);
		std::cout << Form("               User     : %4.1f %\n", procinfo.fCpuUser);
		std::cout << Form("        Memory Resident : %2ld GB %4ld MB\n", memRes / MemSize, memRes % MemSize);
		std::cout << Form("               Virtual  : %2ld GB %4ld MB\n", memVrl / MemSize, memVrl % MemSize);
		std::cout << Form("Info :: AMSRoot Files Processed Successfully Finished.\n");
	}
	else {
		std::cout << Form("Info :: AMSRoot Files Processed Seems Failed.\n");
		std::cout << Form("        Processed %ld in %ld\n", nprocessed, loop_entries);
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

	std::cout << "\n**----------------------------**\n";
	std::cout << "\n**    Loop Event Chain End    **\n";
	std::cout << "\n**----------------------------**\n";
}
#endif // __YiProdNtuple_TCC__
