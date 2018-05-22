#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_C__
#define __TRACKLibs_InterfaceAms_C__


namespace TrackSys {
namespace InterfaceAms {


void Event::init() {
    status_ = false;
    withQ_  = false;
    going_  = 0;
    
    event_  = nullptr;
    trtk_   = nullptr;
    btah_   = nullptr;
    rich_   = nullptr;

    //TRK
    trHitIn_.clear();
    trHitL1_ = HitStTRK();
    trHitL9_ = HitStTRK();

    // TOF
    tfHit_.clear();

    // RICH
    rhHit_ = HitStRICH();
}


Event::Event(AMSEventR* event, Bool_t withQ) {
    init();
    if (event == nullptr) return;
    ParticleR* part = event->pParticle(0);
    if (part == nullptr) return;
	Short_t iTrTrack    = (part->iTrTrack()    < 0) ? -1 : part->iTrTrack()   ;
	Short_t iBetaH      = (part->iBetaH()      < 0) ? -1 : part->iBetaH()     ;
	Short_t iRichRing   = (part->iRichRing()   < 0) ? -1 : part->iRichRing()  ;
    if (iTrTrack < 0) return;
    TrTrackR*  trtk = (iTrTrack  >= 0) ? event->pTrTrack(iTrTrack)   : nullptr;
    BetaHR*    btah = (iBetaH    >= 0) ? event->pBetaH(iBetaH)       : nullptr;
    RichRingR* rich = (iRichRing >= 0) ? event->pRichRing(iRichRing) : nullptr;
    if (trtk == nullptr) return;

    status_ = true;
    withQ_  = withQ;
    going_  = -1;

    event_ = event;
    trtk_  = trtk;
    btah_  = btah;
    rich_  = rich;

    if (!bulid_HitStTRK()) { init(); return; }
    bulid_HitStTOF();
    bulid_HitStRICH();
}


Bool_t Event::bulid_HitStTRK() {
    if (trtk_ == nullptr) return false;
    Short_t cntX = 0;
    Short_t cntY = 0;

    Bool_t hasL1 = trtk_->TestHitLayerJ(1);
    Bool_t hasL9 = trtk_->TestHitLayerJ(9);

    std::vector<HitStTRK> hitIn;
    HitStTRK              hitL1;
    HitStTRK              hitL9;

    for (Short_t layJ = 1; layJ <= 9; ++layJ) {
	    if (!trtk_->TestHitLayerJ(layJ)) continue;
		TrRecHitR* recHit = trtk_->GetHitLJ(layJ);
		if (recHit == nullptr) continue;
        
        AMSPoint    pnt  = ((layJ == 1 || layJ == 9) ? (trtk_->GetHitCooLJ(layJ, 0) + trtk_->GetHitCooLJ(layJ, 1)) * Numc::HALF : trtk_->GetHitCooLJ(layJ)); // (CIEMAT+PG)/2
        TrClusterR* xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : nullptr;
		TrClusterR* ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : nullptr;

		Double_t qx = (xcls == nullptr || !TrCharge::GoodChargeReconHit(recHit, 0)) ? -1.0 : recHit->GetSignalCombination(0, QOptTracker, 1, 0, 0); 
		Double_t qy = (ycls == nullptr || !TrCharge::GoodChargeReconHit(recHit, 1)) ? -1.0 : recHit->GetSignalCombination(1, QOptTracker, 1, 0, 0); 
			
        std::vector<float> xstripSig;
		std::vector<float> xstripSgm;
		std::vector<float> xstripSS;
		for (int it = 0; (xcls!=nullptr) && (it < xcls->GetLength()); ++it) {
			xstripSig.push_back(xcls->GetSignal(it));
			xstripSgm.push_back(xcls->GetNoise(it));
            xstripSS.push_back(xcls->GetSignal(it)/xcls->GetNoise(it));
		}
		
		std::vector<float> ystripSig;
		std::vector<float> ystripSgm;
		std::vector<float> ystripSS;
		for (int it = 0; (ycls!=nullptr) && (it < ycls->GetLength()); ++it) {
			ystripSig.push_back(ycls->GetSignal(it));
			ystripSgm.push_back(ycls->GetNoise(it));
            ystripSS.push_back(ycls->GetSignal(it)/ycls->GetNoise(it));
		}
        
        const float GateTh = 0.3;
        const float LenTh  = 0.082;
        short xseedAddr = (xcls) ? xcls->GetSeedAddress() : -1;
        short xseedIndx = (xcls) ? xcls->GetSeedIndex() : -1;
        float xlenTh    = (xcls) ? LenTh*xstripSS.at(xseedIndx) : 0.;
        short xreg[2] = { xseedIndx, xseedIndx };
		for (int it = xseedIndx-1; (xcls!=nullptr) && it >= 0; --it)
            if (xstripSS.at(it) < xstripSS.at(it+1)+GateTh && xstripSS.at(it) > xlenTh) xreg[0] = it; else break;
		for (int it = xseedIndx+1; (xcls!=nullptr) && it < xstripSS.size(); ++it)
            if (xstripSS.at(it) < xstripSS.at(it-1)+GateTh && xstripSS.at(it) > xlenTh) xreg[1] = it; else break;
        
        short yseedAddr = (ycls) ? ycls->GetSeedAddress() : -1;
        short yseedIndx = (ycls) ? ycls->GetSeedIndex() : -1;
        float ylenTh    = (ycls) ? LenTh*ystripSS.at(yseedIndx) : 0.;
        short yreg[2] = { yseedIndx, yseedIndx };
		for (int it = yseedIndx-1; (ycls!=nullptr) && it >= 0; --it)
            if (ystripSS.at(it) < ystripSS.at(it+1)+GateTh && ystripSS.at(it) > ylenTh) yreg[0] = it; else break;
		for (int it = yseedIndx+1; (ycls!=nullptr) && it < ystripSS.size(); ++it)
            if (ystripSS.at(it) < ystripSS.at(it-1)+GateTh && ystripSS.at(it) > ylenTh) yreg[1] = it; else break;
            
        Short_t nsrx = (xcls != nullptr) ? (xreg[1]-xreg[0]+1) : -1;
        Short_t nsry = (ycls != nullptr) ? (yreg[1]-yreg[0]+1) : -1;
        Bool_t scx = (nsrx > 0);
        Bool_t scy = (nsry > 0);

        HitStTRK hit(scx, scy, layJ);
        hit.set_coo(pnt.x(), pnt.y(), pnt.z());
        hit.set_nsr(nsrx, nsry);
        //if (withQ_) hit.set_q(qx, qy);

        if      (layJ == 1) hitL1 = hit;
        else if (layJ == 9) hitL9 = hit;
        else                hitIn.push_back(hit);

        Bool_t isInn = (layJ != 1 && layJ != 9);
        if (isInn && scx) cntX++;
        if (isInn && scy) cntY++;
    }
    if (cntX <= 2 || cntY <= 3) return false;

    if (hasL1 && (hitL1.scx() || hitL1.scy())) trHitL1_ = hitL1;
    if (hasL9 && (hitL9.scx() || hitL9.scy())) trHitL9_ = hitL9;
    if (hitIn.size() > 0) trHitIn_ = hitIn;
    return true;
}


Bool_t Event::bulid_HitStTOF() {
    if (btah_ == nullptr) return false;
    std::vector<HitStTOF> hits;

    Double_t minT = Numc::ZERO<>;
    std::vector<Short_t>  lay;
    std::vector<SVecD<3>> coo;
    std::vector<Double_t> tme;
    std::vector<Double_t> chg;
    for (Int_t il = 0; il < 4; ++il) {
        if (!btah_->TestExistHL(il)) continue;
        TofClusterHR* cls = btah_->GetClusterHL(il);
        if (cls == nullptr) continue;
        if (!cls->IsGoodTime()) continue;
        lay.push_back(il);
        coo.push_back(SVecD<3>(cls->Coo[0], cls->Coo[1], cls->Coo[2]));
        tme.push_back(btah_->GetTime(il)*HitStTOF::TRANS_NS_TO_CM);
        chg.push_back(btah_->GetQL(il, 2, QOptTOF));
        minT = std::min(minT, tme.back());
    }

    for (UInt_t it = 0; it < lay.size(); ++it) {
        tme.at(it) -= minT;
        HitStTOF hit(lay.at(it));
        hit.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        hit.set_t(tme.at(it));
        //if (withQ_) hit.set_q(chg.at(it));
        
        hits.push_back(hit);
    }

    Double_t beta  = btah_->GetBeta();
    Short_t  going = -Numc::Compare(beta);

    if (hits.size() <= 2) return false;
    else {
        tfHit_ = hits;
        going_ = going;
        return true;
    }
}


Bool_t Event::bulid_HitStRICH() {
    if (rich_ == nullptr) return false;
    HitStRICH::Radiator radiator = (rich_->IsNaF() ? HitStRICH::Radiator::NAF : HitStRICH::Radiator::AGL);
    Double_t ibta = ((Numc::Compare(rich_->getBeta()) > 0) ? (Numc::ONE<> / rich_->getBeta()) : Numc::ZERO<>);
    SVecD<3> coo(rich_->AMSTrPars[0], rich_->AMSTrPars[1], rich_->AMSTrPars[2]);

    HitStRICH hit(radiator);
    hit.set_coo(coo(0), coo(1), coo(2));
    hit.set_ib(ibta);

    rhHit_ = hit;
    return true;
}
        

TrFitPar Event::get(const PartInfo& info, const TrackerPatt& trPatt, Bool_t withTOF, Bool_t withRICH) const {
    TrFitPar::Orientation ortt = ((going_ <= 0) ? TrFitPar::Orientation::kDownward : TrFitPar::Orientation::kUpward);
    if (!status_) return TrFitPar();
    if (withTOF  && (btah_ == nullptr || tfHit_.size() == 0)) return TrFitPar();
    if (withRICH && (rich_ == nullptr || !rhHit_.sib()))      return TrFitPar();

    Bool_t reqL1 = (TrackerPatt::InnerL1 == trPatt || TrackerPatt::FullSpan == trPatt);
    Bool_t reqL9 = (TrackerPatt::InnerL9 == trPatt || TrackerPatt::FullSpan == trPatt);
    if (reqL1 && !(trHitL1_.scx() || trHitL1_.scy())) return TrFitPar();
    if (reqL9 && !(trHitL9_.scx() || trHitL9_.scy())) return TrFitPar();
    
    Bool_t hasL1 = ((trHitL1_.scx() || trHitL1_.scy()) && (TrackerPatt::MaxSpan == trPatt || TrackerPatt::InnerL1 == trPatt || TrackerPatt::FullSpan == trPatt));
    Bool_t hasL9 = ((trHitL9_.scx() || trHitL9_.scy()) && (TrackerPatt::MaxSpan == trPatt || TrackerPatt::InnerL9 == trPatt || TrackerPatt::FullSpan == trPatt));

    TrFitPar fitPar(info, ortt);
    fitPar.add_hit(trHitIn_);
    if (hasL1)    fitPar.add_hit(trHitL1_);
    if (hasL9)    fitPar.add_hit(trHitL9_);
    if (withTOF)  fitPar.add_hit(tfHit_);
    if (withRICH) fitPar.add_hit(rhHit_);
    
    if (fitPar.check()) return fitPar;
    else                return TrFitPar();
}












































/*




static std::vector<HitStTRK> GetHitStTRK(TrTrackR& trtk, const TrackerPatt& patt, Bool_t hasQ) {
    Bool_t hasL1 = trtk.TestHitLayerJ(1);
    Bool_t hasL9 = trtk.TestHitLayerJ(9);
    
    std::vector<HitStTRK> hits;
    if (!hasL1 && (TrackerPatt::InnerL1 == patt || TrackerPatt::FullSpan == patt)) return hits;
    if (!hasL9 && (TrackerPatt::InnerL9 == patt || TrackerPatt::FullSpan == patt)) return hits;

    Int_t satLayJ = (TrackerPatt::MaxSpan == patt || TrackerPatt::InnerL1 == patt || TrackerPatt::FullSpan == patt) ? 1 : 2;
    Int_t endLayJ = (TrackerPatt::MaxSpan == patt || TrackerPatt::InnerL9 == patt || TrackerPatt::FullSpan == patt) ? 9 : 8;

    std::vector<Short_t>  lay;
    std::vector<SVecS<2>> nsr;
    std::vector<SVecD<3>> coo;
    std::vector<SVecD<2>> chg;
    for (Int_t layJ = satLayJ; layJ <= endLayJ; ++layJ) {
	    if (!trtk.TestHitLayerJ(layJ)) continue;
		TrRecHitR* recHit = trtk.GetHitLJ(layJ);
		if (recHit == nullptr) continue;
        
        AMSPoint    pnt  = ((layJ == 1 || layJ == 9) ? (trtk.GetHitCooLJ(layJ, 0) + trtk.GetHitCooLJ(layJ, 1)) * Numc::HALF : trtk.GetHitCooLJ(layJ)); // (CIEMAT+PG)/2
        TrClusterR* xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : nullptr;
		TrClusterR* ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : nullptr;

		Double_t qx = (xcls == nullptr || !TrCharge::GoodChargeReconHit(recHit, 0)) ? -1.0 : recHit->GetSignalCombination(0, QOptTracker, 1, 0, 0); 
		Double_t qy = (ycls == nullptr || !TrCharge::GoodChargeReconHit(recHit, 1)) ? -1.0 : recHit->GetSignalCombination(1, QOptTracker, 1, 0, 0); 
			
        std::vector<float> xstripSig;
		std::vector<float> xstripSgm;
		std::vector<float> xstripSS;
		for (int it = 0; (xcls!=nullptr) && (it < xcls->GetLength()); ++it) {
			xstripSig.push_back(xcls->GetSignal(it));
			xstripSgm.push_back(xcls->GetNoise(it));
            xstripSS.push_back(xcls->GetSignal(it)/xcls->GetNoise(it));
		}
		
		std::vector<float> ystripSig;
		std::vector<float> ystripSgm;
		std::vector<float> ystripSS;
		for (int it = 0; (ycls!=nullptr) && (it < ycls->GetLength()); ++it) {
			ystripSig.push_back(ycls->GetSignal(it));
			ystripSgm.push_back(ycls->GetNoise(it));
            ystripSS.push_back(ycls->GetSignal(it)/ycls->GetNoise(it));
		}
        
        const float GateTh = 0.3;
        const float LenTh  = 0.082;
        short xseedAddr = (xcls) ? xcls->GetSeedAddress() : -1;
        short xseedIndx = (xcls) ? xcls->GetSeedIndex() : -1;
        float xlenTh    = (xcls) ? LenTh*xstripSS.at(xseedIndx) : 0.;
        short xreg[2] = { xseedIndx, xseedIndx };
		for (int it = xseedIndx-1; (xcls!=nullptr) && it >= 0; --it)
            if (xstripSS.at(it) < xstripSS.at(it+1)+GateTh && xstripSS.at(it) > xlenTh) xreg[0] = it; else break;
		for (int it = xseedIndx+1; (xcls!=nullptr) && it < xstripSS.size(); ++it)
            if (xstripSS.at(it) < xstripSS.at(it-1)+GateTh && xstripSS.at(it) > xlenTh) xreg[1] = it; else break;
        
        short yseedAddr = (ycls) ? ycls->GetSeedAddress() : -1;
        short yseedIndx = (ycls) ? ycls->GetSeedIndex() : -1;
        float ylenTh    = (ycls) ? LenTh*ystripSS.at(yseedIndx) : 0.;
        short yreg[2] = { yseedIndx, yseedIndx };
		for (int it = yseedIndx-1; (ycls!=nullptr) && it >= 0; --it)
            if (ystripSS.at(it) < ystripSS.at(it+1)+GateTh && ystripSS.at(it) > ylenTh) yreg[0] = it; else break;
		for (int it = yseedIndx+1; (ycls!=nullptr) && it < ystripSS.size(); ++it)
            if (ystripSS.at(it) < ystripSS.at(it-1)+GateTh && ystripSS.at(it) > ylenTh) yreg[1] = it; else break;
            
        Short_t nsrx = (xcls != nullptr) ? (xreg[1]-xreg[0]+1) : -1;
        Short_t nsry = (ycls != nullptr) ? (yreg[1]-yreg[0]+1) : -1;
          
        lay.push_back(layJ);
        nsr.push_back(SVecS<2>(nsrx, nsry));
        coo.push_back(SVecD<3>(pnt.x(), pnt.y(), pnt.z()));
        chg.push_back(SVecD<2>(qx, qy));
    }
    
    for (UInt_t it = 0; it < lay.size(); ++it) {
        Bool_t scx = (nsr.at(it)(0) > 0);
        Bool_t scy = (nsr.at(it)(1) > 0);
        HitStTRK hit(scx, scy, lay.at(it));
        hit.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        hit.set_nsr(nsr.at(it)(0), nsr.at(it)(1));
        //if (hasQ) hit.set_q(chg.at(it)(0), chg.at(it)(1));
        hits.push_back(hit);
    }

    return hits;
}


static std::vector<HitStTOF> GetHitStTOF(BetaHR& betaH, Bool_t hasQ) {
    std::vector<HitStTOF> hits;

    Double_t minT = Numc::ZERO<>;
    std::vector<Short_t>  lay;
    std::vector<SVecD<3>> coo;
    std::vector<Double_t> tme;
    std::vector<Double_t> chg;
    for (Int_t il = 0; il < 4; ++il) {
        if (!betaH.TestExistHL(il)) continue;
        TofClusterHR* cls = betaH.GetClusterHL(il);
        if (cls == nullptr) continue;
        if (!cls->IsGoodTime()) continue;
        lay.push_back(il);
        coo.push_back(SVecD<3>(cls->Coo[0], cls->Coo[1], cls->Coo[2]));
        tme.push_back(betaH.GetTime(il)*HitStTOF::TRANS_NS_TO_CM);
        chg.push_back(betaH.GetQL(il, 2, QOptTOF));
        minT = std::min(minT, tme.back());
    }

    for (UInt_t it = 0; it < lay.size(); ++it) {
        tme.at(it) -= minT;
        HitStTOF hit(lay.at(it));
        hit.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        hit.set_t(tme.at(it));
        //if (hasQ) hit.set_q(chg.at(it));
        hits.push_back(hit);
    }
    return hits;
}


static HitStRICH GetHitStRICH(RichRingR& rich) {
    HitStRICH::Radiator radiator = (rich.IsNaF() ? HitStRICH::Radiator::NAF : HitStRICH::Radiator::AGL);
    Double_t ibta = ((Numc::Compare(rich.getBeta()) > 0) ? (Numc::ONE<> / rich.getBeta()) : Numc::ZERO<>);
    SVecD<3> coo(rich.AMSTrPars[0], rich.AMSTrPars[1], rich.AMSTrPars[2]);

    HitStRICH hit(radiator);
    hit.set_coo(coo(0), coo(1), coo(2));
    hit.set_ib(ibta);

    return hit;
}
*/

} // namespace InterfaceAms
} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_C__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
