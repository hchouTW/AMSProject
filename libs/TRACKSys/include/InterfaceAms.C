#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_C__
#define __TRACKLibs_InterfaceAms_C__


namespace TrackSys {
namespace InterfaceAms {


static std::vector<HitStTRK> GetHitStTRK(TrTrackR& trtk, Int_t patt, Bool_t hasQ) {
    const int qopt = TrClusterR::kTotSign2017 | TrClusterR::kSimAsym | TrClusterR::kSimSignal | TrClusterR::kLoss | TrClusterR::kAngle;
    Bool_t hasL1 = (patt == 0 || patt == 5 || patt == 7);
    Bool_t hasL9 = (patt == 0 || patt == 6 || patt == 7);
    std::vector<HitStTRK> hits;
	
    std::vector<Short_t>  lay;
    std::vector<SVecS<2>> nsr;
    std::vector<SVecD<3>> coo;
    std::vector<SVecD<2>> chg;
    for (Int_t layJ = 1; layJ <= 9; ++layJ) {
	    if (!trtk.TestHitLayerJ(layJ)) continue;
		TrRecHitR* recHit = trtk.GetHitLJ(layJ);
		if (recHit == nullptr) continue;
        if (layJ == 1 && !hasL1) continue; 
        if (layJ == 9 && !hasL9) continue; 
        
        AMSPoint    pnt  = ((layJ==1 || layJ==9) ? (trtk.GetHitCooLJ(layJ, 0) + trtk.GetHitCooLJ(layJ, 1)) * Numc::HALF : trtk.GetHitCooLJ(layJ)); // (CIEMAT+PG)/2
        TrClusterR* xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : nullptr;
		TrClusterR* ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : nullptr;

		Double_t qx = (xcls == nullptr || !TrCharge::GoodChargeReconHit(recHit, 0)) ? -1.0 : recHit->GetSignalCombination(0, qopt, 1, 0, 0); 
		Double_t qy = (ycls == nullptr || !TrCharge::GoodChargeReconHit(recHit, 1)) ? -1.0 : recHit->GetSignalCombination(1, qopt, 1, 0, 0); 
			
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


static std::vector<HitStTOF> GetHitStTOF(BetaHR& betaH, Bool_t hasT, Bool_t hasQ) {
    const int qopt = TofRecH::kThetaCor|TofRecH::kBirkCor|TofRecH::kReAttCor|TofRecH::kDAWeight|TofRecH::kQ2Q;
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
        tme.push_back(betaH.GetTime(il));
        chg.push_back(betaH.GetQL(il, 2, qopt));
        minT = std::min(minT, tme.back());
    }

    for (UInt_t it = 0; it < lay.size(); ++it) {
        tme.at(it) -= minT;
        HitStTOF hit(lay.at(it));
        hit.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        if (hasT) hit.set_t(tme.at(it));
        if (hasQ) hit.set_q(chg.at(it));
        hits.push_back(hit);
    }
    return hits;
}


} // namespace InterfaceAms
} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_C__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
