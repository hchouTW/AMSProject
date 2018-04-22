#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_C__
#define __TRACKLibs_InterfaceAms_C__


namespace TrackSys {
namespace InterfaceAms {

static std::vector<HitStTRK> GetHitStTRK(const TrTrackR& trtk, Int_t patt, Bool_t hasQ) {
    std::vector<HitStTRK> hits;
	
    std::vector<Short_t>  lay;
    std::vector<SVecI<2>> nsr;
    std::vector<SVecD<3>> coo;
    std::vector<SVecD<2>> chg;
    for (Int_t layJ = 0; layJ <= 9; ++layJ) {
	    if (!trtk.TestHitLayerJ(layJ)) continue;
		TrRecHitR* recHit = trtk.GetHitLJ(layJ);
		if (recHit == nullptr) continue;
        TrClusterR* xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : nullptr;
		TrClusterR* ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : nullptr;

        AMSPoint loc = (layJ==1 || layJ==9) ? (trtk.GetHitCooLJ(layJ, 0) + trtk.GetHitCooLJ(layJ, 1))*Numc::HALF : trtk.GetHitCooLJ(layJ); // (CIEMAT+PG)/2

        lay.push_back(layJ);
        coo.push_back(SVecD<3>(loc.x(), loc.y(), loc.z());

    }
    
    for (UInt_t it = 0; it < lay.size(); ++it) {
        Bool_t scx = (nsr.at(it)(0) > 0);
        Bool_t scy = (nsr.at(it)(1) > 0);
        HitStTOF hit(scx, scy, lay.at(it));
        hit.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        hit.set_nsr(nsr.at(it)(0), nsr.at(it)(1));
        hit.set_q(chg.at(it)(0), chg.at(it)(1));
        hits.push_back(hit);
    }

    return hits;
}

static std::vector<HitStTOF> GetHitStTOF(const BetaHR& betaH, Bool_t hasQ, Bool_t hasT) {
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
        minT = std::min(minT, betaH.GetTime(il));
    }

    for (UInt_t it = 0; it < lay.size(); ++it) {
        Bool_t scx = (lay.at(it) == 1 || lay.at(it) == 2);
        Bool_t scy = (lay.at(it) == 0 || lay.at(it) == 3);
        tme.at(it) -= minT;
        HitStTOF hit(scx, scy, lay.at(it));
        hit.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        hit.set_t(tme.at(it));
        hit.set_q(chg.at(it));
        hits.push_back(hit);
    }
    return hits;
}


} // namespace InterfaceAms
} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_C__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
