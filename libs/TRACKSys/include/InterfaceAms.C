#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_InterfaceAms_C__
#define __TRACKLibs_InterfaceAms_C__


#include "Sys.h"
#include "Math.h"
#include "IonEloss.h"
#include "GmIonEloss.h"
#include "PartInfo.h"
#include "PhySt.h"
#include "MagEnv.h"
#include "MatEnv.h"
#include "Prop.h"
#include "HitSt.h"
#include "PhyFit.h"
#include "InterfaceAms.h"


namespace TrackSys {
namespace InterfaceAms {

 
TrFitPar::Orientation         Event::ArgOrtt       = TrFitPar::Orientation::kDownward;
VirtualHitSt::NoiseController Event::ArgNoiseCtler = VirtualHitSt::NoiseController::ON;
Bool_t                        Event::ArgSwMscat    = PhyArg::OptMscat();
Bool_t                        Event::ArgSwEloss    = PhyArg::OptEloss();

UInt_t     Event::RunID    = 0;
UInt_t     Event::EvID     = 0;
UInt_t     Event::PtID     = 0;
AMSEventR* Event::Ev       = nullptr;
TrTrackR*  Event::Trtk     = nullptr;
BetaHR*    Event::Btah     = nullptr;
RichRingR* Event::Rich     = nullptr;
Bool_t     Event::StatusTk = false;
Bool_t     Event::StatusTf = false;
Bool_t     Event::StatusRh = false;

std::vector<HitStTRK> Event::TkHitIn = std::vector<HitStTRK>();
HitStTRK              Event::TkHitL1 = HitStTRK();
HitStTRK              Event::TkHitL9 = HitStTRK();

std::vector<HitStTRK> Event::TkHitInQ = std::vector<HitStTRK>();
HitStTRK              Event::TkHitL1Q = HitStTRK();
HitStTRK              Event::TkHitL9Q = HitStTRK();

std::vector<HitStTOF> Event::TfHitT  = std::vector<HitStTOF>();
std::vector<HitStTOF> Event::TfHitQ  = std::vector<HitStTOF>();
std::vector<HitStTOF> Event::TfHitTQ = std::vector<HitStTOF>();

HitStRICH             Event::RhHit = HitStRICH();


void Event::Clear() {
    RunID    = 0;
    EvID     = 0;
    PtID     = 0;
    Ev       = nullptr;
    Trtk     = nullptr;
    Btah     = nullptr;
    Rich     = nullptr;
    StatusTk = false;
    StatusTf = false;
    StatusRh = false;
}


Bool_t Event::Load(AMSEventR* event, UInt_t ipart, Bool_t rebuild) {
    MagMgnt::Load();
    MatMgnt::Load();

    Clear(); Init();
    if (event == nullptr || ipart >= event->NParticle()) return false;
    ParticleR* part = event->pParticle(ipart);
    if (part == nullptr) return false;
	Short_t iTrTrack  = (part->iTrTrack()    < 0) ? -1 : part->iTrTrack();
	Short_t iBetaH    = (part->iBetaH()      < 0) ? -1 : part->iBetaH();
	Short_t iRichRing = (part->iRichRing()   < 0) ? -1 : part->iRichRing();
    
    TrTrackR*  trtk = (iTrTrack  >= 0) ? event->pTrTrack(iTrTrack)   : nullptr;
    BetaHR*    btah = (iBetaH    >= 0) ? event->pBetaH(iBetaH)       : nullptr;
    RichRingR* rich = (iRichRing >= 0) ? event->pRichRing(iRichRing) : nullptr;
    if (trtk == nullptr) return false;
    if (rebuild) trtk->iTrTrackPar(1, 0, 23);

    RunID  = event->Run();
    EvID   = event->Event();
    PtID   = ipart;
    Ev     = event;
    Trtk   = trtk;
    Btah   = btah;
    Rich   = rich;

    StatusTk = BulidHitStTRK();
    StatusTf = BulidHitStTOF();
    StatusRh = BulidHitStRICH();
    if (!StatusTk) { Clear(); Init(); return false; }

    return true;
} 


TrFitPar Event::Get(const PartInfo& info, const TkOpt& tkOpt, const TfOpt& tfOpt, const RhOpt& rhOpt) {
    if (Ev == nullptr) return TrFitPar();
    if (Trtk == nullptr || !StatusTk) return TrFitPar();
    Bool_t hasL1 = (tkOpt.dedx() ? (TkHitL1Q.scx() || TkHitL1Q.scy()) : (TkHitL1.scx() || TkHitL1.scy()));
    Bool_t hasL9 = (tkOpt.dedx() ? (TkHitL9Q.scx() || TkHitL9Q.scy()) : (TkHitL9.scx() || TkHitL9.scy()));
    if (tkOpt.reqL1() && !hasL1) return TrFitPar();
    if (tkOpt.reqL9() && !hasL9) return TrFitPar();

    if (tfOpt.used() && !StatusTf) return TrFitPar();
    if (tfOpt.used()) {
        if (Btah == nullptr) return TrFitPar();
        if (( tfOpt.time() && !tfOpt.dedx()) && TfHitT.size()  == 0) return TrFitPar(); 
        if ((!tfOpt.time() &&  tfOpt.dedx()) && TfHitQ.size()  == 0) return TrFitPar(); 
        if (( tfOpt.time() &&  tfOpt.dedx()) && TfHitTQ.size() == 0) return TrFitPar(); 
    }

    if (rhOpt.used() && !StatusRh) return TrFitPar();
    if (rhOpt.used()) {
        if (Rich == nullptr) return TrFitPar();
        if (RhOpt::Rad::AGL == rhOpt.rad())
            if (HitStRICH::Radiator::AGL != RhHit.rad()) return TrFitPar();
        if (RhOpt::Rad::NAF == rhOpt.rad())
            if (HitStRICH::Radiator::NAF != RhHit.rad()) return TrFitPar();
        if (!RhHit.sib()) return TrFitPar(); 
    }

    TrFitPar fitPar(info, ArgOrtt, ArgNoiseCtler, ArgSwMscat, ArgSwEloss);

    fitPar.add_hit(tkOpt.dedx() ? TkHitInQ : TkHitIn);
    if (tkOpt.useL1()) fitPar.add_hit(tkOpt.dedx() ? TkHitL1Q : TkHitL1);
    if (tkOpt.useL9()) fitPar.add_hit(tkOpt.dedx() ? TkHitL9Q : TkHitL9);

    if (tfOpt.used()) {
        if      ( tfOpt.time() && !tfOpt.dedx()) fitPar.add_hit(TfHitT);
        else if (!tfOpt.time() &&  tfOpt.dedx()) fitPar.add_hit(TfHitQ);
        else if ( tfOpt.time() &&  tfOpt.dedx()) fitPar.add_hit(TfHitTQ);
        else return TrFitPar();
    }
    
    if (rhOpt.used()) fitPar.add_hit(RhHit);

    if (fitPar.check()) return fitPar;
    else                return TrFitPar();
}


void Event::Init() {
    TkHitIn.clear();
    TkHitL1 = HitStTRK();
    TkHitL9 = HitStTRK();
    
    TkHitInQ.clear();
    TkHitL1Q = HitStTRK();
    TkHitL9Q = HitStTRK();
    
    TfHitT.clear();
    TfHitQ.clear();
    TfHitTQ.clear();
    
    RhHit = HitStRICH();
}


Bool_t Event::BulidHitStTRK() {
    if (Trtk == nullptr) return false;
    Short_t cntX = 0;
    Short_t cntY = 0;

    Bool_t hasL1 = Trtk->TestHitLayerJ(1);
    Bool_t hasL9 = Trtk->TestHitLayerJ(9);

    for (Short_t layJ = 1; layJ <= 9; ++layJ) {
	    if (!Trtk->TestHitLayerJ(layJ)) continue;
		TrRecHitR* recHit = Trtk->GetHitLJ(layJ);
		if (recHit == nullptr) continue;
        
        //AMSPoint pnt = ((layJ == 1 || layJ == 9) ? (Trtk->GetHitCooLJ(layJ, 0) + Trtk->GetHitCooLJ(layJ, 1)) * Numc::HALF : Trtk->GetHitCooLJ(layJ)); // (CIEMAT+PG)/2
        AMSPoint pnt = TrTrackR::FitCoo[layJ-1]; // (CIEMAT+PG)/2 after TrTrackR maxspan refit 23

        TrClusterR* xcls = (recHit->GetXClusterIndex() >= 0 && recHit->GetXCluster()) ? recHit->GetXCluster() : nullptr;
		TrClusterR* ycls = (recHit->GetYClusterIndex() >= 0 && recHit->GetYCluster()) ? recHit->GetYCluster() : nullptr;

        Int_t    qopt = ((Ev->NMCEventg() > 0 && Ev->Version() >= 1107) ? InterfaceAms::TkOpt::QOptMC : InterfaceAms::TkOpt::QOptISS);
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
        Bool_t scx = (nsrx > 0);
        Bool_t scy = (nsry > 0);

        HitStTRK hit(scx, scy, layJ);
        hit.set_coo(pnt.x(), pnt.y(), pnt.z());
        hit.set_nsr(nsrx, nsry);
        
        if      (layJ == 1) TkHitL1 = hit;
        else if (layJ == 9) TkHitL9 = hit;
        else                TkHitIn.push_back(hit);
        
        HitStTRK hitQ(scx, scy, layJ);
        hitQ.set_coo(pnt.x(), pnt.y(), pnt.z());
        hitQ.set_nsr(nsrx, nsry);
        hitQ.set_q(qx, qy);
        
        if      (layJ == 1) TkHitL1Q = hitQ;
        else if (layJ == 9) TkHitL9Q = hitQ;
        else                TkHitInQ.push_back(hitQ);

        Bool_t isInn = (layJ != 1 && layJ != 9);
        if (isInn && scx) cntX++;
        if (isInn && scy) cntY++;
    }
    if (cntX <= 2 || cntY <= 3) return false;
    
    return true;
}


Bool_t Event::BulidHitStTOF() {
    if (Btah == nullptr) return false;

    Double_t minT = Numc::ZERO<>;
    std::vector<Short_t>  lay;
    std::vector<SVecD<3>> coo;
    std::vector<Double_t> tme;
    std::vector<Double_t> chg;
    for (Int_t il = 0; il < 4; ++il) {
        if (!Btah->TestExistHL(il)) continue;
        TofClusterHR* cls = Btah->GetClusterHL(il);
        if (cls == nullptr) continue;
        if (!cls->IsGoodTime()) continue;
        lay.push_back(il);
        coo.push_back(SVecD<3>(cls->Coo[0], cls->Coo[1], cls->Coo[2]));
        tme.push_back(Btah->GetTime(il)*HitStTOF::TRANS_NS_TO_CM);
        chg.push_back(Btah->GetQL(il, 2, InterfaceAms::TfOpt::QOpt));
        minT = std::min(minT, tme.back());
    }
    for (UInt_t it = 0; it < lay.size(); ++it) tme.at(it) -= minT;

    for (UInt_t it = 0; it < lay.size(); ++it) {
        HitStTOF hitT(lay.at(it));
        hitT.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        hitT.set_t(tme.at(it));
        TfHitT.push_back(hitT);
        
        HitStTOF hitQ(lay.at(it));
        hitQ.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        hitQ.set_q(chg.at(it));
        TfHitQ.push_back(hitQ);
        
        HitStTOF hitTQ(lay.at(it));
        hitTQ.set_coo(coo.at(it)(0), coo.at(it)(1), coo.at(it)(2));
        hitTQ.set_t(tme.at(it));
        hitTQ.set_q(chg.at(it));
        TfHitTQ.push_back(hitTQ);
    }
    if (TfHitT.size()  <= 2) return false;
    if (TfHitQ.size()  <= 2) return false;
    if (TfHitTQ.size() <= 2) return false;

    return true;
}


Bool_t Event::BulidHitStRICH() {
    if (Rich == nullptr) return false;
    HitStRICH::Radiator radiator = (Rich->IsNaF() ? HitStRICH::Radiator::NAF : HitStRICH::Radiator::AGL);
    Double_t ibta = ((Numc::Compare(Rich->getBeta()) > 0) ? (Numc::ONE<> / Rich->getBeta()) : Numc::ZERO<>);
    SVecD<3> coo(Rich->AMSTrPars[0], Rich->AMSTrPars[1], Rich->AMSTrPars[2]);

    HitStRICH hit(radiator);
    hit.set_coo(coo(0), coo(1), coo(2));
    hit.set_ib(ibta);

    RhHit = hit;
    return true;
}


} // namespace InterfaceAms
} // namespace TrackSys


#endif // __TRACKLibs_InterfaceAms_C__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
