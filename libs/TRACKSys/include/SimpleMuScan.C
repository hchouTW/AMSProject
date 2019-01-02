#ifndef __TRACKLibs_SimpleMuScan_C__
#define __TRACKLibs_SimpleMuScan_C__


#include "Sys.h"
#include "Math.h"
#include "CooMeas.h"
#include "TmeMeas.h"
#include "CherenkovMeas.h"
#include "IonEloss.h"
#include "IonTrEloss.h"
#include "PartInfo.h"
#include "PhySt.h"
#include "MagEnv.h"
#include "MatEnv.h"
#include "Prop.h"
#include "HitSt.h"
#include "TrFitPar.h"
#include "SimpleBtaFit.h"
#include "PhyBtaFit.h"
#include "SimpleTrFit.h"
#include "PhyTrFit.h"
#include "SimpleMuScan.h"


namespace TrackSys {

    
SimpleMuScan& SimpleMuScan::operator=(const SimpleMuScan& rhs) {
    if (this != &rhs) {
        succ_ = rhs.succ_;
        part_ = rhs.part_;
        tsft_ = rhs.tsft_;
        args_ = rhs.args_;
      
        ndof_    = rhs.ndof_;
        nchi_    = rhs.nchi_;
        quality_ = rhs.quality_;

        timer_   = rhs.timer_;
    }
    return *this;
}


void SimpleMuScan::clear() {
    succ_ = false;
    part_.reset(part_.info());
    tsft_ = 0;
    
    args_.clear();

    ndof_.fill(0);
    nchi_.fill(0);
    quality_.fill(0);

    timer_.clear();
}


SimpleMuScan::SimpleMuScan(const TrFitPar& fitPar) {
    Short_t chrg = std::abs(fitPar.info().chrg());
    SimpleMuScan::clear();
    
    part_.reset(chrg, fitPar.info().mass());
    part_.arg().reset(fitPar.sw_mscat(), fitPar.sw_eloss());
    if (part_.chrg() == 0 || part_.chrg() >= LIST_MASS_Q.size()) { clear(); return; }
    
    timer_.start();

    // scan
    MuScanObj condMuObj;
    for (auto&& mass : LIST_MASS_Q.at(part_.chrg())) {
        MuScanObj&& obj = scan(fitPar, mass, true);
        if (obj.chrg() == 0) continue;
        condMuObj = obj;
        break;
    }
    if (condMuObj.chrg() == 0) { clear(); return; }

    // re-scan
    for (Short_t iter = 1; iter <= LMT_SCAN; ++iter) {
        Double_t preM = condMuObj.mass();
        condMuObj = scan(fitPar, condMuObj.mass());
        Bool_t succ_scan = (condMuObj.chrg() != 0);
        if (!succ_scan) { clear(); return; }
        part_.reset(part_.chrg(), condMuObj.mass());
        
        Double_t aftM = condMuObj.mass();
        Double_t difM = std::fabs(preM - aftM);
        Double_t ratM = (difM / (preM + aftM));
        Bool_t   conv = (Numc::Compare(difM, CONVG_DM) < 0 || Numc::Compare(ratM, CONVG_RM) < 0);
        if (conv) break;
    }

    // final fit
    TrFitPar trPar(part_.info(), fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    trPar.add_hit( fitPar.hitsTRK()  );
    trPar.add_hit( fitPar.hitsTOF()  );
    trPar.add_hit( fitPar.hitsRICH() );
    trPar.add_hit( fitPar.hitsTRD()  );
    
    PhyTrFit trfit(trPar);
    if (!trfit.status()) { clear(); return; }
    
    succ_ = true;
    part_ = trfit.part();
    tsft_ = trfit.tsft();
    args_ = trfit.args();
    
    ndof_.at(0) = trfit.ndof(0);
    ndof_.at(1) = trfit.ndof(1);
    nchi_.at(0) = trfit.nchi(0);
    nchi_.at(1) = trfit.nchi(1);
    quality_.at(0) = trfit.quality(0);
    quality_.at(1) = trfit.quality(1);
    
    timer_.stop();

    //if (!succ_) CERR("FAILURE === SimpleMuScan\n");
}


SimpleMuScan::MuScanObj SimpleMuScan::scan(const TrFitPar& fitPar, Double_t mass, Bool_t selopt) {
    PartInfo info(part_.chrg(), mass);
        
    // track fit
    TrFitPar trPar(info, fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    for (auto&& hitTRK : fitPar.hitsTRK()) {
        if (!(hitTRK.scx() || hitTRK.scy())) continue;
        HitStTRK hit(hitTRK.scx(), hitTRK.scy(), hitTRK.lay());
        hit.set_coo(hitTRK.cx(), hitTRK.cy(), hitTRK.cz());
        trPar.add_hit(hit);
    }
    PhyTrFit trfit(trPar);
    if (!trfit.status()) return MuScanObj();
    
    Double_t qltr = std::hypot(trfit.quality(0), trfit.quality(1));
    if (trfit.quality(0) < 0 || trfit.quality(1) < 0)
        qltr = std::max(trfit.quality(0), trfit.quality(1));
    if (selopt && qltr > LMT_QLTR) return MuScanObj();

    // beta fit
    TrFitPar btaPar(info, fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    btaPar.add_hit( fitPar.hitsTRK()  );
    btaPar.add_hit( fitPar.hitsTOF()  );
    btaPar.add_hit( fitPar.hitsRICH() );
    btaPar.add_hit( fitPar.hitsTRD()  );
    PhyBtaFit btafit(btaPar, trfit.part());
    if (!btafit.status()) return MuScanObj();
    
    Double_t qltb = btafit.quality();
    if (selopt && qltb > LMT_QLTB) return MuScanObj();
        
    // selection
    Double_t qlt = std::hypot(qltr, qltb);
    if (qltr < 0 || qltb < 0) qlt = std::max(qltr, qltb);
    if (selopt && qlt > LMT_QLT) return MuScanObj();
    
    // estimate
    auto mmz = std::minmax({ trfit.part().cz(), btafit.part().cz() });
    Double_t tagz = ((fitPar.ortt() == TrFitPar::Orientation::kDownward) ? mmz.first : mmz.second);
    PhySt&& sttTr  = trfit .interpolate_to_z(tagz);
    PhySt&& sttBta = btafit.interpolate_to_z(tagz);
    if (Numc::EqualToZero(sttTr .mom())) return MuScanObj();
    if (Numc::EqualToZero(sttBta.mom())) return MuScanObj();

    Short_t  est_chrg = info.chrg();
    Double_t est_mass = (sttTr.mom() * sttBta.igb());
    Double_t est_qltr = qltr;
    Double_t est_qltb = qltb;
    Double_t est_qlt  = qlt;
    
    if (Numc::Compare(est_mass, LMT_MASS) <= 0) est_mass = LMT_MASS;
    MuScanObj muObj(est_chrg, est_mass, est_qltr, est_qltb, est_qlt);
    return muObj;
}


} // namespace TrackSys


#endif // __TRACKLibs_SimpleMuScan_C__
