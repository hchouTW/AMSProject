#ifndef __TRACKLibs_SimpleMuScan_C__
#define __TRACKLibs_SimpleMuScan_C__


#include "Sys.h"
#include "Math.h"
#include "TmeMeas.h"
#include "IonEloss.h"
#include "GmIonEloss.h"
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
    }
    return *this;
}


void SimpleMuScan::clear() {
    succ_ = false;
    part_.reset(part_.info());
    tsft_ = 0;
    
    args_.clear();

    ndof_.at(0) = 0;
    ndof_.at(1) = 0;
    nchi_.at(0) = 0;
    nchi_.at(1) = 0;
    quality_.at(0) = 0;
    quality_.at(1) = 0;
}


SimpleMuScan::SimpleMuScan(const TrFitPar& fitPar) {
    Short_t chrg = std::abs(fitPar.info().chrg());
    SimpleMuScan::clear();
    
    part_.reset(chrg, fitPar.info().mass());
    part_.arg().reset(fitPar.sw_mscat(), fitPar.sw_eloss());
    if (part_.chrg() == 0 || part_.chrg() >= LIST_MASS_Q.size()) { clear(); return; }
    if (part_.chrg() > 2) { clear(); return; } // testcode

    // scan
    MuScanObj condMuObj;
    for (auto&& mass : LIST_MASS_Q.at(part_.chrg())) {
        MuScanObj&& obj = scan(fitPar, mass);
        if (obj.chrg() == 0) continue;
        if (obj.qltx() > LMT_QLTX) continue;
        if (obj.qlty() > LMT_QLTY) continue;
        if (obj.qltb() > LMT_QLTB) continue;
        condMuObj = obj;
        break;
    }
    if (condMuObj.chrg() == 0) { clear(); return; } 
    part_.reset(part_.chrg(), condMuObj.mass());

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

    //if (!succ_) CERR("FAILURE === SimpleMuScan\n");
}


SimpleMuScan::MuScanObj SimpleMuScan::scan(const TrFitPar& fitPar, Double_t mass) {
    PartInfo info(part_.chrg(), mass);
        
    // track fit
    TrFitPar trPar(info, fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    for (auto&& hitTRK : fitPar.hitsTRK()) {
        if (!(hitTRK.scx() || hitTRK.scy())) continue;
        HitStTRK hit(hitTRK.scx(), hitTRK.scy(), hitTRK.lay(), hitTRK.isInnTr());
        hit.set_coo(hitTRK.cx(), hitTRK.cy(), hitTRK.cz());
        hit.set_nsr(hitTRK.nsrx(), hitTRK.nsry());
        trPar.add_hit(hit);
    }
    PhyTrFit trfit(trPar);
    if (!trfit.status()) return MuScanObj();
    
    // beta fit
    TrFitPar btaPar(info, fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    btaPar.add_hit( fitPar.hitsTRK()  );
    btaPar.add_hit( fitPar.hitsTOF()  );
    btaPar.add_hit( fitPar.hitsRICH() );
    btaPar.add_hit( fitPar.hitsTRD()  );
    PhyBtaFit btafit(btaPar, trfit.part());
    if (!btafit.status()) return MuScanObj();
    
    // estimate
    auto mmz = std::minmax({ trfit.part().cz(), btafit.part().cz() });
    Double_t tagz = ((fitPar.ortt() == TrFitPar::Orientation::kDownward) ? mmz.first : mmz.second);
    PhySt&& sttTr  = trfit .interpolate_to_z(tagz);
    PhySt&& sttBta = btafit.interpolate_to_z(tagz);
    if (Numc::EqualToZero(sttTr .mom())) return MuScanObj();
    if (Numc::EqualToZero(sttBta.mom())) return MuScanObj();

    Short_t  est_chrg = info.chrg();
    Double_t est_mass = (sttTr.mom() * sttBta.igb());
    Double_t est_qltx = trfit.quality(0);
    Double_t est_qlty = trfit.quality(1);
    Double_t est_qltb = btafit.quality();
    
    if (Numc::Compare(est_mass, LMT_MASS) <= 0) est_mass = LMT_MASS;
    MuScanObj muObj(est_chrg, est_mass, est_qltx, est_qlty, est_qltb);
    return muObj;
}


} // namespace TrackSys


#endif // __TRACKLibs_SimpleMuScan_C__
