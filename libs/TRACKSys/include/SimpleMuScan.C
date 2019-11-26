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
    
        ndof_all_    = rhs.ndof_all_;
        nchi_all_    = rhs.nchi_all_;
        quality_all_ = rhs.quality_all_;

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
    
    ndof_all_ = 0;
    nchi_all_ = 0;
    quality_all_ = 0;

    timer_.clear();
}


SimpleMuScan::SimpleMuScan(const TrFitPar& fitPar) {
    Short_t chrg = std::abs(fitPar.info().chrg());
    SimpleMuScan::clear();
    
    part_.reset(fitPar.info().chrg(), fitPar.info().mass());
    part_.arg().reset(fitPar.sw_mscat(), fitPar.sw_eloss());
    if (part_.chrg() == 0 || std::abs(part_.chrg()) >= LIST_MASS_Q.size()) { clear(); return; }
    
    timer_.start();

    // scan
    Short_t  condID  = -1;
    Double_t condQLT = -1;
    std::vector<PhyTrFit> vecTrs;
    for (auto&& mass : LIST_MASS_Q.at(chrg)) {
        TrFitPar trpar(PartInfo(part_.chrg(), mass), fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
        trpar.add_hit( fitPar.hitsTRK()  );
        trpar.add_hit( fitPar.hitsTOF()  );
        trpar.add_hit( fitPar.hitsRICH() );
        trpar.add_hit( fitPar.hitsTRD()  );

        PhyTrFit trAll(trpar);
        if (!trAll.status()) continue;
        
        Double_t qlt = trAll.quality_all();
        if (condID < 0 || qlt < condQLT) { condID = vecTrs.size(); condQLT = qlt; }
        vecTrs.push_back(trAll);
    }
    if (vecTrs.size() == 0 || condID < 0) { vecTrs.clear(); clear(); return; }

    MuScanObj condMuObj = scan(static_cast<const TrFitPar&>(vecTrs.at(condID)));
    if (condMuObj.chrg() == 0) { vecTrs.clear(); clear(); return; }

    // final fit
    TrFitPar trPar(PartInfo(part_.chrg(), condMuObj.mass()), fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    trPar.add_hit( fitPar.hitsTRK()  );
    trPar.add_hit( fitPar.hitsTOF()  );
    trPar.add_hit( fitPar.hitsRICH() );
    trPar.add_hit( fitPar.hitsTRD()  );
    
    PhyTrFit trfit(trPar);
    if (!trfit.status()) { vecTrs.clear(); clear(); return; }
    
    succ_ = true;
    part_ = trfit.part();
    ibta_ = (trfit.part().ibta() / condMuObj.ibta()) * condMuObj.paribta();
    tsft_ = trfit.tsft();
    args_ = trfit.args();
    
    ndof_.at(0) = trfit.ndof(0);
    ndof_.at(1) = trfit.ndof(1);
    nchi_.at(0) = trfit.nchi(0);
    nchi_.at(1) = trfit.nchi(1);
    quality_.at(0) = trfit.quality(0);
    quality_.at(1) = trfit.quality(1);
        
    ndof_all_    = trfit.ndof_all();
    nchi_all_    = trfit.nchi_all();
    quality_all_ = trfit.quality_all();
    
    timer_.stop();

    //if (!succ_) CERR("FAILURE === SimpleMuScan\n");
}


SimpleMuScan::MuScanObj SimpleMuScan::scan(const TrFitPar& fitPar) {
    // track fit
    TrFitPar trPar(fitPar.info(), fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    for (auto&& hitTRK : fitPar.hitsTRK()) {
        if (!(hitTRK.scx() || hitTRK.scy())) continue;
        HitStTRK hit(hitTRK.scx(), hitTRK.scy(), hitTRK.lay());
        hit.set_coo(hitTRK.cx(), hitTRK.cy(), hitTRK.cz());
        trPar.add_hit(hit);
    }
    PhyTrFit trfit(trPar);
    if (!trfit.status()) return MuScanObj();
    Double_t qltr = trfit.quality_all();

    // beta fit
    TrFitPar btaPar(fitPar.info(), fitPar.ortt(), fitPar.sw_mscat(), fitPar.sw_eloss());
    btaPar.add_hit( fitPar.hitsTRK()  );
    btaPar.add_hit( fitPar.hitsTOF()  );
    btaPar.add_hit( fitPar.hitsRICH() );
    btaPar.add_hit( fitPar.hitsTRD()  );
    PhyBtaFit btafit(btaPar, trfit.part());
    if (!btafit.status()) return MuScanObj();
    Double_t qltb = btafit.quality();
        
    // estimate
    auto mmz = std::minmax({ trfit.part().cz(), btafit.part().cz() });
    Double_t tagz = ((fitPar.ortt() == TrFitPar::Orientation::kDownward) ? mmz.first : mmz.second);
    PhySt&& sttTr  = trfit .interpolate_to_z(tagz);
    PhySt&& sttBta = btafit.interpolate_to_z(tagz);
    if (Numc::EqualToZero(sttTr .mom())) return MuScanObj();
    if (Numc::EqualToZero(sttBta.mom())) return MuScanObj();
    Double_t init_ibta = btafit.part().ibta();

    Short_t  est_chrg = fitPar.info().chrg();
    Double_t est_qltr = qltr;
    Double_t est_qltb = qltb;

    Double_t est_mom     = sttTr.mom();
    Double_t est_paribta = (sttBta.ibta() / init_ibta) * btafit.ibta();
    Double_t est_ibta    = ((est_paribta > LMTL_IBTA_APPROX_LIGHT) ? est_paribta : LMTL_IBTA_APPROX_LIGHT);
    Double_t est_sqrm    = (est_mom * (est_paribta + Numc::ONE<>)) * (est_mom * (est_ibta - Numc::ONE<>));
    if (!Numc::Valid(est_sqrm)) return MuScanObj();
    
    Double_t est_mass = std::sqrt(std::fabs(est_sqrm));
    //if (Numc::Compare(est_mass, EL_MASS) <= 0) est_mass = EL_MASS;
    if (Numc::Compare(est_mass, EL_MASS) <= 0 || est_sqrm <= 0.0) est_mass = EL_MASS; // testcode
    
    Short_t est_sign = (Numc::Compare(est_sqrm) >= 0 ? 1 : -1);
    est_sqrm = est_sign * (est_mass * est_mass);

    MuScanObj muObj(est_chrg, est_mass, est_qltr, est_qltb, est_mom, est_ibta, est_sqrm, est_paribta);
    return muObj;
}


} // namespace TrackSys


#endif // __TRACKLibs_SimpleMuScan_C__
