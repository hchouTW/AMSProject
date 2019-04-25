#if defined(_PGTRACK_) || defined(__ROOTSHAREDLIBRARY__)
#ifndef __TRACKLibs_RichAms_C__
#define __TRACKLibs_RichAms_C__


#include "Sys.h"
#include "Math.h"
#include "CherenkovMeas.h"
#include "RichAms.h"


namespace TrackSys {


RichHitAms::RichHitAms(RichHitR* hit, double dbta, double rbtaA, double rbtaB, double dist) : RichHitAms() {
    if (hit == nullptr || (dbta <= 0 && rbtaA <= 0 && rbtaB <= 0)) return;
    status_ = true;
    chann_  = hit->Channel;
    pmtid_  = chann_ / NUM_CHANN_IN_2D_PMT;
    pixel_  = chann_ - pmtid_ * NUM_CHANN_IN_2D_PMT;
    locid_[0] = (pixel_ / NUM_CHANN_IN_1D_PMT);
    locid_[1] = (pixel_ % NUM_CHANN_IN_1D_PMT);

    type_  = (dbta > 0) + (rbtaA > 0) * 2 + (rbtaB > 0) * 4;
    dbta_  = (type_&1 == 1) ? dbta  : -1.0;
    rbtaA_ = (type_&2 == 2) ? rbtaA : -1.0;
    rbtaB_ = (type_&4 == 4) ? rbtaB : -1.0;
    npe_   = hit->Npe;
    cx_    = hit->Coo[0];
    cy_    = hit->Coo[1];
    cz_    = hit->Coo[2];
    dist_  = dist;
    hit_   = hit;
    
    cross_ = hit->IsCrossed();
}


CherenkovHit RichHitAms::trans() const {
    if (!status_) return CherenkovHit();
    return CherenkovHit(chann_, pmtid_, dbta_, rbtaA_, rbtaB_, npe_, cx_, cy_);
}


void RichHitAms::clear() {
    status_   = false;
    chann_    = -1; 
    pmtid_    = -1;
    pixel_    = -1;
    locid_[0] = -1;
    locid_[1] = -1;
    type_     = -1;
    dbta_     = -1;
    rbtaA_    = -1;
    rbtaB_    = -1;
    npe_      = -1;
    cx_       = 0;
    cy_       = 0;
    cz_       = 0;
    dist_     = 0;
    hit_      = nullptr;
    cross_    = false;
}


RichAms::RichAms(AMSEventR* event, TrTrackR* trtk) : RichAms() {
    if (event == nullptr || event->NRichHit() == 0) return;
    timer_.start();

    TrTrackR* cand_trtk = nullptr;
    if (trtk == nullptr && event->NParticle() != 0 && event->pParticle(0) != nullptr) {
        cand_trtk = event->pParticle(0)->pTrTrack();
    }
    else { cand_trtk = trtk; }
    if (cand_trtk == nullptr) return;
    
    AMSPoint pmtp;
    AMSDir   pmtd;
    cand_trtk->Interpolate(PMT_CZ, pmtp, pmtd);
    pmtp_ = pmtp;

    AMSPoint pnt;
    AMSDir   dir;
    cand_trtk->Interpolate(RichOffline::RICHDB::RICradpos(), pnt, dir);
    RichOffline::TrTrack track(pnt, dir);

    event_ = event;
    trtk_  = cand_trtk;
    
    status_ = build(track);

    timer_.stop();
    if (!status_) clear();
}


void RichAms::clear() {
    status_  = false;
    kind_    = KIND_EMPTY;
    tile_    = -1;
    index_   = 0;
    dist_    = 0;
    dirp_    = AMSPoint();
    dird_    = AMSDir();
    refp_    = AMSPoint();
    refd_    = AMSDir();
    
    radp_    = AMSPoint();
    radd_    = AMSDir();
    pmtp_    = AMSPoint();
    npe_col_ = 0.0;
    bta_crr_ = 1.0;
    
    is_good_geom_ = false;
    is_bad_tile_ = false;

    is_in_pmt_plane_ = false;

    event_   = nullptr;
    trtk_    = nullptr;
    rhhits_.clear();
    chhits_.clear();
    timer_.clear();
}


CherenkovFit RichAms::fit() {
    if (!status_ || kind_ == KIND_EMPTY || chhits_.size() == 0) return CherenkovFit();
    std::array<double, 2> pmtc({ pmtp_[0], pmtp_[1] });
    if (kind_ == KIND_AGL) {
        return CherenkovFit(chhits_, pmtc, index_, AGL_BETA_WIDTH, AGL_BETA_PDF, bta_crr_);
    }
    else if (kind_ == KIND_NAF) {
        return CherenkovFit(chhits_, pmtc, index_, NAF_BETA_WIDTH, NAF_BETA_PDF, bta_crr_);
    }
    return CherenkovFit();
}
        

RichObjAms RichAms::get_obj_by_fit(int zin) {
    if (!status_ || kind_ == KIND_EMPTY || chhits_.size() == 0 || zin <= 0) return RichObjAms();
    CherenkovFit&& fit = this->fit();

    RichObjAms obj;
    obj.status = status_;
    obj.kind   = kind_;
    obj.tile   = tile_;
    obj.index  = index_;
    obj.dist   = dist_;
    
    obj.is_good_geom = is_good_geom_;
    obj.is_bad_tile  = is_bad_tile_;
   
    obj.bta_crr  = bta_crr_;

    obj.radp[0]  = radp_[0];
    obj.radp[1]  = radp_[1];
    obj.radp[2]  = radp_[2];
    obj.radd[0]  = radd_[0];
    obj.radd[1]  = radd_[1];
    obj.radd[2]  = radd_[2];
    obj.pmtp[0]  = pmtp_[0];
    obj.pmtp[1]  = pmtp_[1];
    obj.pmtp[2]  = pmtp_[2];

    if (!fit.status()) return obj;
    obj.nstn = fit.stns().size(); 
    obj.ncld = fit.clds().size(); 
    obj.ntmr = fit.tmrs().size(); 
    obj.ngst = fit.gsts().size(); 

    obj.nhit_ttl = fit.nhit_total();
    obj.nhit_stn = fit.nhit_stone();
    obj.nhit_cld = fit.nhit_cloud();
    obj.nhit_tmr = fit.nhit_tumor();
    obj.nhit_gst = fit.nhit_ghost();
    obj.nhit_oth = fit.nhit_other();
    
    obj.npe_ttl = fit.npe_total();
    obj.npe_stn = fit.npe_stone();
    obj.npe_cld = fit.npe_cloud();
    obj.npe_tmr = fit.npe_tumor();
    obj.npe_gst = fit.npe_ghost();
    obj.npe_oth = fit.npe_other();

    if (fit.stns().size() != 0) {
        auto& stn = fit.stns().at(0);
        obj.stn_status = stn.status();
        obj.stn_nhit   = stn.nhit();
        obj.stn_npmt   = stn.npmt();
        obj.stn_cx     = stn.cx();
        obj.stn_cy     = stn.cy();
        obj.stn_dist   = stn.dist();
        obj.stn_npe    = stn.npe();
        obj.stn_cnt    = stn.cnt();
        obj.stn_nchi   = stn.nchi();
        obj.stn_chic   = stn.chic();
    }
    
    if (fit.clds().size() != 0) {
        auto& cld = fit.clds().at(0);
        double expnpe = RichRingR::ComputeNpExp(trtk_, cld.cbta(), zin);
        auto&& trace  = cal_trace(cld.cbta(), &cld);
        
        obj.cld_status   = cld.status();
        obj.cld_nhit     = cld.nhit();
        obj.cld_npmt     = cld.npmt();
        obj.cld_border   = trace[0];
        obj.cld_trace    = trace[1];
        obj.cld_accuracy = trace[2];
        obj.cld_uniform  = trace[3];
        obj.cld_crrch    = trace[4];
        obj.cld_beta     = cld.beta();
        obj.cld_cbta     = cld.cbta();
        obj.cld_npe      = cld.npe();
        obj.cld_expnpe   = expnpe;
        obj.cld_cnt      = cld.cnt();
        obj.cld_nchi     = cld.nchi();
        obj.cld_misjudge = cld.misjudge();
    }

    if (fit.hits().size() != 0) {
        for (auto&& hit : fit.hits()) {
            obj.hit_chann.push_back(hit.chann());
            obj.hit_pmtid.push_back(hit.pmtid());
            obj.hit_cls.push_back(hit.cluster());
            obj.hit_mode.push_back(hit.mode());
            obj.hit_beta.push_back(hit.beta());
            obj.hit_npe.push_back(hit.npe());
            obj.hit_cx.push_back(hit.cx());
            obj.hit_cy.push_back(hit.cy());
            obj.nhit++;
        }
    }

    obj.cpuTime = fit.timer().time() * 1.0e+03;
    return obj;
}
        

std::array<double, 5> RichAms::cal_trace(double cbta, const CherenkovCloud* cloud) const {
    const std::array<double, 2> RAD_HEIGHT { 2.5, 0.5 }; // AGL, NAF
    std::array<double, 5> rlt_trace({ 0.0, 0.0, 1.0, 1.0 });
    bool has_cloud = (cloud != nullptr && cloud->status());
    if (!status_ || kind_ == KIND_EMPTY || (has_cloud ? cloud->cbta() : cbta) <= (Numc::ONE<> / index_)) return rlt_trace;

    double cosv = Numc::ONE<> / (index_ * (has_cloud ? cloud->cbta() : cbta));
    double sinv = std::sqrt((Numc::ONE<> - cosv) * (Numc::ONE<> + cosv));
    if (!Numc::Valid(cosv) || !Numc::Valid(sinv)) return rlt_trace;
    int kind = (kind_ == KIND_AGL) ? 0 : 1;

    double lmtrfr = -std::cos(std::asin(Numc::ONE<> / index_));

    SVecD<3> radp(radp_[0], radp_[1], radp_[2]);
    SVecD<3> radd(radd_[0], radd_[1], radd_[2]);
    OrthCoord orth(radd);
    
    double topx = radp[0] + (Numc::ONE_TO_TWO * RAD_HEIGHT[kind]) * (radd[0] / radd[2]);
    double topy = radp[1] + (Numc::ONE_TO_TWO * RAD_HEIGHT[kind]) * (radd[1] / radd[2]);
    double topz = radp[2] + (Numc::ONE_TO_TWO * RAD_HEIGHT[kind]);

    std::array<std::array<short,  TRACE_NSET>, TRACE_NPHI> border;  // 0 no,  1 pass
    std::array<std::array<short,  TRACE_NSET>, TRACE_NPHI> simphs;  // 0 no,  1 dir,  2 rfl
    std::array<std::array<double, TRACE_NSET>, TRACE_NPHI> simphsx;
    std::array<std::array<double, TRACE_NSET>, TRACE_NPHI> simphsy;
    for(auto&& phs : border) phs.fill(0);
    for(auto&& phs : simphs) phs.fill(0);
    for(auto&& phsx : simphsx) phsx.fill(0.0);
    for(auto&& phsy : simphsy) phsy.fill(0.0);

    for (int it = 0; it < TRACE_NPHI; ++it) {
    for (int is = 0; is < TRACE_NSET; ++is) {
        double locwgt = (static_cast<double>(is) + Numc::ONE_TO_TWO) / static_cast<double>(TRACE_NSET);
        
        double phi = Numc::TWO_PI * static_cast<double>(it) / static_cast<double>(TRACE_NPHI);
        SVecD<3>&& phivec = orth.tau() * std::cos(phi) + orth.rho() * std::sin(phi);
        SVecD<3>&& srcvec = radd * cosv + phivec * sinv;
        if (srcvec[2] >= lmtrfr) continue;

        double planex = topx + (-RAD_HEIGHT[kind]) * (locwgt * (radd[0] / radd[2]) + (Numc::ONE<> - locwgt) * (srcvec[0] / srcvec[2]));
        double planey = topy + (-RAD_HEIGHT[kind]) * (locwgt * (radd[1] / radd[2]) + (Numc::ONE<> - locwgt) * (srcvec[1] / srcvec[2]));
        double planez = topz + (-RAD_HEIGHT[kind]);
        if (std::hypot(planex, planey) > MIRROR_TOP_RADIUS) continue;
        
        short tile = RichRingR::getTileIndex(planex, planey);
        if (tile < 0 || tile != tile_) continue;
        border[it][is] = 1; 

        double cosvac = Numc::NEG<> * std::sqrt(Numc::ONE<> - (index_ * index_) * (Numc::ONE<> - srcvec[2]) * (Numc::ONE<> + srcvec[2]));
        double sinvac = std::sqrt((Numc::ONE<> - cosvac) * (Numc::ONE<> + cosvac));
        if (!Numc::Valid(cosvac) || !Numc::Valid(sinvac)) continue;
        SVecD<2>&& vacr = LA::Unit(SVecD<2>(srcvec[0], srcvec[1])) * sinvac;
        if (!Numc::Valid(vacr[0]) || !Numc::Valid(vacr[1])) vacr = std::move(SVecD<2>());

        SVecD<3> vacp(planex, planey, planez);
        SVecD<3> vacd(vacr[0], vacr[1], cosvac);

        double dheight = (PMT_CZ - vacp[2]);
        double pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        double pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);
        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) { 
            if (RichOffline::RichPMTsManager::FindChannel(pmtx, pmty) < 0) continue;
            simphs[it][is]  = 1; 
            simphsx[it][is] = pmtx; 
            simphsy[it][is] = pmty; 
            continue; 
        }

        double vactx = (vacd[0] / vacd[2]);
        double vacty = (vacd[1] / vacd[2]);
        double mirtr = ((MIRROR_TOP_RADIUS - MIRROR_BTM_RADIUS) / MIRROR_HEIGHT);
        double poly0 = (mirtr * mirtr - (vactx * vactx + vacty * vacty));
        double poly1 = Numc::TWO<> * (MIRROR_TOP_RADIUS * mirtr - (vacp[0] * vactx + vacp[1] * vacty));
        double poly2 = (MIRROR_TOP_RADIUS * MIRROR_TOP_RADIUS - (vacp[0] * vacp[0] + vacp[1] * vacp[1]));

        double det  = std::sqrt(poly1 * poly1 - Numc::FOUR<> * poly0 * poly2);
        double sol1 = Numc::ONE_TO_TWO * (-poly1 + det) / poly0;
        double sol2 = Numc::ONE_TO_TWO * (-poly1 - det) / poly0;
        bool is_ok1 = (Numc::Valid(sol1) && sol1 < 0.0 && sol1 > dheight);
        bool is_ok2 = (Numc::Valid(sol2) && sol2 < 0.0 && sol2 > dheight);

        if (!is_ok1 && !is_ok2) continue;

        double sol = 0;
        if (is_ok1) sol = sol1;
        if (is_ok2) sol = sol2;
        vacp[0] = vacp[0] + sol * vactx;
        vacp[1] = vacp[1] + sol * vacty;
        vacp[2] = vacp[2] + sol;
        
        double cosmir = Numc::NEG<> * std::sqrt((mirtr * mirtr) / (Numc::ONE<> + mirtr * mirtr));
        double sinmir = std::sqrt(Numc::ONE<> - cosmir * cosmir);
        SVecD<2>&& mirr = LA::Unit(SVecD<2>(-vacp[0], -vacp[1])) * sinmir;
        SVecD<3> mird(mirr[0], mirr[1], cosmir);

        double cosrfl = std::fabs(LA::Dot(vacd, mird));
        if (!Numc::Valid(cosrfl)) continue;
        
        vacd = std::move(LA::Unit(vacd + Numc::TWO<> * cosrfl * mird));
        if (vacd[2] >= Numc::ZERO<>) continue;

        dheight = (PMT_CZ - vacp[2]);
        pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);

        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) {
            if (RichOffline::RichPMTsManager::FindChannel(pmtx, pmty) < 0) continue;
            simphs[it][is]  = 2;
            simphsx[it][is] = pmtx;
            simphsy[it][is] = pmty;
        }
    }}
    
    short count_bd = 0;
    for(auto&& phs : border) {
    for (auto&& ph : phs) {
        count_bd += (ph != 0);
    }}
    double trace_bd = static_cast<double>(count_bd) / static_cast<double>(TRACE_NPHI * TRACE_NSET);
    
    short count = 0;
    for(auto&& phs : simphs) {
    for (auto&& ph : phs) {
        count += (ph != 0);
    }}
    double trace = static_cast<double>(count) / static_cast<double>(TRACE_NPHI * TRACE_NSET);
    double crrch = (Numc::Compare(trace) > 0) ? (0.1 * ((index_ * index_ - Numc::ONE<>) / (sinv * sinv)) * std::fabs(radd[2]) / trace) : 0.0;

    if (!Numc::Valid(trace_bd) || !Numc::Valid(trace) || !Numc::Valid(crrch)) return rlt_trace;

    if (!has_cloud) {
        rlt_trace = std::move(std::array<double, 5>({ trace_bd, trace, 1.0, 1.0, crrch }));
        return rlt_trace;
    }

    const double thres = 3.4; // pmt width
    const std::vector<CherenkovHit>& hits = cloud->hits();
    std::vector<short>  hphi(hits.size(), -1);
    std::vector<double> hdxy(hits.size(), -1.0);
    std::vector<short>  hpmt(hits.size(), -1);
    for (int is = 0; is < TRACE_NSET; ++is) {
    for (int iphi = 0; iphi < TRACE_NPHI; ++iphi) {
        if (simphs[iphi][is] == 0) continue;
        for (int ih = 0; ih < hphi.size(); ++ih) {
            const CherenkovHit& hit = hits.at(ih);
            if (simphs[iphi][is] == 1 && (hit.mode() != 0)) continue;
            if (simphs[iphi][is] == 2 && (hit.mode() != 1 && hit.mode() != 2)) continue;
            
            double dist = std::hypot(hit.cx() - simphsx[iphi][is], hit.cy() - simphsy[iphi][is]);
            if (dist > thres || (hdxy.at(ih) >= 0.0 && hdxy.at(ih) < dist)) continue;
            
            hphi.at(ih) = iphi;
            hdxy.at(ih) = dist;
            hpmt.at(ih) = hit.pmtid();
        }
    }}

    for (int ih = hphi.size() - 1; ih >= 0; --ih) {
        if (hphi.at(ih) >= 0) continue;
        hphi.erase(hphi.begin() + ih);
        hpmt.erase(hpmt.begin() + ih);
    }
    std::set<short> setpmt;
    for (int ih = 0; ih < hpmt.size(); ++ih) setpmt.insert(hpmt.at(ih));
    
    double accuracy = static_cast<double>(hphi.size()) / static_cast<double>(hits.size());

    if (hphi.size() < 2 || setpmt.size() < 2) {
        rlt_trace = std::move(std::array<double, 5>({ trace_bd, trace, accuracy, 0.0, crrch }));
        return rlt_trace;
    }

    int uniformity_width = static_cast<int>(Numc::ONE_TO_TWO * static_cast<double>(count) / static_cast<double>(setpmt.size()));
    if (uniformity_width <= 0) uniformity_width = 1;

    std::set<int> withphi;
    for (int ih = 0; ih < hphi.size(); ++ih) {
        int init_phi = hphi.at(ih);
        int cnt_lw = 0;
        for (int iphi = init_phi; iphi > (init_phi - TRACE_NPHI); --iphi) {
            if (cnt_lw > uniformity_width) break;
            int idx_phi = (iphi >= 0) ? iphi : (iphi + TRACE_NPHI);
            for (int is = 0; is < TRACE_NSET; ++is) {
                if (simphs[idx_phi][is] == 0) continue;
                int idx = idx_phi * TRACE_NSET + is;
                withphi.insert(idx);
                cnt_lw++;
            }
        }
        int cnt_up = 0;
        for (int iphi = init_phi; iphi < (init_phi + TRACE_NPHI); ++iphi) {
            if (cnt_up > uniformity_width) break;
            int idx_phi = (iphi < TRACE_NPHI) ? iphi : (iphi - TRACE_NPHI);
            for (int is = 0; is < TRACE_NSET; ++is) {
                if (simphs[idx_phi][is] == 0) continue;
                int idx = idx_phi * TRACE_NSET + is;
                withphi.insert(idx);
                cnt_up++;
            }
        }
    }
    
    double uniformity = static_cast<double>(withphi.size()) / static_cast<double>(count);
    if (!Numc::Valid(uniformity)) {
        rlt_trace = std::move(std::array<double, 5>({ trace_bd, trace, accuracy, 1.0, crrch }));
        return rlt_trace;
    }

    rlt_trace = std::move(std::array<double, 5>({ trace_bd, trace, accuracy, uniformity, crrch }));
    return rlt_trace;
}
        

bool RichAms::build(RichOffline::TrTrack track) {
    const std::array<double, 2> RAD_HEIGHT { 2.5, 0.5 }; // AGL, NAF
    
    RichOffline::RichRadiatorTileManager crossed_tile(&track);
    if(crossed_tile.getkind() == KIND_EMPTY) return false;
    
    kind_   = static_cast<short>(crossed_tile.getkind());
    tile_   = static_cast<short>(crossed_tile.getcurrenttile());
    index_  = static_cast<double>(crossed_tile.getindex());
    dist_   = static_cast<double>(crossed_tile.getdistance());
    dirp_   = crossed_tile.getemissionpoint();
    dird_   = crossed_tile.getemissiondir();
    refp_   = crossed_tile.getemissionpoint(1);
    refd_   = crossed_tile.getemissiondir(1);

    // Charge to AMS coordinate
    AMSPoint amsp = RichOffline::RichAlignment::RichToAMS(dirp_);
    AMSDir   amsd = RichOffline::RichAlignment::RichToAMS(dird_);
    if (std::cos(amsd.gettheta()) > 0) amsd = AMSDir(-amsd[0], -amsd[1], -amsd[2]);
    radp_ = amsp;
    radd_ = amsd;

    if (RichBetaUniformityCorrection::getHead() != nullptr) {
        RichBetaUniformityCorrection* ptr = RichBetaUniformityCorrection::getHead();
        double bta_crr = ptr->getCorrection(amsp[0], amsp[1], amsd[0], amsd[1]);
        if (bta_crr > 0) bta_crr_ = bta_crr;
    }

    double radius = std::hypot(radp_[0], radp_[1]);
    is_good_geom_ = (radius < EXTERNAL_RAD_RADIUS);
    if (is_good_geom_ && kind_ == KIND_AGL) {
        is_good_geom_ = (std::max(std::fabs(radp_[0]), std::fabs(radp_[1])) > RAD_BOUNDARY[0]);
    }
    if (is_good_geom_ && kind_ == KIND_NAF) {
        is_good_geom_ = (std::max(std::fabs(radp_[0]), std::fabs(radp_[1])) < RAD_BOUNDARY[1]);
    }
    tile_ = RichRingR::getTileIndex(radp_[0], radp_[1]);
    if (tile_ < 0) return false;
        
    bool is_bad_tile = false;
    std::array<short, 7>  BAD_TILE_INDEX { 3, 7, 12, 20, 87, 100, 108 };
    for (int idx = 0; idx < BAD_TILE_INDEX.size(); ++idx) {
        if (tile_ == BAD_TILE_INDEX[idx]) is_bad_tile = true;
    }
    is_bad_tile_ = is_bad_tile;

    bool is_in_pmt_plane = false;
    {
        int cntin = 0;
        for (int ir = 1; ir <= 10; ++ir) {
            double dr = WIDTH_PMT * 0.15 * static_cast<double>(ir);
            for (int iphi = 0; iphi < 10; ++iphi) {
                double phi = Numc::TWO_PI * 0.1 * static_cast<double>(iphi);
                double dx  = dr * std::cos(phi);
                double dy  = dr * std::sin(phi);
                if (RichOffline::RichPMTsManager::FindChannel(pmtp_[0] + dx, pmtp_[1] + dy) < 0) continue;
                cntin++;
            }
        }
        double cross = 0.01 * static_cast<double>(cntin);
        if (cross > 0.75) is_in_pmt_plane = true;
    }
    is_in_pmt_plane_ = is_in_pmt_plane;

    int bit = event_->nRichRing();

    const double betamin = 0.01;
    const double betamax = 5.00;
    for(RichOffline::RichRawEvent* rawhit = new RichOffline::RichRawEvent(event_); rawhit != nullptr; rawhit = rawhit->next()) {
        if(!rawhit->getbit(RichOffline::ok_status_bit)) continue;
        RichHitR* phit = rawhit->getpointer();

        float recs[3] = { 0., 0., 0. };
        rawhit->reconstruct(dirp_, refp_, dird_, refd_, betamin, betamax, recs, index_, kind_);

        if (recs[0] <= 0 && recs[1] <= 0 && recs[2] <= 0) continue;
        double dbta  = (recs[0] > 0) ? recs[0] : -1.0;
        double rbtaA = (recs[1] > 0) ? recs[1] : -1.0;
        double rbtaB = (recs[2] > 0) ? recs[2] : -1.0;
        
        if (phit == nullptr || phit->Npe <= 0) continue;
        double dist = std::hypot(phit->Coo[0] - pmtp_[0], phit->Coo[1] - pmtp_[1]);

        RichHitAms rhhit(phit, dbta, rbtaA, rbtaB, dist);
        if (!rhhit.status()) continue;
        
        CherenkovHit&& chhit = rhhit.trans();
        if (!chhit.status()) continue;
        
        rhhits_.push_back(rhhit);
        chhits_.push_back(chhit);
    }
    if (rhhits_.size() > 1) std::sort(rhhits_.begin(), rhhits_.end(), RichHitAms_sort());
    if (chhits_.size() > 1) std::sort(chhits_.begin(), chhits_.end(), CherenkovHit_sort());

    for (auto&& hit : rhhits_) {
        if (hit.cross()) continue;
        npe_col_ += hit.npe();
    }

    return true;
}
        

std::vector<std::array<double, 2>> RichAms::RayTrace(const std::array<double, 6>& part, double cbta, double index, double height, short tile) {
    std::vector<std::array<double, 2>> trace;
    if (cbta <= 0.0 || index <= 1.0 || height < 0.0 || tile < 0) return trace;

    double cosv = Numc::ONE<> / (index * cbta);
    double sinv = std::sqrt((Numc::ONE<> - cosv) * (Numc::ONE<> + cosv));
    if (!Numc::Valid(cosv) || !Numc::Valid(sinv)) return trace;

    double lmtrfr = -std::cos(std::asin(Numc::ONE<> / index));

    SVecD<3> partp(part[0], part[1], part[2]);
    SVecD<3> partd(part[3], part[4], part[5]);
    OrthCoord orth(partd);
    
    double topx = partp[0] + (Numc::ONE_TO_TWO * height) * (partd[0] / partd[2]);
    double topy = partp[1] + (Numc::ONE_TO_TWO * height) * (partd[1] / partd[2]);
    double topz = partp[2] + (Numc::ONE_TO_TWO * height);

    const int trace_nphi = 180;
    const int trace_nset = 5;
    for (int it = 0; it < trace_nphi; ++it) {
    for (int is = 0; is < trace_nset; ++is) {
        double locwgt = (static_cast<double>(is) + Numc::ONE_TO_TWO) / static_cast<double>(trace_nset);
        
        double phi = Numc::TWO_PI * static_cast<double>(it) / static_cast<double>(trace_nphi);
        SVecD<3>&& phivec = orth.tau() * std::cos(phi) + orth.rho() * std::sin(phi);
        SVecD<3>&& srcvec = partd * cosv + phivec * sinv;
        if (srcvec[2] >= lmtrfr) continue;

        double planex = topx + (-height) * (locwgt * (partd[0] / partd[2]) + (Numc::ONE<> - locwgt) * (srcvec[0] / srcvec[2]));
        double planey = topy + (-height) * (locwgt * (partd[1] / partd[2]) + (Numc::ONE<> - locwgt) * (srcvec[1] / srcvec[2]));
        double planez = topz + (-height);
        if (std::hypot(planex, planey) > MIRROR_TOP_RADIUS) continue;
        
        short otile = RichRingR::getTileIndex(planex, planey);
        if (otile < 0 || otile != tile) continue;

        double cosvac = Numc::NEG<> * std::sqrt(Numc::ONE<> - (index * index) * (Numc::ONE<> - srcvec[2]) * (Numc::ONE<> + srcvec[2]));
        double sinvac = std::sqrt((Numc::ONE<> - cosvac) * (Numc::ONE<> + cosvac));
        if (!Numc::Valid(cosvac) || !Numc::Valid(sinvac)) continue;
        SVecD<2>&& vacr = LA::Unit(SVecD<2>(srcvec[0], srcvec[1])) * sinvac;
        if (!Numc::Valid(vacr[0]) || !Numc::Valid(vacr[1])) vacr = std::move(SVecD<2>());

        SVecD<3> vacp(planex, planey, planez);
        SVecD<3> vacd(vacr[0], vacr[1], cosvac);

        double dheight = (PMT_CZ - vacp[2]);
        double pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        double pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);
        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) { 
            //if (RichOffline::RichPMTsManager::FindChannel(pmtx, pmty) < 0) continue;
            trace.push_back(std::array<double, 2>({pmtx, pmty})); 
            continue; 
        }

        double vactx = (vacd[0] / vacd[2]);
        double vacty = (vacd[1] / vacd[2]);
        double mirtr = ((MIRROR_TOP_RADIUS - MIRROR_BTM_RADIUS) / MIRROR_HEIGHT);
        double poly0 = (mirtr * mirtr - (vactx * vactx + vacty * vacty));
        double poly1 = Numc::TWO<> * (MIRROR_TOP_RADIUS * mirtr - (vacp[0] * vactx + vacp[1] * vacty));
        double poly2 = (MIRROR_TOP_RADIUS * MIRROR_TOP_RADIUS - (vacp[0] * vacp[0] + vacp[1] * vacp[1]));

        double det  = std::sqrt(poly1 * poly1 - Numc::FOUR<> * poly0 * poly2);
        double sol1 = Numc::ONE_TO_TWO * (-poly1 + det) / poly0;
        double sol2 = Numc::ONE_TO_TWO * (-poly1 - det) / poly0;
        bool is_ok1 = (Numc::Valid(sol1) && sol1 < 0.0 && sol1 > dheight);
        bool is_ok2 = (Numc::Valid(sol2) && sol2 < 0.0 && sol2 > dheight);

        if (!is_ok1 && !is_ok2) continue;

        double sol = 0;
        if (is_ok1) sol = sol1;
        if (is_ok2) sol = sol2;
        vacp[0] = vacp[0] + sol * vactx;
        vacp[1] = vacp[1] + sol * vacty;
        vacp[2] = vacp[2] + sol;
        
        double cosmir = Numc::NEG<> * std::sqrt((mirtr * mirtr) / (Numc::ONE<> + mirtr * mirtr));
        double sinmir = std::sqrt(Numc::ONE<> - cosmir * cosmir);
        SVecD<2>&& mirr = LA::Unit(SVecD<2>(-vacp[0], -vacp[1])) * sinmir;
        SVecD<3> mird(mirr[0], mirr[1], cosmir);

        double cosrfl = std::fabs(LA::Dot(vacd, mird));
        if (!Numc::Valid(cosrfl)) continue;
        
        vacd = std::move(LA::Unit(vacd + Numc::TWO<> * cosrfl * mird));
        if (vacd[2] >= Numc::ZERO<>) continue;

        dheight = (PMT_CZ - vacp[2]);
        pmtx = vacp[0] + dheight * (vacd[0] / vacd[2]);
        pmty = vacp[1] + dheight * (vacd[1] / vacd[2]);

        if (std::fabs(pmtx) < PMT_HOLE[0] && std::fabs(pmty) < PMT_HOLE[1]) continue;
        if (std::hypot(pmtx, pmty) < MIRROR_BTM_RADIUS) {
            //if (RichOffline::RichPMTsManager::FindChannel(pmtx, pmty) < 0) continue;
            trace.push_back(std::array<double, 2>({pmtx, pmty})); 
        }
    }}

    return trace;
}


} // namespace TrackSys


#endif // __TRACKLibs_RichAms_C__
#endif // _PGTRACK_ __ROOTSHAREDLIBRARY__ 
