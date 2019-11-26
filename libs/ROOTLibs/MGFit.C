#ifndef __ROOTLibs_MGFit_C__
#define __ROOTLibs_MGFit_C__

namespace MGROOT {
namespace Fit {

void RooVar::set(const std::string& name, Hist * samp, Hist * sumt, const HistList& temp, Bool_t link, Double_t min, Double_t max) {
    clear();
 
    if (name == "") {
        MGSys::ShowError("<< RooVar::RooVar >> No name.");
        return;
    }
    
    if (samp == nullptr || !samp->exist()) {
        MGSys::ShowError("<< RooVar::RooVar >> No sample histogram.");
        return;
    }
    else {
        if (samp->dim() != HistDim::k1D) {
            MGSys::ShowError("<< RooVar::RooVar >> Now, only for 1-D histogram.");
        }
    }
    
    if (temp.size() == 0) {
        MGSys::ShowError("<< RooVar::RooVar >> Empty template histogram.");
        return;
    }
    else {
        for (auto&& elem : temp) {
            if (elem == nullptr || !elem->exist()) {
                MGSys::ShowError("<< RooVar::RooVar >> Has empty template histogram.");
                return;
            }
            else if (elem->dim() != samp->dim()) {
                MGSys::ShowError(STR("<< RooVar::RooVar >> Template dimension(%d) vs. Sample dimension(%d).", elem->dim(), samp->dim()));
                return;
            }
        }
    }

    Double_t minEdge = samp->axis().x().min();
    Double_t maxEdge = samp->axis().x().max();
    
    Bool_t norm = (MGNumc::Compare(min, max) < 0);
    fMin = (min < minEdge || !norm) ? minEdge : min;
    fMax = (max > maxEdge || !norm) ? maxEdge : max;
    fMinBin = (*samp)()->FindBin(fMin);
    fMaxBin = (*samp)()->FindBin(fMax);

    fName = name;
    fLink = link;
    if (fLink) {
        fSamp = samp;
        if (sumt != nullptr)
            fSumt = sumt;
        fTemp = temp;
    }
    else {
        fSamp = Hist::New(STR("RVAR_%s__%s", fName.c_str(), samp->name().c_str()), samp->title(), (*samp)());
        if (sumt != nullptr)
            fSumt = Hist::New(STR("RVAR_%s__%s", fName.c_str(), sumt->name().c_str()), sumt->title(), (*sumt)());
        for (auto&& elem : temp) {
            Hist * hist = Hist::New(STR("RVAR_%s__%s", fName.c_str(), elem->name().c_str()), elem->title(), (*elem)());
            fTemp.push_back(hist);
        }
    }
}

void RooVar::set(const std::string& name, Hist * samp, const HistList& temp, Bool_t link, Double_t min, Double_t max) {
    set(name, samp, nullptr, temp, link, min, max);
}


void RooVar::clear() {
    fName = "";
    if (!fLink) Hist::Delete(fSamp);
    fSamp = nullptr;
    if (!fLink) Hist::Delete(fSumt);
    fSumt = nullptr;
    if (!fLink) {
        for (auto&& elem : fTemp)
            Hist::Delete(elem);
    }
    fTemp.clear();

    fLink = true;
    fMin = 0;
    fMax = 0;
    fMinBin = 0;
    fMaxBin = 0;
}



RooResult::RooResult(const RooVar& var, Bool_t extended, Bool_t fluc) : exist_(false) {
    RooResult::fCounter++;
    if (!var.exist()) return;
    const Double_t LMTVAL = 1.0e-8;

    //---- RooVar ----//
    Double_t    min  = var.min();
    Double_t    max  = var.max();
    std::string name = var.name();
    Hist *      sumt = Hist::New(STR("SUMT%06d_%s", RooResult::fCounter, var.samp()->name().c_str()), "", (*var.samp())(), true);
    Hist *      samp = Hist::New(STR("SAMP%06d_%s", RooResult::fCounter, var.samp()->name().c_str()), "", (*var.samp())());
    HistList    temp;
    for (auto&& elem : var.temp()) {
        Hist * hist = Hist::New(STR("TEMP%06d_%s", RooResult::fCounter, elem->name().c_str()), "", (*elem)());

        if (fluc) {
            Double_t sum = (*hist)()->Integral();
            Bool_t opt = (sum < 100.);
            for (Int_t ib = 0; ib <= (*hist)()->GetNbinsX()+1; ++ib) {
                if (MGNumc::Compare((*hist)()->GetBinContent(ib), LMTVAL) <= 0) continue;
                Double_t men = (*hist)()->GetBinContent(ib);
                Double_t err = (*hist)()->GetBinError  (ib);
                Double_t prob = men / sum;
                Double_t sgm = std::sqrt(sum * prob * (1. - prob));
                const Int_t max_ntry = 10; Int_t itry = 1;
                Double_t val = (opt) ? Rndm::Binomial(Int_t(sum+0.5), prob) : Rndm::Gaus(men, sgm);
                while (MGNumc::Compare(val) < 0 && itry <= max_ntry) {
                    val = (opt) ? Rndm::Binomial(Int_t(sum+0.5), prob) : Rndm::Gaus(men, sgm);
                    itry++;
                }
                if (itry > max_ntry) continue;
                (*hist)()->SetBinContent(ib, val);
                
                // TODO (hchou): find Correct way
                (*hist)()->SetBinError  (ib, err * (val/men));
            }
        }
        
        for (Int_t ib = 0; ib <= (*hist)()->GetNbinsX()+1; ++ib) {
            if (MGNumc::Compare((*hist)()->GetBinContent(ib), LMTVAL) <= 0) {
                (*hist)()->SetBinContent(ib, LMTVAL);
                (*hist)()->SetBinError(ib, 0.);
            }
        }
 
        hist->normalized(HistNorm::kIntegral);
        temp.push_back(hist);
    }
    var_.set(name, samp, sumt, temp, true, min, max);

    // MsgLevel  DEBUG=0, INFO=1, PROGRESS=2, WARNING=3, ERROR=4, FATAL=5
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    
    //---- Binning ----//
    Int_t lwBin = (*var_.samp())()->FindBin(var_.min());
    Int_t upBin = (*var_.samp())()->FindBin(var_.max());
    Double_t total = (*var_.samp())()->Integral(lwBin, upBin);
    Double_t lwLmt = 0.0 * std::sqrt(total); // testcode
    //Double_t lwLmt = -7.0 * std::sqrt(total);
    Double_t upLmt = total + 7.0 * std::sqrt(total);
    Double_t meanVal = (total / var_.num());
    
    //---- RooFit ----//
    RooRealVar xvar(var_.name().c_str(), "", var_.min(), var_.max());
    std::vector<RooRealVar*>  tempVar(var_.num(), nullptr);
    RooArgList tempVarList;
    for (Int_t it = 0; it < var_.num(); ++it) {
        tempVar.at(it) = new RooRealVar(CSTR("ROOVAR__TEMP%03d", it), "", meanVal, lwLmt, upLmt);
        tempVarList.add(*tempVar.at(it));
    }
    
    //---- Templates and Model ----//
    std::vector<RooDataHist *> tempHist(var_.num(), nullptr);
    std::vector<RooHistPdf *>  tempPdf(var_.num(), nullptr);
    RooArgList tempPdfList;
    for (Int_t it = 0; it < var_.num(); ++it) {
        tempHist.at(it) = new RooDataHist(CSTR("ROOHIST__TEMP%03d", it), "", xvar, (*var_.temp(it))());
        tempPdf.at(it) = new RooHistPdf(CSTR("ROOPDF__TEMP%03d", it), "", xvar, *tempHist.at(it));
        tempPdfList.add(*tempPdf.at(it));
    }
    RooAddPdf model("model", "model", tempPdfList, tempVarList);

    //---- Sample ----//
    RooDataHist sampHist("ROOHIST__SAMP", "", xvar, (*var_.samp())());
  
    //---- Fit ----//
    model.fitTo(sampHist, RooFit::Extended(extended), RooFit::PrintLevel(-1));

    par_.clear();
    for (Int_t it = 0; it < var_.num(); ++it) {
        for (Int_t ib = 0; ib <= (*var_.temp(it))()->GetNbinsX()+1; ++ib)
            if(MGNumc::Compare((*var_.temp(it))()->GetBinContent(ib), LMTVAL) <= 0)
                (*var_.temp(it))()->SetBinContent(ib, 0.0);
        
        par_.push(tempVar.at(it)->getVal(), tempVar.at(it)->getError()); 

        Double_t cntNum = (*var_.temp(it))()->Integral(lwBin, upBin);
        (*var_.temp(it))()->Scale(par_.val(it) / cntNum);
        (*var_.sumt())()->Add((*var_.temp(it))());
    }
    
    Double_t chi = 0;
    Int_t ndf = -(var_.num() - (!extended));
    for (Int_t ibin = lwBin; ibin <= upBin; ++ibin) {
        Double_t sampVal = (*var_.samp())()->GetBinContent(ibin);
        Double_t sampErr = (*var_.samp())()->GetBinError(ibin);
        if (MGNumc::EqualToZero(sampVal)) continue;
        
        Double_t tempVal = 0.;
        Double_t tempErr = 0.;
        for (Int_t it = 0; it < var.num(); ++it) {
            tempVal += (*var_.temp(it))()->GetBinContent(ibin);
            tempErr += (*var_.temp(it))()->GetBinError(ibin) * (*var_.temp(it))()->GetBinError(ibin);
        }
        tempErr = std::sqrt(tempErr);
        
        Double_t totErr = std::sqrt(sampErr * sampErr + tempErr * tempErr);
        Double_t resVal = (sampVal - tempVal) / totErr;
        
        chi += (resVal * resVal);
        ndf++;
    }
    par_.set_ndf_and_chi(ndf, chi);

    exist_ = true;

    for (Int_t it = 0; it < var.num(); ++it) {
    delete tempHist.at(it); delete tempPdf.at(it); delete tempVar.at(it);
        tempHist.at(it) = nullptr; tempPdf.at(it) = nullptr; tempVar.at(it) = nullptr;
    }
    
    RooMsgService::instance().reset();
}


RooSysResult::RooSysResult(const RooVar& var, Bool_t extended, Int_t ntimes) : exist_(false) {
    exist_ = false;
    RooResult result(var, extended);
    if (!var.exist()) return;
    if (!result.exist()) return;
    if (ntimes <= 0) return;

    const Double_t CIlevel95 = 0.95;
    const Double_t CIlevel68 = 0.68;
    const Double_t CIlevel = CIlevel95;
    const Double_t CIalpha = 0.5 * (1.0 - CIlevel);
    Long64_t nreqs = Long64_t(ntimes / CIlevel + 10);

    std::vector<RooPar> perparvec(nreqs);
    Bool_t FlucOpt = true;

    Long64_t ndf = 0;
    Long64_t cntreq = 0;
    Long64_t iter = 0, niter = 10 * nreqs;
    while (iter < niter && cntreq < nreqs) {
        iter++;
        RooResult rlt(var, extended, FlucOpt);
        if (!rlt.exist()) continue;
        perparvec.at(cntreq) = rlt.par();
        if (ndf==0) ndf = rlt.par().ndf();
        cntreq++;
    }
    perparvec.resize(cntreq);
    if (perparvec.size() == 0) return;
    std::sort(perparvec.begin(), perparvec.end(), RooPar::RooParSort());
    Long64_t startsize = Long64_t(perparvec.size() * (      CIalpha)); 
    Long64_t finalsize = Long64_t(perparvec.size() * (1.0 - CIalpha)); 
    std::vector<RooPar> parvec(perparvec.begin()+startsize, perparvec.begin()+finalsize);
    
    //---- Binning ---//
    Double_t lwBin = (*var.samp())()->FindBin(var.min());
    Double_t upBin = (*var.samp())()->FindBin(var.max());
    Double_t total = (*var.samp())()->Integral(lwBin, upBin);

    Double_t avgchi = 0.;
    std::vector<Double_t> sumwgt(var.num(), 0.0);
    std::vector<std::pair<Double_t, Double_t>> avgpar(var.num(), std::make_pair(0.0, 0.0));
    for (Int_t it = 0; it < parvec.size(); ++it) {
        RooPar& par = parvec.at(it);
        for (Int_t ip = 0; ip < var.num(); ++ip) {
            Double_t err2 = (par.err(ip) * par.err(ip));
            sumwgt.at(ip) += (1. / err2);
            avgpar.at(ip).first  += par.val(ip) / err2;
            avgpar.at(ip).second += err2;
        }
        avgchi += par.chi();
    }
    avgchi /= parvec.size();
    
    for (Int_t ip = 0; ip < var.num(); ++ip) {
        avgpar.at(ip).first  = avgpar.at(ip).first / sumwgt.at(ip);
        avgpar.at(ip).second = std::sqrt(avgpar.at(ip).second / parvec.size());
    }

    std::vector<Double_t> avgerr(var.num(), 0.0);
    for (Int_t it = 0; it < parvec.size(); ++it) {
        RooPar& par = parvec.at(it);
        for (Int_t ip = 0; ip < var.num(); ++ip) {
            Double_t dif = (par.val(ip) - avgpar.at(ip).first) / par.err(ip);
            avgerr.at(ip) += (dif * dif);
        }
    }
    for (Int_t ip = 0; ip < var.num(); ++ip) {
        avgerr.at(ip) = std::sqrt(avgerr.at(ip) / sumwgt.at(ip));
    }

    //---- Parameters ----//
    var_ = result.var();
    std_par_ = result.par();
    sys_par_.clear();
    for (Int_t ip = 0; ip < var.num(); ++ip) {
        Double_t val = avgpar.at(ip).first;
        Double_t err = avgpar.at(ip).second;
        Double_t sys = avgerr.at(ip);
        sys_par_.push(val, err, sys);
       
        // TODO (hchou): find Correct way
        //CERR("PAR(%d) VAL %8.2f STAT_ERR %8.2f SYS_ERR %8.2f\n", ip, val, err, sys);
    }
    sys_par_.set_ndf_and_chi(ndf, avgchi * ndf);
    sys_fit_set_ = parvec;

    exist_ = true;
}


} // namespace Fit
} // namespace MGROOT

#endif // __ROOTLibs_MGFit_C__
