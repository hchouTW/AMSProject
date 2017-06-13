#ifndef __MGFit_C__
#define __MGFit_C__

void MGROOT::Fit::RooVar::set(const std::string& name, Hist * samp, Hist * sumt, const HistList& temp, Bool_t link, Double_t min, Double_t max) {
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
				MGSys::ShowError(StrFmt("<< RooVar::RooVar >> Template dimension(%d) vs. Sample dimension(%d).", elem->dim(), samp->dim()));
				return;
			}
		}
	}

	Double_t minEdge = samp->axis().x().min();
	Double_t maxEdge = samp->axis().x().max();
	
	Bool_t norm = (MGNumc::Compare(min, max) < 0);
	fMin = (min < minEdge || !norm) ? minEdge : min;
	fMax = (max > maxEdge || !norm) ? maxEdge : max;

	fName = name;
	fLink = link;
	if (fLink) {
		fSamp = samp;
		if (sumt != nullptr)
			fSumt = sumt;
		fTemp = temp;
	}
	else {
		fSamp = Hist::New(StrFmt("ROOVAR__%s", samp->name().c_str()), samp->title(), (*samp)());
		if (sumt != nullptr)
			fSumt = Hist::New(StrFmt("ROOVAR__%s", sumt->name().c_str()), sumt->title(), (*sumt)());
		for (auto&& elem : temp) {
			Hist * hist = Hist::New(StrFmt("ROOVAR__%s", elem->name().c_str()), elem->title(), (*elem)());
			fTemp.push_back(hist);
		}
	}
}

void MGROOT::Fit::RooVar::set(const std::string& name, Hist * samp, const HistList& temp, Bool_t link, Double_t min, Double_t max) {
	set(name, samp, nullptr, temp, link, min, max);
}


void MGROOT::Fit::RooVar::clear() {
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
}



MGROOT::Fit::RooResult::RooResult(const RooVar& var, Bool_t extended, Bool_t fluc) : fExist(false) {
	RooResult::COUNT++;
	if (!var.exist()) return;
	const Double_t LMTVAL = 1.0e-8;

	//---- RooVar ----//
	Double_t    min  = var.min();
	Double_t    max  = var.max();
	std::string name = var.name();
	Hist *      sumt = Hist::New(StrFmt("SUMT%06d_%s", RooResult::COUNT, var.samp()->name().c_str()), "", (*var.samp())(), true);
	Hist *      samp = Hist::New(StrFmt("SAMP%06d_%s", RooResult::COUNT, var.samp()->name().c_str()), "", (*var.samp())());
	HistList    temp;
	for (auto&& elem : var.temp()) {
		Hist * hist = Hist::New(StrFmt("TEMP%06d_%s", RooResult::COUNT, elem->name().c_str()), "", (*elem)());
		//////////////////////////////////
		Double_t sum = (*hist)()->Integral();
		Bool_t opt = (sum < 100.);
		for (Int_t ib = 0; ib <= (*hist)()->GetNbinsX()+1; ++ib) {
			if (fluc) {
				Double_t men = (*hist)()->GetBinContent(ib);
				Double_t prob = men / sum;
				Double_t sgm = std::sqrt(sum * prob * (1. - prob));
				const Int_t maxNTry = 10; Int_t itry = 1;
				Double_t val = (opt) ? MGROOT::Random::Binomial(Int_t(sum+0.5), prob) : MGROOT::Random::Gaus(men, sgm);
				while (MGNumc::Compare(val, 0.) < 0 && itry <= maxNTry) {
					val = (opt) ? MGROOT::Random::Binomial(Int_t(sum+0.5), prob) : MGROOT::Random::Gaus(men, sgm);
					itry++;
				}
				if (itry > maxNTry) continue;
				(*hist)()->SetBinContent(ib, val);
				(*hist)()->SetBinContent(ib, std::sqrt(val));
			}
			if(MGNumc::Compare((*hist)()->GetBinContent(ib), LMTVAL) <= 0) {
				(*hist)()->SetBinContent(ib, LMTVAL);
				(*hist)()->SetBinError(ib, 0.);
			}
		}
 
		hist->normalized(HistNorm::kIntegral);
		temp.push_back(hist);
	}
	//fVar.set(name, samp);



/*



	//---- Histogram ----//
	fName = var.name();

	fSumT = MGROOT::Hist::New(StrFmt("SUMT%06d_%s", RooResult::COUNT, var.samp()->name().c_str()), "", (*var.samp())(), true);
	fSamp = MGROOT::Hist::New(StrFmt("SAMP%06d_%s", RooResult::COUNT, var.samp()->name().c_str()), "", (*var.samp())());
	TH1 * hSamp = (*fSamp)();
	TH1 * hSumT = (*fSumT)();

	const Double_t LMTVAL = 1.0e-8;
	std::vector<TH1 *> hTemp(var.num(), nullptr);
	for (Int_t it = 0; it < var.num(); ++it) {
		Hist * temp = MGROOT::Hist::New(StrFmt("TEMP%06d_%s", RooResult::COUNT, var.temp(it)->name().c_str()), "", (*var.temp(it))());
		hTemp.at(it) = (*temp)();

		Double_t sum = hTemp.at(it)->Integral();
		Bool_t opt = (sum < 100.);
		for (Int_t ib = 0; ib <= hTemp.at(it)->GetNbinsX()+1; ++ib) {
			if (fluc) {
				Double_t men = hTemp.at(it)->GetBinContent(ib);
				Double_t prob = men / sum;
				Double_t sgm = std::sqrt(sum * prob * (1. - prob));
				const Int_t maxNTry = 10; Int_t itry = 1;
				Double_t val = (opt) ? MGROOT::Random::Binomial(Int_t(sum+0.5), prob) : MGROOT::Random::Gaus(men, sgm);
				while (MGNumc::Compare(val, 0.) < 0 && itry <= maxNTry) {
					val = (opt) ? MGROOT::Random::Binomial(Int_t(sum+0.5), prob) : MGROOT::Random::Gaus(men, sgm);
					itry++;
				}
				if (itry > maxNTry) continue;
				hTemp.at(it)->SetBinContent(ib, val);
			}
			if(MGNumc::Compare(hTemp.at(it)->GetBinContent(ib), LMTVAL) <= 0)
  		   hTemp.at(it)->SetBinContent(ib, LMTVAL);
		}
 
		temp->normalized(MGROOT::Hist::kIntegral);
		fTemp.push_back(temp);
	}
	fParam.resize(var.num(), std::make_pair(0., 0.));

	// MsgLevel  DEBUG=0, INFO=1, PROGRESS=2, WARNING=3, ERROR=4, FATAL=5
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	
	//---- Binning ----//
	Double_t lwBin = hSamp->FindBin(min());
	Double_t upBin = hSamp->FindBin(max());
	Double_t total = hSamp->Integral(lwBin, upBin);
	Double_t lwLmt = -7.0 * std::sqrt(total);
	Double_t upLmt = total + 7.0 * std::sqrt(total);
	Double_t meanVal = (total / var.num());
	
	//---- RooFit ----//
	RooRealVar xvar(fName.c_str(), "", min(), max());
	std::vector<RooRealVar *>  tempVar(var.num(), nullptr);
	RooArgList tempVarList;
	for (Int_t it = 0; it < var.num(); ++it) {
		tempVar.at(it) = new RooRealVar(CStrFmt("ROOVAR__TEMP%03d", it), "", meanVal, lwLmt, upLmt);
		tempVarList.add(*tempVar.at(it));
	}
	
	//---- Templates and Model ----//
	std::vector<RooDataHist *> tempHist(var.num(), nullptr);
	std::vector<RooHistPdf *>  tempPdf(var.num(), nullptr);
	RooArgList tempPdfList;
	for (Int_t it = 0; it < var.num(); ++it) {
		tempHist.at(it) = new RooDataHist(CStrFmt("ROOHIST__TEMP%03d", it), "", xvar, hTemp.at(it));
		tempPdf.at(it) = new RooHistPdf(CStrFmt("ROOPDF__TEMP%03d", it), "", xvar, *tempHist.at(it));
		tempPdfList.add(*tempPdf.at(it));
	}
	RooAddPdf model("model", "model", tempPdfList, tempVarList);

	//---- Sample ----//
	RooDataHist sampHist("ROOHIST__SAMP", "", xvar, hSamp);
  
	//---- Fit ----//
	model.fitTo(sampHist, RooFit::Extended(extended), RooFit::PrintLevel(-1));
	
	for (Int_t it = 0; it < var.num(); ++it) {
		for (Int_t ib = 0; ib <= hTemp.at(it)->GetNbinsX()+1; ++ib)
			if(MGNumc::Compare(hTemp.at(it)->GetBinContent(ib), LMTVAL) <= 0)
				hTemp.at(it)->SetBinContent(ib, 0.0);
		
		fParam.at(it).first  = tempVar.at(it)->getVal();
		fParam.at(it).second = tempVar.at(it)->getError();
		Double_t cntNum = hTemp.at(it)->Integral(lwBin, upBin);
		hTemp.at(it)->Scale(fParam.at(it).first / cntNum);
		hSumT->Add(hTemp.at(it));
	}
	
	fNdf = -(var.num() - (!extended));
	for (Int_t ibin = lwBin; ibin <= upBin; ++ibin) {
		Double_t sampVal = hSamp->GetBinContent(ibin);
		Double_t sampErr = hSamp->GetBinError(ibin);
		if (MGNumc::EqualToZero(sampVal)) continue;
		
		Double_t tempVal = 0.;
		Double_t tempErr = 0.;
		for (Int_t it = 0; it < var.num(); ++it) {
			tempVal += hTemp.at(it)->GetBinContent(ibin);
			tempErr += hTemp.at(it)->GetBinError(ibin) * hTemp.at(it)->GetBinError(ibin);
		}
		tempErr = std::sqrt(tempErr);
		
		Double_t totErr = std::sqrt(sampErr * sampErr + tempErr * tempErr);
		Double_t resVal = (sampVal - tempVal) / totErr;
		
		fChisq += (resVal * resVal);
		fNdf++;
	}
	fValid = true;

	for (Int_t it = 0; it < var.num(); ++it) {
    delete tempHist.at(it); delete tempPdf.at(it); delete tempVar.at(it);
		tempHist.at(it) = nullptr; tempPdf.at(it) = nullptr; tempVar.at(it) = nullptr;
	}
	
	RooMsgService::instance().reset();
	*/
}
		










/////////////////////////////////////////////////////////////////////
/*
MGROOT::Fit::RooResult::RooResult(MGROOT::Fit::RooParam& param, Bool_t extended, Bool_t fluc) : fValid(false), fName(""), fSamp(nullptr), fSumT(nullptr), fChisq(0.), fNdf(0) {
	RooResult::COUNT++;
	fBound.first  = 0.0;
	fBound.second = 0.0;
	if (!param.valid()) return;
	
	//---- Histogram ----//
	fName = param.name();
	fBound.first  = param.min();
	fBound.second = param.max();

	fSumT = MGROOT::Hist::New(StrFmt("SUMT%06d_%s", RooResult::COUNT, param.samp()->name().c_str()), (*param.samp())(), true);
	fSamp = MGROOT::Hist::New(StrFmt("SAMP%06d_%s", RooResult::COUNT, param.samp()->name().c_str()), (*param.samp())());
	TH1 * hSamp = (*fSamp)();
	TH1 * hSumT = (*fSumT)();

	const Double_t LMTVAL = 1.0e-8;
	std::vector<TH1 *> hTemp(param.num(), nullptr);
	for (Int_t it = 0; it < param.num(); ++it) {
		Hist * temp = MGROOT::Hist::New(StrFmt("TEMP%06d_%s", RooResult::COUNT, param.temp(it)->name().c_str()), (*param.temp(it))());
		hTemp.at(it) = (*temp)();

		Double_t sum = hTemp.at(it)->Integral();
		Bool_t opt = (sum < 100.);
		for (Int_t ib = 0; ib <= hTemp.at(it)->GetNbinsX()+1; ++ib) {
			if (fluc) {
				Double_t men = hTemp.at(it)->GetBinContent(ib);
				Double_t prob = men / sum;
				Double_t sgm = std::sqrt(sum * prob * (1. - prob));
				const Int_t maxNTry = 10; Int_t itry = 1;
				Double_t val = (opt) ? MGROOT::Random::Binomial(Int_t(sum+0.5), prob) : MGROOT::Random::Gaus(men, sgm);
				while (MGNumc::Compare(val, 0.) < 0 && itry <= maxNTry) {
					val = (opt) ? MGROOT::Random::Binomial(Int_t(sum+0.5), prob) : MGROOT::Random::Gaus(men, sgm);
					itry++;
				}
				if (itry > maxNTry) continue;
				hTemp.at(it)->SetBinContent(ib, val);
			}
			if(MGNumc::Compare(hTemp.at(it)->GetBinContent(ib), LMTVAL) <= 0)
  		   hTemp.at(it)->SetBinContent(ib, LMTVAL);
		}
 
		temp->normalized(MGROOT::Hist::kIntegral);
		fTemp.push_back(temp);
	}
	fParam.resize(param.num(), std::make_pair(0., 0.));

	// MsgLevel  DEBUG=0, INFO=1, PROGRESS=2, WARNING=3, ERROR=4, FATAL=5
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	
	//---- Binning ----//
	Double_t lwBin = hSamp->FindBin(min());
	Double_t upBin = hSamp->FindBin(max());
	Double_t total = hSamp->Integral(lwBin, upBin);
	Double_t lwLmt = -7.0 * std::sqrt(total);
	Double_t upLmt = total + 7.0 * std::sqrt(total);
	Double_t meanVal = (total / param.num());
	
	//---- RooFit ----//
	RooRealVar xvar(fName.c_str(), "", min(), max());
  std::vector<RooRealVar *>  tempVar(param.num(), nullptr);
  RooArgList tempVarList;
	for (Int_t it = 0; it < param.num(); ++it) {
    tempVar.at(it) = new RooRealVar(CStrFmt("ROOVAR__TEMP%03d", it), "", meanVal, lwLmt, upLmt);
    tempVarList.add(*tempVar.at(it));
  }
	
	//---- Templates and Model ----//
	std::vector<RooDataHist *> tempHist(param.num(), nullptr);
  std::vector<RooHistPdf *>  tempPdf(param.num(), nullptr);
  RooArgList tempPdfList;
	for (Int_t it = 0; it < param.num(); ++it) {
		tempHist.at(it) = new RooDataHist(CStrFmt("ROOHIST__TEMP%03d", it), "", xvar, hTemp.at(it));
    tempPdf.at(it) = new RooHistPdf(CStrFmt("ROOPDF__TEMP%03d", it), "", xvar, *tempHist.at(it));
    tempPdfList.add(*tempPdf.at(it));
	}
	RooAddPdf model("model", "model", tempPdfList, tempVarList);

	//---- Sample ----//
	RooDataHist sampHist("ROOHIST__SAMP", "", xvar, hSamp);
  
	//---- Fit ----//
  model.fitTo(sampHist, RooFit::Extended(extended), RooFit::PrintLevel(-1));

  for (Int_t it = 0; it < param.num(); ++it) {
		for (Int_t ib = 0; ib <= hTemp.at(it)->GetNbinsX()+1; ++ib)
			if(MGNumc::Compare(hTemp.at(it)->GetBinContent(ib), LMTVAL) <= 0)
  		   hTemp.at(it)->SetBinContent(ib, 0.0);
		
		fParam.at(it).first  = tempVar.at(it)->getVal();
		fParam.at(it).second = tempVar.at(it)->getError();
		Double_t cntNum = hTemp.at(it)->Integral(lwBin, upBin);
		hTemp.at(it)->Scale(fParam.at(it).first / cntNum);
		hSumT->Add(hTemp.at(it));
  }
	
	fNdf = -(param.num() - (!extended));
	for (Int_t ibin = lwBin; ibin <= upBin; ++ibin) {
		Double_t sampVal = hSamp->GetBinContent(ibin);
		Double_t sampErr = hSamp->GetBinError(ibin);
		if (MGNumc::EqualToZero(sampVal)) continue;
		
		Double_t tempVal = 0.;
		Double_t tempErr = 0.;
		for (Int_t it = 0; it < param.num(); ++it) {
			tempVal += hTemp.at(it)->GetBinContent(ibin);
			tempErr += hTemp.at(it)->GetBinError(ibin) * hTemp.at(it)->GetBinError(ibin);
		}
		tempErr = std::sqrt(tempErr);
		
		Double_t totErr = std::sqrt(sampErr * sampErr + tempErr * tempErr);
		Double_t resVal = (sampVal - tempVal) / totErr;
		
		fChisq += (resVal * resVal);
		fNdf++;
	}
	fValid = true;

	for (Int_t it = 0; it < param.num(); ++it) {
    delete tempHist.at(it); delete tempPdf.at(it); delete tempVar.at(it);
		tempHist.at(it) = nullptr; tempPdf.at(it) = nullptr; tempVar.at(it) = nullptr;
	}
	
	RooMsgService::instance().reset();
}
		
MGROOT::Fit::RooResult::~RooResult() {
	fValid = false;
	fName = "";
	fBound.first  = 0.0;
	fBound.second = 0.0;
	if (fSamp) { fSamp->clear(); fSamp = nullptr; };
	if (fSumT) { fSumT->clear(); fSumT = nullptr; };
	for (Int_t it = 0; it < fTemp.size(); ++it)
		if (fTemp.at(it)) { fTemp.at(it)->clear(); fTemp.at(it) = nullptr; };
	fTemp.clear();
	fParam.clear();
	fChisq = 0.0;
	fNdf = 0;
}

MGROOT::Fit::RooSysResult::RooSysResult(MGROOT::Fit::RooParam & param, Bool_t extended, Int_t ntimes) : fRooResult(param, extended), fSysValid(false), fSysChisq(0.), fSysNdf(0) {
	if (!param.valid()) return;
	if (!fRooResult.valid()) return;
	if (ntimes <= 0) return;
	const Double_t CIlevel95 = 0.95;
	const Double_t CIlevel68 = 0.68;
	const Double_t CIlevel = CIlevel95;
	const Double_t CIalpha = 0.5 * (1.0 - CIlevel);
	Long64_t nreqs = Long64_t(ntimes / CIlevel + 10);

	std::vector<MGROOT::Fit::RooSysResult::PARAM> perparvec(nreqs);
	Bool_t FlucOpt = true;

	Long64_t ndf = 0;
	Long64_t cntreq = 0;
	Long64_t iter = 0, niter = 10 * nreqs;
	while (iter < niter && cntreq < nreqs) {
		iter++;
		MGROOT::Fit::RooResult rlt(param, extended, FlucOpt);
		if (!rlt.valid()) continue;
		perparvec.at(cntreq) = MGROOT::Fit::RooSysResult::PARAM(rlt);
		if (ndf==0) ndf = rlt.ndf();
		cntreq++;
	}
	perparvec.resize(cntreq);
	if (perparvec.size() == 0) return;
	std::sort(perparvec.begin(), perparvec.end(), RooSysResult::PARAM_sort());
	Long64_t startsize = Long64_t(perparvec.size() * (      CIalpha)); 
	Long64_t finalsize = Long64_t(perparvec.size() * (1.0 - CIalpha)); 
	std::vector<MGROOT::Fit::RooSysResult::PARAM> parvec(perparvec.begin()+startsize, perparvec.begin()+finalsize);
	
	//---- Binning ---//
	Double_t lwBin = (*param.samp())()->FindBin(param.min());
	Double_t upBin = (*param.samp())()->FindBin(param.max());
	Double_t total = (*param.samp())()->Integral(lwBin, upBin);

	Double_t avgchi = 0.;
	Double_t sumpar = 0.;
	std::vector<Double_t> sumwgt(param.num(), 0.0);
	std::vector<std::pair<Double_t, Double_t>> avgpar(param.num(), std::make_pair(0.0, 0.0));
	for (Int_t it = 0; it < parvec.size(); ++it) {
		MGROOT::Fit::RooSysResult::PARAM & par = parvec.at(it);
		for (Int_t ip = 0; ip < param.num(); ++ip) {
			Double_t err2 = (par.err.at(ip) * par.err.at(ip));
			sumwgt.at(ip) += (1. / err2);
			avgpar.at(ip).first  += par.val.at(ip) / err2;
			avgpar.at(ip).second += err2;
		}
		avgchi += par.chi;
	}
	avgchi /= parvec.size();
	
	for (Int_t ip = 0; ip < param.num(); ++ip) {
		avgpar.at(ip).first  = avgpar.at(ip).first / sumwgt.at(ip);
		avgpar.at(ip).second = std::sqrt(avgpar.at(ip).second / parvec.size());
		sumpar += avgpar.at(ip).first;
	}

	std::vector<Double_t> avgerr(param.num(), 0.0);
	for (Int_t it = 0; it < parvec.size(); ++it) {
		MGROOT::Fit::RooSysResult::PARAM & par = parvec.at(it);
		for (Int_t ip = 0; ip < param.num(); ++ip) {
			Double_t dif = (par.val.at(ip) - avgpar.at(ip).first) / par.err.at(ip);
			avgerr.at(ip) += (dif * dif);
		}
	}
	for (Int_t ip = 0; ip < param.num(); ++ip) {
		avgerr.at(ip) = std::sqrt(avgerr.at(ip) / sumwgt.at(ip));
	}

	//---- Parameters ----//
	fSysFitSet = perparvec;
	fSysErr.resize(param.num(), 0.0);
	fSysParam.resize(param.num(), std::make_pair(0.0, 0.0));
	for (Int_t ip = 0; ip < param.num(); ++ip) {
		fSysParam.at(ip).first  = avgpar.at(ip).first;
		fSysParam.at(ip).second = avgpar.at(ip).second;
		fSysErr.at(ip) = avgerr.at(ip);
	}
	fSysChisq = avgchi * ndf;
	fSysNdf   = ndf;

	fSysValid = true;
}

MGROOT::Fit::RooSysResult::~RooSysResult() {
	fSysFitSet.clear();
	fSysValid = false;
	fSysParam.clear();
	fSysErr.clear();
	fSysChisq = 0.0;
	fSysNdf = 0;
}

*/

#endif // __MGFit_C__
