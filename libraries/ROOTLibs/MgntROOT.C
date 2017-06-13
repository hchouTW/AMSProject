#ifndef __MgntROOT_C__
#define __MgntROOT_C__
#include "MgntROOT.h"


/********************/
/****  MgntROOT  ****/
/********************/

/****************/
/****  Math  ****/
/****************/
#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"
#include "Math/GSLMinimizer1D.h"
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

namespace MgntROOT {
	namespace Math {
		// One-Dimensional Minimization
		using ROOT::Math::Functor1D;
		using ROOT::Math::BrentMinimizer1D;
		using ROOT::Math::GSLMinimizer1D;

		// Multi-Dimensional Minimization
		using ROOT::Math::Functor;
		using ROOT::Math::Minimizer;
		using ROOT::Math::GSLMinimizer;
		using ROOT::Math::GSLSimAnMinimizer;
		using ROOT::Minuit2::Minuit2Minimizer;
	}
}


/*****************/
/****  MtxLB  ****/
/*****************/
#include "TArray.h"
#include "TVector.h"
#include "TMatrix.h"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
namespace MgntROOT {
	namespace MtxLB {
		// TArray 
  	using TArrI = TArrayI;
  	using TArrF = TArrayF;
  	using TArrD = TArrayD;
  	
		// TVector TMatrix Package
		using TVecI = TVectorT<Int_t>;
  	using TVecF = TVectorT<Float_t>;
  	using TVecD = TVectorT<Double_t>;

  	using TMtxI = TMatrixT<Int_t>;
  	using TMtxF = TMatrixT<Float_t>;
  	using TMtxD = TMatrixT<Double_t>;
		
		using TMtxSymI = TMatrixTSym<Int_t>;
		using TMtxSymF = TMatrixTSym<Float_t>;
  	using TMtxSymD = TMatrixTSym<Double_t>;
		
		using TMtxSparseI = TMatrixTSparse<Int_t>;
		using TMtxSparseF = TMatrixTSparse<Float_t>;
  	using TMtxSparseD = TMatrixTSparse<Double_t>;

		// SVector SMatrix Package
  	template <unsigned int D>
  	using SVecI = ROOT::Math::SVector<Int_t, D>;

  	template <unsigned int D>
  	using SVecF = ROOT::Math::SVector<Float_t, D>;

		template <unsigned int D>
  	using SVecD = ROOT::Math::SVector<Double_t, D>;

  	template <unsigned int D1, unsigned int D2 = D1>
  	using SMtxI = ROOT::Math::SMatrix<Int_t, D1, D2, ROOT::Math::MatRepStd<Int_t, D1, D2>>;
  	
		template <unsigned int D1, unsigned int D2 = D1>
  	using SMtxF = ROOT::Math::SMatrix<Float_t, D1, D2, ROOT::Math::MatRepStd<Float_t, D1, D2>>;

  	template <unsigned int D1, unsigned int D2 = D1>
  	using SMtxD = ROOT::Math::SMatrix<Double_t, D1, D2, ROOT::Math::MatRepStd<Double_t, D1, D2>>;

		template <unsigned int D>
  	using SMtxSymI = ROOT::Math::SMatrix<Int_t, D, D, ROOT::Math::MatRepSym<Int_t, D>>;
  	
		template <unsigned int D>
  	using SMtxSymF = ROOT::Math::SMatrix<Float_t, D, D, ROOT::Math::MatRepSym<Float_t, D>>;
  
		template <unsigned int D>
  	using SMtxSymD = ROOT::Math::SMatrix<Double_t, D, D, ROOT::Math::MatRepSym<Double_t, D>>;
  
		using SIdMtx = ROOT::Math::SMatrixIdentity;

		// SVector Template Function 
		using ROOT::Math::Dot;
		using ROOT::Math::Mag2;
		using ROOT::Math::Mag;
		using ROOT::Math::Lmag2;
		using ROOT::Math::Lmag;
		using ROOT::Math::Cross;
		using ROOT::Math::Unit;
		using ROOT::Math::TensorProd;
		using ROOT::Math::fabs;
		using ROOT::Math::sqr;
		using ROOT::Math::sqrt;
		
		// SMatrix Template Function 
		using ROOT::Math::Transpose;
		using ROOT::Math::Similarity;
		using ROOT::Math::SimilarityT;
		using ROOT::Math::fabs;
		using ROOT::Math::sqr;
		using ROOT::Math::sqrt;
	}
}

/***************/
/****  Fit  ****/
/***************/
//#define ROOFIT
#ifdef ROOFIT
#include "TFractionFitter.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooBifurGauss.h"
#include "RooMinuit.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"
namespace MgntROOT {
	namespace Fit {
		class RooParam;
		class RooResult;
		class RooSysResult;
	}
}

class MgntROOT::Fit::RooParam {
	public :
		using LIST = std::vector<std::string>;
		using VLIST = std::vector<Hist *>;

	public :
		RooParam(const std::string& name, MgntROOT::Hist * samp, const RooParam::LIST& list, Double_t min = 0., Double_t max = 0.); 
		RooParam(const std::string& name, MgntROOT::Hist * samp, const RooParam::VLIST& vlist, Double_t min = 0., Double_t max = 0.); 
		~RooParam() {}

		Bool_t      valid() { return (fName != "" && fSamp != nullptr && fTemp.size()); }
		std::string name() { return fName; }
		Hist *      samp() { return fSamp; }
		Hist *      temp(Int_t idx) { return ((idx < 0 || idx >= fTemp.size()) ? nullptr : fTemp.at(idx)); }
		Int_t       num() { return fTemp.size(); }
		Double_t    min() { return fBound.first; }
		Double_t    max() { return fBound.second; }

	protected :
		std::string                   fName;
		MgntROOT::Hist *              fSamp;
		std::vector<MgntROOT::Hist *> fTemp;
		std::pair<Double_t, Double_t> fBound;
};

class MgntROOT::Fit::RooResult {
	public :
		RooResult(MgntROOT::Fit::RooParam & param, Bool_t extended = true, Bool_t fluc = false);
		~RooResult();
		
		Bool_t      valid() { return fValid; }
		std::string name() { return fName; }
		Double_t    min()  { return fBound.first; }
		Double_t    max()  { return fBound.second; }
		Hist *      samp() { return fSamp; }
		Hist *      sumt() { return fSumT; }
		Hist *      temp(Int_t idx) { return ((idx < 0 || idx >= fTemp.size()) ? nullptr : fTemp.at(idx)); }
		Double_t    val(Int_t idx) { return ((idx < 0 || idx >= fParam.size()) ? 0.0 : fParam.at(idx).first); }
		Double_t    err(Int_t idx) { return ((idx < 0 || idx >= fParam.size()) ? 0.0 : fParam.at(idx).second); }
		Int_t       num() { return fTemp.size(); }
		Double_t    chisq() { return fChisq; }
		Int_t       ndf() { return fNdf; }

	protected :
		Bool_t                                     fValid;
		std::string                                fName;
		std::pair<Double_t, Double_t>              fBound;
		MgntROOT::Hist *                           fSamp;
		MgntROOT::Hist *                           fSumT;
		std::vector<MgntROOT::Hist *>              fTemp;
		std::vector<std::pair<Double_t, Double_t>> fParam;
		Double_t                                   fChisq;
		Int_t                                      fNdf;

		static Long64_t COUNT;
};
	
Long64_t MgntROOT::Fit::RooResult::COUNT = 0;

class MgntROOT::Fit::RooSysResult {
	public :
		class PARAM {
			public :
				PARAM() : chi(0.0) {}
				PARAM(MgntROOT::Fit::RooResult& rlt)
					{
						val.resize(rlt.num(), 0.0); 
						err.resize(rlt.num(), 0.0);
						for (Int_t it = 0; it < rlt.num(); ++it) {
							val.at(it) = rlt.val(it);
							err.at(it) = rlt.err(it);
						}
						chi = (rlt.num()!=0) ? rlt.chisq() / rlt.ndf() : 0.0;
					}
				~PARAM() {}
			
				std::vector<Double_t> val;
				std::vector<Double_t> err;
				Double_t              chi;
		};

		struct PARAM_sort {
			Bool_t operator() (const RooSysResult::PARAM& param1, const RooSysResult::PARAM& param2) {
				if (param1.chi < param2.chi) return true;
				else return false;
			}
		};

	public :
		RooSysResult(MgntROOT::Fit::RooParam & param, Bool_t extended = true, Int_t ntimes = 400);
		~RooSysResult();

		MgntROOT::Fit::RooResult& operator()() { return fRooResult; }

		Bool_t      valid() { return fRooResult.valid(); }
		std::string name()  { return fRooResult.name(); }
		Double_t    min()   { return fRooResult.min(); }
		Double_t    max()   { return fRooResult.max(); }
		Int_t       num()   { return fRooResult.num(); }

		Hist *      samp()          { return fRooResult.samp(); }
		Hist *      sumt()          { return fRooResult.sumt(); }
		Hist *      temp(Int_t idx) { return ((idx < 0 || idx >= fRooResult.num()) ? nullptr : fRooResult.temp(idx)); }
		Double_t    val(Int_t idx)  { return ((idx < 0 || idx >= fRooResult.num()) ? 0.0 : fRooResult.val(idx)); }
		Double_t    err(Int_t idx)  { return ((idx < 0 || idx >= fRooResult.num()) ? 0.0 : fRooResult.err(idx)); }
		Double_t    chisq()         { return fRooResult.chisq(); }
		Int_t       ndf()           { return fRooResult.ndf(); }
		
		Bool_t      sysValid()        { return fSysValid; }
		Double_t    sysVal(Int_t idx) { return ((idx < 0 || idx >= fSysParam.size()) ? 0.0 : fSysParam.at(idx).first); }
		Double_t    satErr(Int_t idx) { return ((idx < 0 || idx >= fSysParam.size()) ? 0.0 : fSysParam.at(idx).second); }
		Double_t    sysErr(Int_t idx) {  return ((idx < 0 || idx >= fSysErr.size()) ? 0.0 : fSysErr.at(idx)); }
		Double_t    sysChisq()        { return fSysChisq; }
		Int_t       sysNdf()          { return fSysNdf; }

		Int_t                 nTrails() { return fSysFitSet.size(); }
		RooSysResult::PARAM * trail(Int_t idx) { return (idx < 0 || idx > fSysFitSet.size()) ? nullptr : &fSysFitSet.at(idx); } 

	protected :
		MgntROOT::Fit::RooResult fRooResult;

		Bool_t                                     fSysValid;
		std::vector<std::pair<Double_t, Double_t>> fSysParam;
		std::vector<Double_t>                      fSysErr;
		Double_t                                   fSysChisq;
		Int_t                                      fSysNdf;
		
		std::vector<RooSysResult::PARAM> fSysFitSet;
};
		

MgntROOT::Fit::RooParam::RooParam(const std::string& name, MgntROOT::Hist * samp, const RooParam::LIST& list, Double_t min, Double_t max) : fName(""), fSamp(nullptr) {
	fBound.first  = 0.0;
	fBound.second = 0.0;
	if (name == "") return;
	if (samp == nullptr) return;
	for (const std::string& elem : list)
		if (Hist::Head(elem) == nullptr) return;

	Bool_t norm = (MgntNum::Compare(min, max) < 0);
  Double_t minEdge = (*samp)()->GetBinLowEdge(1);
  Double_t maxEdge = (*samp)()->GetBinLowEdge((*samp)()->GetNbinsX()+1);
  fBound.first  = (min < minEdge || !norm) ? minEdge : min;
  fBound.second = (max > maxEdge || !norm) ? maxEdge : max;

	fName = name;
	fSamp = samp;

	for (const std::string& elem : list)
		fTemp.push_back(Hist::Head(elem));
}

MgntROOT::Fit::RooParam::RooParam(const std::string& name, MgntROOT::Hist * samp, const RooParam::VLIST& vlist, Double_t min, Double_t max) : fName(""), fSamp(nullptr) {
	fBound.first  = 0.0;
	fBound.second = 0.0;
	if (name == "") return;
	if (samp == nullptr) return;
	for (MgntROOT::Hist * temp : vlist)
		if (temp == nullptr) return;

	Bool_t norm = (MgntNum::Compare(min, max) < 0);
  Double_t minEdge = (*samp)()->GetBinLowEdge(1);
  Double_t maxEdge = (*samp)()->GetBinLowEdge((*samp)()->GetNbinsX()+1);
  fBound.first  = (min < minEdge || !norm) ? minEdge : min;
  fBound.second = (max > maxEdge || !norm) ? maxEdge : max;

	fName = name;
	fSamp = samp;
	for (MgntROOT::Hist * temp : vlist)
		fTemp.push_back(temp);
}

MgntROOT::Fit::RooResult::RooResult(MgntROOT::Fit::RooParam& param, Bool_t extended, Bool_t fluc) : fValid(false), fName(""), fSamp(nullptr), fSumT(nullptr), fChisq(0.), fNdf(0) {
	RooResult::COUNT++;
	fBound.first  = 0.0;
	fBound.second = 0.0;
	if (!param.valid()) return;
	
	//---- Histogram ----//
	fName = param.name();
	fBound.first  = param.min();
	fBound.second = param.max();

	fSumT = MgntROOT::Hist::New(StrFmt("SUMT%06d_%s", RooResult::COUNT, param.samp()->name().c_str()), (*param.samp())(), true);
	fSamp = MgntROOT::Hist::New(StrFmt("SAMP%06d_%s", RooResult::COUNT, param.samp()->name().c_str()), (*param.samp())());
	TH1 * hSamp = (*fSamp)();
	TH1 * hSumT = (*fSumT)();

	const Double_t LMTVAL = 1.0e-8;
	std::vector<TH1 *> hTemp(param.num(), nullptr);
	for (Int_t it = 0; it < param.num(); ++it) {
		Hist * temp = MgntROOT::Hist::New(StrFmt("TEMP%06d_%s", RooResult::COUNT, param.temp(it)->name().c_str()), (*param.temp(it))());
		hTemp.at(it) = (*temp)();

		Double_t sum = hTemp.at(it)->Integral();
		Bool_t opt = (sum < 100.);
		for (Int_t ib = 0; ib <= hTemp.at(it)->GetNbinsX()+1; ++ib) {
			if (fluc) {
				Double_t men = hTemp.at(it)->GetBinContent(ib);
				Double_t prob = men / sum;
				Double_t sgm = std::sqrt(sum * prob * (1. - prob));
				const Int_t maxNTry = 10; Int_t itry = 1;
				Double_t val = (opt) ? MgntROOT::Random::Binomial(Int_t(sum+0.5), prob) : MgntROOT::Random::Gaus(men, sgm);
				while (MgntNum::Compare(val, 0.) < 0 && itry <= maxNTry) {
					val = (opt) ? MgntROOT::Random::Binomial(Int_t(sum+0.5), prob) : MgntROOT::Random::Gaus(men, sgm);
					itry++;
				}
				if (itry > maxNTry) continue;
				hTemp.at(it)->SetBinContent(ib, val);
			}
			if(MgntNum::Compare(hTemp.at(it)->GetBinContent(ib), LMTVAL) <= 0)
  		   hTemp.at(it)->SetBinContent(ib, LMTVAL);
		}
 
		temp->normalized(MgntROOT::Hist::kIntegral);
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
			if(MgntNum::Compare(hTemp.at(it)->GetBinContent(ib), LMTVAL) <= 0)
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
		if (MgntNum::EqualToZero(sampVal)) continue;
		
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
		
MgntROOT::Fit::RooResult::~RooResult() {
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

MgntROOT::Fit::RooSysResult::RooSysResult(MgntROOT::Fit::RooParam & param, Bool_t extended, Int_t ntimes) : fRooResult(param, extended), fSysValid(false), fSysChisq(0.), fSysNdf(0) {
	if (!param.valid()) return;
	if (!fRooResult.valid()) return;
	if (ntimes <= 0) return;
	const Double_t CIlevel95 = 0.95;
	const Double_t CIlevel68 = 0.68;
	const Double_t CIlevel = CIlevel95;
	const Double_t CIalpha = 0.5 * (1.0 - CIlevel);
	Long64_t nreqs = Long64_t(ntimes / CIlevel + 10);

	std::vector<MgntROOT::Fit::RooSysResult::PARAM> perparvec(nreqs);
	Bool_t FlucOpt = true;

	Long64_t ndf = 0;
	Long64_t cntreq = 0;
	Long64_t iter = 0, niter = 10 * nreqs;
	while (iter < niter && cntreq < nreqs) {
		iter++;
		MgntROOT::Fit::RooResult rlt(param, extended, FlucOpt);
		if (!rlt.valid()) continue;
		perparvec.at(cntreq) = MgntROOT::Fit::RooSysResult::PARAM(rlt);
		if (ndf==0) ndf = rlt.ndf();
		cntreq++;
	}
	perparvec.resize(cntreq);
	if (perparvec.size() == 0) return;
	std::sort(perparvec.begin(), perparvec.end(), RooSysResult::PARAM_sort());
	Long64_t startsize = Long64_t(perparvec.size() * (      CIalpha)); 
	Long64_t finalsize = Long64_t(perparvec.size() * (1.0 - CIalpha)); 
	std::vector<MgntROOT::Fit::RooSysResult::PARAM> parvec(perparvec.begin()+startsize, perparvec.begin()+finalsize);
	
	//---- Binning ---//
	Double_t lwBin = (*param.samp())()->FindBin(param.min());
	Double_t upBin = (*param.samp())()->FindBin(param.max());
	Double_t total = (*param.samp())()->Integral(lwBin, upBin);

	Double_t avgchi = 0.;
	Double_t sumpar = 0.;
	std::vector<Double_t> sumwgt(param.num(), 0.0);
	std::vector<std::pair<Double_t, Double_t>> avgpar(param.num(), std::make_pair(0.0, 0.0));
	for (Int_t it = 0; it < parvec.size(); ++it) {
		MgntROOT::Fit::RooSysResult::PARAM & par = parvec.at(it);
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
		MgntROOT::Fit::RooSysResult::PARAM & par = parvec.at(it);
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

MgntROOT::Fit::RooSysResult::~RooSysResult() {
	fSysFitSet.clear();
	fSysValid = false;
	fSysParam.clear();
	fSysErr.clear();
	fSysChisq = 0.0;
	fSysNdf = 0;
}
#endif


/*****************/
/****  Color  ****/
/*****************/
Color_t MgntROOT::Color::At(UInt_t idx, Int_t nset) {
	Color_t color;
	if      (nset <  0) color = idx;
	else if (nset == 0) color = TColor::GetColorPalette(idx);
	else {
		UInt_t step = TColor::GetNumberOfColors() / nset;
		if (step == 0)	color = TColor::GetColorPalette(idx);
		else            color = TColor::GetColorPalette(idx * step);
	}
	return color;
}


/****************/
/****  Line  ****/
/****************/


/******************/
/****  Marker  ****/
/******************/
Style_t MgntROOT::Marker::Style(UInt_t idx) {
	const UInt_t NumStyle = 14;
	Style_t markerStyle = 20;
	Style_t style = (idx%NumStyle);
	switch (style) {
		// Circle (Full, Open)
		case  0 : markerStyle = 20; break;
		case  1 : markerStyle = 24; break;
		// Square (Full, Open)
		case  2 : markerStyle = 21; break;
		case  3 : markerStyle = 25; break;
		// TriangleUp (Full, Open)
		case  4 : markerStyle = 22; break;
		case  5 : markerStyle = 26; break;
		// TriangleDown (Full, Open)
		case  6 : markerStyle = 23; break;
		case  7 : markerStyle = 32; break;
		// Diamond (Full, Open)
		case  8 : markerStyle = 33; break;
		case  9 : markerStyle = 27; break;
		// Cross (Full, Open)
		case 10 : markerStyle = 34; break;
		case 11 : markerStyle = 28; break;
		// Star (Full, Open)
		case 12 : markerStyle = 29; break;
		case 13 : markerStyle = 30; break;
		default : break;
	}
	return markerStyle;
}


/****************/
/****  Text  ****/
/****************/
void MgntROOT::Text::Draw(const Text::Txt_t& txt, const Text::Att_t& att) {
	TLatex ltx;
	ltx.SetNDC();
	ltx.SetTextColor(txt.color);
	ltx.SetTextFont(txt.font);
	ltx.SetTextSize(txt.size);
	ltx.SetTextAlign(att.align);
	ltx.SetTextAngle(att.angle);
	ltx.DrawLatex(att.x, att.y, txt.text.c_str());
}


/*****************/
/****  Style  ****/
/*****************/
void MgntROOT::Style::LoadDefaultEnvironment() {
	std::cout << "MgntROOT::Style : Load default environment.\n";
  Int_t col = 10;
  fStyle.SetFrameBorderMode(col);
  fStyle.SetFrameFillColor(col);
  fStyle.SetCanvasBorderMode(col);
  fStyle.SetCanvasColor(col);
  fStyle.SetPadBorderMode(col);
  fStyle.SetPadColor(col);
  fStyle.SetStatColor(col);

  // set the paper & margin sizes
  fStyle.SetPaperSize(20, 26);
	    
  // set margin sizes
  fStyle.SetPadTopMargin(0.12);
  fStyle.SetPadBottomMargin(0.12);
  fStyle.SetPadRightMargin(0.15);
  fStyle.SetPadLeftMargin(0.15);
    
  // use large fonts
  Int_t font = 42; // Helvetica

 	fStyle.SetTextFont(font);
  fStyle.SetTextSize(0.03);
	
	fStyle.SetAxisColor(1, "xyz");
	fStyle.SetStripDecimals(false);	

	fStyle.SetLabelColor(1, "xyz");
  fStyle.SetLabelSize(0.03, "xyz");
	fStyle.SetLabelFont(font, "xyz");
	fStyle.SetLabelOffset(0.005, "xyz");

 	fStyle.SetTitleColor(1, "xyz"); 
  fStyle.SetTitleSize(0.04, "xyz");
	fStyle.SetTitleFont(font, "xyz");
	fStyle.SetTitleOffset(1.05, "x");
	fStyle.SetTitleOffset(1.20, "yz");

	// legend
	fStyle.SetLegendBorderSize(0);
	fStyle.SetLegendFillColor(0);
	fStyle.SetLegendTextSize(0.025);
	fStyle.SetLegendFont(42);

  // use bold lines and markers
	fStyle.SetMarkerColor(1);
  fStyle.SetMarkerStyle(20);
  fStyle.SetMarkerSize(1.2);
  fStyle.SetHistLineWidth((Width_t) 2.);
  fStyle.SetLineStyleString(2, "[12 12]"); // postscript dashes
	
	// set the hist title
	fStyle.SetTitleAlign(22);
	fStyle.SetTitleBorderSize(0);
	fStyle.SetTitleFillColor(10);
	fStyle.SetTitleFontSize(0.05);
	fStyle.SetTitleX(0.5);  
	fStyle.SetTitleY(0.95);
	
  // do not display any of the standard histogram decorations
  fStyle.SetOptTitle(0);
  fStyle.SetOptStat(0);
  fStyle.SetOptFit(0);

	// set histogram x-axis error to zero
	fStyle.SetErrorX(0);
  
	// put tick marks on top and RHS of plots
	fStyle.SetPadGridX(0);
	fStyle.SetPadGridY(0);
  fStyle.SetPadTickX(1);
  fStyle.SetPadTickY(1);

	// set axis scale
	fStyle.SetOptLogx(0);
	fStyle.SetOptLogy(0);
	fStyle.SetOptLogz(1);

	// set color palette
	Style::ColorPalette(-1);
	fStyle.SetNumberContours(254);

	// set windows
	Style::WindowsScale(1024, 768);

	// set event status
	if (fStyle.GetShowEventStatus() == 0)
		fStyle.ToggleEventStatus();

	// update
	fStyle.cd();
}

void MgntROOT::Style::ColorPalette(Int_t pattern) {
	if (pattern >= 0) {
		fStyle.SetPalette(pattern, 0, 1);
	}
	else {
		const Int_t Num = 254;
		const Int_t NSet = 6;
		Int_t    palette[Num] = {0};
		Double_t red[NSet]    = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
		Double_t green[NSet]  = { 0.0, 1.0, 1.0, 1.0, 0.0, 0.0 };
		Double_t blue[NSet]   = { 1.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
		Double_t length[NSet] = { 0.0, .20, .40, .60, .80, 1.0 };
		Int_t FI = TColor::CreateGradientColorTable(NSet, length, red, green, blue, Num);
		for (Int_t i = 0; i < Num; ++i) palette[i] = FI + i;
		fStyle.SetPalette(Num, palette, 1);
	}
	fStyle.cd();
}

void MgntROOT::Style::WindowsScale(UInt_t w, UInt_t h) {
	fStyle.SetCanvasDefW(w+20);
	fStyle.SetCanvasDefH(h+40);
	fStyle.cd();
}

void MgntROOT::Style::AxisScale(Bool_t xaxis, Bool_t yaxis, Bool_t zaxis) {
	fStyle.SetOptLogx(xaxis);
	fStyle.SetOptLogy(yaxis);
	fStyle.SetOptLogz(zaxis);
	fStyle.cd();
}


/******************/
/****  Canvas  ****/
/******************/
TVirtualPad * MgntROOT::Canvas::cd(UInt_t idx, const Canvas::AxisScl_t& scl) {
	TVirtualPad * pad = fCanvas.cd(idx);
	if (pad != nullptr) {
		pad->SetLogx(scl.x);
		pad->SetLogy(scl.y);
		pad->SetLogz(scl.z);
    pad->Modified();
    pad->Update();
	}
	return pad;
}

void MgntROOT::Canvas::setTemp(Canvas::Temp temp) {
	switch(temp) {
		case Canvas::kNone :
			setTemp(1, 1);
			break;
		case Canvas::kStat :
			{
				setTemp(2, 1);
				if (checkPad(1)) {
					TVirtualPad * pad = setPad(1, Canvas::Pad_t(0.0, 1.0, 0.4, 1.0));
					pad->SetTopMargin(0.15);
					pad->SetBottomMargin(0.00);
				}
				if (checkPad(2)) {
					TVirtualPad * pad = setPad(2, Canvas::Pad_t(0.0, 1.0, 0.0, 0.4));
					pad->SetTopMargin(0.00);
					pad->SetBottomMargin(0.20);
				}
			}
			break;
		default : break;
	}
	fTemp = temp;
	fCanvas.cd(0);
	update();
}
		
void MgntROOT::Canvas::setTemp(Int_t nx, Int_t ny, Float_t wmargin, Float_t hmargin, Int_t color) {
	Int_t divx = (nx < 1 || ny < 1) ? 1 : nx;
	Int_t divy = (nx < 1 || ny < 1) ? 1 : ny;
	fCanvas.cd(0);
	fCanvas.Divide(divx, divy, wmargin, hmargin, color);
	Int_t npad = divx * divy;
	Double_t dx = 1. / divx;
	Double_t dy = 1. / divy;
	for (Int_t iy = 0; iy < divy; ++iy) {
		for (Int_t ix = 0; ix < divx; ++ix) {
			Int_t it = ix + iy * divx + 1;
			Double_t xi = dx * (ix);
			Double_t xf = dx * (ix + 1);
			Double_t yi = dy * (ny - iy - 1);
			Double_t yf = dy * (ny - iy);
			if (!checkPad(it)) continue;
			TVirtualPad * pad = setPad(2, Canvas::Pad_t(xi, xf, yi, yf));
		}
	}
	fTemp = Canvas::kNone;
	fCanvas.cd(0);
	update();
}

void MgntROOT::Canvas::setSize(Canvas::Size size) {
	Canvas::Size_t type(1024, 768);
	switch(size) {
		case Canvas::kSliceH :     type = Size_t(4096, 3072); break;
		case Canvas::kSliceM :     type = Size_t(2048, 1536); break;
		case Canvas::kSliceL :     type = Size_t(1024,  768); break;
		case Canvas::kA4Vertical : type = Size_t(4961, 7016); break;
		case Canvas::kA4Horizon :  type = Size_t(7016, 4961); break;
		case Canvas::kMacH :       type = Size_t(7680, 4800); break;
		case Canvas::kMacM :       type = Size_t(3840, 2400); break;
		case Canvas::kMacL :       type = Size_t( 960,  600); break;
		case Canvas::kMovie :      type = Size_t(4200, 1800); break;
		default : break;
	}
	fSize = size;
	fCanvas.cd(0);
	fCanvas.SetCanvasSize(type.w, type.h);
	fCanvas.SetTopMargin(0.12);
	fCanvas.SetBottomMargin(0.12);
	fCanvas.SetRightMargin(0.15);
	fCanvas.SetLeftMargin(0.15);
	update();
}

TVirtualPad * MgntROOT::Canvas::setPad(UInt_t idx, const Pad_t& pad) {
	TVirtualPad * vpad = fCanvas.cd(idx);
	if (vpad == nullptr || idx == 0) return nullptr;
	std::string name = Form("%s_%d", fCanvas.GetName(), idx);
	vpad->SetPad(name.c_str(), name.c_str(), pad.xlw, pad.ylw, pad.xup, pad.yup, pad.color, pad.bordersize, pad.bordermode);
	vpad->SetTopMargin(0.12);
	vpad->SetBottomMargin(0.12);
	vpad->SetRightMargin(0.15);
	vpad->SetLeftMargin(0.15);
	return vpad;
}


/*********************/
/****  PdfEditor  ****/
/*********************/
void MgntROOT::PdfEditor::open(PdfEditor::Size size, const std::string& filename, const std::string& filepath, Option_t * option) {
	if (exist()) return;
	if (!MgntSys::TestFile(filepath, 'd')) return;
	fExist = true;
	fFilePage = 0;
	fFilePath = filepath;
	fFileName = (filename != "") ? filename : "PdfPainter";
	MgntROOT::Canvas::Size cvsSize = MgntROOT::Canvas::kSliceH;
	switch (size) {
		case PdfEditor::kSlice      : cvsSize = MgntROOT::Canvas::kSliceH;     break;
		case PdfEditor::kA4Vertical : cvsSize = MgntROOT::Canvas::kA4Vertical; break;
		case PdfEditor::kA4Horizon  : cvsSize = MgntROOT::Canvas::kA4Horizon;  break;
		case PdfEditor::kMac        : cvsSize = MgntROOT::Canvas::kMacM;       break;
		case PdfEditor::kMovie      : cvsSize = MgntROOT::Canvas::kMovie;      break;
		default : break;
	}
	fCanvas.setNameTitle(fFileName, fFileName);
	fCanvas.create(cvsSize);
	std::string fullPathPDF = Form("%s/%s.pdf", fFilePath.c_str(), fFileName.c_str());	
	std::cout << Form("\n\n**************** PdfEditor::OPEN    %-45s  ****************\n", fullPathPDF.c_str());
	create("START");
	MgntROOT::Text::Draw(
		MgntROOT::Text::Txt_t(Form("<< START >>  %s  ", fFileName.c_str())), 
		MgntROOT::Text::Att_t(0.5, 0.5)
	);
	fCanvas.update();
	fCanvas.save(Form("%s(", fullPathPDF.c_str()));
}

void MgntROOT::PdfEditor::create(const std::string& title, PdfEditor::Temp temp) {
	if (!exist()) return;
	++fFilePage;
	std::string statement = (title == "") ? fFileName : Form("%s  >~@~>  %s", fFileName.c_str(), title.c_str());
	MgntROOT::Canvas::Temp cvsTemp = MgntROOT::Canvas::kNone;
	switch(temp) {
		case PdfEditor::kNone : cvsTemp = MgntROOT::Canvas::kNone; break;
		case PdfEditor::kStat : cvsTemp = MgntROOT::Canvas::kStat; break;
		default : break;
	}
	fCanvas.setName(CStrFmt("%s__PAGE%06d", fFileName.c_str(), fFilePage));
	fCanvas.create(fCanvas.getSize(), cvsTemp);
	std::cout << Form("PdfEditor::NEW_PAGE  ~~@ %05d @~~  << %-57s >>\n", fFilePage, statement.c_str());
}
		
void MgntROOT::PdfEditor::create(const std::string& title, Int_t nx, Int_t ny) {
	if (!exist()) return;
	++fFilePage;
	std::string statement = (title == "") ? fFileName : Form("%s  >~@~>  %s", fFileName.c_str(), title.c_str());
	fCanvas.setName(CStrFmt("%s__PAGE%06d", fFileName.c_str(), fFilePage));
	fCanvas.create(fCanvas.getSize(), nx, ny);
	std::cout << Form("PdfEditor::NEW_PAGE  ~~@ %05d @~~  << %-57s >>\n", fFilePage, statement.c_str());
}
		
inline void MgntROOT::PdfEditor::close() {
	if (!exist()) return;
	std::string fullPathPDF = Form("%s/%s.pdf", fFilePath.c_str(), fFileName.c_str());
	create("END");
	MgntROOT::Text::Draw(
		MgntROOT::Text::Txt_t(Form("<< END >>  %s  ", fFileName.c_str())), 
		MgntROOT::Text::Att_t(0.5, 0.5)
	);
	fCanvas.update();
	fCanvas.save(Form("%s)", fullPathPDF.c_str()));
	std::cout << Form("**************** PdfEditor::CLOSE   %-45s  ****************\n\n", fullPathPDF.c_str());
	fExist = false;
	fFileName = "PdfEditor";
	fFilePath = ".";
	fFilePage = 0;
}


/****************/
/****  Axis  ****/
/****************/
MgntROOT::Axis MgntROOT::Axis::Invert(const std::string& title, MgntROOT::Axis& axis) {
	Bool_t result = axis.exist();
	for (Int_t idx = 0; result && idx < axis.nbin()+1; ++idx)
		if (MgntNum::Compare(axis.bins(idx)) <= 0) result = false;
	if (result) {
		std::vector<Double_t> list(2*axis.nbin()+3);
		list.at(axis.nbin()) = 0.0;
		UInt_t nbin = list.size() - 1;
		for (Int_t it = 0; it <= axis.nbin(); ++it) {
			Double_t inv = 1. / axis.bins(it);
			list.at(it)        = -1. * inv;
			list.at(nbin - it) = inv;
		}
		MgntROOT::Axis axis(title, list);
		return axis;
	}
	else {
		MgntSys::Error("<< Axis >>  Invert Axis failure.");
		MgntROOT::Axis axis;
		return axis;
	}
}

TH1D * MgntROOT::Axis::CreateHist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX) {
	if (!axisX.exist()) return nullptr;
	TH1D * hist = new TH1D(name.c_str(), 
		Form("%s;%s", title.c_str(), axisX.title().c_str()), 
		axisX.nbin(), axisX.bins()
	);
	for (Int_t ibin = 1; ibin <= axisX.nbin(); ++ibin)
		hist->SetBinContent(ibin, axisX.width(ibin));
	return hist;
}

TH2D * MgntROOT::Axis::CreateHist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY) {
	if (!axisX.exist() || !axisY.exist()) return nullptr;
	TH2D * hist = new TH2D(name.c_str(), 
		Form("%s;%s;%s", title.c_str(), axisX.title().c_str(), axisY.title().c_str()), 
		axisX.nbin(), axisX.bins(), axisY.nbin(), axisY.bins()
	);
	for (Int_t ibin = 1; ibin <= axisX.nbin(); ++ibin)
		for (Int_t jbin = 1; jbin <= axisY.nbin(); ++jbin)
		hist->SetBinContent(ibin, jbin, (axisX.width(ibin) * axisY.width(jbin)));
	return hist;
}

TH3D * MgntROOT::Axis::CreateHist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, MgntROOT::Axis& axisZ) {
	if (!axisX.exist() || !axisY.exist() || !axisZ.exist()) return nullptr;
	TH3D * hist = new TH3D(name.c_str(), 
		Form("%s;%s;%s;%s", title.c_str(), axisX.title().c_str(), axisY.title().c_str(), axisZ.title().c_str()), 
		axisX.nbin(), axisX.bins(), axisY.nbin(), axisY.bins(), axisZ.nbin(), axisZ.bins()
	);
	for (Int_t ibin = 1; ibin <= axisX.nbin(); ++ibin)
		for (Int_t jbin = 1; jbin <= axisY.nbin(); ++jbin)
			for (Int_t kbin = 1; kbin <= axisZ.nbin(); ++kbin)
				hist->SetBinContent(ibin, jbin, kbin, (axisX.width(ibin) * axisY.width(jbin) * axisZ.width(kbin)));
	return hist;
}

Double_t MgntROOT::Axis::bins(UInt_t idx) {
	Bool_t fine = (exist() && idx < fList.size());
	if (!fine) MgntSys::Error("<< Axis >>  index is failure.");
	return (fine ? fList.at(idx) : 0.0); 
}

Double_t MgntROOT::Axis::width(Int_t ibin) {
	return (check(ibin) ? (fList.at(ibin) - fList.at(ibin-1)) : 0.0);
}

Double_t MgntROOT::Axis::center(UInt_t ibin, Axis::Scl scl) {
	if (exist()) {
		if (ibin == 0)            return fList.at(ibin);
		if (ibin == fList.size()) return fList.at(ibin-1);
	}
	switch (scl) {
		case Axis::kLinear : 
			return (check(ibin) ? 0.5 * (fList.at(ibin-1) + fList.at(ibin)) : 0.0); break;
		case Axis::kLog : 
			return (check(ibin) ? std::sqrt(fList.at(ibin-1) * fList.at(ibin)) : 0.0); break;
		default : return 0.0; break;
	}
}

Int_t MgntROOT::Axis::find(Double_t val) {
	if (!exist()) return std::numeric_limits<Int_t>::max();
	Short_t bdlw = MgntNum::Compare(val, min());
	Short_t bdup = MgntNum::Compare(val, max());
	if (bdlw  < 0) return 0;
	if (bdup >= 0) return fList.size();
	UInt_t iter = 0;
	Int_t range[2] = { 0, nbin() };
	while ((range[1] - range[0]) != 1) {
		Int_t mid = (range[0] + range[1]) /2;
		if      (MgntNum::Compare(val, fList.at(mid))  < 0) range[1] = mid;
		else if (MgntNum::Compare(val, fList.at(mid)) >= 0) range[0] = mid;
		if (iter++ > 1e3) MgntSys::Error("<< Axis >> iter %lu", iter);
	}
	Int_t bin = range[1];
	return bin;
}

void MgntROOT::Axis::print(Axis::Scl scl, std::ostream& out) {
	if (!exist()) return;
	out << ("*****************************************  AxisBin  *****************************************\n");
	out << Form("AxisBin : TITLE(\"%s\") NBIN(%03d) REG[%16.8f %16.8f]\n", fTitle.c_str(), nbin(), min(), max());
	for (Int_t it = 1; it <= nbin(); ++it)
		out << Form("BIN(%03d) CEN(%16.8f) REG[%16.8f %16.8f] WIDTH(%16.8f)\n",
		                  it, center(it, scl), fList.at(it-1), fList.at(it), width(it));
	out << ("*********************************************************************************************\n");
}

MgntROOT::Axis::Axis(const Axis& axis) {
	*this = axis;
}

MgntROOT::Axis::Axis(Axis& axis, UInt_t mergeFT) {
	if (!axis.exist()) { MgntSys::Error("<< Axis >>  Axis is not exist."); return; }
	std::vector<Double_t> list(axis.bins(), axis.bins() + axis.nbin() + 1);
	if (merge(list, mergeFT)) fTitle = axis.title();
}

MgntROOT::Axis::Axis(std::initializer_list<Double_t> list) {
	if (!check(list)) return;
	if (merge(list)) fTitle = "";
}

MgntROOT::Axis::Axis(const std::string& title, const std::vector<Double_t>& list, UInt_t mergeFT) {
	if (!check(list)) return;
	if (merge(list, mergeFT)) fTitle = title;
}

 MgntROOT::Axis::Axis(const std::string& title, const UInt_t nbin, const Double_t * bins, UInt_t mergeFT) {
	if (bins == nullptr) { MgntSys::Error("<< Axis >>  bins is nullptr."); return; }
	std::vector<Double_t> list(bins, bins + nbin + 1);
	if (!check(list)) return;
	if (merge(list, mergeFT)) fTitle = title;
}

MgntROOT::Axis::Axis(const std::string& title, const UInt_t nbin, Double_t lw, Double_t up, Axis::Scl scl) {
	if (MgntSys::Error("<< Axis >>  nbin is less than 1.", (nbin < 1))) return;
	if (MgntSys::Error("<< Axis >>  lw >= up.", (MgntNum::Compare(lw, up) >= 0))) return;
	Bool_t scllog = (scl==Axis::kLog && (MgntNum::Compare(lw) <= 0 || MgntNum::Compare(up) <= 0));
	if (MgntSys::Error("<< Axis >>  lw or up <= 0. (Axis::kLog)", scllog)) return;
	std::vector<Double_t> list(nbin+1);
	switch (scl) {
		case Axis::kLinear :
			{
				Double_t step = (up - lw) / Double_t(nbin);
				UInt_t cnt = 0;
				for (Double_t& elem : list) {
					elem = (lw + cnt * step);
					cnt++;
				}
				break;
			}
		case Axis::kLog :
			{
				Double_t loglw = std::log(lw);
				Double_t logup = std::log(up);
				Double_t step = (logup - loglw) / Double_t(nbin);
				UInt_t cnt = 0;
				for (Double_t& elem : list) {
					elem = std::exp(loglw + cnt * step);
					cnt++;
				}
				break;
			}
		default : break;
	}
	if (merge(list)) fTitle = title;
}

Bool_t MgntROOT::Axis::init_TAxis(TAxis * axis, UInt_t mergeFT) {
	if (axis == nullptr) { MgntSys::Error("<< Axis >>  TAxis is nullptr."); return false; }
	UInt_t nbin = axis->GetNbins();
	std::vector<Double_t> list(nbin+1);
	for (Int_t ibin = 1; ibin <= nbin+1; ++ibin)
		list.at(ibin-1) = axis->GetBinLowEdge(ibin);
	
	if (merge(list, mergeFT)) fTitle = axis->GetTitle();
	return true;
}

Bool_t MgntROOT::Axis::init_TAxis(const std::string& title, TAxis * axis, UInt_t mergeFT) {
	if (init_TAxis(axis, mergeFT)) fTitle = title;
	return exist();
}

Bool_t MgntROOT::Axis::init_TH1(TH1 * hist, Axis::Dim dim, UInt_t mergeFT) {
	if (hist == nullptr) { MgntSys::Error("<< Axis >>  TH1 is nullptr."); return false; }
	Int_t ndim = hist->GetDimension();
	if (dim > ndim) { MgntSys::Error("<< Axis >>  require axis dimension too high."); return false; }
	TAxis * axis = nullptr;
	switch (dim) {
		case Axis::kX : axis = hist->GetXaxis(); break;
		case Axis::kY : axis = hist->GetYaxis(); break;
		case Axis::kZ : axis = hist->GetZaxis(); break;
		default : break;
	}
	return init_TAxis(axis, mergeFT);
}

Bool_t MgntROOT::Axis::init_TH1(const std::string& title, TH1 * hist, Axis::Dim dim, UInt_t mergeFT) {
	if (init_TH1(hist, dim, mergeFT)) fTitle = title;
	return exist();
}

Bool_t MgntROOT::Axis::init_TObject(TObject * obj, Axis::Dim dim, UInt_t mergeFT) {
	if (obj == nullptr) { MgntSys::Error("<< Axis >>  TObject is nullptr."); return false; }
	if (!dynamic_cast<TH1*>(obj)) { MgntSys::Error("<< Axis >>  TObject (%s) isn't a histogram.", obj->ClassName()); return false; }
	TH1 * hist = (TH1*)obj;
	return init_TH1(hist, dim, mergeFT);
}

Bool_t MgntROOT::Axis::init_TObject(const std::string& title, TObject * obj, Axis::Dim dim, UInt_t mergeFT) {
	if (init_TObject(obj, dim, mergeFT)) fTitle = title;
	return exist();
}

MgntROOT::Axis::Axis(const std::string& title, TAxis * axis, UInt_t mergeFT) {
	init_TAxis(title, axis, mergeFT);
}

MgntROOT::Axis::Axis(TAxis * axis, UInt_t mergeFT) {
	init_TAxis(axis, mergeFT);
}

MgntROOT::Axis::Axis(const std::string& title, TH1 * hist, Axis::Dim dim, UInt_t mergeFT) {
	init_TH1(title, hist, dim, mergeFT);
}

MgntROOT::Axis::Axis(TH1 * hist, Axis::Dim dim, UInt_t mergeFT) {
	init_TH1(hist, dim, mergeFT);
}

MgntROOT::Axis::Axis(const std::string& title, TObject * obj, Axis::Dim dim, UInt_t mergeFT) {
	init_TObject(title, obj, dim, mergeFT);
}

MgntROOT::Axis::Axis(TObject * obj, Axis::Dim dim, UInt_t mergeFT) {
	init_TObject(obj, dim, mergeFT);
}
		
MgntROOT::Axis::Axis(const std::string& title, MgntROOT::Hist * hist, Axis::Dim dim, UInt_t mergeFT) {
	if (hist == nullptr || !hist->exist()) return;
	init_TH1(title, (*hist)(), dim, mergeFT);
}

MgntROOT::Axis::Axis(MgntROOT::Hist * hist, Axis::Dim dim, UInt_t mergeFT) {
	if (hist == nullptr || !hist->exist()) return;
	init_TH1((*hist)(), dim, mergeFT);
}

Bool_t MgntROOT::Axis::check(const std::vector<Double_t>& list) {
	Bool_t result = true;
	if (result)
		for (const Double_t& elem : list) {
			if (!MgntNum::Valid(elem)) { result = false; break; }
		}
	if (result)
		for (Int_t it = 1; it < list.size(); ++it)
			if (MgntNum::Compare(list.at(it), list.at(it-1)) <= 0) { result = false; break; }
	if (!result) MgntSys::Error("<< Axis >> check list failure.");
	return result;
}

Bool_t MgntROOT::Axis::check(UInt_t ibin) {
	Bool_t result = (exist() && ibin >= 1 && ibin <= fList.size()-1);
	if (!result) MgntSys::Error("<< Axis >> check bin failure.");
	return result;
}

Bool_t MgntROOT::Axis::merge(const std::vector<Double_t>& list, UInt_t mergeFT) {
	Bool_t result = !(mergeFT < 1 || list.size() < 2 || ((list.size()-1) % mergeFT) != 0);
	if (!result) { MgntSys::Error("<< Axis >> merge failure."); return result; }
	if (mergeFT == 1) fList = list;
	else {
		UInt_t nbins = ((list.size()-1) / mergeFT) + 1;
		UInt_t cnt = 0;
		fList.resize(nbins);
		for (Double_t& elem : fList) {
			elem = list.at(cnt);
			cnt += mergeFT;
		}
	}
	return result;
}


/****************/
/****  Hist  ****/
/****************/
MgntROOT::Hist * MgntROOT::Hist::Head(const std::string& name) {
	if (name == "") return nullptr;
	std::map<std::string, Hist*>::iterator search = fHistMap.find(name);
	MgntROOT::Hist * hist = ((search != fHistMap.end()) ? search->second : nullptr);
	if (Hist::fDebug && hist == nullptr) MgntSys::Error(StrFmt("<< MgntRoot::Hist >> %s is not found.", name.c_str())); 
	return hist;
}

Bool_t MgntROOT::Hist::Exist(const std::string& name) {
	return (Head(name) != nullptr);
}

void MgntROOT::Hist::Delete(const std::string& name) {
	Hist * hist = Hist::Head(name);
	if (hist != nullptr) { hist->clear(); delete hist; }
}

void MgntROOT::Hist::Delete(Hist * hist) {
	if (hist != nullptr) { hist->clear(); delete hist; }
}

void MgntROOT::Hist::Delete(std::vector<Hist *>& histvec) {
	for (auto* hist : histvec)
		Hist::Delete(hist);
}

void MgntROOT::Hist::Load(const std::string& path) {
	if (!MgntSys::TestFile(path, 'f')) {
		MgntSys::Error(StrFmt("<< MgntROOT::Hist >>  %s is not open.", path.c_str()));
		return;
	}

	TFile * file = TFile::Open(path.c_str(), "READ");
	if (file != nullptr && file->IsOpen()) {
		TList * list = file->GetListOfKeys();
		TIter next((TList*)list);
		while (TObject * key = next()) {
			TString keyName = key->GetName();
			TH1 * obj = dynamic_cast<TH1*>(file->Get(keyName.Data()));
			if (obj == nullptr) continue;
			MgntROOT::Hist * hist = new MgntROOT::Hist(obj);
		}
		file->Close();
	}
}

void MgntROOT::Hist::Write() {
	for (auto const& it : fHistMap) (it.second)->write();
}


void MgntROOT::Hist::Save(const std::string& filename, const std::string& filepath, const std::string& dirname) {
	std::string fullpath = Form("%s/%s.root", filepath.c_str(), filename.c_str());
	if (!MgntSys::TestFile(filepath, 'd')) return;
	TDirectory * gdir = gDirectory;
	TFile * file = new TFile(fullpath.c_str(), "RECREATE");
	file->cd();
	TDirectoryFile * dir = nullptr;
	if (dirname != "") {
		dir = new TDirectoryFile(dirname.c_str(), dirname.c_str());
		dir->cd();
	}
	
	MgntROOT::Hist::Write();

	file->Write();
	file->Close();
	if (dir != nullptr) dir->Delete();
	file->Delete();
	gdir->cd();
}

void MgntROOT::Hist::Print(std::ostream& out) {
	out << Form("*********************************  Histogram List  *********************************\n");
	out << Form("%-8s   %12s   %-18s   %-50s\n", "UniqueID", "ClassName", "Name", "Title");
	for (auto const& it : fHistMap) { 
		const TH1 * hist = (*(it.second))();
		std::cout << Form("%08u   %12s   %-18s   %-50s\n", hist->GetUniqueID(), hist->ClassName(), hist->GetName(), hist->GetTitle());
	}
	out << Form("************************************************************************************\n");
}
		
MgntROOT::Hist * MgntROOT::Hist::Calculate(Hist::ArithOpt opt, const std::string& name, const std::string& title, Hist * hNum, Hist * hDen) {
	if (hNum == nullptr || 
	   (*hNum)()->GetEntries() == 0 || 
		 !MgntNum::Valid((*hNum)()->Integral())) { 
		MgntSys::Error("<< MgntROOT::Hist >>  hNum is nullptr or NULL/INF/NAN."); return nullptr;
	}
	if (hDen == nullptr || 
	    (*hDen)()->GetEntries() == 0 || 
			!MgntNum::Valid((*hDen)()->Integral())) { 
		MgntSys::Error("<< MgntROOT::Hist >>  hDen is nullptr or NULL/INF/NAN."); return nullptr; 
	}

	TH1 * hNumDen = (TH1*) ((*hNum)()->Clone(name.c_str()));
	hNumDen->SetNameTitle(name.c_str(), title.c_str());
	switch (opt) {
		case Hist::kAddition : hNumDen->Add((*hDen)(),  1.0); break;
		case Hist::kSubtract : hNumDen->Add((*hDen)(), -1.0); break;
		case Hist::kMultiply : hNumDen->Multiply((*hDen)());  break;
		case Hist::kDivide   : hNumDen->Divide((*hDen)());    break;
		default : break;
	}
	Hist * hist = new Hist(hNumDen);
	hNumDen->Delete();
	return hist;
}
	
MgntROOT::Hist * MgntROOT::Hist::Project(Hist::ProjOpt opt, Hist * hmon, Int_t ibin, Int_t jbin) {
	if (hmon == nullptr || hmon->isProfile() || hmon->ndim() != 2) return nullptr;
	Axis * axis = (opt == Hist::kProjX) ? hmon->axisY() : hmon->axisX();
	Int_t sbin = ((ibin <= axis->nbin()+1) ? ibin : axis->nbin()+1);
	Int_t ebin = ((jbin <= axis->nbin()+1) ? jbin : axis->nbin()+1);
	if (sbin < 0) sbin = ((sbin == -1) ? 1 : 0);
	if (ebin < 0) ebin = ((ebin == -1) ? 0 : 1) + axis->nbin();
	std::string name = StrFmt("%s__PROJ%c%05dTO%05d", hmon->name().c_str(), ((opt==Hist::kProjX)?'X':'Y'), sbin, ebin);

	Hist * hist = nullptr;
	if (opt == Hist::kProjX)
		hist = Hist::New(((TH2D*)(*hmon)())->ProjectionX(name.c_str(), sbin, ebin));
	if (opt == Hist::kProjY)
		hist = Hist::New(((TH2D*)(*hmon)())->ProjectionY(name.c_str(), sbin, ebin));
	
	return hist;
}

MgntROOT::Hist * MgntROOT::Hist::Project(Hist::ProjOpt opt, const std::string& nmon, Int_t ibin, Int_t jbin) {
	Hist * hmon = Hist::Head(nmon);
	return Hist::Project(opt, hmon, ibin, jbin);
}

std::vector<MgntROOT::Hist *> MgntROOT::Hist::Project(Hist::ProjOpt opt, Hist * hmon, UInt_t mergeFT) {
	std::vector<MgntROOT::Hist *> list;
	if (hmon == nullptr || hmon->isProfile() || hmon->ndim() != 2) return list;
	if (opt == Hist::kProjX) {
		Axis * axis = hmon->axisY();
		if (mergeFT == 0 || (mergeFT > axis->nbin()) || (axis->nbin() % mergeFT != 0)) return list;
		Int_t nstep = axis->nbin() / mergeFT;
		list.push_back(Hist::Project(opt, hmon, 0, 0));
		for (Int_t istep = 0; istep < nstep; ++istep) {
			Int_t ibin = (mergeFT * istep) + 1;
			Int_t jbin = (mergeFT * (istep + 1));
			list.push_back(Hist::Project(opt, hmon, ibin, jbin));
		}
		list.push_back(Hist::Project(opt, hmon, axis->nbin()+1, axis->nbin()+1));
	}
	else if (opt == Hist::kProjY) {
		Axis * axis = hmon->axisX();
		if (mergeFT == 0 || (mergeFT > axis->nbin()) || (axis->nbin() % mergeFT != 0)) return list;
		Int_t nstep = axis->nbin() / mergeFT;
		list.push_back(Hist::Project(opt, hmon, 0, 0));
		for (Int_t istep = 0; istep < nstep; ++istep) {
			Int_t ibin = (mergeFT * istep) + 1;
			Int_t jbin = (mergeFT * (istep + 1));
			list.push_back(Hist::Project(opt, hmon, ibin, jbin));
		}
		list.push_back(Hist::Project(opt, hmon, axis->nbin()+1, axis->nbin()+1));
	}
	return list;
}

std::vector<MgntROOT::Hist *> MgntROOT::Hist::Project(Hist::ProjOpt opt, const std::string& nmon, UInt_t mergeFT) {
	Hist * hmon = Hist::Head(nmon.c_str());
	return Hist::Project(opt, hmon, mergeFT);
}
		
MgntROOT::Hist * MgntROOT::Hist::ProjectAll(Hist::ProjOpt opt, Hist * hmon, Bool_t lwflow, Bool_t upflow) {
	return Hist::Project(opt, hmon, (lwflow?-2:-1), (upflow?-2:-1));
}

MgntROOT::Hist * MgntROOT::Hist::ProjectAll(Hist::ProjOpt opt, const std::string& nmon, Bool_t lwflow, Bool_t upflow) {
	Hist * hmon = Hist::Head(nmon.c_str());
	return Hist::Project(opt, hmon, (lwflow?-2:-1), (upflow?-2:-1));
}

MgntROOT::Hist * MgntROOT::Hist::Project3D(Hist::ProjOpt opt, Hist * hmon, Int_t ibin, Int_t jbin) {
	if (hmon == nullptr || hmon->isProfile() || hmon->ndim() != 3) return nullptr;
	Axis * axis = (opt == Hist::kProjYZ) ? hmon->axisX() : hmon->axisX();
	Int_t sbin = ((ibin <= axis->nbin()+1) ? ibin : axis->nbin()+1);
	Int_t ebin = ((jbin <= axis->nbin()+1) ? jbin : axis->nbin()+1);
	if (sbin < 0) sbin = ((sbin == -1) ? 1 : 0);
	if (ebin < 0) ebin = ((ebin == -1) ? 0 : 1) + axis->nbin();
	std::string name = StrFmt("%s__PROJ%c%05dTO%05d", hmon->name().c_str(), ((opt==Hist::kProjX)?'X':'X'), sbin, ebin);

	Hist * hist = nullptr;
	//if (Hist::Exist(name))
	//	hist = Hist::Head(name);
	//else if (opt == Hist::kProjYZ) {
	if (opt == Hist::kProjYZ) {
		(*hmon)()->GetXaxis()->SetRange(sbin, ebin);
		TH2D * hson = (TH2D*)(((TH3D*)((*hmon)()))->Project3D("zy"));
		hson->SetName(name.c_str());
		(*hmon)()->GetXaxis()->UnZoom();
		hist = Hist::New(hson);
	}
	
	return hist;
}

MgntROOT::Hist * MgntROOT::Hist::Project3D(Hist::ProjOpt opt, const std::string& nmon, Int_t ibin, Int_t jbin) {
	Hist * hmon = Hist::Head(nmon);
	return Hist::Project3D(opt, hmon, ibin, jbin);
}

std::vector<MgntROOT::Hist *> MgntROOT::Hist::Project3D(Hist::ProjOpt opt, Hist * hmon, UInt_t mergeFT) {
	std::vector<MgntROOT::Hist *> list;
	if (hmon == nullptr || hmon->isProfile() || hmon->ndim() != 3) return list;
	if (opt == Hist::kProjYZ) {
		Axis * axis = hmon->axisX();
		if (mergeFT == 0 || (mergeFT > axis->nbin()) || (axis->nbin() % mergeFT != 0)) return list;
		Int_t nstep = axis->nbin() / mergeFT;
		list.push_back(Hist::Project(opt, hmon, 0, 0));
		for (Int_t istep = 0; istep < nstep; ++istep) {
			Int_t ibin = (mergeFT * istep) + 1;
			Int_t jbin = (mergeFT * (istep + 1));
			list.push_back(Hist::Project3D(opt, hmon, ibin, jbin));
		}
		list.push_back(Hist::Project3D(opt, hmon, axis->nbin()+1, axis->nbin()+1));
	}
	return list;
}

std::vector<MgntROOT::Hist *> MgntROOT::Hist::Project3D(Hist::ProjOpt opt, const std::string& nmon, UInt_t mergeFT) {
	Hist * hmon = Hist::Head(nmon.c_str());
	return Hist::Project3D(opt, hmon, mergeFT);
}

MgntROOT::Hist * MgntROOT::Hist::Project3DAll(Hist::ProjOpt opt, Hist * hmon, Bool_t lwflow, Bool_t upflow) {
	return Hist::Project3D(opt, hmon, (lwflow?-2:-1), (upflow?-2:-1));
}

MgntROOT::Hist * MgntROOT::Hist::Project3DAll(Hist::ProjOpt opt, const std::string& nmon, Bool_t lwflow, Bool_t upflow) {
	Hist * hmon = Hist::Head(nmon.c_str());
	return Hist::Project3D(opt, hmon, (lwflow?-2:-1), (upflow?-2:-1));
}
		
THStack * MgntROOT::Hist::Collect(const std::string& name, const std::string& title, const Hist::LIST& list) {
	if (list.size() == 0) return nullptr;
	THStack * coll = new THStack();
	coll->SetNameTitle(name.c_str(), title.c_str());
	for (const std::string& elem : list) {
		Hist * hist = Hist::Head(elem);
		if (hist == nullptr) continue;
		coll->Add((*hist)());
	}
	return coll;
}

THStack * MgntROOT::Hist::Collect(const std::string& name, const std::string& title, const Hist::VLIST& vlist) {
	if (vlist.size() == 0) return nullptr;
	THStack * coll = new THStack();
	coll->SetNameTitle(name.c_str(), title.c_str());
	for (Hist * hist : vlist) {
		if (hist == nullptr) continue;
		coll->Add((*hist)());
	}
	return coll;
}
		
MgntROOT::Hist * MgntROOT::Hist::New(TH1 * hist, Bool_t reset) {
	if (hist == nullptr) return nullptr;
	new Hist(hist, reset);
	return Hist::Head(hist->GetName());
}

MgntROOT::Hist * MgntROOT::Hist::New(const std::string& name, TH1 * hist, Bool_t reset) {
	if (hist == nullptr) return nullptr;
	new Hist(name, hist, reset);
	return Hist::Head(name);
}

MgntROOT::Hist * MgntROOT::Hist::New(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, Bool_t isProfile) {
	new Hist(name, title, axisX, isProfile);
	return Hist::Head(name);
}

MgntROOT::Hist * MgntROOT::Hist::New(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, Bool_t isProfile2D) {
	new Hist(name, title, axisX, axisY, isProfile2D);
	return Hist::Head(name);
}

MgntROOT::Hist * MgntROOT::Hist::New(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, MgntROOT::Axis& axisZ, Bool_t isProfile3D) {
	new Hist(name, title, axisX, axisY, axisZ, isProfile3D);
	return Hist::Head(name);
}

MgntROOT::Hist::Hist(TH1 * hist, Bool_t reset) : Hist() {
	if (hist == nullptr) {
		MgntSys::Error("<< MgntROOT::Hist >>  histogram is not exist.");
		return;
	}
	Hist::DebugOFF();
	if (isHistExist(hist->GetName())) 
		{ Hist::DebugON(); return; }
	Hist::DebugON();

	fName = hist->GetName();
	fTitle = hist->GetTitle();
	setType(hist);

	fHist = (TH1*) hist->Clone(hist->GetName());
	if (reset) {
		fHist->GetListOfFunctions()->Clear();
		fHist->Reset();
	}

	if (fHist->GetDimension() >= 1) fAxisX = MgntROOT::Axis(fHist, MgntROOT::Axis::kX);
	if (fHist->GetDimension() >= 2) fAxisY = MgntROOT::Axis(fHist, MgntROOT::Axis::kY);
	if (fHist->GetDimension() >= 3) fAxisZ = MgntROOT::Axis(fHist, MgntROOT::Axis::kZ);

	setStyle();
	pushBack();
}

MgntROOT::Hist::Hist(const std::string& name, TH1 * hist, Bool_t reset) : Hist() {
	if (hist == nullptr) {
		MgntSys::Error("<< MgntROOT::Hist >>  histogram is not exist.");
		return;
	}
	Hist::DebugOFF();
	if (isHistExist(name))
		{ Hist::DebugON(); return; }
	Hist::DebugON();
	
	fName = name;
	fTitle = hist->GetTitle();
	setType(hist);

	fHist = (TH1*) hist->Clone(name.c_str());
	if (reset) {
		fHist->GetListOfFunctions()->Clear();
		fHist->Reset();
	}

	if (fHist->GetDimension() >= 1) fAxisX = MgntROOT::Axis(fHist, MgntROOT::Axis::kX);
	if (fHist->GetDimension() >= 2) fAxisY = MgntROOT::Axis(fHist, MgntROOT::Axis::kY);
	if (fHist->GetDimension() >= 3) fAxisZ = MgntROOT::Axis(fHist, MgntROOT::Axis::kZ);

	setStyle();
	pushBack();
}

MgntROOT::Hist::Hist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, Bool_t isProfile) : Hist() {
	Hist::DebugOFF();
	if (isHistExist(name))
		{ Hist::DebugON(); return; }
	Hist::DebugON();
	if (!axisX.exist()) {
		MgntSys::Error("<< MgntROOT::Hist >>  axis is not exist.");
		return;
	}
	
	std::string tle = "";
	std::vector<std::string>&& strvec = MgntRegex::Split(title, MgntRegex::Formula::Semicolon);
	if      (strvec.size() == 0) tle = Form(";%s;", axisX.title().c_str());
	else if (strvec.size() == 1) tle = Form("%s;%s;", strvec.at(0).c_str(), axisX.title().c_str());
	else if (strvec.size() == 2) tle = Form("%s;%s;%s", strvec.at(0).c_str(), axisX.title().c_str(), strvec.at(1).c_str());
	else { 
		MgntSys::Error("<< MgntROOT::Hist >>  title formula is failure.");
		return;
	}

	TH1 * hist = nullptr;
	if (isProfile)
		hist = new TProfile(name.c_str(), tle.c_str(), axisX.nbin(), axisX.bins());
	else
		hist = new TH1D(name.c_str(), tle.c_str(), axisX.nbin(), axisX.bins());
	fAxisX = axisX;

	fName = name;
	fTitle = hist->GetTitle();
	fType.first = 1;
	fType.second = isProfile;
	fHist = hist;
	setStyle();
	pushBack();
}
		
MgntROOT::Hist::Hist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, Bool_t isProfile2D) : Hist() {
	Hist::DebugOFF();
	if (isHistExist(name))
		{ Hist::DebugON(); return; }
	Hist::DebugON();
	if (!axisX.exist() || !axisY.exist()) {
		MgntSys::Error("<< MgntROOT::Hist >>  axis is not exist.");
		return;
	}
	
	std::string tle = "";
	std::vector<std::string>&& strvec = MgntRegex::Split(title, MgntRegex::Formula::Semicolon);
	if      (strvec.size() == 0) tle = Form(";%s;%s;", axisX.title().c_str(), axisY.title().c_str());
	else if (strvec.size() == 1) tle = Form("%s;%s;%s", strvec.at(0).c_str(), axisX.title().c_str(), axisY.title().c_str());
	else if (strvec.size() == 2) tle = Form("%s;%s;%s;%s", strvec.at(0).c_str(), axisX.title().c_str(), axisY.title().c_str(), strvec.at(1).c_str());
	else { 
		MgntSys::Error("<< MgntROOT::Hist >>  title formula is failure.");
		return;
	}

	TH1 * hist = nullptr;
	if (isProfile2D)
		hist = new TProfile2D(name.c_str(), tle.c_str(), axisX.nbin(), axisX.bins(), axisY.nbin(), axisY.bins());
	else
		hist = new TH2D(name.c_str(), tle.c_str(), axisX.nbin(), axisX.bins(), axisY.nbin(), axisY.bins());
	fAxisX = axisX;
	fAxisY = axisY;

	fName = name;
	fTitle = hist->GetTitle();
	fType.first = 2;
	fType.second = isProfile2D;
	fHist = hist;
	setStyle();
	pushBack();
}

MgntROOT::Hist::Hist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, MgntROOT::Axis& axisZ, Bool_t isProfile3D) : Hist() {
	Hist::DebugOFF();
	if (isHistExist(name))
		{ Hist::DebugON(); return; }
	Hist::DebugON();
	if (!axisX.exist() || !axisY.exist() | !axisZ.exist()) {
		MgntSys::Error("<< MgntROOT::Hist >>  axis is not exist.");
		return;
	}
	
	std::string tle = "";
	std::vector<std::string>&& strvec = MgntRegex::Split(title, MgntRegex::Formula::Semicolon);
	if      (strvec.size() == 0) tle = Form(";%s;%s;%s", axisX.title().c_str(), axisY.title().c_str(), axisZ.title().c_str());
	else if (strvec.size() == 1) tle = Form("%s;%s;%s;%s", strvec.at(0).c_str(), axisX.title().c_str(), axisY.title().c_str(), axisZ.title().c_str());
	else if (strvec.size() == 2) tle = Form("%s;%s;%s;%s;%s", strvec.at(0).c_str(), axisX.title().c_str(), axisY.title().c_str(), axisZ.title().c_str(), strvec.at(1).c_str());
	else { 
		MgntSys::Error("<< MgntROOT::Hist >>  title formula is failure.");
		return;
	}

	TH1 * hist = nullptr;
	if (isProfile3D)
		hist = new TProfile3D(name.c_str(), tle.c_str(), axisX.nbin(), axisX.bins(), axisY.nbin(), axisY.bins(), axisZ.nbin(), axisZ.bins());
	else
		hist = new TH3D(name.c_str(), tle.c_str(), axisX.nbin(), axisX.bins(), axisY.nbin(), axisY.bins(), axisZ.nbin(), axisZ.bins());
	fAxisX = axisX;
	fAxisY = axisY;
	fAxisZ = axisZ;

	fName = name;
	fTitle = hist->GetTitle();
	fType.first = 3;
	fType.second = false;
	fHist = hist;
	setStyle();
	pushBack();
}
		
void MgntROOT::Hist::fill(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e) {
	if (!exist()) return;
	switch (fType.first) {
		case 1 :
			if (fType.second) ((TProfile*)fHist)->Fill(a, b, c);
			else              ((TH1D*)    fHist)->Fill(a, b);
			break;
		case 2 :
			if (fType.second) ((TProfile2D*)fHist)->Fill(a, b, c, d);
			else              ((TH2D*)      fHist)->Fill(a, b, c);
			break;
		case 3 :
			if (fType.second) ((TProfile3D*)fHist)->Fill(a, b, c, d, e);
			else              ((TH3D*)      fHist)->Fill(a, b, c, d);
			break;
		default : break;
	}
}
		
void MgntROOT::Hist::setContent(Double_t content, Int_t ibin, Int_t jbin, Int_t kbin) {
	if (!exist()) return;
	switch (fType.first) {
		case 1 :
			if (fType.second) ((TProfile*)fHist)->SetBinContent(ibin, content);
			else              ((TH1D*)    fHist)->SetBinContent(ibin, content);
			break;
		case 2 :
			if (fType.second) ((TProfile2D*)fHist)->SetBinContent(ibin, jbin, content);
			else              ((TH2D*)      fHist)->SetBinContent(ibin, jbin, content);
			break;
		case 3 :
			if (fType.second) ((TProfile3D*)fHist)->SetBinContent(ibin, jbin, kbin, content);
			else              ((TH3D*)      fHist)->SetBinContent(ibin, jbin, kbin, content);
			break;
		default : break;
	}
}

void MgntROOT::Hist::setError(Double_t error, Int_t ibin, Int_t jbin, Int_t kbin) {
	if (!exist()) return;
	switch (fType.first) {
		case 1 :
			if (fType.second) ((TProfile*)fHist)->SetBinError(ibin, error);
			else              ((TH1D*)    fHist)->SetBinError(ibin, error);
			break;
		case 2 :
			if (fType.second) ((TProfile2D*)fHist)->SetBinError(ibin, jbin, error);
			else              ((TH2D*)      fHist)->SetBinError(ibin, jbin, error);
			break;
		case 3 :
			if (fType.second) ((TProfile3D*)fHist)->SetBinError(ibin, jbin, kbin, error);
			else              ((TH3D*)      fHist)->SetBinError(ibin, jbin, kbin, error);
			break;
		default : break;
	}
}

void MgntROOT::Hist::setContentWithError(Double_t content, Double_t error, Int_t ibin, Int_t jbin, Int_t kbin) {
	if (!exist()) return;
	setContent(content, ibin, jbin, kbin);
	setError(error, ibin, jbin, kbin);
}
		
Double_t MgntROOT::Hist::getContent(Int_t ibin, Int_t jbin, Int_t kbin) {
	if (!exist()) return 0.;
	Double_t content = 0.;
	switch (fType.first) {
		case 1 :
			if (fType.second) content = ((TProfile*)fHist)->GetBinContent(ibin);
			else              content = ((TH1D*)    fHist)->GetBinContent(ibin);
			break;
		case 2 :
			if (fType.second) content = ((TProfile2D*)fHist)->GetBinContent(ibin, jbin);
			else              content = ((TH2D*)      fHist)->GetBinContent(ibin, jbin);
			break;
		case 3 :
			if (fType.second) content = ((TProfile3D*)fHist)->GetBinContent(ibin, jbin, kbin);
			else              content = ((TH3D*)      fHist)->GetBinContent(ibin, jbin, kbin);
			break;
		default : break;
	}
	return content;
}

Double_t MgntROOT::Hist::getError(Int_t ibin, Int_t jbin, Int_t kbin) {
	if (!exist()) return 0.;
	Double_t error = 0.;
	switch (fType.first) {
		case 1 :
			if (fType.second) error = ((TProfile*)fHist)->GetBinError(ibin);
			else              error = ((TH1D*)    fHist)->GetBinError(ibin);
			break;
		case 2 :
			if (fType.second) error = ((TProfile2D*)fHist)->GetBinError(ibin, jbin);
			else              error = ((TH2D*)      fHist)->GetBinError(ibin, jbin);
			break;
		case 3 :
			if (fType.second) error = ((TProfile3D*)fHist)->GetBinError(ibin, jbin, kbin);
			else              error = ((TH3D*)      fHist)->GetBinError(ibin, jbin, kbin);
			break;
		default : break;
	}
	return error;
}

std::pair<Double_t, Double_t> MgntROOT::Hist::getContentWithError(Int_t ibin, Int_t jbin, Int_t kbin) {
	if (!exist()) return std::make_pair(0., 0.);
	return std::make_pair(getContent(ibin, jbin, kbin), getError(ibin, jbin, kbin));
}
		
void MgntROOT::Hist::scale(Double_t scl, Option_t * option) {
	if (!exist()) return;
	fHist->Scale(scl, option);
}

void MgntROOT::Hist::normalized(Hist::NormOpt opt) {
	if (!exist()) return;
	if (!fHist->GetDefaultSumw2()) fHist->Sumw2();
	switch (opt) {
		case Hist::kEntries :
	    {
			  Double_t entries = (Double_t)fHist->GetEntries();
			  if (MgntNum::Compare(entries) > 0) fHist->Scale(1. / entries);
	    }
			break;
		case Hist::kIntegral :
	    {
			  Double_t integral = (Double_t)fHist->Integral();
			  if (MgntNum::Compare(integral) > 0) fHist->Scale(1. / integral);
	    }
			break;
		case Hist::kArea :
	    {
			  Double_t area = (Double_t)fHist->Integral("width");
			  if (MgntNum::Compare(area) > 0) fHist->Scale(1. / area);
	    }
			break;
		default : break;
	}
}

void MgntROOT::Hist::setArea(MgntROOT::Area area) {
	if (!exist()) return;
	fHist->SetFillColor(area().GetFillColor());
	fHist->SetFillStyle(area().GetFillStyle());
}

void MgntROOT::Hist::setLine(MgntROOT::Line line) {
	if (!exist()) return;
	fHist->SetLineColor(line().GetLineColor());
	fHist->SetLineStyle(line().GetLineStyle());
	fHist->SetLineWidth(line().GetLineWidth());
}

void MgntROOT::Hist::setMarker(MgntROOT::Marker marker) {
	if (!exist()) return;
	fHist->SetMarkerColor(marker().GetMarkerColor());
	fHist->SetMarkerStyle(marker().GetMarkerStyle());
	fHist->SetMarkerSize (marker().GetMarkerSize());
}
		
void MgntROOT::Hist::setStyle(MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker) {
	if (!exist()) return;
	setArea(area);
	setLine(line);
	setMarker(marker);
}
		
void MgntROOT::Hist::draw(Option_t * option,  MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker) {
	if (!exist()) return;
	setStyle(area, line, marker);
	fHist->Draw(option);
}

inline Bool_t MgntROOT::Hist::isHistExist(const std::string& name) { 
	if (MgntROOT::Hist::Exist(name)) { 
		MgntSys::Error(Form("<< MgntROOT::Hist >> %s is already exist.", name.c_str())); 
		return true; 
	}
	else return false;
}
	
void MgntROOT::Hist::setType(TH1 * hist) {
	if (hist == nullptr) return;
	fType.first = hist->GetDimension();
	switch (fType.first) {
		case 1 :
			if (hist->ClassName() == TString("TProfile")) fType.second = true;
			else                                          fType.second = false;
			break;
		case 2 :
			if (hist->ClassName() == TString("TProfile2D")) fType.second = true;
			else                                            fType.second = false;
			break;
		case 3 :
			if (hist->ClassName() == TString("TProfile3D")) fType.second = true;
			else                                            fType.second = false;
			break;
		default : break;
	}
}


/*****************/
/****  Graph  ****/
/*****************/
MgntROOT::Graph * MgntROOT::Graph::Head(const std::string& name) {
	if (name == "") return nullptr;
	std::map<std::string, Graph*>::iterator search = fGraphMap.find(name);
	MgntROOT::Graph * graph = ((search != fGraphMap.end()) ? search->second : nullptr);
	return graph;
}

Bool_t MgntROOT::Graph::Exist(const std::string& name) {
	return (Head(name) != nullptr);
}

void MgntROOT::Graph::Delete(const std::string& name) {
	Graph * graph = Graph::Head(name);
	if (graph != nullptr) { graph->clear(); delete graph; }
}

void MgntROOT::Graph::Delete(Graph * graph) {
	if (graph != nullptr) { graph->clear(); delete graph; }
}

void MgntROOT::Graph::Delete(std::vector<Graph *>& graphvec) {
	for (auto* graph : graphvec)
		Graph::Delete(graph);
}

void MgntROOT::Graph::Load(const std::string& path) {
	if (!MgntSys::TestFile(path, 'f')) {
		MgntSys::Error(Form("<< MgntROOT::Graph >>  %s is not open.", path.c_str()));
		return;
	}

	TFile * file = TFile::Open(path.c_str(), "READ");
	if (file != nullptr && file->IsOpen()) {
		TList * list = file->GetListOfKeys();
		TIter next((TList*)list);
		while (TObject * key = next()) {
			TString keyName = key->GetName();
			TGraph * obj = dynamic_cast<TGraph*>(file->Get(keyName.Data()));
			if (obj == nullptr) continue;
			MgntROOT::Graph * graph = new MgntROOT::Graph(obj);
		}
		file->Close();
	}
}

void MgntROOT::Graph::Write() {
	for (auto const& it : fGraphMap) (it.second)->write();
}

void MgntROOT::Graph::Save(const std::string& filename, const std::string& filepath, const std::string& dirname) {
	std::string fullpath = Form("%s/%s.root", filepath.c_str(), filename.c_str());
	if (!MgntSys::TestFile(filepath, 'd')) return;
	TDirectory * gdir = gDirectory;
	TFile * file = new TFile(fullpath.c_str(), "RECREATE");
	file->cd();
	TDirectoryFile * dir = nullptr;
	if (dirname != "") {
		dir = new TDirectoryFile(dirname.c_str(), dirname.c_str());
		dir->cd();
	}
	
	MgntROOT::Graph::Write();

	file->Write();
	file->Close();
}

void MgntROOT::Graph::Print(std::ostream& out) {
	out << Form("***********************************  Graph List  ***********************************\n");
	out << Form("%-8s   %-18s   %-50s\n", "UniqueID", "Name", "Title");
	for (auto const& it : fGraphMap) { 
		const TGraphAsymmErrors * graph = (*(it.second))();
		std::cout << Form("%08u   %-18s   %-50s\n", graph->GetUniqueID(), graph->GetName(), graph->GetTitle());
	}
	out << Form("************************************************************************************\n");
}


TMultiGraph * MgntROOT::Graph::Collect(const std::string& name, const std::string& title, const Graph::LIST& list) {
	if (list.size() == 0) return nullptr;
	TMultiGraph * coll = new TMultiGraph();
	coll->SetName(name.c_str());
	coll->SetTitle(title.c_str());
	for (const std::string& elem : list) {
		Graph * graph = Graph::Head(elem);
		if (graph == nullptr) continue;
		coll->Add((*graph)());
	}
	return coll;
}
		
TMultiGraph * MgntROOT::Graph::Collect(const std::string& name, const std::string& title, const Graph::VLIST& vlist) {
	if (vlist.size() == 0) return nullptr;
	TMultiGraph * coll = new TMultiGraph();
	coll->SetName(name.c_str());
	coll->SetTitle(title.c_str());
	for (Graph * graph : vlist) {
		if (graph == nullptr) continue;
		coll->Add((*graph)());
	}
	return coll;
}

MgntROOT::Graph * MgntROOT::Graph::New(Hist * hist, Bool_t islogx) {
	if (hist == nullptr) return nullptr;
	new Graph(hist, islogx);
	return Graph::Head(hist->name());
}

MgntROOT::Graph * MgntROOT::Graph::New(TGraph * graph, Bool_t reset) {
	if (graph == nullptr) return nullptr;
	new Graph(graph, reset);
	return Graph::Head(graph->GetName());
}

MgntROOT::Graph * MgntROOT::Graph::New(const std::string& name, const std::string& title, const std::string& axtlex, const std::string& axtley) {
	new Graph(name, title, axtlex, axtley);
	return Graph::Head(name);
}
		
MgntROOT::Graph::Graph(Hist * hist, Bool_t islogx) : Graph() {
	if (hist == nullptr || !hist->exist() || hist->ndim() != 1) return;
	if (isGraphExist(hist->name())) return;
	fName = hist->name();
	fTitle = hist->title();
	fGraph = new TGraphAsymmErrors();
	fGraph->SetName(hist->name().c_str());
	fGraph->SetTitle(hist->title().c_str());
	fGraph->GetXaxis()->SetTitle(hist->axisX()->title().c_str());
	fGraph->GetYaxis()->SetTitle(hist->axisY()->title().c_str());


	MgntROOT::Axis * axis = hist->axisX();
	fGraph->GetXaxis()->SetRange(axis->min(), axis->max());	
	for (Int_t ibin = 1; ibin <= axis->nbin(); ++ibin) {
		Double_t cen = axis->center(ibin, (islogx ? MgntROOT::Axis::kLog : MgntROOT::Axis::kLinear));
		Double_t val = hist->getContent(ibin);
		Double_t err = hist->getError(ibin);
		fGraph->SetPoint(ibin-1, cen, val);
		fGraph->SetPointError(ibin-1, 0., 0., err, err);
	}

	setStyle();
	pushBack();
}

MgntROOT::Graph::Graph(TGraph * graph, Bool_t reset) : Graph() {
	if (graph == nullptr) {
		MgntSys::Error("<< MgntROOT::Graph >>  graph is not exist.");
		return;
	}
	if (isGraphExist(graph->GetName())) return;

	fName  = graph->GetName();
	fTitle = graph->GetTitle();
	fGraph = new TGraphAsymmErrors();
	fGraph->SetName(graph->GetName());
	fGraph->SetTitle(graph->GetTitle());
	fGraph->GetXaxis()->SetTitle(graph->GetXaxis()->GetTitle());
	fGraph->GetYaxis()->SetTitle(graph->GetYaxis()->GetTitle());

	if (!reset) {
		UInt_t type = 0;
		if      (std::string(graph->ClassName()) == std::string("TGraph")) type = 1;
		else if (std::string(graph->ClassName()) == std::string("TGraphErrors")) type = 2;
		else if (std::string(graph->ClassName()) == std::string("TGraphAsymmErrors")) type = 3;
		switch (type) {
			case 1 :
				for (Int_t it = 0; it < graph->GetN(); ++it) 
					setPoint(it, Graph::Point_t(graph->GetX()[it], graph->GetY()[it]));
				break;
			case 2 :
				for (Int_t it = 0; it < graph->GetN(); ++it) 
					setPointWithError(it, 
						Graph::Point_t(graph->GetX()[it], graph->GetY()[it]),
						Graph::Error_t(graph->GetEX()[it], graph->GetEY()[it])
					);
				break;
			case 3 :
				for (Int_t it = 0; it < graph->GetN(); ++it) 
					setPointWithError(it, 
						Graph::Point_t(graph->GetX()[it], graph->GetY()[it]),
						Graph::Error_t(graph->GetErrorXlow(it), graph->GetErrorXhigh(it), 
						               graph->GetErrorYlow(it), graph->GetErrorYhigh(it))
					);
				break;
			default : break;
		}
	}
	
	setStyle();
	pushBack();
}

MgntROOT::Graph::Graph(const std::string& name, const std::string& title, const std::string& axtlex, const std::string& axtley) : Graph() {
	if (isGraphExist(name)) return;
	TGraphAsymmErrors * graph = new TGraphAsymmErrors();
	graph->SetName(name.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle(axtlex.c_str());
	graph->GetYaxis()->SetTitle(axtley.c_str());
	
	fName  = name;
	fTitle = title;
	fGraph = graph;
	setStyle();
	pushBack();
}

void MgntROOT::Graph::setPoint(UInt_t it, const Graph::Point_t& pnt) {
	if (!exist()) return;
	fGraph->SetPoint(it, pnt.x, pnt.y);
}

void MgntROOT::Graph::setError(UInt_t it, const Graph::Error_t& err) {
	if (!exist()) return;
	fGraph->SetPointError(it, err.exl, err.exu, err.eyl, err.eyu);
}

void MgntROOT::Graph::setPointWithError(UInt_t it, const Graph::Point_t& pnt, const Graph::Error_t& err) {
	if (!exist()) return;
	fGraph->SetPoint(it, pnt.x, pnt.y);
	fGraph->SetPointError(it, err.exl, err.exu, err.eyl, err.eyu);
}

void MgntROOT::Graph::pushPoint(const Graph::Point_t& pnt) {
	if (!exist()) return;
	UInt_t it = fGraph->GetN();
	fGraph->SetPoint(it, pnt.x, pnt.y);
}

void MgntROOT::Graph::pushError(const Graph::Error_t& err) {
	if (!exist()) return;
	UInt_t it = fGraph->GetN();
	fGraph->SetPointError(it, err.exl, err.exu, err.eyl, err.eyu);
}

void MgntROOT::Graph::pushPointWithError(const Graph::Point_t& pnt, const Graph::Error_t& err) {
	if (!exist()) return;
	UInt_t it = fGraph->GetN();
	fGraph->SetPoint(it, pnt.x, pnt.y);
	fGraph->SetPointError(it, err.exl, err.exu, err.eyl, err.eyu);
}

void MgntROOT::Graph::setArea(MgntROOT::Area area) {
	if (!exist()) return;
	fGraph->SetFillColor(area().GetFillColor());
	fGraph->SetFillStyle(area().GetFillStyle());
}

void MgntROOT::Graph::setLine(MgntROOT::Line line) {
	if (!exist()) return;
	fGraph->SetLineColor(line().GetLineColor());
	fGraph->SetLineStyle(line().GetLineStyle());
	fGraph->SetLineWidth(line().GetLineWidth());
}

void MgntROOT::Graph::setMarker(MgntROOT::Marker marker) {
	if (!exist()) return;
	fGraph->SetMarkerColor(marker().GetMarkerColor());
	fGraph->SetMarkerStyle(marker().GetMarkerStyle());
	fGraph->SetMarkerSize (marker().GetMarkerSize());
}
		
void MgntROOT::Graph::setStyle(MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker) {
	if (!exist()) return;
	setArea(area);
	setLine(line);
	setMarker(marker);
}
		
void MgntROOT::Graph::draw(Option_t * option,  MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker) {
	if (!exist()) return;
	setStyle(area, line, marker);
	fGraph->GetHistogram()->SetLineColor(0);
	fGraph->Draw(option);
}

Bool_t MgntROOT::Graph::isGraphExist(const std::string& name) { 
	if (MgntROOT::Graph::Exist(name)) { 
		MgntSys::Error(Form("<< MgntROOT::Graph >> %s is already exist.", name.c_str())); 
		return true; 
	}
	else return false;
}


/****************/
/****  Func  ****/
/****************/
MgntROOT::Func::FitResult::FitResult(TFitResultPtr& vtrptr) : FitResult() {
	TFitResult * ptr = vtrptr.Get();
	if (ptr == nullptr || ptr->IsEmpty() || !ptr->IsValid() || ptr->NPar() == 0) return;
	params.resize(ptr->NPar());
	covmtx.resize(ptr->NPar(), std::vector<Double_t>(ptr->NPar()));
	for (Int_t i = 0; i < ptr->NPar(); ++i) {
		params.at(i).name  = ptr->ParName(i);
		params.at(i).fixed = ptr->IsParameterFixed(i);
		params.at(i).value = ptr->Value(i);
		params.at(i).error = ptr->Error(i);
		if (ptr->IsParameterBound(i))
			ptr->ParameterBounds(i, params.at(i).lmtlw, params.at(i).lmtup);
		for (Int_t j = 0; j < ptr->NPar(); ++j) {
			covmtx.at(i).at(j) = ptr->CovMatrix(i, j);
		}
	}
	ndf   = ptr->Ndf();
	chisq = ptr->Chi2();
	prob  = ptr->Prob();
	CI95  = ptr->GetConfidenceIntervals(0.95);
	CI68  = ptr->GetConfidenceIntervals(0.68);
	valid = true;
}

MgntROOT::Func * MgntROOT::Func::Head(const std::string& name) {
	if (name == "") return nullptr;
	std::map<std::string, Func*>::iterator search = fFuncMap.find(name);
	MgntROOT::Func * Func = ((search != fFuncMap.end()) ? search->second : nullptr);
	return Func;
}

Bool_t MgntROOT::Func::Exist(const std::string& name) {
	return (Head(name) != nullptr);
}

void MgntROOT::Func::Delete(const std::string& name) {
	Func * func = Func::Head(name);
	if (func != nullptr) { func->clear(); delete func; }
}

void MgntROOT::Func::Delete(Func * func) {
	if (func != nullptr) { func->clear(); delete func; }
}

void MgntROOT::Func::Delete(std::vector<Func *>& funcvec) {
	for (auto* func : funcvec)
		Func::Delete(func);
}

void MgntROOT::Func::Load(const std::string& path) {
	if (!MgntSys::TestFile(path, 'f')) {
		MgntSys::Error(Form("<< MgntROOT::Func >>  %s is not open.", path.c_str()));
		return;
	}

	TFile * file = TFile::Open(path.c_str(), "READ");
	if (file != nullptr && file->IsOpen()) {
		TList * list = file->GetListOfKeys();
		TIter next((TList*)list);
		while (TObject * key = next()) {
			TString keyName = key->GetName();
			TF1 * obj = dynamic_cast<TF1*>(file->Get(keyName.Data()));
			if (obj == nullptr) continue;
			MgntROOT::Func * Func = new MgntROOT::Func(obj);
		}
		file->Close();
	}
}

void MgntROOT::Func::Write() {
	for (auto const& it : fFuncMap) (it.second)->write();
}

void MgntROOT::Func::Save(const std::string& filename, const std::string& filepath, const std::string& dirname) {
	std::string fullpath = Form("%s/%s.root", filepath.c_str(), filename.c_str());
	if (!MgntSys::TestFile(filepath, 'd')) return;
	TDirectory * gdir = gDirectory;
	TFile * file = new TFile(fullpath.c_str(), "RECREATE");
	file->cd();
	TDirectoryFile * dir = nullptr;
	if (dirname != "") {
		dir = new TDirectoryFile(dirname.c_str(), dirname.c_str());
		dir->cd();
	}
	
	MgntROOT::Func::Write();

	file->Write();
	file->Close();
}

void MgntROOT::Func::Print(std::ostream& out) {
	out << Form("***********************************  Func List  ***********************************\n");
	out << Form("%-8s   %-18s   %-50s\n", "UniqueID", "Name", "Title");
	for (auto const& it : fFuncMap) { 
		std::string&& title = (it.second)->title();
		const TF1 * func = (*(it.second))();
		std::cout << Form("%08u   %-18s   %-50s\n", func->GetUniqueID(), func->GetName(), title.c_str());
	}
	out << Form("************************************************************************************\n");
}

MgntROOT::Func * MgntROOT::Func::New(TF1 * func, const std::string& title) {
	if (func == nullptr) return nullptr;
	new Func(func, title);
	return Func::Head(func->GetName());
}

MgntROOT::Func * MgntROOT::Func::New(const std::string& name, const std::string& title, const std::string& formula, const Bound& bound, const LIST& list) {
	new Func(name, title, formula, bound, list);
	return Func::Head(name);
}

MgntROOT::Func * MgntROOT::Func::New(const std::string& name, const std::string& title, const std::string& formula, const Bound& bound, const VLIST& vlist) {
	new Func(name, title, formula, bound, vlist);
	return Func::Head(name);
}

MgntROOT::Func::Func(TF1 * func, const std::string& title) : Func() {
	if (func == nullptr) {
		MgntSys::Error("<< MgntROOT::Func >>  func is not exist.");
		return;
	}
	if (isFuncExist(func->GetName())) return;
	fFunc  = (TF1 *) func->Clone(func->GetName());
	fName  = fFunc->GetName();
	fTitle = title;
	setStyle();
	pushBack();
}

MgntROOT::Func::Func(const std::string& name, const std::string& title, const std::string& formula, const Bound& bound, const LIST& list) : Func() {
	TF1 * func = nullptr;
	Bool_t hasBound = (MgntNum::Compare(bound.lw, bound.up) < 0);
	if (hasBound)
		func = new TF1(name.c_str(), formula.c_str(), bound.lw, bound.up);
	else
		func = new TF1(name.c_str(), formula.c_str());
	func->SetNpx(4000);

	fName  = name;
	fTitle = title;
	fFunc  = func;
	setParam(list);
	setStyle();
	pushBack();
}
		
MgntROOT::Func::Func(const std::string& name, const std::string& title, const std::string& formula, const Bound& bound, const VLIST& vlist) : Func() {
	Func::LIST list(vlist.size());
	for (Int_t it = 0; it < vlist.size(); ++it)
		list.at(it).value = vlist.at(it);
	
	TF1 * func = nullptr;
	Bool_t hasBound = (MgntNum::Compare(bound.lw, bound.up) < 0);
	if (hasBound)
		func = new TF1(name.c_str(), formula.c_str(), bound.lw, bound.up);
	else
		func = new TF1(name.c_str(), formula.c_str());
	func->SetNpx(4000);

	fName  = name;
	fTitle = title;
	fFunc  = func;
	setParam(list);
	setStyle();
	pushBack();
}
		
void MgntROOT::Func::setParam(UInt_t ipar, const Param& par) {
	if (!exist() || ipar >= fFunc->GetNpar()) return;
	fFunc->SetParameter(ipar, par.value);
	if (MgntNum::Compare(par.lmtlw, par.lmtup) < 0)
		fFunc->SetParLimits(ipar, par.lmtlw, par.lmtup);
	if (par.name != "") 
		fFunc->SetParName(ipar, par.name.c_str());
	if (par.fixed) 
		fFunc->FixParameter(ipar, par.value);
}

void MgntROOT::Func::setParam(UInt_t ipar, Double_t vpar, Bool_t fixed) {
	Func::Param par(vpar, fixed);
	setParam(ipar, par);
}

void MgntROOT::Func::setParam(const LIST& list) {
	if (!exist() || list.size() == 0) return;
	UInt_t ipar = 0;
	for (const Func::Param& par : list) {
		setParam(ipar, par);
		ipar++;
	}
}

void MgntROOT::Func::setParam(const VLIST& vlist) {
	Func::LIST list(vlist.size());
	for (Int_t it = 0; it < vlist.size(); ++it)
		list.at(it).value = vlist.at(it);
	setParam(list);
}
		
MgntROOT::Func::Param MgntROOT::Func::getParam(UInt_t ipar) {
	if (!exist() || ipar >= fFunc->GetNpar()) return Param();
	Func::Param par;
	par.name = fFunc->GetParName(ipar);
	par.value = fFunc->GetParameter(ipar);
	par.error = fFunc->GetParError(ipar);
	fFunc->GetParLimits(ipar, par.lmtlw, par.lmtup);
	par.fixed = false;
	return par;
}

MgntROOT::Func::FitResult MgntROOT::Func::fitTo(Hist * hist, Option_t * option, Option_t * goption, Double_t xmin, Double_t xmax) {
	if (!exist() || hist == nullptr || hist->ndim() != 1 || (*hist)()->GetEntries() == 0) return Func::FitResult();
	if (MgntNum::Compare(xmin, xmax) < 0) 
		(*this)()->SetRange(xmin, xmax);
	else {
		Double_t min = (*hist)()->GetBinLowEdge(1);
		Double_t max = (*hist)()->GetBinLowEdge((*hist)()->GetNbinsX());
		(*this)()->SetRange(min, max);
	}
	TFitResultPtr rltptr = (*hist)()->Fit((*this)(), CStrFmt("%sS", option), goption, xmin, xmax);
	return Func::FitResult(rltptr);
}

MgntROOT::Func::FitResult MgntROOT::Func::fitTo(Graph * graph, Option_t * option, Option_t * goption, Double_t xmin, Double_t xmax) {
	if (!exist() || graph == nullptr || (*graph)()->GetN() == 0) return Func::FitResult();
	if (MgntNum::Compare(xmin, xmax) < 0)
		(*this)()->SetRange(xmin, xmax);
	else {
		Int_t      num  = (*graph)()->GetN();
		Double_t * xarr = (*graph)()->GetX();
		Double_t min = *std::min_element(xarr, xarr + num);
		Double_t max = *std::max_element(xarr, xarr + num);
		(*this)()->SetRange(min, max);
	}
	TFitResultPtr rltptr = (*graph)()->Fit((*this)(), CStrFmt("%sS", option), goption, xmin, xmax);
	return Func::FitResult(rltptr);
}

void MgntROOT::Func::setArea(MgntROOT::Area area) {
	if (!exist()) return;
	fFunc->SetFillColor(area().GetFillColor());
	fFunc->SetFillStyle(area().GetFillStyle());
}

void MgntROOT::Func::setLine(MgntROOT::Line line) {
	if (!exist()) return;
	fFunc->SetLineColor(line().GetLineColor());
	fFunc->SetLineStyle(line().GetLineStyle());
	fFunc->SetLineWidth(line().GetLineWidth());
}

void MgntROOT::Func::setMarker(MgntROOT::Marker marker) {
	if (!exist()) return;
	fFunc->SetMarkerColor(marker().GetMarkerColor());
	fFunc->SetMarkerStyle(marker().GetMarkerStyle());
	fFunc->SetMarkerSize (marker().GetMarkerSize());
}
		
void MgntROOT::Func::setStyle(MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker) {
	if (!exist()) return;
	setArea(area);
	setLine(line);
	setMarker(marker);
}
		
void MgntROOT::Func::draw(Option_t * option,  MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker) {
	if (!exist()) return;
	setStyle(area, line, marker);
	fFunc->Draw(option);
}

Bool_t MgntROOT::Func::isFuncExist(const std::string& name) { 
	if (MgntROOT::Func::Exist(name)) { 
		MgntSys::Error(Form("<< MgntROOT::Func >> %s is already exist.", name.c_str())); 
		return true; 
	}
	else return false;
}

#endif // __MgntROOT_C__
