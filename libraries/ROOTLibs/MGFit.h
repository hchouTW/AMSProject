#ifndef __MGFit_H__
#define __MGFit_H__

#include <TFractionFitter.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooBifurGauss.h>
#include <RooMinuit.h>
#include <RooChi2Var.h>
#include <RooMsgService.h>


namespace MGROOT {
	namespace Fit {
		class RooPar;
		class RooVar;
		class RooResult;
		//class RooSysResult;
	}
}

class MGROOT::Fit::RooPar {
	public :
		RooPar() : fExist(false), fChi(0), fNdf(0) {}
		~RooPar() {}

		inline const Bool_t&   exist()        const { return fExist; }
		inline Int_t           num()          const { return fVal.size(); }
		
		inline const Double_t& val(Int_t idx) const { return fVal.at(idx); }
		inline const Double_t& err(Int_t idx) const { return fErr.at(idx); }
		inline const Double_t& chi()          const { return fChi; }
		inline const Int_t&    ndf()          const { return fNdf; }

		inline void push(Double_t val, Double_t err = 0.) 
			{ fExist = true; fVal.push_back(val); fErr.push_back(err); }

		inline void clear() 
			{ fExist = false; fVal.clear(); fErr.clear(); fChi = 0; fNdf = 0; }

	protected :
		Bool_t                fExist;
		std::vector<Double_t> fVal;
		std::vector<Double_t> fErr;
		Double_t              fChi;
		Int_t                 fNdf;
};


class MGROOT::Fit::RooVar {
	public :
		RooVar() { clear(); }
		~RooVar() { clear(); }

		RooVar(const std::string& name, Hist * samp, Hist * sumt, const HistList& temp, Bool_t link = true, Double_t min = 0., Double_t max = 0.)
			{ set(name, samp, sumt, temp, min, max); } 
		RooVar(const std::string& name, Hist * samp, const HistList& temp, Bool_t link = true, Double_t min = 0., Double_t max = 0.)
			{ set(name, samp, nullptr, temp, min, max); } 

		inline Bool_t exist() const { return (fName != "" && fSamp != nullptr && fTemp.size() != 0); }
		inline Int_t    num() const { return fTemp.size(); }

		inline const std::string& name() const { return fName; }
		inline const Hist *       samp() const { return fSamp; }
		inline const Hist *       sumt() const { return fSumt; }
		inline const HistList&    temp() const { return fTemp; }
		inline const Hist *       temp(Int_t idx) const { return fTemp.at(idx); }
		inline const Double_t&    min() const { return fMin; }
		inline const Double_t&    max() const { return fMax; }

		void set(const std::string& name, Hist * samp, Hist * sumt, const HistList& temp, Bool_t link = true, Double_t min = 0., Double_t max = 0.);
		void set(const std::string& name, Hist * samp, const HistList& temp, Bool_t link = true, Double_t min = 0., Double_t max = 0.);
		void clear();

	protected :
		Bool_t       fLink;
		std::string  fName;
		Hist *       fSamp;
		Hist *       fSumt;
		HistList     fTemp;
		Double_t     fMin;
		Double_t     fMax;
};


class MGROOT::Fit::RooResult {
	public :
		RooResult(const RooVar& var, Bool_t extended = true, Bool_t fluc = false);
		~RooResult();

		inline const Bool_t& exist() const { return fExist; }
		inline const RooVar& var()   const { return fVar; }
		inline const RooPar& par()   const { return fPar; }

	protected :
		Bool_t fExist;
		RooVar fVar;
		RooPar fPar;

	protected :
		static Long64_t COUNT;
};
	
Long64_t MGROOT::Fit::RooResult::COUNT = 0;





/*
class MGROOT::Fit::RooSysResult {
	public :
		class PARAM {
			public :
				PARAM() : chi(0.0) {}
				PARAM(MGROOT::Fit::RooResult& rlt)
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
		RooSysResult(MGROOT::Fit::RooParam & param, Bool_t extended = true, Int_t ntimes = 400);
		~RooSysResult();

		MGROOT::Fit::RooResult& operator()() { return fRooResult; }

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
		MGROOT::Fit::RooResult fRooResult;

		Bool_t                                     fSysValid;
		std::vector<std::pair<Double_t, Double_t>> fSysParam;
		std::vector<Double_t>                      fSysErr;
		Double_t                                   fSysChisq;
		Int_t                                      fSysNdf;
		
		std::vector<RooSysResult::PARAM> fSysFitSet;
};
*/

#endif // __MGFit_H__
