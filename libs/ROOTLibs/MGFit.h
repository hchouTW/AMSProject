#ifndef __ROOTLibs_MGFit_H__
#define __ROOTLibs_MGFit_H__

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

//using namespace RooFit;

namespace MGROOT {
namespace Fit {

class RooPar {
    public :
        RooPar() : exist_(false), fChi(0), fNdf(0) {}
        ~RooPar() { clear(); }

        inline const Bool_t&   exist()        const { return exist_; }
        inline Int_t           num()          const { return fVal.size(); }
        
        inline const Double_t& val(Int_t idx) const { return fVal.at(idx); }
        inline const Double_t& err(Int_t idx) const { return fErr.at(idx); }
        inline const Double_t& sys(Int_t idx) const { return fSys.at(idx); }
        inline const Double_t& chi()          const { return fChi; }
        inline const Int_t&    ndf()          const { return fNdf; }

        inline void push(Double_t val, Double_t err = 0., Double_t sys = 0.) { exist_ = true; fVal.push_back(val); fErr.push_back(err); fSys.push_back(sys); }

        inline void set_ndf_and_chi(Int_t ndf, Double_t chi) { fNdf = ndf; fChi = chi; exist_ = true; }

        inline void clear() { exist_ = false; fVal.clear(); fErr.clear(); fSys.clear(); fChi = 0; fNdf = 0; }

    protected :
        Bool_t                exist_;
        std::vector<Double_t> fVal;
        std::vector<Double_t> fErr;
        std::vector<Double_t> fSys;
        Double_t              fChi;
        Int_t                 fNdf;
        
    public :
        struct RooParSort {
            Bool_t operator() (const RooPar& par1, const RooPar& par2) {
                if (par1.chi() < par2.chi()) return true;
                else return false;
            }
        };
};


class RooVar {
    public :
        RooVar() : fSamp(nullptr), fSumt(nullptr) { clear(); }
        ~RooVar() { clear(); }

        RooVar(const std::string& name, Hist * samp, Hist * sumt, const HistList& temp, Bool_t link = true, Double_t min = 0., Double_t max = 0.) { set(name, samp, sumt, temp, min, max); } 
        RooVar(const std::string& name, Hist * samp, const HistList& temp, Bool_t link = true, Double_t min = 0., Double_t max = 0.) { set(name, samp, nullptr, temp, min, max); } 

        inline Bool_t exist() const { return (fName != "" && fSamp != nullptr && fTemp.size() != 0); }
        inline Int_t    num() const { return fTemp.size(); }

        // TODO (hchou): add const to each function
        inline const std::string& name() const { return fName; }
        inline  Hist *       samp() const { return fSamp; }
        inline  Hist *       sumt() const { return fSumt; }
        inline const HistList&    temp() const { return fTemp; }
        inline  Hist *       temp(Int_t idx) const { return fTemp.at(idx); }
        inline const Double_t&    min() const { return fMin; }
        inline const Double_t&    max() const { return fMax; }
        inline const Int_t&       min_bin() const { return fMinBin; }
        inline const Int_t&       max_bin() const { return fMaxBin; }

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
        Int_t        fMinBin;
        Int_t        fMaxBin;
};


class RooResult {
    public :
        RooResult(const RooVar& var, Bool_t extended = true, Bool_t fluc = false);
        ~RooResult() { exist_ = false; var_.clear(); par_.clear(); }

        inline const Bool_t& exist() const { return exist_; }
        inline const RooVar& var()   const { return var_; }
        inline const RooPar& par()   const { return par_; }

    protected :
        Bool_t exist_;
        RooVar var_;
        RooPar par_;

    protected :
        static Long64_t fCounter;
};
    
Long64_t RooResult::fCounter = 0;


class RooSysResult {
    public :
        RooSysResult(const RooVar& var, Bool_t extended = true, Int_t ntimes = 400);
        ~RooSysResult() { exist_ = false; var_.clear(); std_par_.clear(); sys_par_.clear(); sys_fit_set_.clear(); }

        inline const Bool_t& exist()   const { return exist_; }
        inline const RooVar& var()     const { return var_; }
        inline const RooPar& std_par() const { return std_par_; }
        inline const RooPar& sys_par() const { return sys_par_; }

    protected :
        Bool_t exist_;
        RooVar var_;
        RooPar std_par_;
        RooPar sys_par_;

    private :
        std::vector<RooPar> sys_fit_set_;
};
    
} // namespace Fit
} // namespace MGROOT

#endif // __ROOTLibs_MGFit_H__
