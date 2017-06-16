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


namespace MGROOT {
namespace Fit {

class RooPar {
    public :
        RooPar() : fExist(false), fChi(0), fNdf(0) {}
        ~RooPar() { clear(); }

        inline const Bool_t&   exist()        const { return fExist; }
        inline Int_t           num()          const { return fVal.size(); }
        
        inline const Double_t& val(Int_t idx) const { return fVal.at(idx); }
        inline const Double_t& err(Int_t idx) const { return fErr.at(idx); }
        inline const Double_t& sys(Int_t idx) const { return fSys.at(idx); }
        inline const Double_t& chi()          const { return fChi; }
        inline const Int_t&    ndf()          const { return fNdf; }

        inline void push(Double_t val, Double_t err = 0., Double_t sys = 0.) { fExist = true; fVal.push_back(val); fErr.push_back(err); fSys.push_back(sys); }

        inline void set_ndf_and_chi(Int_t ndf, Double_t chi) { fNdf = ndf; fChi = chi; fExist = true; }

        inline void clear() { fExist = false; fVal.clear(); fErr.clear(); fSys.clear(); fChi = 0; fNdf = 0; }

    protected :
        Bool_t                fExist;
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


class RooResult {
    public :
        RooResult(const RooVar& var, Bool_t extended = true, Bool_t fluc = false);
        ~RooResult() { fExist = false; fVar.clear(); fPar.clear(); }

        inline const Bool_t& exist() const { return fExist; }
        inline const RooVar& var()   const { return fVar; }
        inline const RooPar& par()   const { return fPar; }

    protected :
        Bool_t fExist;
        RooVar fVar;
        RooPar fPar;

    protected :
        static Long64_t fCounter;
};
    
Long64_t RooResult::fCounter = 0;


class RooSysResult {
    public :
        RooSysResult(const RooVar& var, Bool_t extended = true, Int_t ntimes = 400);
        ~RooSysResult() { fExist = false; fVar.clear(); fStdPar.clear(); fSysPar.clear(); fSysFitSet.clear(); }

        inline const Bool_t& exist()   const { return fExist; }
        inline const RooVar& var()     const { return fVar; }
        inline const RooPar& std_par() const { return fStdPar; }
        inline const RooPar& sys_par() const { return fSysPar; }

    protected :
        Bool_t fExist;
        RooVar fVar;
        RooPar fStdPar;
        RooPar fSysPar;

    private :
        std::vector<RooPar> fSysFitSet;
};
    
} // namespace Fit
} // namespace MGROOT

#endif // __ROOTLibs_MGFit_H__
