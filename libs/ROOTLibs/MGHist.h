#ifndef __ROOTLibs_MGHist_H__
#define __ROOTLibs_MGHist_H__

#include <TMath.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TList.h>
#include <TKey.h>
#include <TObjArray.h>
#include <THStack.h>

#include <unordered_map>

/******************/
/****  MGROOT  ****/
/******************/
namespace MGROOT {

enum class HistType : Int_t { kNone, kHist, kProfile };
enum class HistDim  : Int_t { kNone = 0, k1D = 1, k2D = 2, k3D = 3 };
enum class HistProj : Int_t { kX = 1, kY = 2, kZ = 3, kYZ = 4, kZY = 7, kXZ = 5, kZX = 8, kXY = 6, kYX = 9 };
enum class HistNorm : Int_t { kEntries, kIntegral, kArea };
enum class HistArith : Int_t { kAddition, kSubtract, kMultiply, kDivide };


class HistBin {
	public :
		HistBin() : dim_(HistDim::kNone), ibin_(-1), jbin_(-1), kbin_(-1) {}
		HistBin(Int_t ibin) : HistBin() { dim_ = HistDim::k1D; ibin_ = ibin; }
		HistBin(Int_t ibin, Int_t jbin) : HistBin() { dim_ = HistDim::k2D; ibin_ = ibin; jbin_ = jbin; }
		HistBin(Int_t ibin, Int_t jbin, Int_t kbin) : HistBin() { dim_ = HistDim::k3D; ibin_ = ibin; jbin_ = jbin; kbin_ = kbin; }
		~HistBin() {}

		inline const HistDim& dim() const { return dim_; }
		inline const Int_t& x() const { return ibin_; }
		inline const Int_t& y() const { return jbin_; }
		inline const Int_t& z() const { return kbin_; }

	protected :
		HistDim dim_;
		Int_t   ibin_;
		Int_t   jbin_;
		Int_t   kbin_;
};


class HistAxis {
	public :
		HistAxis() : dim_(HistDim::kNone) {}
		HistAxis(const Axis& xaxis, const std::string& ytitle = "") : HistAxis() { if (!xaxis.exist()) return; dim_ = HistDim::k1D; xaxis_ = xaxis; yaxis_ = Axis(ytitle); }
		HistAxis(const Axis& xaxis, const Axis& yaxis, const std::string& ztitle = "") : HistAxis() { if (!xaxis.exist() || !yaxis.exist()) return; dim_ = HistDim::k2D; xaxis_ = xaxis; yaxis_ = yaxis; zaxis_ = Axis(ztitle); }
		HistAxis(const Axis& xaxis, const Axis& yaxis, const Axis& zaxis) : HistAxis() { if (!xaxis.exist() || !yaxis.exist() || !zaxis.exist()) return; dim_ = HistDim::k3D; xaxis_ = xaxis; yaxis_ = yaxis; zaxis_ = zaxis; }
		~HistAxis() {}
		
        HistAxis(const std::string& xtitle) : HistAxis() { xaxis_ = Axis(xtitle); }
        HistAxis(const std::string& xtitle, const std::string& ytitle) : HistAxis() { xaxis_ = Axis(xtitle); yaxis_ = Axis(ytitle); }
        HistAxis(const std::string& xtitle, const std::string& ytitle, const std::string& ztitle) : HistAxis() { xaxis_ = Axis(xtitle); yaxis_ = Axis(ytitle); zaxis_ = Axis(ztitle); }

		inline const HistDim& dim() const { return dim_; }
		
		inline const Axis& x() const { return xaxis_; }
		inline const Axis& y() const { return yaxis_; }
		inline const Axis& z() const { return zaxis_; }

		inline const Double_t& x(Int_t idx) const { return xaxis_(idx); }
		inline const Double_t& y(Int_t idx) const { return yaxis_(idx); }
		inline const Double_t& z(Int_t idx) const { return zaxis_(idx); }

	protected :
		HistDim dim_;
		Axis xaxis_;
		Axis yaxis_;
		Axis zaxis_;
};


class Hist {
    public :
        static void LoadDefaultEnvironment() { COUT("MGROOT::Hist : Load default environment.\n"); TH1::SetDefaultSumw2(true); TH1::AddDirectory(false); }
        static void AddDirectory(Bool_t opt = true) { TH1::AddDirectory(opt); }

	protected :
		Hist() { unique_id_ = -1; info_.first = HistType::kNone; info_.second = HistDim::kNone; hist_ = nullptr; }

	public :
		Hist(const TH1 * hist, Bool_t reset = false) : Hist("", "", hist, reset) {}
		Hist(const std::string& name, const TH1 * hist, Bool_t reset = false) : Hist(name, "", hist, reset) {}
		Hist(const std::string& name, const std::string& title, const TH1 * hist, Bool_t reset = false);
        Hist(const std::string& name, const std::string& title, const HistAxis& axis, HistType type = HistType::kHist);
        Hist(const std::string& name, const HistAxis& axis, HistType type = HistType::kHist) : Hist(name, "", axis, type) {}
		~Hist() { clear(); }

		inline Bool_t exist() const { return (hist_ != nullptr && info_.first != HistType::kNone && info_.second != HistDim::kNone); }
		inline Bool_t exist(const HistType& type) const { return (hist_ != nullptr && info_.first == type); }
		inline Bool_t exist(const HistType& type, const HistDim& dim) const { return (hist_ != nullptr && info_.first == type && info_.second == dim); }
		
		inline const std::string& name()       const { return name_; }
		inline const std::string& title()      const { return title_; }
		inline const HistType&    type()       const { return info_.first; }
		inline const HistDim&     dim()        const { return info_.second; }
		inline const HistAxis&    axis()       const { return axis_; }
		inline TH1 *              operator()() const { return hist_; }
		
        inline const Axis&        xaxis()      const { return axis_.x(); }
		inline const Axis&        yaxis()      const { return axis_.y(); }
		inline const Axis&        zaxis()      const { return axis_.z(); }
		
		void                          set_bin(const HistBin& bin, Double_t content, Double_t error = -1);
		std::pair<Double_t, Double_t> get_bin(const HistBin& bin) const;

		void scale(Double_t scl = 1., Option_t * option = "");
		void normalized(const HistNorm& norm = HistNorm::kEntries);
		void style(const TAttLine& line = Line(), const TAttMarker& marker = Marker(), const TAttFill& fill = Fill());

		void fill(Double_t a, Double_t b = 1.0, Double_t c = 1.0, Double_t d = 1.0, Double_t e = 1.0);
        
		void fillH(Double_t a, Double_t b = 1.0, Double_t c = 1.0, Double_t d = 1.0);
        inline void fillH1D(Double_t a, Double_t b = 1.0) { if (exist(HistType::kHist, HistDim::k1D)) dynamic_cast<TH1D*>(hist_)->Fill(a, b); }
        inline void fillH2D(Double_t a, Double_t b, Double_t c = 1.0) { if (exist(HistType::kHist, HistDim::k2D)) dynamic_cast<TH2D*>(hist_)->Fill(a, b, c); }
        inline void fillH3D(Double_t a, Double_t b, Double_t c, Double_t d = 1.0) { if (exist(HistType::kHist, HistDim::k3D)) dynamic_cast<TH3D*>(hist_)->Fill(a, b, c, d); }
        
		void fillP(Double_t a, Double_t b = 1.0, Double_t c = 1.0, Double_t d = 1.0, Double_t e = 1.0);
        inline void fillP1D(Double_t a, Double_t b, Double_t c = 1.0) { if (exist(HistType::kProfile, HistDim::k1D)) dynamic_cast<TProfile*>(hist_)->Fill(a, b, c); }
        inline void fillP2D(Double_t a, Double_t b, Double_t c, Double_t d = 1.0) { if (exist(HistType::kProfile, HistDim::k2D)) dynamic_cast<TProfile2D*>(hist_)->Fill(a, b, c, d); }
        inline void fillP3D(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e = 1.0) { if (exist(HistType::kProfile, HistDim::k3D)) dynamic_cast<TProfile3D*>(hist_)->Fill(a, b, c, d, e); }

		void draw(Option_t * option = "") { if (exist()) hist_->Draw(option); }
		void write(const std::string& name = "") { if (exist()) hist_->Write(name.c_str(), 0, 0); }
		
		Hist(const Hist&) = delete;
		Hist& operator= (const Hist&) = delete;

	protected :
		void push();
		void clear();

	protected :
		Long64_t                     unique_id_;
		std::string                  name_;
		std::string                  title_;
		std::pair<HistType, HistDim> info_;
		TH1*                         hist_;
		HistAxis                     axis_;

	public :
		static Hist * Head(const std::string& name, Bool_t show = true);
		
		static Hist * New(const TH1 * hist, Bool_t reset = false) { return (new Hist("", "", hist, reset)); }
		static Hist * New(const std::string& name, const TH1 * hist, Bool_t reset = false) { return (new Hist(name, "", hist, reset)); }
		static Hist * New(const std::string& name, const std::string& title, const TH1 * hist, Bool_t reset = false) { return (new Hist(name, title, hist, reset)); }
        static Hist * New(const std::string& name, const std::string& title, const HistAxis& axis, HistType type = HistType::kHist) { return (new Hist(name, title, axis, type)); }
        static Hist * New(const std::string& name, const HistAxis& axis, HistType type = HistType::kHist) { return (new Hist(name, axis, type)); }
		
		static void Delete(Hist* hist) { if (hist != nullptr) { hist->clear(); delete hist; hist = nullptr; } }
		static void Delete(Hist& hist) { hist.clear(); }
		static void Delete(const std::string& name);
		static void Delete(std::vector<Hist*>& histvec);
		
		static void Write();
		static Bool_t Load(const std::string& filename, const std::string& filepath = ".");
		static Bool_t Save(const std::string& filename, const std::string& filepath = ".");

		static void Print();
	
	public :
		static Hist * Calculate(const HistArith& arith, const std::string& name, const std::string& title, Hist * hSon, Hist * hMom);
		static Hist * Project(const HistProj& proj, Hist * hMom, Int_t isb = -1, Int_t ieb = -1, Int_t jsb = -1, Int_t jeb = -1);
		static std::vector<Hist*> ProjectAll(const HistProj& proj, Hist * hMom);
		
		//static THStack * Collect(const std::string& name, const std::string& title, const std::vector<std::string>& list);
		//static THStack * Collect(const std::string& name, const std::string& title, const std::vector<Hist*>& list);
		
        static THStack * Collect(const std::string& name, const std::vector<std::string>& list);
		static THStack * Collect(const std::string& name, const std::vector<Hist*>& list);

	protected :
		static Long64_t                               counter_;
		static std::unordered_map<std::string, Hist*> hist_map_;
};

Long64_t                               Hist::counter_ = 0;
std::unordered_map<std::string, Hist*> Hist::hist_map_;


using HistList = std::vector<Hist*>;


} // namespace MGROOT

#endif // __ROOTLibs_MGHist_H__
