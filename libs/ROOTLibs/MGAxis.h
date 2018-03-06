#ifndef __ROOTLibs_MGAxis_H__
#define __ROOTLibs_MGAxis_H__

#include <TAxis.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>


namespace MGROOT {

enum class AxisScale : Int_t { kLinear = 0, kLog = 1 };
enum class AxisDim : Int_t { kX = 1, kY = 2, kZ = 3 };

using BinList = std::vector<Double_t>;

/******************/
/****  BinSet  ****/
/******************/
/*  MGROOT::BinSet : Example
 *  nbin = 6, bins = {0, 1, 2, 3, 4, 5, 6}
 *
 *  DEMO :
 *         -----|---|---|---|---|---|---|-----
 *  bound       0   1   2   3   4   5   6
 *               \ / \ / \ / \ / \ / \ /
 *  bin      0    1   2   3   4   5   6    7
 */

class BinSet {
	public :
		BinSet() {}
		BinSet(const BinSet& binset, UInt_t merge = 1, Bool_t invert = false);
		BinSet(const std::vector<Double_t>& bins, UInt_t merge = 1, Bool_t invert = false);
		BinSet(const UInt_t nbin, const Double_t * bins, UInt_t merge = 1, Bool_t invert = false);
		~BinSet() {}

		inline Bool_t   exist() const { return (bins_.size() != 0); }
		inline Int_t    nbin()  const { return (bins_.size()-1); }
		inline Double_t min()   const { return bins_.front(); }
		inline Double_t max()   const { return bins_.back(); }
	
		inline Double_t width(UInt_t ibin) const { return ((ibin>0&&ibin<bins_.size())?(bins_.at(ibin)-bins_.at(ibin-1)):-1); }
		inline Int_t find(Double_t val) const;

		inline const std::vector<Double_t>& operator()() const { return bins_; }
		inline const Double_t& operator()(UInt_t idx) const { return bins_.at(idx); }

	protected :
		Bool_t check_sequence(const std::vector<Double_t>& list);
		std::vector<Double_t> merge_bins(const std::vector<Double_t>& bins, UInt_t merge = 1);
        Bool_t invert_bins();

	protected :
		std::vector<Double_t> bins_;
};


/****************/
/****  Axis  ****/
/****************/
class Axis {
	public :
		Axis(const std::string& title = "") : title_(title) {}
		
        Axis(const Axis& axis, UInt_t mergeFT = 1, Bool_t invert = false);
		Axis(std::initializer_list<Double_t> list) { binset_ = BinSet(list); }
		Axis(const std::string& title, const Axis& axis, UInt_t mergeFT = 1, Bool_t invert = false);
		Axis(const std::string& title, const std::vector<Double_t>& list, UInt_t mergeFT = 1, Bool_t invert = false);
		Axis(const std::string& title, UInt_t nbin, Double_t lw, Double_t up, AxisScale scl = AxisScale::kLinear);

		Axis(const TAxis * axis, UInt_t merge = 1) { init_TAxis(axis, merge, true); }
		Axis(const TH1 * hist, AxisDim dim = AxisDim::kX, UInt_t merge = 1) { init_TH1(hist, dim, merge, true); }
		Axis(const TObject * obj, AxisDim dim = AxisDim::kX, UInt_t merge = 1) { init_TObject(obj, dim, merge, true); }
		
		Axis(const std::string& title, const TAxis * axis, UInt_t merge = 1) { if (init_TAxis(axis, merge)) set_name_and_title(VAR_NAME(this), title); }
		Axis(const std::string& title, const TH1 * hist, AxisDim dim = AxisDim::kX, UInt_t merge = 1) { if (init_TH1(hist, dim, merge)) set_name_and_title(VAR_NAME(this), title); }
		Axis(const std::string& title, const TObject * obj, AxisDim dim = AxisDim::kX, UInt_t merge = 1) { if (init_TObject(obj, dim, merge)) set_name_and_title(VAR_NAME(this), title); }

		~Axis() {}

	public :
		inline const std::string& name()  const { return name_; }
		inline const std::string& title() const { return title_; }
		inline const BinSet& operator()() const { return binset_; }     
		
		inline Bool_t   exist() const { return binset_.exist(); }
		inline Int_t    nbin()  const { return binset_.nbin(); }
		inline Double_t min()   const { return binset_.min(); }
		inline Double_t max()   const { return binset_.max(); }

		Double_t center(UInt_t ibin, AxisScale scl = AxisScale::kLinear) const;
		Double_t center(UInt_t ibin, Double_t power) const;
		inline Double_t width(UInt_t ibin) const { return binset_.width(ibin); }
		inline Int_t find(Double_t val) const { return binset_.find(val); }
		inline const Double_t& operator()(UInt_t idx) const { return binset_().at(idx); }
		
		void print() const;

	protected :
		void   set_name_and_title(const std::string& name, const std::string& title) { name_ = name; title_ = title; }
		void   set_name_and_title(const TAxis * axis) { if (axis) { name_ = axis->GetName(); title_ = axis->GetTitle(); } }
		Bool_t init_TAxis(const TAxis * axis, UInt_t merge = 1, Bool_t copyNT = false);
		Bool_t init_TH1(const TH1 * hist, AxisDim dim = AxisDim::kX, UInt_t merge = 1, Bool_t copyNT = false);
		Bool_t init_TObject(const TObject * obj, AxisDim dim = AxisDim::kX, UInt_t merge = 1, Bool_t copyNT = false);

	protected :
		std::string name_;
		std::string title_;
		BinSet      binset_;
};

} // namespace MGROOT

#endif // __ROOTLibs_MGAxis_H__
