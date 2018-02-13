#ifndef __ROOTLibs_MGAxis_C__
#define __ROOTLibs_MGAxis_C__

namespace MGROOT {

/******************/
/****  BinSet  ****/
/******************/
BinSet::BinSet(const BinSet& binset, UInt_t merge, Bool_t invert) {
	*this = binset;
	if (merge > 1) 
		bins_ = merge_bins(bins_, merge);
    if (invert) invert_bins();
}


BinSet::BinSet(const std::vector<Double_t>& bins, UInt_t merge, Bool_t invert) {
	bins_ = merge_bins(bins, merge);
    if (invert) invert_bins();
}


BinSet::BinSet(const UInt_t nbin, const Double_t * bins, UInt_t merge, Bool_t invert) {
	std::vector<Double_t> vbins(bins, bins + nbin + 1);
	bins_ = merge_bins(vbins, merge);
    if (invert) invert_bins();
}


Bool_t BinSet::check_sequence(const std::vector<Double_t>& list) {
	if (list.size() == 0) return false;
	for (auto&& elem : list)
		if (!MGNumc::Valid(elem))
			return false;
	for (UInt_t it = 1; it < list.size(); ++it)
		if (MGNumc::Compare(list.at(it), list.at(it-1)) <= 0)
			return false;
	return true;
}


Int_t BinSet::find(Double_t val) const {
	if (!exist()) return -1;
	Bool_t bdl = (MGNumc::Compare(val, min()) < 0);
	Bool_t bdu = (MGNumc::Compare(val, max()) >= 0);
	if (bdl) return 0;
	if (bdu) return bins_.size();

	Int_t range[2] = { 0, nbin() };
	while ((range[1] - range[0]) != 1) {
		Int_t mid = (range[0] + range[1]) /2;
		if      (MGNumc::Compare(val, bins_.at(mid))  < 0) range[1] = mid;
		else if (MGNumc::Compare(val, bins_.at(mid)) >= 0) range[0] = mid;
	}
	Int_t bin = range[1];
	return bin;
}


std::vector<Double_t> BinSet::merge_bins(const std::vector<Double_t>& bins, UInt_t merge) {
	std::vector<Double_t> vbins;
	if (merge < 1 || bins.size() == 0 || ((bins.size()-1)%merge)!=0) return vbins;
	if (!check_sequence(bins)) return vbins;
	if (merge == 1) vbins = bins;
	else {
		Int_t cnt = 0;
		Int_t num = (bins.size()-1)/merge + 1;
		vbins.resize(num);
		for (auto&& elem : vbins) {
			elem = bins.at(cnt);
			cnt += merge;
		}
	}
	return vbins;
}


Bool_t BinSet::invert_bins() {
    Bool_t result = (bins_.size() != 0);
	for (UInt_t idx = 0; result && idx < bins_.size(); ++idx)
		if (MGNumc::Compare(bins_.at(idx)) <= 0) result = false;

    if (result) {
		std::vector<Double_t> list(2*bins_.size()+1);
		list.at(bins_.size()) = 0.0;
		UInt_t nbin = list.size() - 1;
		for (UInt_t it = 0; it < bins_.size(); ++it) {
			Double_t inv = 1. / bins_.at(it);
			list.at(it)                                                = -1. * inv;
			list.at(static_cast<Int_t>(nbin) - static_cast<Int_t>(it)) = inv;
		}
        bins_ = list;
        return true;
	}
    else {
        bins_.clear();
	    MGSys::ShowError("<< BinSet >>  Invert bins failure.");
        return false;
    }
}


/****************/
/****  Axis  ****/
/****************/
Axis::Axis(const Axis& axis, UInt_t mergeFT, Bool_t invert) {
	binset_ = BinSet(axis(), mergeFT, invert);
	if (binset_().size() != 0)
		set_name_and_title(VAR_NAME(this), axis.title());
}


Axis::Axis(const std::string& title, const std::vector<Double_t>& list, UInt_t mergeFT, Bool_t invert) {
	binset_ = BinSet(list, mergeFT, invert);
	if (binset_().size() != 0)
		set_name_and_title(VAR_NAME(this), title);
}


Axis::Axis(const std::string& title, const Axis& axis, UInt_t mergeFT, Bool_t invert) {
	binset_ = BinSet(axis(), mergeFT, invert);
	if (binset_().size() != 0)
		set_name_and_title(VAR_NAME(this), title);
}


Axis::Axis(const std::string& title, UInt_t nbin, Double_t lw, Double_t up, AxisScale scl) {
	if (nbin < 1) return;
	if (MGNumc::Compare(lw, up) >= 0) return;
	if (scl == AxisScale::kLog && (MGNumc::Compare(lw) <= 0 || MGNumc::Compare(up) <= 0)) return;
	
	std::vector<Double_t> list(nbin+1);
	switch (scl) {
		case AxisScale::kLinear :
		{
			Double_t step = (up - lw) / Double_t(nbin);UInt_t cnt = 0;
			for (Double_t& elem : list) {
				elem = (lw + cnt * step);
				cnt++;
			}
			break;
		}
		case AxisScale::kLog :
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
	binset_ = BinSet(list);
	
	if (binset_().size() != 0)
		set_name_and_title(VAR_NAME(this), title);
}


Bool_t Axis::init_TAxis(const TAxis * axis, UInt_t mergeFT, Bool_t copyNT) {
	if (axis == nullptr) return false;
	UInt_t nbin = axis->GetNbins();
	std::vector<Double_t> list(nbin+1);
	for (UInt_t ibin = 1; ibin <= nbin+1; ++ibin)
		list.at(ibin-1) = axis->GetBinLowEdge(ibin);
	binset_ = BinSet(list, mergeFT);	
	if (copyNT) set_name_and_title(axis);
	return true;
}


Bool_t Axis::init_TH1(const TH1 * hist, AxisDim dim, UInt_t mergeFT, Bool_t copyNT) {
	if (hist == nullptr) return false;
	Int_t ndim = hist->GetDimension();
	if (static_cast<Int_t>(dim) > ndim) return false;
	const TAxis * axis = nullptr;
	switch (dim) {
		case AxisDim::kX : axis = hist->GetXaxis(); break;
		case AxisDim::kY : axis = hist->GetYaxis(); break;
		case AxisDim::kZ : axis = hist->GetZaxis(); break;
		default : break;
	}
	return init_TAxis(axis, mergeFT, copyNT);
}


Bool_t Axis::init_TObject(const TObject * obj, AxisDim dim, UInt_t mergeFT, Bool_t copyNT) {
	if (obj == nullptr) return false;
	if (!dynamic_cast<const TH1*>(obj)) return false;
	TH1 * hist = (TH1*)obj;
	return init_TH1(hist, dim, mergeFT, copyNT);
}


Double_t Axis::center(UInt_t ibin, AxisScale scl) const {
	if (!exist()) return 0.0;
	if      (ibin ==                             0) return binset_.min();
	else if (ibin == static_cast<UInt_t>(nbin()+1)) return binset_.max();
	else {
		switch (scl) {
			case AxisScale::kLinear : 
				return (0.5 * (binset_(ibin-1) + binset_(ibin))); break;
			case AxisScale::kLog : 
				return (std::sqrt(binset_(ibin-1) * binset_(ibin))); break;
			default : break;
		}
	}
	return 0.0;
}


Double_t Axis::center(UInt_t ibin, Double_t power) const {
	if (!exist()) return 0.0;
	if      (ibin ==                             0) return binset_.min();
	else if (ibin == static_cast<UInt_t>(nbin()+1)) return binset_.max();
    else if (MGNumc::Compare(power)       == 0) return center(ibin, AxisScale::kLinear);
    else if (MGNumc::Compare(power, -1.0) == 0) return center(ibin, AxisScale::kLog);
	else {
		Double_t lw = binset_(ibin-1);
        Double_t up = binset_(ibin);
        if (MGNumc::Compare(lw) <= 0 || MGNumc::Compare(up) <= 0) return 0.0;
        return std::sqrt(lw * up);
    }
    return 0.0;
}


void Axis::print() const {
	if (!exist()) return;
	COUT("************************************  Axis  ************************************\n");
	COUT("Axis : NAME(\"%s\")\n", name_.c_str());
	COUT("Axis : TITLE(\"%s\")\n", title_.c_str());
	COUT("Axis : NBIN(%d)\n", nbin());
	COUT("Axis : REG[%16.8f %16.8f]\n", min(), max());
	COUT("--------------------------------------------------------------------------------\n");
	for (Int_t it = 1; it <= nbin(); ++it)
		COUT("BIN(%4d)     REG[%16.6f  %16.6f]    WIDTH(%16.6f)\n", it, binset_(it-1), binset_(it), binset_.width(it));
	COUT("********************************************************************************\n");
}


} // namespace MGROOT

#endif // __ROOTLibs_MGAxis_C__
