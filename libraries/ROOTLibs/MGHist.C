#ifndef __ROOTLibs_MGHist_C__
#define __ROOTLibs_MGHist_C__

namespace MGROOT {

/****************/
/****  Hist  ****/
/****************/
Hist::Hist(const std::string& name, const std::string& title, const TH1 * hist, Bool_t reset) : Hist() {
	if (hist == nullptr) {
		MGSys::ShowError("<< Hist::Hist >> Histogram is not existence..");
		return;
	}
	else if (name != "" && Head(name, false) != nullptr) {
		MGSys::ShowError(STR_FMT("<< Hist::Hist >> Histogram (%s) is existence.", name.c_str()));
		return;	
	}
	else if (name == "" && Head(hist->GetName(), false) != nullptr) {
		MGSys::ShowError(STR_FMT("<< Hist::Hist >> Histogram (%s) is existence.", name.c_str()));
		return;
	}
	
	std::string classname = hist->ClassName();
	HistType type = (classname.find("TProfile") != std::string::npos) ? HistType::kProfile : HistType::kHist;
	
	Axis xaxis(hist, AxisDim::kX);
	Axis yaxis(hist, AxisDim::kY);
	Axis zaxis(hist, AxisDim::kZ);
	HistDim dim = static_cast<HistDim>(hist->GetDimension());
	switch (dim) {
		case HistDim::k1D : axis_ = HistAxis(xaxis); break;
		case HistDim::k2D : axis_ = HistAxis(xaxis, yaxis); break;
		case HistDim::k3D : axis_ = HistAxis(xaxis, yaxis, zaxis); break;
		default : break;
	}

	name_ =  (name == "")  ? hist->GetName() : name;
	title_ = (title == "") ? hist->GetTitle() : title;
	info_.first  = type;
	info_.second = dim;

	hist_ = (TH1*) hist->Clone(name_.c_str());
	hist_->SetTitle(title_.c_str());
	if (reset) {
		hist_->GetListOfFunctions()->Clear();
		hist_->Reset();
	}
	
	style();
	push();
}


Hist::Hist(const std::string& name, const std::string& title, const HistAxis& axis, HistType type) : Hist() {
	if (name == "") {
		MGSys::ShowError("<< Hist::Hist >> Name is empty.");
		return;
	}
	if (axis.dim() == HistDim::kNone) {
		MGSys::ShowError("<< Hist::Hist >> Axis is empty.");
		return;
	}
	if (Head(name, false) != nullptr) {
		MGSys::ShowError(STR_FMT("<< Hist::Hist >> Histogram (%s) is existence.", name.c_str()));
		return;	
	}

	switch(axis.dim()) {
		case HistDim::k1D : 
			{
				if (type == HistType::kHist)
					hist_ = new TH1D(name.c_str(), title.c_str(), axis.x().nbin(), &axis.x(0));
				else 
					hist_ = new TProfile(name.c_str(), title.c_str(), axis.x().nbin(), &axis.x(0));
				
				hist_->GetXaxis()->SetNameTitle(axis.x().name().c_str(), axis.x().title().c_str());
				break;
			}
		case HistDim::k2D :
			{
				if (type == HistType::kHist)
					hist_ = new TH2D(name.c_str(), title.c_str(), axis.x().nbin(), &axis.x(0), axis.y().nbin(), &axis.y(0));
				else
					hist_ = new TProfile2D(name.c_str(), title.c_str(), axis.x().nbin(), &axis.x(0), axis.y().nbin(), &axis.y(0));
				
				hist_->GetXaxis()->SetNameTitle(axis.x().name().c_str(), axis.x().title().c_str());
				hist_->GetYaxis()->SetNameTitle(axis.y().name().c_str(), axis.y().title().c_str());
				break;
			}
		case HistDim::k3D : 
			{
				if (type == HistType::kHist)
					hist_ = new TH3D(name.c_str(), title.c_str(), axis.x().nbin(), &axis.x(0), axis.y().nbin(), &axis.y(0), axis.z().nbin(), &axis.z(0));
				else
					hist_ = new TProfile3D(name.c_str(), title.c_str(), axis.x().nbin(), &axis.x(0), axis.y().nbin(), &axis.y(0), axis.z().nbin(), &axis.z(0));
				
				hist_->GetXaxis()->SetNameTitle(axis.x().name().c_str(), axis.x().title().c_str());
				hist_->GetYaxis()->SetNameTitle(axis.y().name().c_str(), axis.y().title().c_str());
				hist_->GetZaxis()->SetNameTitle(axis.z().name().c_str(), axis.z().title().c_str());
				break;
			}
		default : break;
	}
	
	if (hist_ == nullptr) return;
	
	info_.first  = type;
	info_.second = axis.dim();

	axis_ = axis;
	name_ = name;
	title_ = title;
	
	style();
	push();
}


void Hist::setBin(const HistBin& bin, Double_t content, Double_t error) {
	if (!exist()) return;
	if (bin.dim() != info_.second) {
		MGSys::ShowError(STR_FMT("<< Hist::setBin >> Histogram dimension(%d) vs. HistBin dimension(%d).", static_cast<Int_t>(info_.second), static_cast<Int_t>(bin.dim())));
		return;
	}

	Bool_t isSetError = (MGNumc::Compare(error, 0.0) >= 0);
	switch (bin.dim()) {
		case HistDim::k1D : 
			{
				hist_->SetBinContent(bin.x(), content);
				if (isSetError)
					hist_->SetBinError  (bin.x(), error);
				break;
			}
		case HistDim::k2D : 
			{
				hist_->SetBinContent(bin.x(), bin.y(), content);
				if (isSetError)
					hist_->SetBinError  (bin.x(), bin.y(), error);
				break;
			}
		case HistDim::k3D : 
			{
				hist_->SetBinContent(bin.x(), bin.y(), bin.z(), content);
				if (isSetError)
					hist_->SetBinError  (bin.x(), bin.y(), bin.z(), error);
				break;
			}
		default : break;
	}
}


std::pair<Double_t, Double_t> Hist::getBin(const HistBin& bin) const {
	if (!exist()) return std::make_pair(0., 0.);
	if (bin.dim() != info_.second) {
		MGSys::ShowError(STR_FMT("<< Hist::getBin >> Histogram dimension(%d) vs. HistBin dimension(%d).", static_cast<Int_t>(info_.second), static_cast<Int_t>(bin.dim())));
		return std::make_pair(0., 0.);
	}
	
	Double_t val = 0., err = 0.;
	switch (bin.dim()) {
		case HistDim::k1D : 
			{
				val = hist_->GetBinContent(bin.x());
				err = hist_->GetBinError  (bin.x());
				break;
			}
		case HistDim::k2D : 
			{
				val = hist_->GetBinContent(bin.x(), bin.y());
				err = hist_->GetBinError  (bin.x(), bin.y());
				break;
			}
		case HistDim::k3D : 
			{
				val = hist_->GetBinContent(bin.x(), bin.y(), bin.z());
				err = hist_->GetBinError  (bin.x(), bin.y(), bin.z());
				break;
			}
		default : break;
	}
	return std::make_pair(val, err);
}


void Hist::scale(Double_t scl, Option_t * option) {
	if (!exist()) return;
	hist_->Scale(scl, option);
}


void Hist::normalized(const HistNorm& norm) {
	if (!exist()) return;
	if (!hist_->GetDefaultSumw2()) hist_->Sumw2();
	switch (norm) {
		case HistNorm::kEntries :
			{
				Double_t entries = (Double_t)hist_->GetEntries();
				if (MGNumc::Compare(entries) > 0) hist_->Scale(1. / entries);
			}
			break;
		case HistNorm::kIntegral :
			{
				Double_t integral = (Double_t)hist_->Integral();
				if (MGNumc::Compare(integral) > 0) hist_->Scale(1. / integral);
			}
			break;
		case HistNorm::kArea :
			{
				Double_t area = (Double_t)hist_->Integral("width");
				if (MGNumc::Compare(area) > 0) hist_->Scale(1. / area);
			}
			break;
		default : break;
	}
}


void Hist::style(const TAttFill& fill, const TAttLine& line, const TAttMarker& marker) {
	if (!exist()) return;
	hist_->SetFillColor(fill.GetFillColor());
	hist_->SetFillStyle(fill.GetFillStyle());
	hist_->SetLineColor(line.GetLineColor());
	hist_->SetLineStyle(line.GetLineStyle());
	hist_->SetLineWidth(line.GetLineWidth());
	hist_->SetMarkerColor(marker.GetMarkerColor());
	hist_->SetMarkerStyle(marker.GetMarkerStyle());
	hist_->SetMarkerSize (marker.GetMarkerSize());
}
		

void Hist::fill(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e) {
	if (!exist()) return;
	switch (info_.second) {
		case HistDim::k1D :
			if      (info_.first == HistType::kHist   ) dynamic_cast<TH1D*>    (hist_)->Fill(a, b); 
			else if (info_.first == HistType::kProfile) dynamic_cast<TProfile*>(hist_)->Fill(a, b, c); 
			break;
		case HistDim::k2D :
			if      (info_.first == HistType::kHist   ) dynamic_cast<TH2D*>      (hist_)->Fill(a, b, c); 
			else if (info_.first == HistType::kProfile) dynamic_cast<TProfile2D*>(hist_)->Fill(a, b, c, d); 
			break;
		case HistDim::k3D :
			if      (info_.first == HistType::kHist   ) dynamic_cast<TH3D*>      (hist_)->Fill(a, b, c, d); 
			else if (info_.first == HistType::kProfile) dynamic_cast<TProfile3D*>(hist_)->Fill(a, b, c, d, e); 
			break;
		default : break;
	}
}


void Hist::push() {
	hist_map_[name_] = this;
	unique_id_ = ++counter_;
	hist_->SetUniqueID(unique_id_);
}


void Hist::clear() {
	if (hist_ != nullptr) { 
		hist_->Delete();
		hist_ = nullptr;
		hist_map_.erase(name_); 
	} 
	axis_ = HistAxis(); 
	info_.first  = HistType::kNone;
	info_.second = HistDim::kNone;
	name_ = ""; 
	title_ = ""; 
	unique_id_ = -1; 
}


Hist * Hist::Head(const std::string& name, Bool_t show) {
	if (name == "") return nullptr;
	std::map<std::string, Hist*>::iterator search = hist_map_.find(name);
	if (search == hist_map_.end()) {
		if (show) MGSys::ShowError(STR_FMT("<< Hist::Head >> %s is not found.", name.c_str()));
		return nullptr;
	}
	return search->second;
}


void Hist::Delete(const std::string& name) {
	Hist * hist = Hist::Head(name, false);
	Delete(hist);
}


void Hist::Delete(std::vector<Hist*>& histvec) {
	for (auto* hist : histvec)
		Hist::Delete(hist);
}


void Hist::Write() {
	for (auto const& it : hist_map_) (it.second)->write();
}


Bool_t Hist::Load(const std::string& filename, const std::string& filepath) {
	std::string fullPath = STR_FMT("%s/%s.root", filepath.c_str(), filename.c_str());
	if (!MGSys::TestFile(fullPath, 'f')) {
		MGSys::ShowError(STR_FMT("<< Hist::Load >>  %s is not open.", fullPath.c_str()));
		return false;
	}

	TFile * file = TFile::Open(fullPath.c_str(), "READ");
	if (file != nullptr && file->IsOpen()) {
		TList * list = file->GetListOfKeys();
		TIter next((TList*)list);
		while (TObject * key = next()) {
			TString keyName = key->GetName();
			TH1 * obj = dynamic_cast<TH1*>(file->Get(keyName.Data()));
			if (obj == nullptr) continue;
			Hist::New(obj);
		}
		file->Close();
	}

	return true;
}


Bool_t Hist::Save(const std::string& filename, const std::string& filepath) {
	std::string fullPath = STR_FMT("%s/%s.root", filepath.c_str(), filename.c_str());
	if (!MGSys::TestFile(filepath, 'd')) return false;
	
	TFile * file = new TFile(fullPath.c_str(), "RECREATE");
	file->cd();
	
	Hist::Write();

	file->Write();
	file->Close();
	file->Delete();
	file = nullptr;
	
	return true;
}


void Hist::Print() {
	COUT("*********************************  Histogram List  *********************************\n");
	COUT("%-8s   %12s   %-18s   %-50s\n", "UniqueID", "ClassName", "Name", "Title");
	for (auto const& it : hist_map_) { 
		const TH1 * hist = (*(it.second))();
		COUT("%08u   %12s   %-18s   %-50s\n", hist->GetUniqueID(), hist->ClassName(), hist->GetName(), hist->GetTitle());
	}
	COUT("************************************************************************************\n");
}


Hist * Hist::Calculate(const HistArith& arith, const std::string& name, const std::string& title, Hist * hSon, Hist * hMom) {
	if (hSon == nullptr || !MGNumc::Valid((*hSon)()->Integral())) {
		MGSys::ShowError("<< Hist::Calculate >> hSon is nullptr or INF/NAN.");
		return nullptr; 
	}
	if (hMom == nullptr || !MGNumc::Valid((*hMom)()->Integral())) {
		MGSys::ShowError("<< Hist::Calculate >> hMom is nullptr or INF/NAN.");
		return nullptr; 
	}

	Hist * hist = Hist::New(name, title, (*hSon)());
	if (hist == nullptr) {
		MGSys::ShowError("<< Hist::Calculate >> Create histogram failure.");
		return nullptr;
	}

	TH1 * hSonMom = (*hist)();
	switch (arith) {
		case HistArith::kAddition : hSonMom->Add((*hMom)(),  1.0); break;
		case HistArith::kSubtract : hSonMom->Add((*hMom)(), -1.0); break;
		case HistArith::kMultiply : hSonMom->Multiply((*hMom)());  break;
		case HistArith::kDivide   : hSonMom->Divide((*hMom)());    break;
		default : break;
	}

	return hist;
}



Hist * Hist::Project(const HistProj& proj, Hist * hMom, Int_t isb, Int_t ieb, Int_t jsb, Int_t jeb) {
	if (hMom == nullptr || !hMom->exist()) {
		MGSys::ShowError("<< Hist::Project >> hMom is not existence.");
		return nullptr;
	}
	else if (hMom->dim() == HistDim::k1D) {
		MGSys::ShowError("<< Hist::Project >> hMom is 1D.");
		return nullptr;
	}
	else if (hMom->type() == HistType::kProfile) {
		MGSys::ShowError("<< Hist::Project >> hMom is Profile.");
		return nullptr;
	}

	Int_t  ndim  = static_cast<Int_t>(hMom->dim());
	Int_t  vProj = static_cast<Int_t>(proj);
	Bool_t optProj = (vProj <= 3);
	if (optProj && vProj > static_cast<Int_t>(hMom->dim())) {
		MGSys::ShowError("<< Hist::Project >> HistProj is higher than hMom dimension.");
		return nullptr;
	}
	if (!optProj && hMom->dim() != HistDim::k3D) {
		MGSys::ShowError("<< Hist::Project >> hMom is not 3D.");
		return nullptr;
	}
	if (isb > ieb || jsb > jeb) {
		MGSys::ShowError("<< Hist::Project >> Input bins is failure.");
		return nullptr;
	}

	Hist * hist = nullptr;
	if (hMom->dim() == HistDim::k2D) {
		switch (proj) {
			case HistProj::kX : 
				{
					const Axis& axis = hMom->axis().y();
					Int_t sb = ((isb < 0 || isb > axis.nbin()+1) ?           1 : isb);
					Int_t eb = ((ieb < 0 || ieb > axis.nbin()+1) ? axis.nbin() : ieb);
					std::string name = STR_FMT("%s__PROJX%05dTO%05d", hMom->name().c_str(), sb, eb);
					TH1D * projHist = ((TH2D*)(*hMom)())->ProjectionX(name.c_str(), sb, eb);
					hist = Hist::New(name, "", projHist);
					projHist->Delete();
					break;
				}
			case HistProj::kY :
				{
					const Axis& axis = hMom->axis().x();
					Int_t sb = ((isb < 0 || isb > axis.nbin()+1) ?           1 : isb);
					Int_t eb = ((ieb < 0 || ieb > axis.nbin()+1) ? axis.nbin() : ieb);
					std::string name = STR_FMT("%s__PROJY%05dTO%05d", hMom->name().c_str(), sb, eb);
					TH1D * projHist = ((TH2D*)(*hMom)())->ProjectionY(name.c_str(), sb, eb);
					hist = Hist::New(name, "", projHist);
					projHist->Delete();
					break;
				}
			default : break;
		}
	}
	else if (hMom->dim() == HistDim::k3D) {
		switch (proj) {
			case HistProj::kX :
				{
					const Axis& iaxis = hMom->axis().y();
					Int_t _isb = ((isb < 0 || isb > iaxis.nbin()+1) ?            1 : isb);
					Int_t _ieb = ((ieb < 0 || ieb > iaxis.nbin()+1) ? iaxis.nbin() : ieb);
					const Axis& jaxis = hMom->axis().z();
					Int_t _jsb = ((jsb < 0 || jsb > jaxis.nbin()+1) ?            1 : jsb);
					Int_t _jeb = ((jeb < 0 || jeb > jaxis.nbin()+1) ? jaxis.nbin() : jeb);
					std::string name = STR_FMT("%s__PROJY%05dTO%05d__PROJZ%05dTO%05d", hMom->name().c_str(), _isb, _ieb, _jsb, _jeb);
					(*hMom)()->GetYaxis()->SetRange(_isb, _ieb);
					(*hMom)()->GetZaxis()->SetRange(_jsb, _jeb);
					TH1D * projHist = (TH1D*)(((TH3D*)((*hMom)()))->Project3D("X"));
					hist = Hist::New(name, "", projHist);
					projHist->Delete();
					(*hMom)()->GetYaxis()->UnZoom();
					(*hMom)()->GetZaxis()->UnZoom();
					break;
				}
			case HistProj::kY :
				{
					const Axis& iaxis = hMom->axis().z();
					Int_t _isb = ((isb < 0 || isb > iaxis.nbin()+1) ?            1 : isb);
					Int_t _ieb = ((ieb < 0 || ieb > iaxis.nbin()+1) ? iaxis.nbin() : ieb);
					const Axis& jaxis = hMom->axis().x();
					Int_t _jsb = ((jsb < 0 || jsb > jaxis.nbin()+1) ?            1 : jsb);
					Int_t _jeb = ((jeb < 0 || jeb > jaxis.nbin()+1) ? jaxis.nbin() : jeb);
					std::string name = STR_FMT("%s__PROJZ%05dTO%05d__PROJX%05dTO%05d", hMom->name().c_str(), _isb, _ieb, _jsb, _jeb);
					(*hMom)()->GetZaxis()->SetRange(_isb, _ieb);
					(*hMom)()->GetXaxis()->SetRange(_jsb, _jeb);
					TH1D * projHist = (TH1D*)(((TH3D*)((*hMom)()))->Project3D("Y"));
					hist = Hist::New(name, "", projHist);
					projHist->Delete();
					(*hMom)()->GetZaxis()->UnZoom();
					(*hMom)()->GetXaxis()->UnZoom();
					break;
				}
			case HistProj::kZ :
				{
					const Axis& iaxis = hMom->axis().x();
					Int_t _isb = ((isb < 0 || isb > iaxis.nbin()+1) ?            1 : isb);
					Int_t _ieb = ((ieb < 0 || ieb > iaxis.nbin()+1) ? iaxis.nbin() : ieb);
					const Axis& jaxis = hMom->axis().y();
					Int_t _jsb = ((jsb < 0 || jsb > jaxis.nbin()+1) ?            1 : jsb);
					Int_t _jeb = ((jeb < 0 || jeb > jaxis.nbin()+1) ? jaxis.nbin() : jeb);
					std::string name = STR_FMT("%s__PROJX%05dTO%05d__PROJY%05dTO%05d", hMom->name().c_str(), _isb, _ieb, _jsb, _jeb);
					(*hMom)()->GetXaxis()->SetRange(_isb, _ieb);
					(*hMom)()->GetYaxis()->SetRange(_jsb, _jeb);
					TH1D * projHist = (TH1D*)(((TH3D*)((*hMom)()))->Project3D("Z"));
					hist = Hist::New(name, "", projHist);
					projHist->Delete();
					(*hMom)()->GetXaxis()->UnZoom();
					(*hMom)()->GetYaxis()->UnZoom();
					break;
				}
			case HistProj::kYZ :
			case HistProj::kZY :
				{
					std::string opt = ((proj == HistProj::kYZ) ? "YZ" : "ZY");
					const Axis& axis = hMom->axis().x();
					Int_t sb = ((isb < 0 || isb > axis.nbin()+1) ?           1 : isb);
					Int_t eb = ((ieb < 0 || ieb > axis.nbin()+1) ? axis.nbin() : ieb);
					std::string name = STR_FMT("%s__PROJ%s%05dTO%05d", hMom->name().c_str(), opt.c_str(), sb, eb);
					(*hMom)()->GetXaxis()->SetRange(sb, eb);
					TH2D * projHist = (TH2D*)(((TH3D*)((*hMom)()))->Project3D(opt.c_str()));
					hist = Hist::New(name, "", projHist);
					projHist->Delete();
					(*hMom)()->GetXaxis()->UnZoom();
					break;
				}
			case HistProj::kXZ :
			case HistProj::kZX :
				{
					std::string opt = ((proj == HistProj::kXZ) ? "XZ" : "ZX");
					const Axis& axis = hMom->axis().y();
					Int_t sb = ((isb < 0 || isb > axis.nbin()+1) ?           1 : isb);
					Int_t eb = ((ieb < 0 || ieb > axis.nbin()+1) ? axis.nbin() : ieb);
					std::string name = STR_FMT("%s__PROJ%s%05dTO%05d", hMom->name().c_str(), opt.c_str(), sb, eb);
					(*hMom)()->GetYaxis()->SetRange(sb, eb);
					TH2D * projHist = (TH2D*)(((TH3D*)((*hMom)()))->Project3D(opt.c_str()));
					hist = Hist::New(name, "", projHist);
					projHist->Delete();
					(*hMom)()->GetYaxis()->UnZoom();
					break;
				}
			case HistProj::kXY :
			case HistProj::kYX :
				{
					std::string opt = ((proj == HistProj::kXY) ? "XY" : "YX");
					const Axis& axis = hMom->axis().z();
					Int_t sb = ((isb < 0 || isb > axis.nbin()+1) ?           1 : isb);
					Int_t eb = ((ieb < 0 || ieb > axis.nbin()+1) ? axis.nbin() : ieb);
					std::string name = STR_FMT("%s__PROJ%s%05dTO%05d", hMom->name().c_str(), opt.c_str(), sb, eb);
					(*hMom)()->GetZaxis()->SetRange(sb, eb);
					TH2D * projHist = (TH2D*)(((TH3D*)((*hMom)()))->Project3D(opt.c_str()));
					hist = Hist::New(name, "", projHist);
					projHist->Delete();
					(*hMom)()->GetZaxis()->UnZoom();
					break;
				}
			default : break;
		}
	}

	return hist;
}
	

std::vector<Hist*> Hist::ProjectAll(const HistProj& proj, Hist * hMom) {
	std::vector<Hist *> histvec;
	if (hMom == nullptr || !hMom->exist()) return histvec;
	if (hMom->dim() == HistDim::k1D) return histvec;

	if (hMom->dim() == HistDim::k2D) {
		switch (proj) {
			case HistProj::kX :
			case HistProj::kY :
				{
					Int_t nbin = (proj == HistProj::kX) ? hMom->axis().y().nbin() : hMom->axis().x().nbin();
					if (nbin < 0) break;
					for (Int_t it = 0; it <= nbin+1; ++it) {
						Hist * hist = Hist::Project(proj, hMom, it, it);
						if (hist == nullptr || !hist->exist()) continue;
						histvec.push_back(hist);
					}
					break;
				}
			default : return histvec;
		}
	}
	
	if (hMom->dim() == HistDim::k3D) {
		switch (proj) {
			case HistProj::kX :
			case HistProj::kY :
			case HistProj::kZ :
				{
					Int_t nbinA = -1, nbinB = -1;
					if (proj == HistProj::kX) { nbinA = hMom->axis().y().nbin(); nbinB = hMom->axis().z().nbin(); }
					if (proj == HistProj::kY) { nbinA = hMom->axis().z().nbin(); nbinB = hMom->axis().x().nbin(); }
					if (proj == HistProj::kZ) { nbinA = hMom->axis().x().nbin(); nbinB = hMom->axis().y().nbin(); }
					if (nbinA < 0 || nbinB < 0) break;
					for (Int_t itA = 0; itA <= nbinA+1; ++itA) {
						for (Int_t itB = 0; itB <= nbinB+1; ++itB) {
							Hist * hist = Hist::Project(proj, hMom, itA, itA, itB, itB);
							if (hist == nullptr || !hist->exist()) continue;
							histvec.push_back(hist);
						}
					}
					break;
				}
			case HistProj::kYZ :
			case HistProj::kZY :
			case HistProj::kXZ :
			case HistProj::kZX :
			case HistProj::kXY :
			case HistProj::kYX :
				{
					Int_t nbin = -1;
					if (proj == HistProj::kYZ || proj == HistProj::kZY) nbin = hMom->axis().x().nbin();
					if (proj == HistProj::kXZ || proj == HistProj::kZX) nbin = hMom->axis().y().nbin();
					if (proj == HistProj::kXY || proj == HistProj::kYX) nbin = hMom->axis().z().nbin();
					if (nbin < 0) break;
					for (Int_t it = 0; it <= nbin+1; ++it) {
						Hist * hist = Hist::Project(proj, hMom, it, it);
						if (hist == nullptr || !hist->exist()) continue;
						histvec.push_back(hist);
					}
					break;
				}
			default : return histvec;
		}
	}
	
	return histvec;
}


THStack * Hist::Collect(const std::string& name, const std::string& title, const std::vector<std::string>& list) {
	if (list.size() == 0) return nullptr;
	THStack * coll = new THStack();
	coll->SetNameTitle(name.c_str(), title.c_str());
	for (auto&& elem : list) {
		Hist * hist = Hist::Head(elem);
		if (hist == nullptr) continue;
		coll->Add((*hist)());
	}
	return coll;
}


THStack * Hist::Collect(const std::string& name, const std::string& title, const std::vector<Hist*>& list) {
	if (list.size() == 0) return nullptr;
	THStack * coll = new THStack();
	coll->SetNameTitle(name.c_str(), title.c_str());
	for (auto&& hist : list) {
		if (hist == nullptr) continue;
		coll->Add((*hist)());
	}
	return coll;
}


} // namesapce MGROOT

#endif // __ROOTLibs_MGHist_C__
