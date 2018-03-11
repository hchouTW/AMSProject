#ifndef __ROOTLibs_MGCanvas_C__
#define __ROOTLibs_MGCanvas_C__

#include "MGCanvas.h"

namespace MGROOT {

TVirtualPad * SetPadWindow(TCanvas * canvas, UInt_t idx, const PadWindow& window) {
	if (canvas == nullptr) return nullptr;
	TVirtualPad * pad = canvas->cd(idx);
	if (pad == nullptr) return nullptr;
	
	pad->SetPad(window.xl, window.yl, window.xu, window.yu);
	
	pad->Modified();
    pad->Update();
	
	return pad;
}


TVirtualPad * SetPadMargin(TCanvas * canvas, UInt_t idx, const PadMargin& margin) {
	if (canvas == nullptr) return nullptr;
	TVirtualPad * pad = canvas->cd(idx);
	if (pad == nullptr) return nullptr;
	
	pad->SetTopMargin   (margin.top);
	pad->SetBottomMargin(margin.bottom);
	pad->SetLeftMargin  (margin.left);
	pad->SetRightMargin (margin.right);
	
	pad->Modified();
    pad->Update();
	
	return pad;
}


TVirtualPad * SetPadBorder(TCanvas * canvas, UInt_t idx, const PadBorder& border) {
	if (canvas == nullptr) return nullptr;
	TVirtualPad * pad = canvas->cd(idx);
	if (pad == nullptr) return nullptr;
	
	pad->SetFillColor(border.color);
	pad->SetFillStyle(1001);
	
	pad->Modified();
    pad->Update();
	
	pad->SetBorderMode(border.mode);
	pad->SetBorderSize(border.size);
	
	return pad;
}


TVirtualPad * SetPadAxis(TCanvas * canvas, UInt_t idx, const PadAxis& axis) {
	if (canvas == nullptr) return nullptr;
	TVirtualPad * pad = canvas->cd(idx);
	if (pad == nullptr) return nullptr;

	pad->SetLogx(axis.logx);
	pad->SetLogy(axis.logy);
	pad->SetLogz(axis.logz);
    
	pad->Modified();
	pad->Update();

	return pad;
}


//---- Canvas  ----//
void Canvas::create(const Window& window, const PadMargin& margin, const PadBorder& border) {
	window_ = window;
	
	canvas_.cd(0);
	canvas_.SetCanvasSize(window.width, window.height);
	SetPadMargin(&canvas_, 0, margin);
	SetPadBorder(&canvas_, 0, border);
	
	canvas_.Modified(); 
	canvas_.Update();
}


void Canvas::create(UInt_t ndivx, UInt_t ndivy, const Window& window, const PadMargin& margin, const PadBorder& border) {
	window_ = window;

	UInt_t  npad = (ndivx<1?1:ndivx) * (ndivy<1?1:ndivy);
	Float_t wmargin = 0.01;
	Float_t hmargin = 0.01;

	canvas_.cd(0);
	canvas_.SetCanvasSize(window.width, window.height);
	canvas_.Divide((ndivx<1?1:ndivx), (ndivy<1?1:ndivy), wmargin, hmargin, border.color);
	
	for (UInt_t it = 1; it <= npad; ++it) {
		SetPadMargin(&canvas_, it, margin);
		SetPadBorder(&canvas_, it, border);
	}

	canvas_.Modified(); 
	canvas_.Update();
}


//---- PdfEditor ----//
PdfEditor::PdfEditor(const Window& window, const std::string& filename, const std::string& filepath) : exist_(false) {
	if (!MGSys::TestFile(filepath, 'd')) return;
	file_page_ = 0;
	file_path_ = filepath;
	file_name_ = (filename != "") ? filename : "PdfPainter";
	canvas_.create(window);
	canvas_.setNameTitle("START", "START");

	std::string fullPathPDF = STR("%s/%s.pdf", file_path_.c_str(), file_name_.c_str());	
	COUT("\n\n**************** PdfEditor::OPEN    %-45s  ****************\n", fullPathPDF.c_str());
	TextDraw(STR("<< START >>  %s  ", file_name_.c_str()));
	canvas_.save(STR("%s(", fullPathPDF.c_str()));
	canvas_.clear();
	exist_ = true;
}


void PdfEditor::create(const std::string& title, const PadMargin& margin, const PadBorder& border) {
	if (!exist_) return;
	++file_page_;
	canvas_.create(canvas_.window(), margin, border);
	canvas_.setNameTitle(STR("%s__PAGE%06d", file_name_.c_str(), file_page_), title);
	
	std::string statement = STR("%s  >~@~>  %s", file_name_.c_str(), title.c_str());
	COUT("PdfEditor::NEW_PAGE  ~~@ %05d @~~  << %-57s >>\n", file_page_, statement.c_str());
}


void PdfEditor::create(const std::string& title, UInt_t ndivx, UInt_t ndivy, const PadMargin& margin, const PadBorder& border) { 
	if (!exist_) return;
	++file_page_;
	canvas_.create(ndivx, ndivy, canvas_.window(), margin, border); 
	canvas_.setNameTitle(STR("%s__PAGE%06d", file_name_.c_str(), file_page_), title);
	
	std::string statement = STR("%s  >~@~>  %s", file_name_.c_str(), title.c_str());
	COUT("PdfEditor::NEW_PAGE  ~~@ %05d @~~  << %-57s >>\n", file_page_, statement.c_str());
}


void PdfEditor::close() {
	if (!exist_) return;
	++file_page_;
	canvas_.create(canvas_.window());
	canvas_.setNameTitle("END", "END");
	
	std::string fullPathPDF = STR("%s/%s.pdf", file_path_.c_str(), file_name_.c_str());
	COUT("**************** PdfEditor::CLOSE   %-45s  ****************\n\n", fullPathPDF.c_str());
	TextDraw(STR("<< END >>  %s  ", file_name_.c_str()));
	canvas_.save(STR("%s)", fullPathPDF.c_str()));
	canvas_.clear();
	exist_ = false;
	file_name_ = "PdfEditor";
	file_path_ = ".";
	file_page_ = 0;
}
		

void PdfEditor::save() {
	if (!exist_) return;
	std::string fullPathPDF = STR("%s/%s.pdf", file_path_.c_str(), file_name_.c_str());
	canvas_.save(fullPathPDF);
	canvas_.clear();
}

} // namespace MGROOT

#endif // __ROOTLibs_MGCanvas_C__
