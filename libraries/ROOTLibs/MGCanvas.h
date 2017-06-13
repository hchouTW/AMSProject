#ifndef __ROOTLibs_MGCanvas_H__
#define __ROOTLibs_MGCanvas_H__

#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>

/********************/
/****  MGROOT  ****/
/********************/
namespace MGROOT {

//---- Pad ----//
struct PadWindow {
	PadWindow(Double_t _xl = 0.0, Double_t _xu = 1.0, Double_t _yl = 0.0, Double_t _yu = 1.0) : xl(_xl), xu(_xu), yl(_yl), yu(_yu) {}
	Double_t xl, xu;
	Double_t yl, yu;
};
inline TVirtualPad * SetPadWindow(TCanvas * canvas, UInt_t idx = 0, const PadWindow& window = PadWindow());

struct PadMargin {
	PadMargin(Double_t _top = 0.12, Double_t _bottom = 0.12, Double_t _left = 0.15, Double_t _right = 0.15) : top(_top), bottom(_bottom), left(_left), right(_right) {}
	Double_t top, bottom;
	Double_t left, right;
};
inline TVirtualPad * SetPadMargin(TCanvas * canvas, UInt_t idx = 0, const PadMargin& margin = PadMargin());

struct PadBorder {
	PadBorder(Color_t _color = kWhite, Short_t _mode = 0, Short_t _size = 0) : color(_color), mode(_mode), size(_size) {}
	Color_t color;
	Short_t mode, size;
};
inline TVirtualPad * SetPadBorder(TCanvas * canvas, UInt_t idx = 0, const PadBorder& border = PadBorder());

struct PadAxis {
	PadAxis(Bool_t _logx = false, Bool_t _logy = false, Bool_t _logz = true) : logx(_logx), logy(_logy), logz(_logz) {}
	Bool_t logx, logy, logz;
};
inline TVirtualPad * SetPadAxis(TCanvas * canvas, UInt_t idx = 0, const PadAxis& axis = PadAxis());


//---- Canvas ----//
class Canvas {
	public :
		Canvas(const std::string& name = "canvas", const std::string& title = "") { setNameTitle(name, title); }
		~Canvas() {}

		void create(const Window& window = Window(WindowSize::kSlice), const PadMargin& margin = PadMargin(), const PadBorder& border = PadBorder());
		void create(UInt_t ndivx, UInt_t ndivy, const Window& window = Window(WindowSize::kSlice), const PadMargin& margin = PadMargin(), const PadBorder& border = PadBorder());

		TVirtualPad * cd(UInt_t idx = 0, const PadAxis& axis = PadAxis()) { return SetPadAxis(&canvas_, idx, axis); }
		void save(const std::string& fullpath = "", Option_t * option = "") { canvas_.Modified(); canvas_.Update(); canvas_.SaveAs(fullpath.c_str(), option); }
		void write(const std::string& name = "", Int_t option = 0, Int_t bufsize = 0) { canvas_.Modified(); canvas_.Update(); canvas_.Write(name.c_str(), option, bufsize); }
		void clear(Option_t * option = "") { canvas_.Clear(option); canvas_.Modified(); canvas_.Update(); }
		void update() { canvas_.Modified(); canvas_.Update(); }
		void setNameTitle(const std::string& name, const std::string& title = "") { canvas_.SetName(name.c_str()); canvas_.SetTitle(title.c_str()); } 
		
		inline TCanvas& operator()() { return canvas_; }
		inline Window& window() { return window_; }

	protected :
		TCanvas canvas_;
		Window  window_;
};


//---- PdfEditor ----//
class PdfEditor {
	public :
		PdfEditor(const Window& window = Window(WindowSize::kSlice), const std::string& filename = "", const std::string& filepath = ".");
		~PdfEditor() { close(); }

		inline Bool_t exist() { return exist_; }

		inline TVirtualPad * cd(UInt_t idx = 0, const PadAxis& axis = PadAxis()) { return ((exist_) ? canvas_.cd(idx, axis) : nullptr); }
		
		inline void create(const std::string& title, const PadMargin& margin = PadMargin(), const PadBorder& border = PadBorder());
		inline void create(const std::string& title, UInt_t ndivx, UInt_t ndivy, const PadMargin& margin = PadMargin(), const PadBorder& border = PadBorder());
		
		void close();

		void write() { if (exist_) canvas_.write(canvas_().GetName()); }
		void save();

		inline Canvas& operator()() { return canvas_; }

	protected :
		Bool_t      exist_;
		Canvas      canvas_;
		std::string file_name_;
		std::string file_path_;
		UInt_t      file_page_;
};

}

#endif // __ROOTLibs_MGCanvas_H__
