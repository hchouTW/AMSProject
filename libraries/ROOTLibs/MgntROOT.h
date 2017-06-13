#ifndef __MgntROOT_H__
#define __MgntROOT_H__
#include "/afs/cern.ch/user/h/hchou/private/AMSProject/libraries/CPPLibs/CPPLibs.h" 

#include "TMath.h"

/********************/
/****  MgntROOT  ****/
/********************/
namespace MgntROOT {
	namespace Math {}
	namespace MtxLB {}
	namespace Fit {}

	class Color;
	class Area;
	class Line;
	class Marker;
	class Text;
	class Style;
	class Canvas;
	class PdfEditor;
	class Random;
	class Title;
	class Axis;
	class Hist;
	class Graph;
	class Func;
	//class Stat;
	//class MVA;
}


/*****************/
/****  Color  ****/
/*****************/
#include "TColor.h"
// kBlock, kBlue, kMagenta, kRed, kGreen, kCyan, kAzure, kViolet, kPink, kOrange, kSpring, kTeal, kYellow
class MgntROOT::Color {
	public :
		static Color_t At(UInt_t idx = 1, Int_t nset = -1);

	public :
		Color(UInt_t idx = 1, Int_t nset = -1) { fColor = Color::At(idx, nset); }
		~Color() {}
		
		Color_t& operator()() { return fColor; }

	protected :
		Color_t fColor;
};


/****************/
/****  Area  ****/
/****************/
#include "TAttFill.h"
class MgntROOT::Area {
	public :
		Area(Color_t color = MgntROOT::Color::At(), Style_t style = 0, Float_t alpha = 1.0) { fArea.SetFillColorAlpha(color, alpha); fArea.SetFillStyle(style); }
		Area(MgntROOT::Color color, Style_t style = 0, Float_t alpha = 1.0) { fArea.SetFillColorAlpha(color(), alpha); fArea.SetFillStyle(style); }
		~Area() {}

		TAttFill& operator()() { return fArea; }

	protected :
		TAttFill fArea;
};


/****************/
/****  Line  ****/
/****************/
#include "TAttLine.h"
class MgntROOT::Line {
	public :
		Line(Color_t color = MgntROOT::Color::At(), UInt_t style = 0, Width_t width = 2) { fLine = TAttLine(color, (style%10), width); }
		Line(MgntROOT::Color color, UInt_t style = 0, Width_t width = 2) { fLine = TAttLine(color(), (style%10), width); }
		~Line() {}

		TAttLine& operator()() { return fLine; }

	protected :
		TAttLine fLine;
};


/******************/
/****  Marker  ****/
/******************/
#include "TAttMarker.h"
class MgntROOT::Marker {
	public :
		enum Opt {
			kNumType = 2,
			kNumShape = 7,
			kFull = 0, kOpen = 1,
			kCircle = 0, kSquare = 1, kTriangleUp = 2, kTriangleDown = 3, 
			kDiamond = 4, kCross = 5, kStar = 6,
		};
		
		static Style_t Style(UInt_t idx = 0);
		static Style_t Style(Marker::Opt shape, Marker::Opt type) { return Marker::Style(shape * 2 + type); }
	
	public :
		Marker(Color_t color = MgntROOT::Color::At(), UInt_t style = 0, Size_t size = 3.0) { fMarker = TAttMarker(color, Marker::Style(style), size); }
		Marker(MgntROOT::Color color, UInt_t style = 0, Size_t size = 3.0) { fMarker = TAttMarker(color(), Marker::Style(style), size); }
		~Marker() {}

		TAttMarker& operator()() { return fMarker; }

	protected :
		TAttMarker fMarker;
};	


/****************/
/****  Text  ****/
/****************/
#include "TText.h"
#include "TLatex.h"
class MgntROOT::Text {
	public :
		struct Txt_t {
			Txt_t(const std::string& _text = "", Color_t _color = MgntROOT::Color::At(), Double_t _size = 0.04, Int_t _font = 42) :
				text(_text), color(_color), size(_size), font(_font) {}
			Txt_t(const std::string& _text, MgntROOT::Color _color, Double_t _size = 0.04, Int_t _font = 42) :
				text(_text), color(_color()), size(_size), font(_font) {}
			std::string text;
			Color_t     color;
			Double_t    size;
			Int_t       font;
		};

		struct Att_t {
			Att_t(Double_t _x = 0.5, Double_t _y = 0.5, Int_t _align = 22, Double_t _angle = 0.0) :
				x(_x), y(_y), align(_align), angle(_angle) {}
			Double_t x, y;
			Int_t    align;
			Double_t angle;
		};

	public :
		static void Draw(const Text::Txt_t& txt, const Text::Att_t& att);

	public :
		Text() {}
		~Text() {}
};


/*****************/
/****  Style  ****/
/*****************/
#include "TStyle.h"
class MgntROOT::Style {
	public :
		static TStyle& Head() { return fStyle; }
		static void Update() { fStyle.cd(); }
		static void LoadDefaultEnvironment();

		// kGreyScale, kRainBow, kBird
		static void ColorPalette(Int_t pattern = kRainBow);

		static void WindowsScale(UInt_t w = 1024, UInt_t h = 768);

		static void AxisScale(Bool_t xaxis = false, Bool_t yaxis = false, Bool_t zaxis = true);

	public :
		Style() {}
		~Style() {}

	protected :
		static TStyle fStyle;
};

TStyle MgntROOT::Style::fStyle("MainStyle", "MainStyle");


/******************/
/****  Canvas  ****/
/******************/
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
class MgntROOT::Canvas {
	public :
		struct Size_t { 
			Size_t(UInt_t _w, UInt_t _h) : 
				w(_w), h(_h) {} 
			UInt_t w; UInt_t h; 
		};

		struct Pad_t { 
			Pad_t(Double_t _xlw = 0.0, Double_t _xup = 1.0, Double_t _ylw = 0.0, Double_t _yup = 1.0, Color_t _color = 10, Short_t _bordermode = 0, Short_t _bordersize = 0):
				xlw(_xlw), xup(_xup), ylw(_ylw), yup(_yup), color(_color), bordermode(_bordermode), bordersize(_bordersize) {}
			Double_t xlw, xup, ylw, yup;
			Color_t color;
			Short_t bordermode, bordersize;
		};
	
		struct AxisScl_t {
			AxisScl_t(Bool_t _x = false, Bool_t _y = false, Bool_t _z = true) :
				x(_x), y(_y), z(_z) {}
			Bool_t x, y, z;
		};

		enum Size { kSliceH, kSliceM, kSliceL, kA4Vertical, kA4Horizon, kMacH, kMacM, kMacL, kMovie };
		enum Temp { kNone, kStat };

	public :
		Canvas(const std::string& name = "Canvas", const std::string& title = "", Canvas::Size size = Canvas::kSliceL, Canvas::Temp temp = Canvas::kNone) 
			{ setNameTitle(name, title); setSize(size); setTemp(temp); }
		~Canvas() {}
	
		inline void setNameTitle(const std::string& name = "Canvas", const std::string& title = "")
			{ fCanvas.SetName(name.c_str()); fCanvas.SetTitle(title.c_str()); }
		
		inline void setName(const std::string& name)
			{ fCanvas.SetName(name.c_str()); }

		inline void setTitle(const std::string& title)
			{ fCanvas.SetTitle(title.c_str()); }
		
		inline void create(Canvas::Size size = Canvas::kSliceL, Canvas::Temp temp = Canvas::kNone) 
			{ clear(); setSize(size); setTemp(temp); update(); }

		inline void create(Canvas::Size size, Int_t nx, Int_t ny, Float_t wmargin = 0.01, Float_t hmargin = 0.01, Int_t color = 0) 
			{ clear(); setSize(size); setTemp(nx, ny, wmargin, hmargin, color); update(); }

		TVirtualPad * cd(UInt_t idx, const Canvas::AxisScl_t& scl = Canvas::AxisScl_t());
		
		inline void update() 
			{ fCanvas.Modified(); fCanvas.Update(); }
		
		inline void write(const std::string& name = "", Int_t option = 0, Int_t bufsize = 0) 
			{ update(); fCanvas.Write(name.c_str(), option, bufsize); }

		inline void save(const std::string& fullpath = "", Option_t * option = "") 
			{ update(); fCanvas.SaveAs(fullpath.c_str(), option); fCanvas.Clear(); }

		inline void clear(Option_t * option = "") 
			{ fCanvas.Clear(option); update(); }

		inline std::string name() { return std::string(fCanvas.GetName()); }
		inline std::string title() { return std::string(fCanvas.GetTitle()); }
		inline TCanvas& operator()() { return fCanvas; }

		inline Canvas::Size getSize() { return fSize; }

		inline Canvas::Temp getTemp() { return fTemp; }

	protected :
		Bool_t checkPad(UInt_t idx) { return (fCanvas.GetPad(idx) != nullptr); }
		TVirtualPad * setPad(UInt_t idx, const Pad_t& pad);
		void setSize(Canvas::Size size = Canvas::kSliceL);
		void setTemp(Canvas::Temp temp = Canvas::kNone);
		void setTemp(Int_t nx, Int_t ny, Float_t wmargin = 0.01, Float_t hmargin = 0.01, Int_t color = 0);
	
	protected :
		TCanvas      fCanvas;
		Canvas::Size fSize;
		Canvas::Temp fTemp;
};


/*********************/
/****  PdfEditor  ****/
/*********************/
class MgntROOT::PdfEditor {
	public :
		enum Size { kSlice, kA4Vertical, kA4Horizon, kMac, kMovie };
		enum Temp { kNone, kStat };
	
	public :
		PdfEditor(PdfEditor::Size size = PdfEditor::kSlice, const std::string& filename = "", const std::string& filepath = ".", Option_t * option = "") : fExist(false) { open(size, filename, filepath, option); }
		~PdfEditor() { close(); }

		inline Bool_t exist() { return fExist; }

		void open(PdfEditor::Size size = PdfEditor::kSlice, const std::string& filename = "", const std::string& filepath = ".", Option_t * option = "");

		void create(const std::string& title = "", PdfEditor::Temp temp = PdfEditor::kNone);
		
		void create(const std::string& title, Int_t nx, Int_t ny);
		
		inline TVirtualPad * cd(UInt_t idx, const Canvas::AxisScl_t& scl = Canvas::AxisScl_t()) 
			{ return (exist() ? fCanvas.cd(idx, scl) : nullptr); }

		inline void write()
			{ if (exist()) fCanvas.write(fCanvas.name()); }
		
		inline void save()
			{ if (exist()) fCanvas.save(Form("%s/%s.pdf", fFilePath.c_str(), fFileName.c_str())); }
		
		inline void close();

		inline MgntROOT::Canvas& operator()() { return fCanvas; }

	protected :
		Bool_t           fExist;
		MgntROOT::Canvas fCanvas;
		std::string      fFileName;
		std::string      fFilePath;
		UInt_t           fFilePage;
};


/******************/
/****  Random  ****/
/******************/
#include "TRandom3.h"
class MgntROOT::Random {
	public :
		static inline Double_t Uniform  (Double_t x1 = 0.0, Double_t x2 = 1.0)      { return fRndmEngTR3.Uniform(x1, x2); }
		static inline Double_t Gaus     (Double_t mean = 0.0, Double_t sigma = 1.0) { return fRndmEngTR3.Gaus(mean, sigma); }
		static inline Double_t Landau   (Double_t mpv = 0.0, Double_t sigma = 1.0)  { return fRndmEngTR3.Landau(mpv, sigma); }
		static inline Double_t Exp      (Double_t tau = 1.0)                        { return fRndmEngTR3.Exp(tau); }
		static inline Int_t    Binomial (Int_t ntot = 1, Double_t prob = 0.5)       { return fRndmEngTR3.Binomial(ntot, prob); }
		static inline Int_t    Poisson  (Double_t mean = 1.0)                       { return fRndmEngTR3.Poisson(mean); }
	
	public :
		Random() {}
		~Random() {}

	protected :
		static TRandom3 fRndmEngTR3;
};

TRandom3 MgntROOT::Random::fRndmEngTR3(0);


/*****************/
/****  Title  ****/
/*****************/
class MgntROOT::Title {
	public :
		Title(std::string _title = "", std::string _x = "", std::string _y = "", std::string _z = "") : fTitle(_title), fX(_x), fY(_y), fZ(_z) {}
		~Title() {}
		
		std::string & title() { return fTitle; }
		std::string & x() { return fX; }
		std::string & y() { return fY; }
		std::string & z() { return fZ; }

	protected :	
		std::string fTitle;
		std::string fX;
		std::string fY;
		std::string fZ;
};


/****************/
/****  Axis  ****/
/****************/
/*  MgntROOT::Axis : Example
 *  nbin = 6, bins = {0, 1, 2, 3, 4, 5, 6}
 *
 *  DEMO :
 *         -----|---|---|---|---|---|---|-----
 *  bound       0   1   2   3   4   5   6
 *               \ / \ / \ / \ / \ / \ /
 *  bin      0    1   2   3   4   5   6    7
 */
#include "TAxis.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
class MgntROOT::Axis {
	public :
		using LIST = std::vector<Double_t>;
		enum Scl { kLinear = 0, kLog = 1 };
		enum Dim { kX = 1, kY = 2, kZ = 3 };
	
	public :
		static Axis Invert(const std::string& title, Axis& axis);
		static TH1D * CreateHist(const std::string& name, const std::string& title, Axis& axisX);
		static TH2D * CreateHist(const std::string& name, const std::string& title, Axis& axisX, Axis& axisY);
		static TH3D * CreateHist(const std::string& name, const std::string& title, Axis& axisX, Axis& axisY, Axis& axisZ); 
	
	public :
		inline const std::string& title() { return fTitle; }
		inline Bool_t             exist() { return (fList.size() > 1); }
		inline Int_t              nbin()  { return (exist() ? (fList.size()-1) : 0); }
		inline const Double_t *   bins()  { return (exist() ? (&fList.at(0)) : 0); }
		inline Double_t           min()   { return (exist() ? fList.at(0) : 0.0); }
		inline Double_t           max()   { return (exist() ? fList.at(nbin()) : 0.0); }

		inline Double_t bins(UInt_t idx);
		inline Double_t width(Int_t ibin);
		inline Double_t center(UInt_t ibin, Axis::Scl scl = Axis::kLinear);
		
		Int_t find(Double_t val);	
		
		void print(Axis::Scl scl = Axis::kLinear, std::ostream& out = std::cout);

	public :
		Axis() {}
		~Axis() {}

		Axis(const Axis& axis);
		Axis(Axis& axis, UInt_t mergeFT = 1);
		Axis(std::initializer_list<Double_t> list);
		Axis(const std::string& title, const std::vector<Double_t>& list, UInt_t mergeFT = 1);
		Axis(const std::string& title, const UInt_t nbin, const Double_t * bins, UInt_t mergeFT = 1);
		Axis(const std::string& title, const UInt_t nbin, Double_t lw, Double_t up, Axis::Scl scl = Axis::kLinear);
		Axis(const std::string& title, TAxis * axis, UInt_t mergeFT = 1);
		Axis(TAxis * axis, UInt_t mergeFT = 1);
		Axis(const std::string& title, TH1 * hist, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
		Axis(TH1 * hist, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
		Axis(const std::string& title, TObject * obj, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
		Axis(TObject * obj, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
		Axis(const std::string& title, MgntROOT::Hist * hist, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
		Axis(MgntROOT::Hist * hist, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);

	protected :
		Bool_t init_TAxis(TAxis * axis, UInt_t mergeFT = 1);
		Bool_t init_TAxis(const std::string& title, TAxis * axis, UInt_t mergeFT = 1);
		Bool_t init_TH1(TH1 * hist, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
		Bool_t init_TH1(const std::string& title, TH1 * hist, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
		Bool_t init_TObject(TObject * obj, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
		Bool_t init_TObject(const std::string& title, TObject * obj, Axis::Dim dim = Axis::kX, UInt_t mergeFT = 1);
	
		Bool_t check(const std::vector<Double_t>& list);
		Bool_t check(UInt_t ibin);
		Bool_t merge(const std::vector<Double_t>& list, UInt_t mergeFT = 1); 
	
	protected :
		std::string           fTitle;
		std::vector<Double_t> fList;
};



/****************/
/****  Hist  ****/
/****************/
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TKey.h"
#include "TObjArray.h"
#include "THStack.h"
class MgntROOT::Hist {
	public :
		enum NormOpt { kEntries, kIntegral, kArea };
		enum ArithOpt { kAddition, kSubtract, kMultiply, kDivide };
		enum ProjOpt { kProjX, kProjY, kProjYZ };

		using LIST = std::vector<std::string>;
		using VLIST = std::vector<Hist *>;

	public :
		static Hist * Head(const std::string& name);
		static Bool_t Exist(const std::string& name);
		static void   Delete(const std::string& name);
		static void   Delete(Hist * hist);
		static void   Delete(std::vector<Hist *>& histvec);

		static void Load(const std::string& path);
		static void Write();
		static void Save(const std::string& filename, const std::string& filepath = ".", const std::string& dirname = "");
		static void Print(std::ostream& out = std::cout);
		
		static Hist * Calculate(Hist::ArithOpt opt, const std::string& name, const std::string& title, Hist * hNum, Hist * hDen);

		static Hist * Addition(const std::string& name, const std::string& title, Hist * hNum, Hist * hDen) 
			{ return Calculate(Hist::kAddition, name, title, hNum, hDen); }
		static Hist * Subtract(const std::string& name, const std::string& title, Hist * hNum, Hist * hDen) 
			{ return Calculate(Hist::kSubtract, name, title, hNum, hDen); }
		static Hist * Multiply(const std::string& name, const std::string& title, Hist * hNum, Hist * hDen)
			{ return Calculate(Hist::kMultiply, name, title, hNum, hDen); }
		static Hist * Divide  (const std::string& name, const std::string& title, Hist * hNum, Hist * hDen)
			{ return Calculate(Hist::kDivide, name, title, hNum, hDen); }

		static Hist *              Project(Hist::ProjOpt opt, Hist * hmon, Int_t ibin, Int_t jbin);
		static Hist *              Project(Hist::ProjOpt opt, const std::string& nmon, Int_t ibin, Int_t jbin);
		static std::vector<Hist *> Project(Hist::ProjOpt opt, Hist * hmon, UInt_t mergeFT = 1);
		static std::vector<Hist *> Project(Hist::ProjOpt opt, const std::string& nmon, UInt_t mergeFT = 1);
		static Hist *              ProjectAll(Hist::ProjOpt opt, Hist * hmon, Bool_t lwflow = false, Bool_t upflow = false);
		static Hist *              ProjectAll(Hist::ProjOpt opt, const std::string& nmon, Bool_t lwflow = false, Bool_t upflow = false);
		
		static Hist *              Project3D(Hist::ProjOpt opt, Hist * hmon, Int_t ibin, Int_t jbin);
		static Hist *              Project3D(Hist::ProjOpt opt, const std::string& nmon, Int_t ibin, Int_t jbin);
		static std::vector<Hist *> Project3D(Hist::ProjOpt opt, Hist * hmon, UInt_t mergeFT = 1);
		static std::vector<Hist *> Project3D(Hist::ProjOpt opt, const std::string& nmon, UInt_t mergeFT = 1);
		static Hist *              Project3DAll(Hist::ProjOpt opt, Hist * hmon, Bool_t lwflow = false, Bool_t upflow = false);
		static Hist *              Project3DAll(Hist::ProjOpt opt, const std::string& nmon, Bool_t lwflow = false, Bool_t upflow = false);

		static THStack * Collect(const std::string& name, const std::string& title, const Hist::LIST& list);
		static THStack * Collect(const std::string& name, const std::string& title, const Hist::VLIST& vlist);
		
	public :
		static Hist * New(TH1 * hist, Bool_t reset = false);
		static Hist * New(const std::string& name, TH1 * hist, Bool_t reset = false);
		static Hist * New(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, Bool_t isProfile = false);
		static Hist * New(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, Bool_t isProfile2D = false);
		static Hist * New(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, MgntROOT::Axis& axisZ, Bool_t isProfile3D = false);

	protected :
		static Long64_t                               fCounter;
		static std::map<std::string, MgntROOT::Hist*> fHistMap;

	protected :
		Hist() { TH1::SetDefaultSumw2(true); TH1::AddDirectory(false); fUniqueID = -1; fName = ""; fTitle = ""; fType.first = 0; fType.second = false; fHist = nullptr; }

	public :
		Hist(TH1 * hist, Bool_t reset = false);
		Hist(const std::string& name, TH1 * hist, Bool_t reset = false);
		Hist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, Bool_t isProfile = false);
		Hist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, Bool_t isProfile2D = false);
		Hist(const std::string& name, const std::string& title, MgntROOT::Axis& axisX, MgntROOT::Axis& axisY, MgntROOT::Axis& axisZ, Bool_t isProfile3D = false);

		~Hist() { clear(); } 
	
		inline Bool_t exist() { return (fHist != nullptr); }

		void fill(Double_t a, Double_t b = 1.0, Double_t c = 1.0, Double_t d = 1.0, Double_t e = 1.0);
		void setContent(Double_t content, Int_t ibin, Int_t jbin = -1, Int_t kbin = -1);
		void setError(Double_t error, Int_t ibin, Int_t jbin = -1, Int_t kbin = -1);
		void setContentWithError(Double_t content, Double_t error, Int_t ibin, Int_t jbin = -1, Int_t kbin = -1);

		Double_t getContent(Int_t ibin, Int_t jbin = -1, Int_t kbin = -1);
		Double_t getError(Int_t ibin, Int_t jbin = -1, Int_t kbin = -1);
		std::pair<Double_t, Double_t> getContentWithError(Int_t ibin, Int_t jbin = -1, Int_t kbin = -1);

		void scale(Double_t scl = 1., Option_t * option = "");
		void normalized(Hist::NormOpt opt = Hist::kEntries);

		void setArea(MgntROOT::Area area);
		void setLine(MgntROOT::Line line);
		void setMarker(MgntROOT::Marker marker);
		void setStyle(MgntROOT::Area area = MgntROOT::Area(), MgntROOT::Line line = MgntROOT::Line(), MgntROOT::Marker marker = MgntROOT::Marker());

		void draw(Option_t * option = "") 
			{ if (exist()) fHist->Draw(option); }

		void draw(Option_t * option,  MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker);
		
		inline void write(const std::string& name = "", Int_t option = 0, Int_t bufsize = 0) 
			{ if (exist()) fHist->Write(name.c_str(), option, bufsize); }
	
		inline void clear()
			{ if (fHist != nullptr) { fHistMap.erase(fName); fHist->Delete(); fHist = nullptr; };
			  fUniqueID = -1; fName = ""; fTitle = ""; fType.first = 0; fType.second = false; }

		inline std::string name()       { return fName; }
		inline std::string title()      { return fTitle; }
		inline Short_t     ndim()       { return fType.first; }
		inline Bool_t      isProfile()  { return fType.second; }
		inline Long64_t    uniqueID()   { return fUniqueID; }
		inline TH1 *       operator()() { return fHist; }
		
		inline MgntROOT::Axis * axisX() { return &fAxisX; }
		inline MgntROOT::Axis * axisY() { return &fAxisY; }
		inline MgntROOT::Axis * axisZ() { return &fAxisZ; }

		Hist(const Hist&) = delete;
		Hist& operator= (const Hist&) = delete;

	protected :
		inline Bool_t isHistExist(const std::string& name);
		inline void setType(TH1 * hist);
		inline void pushBack() 
			{ fHistMap[fName] = this; fUniqueID = fCounter++;  (*this)()->SetUniqueID(fUniqueID); }

	protected :
		Long64_t                   fUniqueID;
		std::string                fName;
		std::string                fTitle;
		std::pair<Short_t, Bool_t> fType; // (dim, profile)
		TH1 *                      fHist;
		MgntROOT::Axis             fAxisX;
		MgntROOT::Axis             fAxisY;
		MgntROOT::Axis             fAxisZ;
	
	protected :
		static void   DebugON() { fDebug = true; }
		static void   DebugOFF() { fDebug = false; }
		static Bool_t fDebug;
};

Long64_t                               MgntROOT::Hist::fCounter = 1;
std::map<std::string, MgntROOT::Hist*> MgntROOT::Hist::fHistMap;
Bool_t                                 MgntROOT::Hist::fDebug = true;

/*****************/
/****  Graph  ****/
/*****************/
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
class MgntROOT::Graph {
	public :
		struct Point_t {
			Point_t(Double_t _x, Double_t _y) : x(_x), y(_y) {}
			Double_t x, y;
		};
		
		struct Error_t {
			Error_t(Double_t _ex, Double_t _ey) : exl(_ex), exu(_ex), eyl(_ey), eyu(_ey) {}
			Error_t(Double_t _exl, Double_t _exu, Double_t _eyl, Double_t _eyu) : exl(_exl), exu(_exu), eyl(_eyl), eyu(_eyu) {}
			Double_t exl, exu;
			Double_t eyl, eyu;
		};

		using LIST = std::vector<std::string>;
		using VLIST = std::vector<Graph *>;
	
	public :
		static Graph * Head(const std::string& name);
		static Bool_t  Exist(const std::string& name);
		static void    Delete(const std::string& name);
		static void    Delete(Graph * graph);
		static void    Delete(std::vector<Graph *>& graphvec);

		static void Load(const std::string& path);
		static void Write();
		static void Save(const std::string& filename, const std::string& filepath = ".", const std::string& dirname = "");
		static void Print(std::ostream& out = std::cout);
		
		static TMultiGraph * Collect(const std::string& name, const std::string& title, const Graph::LIST& list);
		static TMultiGraph * Collect(const std::string& name, const std::string& title, const Graph::VLIST& vlist);
	
	public :
		static Graph * New(Hist * hist, Bool_t islogx = false);
		static Graph * New(TGraph * graph, Bool_t reset = false);
		static Graph * New(const std::string& name, const std::string& title, const std::string& axtlex = "", const std::string& axtley = "");
	
	protected :
		static Long64_t                                fCounter;
		static std::map<std::string, MgntROOT::Graph*> fGraphMap;

	protected :
		Graph() { fUniqueID = -1; fName = ""; fTitle = ""; fGraph = nullptr; }

	public :
		Graph(Hist * hist, Bool_t islogx = false);
		Graph(TGraph * graph, Bool_t reset = false);
		Graph(const std::string& name, const std::string& title, const std::string& axtlex = "", const std::string& axtley = "");

		~Graph() { clear(); }
		
		inline Bool_t exist() { return (fGraph != nullptr); }
		
		inline void setPoint(UInt_t it, const Graph::Point_t& pnt);
		inline void setError(UInt_t it, const Graph::Error_t& err);
		inline void setPointWithError(UInt_t it, const Graph::Point_t& pnt, const Graph::Error_t& err);
		
		inline void pushPoint(const Graph::Point_t& pnt);
		inline void pushError(const Graph::Error_t& err);
		inline void pushPointWithError(const Graph::Point_t& pnt, const Graph::Error_t& err);
		
		void setArea(MgntROOT::Area area);
		void setLine(MgntROOT::Line line);
		void setMarker(MgntROOT::Marker marker);
		void setStyle(MgntROOT::Area area = MgntROOT::Area(), MgntROOT::Line line = MgntROOT::Line(), MgntROOT::Marker marker = MgntROOT::Marker());

		void draw(Option_t * option = "") 
			{ if (exist()) { fGraph->GetHistogram()->SetLineColor(0); fGraph->Draw(option); } }

		void draw(Option_t * option,  MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker);
		
		inline void write(const std::string& name = "", Int_t option = 0, Int_t bufsize = 0) 
			{ if (exist()) fGraph->Write(name.c_str(), option, bufsize); }
	
		inline void clear()
			{ if (fGraph != nullptr) { fGraphMap.erase(fName); fGraph->Delete(); fGraph = nullptr; };
			  fUniqueID = -1; fName = ""; fTitle = ""; }
		
		inline std::string         name()       { return fName; }
		inline std::string         title()      { return fTitle; }
		inline Long64_t            uniqueID()   { return fUniqueID; }
		inline TGraphAsymmErrors * operator()() { return fGraph; }

		Graph(const Graph&) = delete;
		Graph& operator= (const Graph&) = delete;

	protected :
		inline Bool_t isGraphExist(const std::string& name);
		inline void pushBack() 
			{ fGraphMap[fName] = this; fUniqueID = fCounter++;  (*this)()->SetUniqueID(fUniqueID); }

	protected :
		Long64_t            fUniqueID;
		std::string         fName;
		std::string         fTitle;
		TGraphAsymmErrors * fGraph;
};

Long64_t                                MgntROOT::Graph::fCounter = 1;
std::map<std::string, MgntROOT::Graph*> MgntROOT::Graph::fGraphMap;


/****************/
/****  Func  ****/
/****************/
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
class MgntROOT::Func {
	public :
		struct Bound {
			Bound(Double_t _lw = 0., Double_t _up = 0.) : lw(_lw), up(_up) {}
			Double_t lw, up;
		};

		struct Param {
			Param() : name(""), value(0.), error(0.), lmtlw(0.), lmtup(0.), fixed(false) {}
			Param(Double_t _value, Bool_t _fixed = false) : name(""), value(_value), error(0.), lmtlw(0.), lmtup(0.), fixed(_fixed) {}
			Param(Double_t _value, Double_t _lmtlw, Double_t _lmtup) : name(""), value(_value), error(0.), lmtlw(_lmtlw), lmtup(_lmtup), fixed(false) {}
			Param(const std::string& _name, Double_t _value = 0., Bool_t _fixed = false) : name(_name), value(_value), error(0.), lmtlw(0.), lmtup(0.), fixed(_fixed) {}
			Param(const std::string& _name, Double_t _value, Double_t _lmtlw, Double_t _lmtup) : name(_name), value(_value), error(0.), lmtlw(_lmtlw), lmtup(_lmtup), fixed(false) {}
			std::string name;
			Double_t    value;
			Double_t    error;
			Double_t    lmtlw, lmtup;
			Bool_t      fixed;
		};

		using LIST  = std::vector<Param>;
		using VLIST = std::vector<Double_t>;

		struct FitResult {
			FitResult() : valid(false), ndf(0), chisq(0.), prob(0.) {}
			FitResult(TFitResultPtr& ptr);

			Bool_t                             valid;
			Int_t                              ndf;
			Double_t                           chisq;
			Double_t                           prob;
			std::vector<Func::Param>           params;
			std::vector<std::vector<Double_t>> covmtx;
			std::vector<Double_t>              CI95;
			std::vector<Double_t>              CI68;

			inline std::string& name(Int_t i)  { return params.at(i).name; }
			inline Double_t&    value(Int_t i) { return params.at(i).value; }
			inline Double_t&    error(Int_t i) { return params.at(i).error; }
			inline Double_t&    lmtlw(Int_t i) { return params.at(i).lmtlw; }
			inline Double_t&    lmtup(Int_t i) { return params.at(i).lmtup; }
			inline Bool_t&      fixed(Int_t i) { return params.at(i).fixed; }
		};

	public :
		static Func * Head(const std::string& name);
		static Bool_t Exist(const std::string& name);
		static void   Delete(const std::string& name);
		static void   Delete(Func * func);
		static void   Delete(std::vector<Func *>& funcvec);
		
		static void Load(const std::string& path);
		static void Write();
		static void Save(const std::string& filename, const std::string& filepath = ".", const std::string& dirname = "");
		static void Print(std::ostream& out = std::cout);
		
	public :
		static Func * New(TF1 * func, const std::string& title = "");
		static Func * New(const std::string& name, const std::string& title, const std::string& formula, const Bound& bound = Bound(), const LIST& list = LIST());
		static Func * New(const std::string& name, const std::string& title, const std::string& formula, const Bound& bound, const VLIST& vlist);

	protected :
		static Long64_t                               fCounter;
		static std::map<std::string, MgntROOT::Func*> fFuncMap;
	
	protected :
		Func() { fUniqueID = -1; fName = ""; fTitle = ""; fFunc = nullptr; }

	public :
		Func(TF1 * func, const std::string& title = "");
		Func(const std::string& name, const std::string& title, const std::string& formula, const Bound& bound = Bound(), const LIST& list = LIST());
		Func(const std::string& name, const std::string& title, const std::string& formula, const Bound& bound, const VLIST& vlist);
		
		~Func() { clear(); }
		
		inline Bool_t exist() { return (fFunc != nullptr); }

		void setParam(UInt_t ipar, const Param& par);
		void setParam(UInt_t ipar, Double_t vpar, Bool_t fixed = false);
		void setParam(const LIST& list);
		void setParam(const VLIST& vlist);

		Func::Param getParam(UInt_t ipar); 

		Func::FitResult fitTo(Hist * hist, Option_t * option = "", Option_t * goption = "", Double_t xmin = 0., Double_t xmax = 0.);
		Func::FitResult fitTo(Graph * graph, Option_t * option = "", Option_t * goption = "", Double_t xmin = 0., Double_t xmax = 0.);

		void setArea(MgntROOT::Area area);
		void setLine(MgntROOT::Line line);
		void setMarker(MgntROOT::Marker marker);
		void setStyle(MgntROOT::Area area = MgntROOT::Area(), MgntROOT::Line line = MgntROOT::Line(), MgntROOT::Marker marker = MgntROOT::Marker());
		
		void draw(Option_t * option = "") 
			{ if (exist()) fFunc->Draw(option); }

		void draw(Option_t * option,  MgntROOT::Area area, MgntROOT::Line line, MgntROOT::Marker marker);
		
		inline void write(const std::string& name = "", Int_t option = 0, Int_t bufsize = 0) 
			{ if (exist()) fFunc->Write(name.c_str(), option, bufsize); }
	
		inline void clear()
			{ if (fFunc != nullptr) { fFuncMap.erase(fName); fFunc->Delete(); fFunc = nullptr; };
			  fUniqueID = -1; fName = ""; fTitle = ""; }
		
		inline std::string name()       { return fName; }
		inline std::string title()      { return fTitle; }
		inline std::string formula()    { return (exist() ? fFunc->GetTitle() : ""); }
		inline TF1 *       operator()() { return fFunc; }

		Func(const Func&) = delete;
		Func& operator= (const Func&) = delete;
	
	protected :
		inline Bool_t isFuncExist(const std::string& name);
		inline void pushBack() 
			{ fFuncMap[fName] = this; fUniqueID = fCounter++;  (*this)()->SetUniqueID(fUniqueID); }

	protected :
		Long64_t    fUniqueID;
		std::string fName;
		std::string fTitle;
		TF1 *       fFunc;
};

Long64_t                               MgntROOT::Func::fCounter = 1;
std::map<std::string, MgntROOT::Func*> MgntROOT::Func::fFuncMap;


#endif // __MgntROOT_H__
