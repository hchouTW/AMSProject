#ifndef __ROOTLibs_MGStyle_H__
#define __ROOTLibs_MGStyle_H__

#include <TColor.h>
#include <TAttFill.h>
#include <TAttLine.h>
#include <TAttMarker.h>
#include <TLatex.h>
#include <TStyle.h>


namespace MGROOT {

//---- Color ----//
// (Color_t) kBlack, kBlue, kMagenta, kRed, kGreen, kCyan, kAzure, kViolet, kPink, kOrange, kSpring, kTeal, kYellow
inline Color_t Color(Color_t col = kBlack) { return col; }
inline Color_t Color(UInt_t idx, UInt_t nset = TColor::GetNumberOfColors()) { return TColor::GetColorPalette( idx * (TColor::GetNumberOfColors()/nset) ); }


//---- Fill ----//
inline TAttFill Fill(Color_t color = kBlack, Style_t style = 0) { return TAttFill(color, style); }


//---- Line ----//
inline TAttLine Line(Color_t color = kBlack, Style_t style = 0, Width_t width = 1) { return TAttLine(color, (style%10), width); }


//---- Marker ----//
enum class MarkerShape : UInt_t { kCircle = 0, kSquare = 1, kTriangleUp = 2, kTriangleDown = 3, kDiamond = 4, kCross = 5, kStar = 6 };
enum class MarkerType : UInt_t { kFull = 0, kOpen = 1 };

inline Style_t MarkerStyle(UInt_t idx = 0);
inline Style_t MarkerStyle(MarkerShape shape, MarkerType type = MarkerType::kFull) { return MarkerStyle(static_cast<UInt_t>(shape) * 2 + static_cast<UInt_t>(type)); }
	
inline TAttMarker Marker(Color_t color = kBlack, Style_t style = MarkerStyle(), Size_t size = 1.0) { return TAttMarker(color, style, size); } 


//---- Text ----//
struct TextStyle {
	TextStyle(Color_t _color = kBlack, Double_t _size = 0.035, Int_t _font = 42) : color(_color), size(_size), font(_font) {}
	Color_t     color;
	Double_t    size;
	Int_t       font;
};
struct TextAlign {
	TextAlign(Double_t _x = 0.5, Double_t _y = 0.5, Int_t _align = 22, Double_t _angle = 0.0) : x(_x), y(_y), align(_align), angle(_angle) {}
	Double_t x;
	Double_t y;
	Int_t    align;
	Double_t angle;
};
struct Text {
    Text(const std::string& _text = "", const TextStyle& _style = TextStyle()) : text(_text), style(_style) {}
    Text(const std::string& _text, Color_t _color, Double_t _size = 0.035, Int_t _font = 42) : text(_text), style(_color, _size, _font) {}
    std::string text;
    TextStyle   style;
};

inline void TextDraw(const std::string& text, const TextStyle& style = TextStyle(), const TextAlign& align = TextAlign());
inline void TextDraw(const std::vector<Text>& text, const TextAlign& align = TextAlign(0.15, 0.85, 12));
inline void TitleDraw(const std::string& text) { TextDraw(text, TextStyle(kBlack, 0.04), TextAlign(0.5, 0.97)); }

using TextList = std::vector<Text>;


//---- Window ----//
enum class WindowSize : Int_t { kSliceHR, kSliceMR, kSliceLR, kWideSliceHR, kWideSliceMR, kWideSliceLR, kNarrowSliceHR, kNarrowSliceMR, kNarrowSliceLR, kA4Vertical, kA4Horizon, kMac, kMovieHR, kMovieMR, kMovieLR };
struct Window {
	Window(UInt_t _width = 800, UInt_t _height = 600) : width(_width), height(_height) {}
	Window(WindowSize size);
	UInt_t width, height;
};


//---- Style ----//
namespace Style {
	TStyle style("MainStyle", "MainStyle");
	
    TStyle& Head() { return style; }
	void    Update() { style.cd(); }
	void    LoadDefaultEnvironment();

	// kGreyScale, kRainBow, kBird
    void  SetColorPalette(Int_t pattern = -1);

	void SetWindow(const Window& window = Window());
	void SetWindow(WindowSize size);

	void SetAxis(Bool_t xaxis = false, Bool_t yaxis = false, Bool_t zaxis = true);
}


} // namespace MGROOT


#endif // __ROOTLibs_MGStyle_H__
