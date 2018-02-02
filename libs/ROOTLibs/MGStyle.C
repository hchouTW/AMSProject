#ifndef __ROOTLibs_MGStyle_C__
#define __ROOTLibs_MGStyle_C__

#include "MGStyle.h"

namespace MGROOT {

// Marker
Style_t MarkerStyle(UInt_t idx) {
	const UInt_t NumStyle = 14;
	Style_t markerStyle = 20;
	Style_t style = (idx%NumStyle);
	switch (style) {
		// Circle (Full, Open)
		case  0 : markerStyle = 20; break;
		case  1 : markerStyle = 24; break;
		// Square (Full, Open)
		case  2 : markerStyle = 21; break;
		case  3 : markerStyle = 25; break;
		// TriangleUp (Full, Open)
		case  4 : markerStyle = 22; break;
		case  5 : markerStyle = 26; break;
		// TriangleDown (Full, Open)
		case  6 : markerStyle = 23; break;
		case  7 : markerStyle = 32; break;
		// Diamond (Full, Open)
		case  8 : markerStyle = 33; break;
		case  9 : markerStyle = 27; break;
		// Cross (Full, Open)
		case 10 : markerStyle = 34; break;
		case 11 : markerStyle = 28; break;
		// Star (Full, Open)
		case 12 : markerStyle = 29; break;
		case 13 : markerStyle = 30; break;
		default : break;
	}
	return markerStyle;
}
	
// Text
void TextDraw(const std::string& text, const TextStyle& style, const TextAlign& align) {
	TLatex ltx;
	ltx.SetNDC();
	ltx.SetTextColor(style.color);
	ltx.SetTextFont(style.font);
	ltx.SetTextSize(style.size);
	ltx.SetTextAlign(align.align);
	ltx.SetTextAngle(align.angle);
	ltx.DrawLatex(align.x, align.y, text.c_str());
}


// Window
Window::Window(WindowSize size) {
	switch(size) {
		case WindowSize::kSliceHR :    width = 2048; height = 1536; break;
		case WindowSize::kSliceLR :    width = 1024; height =  768; break;
		case WindowSize::kA4Vertical : width = 2480; height = 3508; break;
		case WindowSize::kA4Horizon :  width = 3508; height = 2480; break;
		case WindowSize::kMac :        width = 3840; height = 2400; break;
		case WindowSize::kMovie :      width = 4200; height = 1800; break;
		default : break;
	}
}


// Style
namespace Style {

void LoadDefaultEnvironment() {
	COUT("MGROOT::Style : Load default environment.\n");

    // Set Canvas
	Int_t col = 10;
	style.SetFrameBorderMode(col);
	style.SetFrameFillColor(col);
	style.SetCanvasBorderMode(col);
	style.SetCanvasColor(col);
	style.SetPadBorderMode(col);
	style.SetPadColor(col);
	style.SetStatColor(col);
	
	// Set the paper & margin sizes
	style.SetPaperSize(20, 26);
	      
	// Set margin sizes
	style.SetPadTopMargin(0.12);
	style.SetPadBottomMargin(0.12);
	style.SetPadRightMargin(0.15);
	style.SetPadLeftMargin(0.15);
	  
	// Use large fonts
	Int_t font = 42; // Helvetica
	
	style.SetTextFont(font);
	style.SetTextSize(0.03);
	
	style.SetAxisColor(1, "xyz");
	style.SetStripDecimals(false);	
	
	style.SetLabelColor(1, "xyz");
	style.SetLabelSize(0.03, "xyz");
	style.SetLabelFont(font, "xyz");
	style.SetLabelOffset(0.005, "xyz");
	
	style.SetTitleColor(1, "xyz"); 
	style.SetTitleSize(0.04, "xyz");
	style.SetTitleFont(font, "xyz");
	style.SetTitleOffset(1.05, "x");
	style.SetTitleOffset(1.20, "yz");
	
	// Legend
	style.SetLegendBorderSize(0);
	style.SetLegendFillColor(0);
	//style.SetLegendTextSize(0.025); // only for ROOT6 or later
	style.SetLegendFont(42);
	
	// Use bold lines and markers
	style.SetMarkerColor(1);
	style.SetMarkerStyle(20);
	style.SetMarkerSize(1.0);
	style.SetHistLineWidth((Width_t) 1.);
	style.SetLineStyleString(2, "[12 12]"); // postscript dashes
	
	// Set the hist title
	style.SetTitleAlign(22);
	style.SetTitleBorderSize(0);
	style.SetTitleFillColor(10);
	style.SetTitleFontSize(0.05);
	style.SetTitleX(0.5);  
	style.SetTitleY(0.95);
	
	// Do not display any of the standard histogram decorations
	style.SetOptTitle(0);
	style.SetOptStat(0);
	style.SetOptFit(0);
	
	// Set histogram x-axis error to zero
	style.SetErrorX(0);
	
	// Put tick marks on top and RHS of plots
	style.SetPadGridX(0);
	style.SetPadGridY(0);
	style.SetPadTickX(1);
	style.SetPadTickY(1);
	
	// Set axis scale
	style.SetOptLogx(0);
	style.SetOptLogy(0);
	style.SetOptLogz(1);
	
	// Set color palette
	SetColorPalette(-1);
	style.SetNumberContours(254);
	
	// Set windows
	Style::SetWindow();
	
	// Set event status
	if (style.GetShowEventStatus() == 0) style.ToggleEventStatus();
	
	// Update
	style.cd();
}


void SetColorPalette(Int_t pattern) {
    Printf("MGROOT::Style : Set New Color Palette.");
	if (pattern >= 0) {
		style.SetPalette(pattern, 0, 1);
	}
	else {
		const Int_t Num = 254;
		const Int_t NSet = 6;
		Int_t    palette[Num] = {0};
		Double_t red[NSet]    = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
		Double_t green[NSet]  = { 0.0, 1.0, 1.0, 1.0, 0.0, 0.0 };
		Double_t blue[NSet]   = { 1.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
		Double_t length[NSet] = { 0.0, .20, .40, .60, .80, 1.0 };
		Int_t FI = TColor::CreateGradientColorTable(NSet, length, red, green, blue, Num);
		for (Int_t i = 0; i < Num; ++i) palette[i] = FI + i;
		style.SetPalette(Num, palette, 1);
	}
	style.cd();
}


void SetWindow(const Window& window) {
	style.SetCanvasDefW(window.width+20);
	style.SetCanvasDefH(window.height+40);
	style.cd();
}


void SetWindow(WindowSize size) {
	Window window(size);
	SetWindow(window);
}


void SetAxis(Bool_t xaxis, Bool_t yaxis, Bool_t zaxis) {
	style.SetOptLogx(xaxis);
	style.SetOptLogy(yaxis);
	style.SetOptLogz(zaxis);
	style.cd();
}

} // namespace Style

} // namespace MGROOT

#endif // __ROOTLibs_MGStyle_C__
