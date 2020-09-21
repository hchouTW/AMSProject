#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "DataFit2D.h"
#include "DataFit2D.C"

void stdfmt(TH1* hist) {
    if (hist == nullptr) return;
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleFont(43);
    hist->GetXaxis()->SetTitleSize(20);
    hist->GetXaxis()->SetTitleOffset(2.5);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelSize(15);
    hist->GetYaxis()->SetTitleFont(43);
    hist->GetYaxis()->SetTitleSize(20);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelSize(15);
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    
    TGraphAsymmErrors* offap = (TGraphAsymmErrors*) (TFile::Open("others/ams02_ap2pr.root")->Get("gr_exp1"));
    
    Hist* hLPcnt = Hist::New("hLPcnt", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hPcnt"));
    Hist* hLNcnt = Hist::New("hLNcnt", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hNcnt"));
    Hist* hLStat = Hist::New("hLStat", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hStat"));
    Hist* hLSyst = Hist::New("hLSyst", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hSyst"));
    Hist* hLRate = Hist::New("hLRate", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hRate"));
    Hist* hLCrrR = Hist::New("hLCrrR", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hCrrR"));

    Hist* hMPcnt = Hist::New("hMPcnt", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hPcnt"));
    Hist* hMNcnt = Hist::New("hMNcnt", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hNcnt"));
    Hist* hMStat = Hist::New("hMStat", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hStat"));
    Hist* hMSyst = Hist::New("hMSyst", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hSyst"));
    Hist* hMRate = Hist::New("hMRate", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hRate"));
    Hist* hMCrrR = Hist::New("hMCrrR", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hCrrR"));

    Hist* hIPcnt = Hist::New("hIPcnt", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hPcnt"));
    Hist* hINcnt = Hist::New("hINcnt", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hNcnt"));
    Hist* hIStat = Hist::New("hIStat", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hStat"));
    Hist* hISyst = Hist::New("hISyst", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hSyst"));
    Hist* hIRate = Hist::New("hIRate", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hRate"));
    Hist* hICrrR = Hist::New("hICrrR", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hCrrR"));

    Hist* hHiePcnt = Hist::New("hHiePcnt", (TH1*)TFile::Open("out/apflux_iex.root")->Get("hPcnt"));
    Hist* hHieNcnt = Hist::New("hHieNcnt", (TH1*)TFile::Open("out/apflux_iex.root")->Get("hNcnt"));
    Hist* hHieStat = Hist::New("hHieStat", (TH1*)TFile::Open("out/apflux_iex.root")->Get("hStat"));
    Hist* hHieSyst = Hist::New("hHieSyst", (TH1*)TFile::Open("out/apflux_iex.root")->Get("hSyst"));
    Hist* hHieRate = Hist::New("hHieRate", (TH1*)TFile::Open("out/apflux_iex.root")->Get("hRate"));
    Hist* hHieCrrR = Hist::New("hHieCrrR", (TH1*)TFile::Open("out/apflux_iex.root")->Get("hCrrR"));
    
    Hist* hHexPcnt = Hist::New("hHexPcnt", (TH1*)TFile::Open("out/apflux_hex.root")->Get("hPcnt"));
    Hist* hHexNcnt = Hist::New("hHexNcnt", (TH1*)TFile::Open("out/apflux_hex.root")->Get("hNcnt"));
    Hist* hHexStat = Hist::New("hHexStat", (TH1*)TFile::Open("out/apflux_hex.root")->Get("hStat"));
    Hist* hHexSyst = Hist::New("hHexSyst", (TH1*)TFile::Open("out/apflux_hex.root")->Get("hSyst"));
    Hist* hHexRate = Hist::New("hHexRate", (TH1*)TFile::Open("out/apflux_hex.root")->Get("hRate"));
    Hist* hHexCrrR = Hist::New("hHexCrrR", (TH1*)TFile::Open("out/apflux_hex.root")->Get("hCrrR"));
    
    Hist* hHfsPcnt = Hist::New("hHfsPcnt", (TH1*)TFile::Open("out/apflux_hfs.root")->Get("hPcnt"));
    Hist* hHfsNcnt = Hist::New("hHfsNcnt", (TH1*)TFile::Open("out/apflux_hfs.root")->Get("hNcnt"));
    Hist* hHfsStat = Hist::New("hHfsStat", (TH1*)TFile::Open("out/apflux_hfs.root")->Get("hStat"));
    Hist* hHfsSyst = Hist::New("hHfsSyst", (TH1*)TFile::Open("out/apflux_hfs.root")->Get("hSyst"));
    Hist* hHfsRate = Hist::New("hHfsRate", (TH1*)TFile::Open("out/apflux_hfs.root")->Get("hRate"));
    Hist* hHfsCrrR = Hist::New("hHfsCrrR", (TH1*)TFile::Open("out/apflux_hfs.root")->Get("hCrrR"));
    
    Hist* hTLPcnt = Hist::New("hTLPcnt", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hTPcnt"));
    Hist* hTLNcnt = Hist::New("hTLNcnt", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hTNcnt"));
    Hist* hTLStat = Hist::New("hTLStat", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hTStat"));
    Hist* hTLSyst = Hist::New("hTLSyst", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hTSyst"));
    Hist* hTLRate = Hist::New("hTLRate", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hTRate"));
    Hist* hTLCrrR = Hist::New("hTLCrrR", (TH1*)TFile::Open("out/apflux_ltf.root")->Get("hTCrrR"));

    Hist* hTMPcnt = Hist::New("hTMPcnt", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hTPcnt"));
    Hist* hTMNcnt = Hist::New("hTMNcnt", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hTNcnt"));
    Hist* hTMStat = Hist::New("hTMStat", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hTStat"));
    Hist* hTMSyst = Hist::New("hTMSyst", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hTSyst"));
    Hist* hTMRate = Hist::New("hTMRate", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hTRate"));
    Hist* hTMCrrR = Hist::New("hTMCrrR", (TH1*)TFile::Open("out/apflux_lrh.root")->Get("hTCrrR"));
    
    Hist* hTIPcnt = Hist::New("hTIPcnt", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hTPcnt"));
    Hist* hTINcnt = Hist::New("hTINcnt", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hTNcnt"));
    Hist* hTIStat = Hist::New("hTIStat", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hTStat"));
    Hist* hTISyst = Hist::New("hTISyst", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hTSyst"));
    Hist* hTIRate = Hist::New("hTIRate", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hTRate"));
    Hist* hTICrrR = Hist::New("hTICrrR", (TH1*)TFile::Open("out/apflux_iin.root")->Get("hTCrrR"));

    PdfEditor editor(Window(), "apflux_all", "out");

    const Axis& AXrig = Hist::Head("hLCrrR")->xaxis();
    
    const Axis& AXTtme = Hist::Head("hTLCrrR")->xaxis();
    const Axis& AXTrig = Hist::Head("hTLCrrR")->yaxis();
 
    Hist* hPcnt = Hist::New("hPcnt", HistAxis(AXrig));
    Hist* hNcnt = Hist::New("hNcnt", HistAxis(AXrig));
    Hist* hStat = Hist::New("hStat", HistAxis(AXrig));
    Hist* hSyst = Hist::New("hSyst", HistAxis(AXrig));
    Hist* hRate = Hist::New("hRate", HistAxis(AXrig));
    Hist* hCrrR = Hist::New("hCrrR", HistAxis(AXrig));
    
    for (int ir = 1; ir <= 14; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        
        (*hPcnt)()->SetBinContent(ir, (*hLPcnt)()->GetBinContent(ir)); 
        (*hNcnt)()->SetBinContent(ir, (*hLNcnt)()->GetBinContent(ir)); 
                                      
        (*hStat)()->SetBinContent(ir, (*hLStat)()->GetBinContent(ir)); 
        (*hSyst)()->SetBinContent(ir, (*hLSyst)()->GetBinContent(ir));
        
        (*hRate)()->SetBinContent(ir, (*hLRate)()->GetBinContent(ir)); 
        (*hRate)()->SetBinError  (ir, (*hLRate)()->GetBinError  (ir)); 
                                      
        (*hCrrR)()->SetBinContent(ir, (*hLCrrR)()->GetBinContent(ir)); 
        (*hCrrR)()->SetBinError  (ir, (*hLCrrR)()->GetBinError  (ir)); 
    }
    
    for (int ir = 15; ir <= 24; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        
        (*hPcnt)()->SetBinContent(ir, (*hMPcnt)()->GetBinContent(ir)); 
        (*hNcnt)()->SetBinContent(ir, (*hMNcnt)()->GetBinContent(ir)); 
                                      
        (*hStat)()->SetBinContent(ir, (*hMStat)()->GetBinContent(ir)); 
        (*hSyst)()->SetBinContent(ir, (*hMSyst)()->GetBinContent(ir));
        
        (*hRate)()->SetBinContent(ir, (*hMRate)()->GetBinContent(ir)); 
        (*hRate)()->SetBinError  (ir, (*hMRate)()->GetBinError  (ir)); 
                                      
        (*hCrrR)()->SetBinContent(ir, (*hMCrrR)()->GetBinContent(ir)); 
        (*hCrrR)()->SetBinError  (ir, (*hMCrrR)()->GetBinError  (ir)); 
    }
    
    //for (int ir = 25; ir <= 37; ++ir) {
    for (int ir = 25; ir <= 39; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        
        (*hPcnt)()->SetBinContent(ir, (*hIPcnt)()->GetBinContent(ir)); 
        (*hNcnt)()->SetBinContent(ir, (*hINcnt)()->GetBinContent(ir)); 
                                      
        (*hStat)()->SetBinContent(ir, (*hIStat)()->GetBinContent(ir)); 
        (*hSyst)()->SetBinContent(ir, (*hISyst)()->GetBinContent(ir));
        
        (*hRate)()->SetBinContent(ir, (*hIRate)()->GetBinContent(ir)); 
        (*hRate)()->SetBinError  (ir, (*hIRate)()->GetBinError  (ir)); 
                                      
        (*hCrrR)()->SetBinContent(ir, (*hICrrR)()->GetBinContent(ir)); 
        (*hCrrR)()->SetBinError  (ir, (*hICrrR)()->GetBinError  (ir)); 
    }
    
    //for (int ir = 38; ir <= AXrig.nbin()-6; ++ir) {
    //for (int ir = 40; ir <= AXrig.nbin()-6; ++ir) {
    for (int ir = 30; ir <= AXrig.nbin()-12; ++ir) {
    //for (int ir = 40; ir <= AXrig.nbin()-12; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        
        (*hPcnt)()->SetBinContent(ir, (*hHiePcnt)()->GetBinContent(ir)); 
        (*hNcnt)()->SetBinContent(ir, (*hHieNcnt)()->GetBinContent(ir)); 
                                      
        (*hStat)()->SetBinContent(ir, (*hHieStat)()->GetBinContent(ir)); 
        (*hSyst)()->SetBinContent(ir, (*hHieSyst)()->GetBinContent(ir));
        
        (*hRate)()->SetBinContent(ir, (*hHieRate)()->GetBinContent(ir)); 
        (*hRate)()->SetBinError  (ir, (*hHieRate)()->GetBinError  (ir)); 
                                      
        (*hCrrR)()->SetBinContent(ir, (*hHieCrrR)()->GetBinContent(ir)); 
        (*hCrrR)()->SetBinError  (ir, (*hHieCrrR)()->GetBinError  (ir)); 
    }
    
    //for (int ir = 38; ir <= AXrig.nbin()-6; ++ir) {
    //for (int ir = 40; ir <= AXrig.nbin()-6; ++ir) {
    //for (int ir = 30; ir <= AXrig.nbin()-6; ++ir) {
    //for (int ir = 50; ir <= AXrig.nbin()-6; ++ir) {
    for (int ir = AXrig.nbin()-16; ir <= AXrig.nbin()-6; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        
        (*hPcnt)()->SetBinContent(ir, (*hHexPcnt)()->GetBinContent(ir)); 
        (*hNcnt)()->SetBinContent(ir, (*hHexNcnt)()->GetBinContent(ir)); 
                                      
        (*hStat)()->SetBinContent(ir, (*hHexStat)()->GetBinContent(ir)); 
        (*hSyst)()->SetBinContent(ir, (*hHexSyst)()->GetBinContent(ir));
        
        (*hRate)()->SetBinContent(ir, (*hHexRate)()->GetBinContent(ir)); 
        (*hRate)()->SetBinError  (ir, (*hHexRate)()->GetBinError  (ir)); 
                                      
        (*hCrrR)()->SetBinContent(ir, (*hHexCrrR)()->GetBinContent(ir)); 
        (*hCrrR)()->SetBinError  (ir, (*hHexCrrR)()->GetBinError  (ir)); 
    }

    //for (int ir = 38; ir <= AXrig.nbin()-3; ++ir) {
    for (int ir = AXrig.nbin()-5; ir <= AXrig.nbin()-2; ++ir) {
    //for (int ir = AXrig.nbin()-16; ir <= AXrig.nbin()-2; ++ir) {
    //for (int ir = AXrig.nbin()-11; ir <= AXrig.nbin()-2; ++ir) {
    //for (int ir = 30; ir <= AXrig.nbin()-2; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        
        (*hPcnt)()->SetBinContent(ir, (*hHfsPcnt)()->GetBinContent(ir)); 
        (*hNcnt)()->SetBinContent(ir, (*hHfsNcnt)()->GetBinContent(ir)); 
                                      
        (*hStat)()->SetBinContent(ir, (*hHfsStat)()->GetBinContent(ir)); 
        (*hSyst)()->SetBinContent(ir, (*hHfsSyst)()->GetBinContent(ir));
        
        (*hRate)()->SetBinContent(ir, (*hHfsRate)()->GetBinContent(ir)); 
        (*hRate)()->SetBinError  (ir, (*hHfsRate)()->GetBinError  (ir)); 
                                      
        (*hCrrR)()->SetBinContent(ir, (*hHfsCrrR)()->GetBinContent(ir)); 
        (*hCrrR)()->SetBinError  (ir, (*hHfsCrrR)()->GetBinError  (ir)); 
    }
    
    Hist* hTPcnt = Hist::New("hTPcnt", HistAxis(AXTtme, AXTrig));
    Hist* hTNcnt = Hist::New("hTNcnt", HistAxis(AXTtme, AXTrig));
    Hist* hTStat = Hist::New("hTStat", HistAxis(AXTtme, AXTrig));
    Hist* hTSyst = Hist::New("hTSyst", HistAxis(AXTtme, AXTrig));
    Hist* hTRate = Hist::New("hTRate", HistAxis(AXTtme, AXTrig));
    Hist* hTCrrR = Hist::New("hTCrrR", HistAxis(AXTtme, AXTrig));
    
    (*hTRate)()->GetXaxis()->SetTitle("Date");
    (*hTRate)()->GetXaxis()->SetTimeDisplay(true);
    (*hTRate)()->GetXaxis()->SetTimeOffset(0, "GMT");
    (*hTRate)()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
    (*hTRate)()->GetZaxis()->SetTitle("#bar{p}/p");
    
    (*hTCrrR)()->GetXaxis()->SetTitle("Date");
    (*hTCrrR)()->GetXaxis()->SetTimeDisplay(true);
    (*hTCrrR)()->GetXaxis()->SetTimeOffset(0, "GMT");
    (*hTCrrR)()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
    (*hTCrrR)()->GetZaxis()->SetTitle("#bar{p}/p");
    
    for (int it = 1; it <= AXTtme.nbin(); ++it) {
    for (int ir = 1; ir <= 5; ++ir) {
        //if (it == 16) continue; // TTCS Off 3m-binning
        if (it ==  48) continue; // 1m-binning
        if (it ==  99) continue; // 1m-binning
        if (it == 101) continue; // 1m-binning
        double tme = AXTtme.center(it, AxisScale::kLinear);
        double rig = AXTrig.center(ir, AxisScale::kLog);
        
        (*hTPcnt)()->SetBinContent(it, ir, (*hTLPcnt)()->GetBinContent(it, ir)); 
        (*hTNcnt)()->SetBinContent(it, ir, (*hTLNcnt)()->GetBinContent(it, ir)); 
        
        (*hTStat)()->SetBinContent(it, ir, (*hTLStat)()->GetBinContent(it, ir)); 
        (*hTSyst)()->SetBinContent(it, ir, (*hTLSyst)()->GetBinContent(it, ir));
        
        (*hTRate)()->SetBinContent(it, ir, (*hTLRate)()->GetBinContent(it, ir)); 
        (*hTRate)()->SetBinError  (it, ir, (*hTLRate)()->GetBinError  (it, ir)); 
        
        (*hTCrrR)()->SetBinContent(it, ir, (*hTLCrrR)()->GetBinContent(it, ir)); 
        (*hTCrrR)()->SetBinError  (it, ir, (*hTLCrrR)()->GetBinError  (it, ir)); 
    }}
    
    for (int it = 1; it <= AXTtme.nbin(); ++it) {
    for (int ir = 6; ir <= 13; ++ir) {
        //if (it == 16) continue; // TTCS Off 3m-binning
        if (it ==  48) continue; // 1m-binning
        if (it ==  99) continue; // 1m-binning
        if (it == 101) continue; // 1m-binning
        double tme = AXTtme.center(it, AxisScale::kLinear);
        double rig = AXTrig.center(ir, AxisScale::kLog);
        
        (*hTPcnt)()->SetBinContent(it, ir, (*hTMPcnt)()->GetBinContent(it, ir)); 
        (*hTNcnt)()->SetBinContent(it, ir, (*hTMNcnt)()->GetBinContent(it, ir)); 
        
        (*hTStat)()->SetBinContent(it, ir, (*hTMStat)()->GetBinContent(it, ir)); 
        (*hTSyst)()->SetBinContent(it, ir, (*hTMSyst)()->GetBinContent(it, ir));
        
        (*hTRate)()->SetBinContent(it, ir, (*hTMRate)()->GetBinContent(it, ir)); 
        (*hTRate)()->SetBinError  (it, ir, (*hTMRate)()->GetBinError  (it, ir)); 
        
        (*hTCrrR)()->SetBinContent(it, ir, (*hTMCrrR)()->GetBinContent(it, ir)); 
        (*hTCrrR)()->SetBinError  (it, ir, (*hTMCrrR)()->GetBinError  (it, ir)); 
    }}
    
    for (int it = 1; it <= AXTtme.nbin(); ++it) {
    for (int ir = 14; ir <= AXTrig.nbin(); ++ir) {
        //if (it == 16) continue; // TTCS Off 3m-binning
        if (it ==  48) continue; // 1m-binning
        if (it ==  99) continue; // 1m-binning
        if (it == 101) continue; // 1m-binning
        double tme = AXTtme.center(it, AxisScale::kLinear);
        double rig = AXTrig.center(ir, AxisScale::kLog);
        
        (*hTPcnt)()->SetBinContent(it, ir, (*hTIPcnt)()->GetBinContent(it, ir)); 
        (*hTNcnt)()->SetBinContent(it, ir, (*hTINcnt)()->GetBinContent(it, ir)); 
        
        (*hTStat)()->SetBinContent(it, ir, (*hTIStat)()->GetBinContent(it, ir)); 
        (*hTSyst)()->SetBinContent(it, ir, (*hTISyst)()->GetBinContent(it, ir));
        
        (*hTRate)()->SetBinContent(it, ir, (*hTIRate)()->GetBinContent(it, ir)); 
        (*hTRate)()->SetBinError  (it, ir, (*hTIRate)()->GetBinError  (it, ir)); 
        
        (*hTCrrR)()->SetBinContent(it, ir, (*hTICrrR)()->GetBinContent(it, ir)); 
        (*hTCrrR)()->SetBinError  (it, ir, (*hTICrrR)()->GetBinError  (it, ir)); 
    }}
   

    editor.create();
    editor.cd(0, PadAxis(1, 0));
    hCrrR->style(Line(kRed, 0, 1.5), Marker(kRed, MarkerStyle(MarkerShape::kCross)));
    //(*hCrrR)()->GetXaxis()->SetRangeUser(0.5, 800);
    (*hCrrR)()->GetXaxis()->SetRangeUser(0.5, 500);
    //(*hCrrR)()->GetXaxis()->SetRangeUser(0.5, 33);
    (*hCrrR)()->GetXaxis()->SetMoreLogLabels();
    (*hCrrR)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hCrrR)()->GetYaxis()->SetTitle("#bar{p}/p");
    (*hCrrR)()->Draw("pe");
    offap->Draw("p");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    hCrrR->style(Line(kRed, 0, 1.5), Marker(kRed, MarkerStyle(MarkerShape::kCross)));
    //(*hCrrR)()->GetXaxis()->SetRangeUser(0.5, 800);
    (*hCrrR)()->GetXaxis()->SetRangeUser(0.5, 500);
    //(*hCrrR)()->GetXaxis()->SetRangeUser(0.5, 33);
    (*hCrrR)()->GetXaxis()->SetMoreLogLabels();
    (*hCrrR)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
    (*hCrrR)()->GetYaxis()->SetTitle("#bar{p}/p");
    (*hCrrR)()->Draw("pe");
    offap->Draw("p");
    editor.save();
  
    Hist* hRatio = Hist::New("hh_Ratio_StatErr", (TH2D*)(TFile::Open("others/AntiprotonFlux_AntiPBin_3Months_new.root")->Get("hh_Ratio_StatErr")));
    std::vector<Hist*> vhR1 = Hist::ProjectAll(HistProj::kX, Hist::Head("hh_Ratio_StatErr"));
    std::vector<Hist*> vhR2 = Hist::ProjectAll(HistProj::kX, Hist::Head("hTCrrR"));

    for (int ir = 1; ir <= AXTrig.nbin(); ++ir) {
        vhR1.at(ir+1)->style(Line(kBlue, 0, 1.5),  Marker(kBlue,  MarkerStyle(MarkerShape::kCircle  )));
        vhR2.at(ir)->style(Line(kRed, 0, 1.5),  Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));

        THStack* stack = Hist::Collect(Form("hR%03d", ir), HistList({ vhR2.at(ir), vhR1.at(ir+1) }));
        editor.create();
        editor.cd(0, PadAxis(0, 0));
        stack->Draw("nostack hist");
        stack->GetHistogram()->SetLineColor(0);
        stack->GetHistogram()->GetXaxis()->SetTitle("Date");
        stack->GetHistogram()->GetXaxis()->SetTimeDisplay(true);
        stack->GetHistogram()->GetXaxis()->SetTimeOffset(0, "GMT");
        stack->GetHistogram()->GetXaxis()->SetTimeFormat("%y/%m/%d");
        stack->GetHistogram()->GetYaxis()->SetTitle("#bar{p}/p");
        stack->Draw("nostack pe");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.20, 0.45, 0.15, 0.25));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXTrig()(ir-1), AXTrig()(ir)));
        leg_table()->AddEntry((*vhR1.at(ir+1))(), "IHEP (pass6)", "lp");
        leg_table()->AddEntry((*vhR2.at(ir))(), "Result (pass7)", "lp");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save(); 
    }
    
    editor.close();

    TFile * ofle = new TFile("out/apflux_all.root", "RECREATE");
    ofle->cd();

    (*hPcnt)()->Write();
    (*hNcnt)()->Write();
    (*hStat)()->Write();
    (*hSyst)()->Write();
    (*hRate)()->Write();
    (*hCrrR)()->Write();
    
    (*hTPcnt)()->Write();
    (*hTNcnt)()->Write();
    (*hTStat)()->Write();
    (*hTSyst)()->Write();
    (*hTRate)()->Write();
    (*hTCrrR)()->Write();

    offap->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
