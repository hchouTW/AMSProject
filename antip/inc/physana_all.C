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
    
    TGraphAsymmErrors* offap = (TGraphAsymmErrors*) (TFile::Open("inc/database_plot.root")->Get("gr_exp1"));
    
    Hist* hLPcnt = Hist::New("hLPcnt", (TH1*)TFile::Open("out/antip_l.root")->Get("hPcnt"));
    Hist* hLNcnt = Hist::New("hLNcnt", (TH1*)TFile::Open("out/antip_l.root")->Get("hNcnt"));
    Hist* hLStat = Hist::New("hLStat", (TH1*)TFile::Open("out/antip_l.root")->Get("hStat"));
    Hist* hLSyst = Hist::New("hLSyst", (TH1*)TFile::Open("out/antip_l.root")->Get("hSyst"));
    Hist* hLRate = Hist::New("hLRate", (TH1*)TFile::Open("out/antip_l.root")->Get("hRate"));
    Hist* hLCrrR = Hist::New("hLCrrR", (TH1*)TFile::Open("out/antip_l.root")->Get("hCrrR"));

    Hist* hMPcnt = Hist::New("hMPcnt", (TH1*)TFile::Open("out/antip_m.root")->Get("hPcnt"));
    Hist* hMNcnt = Hist::New("hMNcnt", (TH1*)TFile::Open("out/antip_m.root")->Get("hNcnt"));
    Hist* hMStat = Hist::New("hMStat", (TH1*)TFile::Open("out/antip_m.root")->Get("hStat"));
    Hist* hMSyst = Hist::New("hMSyst", (TH1*)TFile::Open("out/antip_m.root")->Get("hSyst"));
    Hist* hMRate = Hist::New("hMRate", (TH1*)TFile::Open("out/antip_m.root")->Get("hRate"));
    Hist* hMCrrR = Hist::New("hMCrrR", (TH1*)TFile::Open("out/antip_m.root")->Get("hCrrR"));

    Hist* hIPcnt = Hist::New("hIPcnt", (TH1*)TFile::Open("out/antip_i.root")->Get("hPcnt"));
    Hist* hINcnt = Hist::New("hINcnt", (TH1*)TFile::Open("out/antip_i.root")->Get("hNcnt"));
    Hist* hIStat = Hist::New("hIStat", (TH1*)TFile::Open("out/antip_i.root")->Get("hStat"));
    Hist* hISyst = Hist::New("hISyst", (TH1*)TFile::Open("out/antip_i.root")->Get("hSyst"));
    Hist* hIRate = Hist::New("hIRate", (TH1*)TFile::Open("out/antip_i.root")->Get("hRate"));
    Hist* hICrrR = Hist::New("hICrrR", (TH1*)TFile::Open("out/antip_i.root")->Get("hCrrR"));

    Hist* hHl1Pcnt = Hist::New("hHl1Pcnt", (TH1*)TFile::Open("out/antip_hl1.root")->Get("hPcnt"));
    Hist* hHl1Ncnt = Hist::New("hHl1Ncnt", (TH1*)TFile::Open("out/antip_hl1.root")->Get("hNcnt"));
    Hist* hHl1Stat = Hist::New("hHl1Stat", (TH1*)TFile::Open("out/antip_hl1.root")->Get("hStat"));
    Hist* hHl1Syst = Hist::New("hHl1Syst", (TH1*)TFile::Open("out/antip_hl1.root")->Get("hSyst"));
    Hist* hHl1Rate = Hist::New("hHl1Rate", (TH1*)TFile::Open("out/antip_hl1.root")->Get("hRate"));
    Hist* hHl1CrrR = Hist::New("hHl1CrrR", (TH1*)TFile::Open("out/antip_hl1.root")->Get("hCrrR"));
    
    Hist* hHl9Pcnt = Hist::New("hHl9Pcnt", (TH1*)TFile::Open("out/antip_hl9.root")->Get("hPcnt"));
    Hist* hHl9Ncnt = Hist::New("hHl9Ncnt", (TH1*)TFile::Open("out/antip_hl9.root")->Get("hNcnt"));
    Hist* hHl9Stat = Hist::New("hHl9Stat", (TH1*)TFile::Open("out/antip_hl9.root")->Get("hStat"));
    Hist* hHl9Syst = Hist::New("hHl9Syst", (TH1*)TFile::Open("out/antip_hl9.root")->Get("hSyst"));
    Hist* hHl9Rate = Hist::New("hHl9Rate", (TH1*)TFile::Open("out/antip_hl9.root")->Get("hRate"));
    Hist* hHl9CrrR = Hist::New("hHl9CrrR", (TH1*)TFile::Open("out/antip_hl9.root")->Get("hCrrR"));
    
    Hist* hHfsPcnt = Hist::New("hHfsPcnt", (TH1*)TFile::Open("out/antip_hfs.root")->Get("hPcnt"));
    Hist* hHfsNcnt = Hist::New("hHfsNcnt", (TH1*)TFile::Open("out/antip_hfs.root")->Get("hNcnt"));
    Hist* hHfsStat = Hist::New("hHfsStat", (TH1*)TFile::Open("out/antip_hfs.root")->Get("hStat"));
    Hist* hHfsSyst = Hist::New("hHfsSyst", (TH1*)TFile::Open("out/antip_hfs.root")->Get("hSyst"));
    Hist* hHfsRate = Hist::New("hHfsRate", (TH1*)TFile::Open("out/antip_hfs.root")->Get("hRate"));
    Hist* hHfsCrrR = Hist::New("hHfsCrrR", (TH1*)TFile::Open("out/antip_hfs.root")->Get("hCrrR"));
    
    Hist* hTLPcnt = Hist::New("hTLPcnt", (TH1*)TFile::Open("out/antip_l.root")->Get("hTPcnt"));
    Hist* hTLNcnt = Hist::New("hTLNcnt", (TH1*)TFile::Open("out/antip_l.root")->Get("hTNcnt"));
    Hist* hTLStat = Hist::New("hTLStat", (TH1*)TFile::Open("out/antip_l.root")->Get("hTStat"));
    Hist* hTLSyst = Hist::New("hTLSyst", (TH1*)TFile::Open("out/antip_l.root")->Get("hTSyst"));
    Hist* hTLRate = Hist::New("hTLRate", (TH1*)TFile::Open("out/antip_l.root")->Get("hTRate"));
    Hist* hTLCrrR = Hist::New("hTLCrrR", (TH1*)TFile::Open("out/antip_l.root")->Get("hTCrrR"));

    Hist* hTMPcnt = Hist::New("hTMPcnt", (TH1*)TFile::Open("out/antip_m.root")->Get("hTPcnt"));
    Hist* hTMNcnt = Hist::New("hTMNcnt", (TH1*)TFile::Open("out/antip_m.root")->Get("hTNcnt"));
    Hist* hTMStat = Hist::New("hTMStat", (TH1*)TFile::Open("out/antip_m.root")->Get("hTStat"));
    Hist* hTMSyst = Hist::New("hTMSyst", (TH1*)TFile::Open("out/antip_m.root")->Get("hTSyst"));
    Hist* hTMRate = Hist::New("hTMRate", (TH1*)TFile::Open("out/antip_m.root")->Get("hTRate"));
    Hist* hTMCrrR = Hist::New("hTMCrrR", (TH1*)TFile::Open("out/antip_m.root")->Get("hTCrrR"));
    
    Hist* hTIPcnt = Hist::New("hTIPcnt", (TH1*)TFile::Open("out/antip_i.root")->Get("hTPcnt"));
    Hist* hTINcnt = Hist::New("hTINcnt", (TH1*)TFile::Open("out/antip_i.root")->Get("hTNcnt"));
    Hist* hTIStat = Hist::New("hTIStat", (TH1*)TFile::Open("out/antip_i.root")->Get("hTStat"));
    Hist* hTISyst = Hist::New("hTISyst", (TH1*)TFile::Open("out/antip_i.root")->Get("hTSyst"));
    Hist* hTIRate = Hist::New("hTIRate", (TH1*)TFile::Open("out/antip_i.root")->Get("hTRate"));
    Hist* hTICrrR = Hist::New("hTICrrR", (TH1*)TFile::Open("out/antip_i.root")->Get("hTCrrR"));

    PdfEditor editor(Window(), "antip_all", "out");

    const Axis& AXrig = Hist::Head("hLCrrR")->xaxis();
    
    const Axis& AXTtme = Hist::Head("hTLCrrR")->xaxis();
    const Axis& AXTrig = Hist::Head("hTLCrrR")->yaxis();
 
    Hist* hPcnt = Hist::New("hPcnt", HistAxis(AXrig));
    Hist* hNcnt = Hist::New("hNcnt", HistAxis(AXrig));
    Hist* hStat = Hist::New("hStat", HistAxis(AXrig));
    Hist* hSyst = Hist::New("hSyst", HistAxis(AXrig));
    Hist* hRate = Hist::New("hRate", HistAxis(AXrig));
    Hist* hCrrR = Hist::New("hCrrR", HistAxis(AXrig));
    
    for (int ir = 1; ir <= 13; ++ir) {
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
    
    for (int ir = 14; ir <= 25; ++ir) {
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
    
    for (int ir = 26; ir <= 36; ++ir) {
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
    
    for (int ir = 37; ir <= 52; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        
        (*hPcnt)()->SetBinContent(ir, (*hHl1Pcnt)()->GetBinContent(ir) + (*hHl9Pcnt)()->GetBinContent(ir) + (*hHfsPcnt)()->GetBinContent(ir)); 
        (*hNcnt)()->SetBinContent(ir, (*hHl1Ncnt)()->GetBinContent(ir) + (*hHl9Ncnt)()->GetBinContent(ir) + (*hHfsNcnt)()->GetBinContent(ir)); 
                                     
        std::array<double, 3> rate({ (*hHl1CrrR)()->GetBinContent(ir), (*hHl9CrrR)()->GetBinContent(ir), (*hHfsCrrR)()->GetBinContent(ir) });
        std::array<double, 3> ierr({ 1.0/(*hHl1Stat)()->GetBinContent(ir), 1.0/(*hHl9Stat)()->GetBinContent(ir), 1.0/(*hHfsStat)()->GetBinContent(ir) });

        //double err = 1.0 / std::sqrt(ierr[0]*ierr[0] + ierr[1]*ierr[1] + ierr[2]*ierr[2]);
        //double val = (rate[0]*ierr[0]*ierr[0] + rate[1]*ierr[1]*ierr[1] + rate[2]*ierr[2]*ierr[2]) * (err*err);
        
        double err = 1.0 / std::sqrt(ierr[1]*ierr[1] + ierr[2]*ierr[2]);
        double val = (rate[1]*ierr[1]*ierr[1] + rate[2]*ierr[2]*ierr[2]) * (err*err);

        (*hStat)()->SetBinContent(ir, (*hHfsStat)()->GetBinContent(ir)); 
        (*hSyst)()->SetBinContent(ir, (*hHfsSyst)()->GetBinContent(ir));
        
        (*hRate)()->SetBinContent(ir, (*hHfsRate)()->GetBinContent(ir)); 
        (*hRate)()->SetBinError  (ir, (*hHfsRate)()->GetBinError  (ir)); 
                                      
        (*hCrrR)()->SetBinContent(ir, val); 
        (*hCrrR)()->SetBinError  (ir, err); 
    }

    for (int ir = 53; ir <= AXrig.nbin(); ++ir) {
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
        if (it == 16) continue; // TTCS Off
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
        if (it == 16) continue; // TTCS Off
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
        if (it == 16) continue; // TTCS Off
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
    (*hCrrR)()->GetXaxis()->SetRangeUser(1.0, 500);
    (*hCrrR)()->Draw("pe");
    offap->Draw("p");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(1, 1));
    hCrrR->style(Line(kRed, 0, 1.5), Marker(kRed, MarkerStyle(MarkerShape::kCross)));
    (*hCrrR)()->GetXaxis()->SetRangeUser(1.0, 500);
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
        
        Legend leg_table("", PadWindow(0.20, 0.45, 0.15, 0.25));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXTrig()(ir-1), AXTrig()(ir)));
        leg_table()->AddEntry((*vhR1.at(ir+1))(), "IHEP (pass6)", "lp");
        leg_table()->AddEntry((*vhR2.at(ir))(), "Result (pass7)", "lp");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save(); 
    }
    
    std::vector<Hist*> vhT1 = Hist::ProjectAll(HistProj::kY, Hist::Head("hh_Ratio_StatErr"));
    std::vector<Hist*> vhT2 = Hist::ProjectAll(HistProj::kY, Hist::Head("hTCrrR"));
    
    for (int it = 1; it <= AXTtme.nbin()-6; ++it) {
        vhT1.at(it)->style(Line(kBlue, 0, 1.5),  Marker(kBlue,  MarkerStyle(MarkerShape::kCircle  )));
        vhT2.at(it)->style(Line(kRed, 0, 1.5),  Marker(kRed,  MarkerStyle(MarkerShape::kCross  )));

        THStack* stack = Hist::Collect(Form("hT%03d", it), HistList({ vhT2.at(it), vhT1.at(it) }));
        editor.create();
        editor.cd(0, PadAxis(0, 0));
        stack->Draw("nostack hist");
        stack->GetHistogram()->SetLineColor(0);
        stack->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        //stack->GetHistogram()->GetXaxis()->SetTimeDisplay(true);
        //stack->GetHistogram()->GetXaxis()->SetTimeOffset(0, "GMT");
        //stack->GetHistogram()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
        stack->GetHistogram()->GetYaxis()->SetTitle("#bar{p}/p");
        stack->Draw("nostack pe");
        
        Legend leg_table("", PadWindow(0.55, 0.80, 0.15, 0.25));
        leg_table()->AddEntry((TObject*)nullptr, Form("%s", MGClock::ConvertFromUTimeToCTime(AXTtme()(it-1)).c_str()), "");
        leg_table()->AddEntry((TObject*)nullptr, Form("%s", MGClock::ConvertFromUTimeToCTime(AXTtme()(it)).c_str()), "");
        leg_table()->AddEntry((*vhT1.at(it))(), "IHEP (pass6)", "lp");
        leg_table()->AddEntry((*vhT2.at(it))(), "Result (pass7)", "lp");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save(); 
    }

    editor.close();

    TFile * ofle = new TFile("out/antip_all.root", "RECREATE");
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
