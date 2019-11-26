#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "TMultiGraph.h"

void style(TGraphErrors* gr, const TAttLine& line = MGROOT::Line(), const TAttMarker& marker = MGROOT::Marker(), const TAttFill& fill = MGROOT::Fill()) {
    if (gr == nullptr) return;
	gr->SetLineColor(line.GetLineColor());
	gr->SetLineStyle(line.GetLineStyle());
	gr->SetLineWidth(line.GetLineWidth());
	gr->SetMarkerColor(marker.GetMarkerColor());
	gr->SetMarkerStyle(marker.GetMarkerStyle());
	gr->SetMarkerSize (marker.GetMarkerSize());
	gr->SetFillColor(fill.GetFillColor());
	gr->SetFillStyle(fill.GetFillStyle());
}

void stdfmt(TMultiGraph* mg) {
    if (mg == nullptr) return;
    mg->GetHistogram()->GetXaxis()->CenterTitle();
    mg->GetHistogram()->GetXaxis()->SetTitleFont(43);
    mg->GetHistogram()->GetXaxis()->SetTitleSize(20);
    mg->GetHistogram()->GetXaxis()->SetTitleOffset(2.5);
    mg->GetHistogram()->GetXaxis()->SetLabelFont(43);
    mg->GetHistogram()->GetXaxis()->SetLabelSize(15);
    mg->GetHistogram()->GetYaxis()->SetTitleFont(43);
    mg->GetHistogram()->GetYaxis()->SetTitleSize(20);
    mg->GetHistogram()->GetYaxis()->SetLabelFont(43);
    mg->GetHistogram()->GetYaxis()->SetLabelSize(15);
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);

    if (argc != 3) return 0;
    std::string dir_pr = Form("/afs/cern.ch/user/h/hchou/AMSProject/trsys/%s", argv[1]);
    std::string dir_he = Form("/afs/cern.ch/user/h/hchou/AMSProject/trsys/%s", argv[2]);

    PdfEditor editor(Window(), "fit_merge", "out");

    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_mscat.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_MCmen = (TGraphErrors*) file_pr->Get("gr_hMCmen");
        TGraphErrors* gr_pr_SMmen = (TGraphErrors*) file_pr->Get("gr_hSMmen");
       
        style(gr_pr_MCmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_SMmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_mscat.root", dir_he.c_str()));
        TGraphErrors* gr_he_MCmen = (TGraphErrors*) file_he->Get("gr_hMCmen");
        TGraphErrors* gr_he_SMmen = (TGraphErrors*) file_he->Get("gr_hSMmen");
        
        style(gr_he_MCmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_SMmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_MCmen);
        mg_men->Add(gr_pr_SMmen);
        mg_men->Add(gr_he_MCmen);
        mg_men->Add(gr_he_SMmen);
        
        editor.create();
        editor.cd(1, PadAxis(1));
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->SetMoreLogLabels();
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Mean of Cos-Angle Difference");
        mg_men->GetHistogram()->SetMinimum(-0.0008);
        mg_men->GetHistogram()->SetMaximum( 0.0008);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.15, 0.40, 0.70, 0.85));
        leg_men()->AddEntry(gr_pr_MCmen, "Proton: GEANT4 Simulation", "lp");
        leg_men()->AddEntry(gr_pr_SMmen, "Proton: FPM Simulation", "lp");
        leg_men()->AddEntry(gr_he_MCmen, "Helium-4: GEANT4 Simulation", "lp");
        leg_men()->AddEntry(gr_he_SMmen, "Helium-4: FPM Simulation", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_mscat.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_MCrso = (TGraphErrors*) file_pr->Get("gr_hMCrso");
        TGraphErrors* gr_pr_SMrso = (TGraphErrors*) file_pr->Get("gr_hSMrso");
       
        style(gr_pr_MCrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_SMrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TGraphErrors* gr_pr_relrso = (TGraphErrors*) file_pr->Get("gr_hSMrelrso");
        style(gr_pr_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_mscat.root", dir_he.c_str()));
        TGraphErrors* gr_he_MCrso = (TGraphErrors*) file_he->Get("gr_hMCrso");
        TGraphErrors* gr_he_SMrso = (TGraphErrors*) file_he->Get("gr_hSMrso");
        
        style(gr_he_MCrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_SMrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        TGraphErrors* gr_he_relrso = (TGraphErrors*) file_he->Get("gr_hSMrelrso");
        style(gr_he_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_rso = new TMultiGraph();
        mg_rso->Add(gr_pr_MCrso);
        mg_rso->Add(gr_pr_SMrso);
        mg_rso->Add(gr_he_MCrso);
        mg_rso->Add(gr_he_SMrso);
        
        TMultiGraph* mg_relrso = new TMultiGraph();
        mg_relrso->Add(gr_pr_relrso);
        mg_relrso->Add(gr_he_relrso);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(1));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_rso->Draw("ap");
        mg_rso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_rso->GetHistogram()->SetLineColor(0);
        mg_rso->GetHistogram()->SetMarkerColor(0);
        mg_rso->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
        mg_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Cos-Angle Difference");
        stdfmt(mg_rso);
        mg_rso->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_rso->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_rso->Draw("ap");
        Legend leg_rso("", PadWindow(0.15, 0.40, 0.60, 0.85));
        leg_rso()->AddEntry(gr_pr_MCrso, "Proton: GEANT4 Simulation", "lp");
        leg_rso()->AddEntry(gr_pr_SMrso, "Proton: FPM Simulation", "lp");
        leg_rso()->AddEntry(gr_he_MCrso, "Helium-4: GEANT4 Simulation", "lp");
        leg_rso()->AddEntry(gr_he_SMrso, "Helium-4: FPM Simulation", "lp");
        leg_rso()->SetTextFont(43);
        leg_rso()->SetTextSize(15);
        leg_rso()->SetFillColor(0);
        leg_rso.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(1));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relrso->Draw("ap");
        mg_relrso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relrso->GetHistogram()->SetLineColor(0);
        mg_relrso->GetHistogram()->SetMarkerColor(0);
        mg_relrso->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
        mg_relrso->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relrso->GetHistogram()->SetMinimum(0.915);
        mg_relrso->GetHistogram()->SetMaximum(1.085);
        stdfmt(mg_relrso);
        mg_relrso->Draw("ap");
        Legend leg_relrso("", PadWindow(0.15, 0.40, 0.75, 0.95));
        leg_relrso()->AddEntry(gr_pr_relrso, "Proton: FPM / GEANT4", "lp");
        leg_relrso()->AddEntry(gr_he_relrso, "Helium-4: FPM / GEANT4", "lp");
        leg_relrso()->SetTextFont(43);
        leg_relrso()->SetTextSize(15);
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_eloss.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_MCmen = (TGraphErrors*) file_pr->Get("gr_hMCmen");
        TGraphErrors* gr_pr_SMmen = (TGraphErrors*) file_pr->Get("gr_hSMmen");
       
        style(gr_pr_MCmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_SMmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TGraphErrors* gr_pr_relmen = (TGraphErrors*) file_pr->Get("gr_hSMrelmen");
        style(gr_pr_relmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_eloss.root", dir_he.c_str()));
        TGraphErrors* gr_he_MCmen = (TGraphErrors*) file_he->Get("gr_hMCmen");
        TGraphErrors* gr_he_SMmen = (TGraphErrors*) file_he->Get("gr_hSMmen");
        
        style(gr_he_MCmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_SMmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        TGraphErrors* gr_he_relmen = (TGraphErrors*) file_he->Get("gr_hSMrelmen");
        style(gr_he_relmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_MCmen);
        mg_men->Add(gr_pr_SMmen);
        mg_men->Add(gr_he_MCmen);
        mg_men->Add(gr_he_SMmen);
        
        TMultiGraph* mg_relmen = new TMultiGraph();
        mg_relmen->Add(gr_pr_relmen);
        mg_relmen->Add(gr_he_relmen);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(1));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Energy Loss [MeV]");
        stdfmt(mg_men);
        mg_men->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_men->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.15, 0.40, 0.60, 0.85));
        leg_men()->AddEntry(gr_pr_MCmen, "Proton: GEANT4 Simulation", "lp");
        leg_men()->AddEntry(gr_pr_SMmen, "Proton: FPM Simulation", "lp");
        leg_men()->AddEntry(gr_he_MCmen, "Helium-4: GEANT4 Simulation", "lp");
        leg_men()->AddEntry(gr_he_SMmen, "Helium-4: FPM Simulation", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(1));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relmen->Draw("ap");
        mg_relmen->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relmen->GetHistogram()->SetLineColor(0);
        mg_relmen->GetHistogram()->SetMarkerColor(0);
        mg_relmen->GetHistogram()->GetXaxis()->SetTitle("1/(#gamma#beta)");
        mg_relmen->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relmen->GetHistogram()->SetMinimum(0.915);
        mg_relmen->GetHistogram()->SetMaximum(1.085);
        stdfmt(mg_relmen);
        mg_relmen->Draw("ap");
        Legend leg_relmen("", PadWindow(0.15, 0.40, 0.75, 0.95));
        leg_relmen()->AddEntry(gr_pr_relmen, "Proton: FPM / GEANT4", "lp");
        leg_relmen()->AddEntry(gr_he_relmen, "Helium-4: FPM / GEANT4", "lp");
        leg_relmen()->SetTextFont(43);
        leg_relmen()->SetTextSize(15);
        leg_relmen()->SetFillColor(0);
        leg_relmen.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_vel_tf.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get("gr_hOFmen");
        TGraphErrors* gr_pr_HCmen = (TGraphErrors*) file_pr->Get("gr_hTQmen");
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TFile* file_he = TFile::Open(Form("%s/fit_vel_tf.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get("gr_hOFmen");
        TGraphErrors* gr_he_HCmen = (TGraphErrors*) file_he->Get("gr_hTQmen");
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_HCmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_HCmen);
        
        editor.create();
        editor.cd(1, PadAxis(0));
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Bias of Beta [%]");
        mg_men->GetHistogram()->SetMinimum(-15.0);
        mg_men->GetHistogram()->SetMaximum(  3.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.50, 0.85, 0.20, 0.40));
        leg_men()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_HCmen, "Proton: FPM fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_HCmen, "Helium-4: FPM fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_vel_tf.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFrso = (TGraphErrors*) file_pr->Get("gr_hOFrso");
        TGraphErrors* gr_pr_HCrso = (TGraphErrors*) file_pr->Get("gr_hTQrso");
       
        style(gr_pr_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TGraphErrors* gr_pr_relrso = (TGraphErrors*) file_pr->Get("gr_hTQrelrso");
        style(gr_pr_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_vel_tf.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFrso = (TGraphErrors*) file_he->Get("gr_hOFrso");
        TGraphErrors* gr_he_HCrso = (TGraphErrors*) file_he->Get("gr_hTQrso");
        
        style(gr_he_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        TGraphErrors* gr_he_relrso = (TGraphErrors*) file_he->Get("gr_hTQrelrso");
        style(gr_he_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_rso = new TMultiGraph();
        mg_rso->Add(gr_pr_OFrso);
        mg_rso->Add(gr_pr_HCrso);
        mg_rso->Add(gr_he_OFrso);
        mg_rso->Add(gr_he_HCrso);
        
        TMultiGraph* mg_relrso = new TMultiGraph();
        mg_relrso->Add(gr_pr_relrso);
        mg_relrso->Add(gr_he_relrso);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(0));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_rso->Draw("ap");
        mg_rso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_rso->GetHistogram()->SetLineColor(0);
        mg_rso->GetHistogram()->SetMarkerColor(0);
        mg_rso->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Beta [%]");
        stdfmt(mg_rso);
        mg_rso->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_rso->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_rso->GetHistogram()->SetMinimum(0.0);
        mg_rso->GetHistogram()->SetMaximum(4.0);
        mg_rso->Draw("ap");
        Legend leg_rso("", PadWindow(0.15, 0.40, 0.60, 0.85));
        leg_rso()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_rso()->AddEntry(gr_pr_OFrso, "Proton: Official fitting", "lp");
        leg_rso()->AddEntry(gr_pr_HCrso, "Proton: FPM fitting", "lp");
        leg_rso()->AddEntry(gr_he_OFrso, "Helium-4: Official fitting", "lp");
        leg_rso()->AddEntry(gr_he_HCrso, "Helium-4: FPM fitting", "lp");
        leg_rso()->SetTextFont(43);
        leg_rso()->SetTextSize(15);
        leg_rso()->SetFillColor(0);
        leg_rso.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(0));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relrso->Draw("ap");
        mg_relrso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relrso->GetHistogram()->SetLineColor(0);
        mg_relrso->GetHistogram()->SetMarkerColor(0);
        mg_relrso->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_relrso->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relrso->GetHistogram()->SetMinimum(0.45);
        mg_relrso->GetHistogram()->SetMaximum(1.05);
        stdfmt(mg_relrso);
        mg_relrso->Draw("ap");
        Legend leg_relrso("", PadWindow(0.15, 0.40, 0.75, 0.95));
        leg_relrso()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_relrso()->AddEntry(gr_pr_relrso, "Proton: FPM / Official", "lp");
        leg_relrso()->AddEntry(gr_he_relrso, "Helium-4: FPM / Official", "lp");
        leg_relrso()->SetTextFont(43);
        leg_relrso()->SetTextSize(15);
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_vel_rh.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get("gr_hOFmen");
        TGraphErrors* gr_pr_HCmen = (TGraphErrors*) file_pr->Get("gr_hHCmen");
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TFile* file_he = TFile::Open(Form("%s/fit_vel_rh.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get("gr_hOFmen");
        TGraphErrors* gr_he_HCmen = (TGraphErrors*) file_he->Get("gr_hHCmen");
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_HCmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_HCmen);
        
        editor.create();
        editor.cd(1, PadAxis(0));
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Bias of Beta [%]");
        mg_men->GetHistogram()->SetMinimum(-0.15);
        mg_men->GetHistogram()->SetMaximum( 0.05);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.55, 0.85, 0.20, 0.40));
        leg_men()->AddEntry((TObject*)0, "RICH", "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_HCmen, "Proton: FPM fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_HCmen, "Helium-4: FPM fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_vel_rh.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFrso = (TGraphErrors*) file_pr->Get("gr_hOFrso");
        TGraphErrors* gr_pr_HCrso = (TGraphErrors*) file_pr->Get("gr_hHCrso");
       
        style(gr_pr_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TGraphErrors* gr_pr_relrso = (TGraphErrors*) file_pr->Get("gr_hHCrelrso");
        style(gr_pr_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_vel_rh.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFrso = (TGraphErrors*) file_he->Get("gr_hOFrso");
        TGraphErrors* gr_he_HCrso = (TGraphErrors*) file_he->Get("gr_hHCrso");
        
        style(gr_he_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        TGraphErrors* gr_he_relrso = (TGraphErrors*) file_he->Get("gr_hHCrelrso");
        style(gr_he_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_rso = new TMultiGraph();
        mg_rso->Add(gr_pr_OFrso);
        mg_rso->Add(gr_pr_HCrso);
        mg_rso->Add(gr_he_OFrso);
        mg_rso->Add(gr_he_HCrso);
        
        TMultiGraph* mg_relrso = new TMultiGraph();
        mg_relrso->Add(gr_pr_relrso);
        mg_relrso->Add(gr_he_relrso);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(0));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_rso->Draw("ap");
        mg_rso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_rso->GetHistogram()->SetLineColor(0);
        mg_rso->GetHistogram()->SetMarkerColor(0);
        mg_rso->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Beta [%]");
        stdfmt(mg_rso);
        mg_rso->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_rso->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_rso->GetHistogram()->SetMinimum(0.0);
        mg_rso->GetHistogram()->SetMaximum(0.18);
        mg_rso->Draw("ap");
        Legend leg_rso("", PadWindow(0.15, 0.40, 0.60, 0.85));
        leg_rso()->AddEntry((TObject*)0, "RICH", "");
        leg_rso()->AddEntry(gr_pr_OFrso, "Proton: Official fitting", "lp");
        leg_rso()->AddEntry(gr_pr_HCrso, "Proton: FPM fitting", "lp");
        leg_rso()->AddEntry(gr_he_OFrso, "Helium-4: Official fitting", "lp");
        leg_rso()->AddEntry(gr_he_HCrso, "Helium-4: FPM fitting", "lp");
        leg_rso()->SetTextFont(43);
        leg_rso()->SetTextSize(15);
        leg_rso()->SetFillColor(0);
        leg_rso.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(0));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relrso->Draw("ap");
        mg_relrso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relrso->GetHistogram()->SetLineColor(0);
        mg_relrso->GetHistogram()->SetMarkerColor(0);
        mg_relrso->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_relrso->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relrso->GetHistogram()->SetMinimum(0.85);
        mg_relrso->GetHistogram()->SetMaximum(1.05);
        stdfmt(mg_relrso);
        mg_relrso->Draw("ap");
        Legend leg_relrso("", PadWindow(0.15, 0.40, 0.75, 0.95));
        leg_relrso()->AddEntry((TObject*)0, "RICH", "");
        leg_relrso()->AddEntry(gr_pr_relrso, "Proton: FPM / Official", "lp");
        leg_relrso()->AddEntry(gr_he_relrso, "Helium-4: FPM / Official", "lp");
        leg_relrso()->SetTextFont(43);
        leg_relrso()->SetTextSize(15);
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_mutr_tf.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get("gr_hOFmen");
        TGraphErrors* gr_pr_HCmen = (TGraphErrors*) file_pr->Get("gr_hHCmen");
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TFile* file_he = TFile::Open(Form("%s/fit_mutr_tf.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get("gr_hOFmen");
        TGraphErrors* gr_he_HCmen = (TGraphErrors*) file_he->Get("gr_hHCmen");
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_HCmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_HCmen);
        
        editor.create();
        editor.cd(1, PadAxis(0));
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Mean of Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        mg_men->GetHistogram()->SetMinimum(0.0);
        mg_men->GetHistogram()->SetMaximum(4.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.15, 0.40, 0.40, 0.60));
        leg_men()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_HCmen, "Proton: FPM fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_HCmen, "Helium-4: FPM fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_mutr_tf.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get("gr_hOFmen");
        TGraphErrors* gr_pr_HCmen = (TGraphErrors*) file_pr->Get("gr_hHCmen");
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TGraphErrors* gr_pr_OFrelmen = (TGraphErrors*) file_pr->Get("gr_hOFrelmen");
        TGraphErrors* gr_pr_HCrelmen = (TGraphErrors*) file_pr->Get("gr_hHCrelmen");
       
        style(gr_pr_OFrelmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCrelmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_mutr_tf.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get("gr_hOFmen");
        TGraphErrors* gr_he_HCmen = (TGraphErrors*) file_he->Get("gr_hHCmen");
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        TGraphErrors* gr_he_OFrelmen = (TGraphErrors*) file_he->Get("gr_hOFrelmen");
        TGraphErrors* gr_he_HCrelmen = (TGraphErrors*) file_he->Get("gr_hHCrelmen");
        
        style(gr_he_OFrelmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCrelmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_HCmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_HCmen);
        
        TMultiGraph* mg_relmen = new TMultiGraph();
        mg_relmen->Add(gr_pr_OFrelmen);
        mg_relmen->Add(gr_pr_HCrelmen);
        mg_relmen->Add(gr_he_OFrelmen);
        mg_relmen->Add(gr_he_HCrelmen);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(0));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Mean of Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        stdfmt(mg_men);
        mg_men->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_men->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_men->GetHistogram()->SetMinimum(0.0);
        mg_men->GetHistogram()->SetMaximum(4.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.15, 0.40, 0.30, 0.60));
        leg_men()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_HCmen, "Proton: FPM fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_HCmen, "Helium-4: FPM fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(0));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relmen->Draw("ap");
        mg_relmen->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relmen->GetHistogram()->SetLineColor(0);
        mg_relmen->GetHistogram()->SetMarkerColor(0);
        mg_relmen->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_relmen->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relmen->GetHistogram()->SetMinimum(0.645);
        mg_relmen->GetHistogram()->SetMaximum(1.405);
        stdfmt(mg_relmen);
        mg_relmen->Draw("ap");
        Legend leg_relmen("", PadWindow(0.15, 0.40, 0.60, 0.95));
        leg_relmen()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_relmen()->AddEntry(gr_pr_OFrelmen, "Proton: Official / Theory", "lp");
        leg_relmen()->AddEntry(gr_pr_HCrelmen, "Proton: FPM / Theory", "lp");
        leg_relmen()->AddEntry(gr_he_OFrelmen, "Helium-4: Official / Theory", "lp");
        leg_relmen()->AddEntry(gr_he_HCrelmen, "Helium-4: FPM / Theory", "lp");
        leg_relmen()->SetTextFont(43);
        leg_relmen()->SetTextSize(15);
        leg_relmen()->SetFillColor(0);
        leg_relmen.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }

    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_mutr_tf.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFrso = (TGraphErrors*) file_pr->Get("gr_hOFrso");
        TGraphErrors* gr_pr_HCrso = (TGraphErrors*) file_pr->Get("gr_hHCrso");
       
        style(gr_pr_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TGraphErrors* gr_pr_relrso = (TGraphErrors*) file_pr->Get("gr_hHCrelrso");
        style(gr_pr_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_mutr_tf.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFrso = (TGraphErrors*) file_he->Get("gr_hOFrso");
        TGraphErrors* gr_he_HCrso = (TGraphErrors*) file_he->Get("gr_hHCrso");
        
        style(gr_he_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        TGraphErrors* gr_he_relrso = (TGraphErrors*) file_he->Get("gr_hHCrelrso");
        style(gr_he_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_rso = new TMultiGraph();
        mg_rso->Add(gr_pr_OFrso);
        mg_rso->Add(gr_pr_HCrso);
        mg_rso->Add(gr_he_OFrso);
        mg_rso->Add(gr_he_HCrso);
        
        TMultiGraph* mg_relrso = new TMultiGraph();
        mg_relrso->Add(gr_pr_relrso);
        mg_relrso->Add(gr_he_relrso);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(0));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_rso->Draw("ap");
        mg_rso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_rso->GetHistogram()->SetLineColor(0);
        mg_rso->GetHistogram()->SetMarkerColor(0);
        mg_rso->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        stdfmt(mg_rso);
        mg_rso->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_rso->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_rso->GetHistogram()->SetMinimum(0.0);
        mg_rso->GetHistogram()->SetMaximum(5.0);
        mg_rso->Draw("ap");
        Legend leg_rso("", PadWindow(0.15, 0.40, 0.50, 0.85));
        leg_rso()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_rso()->AddEntry(gr_pr_OFrso, "Proton: Official fitting", "lp");
        leg_rso()->AddEntry(gr_pr_HCrso, "Proton: FPM fitting", "lp");
        leg_rso()->AddEntry(gr_he_OFrso, "Helium-4: Official fitting", "lp");
        leg_rso()->AddEntry(gr_he_HCrso, "Helium-4: FPM fitting", "lp");
        leg_rso()->SetTextFont(43);
        leg_rso()->SetTextSize(15);
        leg_rso()->SetFillColor(0);
        leg_rso.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(0));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relrso->Draw("ap");
        mg_relrso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relrso->GetHistogram()->SetLineColor(0);
        mg_relrso->GetHistogram()->SetMarkerColor(0);
        mg_relrso->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_relrso->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relrso->GetHistogram()->SetMinimum(0.795);
        mg_relrso->GetHistogram()->SetMaximum(1.105);
        stdfmt(mg_relrso);
        mg_relrso->Draw("ap");
        Legend leg_relrso("", PadWindow(0.15, 0.40, 0.75, 0.95));
        leg_relrso()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_relrso()->AddEntry(gr_pr_relrso, "Proton: FPM / Official", "lp");
        leg_relrso()->AddEntry(gr_he_relrso, "Helium-4: FPM / Official", "lp");
        leg_relrso()->SetTextFont(43);
        leg_relrso()->SetTextSize(15);
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_mutr_rh.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get("gr_hOFmen");
        TGraphErrors* gr_pr_HCmen = (TGraphErrors*) file_pr->Get("gr_hHCmen");
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TFile* file_he = TFile::Open(Form("%s/fit_mutr_rh.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get("gr_hOFmen");
        TGraphErrors* gr_he_HCmen = (TGraphErrors*) file_he->Get("gr_hHCmen");
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_HCmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_HCmen);
        
        editor.create();
        editor.cd(1, PadAxis(0));
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Mean of Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        mg_men->GetHistogram()->SetMinimum(0.0);
        mg_men->GetHistogram()->SetMaximum(4.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.15, 0.40, 0.40, 0.60));
        leg_men()->AddEntry((TObject*)0, "Tracker L2-L8 & RICH", "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_HCmen, "Proton: FPM fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_HCmen, "Helium-4: FPM fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_mutr_rh.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get("gr_hOFmen");
        TGraphErrors* gr_pr_HCmen = (TGraphErrors*) file_pr->Get("gr_hHCmen");
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TGraphErrors* gr_pr_OFrelmen = (TGraphErrors*) file_pr->Get("gr_hOFrelmen");
        TGraphErrors* gr_pr_HCrelmen = (TGraphErrors*) file_pr->Get("gr_hHCrelmen");
       
        style(gr_pr_OFrelmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCrelmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_mutr_rh.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get("gr_hOFmen");
        TGraphErrors* gr_he_HCmen = (TGraphErrors*) file_he->Get("gr_hHCmen");
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        TGraphErrors* gr_he_OFrelmen = (TGraphErrors*) file_he->Get("gr_hOFrelmen");
        TGraphErrors* gr_he_HCrelmen = (TGraphErrors*) file_he->Get("gr_hHCrelmen");
        
        style(gr_he_OFrelmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCrelmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_HCmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_HCmen);
        
        TMultiGraph* mg_relmen = new TMultiGraph();
        mg_relmen->Add(gr_pr_OFrelmen);
        mg_relmen->Add(gr_pr_HCrelmen);
        mg_relmen->Add(gr_he_OFrelmen);
        mg_relmen->Add(gr_he_HCrelmen);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(0));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Mean of Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        stdfmt(mg_men);
        mg_men->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_men->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_men->GetHistogram()->SetMinimum(0.0);
        mg_men->GetHistogram()->SetMaximum(4.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.15, 0.40, 0.30, 0.60));
        leg_men()->AddEntry((TObject*)0, "RICH", "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_HCmen, "Proton: FPM fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_HCmen, "Helium-4: FPM fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(0));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relmen->Draw("ap");
        mg_relmen->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relmen->GetHistogram()->SetLineColor(0);
        mg_relmen->GetHistogram()->SetMarkerColor(0);
        mg_relmen->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_relmen->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relmen->GetHistogram()->SetMinimum(0.645);
        mg_relmen->GetHistogram()->SetMaximum(1.055);
        stdfmt(mg_relmen);
        mg_relmen->Draw("ap");
        Legend leg_relmen("", PadWindow(0.15, 0.40, 0.30, 0.75));
        leg_relmen()->AddEntry((TObject*)0, "RICH", "");
        leg_relmen()->AddEntry(gr_pr_OFrelmen, "Proton: Official / Theory", "lp");
        leg_relmen()->AddEntry(gr_pr_HCrelmen, "Proton: FPM / Theory", "lp");
        leg_relmen()->AddEntry(gr_he_OFrelmen, "Helium-4: Official / Theory", "lp");
        leg_relmen()->AddEntry(gr_he_HCrelmen, "Helium-4: FPM / Theory", "lp");
        leg_relmen()->SetTextFont(43);
        leg_relmen()->SetTextSize(15);
        leg_relmen()->SetFillColor(0);
        leg_relmen.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }

    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_mutr_rh.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFrso = (TGraphErrors*) file_pr->Get("gr_hOFrso");
        TGraphErrors* gr_pr_HCrso = (TGraphErrors*) file_pr->Get("gr_hHCrso");
       
        style(gr_pr_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        TGraphErrors* gr_pr_relrso = (TGraphErrors*) file_pr->Get("gr_hHCrelrso");
        style(gr_pr_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));

        TFile* file_he = TFile::Open(Form("%s/fit_mutr_rh.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFrso = (TGraphErrors*) file_he->Get("gr_hOFrso");
        TGraphErrors* gr_he_HCrso = (TGraphErrors*) file_he->Get("gr_hHCrso");
        
        style(gr_he_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        TGraphErrors* gr_he_relrso = (TGraphErrors*) file_he->Get("gr_hHCrelrso");
        style(gr_he_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));

        TMultiGraph* mg_rso = new TMultiGraph();
        mg_rso->Add(gr_pr_OFrso);
        mg_rso->Add(gr_pr_HCrso);
        mg_rso->Add(gr_he_OFrso);
        mg_rso->Add(gr_he_HCrso);
        
        TMultiGraph* mg_relrso = new TMultiGraph();
        mg_relrso->Add(gr_pr_relrso);
        mg_relrso->Add(gr_he_relrso);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(0));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_rso->Draw("ap");
        mg_rso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_rso->GetHistogram()->SetLineColor(0);
        mg_rso->GetHistogram()->SetMarkerColor(0);
        mg_rso->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        stdfmt(mg_rso);
        mg_rso->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_rso->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_rso->GetHistogram()->SetMinimum(0.0);
        mg_rso->GetHistogram()->SetMaximum(2.5);
        mg_rso->Draw("ap");
        Legend leg_rso("", PadWindow(0.15, 0.40, 0.50, 0.85));
        leg_rso()->AddEntry((TObject*)0, "Tracker L2-L8 & RICH", "");
        leg_rso()->AddEntry(gr_pr_OFrso, "Proton: Official fitting", "lp");
        leg_rso()->AddEntry(gr_pr_HCrso, "Proton: FPM fitting", "lp");
        leg_rso()->AddEntry(gr_he_OFrso, "Helium-4: Official fitting", "lp");
        leg_rso()->AddEntry(gr_he_HCrso, "Helium-4: FPM fitting", "lp");
        leg_rso()->SetTextFont(43);
        leg_rso()->SetTextSize(15);
        leg_rso()->SetFillColor(0);
        leg_rso.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(0));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relrso->Draw("ap");
        mg_relrso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relrso->GetHistogram()->SetLineColor(0);
        mg_relrso->GetHistogram()->SetMarkerColor(0);
        mg_relrso->GetHistogram()->GetXaxis()->SetTitle("#beta");
        mg_relrso->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relrso->GetHistogram()->SetMinimum(0.835);
        mg_relrso->GetHistogram()->SetMaximum(1.105);
        stdfmt(mg_relrso);
        mg_relrso->Draw("ap");
        Legend leg_relrso("", PadWindow(0.15, 0.40, 0.75, 0.95));
        leg_relrso()->AddEntry((TObject*)0, "Tracker L2-L8 & RICH", "");
        leg_relrso()->AddEntry(gr_pr_relrso, "Proton: FPM / Official", "lp");
        leg_relrso()->AddEntry(gr_he_relrso, "Helium-4: FPM / Official", "lp");
        leg_relrso()->SetTextFont(43);
        leg_relrso()->SetTextSize(15);
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    
    std::array<std::string, 4> trpts({ "in", "l1", "l9", "fs" });
    std::array<std::string, 4> trpts_name({ "Tracker L2-L8", "Tracker L1-L8", "Tracker L2-L9", "Tracker L1-L9" });
    std::array<int, 4> trpts_mdr_pr_ck({ 274,  664,  873, 1830 });
    std::array<int, 4> trpts_mdr_pr_hc({ 297,  729,  953, 1928 });
    std::array<int, 4> trpts_mdr_he_ck({ 600, 1239, 1576, 2708 });
    std::array<int, 4> trpts_mdr_he_hc({ 654, 1371, 1696, 2807 });
    for (int ipt = 0; ipt < 4; ++ipt) {
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_geom.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get(Form("gr_hCKmen_%s", trpts[ipt].c_str()));
        TGraphErrors* gr_pr_HCmen = (TGraphErrors*) file_pr->Get(Form("gr_hHCmen_%s", trpts[ipt].c_str()));
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
       
        for (int ip = gr_pr_OFmen->GetN()-1; ip >= 0; --ip) { if (gr_pr_OFmen->GetY()[ip] == 0.0) gr_pr_OFmen->RemovePoint(ip); }
        for (int ip = gr_pr_HCmen->GetN()-1; ip >= 0; --ip) { if (gr_pr_HCmen->GetY()[ip] == 0.0) gr_pr_HCmen->RemovePoint(ip); }

        TFile* file_he = TFile::Open(Form("%s/fit_geom.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get(Form("gr_hCKmen_%s", trpts[ipt].c_str()));
        TGraphErrors* gr_he_HCmen = (TGraphErrors*) file_he->Get(Form("gr_hHCmen_%s", trpts[ipt].c_str()));
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        for (int ip = gr_he_OFmen->GetN()-1; ip >= 0; --ip) { if (gr_he_OFmen->GetY()[ip] == 0.0) gr_he_OFmen->RemovePoint(ip); }
        for (int ip = gr_he_HCmen->GetN()-1; ip >= 0; --ip) { if (gr_he_HCmen->GetY()[ip] == 0.0) gr_he_HCmen->RemovePoint(ip); }

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_HCmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_HCmen);
        
        editor.create();
        editor.cd(1, PadAxis(1));
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->SetMoreLogLabels();
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Bias of Rigidity [%]");
        mg_men->GetHistogram()->SetMinimum(-3.0);
        mg_men->GetHistogram()->SetMaximum( 25.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.55, 0.80, 0.65, 0.85));
        leg_men()->AddEntry((TObject*)0, trpts_name[ipt].c_str(), "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_HCmen, "Proton: FPM fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_HCmen, "Helium-4: FPM fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        editor.save();

        file_pr->Close();
        file_he->Close();
    }

    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_geom.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFrso = (TGraphErrors*) file_pr->Get(Form("gr_hCKrso_%s", trpts[ipt].c_str()));
        TGraphErrors* gr_pr_HCrso = (TGraphErrors*) file_pr->Get(Form("gr_hHCrso_%s", trpts[ipt].c_str()));
       
        style(gr_pr_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        for (int ip = gr_pr_OFrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_OFrso->GetY()[ip] == 0.0) gr_pr_OFrso->RemovePoint(ip); }
        for (int ip = gr_pr_HCrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_HCrso->GetY()[ip] == 0.0) gr_pr_HCrso->RemovePoint(ip); }
        
        TGraphErrors* gr_pr_relrso = (TGraphErrors*) file_pr->Get(Form("gr_hHCrelrso_%s", trpts[ipt].c_str()));
        style(gr_pr_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        
        for (int ip = gr_pr_relrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_relrso->GetY()[ip] == 0.0) gr_pr_relrso->RemovePoint(ip); }

        TFile* file_he = TFile::Open(Form("%s/fit_geom.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFrso = (TGraphErrors*) file_he->Get(Form("gr_hCKrso_%s", trpts[ipt].c_str()));
        TGraphErrors* gr_he_HCrso = (TGraphErrors*) file_he->Get(Form("gr_hHCrso_%s", trpts[ipt].c_str()));
        
        style(gr_he_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_HCrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        for (int ip = gr_he_OFrso->GetN()-1; ip >= 0; --ip) { if (gr_he_OFrso->GetY()[ip] == 0.0) gr_he_OFrso->RemovePoint(ip); }
        for (int ip = gr_he_HCrso->GetN()-1; ip >= 0; --ip) { if (gr_he_HCrso->GetY()[ip] == 0.0) gr_he_HCrso->RemovePoint(ip); }
        
        TGraphErrors* gr_he_relrso = (TGraphErrors*) file_he->Get(Form("gr_hHCrelrso_%s", trpts[ipt].c_str()));
        style(gr_he_relrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        
        for (int ip = gr_he_relrso->GetN()-1; ip >= 0; --ip) { if (gr_he_relrso->GetY()[ip] == 0.0) gr_he_relrso->RemovePoint(ip); }

        TMultiGraph* mg_rso = new TMultiGraph();
        mg_rso->Add(gr_pr_OFrso);
        mg_rso->Add(gr_pr_HCrso);
        mg_rso->Add(gr_he_OFrso);
        mg_rso->Add(gr_he_HCrso);
        
        TMultiGraph* mg_relrso = new TMultiGraph();
        mg_relrso->Add(gr_pr_relrso);
        mg_relrso->Add(gr_he_relrso);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(1));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_rso->Draw("ap");
        mg_rso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_rso->GetHistogram()->SetLineColor(0);
        mg_rso->GetHistogram()->SetMarkerColor(0);
        mg_rso->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Rigidity [%]");
        stdfmt(mg_rso);
        mg_rso->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_rso->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_rso->GetHistogram()->SetMinimum(0.0);
        //mg_rso->GetHistogram()->SetMaximum(2.5);
        mg_rso->Draw("ap");
        Legend leg_rso("", PadWindow(0.15, 0.60, 0.55, 0.85));
        leg_rso()->AddEntry((TObject*)0, trpts_name[ipt].c_str(), "");
        leg_rso()->AddEntry(gr_pr_OFrso, Form("Proton: Official fitting (MDR = %d GV/c)", trpts_mdr_pr_ck[ipt]), "lp");
        leg_rso()->AddEntry(gr_pr_HCrso, Form("Proton: FPM fitting (MDR = %d GV/c)", trpts_mdr_pr_hc[ipt]), "lp");
        leg_rso()->AddEntry(gr_he_OFrso, Form("Helium-4: Official fitting (MDR = %d GV/c)", trpts_mdr_he_ck[ipt]), "lp");
        leg_rso()->AddEntry(gr_he_HCrso, Form("Helium-4: FPM fitting (MDR = %d GV/c)", trpts_mdr_he_hc[ipt]), "lp");
        leg_rso()->SetTextFont(43);
        leg_rso()->SetTextSize(15);
        leg_rso()->SetFillColor(0);
        leg_rso.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(1));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relrso->Draw("ap");
        mg_relrso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relrso->GetHistogram()->SetLineColor(0);
        mg_relrso->GetHistogram()->SetMarkerColor(0);
        mg_relrso->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_relrso->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relrso->GetHistogram()->SetMinimum(0.3);
        mg_relrso->GetHistogram()->SetMaximum(1.01);
        stdfmt(mg_relrso);
        mg_relrso->Draw("ap");
        Legend leg_relrso("", PadWindow(0.60, 0.85, 0.30, 0.55));
        leg_relrso()->AddEntry((TObject*)0, trpts_name[ipt].c_str(), "");
        leg_relrso()->AddEntry(gr_pr_relrso, "Proton: FPM / Official", "lp");
        leg_relrso()->AddEntry(gr_he_relrso, "Helium-4: FPM / Official", "lp");
        leg_relrso()->SetTextFont(43);
        leg_relrso()->SetTextSize(15);
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
    }
    
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_phys_tf.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get("gr_hOFmen");
        TGraphErrors* gr_pr_GMmen = (TGraphErrors*) file_pr->Get("gr_hGMmen");
        TGraphErrors* gr_pr_PHmen = (TGraphErrors*) file_pr->Get("gr_hPHmen");
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_GMmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        style(gr_pr_PHmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kFull)));
       
        for (int ip = gr_pr_OFmen->GetN()-1; ip >= 0; --ip) { if (gr_pr_OFmen->GetY()[ip] == 0.0) gr_pr_OFmen->RemovePoint(ip); }
        for (int ip = gr_pr_GMmen->GetN()-1; ip >= 0; --ip) { if (gr_pr_GMmen->GetY()[ip] == 0.0) gr_pr_GMmen->RemovePoint(ip); }
        for (int ip = gr_pr_PHmen->GetN()-1; ip >= 0; --ip) { if (gr_pr_PHmen->GetY()[ip] == 0.0) gr_pr_PHmen->RemovePoint(ip); }

        TFile* file_he = TFile::Open(Form("%s/fit_phys_tf.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get("gr_hOFmen");
        TGraphErrors* gr_he_GMmen = (TGraphErrors*) file_he->Get("gr_hGMmen");
        TGraphErrors* gr_he_PHmen = (TGraphErrors*) file_he->Get("gr_hPHmen");
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_GMmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        style(gr_he_PHmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kOpen)));
        
        for (int ip = gr_he_OFmen->GetN()-1; ip >= 0; --ip) { if (gr_he_OFmen->GetY()[ip] == 0.0) gr_he_OFmen->RemovePoint(ip); }
        for (int ip = gr_he_GMmen->GetN()-1; ip >= 0; --ip) { if (gr_he_GMmen->GetY()[ip] == 0.0) gr_he_GMmen->RemovePoint(ip); }
        for (int ip = gr_he_PHmen->GetN()-1; ip >= 0; --ip) { if (gr_he_PHmen->GetY()[ip] == 0.0) gr_he_PHmen->RemovePoint(ip); }

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_GMmen);
        mg_men->Add(gr_pr_PHmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_GMmen);
        mg_men->Add(gr_he_PHmen);
        
        editor.create();
        editor.cd(1, PadAxis(1));
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->SetMoreLogLabels();
        mg_men->GetHistogram()->GetXaxis()->SetNoExponent();
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Bias of Rigidity [%]");
        mg_men->GetHistogram()->SetMinimum(-3.0);
        mg_men->GetHistogram()->SetMaximum( 25.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.50, 0.80, 0.55, 0.85));
        leg_men()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_GMmen, "Proton: FPM rigidity fitting", "lp");
        leg_men()->AddEntry(gr_pr_PHmen, "Proton: FPM track fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_GMmen, "Helium-4: FPM rigidity fitting", "lp");
        leg_men()->AddEntry(gr_he_PHmen, "Helium-4: FPM track fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        editor.save();

        file_pr->Close();
        file_he->Close();
    }

    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_phys_tf.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFrso = (TGraphErrors*) file_pr->Get("gr_hOFrso");
        TGraphErrors* gr_pr_GMrso = (TGraphErrors*) file_pr->Get("gr_hGMrso");
        TGraphErrors* gr_pr_PHrso = (TGraphErrors*) file_pr->Get("gr_hPHrso");
       
        style(gr_pr_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_GMrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        style(gr_pr_PHrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kFull)));
        
        for (int ip = gr_pr_OFrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_OFrso->GetY()[ip] == 0.0) gr_pr_OFrso->RemovePoint(ip); }
        for (int ip = gr_pr_GMrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_GMrso->GetY()[ip] == 0.0) gr_pr_GMrso->RemovePoint(ip); }
        for (int ip = gr_pr_PHrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_PHrso->GetY()[ip] == 0.0) gr_pr_PHrso->RemovePoint(ip); }
        
        TGraphErrors* gr_pr_GMrelrso = (TGraphErrors*) file_pr->Get("gr_hGMrelrso");
        TGraphErrors* gr_pr_PHrelrso = (TGraphErrors*) file_pr->Get("gr_hPHrelrso");

        style(gr_pr_GMrelrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        style(gr_pr_PHrelrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kFull)));
        
        for (int ip = gr_pr_GMrelrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_GMrelrso->GetY()[ip] == 0.0) gr_pr_GMrelrso->RemovePoint(ip); }
        for (int ip = gr_pr_PHrelrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_PHrelrso->GetY()[ip] == 0.0) gr_pr_PHrelrso->RemovePoint(ip); }

        TFile* file_he = TFile::Open(Form("%s/fit_phys_tf.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFrso = (TGraphErrors*) file_he->Get("gr_hOFrso");
        TGraphErrors* gr_he_GMrso = (TGraphErrors*) file_he->Get("gr_hGMrso");
        TGraphErrors* gr_he_PHrso = (TGraphErrors*) file_he->Get("gr_hPHrso");
        
        style(gr_he_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_GMrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        style(gr_he_PHrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kOpen)));
        
        for (int ip = gr_he_OFrso->GetN()-1; ip >= 0; --ip) { if (gr_he_OFrso->GetY()[ip] == 0.0) gr_he_OFrso->RemovePoint(ip); }
        for (int ip = gr_he_GMrso->GetN()-1; ip >= 0; --ip) { if (gr_he_GMrso->GetY()[ip] == 0.0) gr_he_GMrso->RemovePoint(ip); }
        for (int ip = gr_he_PHrso->GetN()-1; ip >= 0; --ip) { if (gr_he_PHrso->GetY()[ip] == 0.0) gr_he_PHrso->RemovePoint(ip); }
        
        TGraphErrors* gr_he_GMrelrso = (TGraphErrors*) file_he->Get("gr_hGMrelrso");
        TGraphErrors* gr_he_PHrelrso = (TGraphErrors*) file_he->Get("gr_hPHrelrso");
        
        style(gr_he_GMrelrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        style(gr_he_PHrelrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kOpen)));
        
        for (int ip = gr_he_GMrelrso->GetN()-1; ip >= 0; --ip) { if (gr_he_GMrelrso->GetY()[ip] == 0.0) gr_he_GMrelrso->RemovePoint(ip); }
        for (int ip = gr_he_PHrelrso->GetN()-1; ip >= 0; --ip) { if (gr_he_PHrelrso->GetY()[ip] == 0.0) gr_he_PHrelrso->RemovePoint(ip); }

        TMultiGraph* mg_rso = new TMultiGraph();
        mg_rso->Add(gr_pr_OFrso);
        mg_rso->Add(gr_pr_GMrso);
        mg_rso->Add(gr_pr_PHrso);
        mg_rso->Add(gr_he_OFrso);
        mg_rso->Add(gr_he_GMrso);
        mg_rso->Add(gr_he_PHrso);
        
        TMultiGraph* mg_relrso = new TMultiGraph();
        mg_relrso->Add(gr_pr_GMrelrso);
        mg_relrso->Add(gr_pr_PHrelrso);
        mg_relrso->Add(gr_he_GMrelrso);
        mg_relrso->Add(gr_he_PHrelrso);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(1));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_rso->Draw("ap");
        mg_rso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_rso->GetHistogram()->SetLineColor(0);
        mg_rso->GetHistogram()->SetMarkerColor(0);
        mg_rso->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Rigidity [%]");
        stdfmt(mg_rso);
        mg_rso->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_rso->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_rso->GetHistogram()->SetMinimum(0.0);
        mg_rso->GetHistogram()->SetMaximum(30.0);
        mg_rso->Draw("ap");
        Legend leg_rso("", PadWindow(0.50, 0.80, 0.45, 0.85));
        leg_rso()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_rso()->AddEntry(gr_pr_OFrso, "Proton: Official fitting", "lp");
        leg_rso()->AddEntry(gr_pr_GMrso, "Proton: FPM rigidity fitting", "lp");
        leg_rso()->AddEntry(gr_pr_PHrso, "Proton: FPM track fitting", "lp");
        leg_rso()->AddEntry(gr_he_OFrso, "Helium-4: Official fitting", "lp");
        leg_rso()->AddEntry(gr_he_GMrso, "Helium-4: FPM rigidity fitting", "lp");
        leg_rso()->AddEntry(gr_he_PHrso, "Helium-4: FPM track fitting", "lp");
        leg_rso()->SetTextFont(43);
        leg_rso()->SetTextSize(15);
        leg_rso()->SetFillColor(0);
        leg_rso.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(1));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relrso->Draw("ap");
        mg_relrso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relrso->GetHistogram()->SetLineColor(0);
        mg_relrso->GetHistogram()->SetMarkerColor(0);
        mg_relrso->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_relrso->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relrso->GetHistogram()->SetMinimum(0.0);
        mg_relrso->GetHistogram()->SetMaximum(1.01);
        stdfmt(mg_relrso);
        mg_relrso->Draw("ap");
        Legend leg_relrso("", PadWindow(0.50, 0.80, 0.30, 0.70));
        leg_relrso()->AddEntry((TObject*)0, "Tracker L2-L8 & TOF", "");
        leg_relrso()->AddEntry(gr_pr_GMrelrso, "Proton: FPM rigidity / Official", "lp");
        leg_relrso()->AddEntry(gr_pr_PHrelrso, "Proton: FPM track / Official", "lp");
        leg_relrso()->AddEntry(gr_he_GMrelrso, "Helium-4: FPM rigidity / Official", "lp");
        leg_relrso()->AddEntry(gr_he_PHrelrso, "Helium-4: FPM track / Official", "lp");
        leg_relrso()->SetTextFont(43);
        leg_relrso()->SetTextSize(15);
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
  
    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_phys_rh.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFmen = (TGraphErrors*) file_pr->Get("gr_hOFmen");
        TGraphErrors* gr_pr_GMmen = (TGraphErrors*) file_pr->Get("gr_hGMmen");
        TGraphErrors* gr_pr_PHmen = (TGraphErrors*) file_pr->Get("gr_hPHmen");
       
        style(gr_pr_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_GMmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        style(gr_pr_PHmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kFull)));
       
        for (int ip = gr_pr_OFmen->GetN()-1; ip >= 0; --ip) { if (gr_pr_OFmen->GetY()[ip] == 0.0) gr_pr_OFmen->RemovePoint(ip); }
        for (int ip = gr_pr_GMmen->GetN()-1; ip >= 0; --ip) { if (gr_pr_GMmen->GetY()[ip] == 0.0) gr_pr_GMmen->RemovePoint(ip); }
        for (int ip = gr_pr_PHmen->GetN()-1; ip >= 0; --ip) { if (gr_pr_PHmen->GetY()[ip] == 0.0) gr_pr_PHmen->RemovePoint(ip); }

        TFile* file_he = TFile::Open(Form("%s/fit_phys_rh.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFmen = (TGraphErrors*) file_he->Get("gr_hOFmen");
        TGraphErrors* gr_he_GMmen = (TGraphErrors*) file_he->Get("gr_hGMmen");
        TGraphErrors* gr_he_PHmen = (TGraphErrors*) file_he->Get("gr_hPHmen");
        
        style(gr_he_OFmen, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_GMmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        style(gr_he_PHmen, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kOpen)));
        
        for (int ip = gr_he_OFmen->GetN()-1; ip >= 0; --ip) { if (gr_he_OFmen->GetY()[ip] == 0.0) gr_he_OFmen->RemovePoint(ip); }
        for (int ip = gr_he_GMmen->GetN()-1; ip >= 0; --ip) { if (gr_he_GMmen->GetY()[ip] == 0.0) gr_he_GMmen->RemovePoint(ip); }
        for (int ip = gr_he_PHmen->GetN()-1; ip >= 0; --ip) { if (gr_he_PHmen->GetY()[ip] == 0.0) gr_he_PHmen->RemovePoint(ip); }

        TMultiGraph* mg_men = new TMultiGraph();
        mg_men->Add(gr_pr_OFmen);
        mg_men->Add(gr_pr_GMmen);
        mg_men->Add(gr_pr_PHmen);
        mg_men->Add(gr_he_OFmen);
        mg_men->Add(gr_he_GMmen);
        mg_men->Add(gr_he_PHmen);
        
        editor.create();
        editor.cd(1, PadAxis(1));
        mg_men->Draw("ap");
        mg_men->GetHistogram()->GetXaxis()->SetRangeUser(3.5, 80.0);
        mg_men->GetHistogram()->GetXaxis()->SetMoreLogLabels();
        mg_men->GetHistogram()->GetXaxis()->CenterTitle();
        mg_men->GetHistogram()->SetLineColor(0);
        mg_men->GetHistogram()->SetMarkerColor(0);
        mg_men->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_men->GetHistogram()->GetYaxis()->SetTitle("Bias of Rigidity [%]");
        mg_men->GetHistogram()->SetMinimum(-2.0);
        mg_men->GetHistogram()->SetMaximum( 5.0);
        mg_men->Draw("ap");
        Legend leg_men("", PadWindow(0.50, 0.80, 0.55, 0.85));
        leg_men()->AddEntry((TObject*)0, "Tracker L2-L8 & RICH", "");
        leg_men()->AddEntry(gr_pr_OFmen, "Proton: Official fitting", "lp");
        leg_men()->AddEntry(gr_pr_GMmen, "Proton: FPM rigidity fitting", "lp");
        leg_men()->AddEntry(gr_pr_PHmen, "Proton: FPM track fitting", "lp");
        leg_men()->AddEntry(gr_he_OFmen, "Helium-4: Official fitting", "lp");
        leg_men()->AddEntry(gr_he_GMmen, "Helium-4: FPM rigidity fitting", "lp");
        leg_men()->AddEntry(gr_he_PHmen, "Helium-4: FPM track fitting", "lp");
        leg_men()->SetTextFont(43);
        leg_men()->SetTextSize(15);
        leg_men()->SetFillColor(0);
        leg_men.draw();
        editor.save();

        file_pr->Close();
        file_he->Close();
    }

    if (true) {
        TFile* file_pr = TFile::Open(Form("%s/fit_phys_rh.root", dir_pr.c_str()));
        TGraphErrors* gr_pr_OFrso = (TGraphErrors*) file_pr->Get("gr_hOFrso");
        TGraphErrors* gr_pr_GMrso = (TGraphErrors*) file_pr->Get("gr_hGMrso");
        TGraphErrors* gr_pr_PHrso = (TGraphErrors*) file_pr->Get("gr_hPHrso");
       
        style(gr_pr_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kFull)));
        style(gr_pr_GMrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        style(gr_pr_PHrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kFull)));
        
        for (int ip = gr_pr_OFrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_OFrso->GetY()[ip] == 0.0) gr_pr_OFrso->RemovePoint(ip); }
        for (int ip = gr_pr_GMrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_GMrso->GetY()[ip] == 0.0) gr_pr_GMrso->RemovePoint(ip); }
        for (int ip = gr_pr_PHrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_PHrso->GetY()[ip] == 0.0) gr_pr_PHrso->RemovePoint(ip); }
        
        TGraphErrors* gr_pr_GMrelrso = (TGraphErrors*) file_pr->Get("gr_hGMrelrso");
        TGraphErrors* gr_pr_PHrelrso = (TGraphErrors*) file_pr->Get("gr_hPHrelrso");

        style(gr_pr_GMrelrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kFull)));
        style(gr_pr_PHrelrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kFull)));
        
        for (int ip = gr_pr_GMrelrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_GMrelrso->GetY()[ip] == 0.0) gr_pr_GMrelrso->RemovePoint(ip); }
        for (int ip = gr_pr_PHrelrso->GetN()-1; ip >= 0; --ip) { if (gr_pr_PHrelrso->GetY()[ip] == 0.0) gr_pr_PHrelrso->RemovePoint(ip); }

        TFile* file_he = TFile::Open(Form("%s/fit_phys_rh.root", dir_he.c_str()));
        TGraphErrors* gr_he_OFrso = (TGraphErrors*) file_he->Get("gr_hOFrso");
        TGraphErrors* gr_he_GMrso = (TGraphErrors*) file_he->Get("gr_hGMrso");
        TGraphErrors* gr_he_PHrso = (TGraphErrors*) file_he->Get("gr_hPHrso");
        
        style(gr_he_OFrso, Line(kBlue), Marker(kBlue, MarkerStyle(MarkerShape::kCircle,  MarkerType::kOpen)));
        style(gr_he_GMrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        style(gr_he_PHrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kOpen)));
        
        for (int ip = gr_he_OFrso->GetN()-1; ip >= 0; --ip) { if (gr_he_OFrso->GetY()[ip] == 0.0) gr_he_OFrso->RemovePoint(ip); }
        for (int ip = gr_he_GMrso->GetN()-1; ip >= 0; --ip) { if (gr_he_GMrso->GetY()[ip] == 0.0) gr_he_GMrso->RemovePoint(ip); }
        for (int ip = gr_he_PHrso->GetN()-1; ip >= 0; --ip) { if (gr_he_PHrso->GetY()[ip] == 0.0) gr_he_PHrso->RemovePoint(ip); }
        
        TGraphErrors* gr_he_GMrelrso = (TGraphErrors*) file_he->Get("gr_hGMrelrso");
        TGraphErrors* gr_he_PHrelrso = (TGraphErrors*) file_he->Get("gr_hPHrelrso");
        
        style(gr_he_GMrelrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kDiamond, MarkerType::kOpen)));
        style(gr_he_PHrelrso, Line(kRed),  Marker(kRed,  MarkerStyle(MarkerShape::kCross,   MarkerType::kOpen)));
        
        for (int ip = gr_he_GMrelrso->GetN()-1; ip >= 0; --ip) { if (gr_he_GMrelrso->GetY()[ip] == 0.0) gr_he_GMrelrso->RemovePoint(ip); }
        for (int ip = gr_he_PHrelrso->GetN()-1; ip >= 0; --ip) { if (gr_he_PHrelrso->GetY()[ip] == 0.0) gr_he_PHrelrso->RemovePoint(ip); }

        TMultiGraph* mg_rso = new TMultiGraph();
        mg_rso->Add(gr_pr_OFrso);
        mg_rso->Add(gr_pr_GMrso);
        mg_rso->Add(gr_pr_PHrso);
        mg_rso->Add(gr_he_OFrso);
        mg_rso->Add(gr_he_GMrso);
        mg_rso->Add(gr_he_PHrso);
        
        TMultiGraph* mg_relrso = new TMultiGraph();
        mg_relrso->Add(gr_pr_GMrelrso);
        mg_relrso->Add(gr_pr_PHrelrso);
        mg_relrso->Add(gr_he_GMrelrso);
        mg_relrso->Add(gr_he_PHrelrso);

        editor.create("", 1, 2, PadMargin(0.1, 0.25));

        TVirtualPad* pad1 = editor.cd(1, PadAxis(1));
        pad1->SetPad(0.0, 0.4, 1.0, 1.0);
        pad1->SetBottomMargin(0.015);
        //pad1->SetGridx();
        mg_rso->Draw("ap");
        mg_rso->GetHistogram()->GetXaxis()->SetRangeUser(3.5, 80.0);
        mg_rso->GetHistogram()->GetXaxis()->SetMoreLogLabels();
        mg_rso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_rso->GetHistogram()->SetLineColor(0);
        mg_rso->GetHistogram()->SetMarkerColor(0);
        mg_rso->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_rso->GetHistogram()->GetYaxis()->SetTitle("Resolution of Rigidity [%]");
        stdfmt(mg_rso);
        mg_rso->GetHistogram()->GetXaxis()->SetTitleSize(0.0);
        mg_rso->GetHistogram()->GetXaxis()->SetLabelSize(0.0);
        mg_rso->GetHistogram()->SetMinimum(0.0);
        mg_rso->GetHistogram()->SetMaximum(30.0);
        mg_rso->Draw("ap");
        Legend leg_rso("", PadWindow(0.15, 0.45, 0.45, 0.85));
        leg_rso()->AddEntry((TObject*)0, "Tracker L2-L8 & RICH", "");
        leg_rso()->AddEntry(gr_pr_OFrso, "Proton: Official fitting", "lp");
        leg_rso()->AddEntry(gr_pr_GMrso, "Proton: FPM rigidity fitting", "lp");
        leg_rso()->AddEntry(gr_pr_PHrso, "Proton: FPM track fitting", "lp");
        leg_rso()->AddEntry(gr_he_OFrso, "Helium-4: Official fitting", "lp");
        leg_rso()->AddEntry(gr_he_GMrso, "Helium-4: FPM rigidity fitting", "lp");
        leg_rso()->AddEntry(gr_he_PHrso, "Helium-4: FPM track fitting", "lp");
        leg_rso()->SetTextFont(43);
        leg_rso()->SetTextSize(15);
        leg_rso()->SetFillColor(0);
        leg_rso.draw();
        
        TVirtualPad* pad2 = editor.cd(2, PadAxis(1));
        pad2->SetPad(0.0, 0.0, 1.0, 0.4);
        pad2->SetTopMargin(0.015);
        //pad2->SetGridx();
        mg_relrso->Draw("ap");
        mg_relrso->GetHistogram()->GetXaxis()->SetRangeUser(3.5, 80.0);
        mg_relrso->GetHistogram()->GetXaxis()->SetMoreLogLabels();
        mg_relrso->GetHistogram()->GetXaxis()->CenterTitle();
        mg_relrso->GetHistogram()->SetLineColor(0);
        mg_relrso->GetHistogram()->SetMarkerColor(0);
        mg_relrso->GetHistogram()->GetXaxis()->SetTitle("Rigidity [GV/c]");
        mg_relrso->GetHistogram()->GetYaxis()->SetTitle("Ratio");
        mg_relrso->GetHistogram()->SetMinimum(0.0);
        mg_relrso->GetHistogram()->SetMaximum(1.01);
        stdfmt(mg_relrso);
        mg_relrso->Draw("ap");
        Legend leg_relrso("", PadWindow(0.55, 0.85, 0.30, 0.70));
        leg_relrso()->AddEntry((TObject*)0, "Tracker L2-L8 & RICH", "");
        leg_relrso()->AddEntry(gr_pr_GMrelrso, "Proton: FPM rigidity / Official", "lp");
        leg_relrso()->AddEntry(gr_pr_PHrelrso, "Proton: FPM track / Official", "lp");
        leg_relrso()->AddEntry(gr_he_GMrelrso, "Helium-4: FPM rigidity / Official", "lp");
        leg_relrso()->AddEntry(gr_he_PHrelrso, "Helium-4: FPM track / Official", "lp");
        leg_relrso()->SetTextFont(43);
        leg_relrso()->SetTextSize(15);
        leg_relrso()->SetFillColor(0);
        leg_relrso.draw();
        
        editor.save();

        file_pr->Close();
        file_he->Close();
    }
  
    editor.close();
    
    return 1;
}
