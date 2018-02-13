//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__

#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

using namespace std;
using namespace MGROOT;
using namespace TrackSys;

int main(int argc, char * argv[]) {
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
   
    Hist::Load("cc_fill.root", "dat");
    Hist* hMCxtk[9];
    for (Int_t il = 0; il < 9; ++il)
        hMCxtk[il] = Hist::Head(STR_FMT("hMCxtk%d", il+1));
    
    Hist* hMCxtkh[9];
    for (Int_t il = 0; il < 9; ++il)
        hMCxtkh[il] = Hist::Head(STR_FMT("hMCxtkh%d", il+1));
    
    Hist* hMCytk[9];
    for (Int_t il = 0; il < 9; ++il)
        hMCytk[il] = Hist::Head(STR_FMT("hMCytk%d", il+1));
    
    Hist* hMCytkh[9];
    for (Int_t il = 0; il < 9; ++il)
        hMCytkh[il] = Hist::Head(STR_FMT("hMCytkh%d", il+1));
    
    const Axis& AXmom = hMCxtkh[0]->xaxis();
   
    PdfEditor editor(Window(1800, 900*12));

    TFile * ofle = new TFile("cc_fit.root", "RECREATE");
    ofle->cd();
    
    Hist::AddDirectory();

    Hist* hMCxtkhS[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il)
        hMCxtkhS[il] = Hist::New(STR_FMT("hMCxtkh%dS", il+1), HistAxis(AXmom));
    
    Hist* hMCxtkhE99[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il)
        hMCxtkhE99[il] = Hist::New(STR_FMT("hMCxtkh%dE99", il+1), HistAxis(AXmom));
    
    Hist* hMCxtkhE999[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il)
        hMCxtkhE999[il] = Hist::New(STR_FMT("hMCxtkh%dE999", il+1), HistAxis(AXmom));
    
    Hist* hMCxtkhE9999[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il)
        hMCxtkhE9999[il] = Hist::New(STR_FMT("hMCxtkh%dE9999", il+1), HistAxis(AXmom));
    
    Hist* hMCytkhS[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il)
        hMCytkhS[il] = Hist::New(STR_FMT("hMCytkh%dS", il+1), HistAxis(AXmom));
    
    Hist* hMCytkhE99[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il)
        hMCytkhE99[il] = Hist::New(STR_FMT("hMCytkh%dE99", il+1), HistAxis(AXmom));
    
    Hist* hMCytkhE999[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il)
        hMCytkhE999[il] = Hist::New(STR_FMT("hMCytkh%dE999", il+1), HistAxis(AXmom));
    
    Hist* hMCytkhE9999[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il)
        hMCytkhE9999[il] = Hist::New(STR_FMT("hMCytkh%dE9999", il+1), HistAxis(AXmom));
    
    
    const Double_t stable = 1.5;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    std::vector<Hist*> vhMCxtkh[9];
    for (Int_t il = 0; il < 9; ++il)
        vhMCxtkh[il] = Hist::ProjectAll(HistProj::kY, hMCxtkh[il]);
    
    std::vector<Hist*> vhMCytkh[9];
    for (Int_t il = 0; il < 9; ++il)
        vhMCytkh[il] = Hist::ProjectAll(HistProj::kY, hMCytkh[il]);

    for (int it = 1; it <= AXmom.nbin(); ++it) {
        double mom = AXmom.center(it, AxisScale::kLog);
        Double_t norm = (30./mom+1.0);
      
        for (Int_t il = 0; il < 9; ++il) {
            Double_t MCmax = (*vhMCxtkh[il].at(it))()->GetBinCenter((*vhMCxtkh[il].at(it))()->GetMaximumBin());
            Double_t MCrms = 0.5 * (*vhMCxtkh[il].at(it))()->GetRMS();
            gaus->SetParameters(1000, MCmax, MCrms);
            (*vhMCxtkh[il].at(it))()->Fit(gaus, "q0", "");
            (*vhMCxtkh[il].at(it))()->Fit(gaus, "q0", "", MCmax-stable*MCrms, MCmax+stable*MCrms);
            (*vhMCxtkh[il].at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhMCxtkh[il].at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhMCxtkh[il].at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhMCxtkh[il].at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
            (*hMCxtkhS[il])()->SetBinContent(it, gaus->GetParameter(2) * norm);
            (*hMCxtkhS[il])()->SetBinError  (it, gaus->GetParError(2) * norm);
         
            Double_t sel = 0;
            Double_t ttl = (*hMCxtk[il])()->GetBinContent(it);
            TGraph grMCeft;
            for (Int_t ib = 0; ib <= (*vhMCxtkh[il].at(it))()->GetNbinsX(); ++ib) {
                Double_t val = (*vhMCxtkh[il].at(it))()->GetXaxis()->GetBinUpEdge(ib);
                sel += (*vhMCxtkh[il].at(it))()->GetBinContent(ib);
                Double_t eft = (sel / ttl);
                grMCeft.SetPoint(ib, eft, val * norm);
            }
            (*hMCxtkhE99[il])()->SetBinContent(it, mom, grMCeft.Eval(0.99));
            (*hMCxtkhE999[il])()->SetBinContent(it, mom, grMCeft.Eval(0.999));
            (*hMCxtkhE9999[il])()->SetBinContent(it, mom, grMCeft.Eval(0.9999));
        }
        
        for (Int_t il = 0; il < 9; ++il) {
            Double_t MCmax = (*vhMCytkh[il].at(it))()->GetBinCenter((*vhMCytkh[il].at(it))()->GetMaximumBin());
            Double_t MCrms = 0.5 * (*vhMCytkh[il].at(it))()->GetRMS();
            gaus->SetParameters(1000, MCmax, MCrms);
            (*vhMCytkh[il].at(it))()->Fit(gaus, "q0", "");
            (*vhMCytkh[il].at(it))()->Fit(gaus, "q0", "", MCmax-stable*MCrms, MCmax+stable*MCrms);
            (*vhMCytkh[il].at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhMCytkh[il].at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhMCytkh[il].at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
            (*vhMCytkh[il].at(it))()->Fit(gaus, "q0", "", gaus->GetParameter(1)-stable*gaus->GetParameter(2), gaus->GetParameter(1)+stable*gaus->GetParameter(2));
    
            (*hMCytkhS[il])()->SetBinContent(it, gaus->GetParameter(2) * norm);
            (*hMCytkhS[il])()->SetBinError  (it, gaus->GetParError(2) * norm);
         
            Double_t sel = 0;
            Double_t ttl = (*hMCytk[il])()->GetBinContent(it);
            TGraph grMCeft;
            for (Int_t ib = 0; ib <= (*vhMCytkh[il].at(it))()->GetNbinsX(); ++ib) {
                Double_t val = (*vhMCytkh[il].at(it))()->GetXaxis()->GetBinUpEdge(ib);
                sel += (*vhMCytkh[il].at(it))()->GetBinContent(ib);
                Double_t eft = (sel / ttl);
                grMCeft.SetPoint(ib, eft, val * norm);
            }
            (*hMCytkhE99[il])()->SetBinContent(it, mom, grMCeft.Eval(0.99));
            (*hMCytkhE999[il])()->SetBinContent(it, mom, grMCeft.Eval(0.999));
            (*hMCytkhE9999[il])()->SetBinContent(it, mom, grMCeft.Eval(0.9999));
        }

        editor.create("", 1, 12);
        for (Int_t il = 2; il < 8; ++il) {
            Double_t sigma = std::max((*hMCxtkhS[il])()->GetBinContent(it), (*hMCytkhS[il])()->GetBinContent(it)) / norm;
            
            vhMCxtkh[il].at(it)->style(Fill(), Line(kBlue), Marker(kBlue));
            vhMCytkh[il].at(it)->style(Fill(), Line(kRed), Marker(kRed));
            TextList list({ 
                Text(STR_FMT("X With NormChisqY<8 (Events %10.0f)", (*hMCxtk[il])()->GetBinContent(it)), kBlue),
                Text(STR_FMT("Central Sigma %8.5f [cm] (%8.5f)", (*hMCxtkhS[il])()->GetBinContent(it), (*hMCxtkhS[il])()->GetBinContent(it)/norm), kBlue),
                Text(STR_FMT("[1.00%] %8.5f [cm] (%8.5f)", (*hMCxtkhE99[il])()->GetBinContent(it), (*hMCxtkhE99[il])()->GetBinContent(it)/norm), kBlue),
                Text(STR_FMT("[0.10%] %8.5f [cm] (%8.5f)", (*hMCxtkhE999[il])()->GetBinContent(it), (*hMCxtkhE999[il])()->GetBinContent(it)/norm), kBlue),
                Text(STR_FMT("[0.01%] %8.5f [cm] (%8.5f)", (*hMCxtkhE9999[il])()->GetBinContent(it), (*hMCxtkhE9999[il])()->GetBinContent(it)/norm), kBlue),
                Text(STR_FMT("Y With NormChisqX<8 (Events %10.0f)", (*hMCytk[il])()->GetBinContent(it)), kRed),
                Text(STR_FMT("Central Sigma %8.5f [cm] (%8.5f)", (*hMCytkhS[il])()->GetBinContent(it), (*hMCytkhS[il])()->GetBinContent(it)/norm), kRed),
                Text(STR_FMT("[1.00%] %8.5f [cm] (%8.5f)", (*hMCytkhE99[il])()->GetBinContent(it), (*hMCytkhE99[il])()->GetBinContent(it)/norm), kRed),
                Text(STR_FMT("[0.10%] %8.5f [cm] (%8.5f)", (*hMCytkhE999[il])()->GetBinContent(it), (*hMCytkhE999[il])()->GetBinContent(it)/norm), kRed),
                Text(STR_FMT("[0.01%] %8.5f [cm] (%8.5f)", (*hMCytkhE9999[il])()->GetBinContent(it), (*hMCytkhE9999[il])()->GetBinContent(it)/norm), kRed),
            });
           
            editor.cd(il*2-3, PadAxis(0, 1));
            THStack* plotMCxtkh = Hist::Collect(STR_FMT("plot%03dMCxtkh%d", it, il), HistList({ vhMCxtkh[il].at(it), vhMCytkh[il].at(it) }));
            plotMCxtkh->Draw("nostack hist");
            plotMCxtkh->GetXaxis()->SetTitle("Residual [cm * (30/Momentum+1)^-1]");
            plotMCxtkh->GetYaxis()->SetTitle("Events/Bin");
            plotMCxtkh->GetYaxis()->SetMoreLogLabels();
            plotMCxtkh->Draw("nostack hist");
            TitleDraw(STR_FMT("Momentum %6.2f~%6.2f GeV   LayerJ %d", AXmom(it-1), AXmom(it), il+1));
            TextDraw(list);

            editor.cd(il*2-2, PadAxis(0, 1));
            THStack* plotMCxtkhScan = Hist::Collect(STR_FMT("plot%03dMCxtkhScan%d", it, il), HistList({ vhMCxtkh[il].at(it), vhMCytkh[il].at(it) }));
            plotMCxtkhScan->Draw("nostack hist");
            plotMCxtkhScan->GetXaxis()->SetRangeUser(-8. * sigma, 8. * sigma);
            plotMCxtkhScan->GetXaxis()->SetTitle("Residual [cm * (30/Momentum+1)^-1]");
            plotMCxtkhScan->GetYaxis()->SetTitle("Events/Bin");
            plotMCxtkhScan->GetYaxis()->SetMoreLogLabels();
            plotMCxtkhScan->Draw("nostack hist");
            TitleDraw(STR_FMT("Momentum %6.2f~%6.2f GeV   LayerJ %d", AXmom(it-1), AXmom(it), il+1));
            TextDraw(list);
        }
        editor.save();
    } 
   
    THStack* chMCxtkh[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il) {
        hMCxtkhS[il]->style(Fill(), Line(kBlue), Marker(kBlue));
        hMCxtkhE99[il]->style(Fill(), Line(kRed), Marker(kRed));
        hMCxtkhE999[il]->style(Fill(), Line(kRed+1), Marker(kRed+1));
        hMCxtkhE9999[il]->style(Fill(), Line(kRed+2), Marker(kRed+2));
        chMCxtkh[il] = Hist::Collect(STR_FMT("chMCxtkh%d", il), HistList({ hMCxtkhS[il], hMCxtkhE99[il], hMCxtkhE999[il], hMCxtkhE9999[il] }));
        chMCxtkh[il]->Write();
    }
    
    THStack* chMCytkh[9] = { nullptr };
    for (Int_t il = 0; il < 9; ++il) {
        hMCytkhS[il]->style(Fill(), Line(kBlue), Marker(kBlue));
        hMCytkhE99[il]->style(Fill(), Line(kRed), Marker(kRed));
        hMCytkhE999[il]->style(Fill(), Line(kRed+1), Marker(kRed+1));
        hMCytkhE9999[il]->style(Fill(), Line(kRed+2), Marker(kRed+2));
        chMCytkh[il] = Hist::Collect(STR_FMT("chMCytkh%d", il), HistList({ hMCytkhS[il], hMCytkhE99[il], hMCytkhE999[il], hMCytkhE9999[il] }));
        chMCytkh[il]->Write();
    }
        
    editor.create("", 1, 12);
    for (Int_t il = 2; il < 8; ++il) {
        editor.cd(il*2-3, PadAxis(1, 1));
        chMCxtkh[il]->Draw("nostack hist");
        chMCxtkh[il]->GetXaxis()->SetTitle("Momentum [GeV]");
        chMCxtkh[il]->GetYaxis()->SetTitle("Residual [cm]");
        chMCxtkh[il]->GetYaxis()->SetMoreLogLabels();
        chMCxtkh[il]->Draw("nostack hist");
        TitleDraw(STR_FMT("LayerJ %d (X With NormChisqY<8)", il+1));
        TextDraw(TextList({ 
            Text("Sigma", kBlue),
            Text("[1.00%]", kRed),
            Text("[0.10%]", kRed+1),
            Text("[0.01%]", kRed+2)
        }), TextAlign(0.85, 0.80, 32));
        
        editor.cd(il*2-2, PadAxis(1, 1));
        chMCytkh[il]->Draw("nostack hist");
        chMCytkh[il]->GetXaxis()->SetTitle("Momentum [GeV]");
        chMCytkh[il]->GetYaxis()->SetTitle("Residual [cm]");
        chMCytkh[il]->GetYaxis()->SetMoreLogLabels();
        chMCytkh[il]->Draw("nostack hist");
        TitleDraw(STR_FMT("LayerJ %d (Y With NormChisqX<8)", il+1));
        TextDraw(TextList({ 
            Text("Sigma", kBlue),
            Text("[1.00%]", kRed),
            Text("[0.10%]", kRed+1),
            Text("[0.01%]", kRed+2)
        }), TextAlign(0.85, 0.80, 32));
    }
    editor.save();
    
    ofle->Write();
    ofle->Close();

    editor.close();

    return 0;
}
