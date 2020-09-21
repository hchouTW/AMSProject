#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "DataFit1D.h"
#include "DataFit1D.C"

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

    TH1D* hresy = (TH1D*) TFile::Open("out_mes/mes_resy")->Get("hresy");
    TH1D* htfql = (TH1D*) TFile::Open("out_mes/mes_tfql")->Get("htfql");
    TH1D* htfqm = (TH1D*) TFile::Open("out_mes/mes_tfqm")->Get("htfqm");
    TH1D* htfqh = (TH1D*) TFile::Open("out_mes/mes_tfqh")->Get("htfqh");


    PdfEditor editor(Window(), "mes", "out_Dec04");

    editor.close();
    
    TFile * ofle = new TFile("out_Dec04/mes.root", "RECREATE");
    ofle->cd();

    ofle->Write();
    ofle->Close();

    return 1;
}
