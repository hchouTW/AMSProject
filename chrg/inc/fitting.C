#include <CPPLibs.h>
#include <ROOTLibs.h>
#include <TrSys.h>

double LinkSMG3func(double *x, double *p) {
    return TrSys::MultiGausFunc(x[0], 0.0, std::vector<std::array<double, 2>>{ {p[0], p[1]}, {p[2], p[3]}, {p[4], p[5]} });
}
double LinkSMG4func(double *x, double *p) {
    return TrSys::MultiGausFunc(x[0], 0.0, std::vector<std::array<double, 2>>{ {p[0], p[1]}, {p[2], p[3]}, {p[4], p[5]}, {p[6], p[7]} });
}
double LinkALGfunc(double *x, double *p) {
    return (p[4] * TrSys::LandauGausFunc(x[0], p[0], p[1], p[2], p[3]));
}

int main() {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    
    Hist::Load("YiMdst.root", "/eos/ams/user/h/hchou/AMSData/subj/chrg/20Jul08/iss16");
    Hist* hZtkEXtf = Hist::Head("hZtkEXtf");
    Hist* hZtfEXtk = Hist::Head("hZtfEXtk");
   
    Axis AXz  = hZtkEXtf->xaxis();
    Axis AXex = hZtkEXtf->yaxis();

    hZtkEXtf->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    hZtfEXtk->style(Line(kRed, 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle)));
    
    Hist* hZKpa = Hist::New("hZKpa", HistAxis(AXz, "#kappa"));
    Hist* hZMpv = Hist::New("hZMpv", HistAxis(AXz, "MPV"));
    Hist* hZSgm = Hist::New("hZSgm", HistAxis(AXz, "#sigma"));
    Hist* hZFlc = Hist::New("hZFlc", HistAxis(AXz, "Fluc"));

    // Sgm
    // 3.07653e-02 + 7.98278e-02 * exp(-6.97023e-02 * Z * Z)

    TF1* ALGfunc = new TF1("ALGfunc", LinkALGfunc, 0.0, 10000.0, 5);
    ALGfunc->SetLineColor(kBlue);
    ALGfunc->SetNpx(100000);

    ALGfunc->SetParameter(0, 0.05);
    ALGfunc->SetParameter(1, 1.0);
    ALGfunc->SetParameter(2, 0.10);
    ALGfunc->SetParameter(3, 0.05);
    ALGfunc->SetParameter(4, 10000.0);
    ALGfunc->SetParLimits(0, 0.0, 1.0);
    ALGfunc->SetParLimits(1, 0.0, 5.0);
    ALGfunc->SetParLimits(2, 0.001, 10.0);
    ALGfunc->SetParLimits(3, 0.001, 10.0);

    std::vector<Hist*> vhZtkEXtf = Hist::ProjectAll(HistProj::kY, hZtkEXtf);
    std::vector<Hist*> vhZtfEXtk = Hist::ProjectAll(HistProj::kY, hZtfEXtk);
    std::vector<Hist*>& vhZEX = vhZtkEXtf;
    //std::vector<Hist*>& vhZEX = vhZtfEXtk;
   
    PdfEditor editor(Window(WindowSize::kWideSliceLR), "fitting", "out");
    
    for (int iz = 1; iz <= AXz.nbin(); ++iz) {
        if (iz > 2) continue;
        Hist* hEX = vhZEX.at(iz);
        
        ALGfunc->SetParameter(1, 1.0);
       
        ALGfunc->FixParameter(0, 0.001);
        //ALGfunc->FixParameter(2, 6.80017e-02 * std::exp(-3.90395e-02 * iz * iz) + 0.02);
        //ALGfunc->FixParameter(3, 0.05);

        (*hEX)()->Fit(ALGfunc, "q0", "", 0.85, 1.5);
        (*hEX)()->Fit(ALGfunc, "q0", "", 0.85, 1.5);
        (*hEX)()->Fit(ALGfunc, "q0", "", 0.85, 1.5);
        
        std::cerr << Form("PARAM[0]  %14.8f %14.8f\n", ALGfunc->GetParameter(0), ALGfunc->GetParError(0)); 
        std::cerr << Form("PARAM[1]  %14.8f %14.8f\n", ALGfunc->GetParameter(1), ALGfunc->GetParError(1)); 
        std::cerr << Form("PARAM[2]  %14.8f %14.8f\n", ALGfunc->GetParameter(2), ALGfunc->GetParError(2)); 
        std::cerr << Form("PARAM[3]  %14.8f %14.8f\n", ALGfunc->GetParameter(3), ALGfunc->GetParError(3)); 
    
        (*hZKpa)()->SetBinContent(iz, ALGfunc->GetParameter(0));
        (*hZKpa)()->SetBinError  (iz, ALGfunc->GetParError (0));
        (*hZMpv)()->SetBinContent(iz, ALGfunc->GetParameter(1));
        (*hZMpv)()->SetBinError  (iz, ALGfunc->GetParError (1));
        (*hZSgm)()->SetBinContent(iz, ALGfunc->GetParameter(2));
        (*hZSgm)()->SetBinError  (iz, ALGfunc->GetParError (2));
        (*hZFlc)()->SetBinContent(iz, ALGfunc->GetParameter(3));
        (*hZFlc)()->SetBinError  (iz, ALGfunc->GetParError (3));

        editor.create();
        editor.cd(0, PadAxis(0, 1));
        (*hEX)()->Draw("hist");
        ALGfunc->Draw("l same");
        editor.save();
    }
        
    editor.create();
    editor.cd(0, PadAxis(0, 0));
    hZKpa->draw("pe");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 0));
    hZMpv->draw("pe");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 0));
    hZSgm->draw("pe");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 0));
    hZFlc->draw("pe");
    editor.save();

    editor.close();

    TFile* ofle = new TFile("out/fitting.root", "RECREATE");
    
    hZKpa->write();
    hZMpv->write();
    hZSgm->write();
    hZFlc->write();

    ofle->Write();
    ofle->Close();

    return 1;
}
