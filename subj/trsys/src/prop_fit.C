#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>
    

static constexpr Double_t stable = 2.0;    
static TF1* fgaus = new TF1("fgaus", "gaus", -3.0, 3.0);
static TF1* flg = new TF1("flg", "[0] * TMath::Exp( (1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) )");

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
    
    //Hist::Load("prop_fill.root", "dat");
    Hist::Load("prop_fill.root", "/afs/cern.ch/work/h/hchou/AMSData/test18");

    // Prop
    Hist* hMcx = Hist::Head("hMcx");
    Hist* hMcy = Hist::Head("hMcy");
    Hist* hTcx = Hist::Head("hTcx");
    Hist* hTcy = Hist::Head("hTcy");
    
    Hist* hMux = Hist::Head("hMux");
    Hist* hMuy = Hist::Head("hMuy");
    Hist* hTux = Hist::Head("hTux");
    Hist* hTuy = Hist::Head("hTuy");
    
    Hist* hMee = Hist::Head("hMee");
    Hist* hTee = Hist::Head("hTee");
    
    const Axis& AXeta = hMee->xaxis();
    const Axis& AXrcx = hMcx->yaxis();
    const Axis& AXrcy = hMcy->yaxis();
    const Axis& AXrux = hMux->yaxis();
    const Axis& AXruy = hMuy->yaxis();
    const Axis& AXels = hMee->yaxis();
    
    TFile * ofle = new TFile("prop_fit.root", "RECREATE");
    ofle->cd();
    
    //Hist::AddDirectory();

    std::vector<Hist*> vhMcx = Hist::ProjectAll(HistProj::kY, hMcx);
    Hist* hMcxM = Hist::New("hMcxM", HistAxis(AXeta, "Mean"));
    Hist* hMcxS = Hist::New("hMcxS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Mcx ITER %d\n", it);
        Double_t men = (*vhMcx.at(it))()->GetBinCenter((*vhMcx.at(it))()->GetMaximumBin());
        Double_t rms = (*vhMcx.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhMcx.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhMcx.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhMcx.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhMcx.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hMcxM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hMcxM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hMcxS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hMcxS)()->SetBinError  (it, fgaus->GetParError(2));
        
        Hist* tmpl = Hist::New(Form("hMcx_tmpl%03d", it), HistAxis(AXrcx));
        for (int jt = 1; jt <= AXrcx.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, fgaus->Eval(AXrcx.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    std::vector<Hist*> vhTcx = Hist::ProjectAll(HistProj::kY, hTcx);
    Hist* hTcxM = Hist::New("hTcxM", HistAxis(AXeta, "Mean"));
    Hist* hTcxS = Hist::New("hTcxS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Tcx ITER %d\n", it);
        Double_t men = (*vhTcx.at(it))()->GetBinCenter((*vhTcx.at(it))()->GetMaximumBin());
        Double_t rms = (*vhTcx.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhTcx.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhTcx.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTcx.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTcx.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hTcxM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hTcxM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hTcxS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hTcxS)()->SetBinError  (it, fgaus->GetParError(2));
        
        Hist* tmpl = Hist::New(Form("hTcx_tmpl%03d", it), HistAxis(AXrcx));
        for (int jt = 1; jt <= AXrcx.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, fgaus->Eval(AXrcx.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        Hist* hMs = vhMcx.at(it); hMs->style(Fill(), Line(kBlue, 5, 2), Marker(kBlue));
        Hist* hTs = vhTcx.at(it); hTs->style(Fill(), Line(kRed, 5, 2), Marker(kRed));
        Hist* hM = Hist::Head(Form("hMcx_tmpl%03d", it));
        Hist* hT = Hist::Head(Form("hTcx_tmpl%03d", it));
        hM->style(Fill(), Line(kBlue), Marker(kBlue));
        hT->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chMT = Hist::Collect(Form("chMTcx%03d", it), HistList({ hMs, hTs, hM, hT }));
        chMT->Write();
    }

    hMcxM->style(Fill(), Line(kBlue), Marker(kBlue));
    hTcxM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTcxM = Hist::Collect("chMTcxM", HistList({ hMcxM, hTcxM }));
    chMTcxM->Write();
    
    hMcxS->style(Fill(), Line(kBlue), Marker(kBlue));
    hTcxS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTcxS = Hist::Collect("chMTcxS", HistList({ hMcxS, hTcxS }));
    chMTcxS->Write();


    std::vector<Hist*> vhMcy = Hist::ProjectAll(HistProj::kY, hMcy);
    Hist* hMcyM = Hist::New("hMcyM", HistAxis(AXeta, "Mean"));
    Hist* hMcyS = Hist::New("hMcyS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Mcy ITER %d\n", it);
        Double_t men = (*vhMcy.at(it))()->GetBinCenter((*vhMcy.at(it))()->GetMaximumBin());
        Double_t rms = (*vhMcy.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhMcy.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhMcy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhMcy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhMcy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hMcyM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hMcyM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hMcyS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hMcyS)()->SetBinError  (it, fgaus->GetParError(2));
        
        Hist* tmpl = Hist::New(Form("hMcy_tmpl%03d", it), HistAxis(AXrcy));
        for (int jt = 1; jt <= AXrcy.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, fgaus->Eval(AXrcy.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    std::vector<Hist*> vhTcy = Hist::ProjectAll(HistProj::kY, hTcy);
    Hist* hTcyM = Hist::New("hTcyM", HistAxis(AXeta, "Mean"));
    Hist* hTcyS = Hist::New("hTcyS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Tcy ITER %d\n", it);
        Double_t men = (*vhTcy.at(it))()->GetBinCenter((*vhTcy.at(it))()->GetMaximumBin());
        Double_t rms = (*vhTcy.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhTcy.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhTcy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTcy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTcy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hTcyM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hTcyM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hTcyS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hTcyS)()->SetBinError  (it, fgaus->GetParError(2));
        
        Hist* tmpl = Hist::New(Form("hTcy_tmpl%03d", it), HistAxis(AXrcy));
        for (int jt = 1; jt <= AXrcy.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, fgaus->Eval(AXrcy.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        Hist* hMs = vhMcy.at(it); hMs->style(Fill(), Line(kBlue, 5, 2), Marker(kBlue));
        Hist* hTs = vhTcy.at(it); hTs->style(Fill(), Line(kRed, 5, 2), Marker(kRed));
        Hist* hM = Hist::Head(Form("hMcy_tmpl%03d", it));
        Hist* hT = Hist::Head(Form("hTcy_tmpl%03d", it));
        hM->style(Fill(), Line(kBlue), Marker(kBlue));
        hT->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chMT = Hist::Collect(Form("chMTcy%03d", it), HistList({ hMs, hTs, hM, hT }));
        chMT->Write();
    }

    hMcyM->style(Fill(), Line(kBlue), Marker(kBlue));
    hTcyM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTcyM = Hist::Collect("chMTcyM", HistList({ hMcyM, hTcyM }));
    chMTcyM->Write();
    
    hMcyS->style(Fill(), Line(kBlue), Marker(kBlue));
    hTcyS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTcyS = Hist::Collect("chMTcyS", HistList({ hMcyS, hTcyS }));
    chMTcyS->Write();

    
    std::vector<Hist*> vhMux = Hist::ProjectAll(HistProj::kY, hMux);
    Hist* hMuxM = Hist::New("hMuxM", HistAxis(AXeta, "Mean"));
    Hist* hMuxS = Hist::New("hMuxS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Mux ITER %d\n", it);
        Double_t men = (*vhMux.at(it))()->GetBinCenter((*vhMux.at(it))()->GetMaximumBin());
        Double_t rms = (*vhMux.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhMux.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhMux.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhMux.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhMux.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hMuxM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hMuxM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hMuxS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hMuxS)()->SetBinError  (it, fgaus->GetParError(2));
        
        Hist* tmpl = Hist::New(Form("hMux_tmpl%03d", it), HistAxis(AXrux));
        for (int jt = 1; jt <= AXrux.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, fgaus->Eval(AXrux.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    std::vector<Hist*> vhTux = Hist::ProjectAll(HistProj::kY, hTux);
    Hist* hTuxM = Hist::New("hTuxM", HistAxis(AXeta, "Mean"));
    Hist* hTuxS = Hist::New("hTuxS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Tux ITER %d\n", it);
        Double_t men = (*vhTux.at(it))()->GetBinCenter((*vhTux.at(it))()->GetMaximumBin());
        Double_t rms = (*vhTux.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhTux.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhTux.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTux.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTux.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hTuxM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hTuxM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hTuxS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hTuxS)()->SetBinError  (it, fgaus->GetParError(2));
        
        Hist* tmpl = Hist::New(Form("hTux_tmpl%03d", it), HistAxis(AXrux));
        for (int jt = 1; jt <= AXrux.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, fgaus->Eval(AXrux.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        Hist* hMs = vhMux.at(it); hMs->style(Fill(), Line(kBlue, 5, 2), Marker(kBlue));
        Hist* hTs = vhTux.at(it); hTs->style(Fill(), Line(kRed, 5, 2), Marker(kRed));
        Hist* hM = Hist::Head(Form("hMux_tmpl%03d", it));
        Hist* hT = Hist::Head(Form("hTux_tmpl%03d", it));
        hM->style(Fill(), Line(kBlue), Marker(kBlue));
        hT->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chMT = Hist::Collect(Form("chMTux%03d", it), HistList({ hMs, hTs, hM, hT }));
        chMT->Write();
    }

    hMuxM->style(Fill(), Line(kBlue), Marker(kBlue));
    hTuxM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTuxM = Hist::Collect("chMTuxM", HistList({ hMuxM, hTuxM }));
    chMTuxM->Write();
    
    hMuxS->style(Fill(), Line(kBlue), Marker(kBlue));
    hTuxS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTuxS = Hist::Collect("chMTuxS", HistList({ hMuxS, hTuxS }));
    chMTuxS->Write();


    std::vector<Hist*> vhMuy = Hist::ProjectAll(HistProj::kY, hMuy);
    Hist* hMuyM = Hist::New("hMuyM", HistAxis(AXeta, "Mean"));
    Hist* hMuyS = Hist::New("hMuyS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Muy ITER %d\n", it);
        Double_t men = (*vhMuy.at(it))()->GetBinCenter((*vhMuy.at(it))()->GetMaximumBin());
        Double_t rms = (*vhMuy.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhMuy.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhMuy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhMuy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhMuy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hMuyM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hMuyM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hMuyS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hMuyS)()->SetBinError  (it, fgaus->GetParError(2));
        
        Hist* tmpl = Hist::New(Form("hMuy_tmpl%03d", it), HistAxis(AXruy));
        for (int jt = 1; jt <= AXruy.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, fgaus->Eval(AXruy.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    std::vector<Hist*> vhTuy = Hist::ProjectAll(HistProj::kY, hTuy);
    Hist* hTuyM = Hist::New("hTuyM", HistAxis(AXeta, "Mean"));
    Hist* hTuyS = Hist::New("hTuyS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Tuy ITER %d\n", it);
        Double_t men = (*vhTuy.at(it))()->GetBinCenter((*vhTuy.at(it))()->GetMaximumBin());
        Double_t rms = (*vhTuy.at(it))()->GetRMS();
        fgaus->SetParameters(1000, men, rms);
        
        (*vhTuy.at(it))()->Fit(fgaus, "q0", "", men-stable*rms, men+stable*rms);
        (*vhTuy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTuy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));
        (*vhTuy.at(it))()->Fit(fgaus, "q0", "", fgaus->GetParameter(1)-stable*fgaus->GetParameter(2), fgaus->GetParameter(1)+stable*fgaus->GetParameter(2));

        (*hTuyM)()->SetBinContent(it, fgaus->GetParameter(1));
        (*hTuyM)()->SetBinError  (it, fgaus->GetParError(1));
        (*hTuyS)()->SetBinContent(it, fgaus->GetParameter(2));
        (*hTuyS)()->SetBinError  (it, fgaus->GetParError(2));
        
        Hist* tmpl = Hist::New(Form("hTuy_tmpl%03d", it), HistAxis(AXruy));
        for (int jt = 1; jt <= AXruy.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, fgaus->Eval(AXruy.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        Hist* hMs = vhMuy.at(it); hMs->style(Fill(), Line(kBlue, 5, 2), Marker(kBlue));
        Hist* hTs = vhTuy.at(it); hTs->style(Fill(), Line(kRed, 5, 2), Marker(kRed));
        Hist* hM = Hist::Head(Form("hMuy_tmpl%03d", it));
        Hist* hT = Hist::Head(Form("hTuy_tmpl%03d", it));
        hM->style(Fill(), Line(kBlue), Marker(kBlue));
        hT->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chMT = Hist::Collect(Form("chMTuy%03d", it), HistList({ hMs, hTs, hM, hT }));
        chMT->Write();
    }

    hMuyM->style(Fill(), Line(kBlue), Marker(kBlue));
    hTuyM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTuyM = Hist::Collect("chMTuyM", HistList({ hMuyM, hTuyM }));
    chMTuyM->Write();
    
    hMuyS->style(Fill(), Line(kBlue), Marker(kBlue));
    hTuyS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTuyS = Hist::Collect("chMTuyS", HistList({ hMuyS, hTuyS }));
    chMTuyS->Write();


    std::vector<Hist*> vhMee = Hist::ProjectAll(HistProj::kY, hMee);
    Hist* hMeeK = Hist::New("hMeeK", HistAxis(AXeta, "Kappa"));
    Hist* hMeeM = Hist::New("hMeeM", HistAxis(AXeta, "Mpv"));
    Hist* hMeeS = Hist::New("hMeeS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Mee ITER %d\n", it);
        Double_t mpv = (*vhMee.at(it))()->GetBinCenter((*vhMee.at(it))()->GetMaximumBin());
        Double_t rms = (*vhMee.at(it))()->GetRMS();
        flg->SetParameters(1000, 0.1, mpv, rms);
        flg->SetParLimits(1, 0.0, 1.0);
        flg->SetParLimits(2, 0.0, 10.0*mpv);
        flg->SetParLimits(3, 0.0, 10.0*rms);
        
        (*vhMee.at(it))()->Fit(flg, "q0", "");
        (*vhMee.at(it))()->Fit(flg, "q0", "");
        (*vhMee.at(it))()->Fit(flg, "q0", "");

        (*hMeeK)()->SetBinContent(it, flg->GetParameter(1));
        (*hMeeK)()->SetBinError  (it, flg->GetParError(1));
        (*hMeeM)()->SetBinContent(it, flg->GetParameter(2));
        (*hMeeM)()->SetBinError  (it, flg->GetParError(2));
        (*hMeeS)()->SetBinContent(it, flg->GetParameter(3));
        (*hMeeS)()->SetBinError  (it, flg->GetParError(3));
        
        Hist* tmpl = Hist::New(Form("hMee_tmpl%03d", it), HistAxis(AXels));
        for (int jt = 1; jt <= AXels.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, flg->Eval(AXels.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    std::vector<Hist*> vhTee = Hist::ProjectAll(HistProj::kY, hTee);
    Hist* hTeeK = Hist::New("hTeeK", HistAxis(AXeta, "Kappa"));
    Hist* hTeeM = Hist::New("hTeeM", HistAxis(AXeta, "Mpv"));
    Hist* hTeeS = Hist::New("hTeeS", HistAxis(AXeta, "Sigma"));
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        COUT("Tee ITER %d\n", it);
        Double_t mpv = (*vhTee.at(it))()->GetBinCenter((*vhTee.at(it))()->GetMaximumBin());
        Double_t rms = (*vhTee.at(it))()->GetRMS();
        flg->SetParameters(1000, 0.1, mpv, rms);
        flg->SetParLimits(1, 0.0, 1.0);
        flg->SetParLimits(2, 0.0, 10.0*mpv);
        flg->SetParLimits(3, 0.0, 10.0*rms);
        
        (*vhTee.at(it))()->Fit(flg, "q0", "");
        (*vhTee.at(it))()->Fit(flg, "q0", "");
        (*vhTee.at(it))()->Fit(flg, "q0", "");

        (*hTeeK)()->SetBinContent(it, flg->GetParameter(1));
        (*hTeeK)()->SetBinError  (it, flg->GetParError(1));
        (*hTeeM)()->SetBinContent(it, flg->GetParameter(2));
        (*hTeeM)()->SetBinError  (it, flg->GetParError(2));
        (*hTeeS)()->SetBinContent(it, flg->GetParameter(3));
        (*hTeeS)()->SetBinError  (it, flg->GetParError(3));
        
        Hist* tmpl = Hist::New(Form("hTee_tmpl%03d", it), HistAxis(AXels));
        for (int jt = 1; jt <= AXels.nbin(); ++jt) {
            (*tmpl)()->SetBinContent(jt, flg->Eval(AXels.center(jt)));
            //(*tmpl)()->SetBinError(jt, 1.0e-6);
        }
    }
    
    for (int it = 1; it <= AXeta.nbin(); ++it) {
        Hist* hMs = vhMee.at(it); hMs->style(Fill(), Line(kBlue, 5, 2), Marker(kBlue));
        Hist* hTs = vhTee.at(it); hTs->style(Fill(), Line(kRed, 5, 2), Marker(kRed));
        Hist* hM = Hist::Head(Form("hMee_tmpl%03d", it));
        Hist* hT = Hist::Head(Form("hTee_tmpl%03d", it));
        hM->style(Fill(), Line(kBlue), Marker(kBlue));
        hT->style(Fill(), Line(kRed), Marker(kRed));
        THStack* chMT = Hist::Collect(Form("chMTee%03d", it), HistList({ hMs, hTs, hM, hT }));
        chMT->Write();
    }
    
    hMeeK->style(Fill(), Line(kBlue), Marker(kBlue));
    hTeeK->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTeeK = Hist::Collect("chMTeeK", HistList({ hMeeK, hTeeK }));
    chMTeeK->Write();
    
    hMeeM->style(Fill(), Line(kBlue), Marker(kBlue));
    hTeeM->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTeeM = Hist::Collect("chMTeeM", HistList({ hMeeM, hTeeM }));
    chMTeeM->Write();
    
    hMeeS->style(Fill(), Line(kBlue), Marker(kBlue));
    hTeeS->style(Fill(), Line(kRed), Marker(kRed));
    THStack* chMTeeS = Hist::Collect("chMTeeS", HistList({ hMeeS, hTeeS }));
    chMTeeS->Write();

    ofle->Write();
    ofle->Close();

    return 0;
}
