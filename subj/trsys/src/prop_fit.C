#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

static TF1* flg = new TF1("flg", "[0] * TMath::Exp( (1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/1.78854160900000003e-01) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) )");

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    //Hist::AddDirectory();
    
    Hist::Load("prop_fill.root", "dat");

    // Prop
    //Hist * hMcx = Hist::Head("hMcx");
    //Hist * hMcy = Hist::Head("hMcy");
    //Hist * hTcx = Hist::Head("hTcx");
    //Hist * hTcy = Hist::Head("hTcy");
    //
    //Hist * hMux = Hist::Head("hMux");
    //Hist * hMuy = Hist::Head("hMuy");
    //Hist * hTux = Hist::Head("hTux");
    //Hist * hTuy = Hist::Head("hTuy");
    
    Hist * hMee = Hist::Head("hMee");
    Hist * hTee = Hist::Head("hTee");
    
    const Axis& AXeta = hMee->xaxis();
    const Axis& AXels = hMee->yaxis();
    
    TFile * ofle = new TFile("prop_fit.root", "RECREATE");
    ofle->cd();
    
    Hist::AddDirectory();

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
            (*tmpl)()->SetBinError(jt, 1.0e-6);
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
            (*tmpl)()->SetBinError(jt, 1.0e-6);
        }
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
    







/*
    TGraphErrors* gMcx_men = new TGraphErrors(); gMcx_men->SetNameTitle("gMcx_men", "");
    TGraphErrors* gMcx_sgm = new TGraphErrors(); gMcx_sgm->SetNameTitle("gMcx_sgm", "");
    TGraphErrors* gMcy_men = new TGraphErrors(); gMcy_men->SetNameTitle("gMcy_men", "");
    TGraphErrors* gMcy_sgm = new TGraphErrors(); gMcy_sgm->SetNameTitle("gMcy_sgm", "");
    TGraphErrors* gMux_men = new TGraphErrors(); gMux_men->SetNameTitle("gMux_men", "");
    TGraphErrors* gMux_sgm = new TGraphErrors(); gMux_sgm->SetNameTitle("gMux_sgm", "");
    TGraphErrors* gMuy_men = new TGraphErrors(); gMuy_men->SetNameTitle("gMuy_men", "");
    TGraphErrors* gMuy_sgm = new TGraphErrors(); gMuy_sgm->SetNameTitle("gMuy_sgm", "");
    TGraphErrors* gTcx_men = new TGraphErrors(); gTcx_men->SetNameTitle("gTcx_men", "");
    TGraphErrors* gTcx_sgm = new TGraphErrors(); gTcx_sgm->SetNameTitle("gTcx_sgm", "");
    TGraphErrors* gTcy_men = new TGraphErrors(); gTcy_men->SetNameTitle("gTcy_men", "");
    TGraphErrors* gTcy_sgm = new TGraphErrors(); gTcy_sgm->SetNameTitle("gTcy_sgm", "");
    TGraphErrors* gTux_men = new TGraphErrors(); gTux_men->SetNameTitle("gTux_men", "");
    TGraphErrors* gTux_sgm = new TGraphErrors(); gTux_sgm->SetNameTitle("gTux_sgm", "");
    TGraphErrors* gTuy_men = new TGraphErrors(); gTuy_men->SetNameTitle("gTuy_men", "");
    TGraphErrors* gTuy_sgm = new TGraphErrors(); gTuy_sgm->SetNameTitle("gTuy_sgm", "");
   
    const Double_t stable = 2.0;
    TF1 * gaus = new TF1("gaus", "gaus", -3.0, 3.0);
    //TF1 * f2gs = new TF1("f2gs", "([1]>0)*([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([3]>0)*([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3])");
    //f2gs->SetParLimits(1, 0, 100);
    //f2gs->SetParLimits(3, 0, 100);
    TF1 * f3gs = new TF1("f3gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5])");
    f3gs->SetParLimits(1, 0, 10);
    f3gs->SetParLimits(3, 0, 10);
    f3gs->SetParLimits(5, 0, 10);
    std::vector<Hist*> vhMcx = Hist::ProjectAll(HistProj::kY, hMcx);
    std::vector<Hist*> vhMcy = Hist::ProjectAll(HistProj::kY, hMcy);
    std::vector<Hist*> vhMux = Hist::ProjectAll(HistProj::kY, hMux);
    std::vector<Hist*> vhMuy = Hist::ProjectAll(HistProj::kY, hMuy);
    std::vector<Hist*> vhTcx = Hist::ProjectAll(HistProj::kY, hTcx);
    std::vector<Hist*> vhTcy = Hist::ProjectAll(HistProj::kY, hTcy);
    std::vector<Hist*> vhTux = Hist::ProjectAll(HistProj::kY, hTux);
    std::vector<Hist*> vhTuy = Hist::ProjectAll(HistProj::kY, hTuy);
    for (int it = 1; it <= AXmom.nbin(); ++it) {
        double mom = AXmom.center(it+1, AxisScale::kLog);
        PhySt part(PartType::Proton);
        part.set_mom(mom);
        //Double_t val = part.gmbta();
        //Double_t val = part.gm();
        Double_t val = part.bta();
        
        gaus->SetParameters(1000, 0, (*vhMcx.at(it))()->GetRMS());
        (*vhMcx.at(it))()->Fit(gaus, "q0", "");
        (*vhMcx.at(it))()->Fit(gaus, "q0", "", -stable*gaus->GetParameter(2), stable*gaus->GetParameter(2));
    
        //f2gs->SetParameters(10., (*vhMcx.at(it))()->GetRMS(), 10., 3*(*vhMcx.at(it))()->GetRMS());
        //(*vhMcx.at(it))()->Fit(f2gs, "q0", "");
        //(*vhMcx.at(it))()->Fit(f2gs, "q0", "");

        gMcx_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gMcx_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gMcx_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gMcx_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        //gMcx_sgm->SetPoint     (it-1, val, std::min(f2gs->GetParameter(1), f2gs->GetParameter(3)));
        //gMcx_sgm->SetPointError(it-1,  0., 0);
        
        gaus->SetParameters(1000, 0, (*vhMcy.at(it))()->GetRMS());
        (*vhMcy.at(it))()->Fit(gaus, "q0", "");
        (*vhMcy.at(it))()->Fit(gaus, "q0", "", -stable*gaus->GetParameter(2), stable*gaus->GetParameter(2));
        
        gMcy_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gMcy_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gMcy_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gMcy_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhMux.at(it))()->GetRMS());
        (*vhMux.at(it))()->Fit(gaus, "q0", "");
        (*vhMux.at(it))()->Fit(gaus, "q0", "", -stable*gaus->GetParameter(2), stable*gaus->GetParameter(2));
        
        //f3gs->SetParameters(100., (*vhMux.at(it))()->GetRMS(), 30., 2*(*vhMux.at(it))()->GetRMS(), 10., 4*(*vhMux.at(it))()->GetRMS());
        //(*vhMux.at(it))()->Fit(f3gs, "q0", "");
        //(*vhMux.at(it))()->Fit(f3gs, "q0", "");
        //(*vhMux.at(it))()->Fit(f3gs, "q0", "");
        //(*vhMux.at(it))()->Fit(f3gs, "q0", "");
        //(*vhMux.at(it))()->Fit(f3gs, "q0", "");
        
        gMux_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gMux_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gMux_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gMux_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        //gMux_sgm->SetPoint     (it-1, val, std::min(std::min(f3gs->GetParameter(1), f3gs->GetParameter(3)), f3gs->GetParameter(5)));
        //gMux_sgm->SetPointError(it-1,  0., 0);
        
        gaus->SetParameters(1000, 0, (*vhMuy.at(it))()->GetRMS());
        (*vhMuy.at(it))()->Fit(gaus, "q0", "");
        (*vhMuy.at(it))()->Fit(gaus, "q0", "", -stable*gaus->GetParameter(2), stable*gaus->GetParameter(2));
        gMuy_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gMuy_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gMuy_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gMuy_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhTcx.at(it))()->GetRMS());
        (*vhTcx.at(it))()->Fit(gaus, "q0", "");
        (*vhTcx.at(it))()->Fit(gaus, "q0", "", -stable*gaus->GetParameter(2), stable*gaus->GetParameter(2));
        gTcx_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gTcx_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gTcx_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gTcx_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhTcy.at(it))()->GetRMS());
        (*vhTcy.at(it))()->Fit(gaus, "q0", "");
        (*vhTcy.at(it))()->Fit(gaus, "q0", "", -stable*gaus->GetParameter(2), stable*gaus->GetParameter(2));
        gTcy_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gTcy_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gTcy_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gTcy_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhTux.at(it))()->GetRMS());
        (*vhTux.at(it))()->Fit(gaus, "q0", "");
        (*vhTux.at(it))()->Fit(gaus, "q0", "", -stable*gaus->GetParameter(2), stable*gaus->GetParameter(2));
        gTux_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gTux_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gTux_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gTux_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
        
        gaus->SetParameters(1000, 0, (*vhTuy.at(it))()->GetRMS());
        (*vhTuy.at(it))()->Fit(gaus, "q0", "");
        (*vhTuy.at(it))()->Fit(gaus, "q0", "", -stable*gaus->GetParameter(2), stable*gaus->GetParameter(2));
        gTuy_men->SetPoint     (it-1, val, gaus->GetParameter(1));
        gTuy_men->SetPointError(it-1,  0., gaus->GetParError(1));
        gTuy_sgm->SetPoint     (it-1, val, gaus->GetParameter(2));
        gTuy_sgm->SetPointError(it-1,  0., gaus->GetParError(2));
    } 

    gMcx_men->Write();
    gMcx_sgm->Write();
    gMcy_men->Write();
    gMcy_sgm->Write();
    gMux_men->Write();
    gMux_sgm->Write();
    gMuy_men->Write();
    gMuy_sgm->Write();
    gTcx_men->Write();
    gTcx_sgm->Write();
    gTcy_men->Write();
    gTcy_sgm->Write();
    gTux_men->Write();
    gTux_sgm->Write();
    gTuy_men->Write();
    gTuy_sgm->Write();

    TGraphErrors* gMee_kpa = new TGraphErrors(); gMee_kpa->SetNameTitle("gMee_kpa",  "");
    TGraphErrors* gMee_mpv = new TGraphErrors(); gMee_mpv->SetNameTitle("gMee_mpv",  "");
    TGraphErrors* gMee_sgm = new TGraphErrors(); gMee_sgm->SetNameTitle("gMee_sgm",  "");
    TGraphErrors* gMee_mos = new TGraphErrors(); gMee_mos->SetNameTitle("gMee_mos",  "");
    TGraphErrors* gTee_kpa = new TGraphErrors(); gTee_kpa->SetNameTitle("gTee_kpa",  "");
    TGraphErrors* gTee_mpv = new TGraphErrors(); gTee_mpv->SetNameTitle("gTee_mpv",  "");
    TGraphErrors* gTee_sgm = new TGraphErrors(); gTee_sgm->SetNameTitle("gTee_sgm",  "");
    TGraphErrors* gTee_mos = new TGraphErrors(); gTee_mos->SetNameTitle("gTee_mos",  "");

    TF1 * feloss = new TF1("feloss", "[0] * TMath::Power( ([2]/x)/[1]/[1], ([2]/x)/[1]/[1] ) / TMath::Gamma( ([2]/x)/[1]/[1] ) * TMath::Exp(-(([2]/x)/[1]/[1]) * ((x-[2])/[3] + TMath::Exp(-(x-[2])/[3])) )");

    std::vector<Hist*> vhMee = Hist::ProjectAll(HistProj::kY, hMee);
    std::vector<Hist*> vhTee = Hist::ProjectAll(HistProj::kY, hTee);
    for (int it = 1; it <= AXmom.nbin(); ++it) {
        double mom = AXmom.center(it+1, AxisScale::kLog);
        PhySt part(PartType::Proton);
        part.set_mom(mom);
        Double_t val = part.gmbta();
        //Double_t val = part.gm();
        ///Double_t val = part.bta();
        
        Double_t sqr_bta = (part.bta() * part.bta());
        Double_t eloss_ion_kpa  = MGMath::SQRT_TWO / (MGMath::HALF * (sqr_bta - MGMath::ONE) + MGMath::ONE/sqr_bta);
        Double_t eloss_ion_kpa2 = MGMath::SQRT_TWO / (-2.83667e+00 * std::pow(std::fabs(MGMath::ONE - sqr_bta), 1.37299e+00) + std::pow(MGMath::ONE/sqr_bta, 1.56239e+00)) + 0.0859215;
       
        int    mpbin = (*vhMee.at(it))()->GetMaximumBin();
        int    mpval = (*vhMee.at(it))()->GetBinContent(mpbin);
        double mpeak = (*vhMee.at(it))()->GetXaxis()->GetBinCenter(mpbin);
        
        int mmin = mpbin, mmax = mpbin;
        for (Int_t bin = mpbin; bin >= 1;                              --bin) { if ((*vhMee.at(it))()->GetBinContent(bin) < 1.5e-1*mpval) { mmin = bin; break; } }
        for (Int_t bin = mpbin; bin <= (*vhMee.at(it))()->GetNbinsX(); ++bin) { if ((*vhMee.at(it))()->GetBinContent(bin) < 5.0e-3*mpval) { mmax = bin; break; } }
        double mbdmin = vhMee.at(it)->xaxis().center(mmin);
        double mbdmax = vhMee.at(it)->xaxis().center(mmax);
        
        feloss->SetParameters(1000., 1.0, mpeak, 0.1*mpeak, 0.5);
        feloss->FixParameter(1, eloss_ion_kpa2);
        (*vhMee.at(it))()->Fit(feloss, "q0", "", mbdmin, mbdmax);

        if (feloss->GetParameter(1) < 0) feloss->SetParameter(1, -feloss->GetParameter(1));
        if (feloss->GetParameter(3) < 0) feloss->SetParameter(3, -feloss->GetParameter(3));
        feloss->ReleaseParameter(1);
        (*vhMee.at(it))()->Fit(feloss, "q0", "", mbdmin, mbdmax);

        if (feloss->GetParameter(1) < 0) feloss->SetParameter(1, -feloss->GetParameter(1));
        if (feloss->GetParameter(3) < 0) feloss->SetParameter(3, -feloss->GetParameter(3));
        (*vhMee.at(it))()->Fit(feloss, "q0", "", mbdmin, mbdmax);

        gMee_kpa->SetPoint     (it-1, val, feloss->GetParameter(1));
        gMee_kpa->SetPointError(it-1,  0., feloss->GetParError(1));
        gMee_mpv->SetPoint     (it-1, val, feloss->GetParameter(2));
        gMee_mpv->SetPointError(it-1,  0., feloss->GetParError(2));
        gMee_sgm->SetPoint     (it-1, val, feloss->GetParameter(3));
        gMee_sgm->SetPointError(it-1,  0., feloss->GetParError(3));
        gMee_mos->SetPoint     (it-1, val, feloss->GetParameter(2)/feloss->GetParameter(3));
        gMee_mos->SetPointError(it-1,  0., 0.);
       
        int    tpbin = (*vhTee.at(it))()->GetMaximumBin();
        int    tpval = (*vhTee.at(it))()->GetBinContent(tpbin);
        double tpeak = (*vhTee.at(it))()->GetXaxis()->GetBinCenter(tpbin);
        
        int tmin = tpbin, tmax = tpbin;
        for (Int_t bin = tpbin; bin >= 1;                              --bin) { if ((*vhTee.at(it))()->GetBinContent(bin) < 2.0e-1*tpval) { tmin = bin; break; } }
        for (Int_t bin = tpbin; bin <= (*vhTee.at(it))()->GetNbinsX(); ++bin) { if ((*vhTee.at(it))()->GetBinContent(bin) < 5.0e-2*tpval) { tmax = bin; break; } }
        double tbdmin = vhTee.at(it)->xaxis().center(tmin);
        double tbdmax = vhTee.at(it)->xaxis().center(tmax);
        
        feloss->SetParameters(1000., 1.0, tpeak, 0.1*tpeak, 0.5);
        feloss->FixParameter(1, eloss_ion_kpa2);
        (*vhTee.at(it))()->Fit(feloss, "q0", "", tbdmin, tbdmax);

        if (feloss->GetParameter(1) < 0) feloss->SetParameter(1, -feloss->GetParameter(1));
        if (feloss->GetParameter(3) < 0) feloss->SetParameter(3, -feloss->GetParameter(3));
        feloss->ReleaseParameter(1);
        (*vhTee.at(it))()->Fit(feloss, "q0", "", tbdmin, tbdmax);
        
        if (feloss->GetParameter(1) < 0) feloss->SetParameter(1, -feloss->GetParameter(1));
        if (feloss->GetParameter(3) < 0) feloss->SetParameter(3, -feloss->GetParameter(3));
        (*vhTee.at(it))()->Fit(feloss, "q0", "", tbdmin, tbdmax);
        
        gTee_kpa->SetPoint     (it-1, val, feloss->GetParameter(1));
        gTee_kpa->SetPointError(it-1,  0., feloss->GetParError(1));
        gTee_mpv->SetPoint     (it-1, val, feloss->GetParameter(2));
        gTee_mpv->SetPointError(it-1,  0., feloss->GetParError(2));
        gTee_sgm->SetPoint     (it-1, val, feloss->GetParameter(3));
        gTee_sgm->SetPointError(it-1,  0., feloss->GetParError(3));
        gTee_mos->SetPoint     (it-1, val, feloss->GetParameter(2)/feloss->GetParameter(3));
        gTee_mos->SetPointError(it-1,  0., 0.);
    }

    gMee_kpa->Write();
    gMee_mpv->Write();
    gMee_sgm->Write();
    gMee_mos->Write();
    gTee_kpa->Write();
    gTee_mpv->Write();
    gTee_sgm->Write();
    gTee_mos->Write();
    
    TGraphErrors* gMTcx = new TGraphErrors(); gMTcx->SetNameTitle("gMTcx", "");
    TGraphErrors* gMTcy = new TGraphErrors(); gMTcy->SetNameTitle("gMTcy", "");
    TGraphErrors* gMTux = new TGraphErrors(); gMTux->SetNameTitle("gMTux", "");
    TGraphErrors* gMTuy = new TGraphErrors(); gMTuy->SetNameTitle("gMTuy", "");
    TGraphErrors* gMTkpa = new TGraphErrors(); gMTkpa->SetNameTitle("gMTkpa", "");
    TGraphErrors* gMTmpv = new TGraphErrors(); gMTmpv->SetNameTitle("gMTmpv", "");
    TGraphErrors* gMTsgm = new TGraphErrors(); gMTsgm->SetNameTitle("gMTsgm", "");
    TGraphErrors* gMTmos = new TGraphErrors(); gMTmos->SetNameTitle("gMTmos", "");
    for (int it = 0; it < AXmom.nbin(); ++it) {
        gMTcx->SetPoint(it, gMcx_sgm->GetX()[it], gMcx_sgm->GetY()[it]/gTcx_sgm->GetY()[it]);
        gMTcy->SetPoint(it, gMcy_sgm->GetX()[it], gMcy_sgm->GetY()[it]/gTcy_sgm->GetY()[it]);
        gMTux->SetPoint(it, gMux_sgm->GetX()[it], gMux_sgm->GetY()[it]/gTux_sgm->GetY()[it]);
        gMTuy->SetPoint(it, gMuy_sgm->GetX()[it], gMuy_sgm->GetY()[it]/gTuy_sgm->GetY()[it]);
        gMTkpa->SetPoint(it, gMee_kpa->GetX()[it], gMee_kpa->GetY()[it]/gTee_kpa->GetY()[it]);
        gMTmpv->SetPoint(it, gMee_mpv->GetX()[it], gMee_mpv->GetY()[it]/gTee_mpv->GetY()[it]);
        gMTsgm->SetPoint(it, gMee_sgm->GetX()[it], gMee_sgm->GetY()[it]/gTee_sgm->GetY()[it]);
        gMTmos->SetPoint(it, gMee_mos->GetX()[it], gMee_mos->GetY()[it]/gTee_mos->GetY()[it]);
    }
    
    gMTcx ->Write(); 
    gMTcy ->Write(); 
    gMTux ->Write(); 
    gMTuy ->Write(); 
    gMTkpa->Write();
    gMTmpv->Write();
    gMTsgm->Write();
    gMTmos->Write();
*/
    ofle->Write();
    ofle->Close();

    return 0;
}
