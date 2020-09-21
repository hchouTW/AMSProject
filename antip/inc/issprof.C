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
    std::string subv = "38";
   
    // ref. li yang p.75
    const double acc_geom_factor = 47.78; // (m^2 sr)
    UInt_t cntev = 0;
    
    UInt_t cntp = 0;
    TFile* fmcp = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcpr%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcp = (TH1D*)fmcp->Get("hPROF_HC_mass_cnt_MC");
    TTree* tmcp = (TTree*)fmcp->Get("ana");
    tmcp->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcp->GetEntries(); ++it) { tmcp->GetEntry(it); cntp+=cntev; }
    double num_fact_p = static_cast<double>(cntp) / (std::log(1000.0) - std::log(0.2));
    
    UInt_t cntd = 0;
    TFile* fmcd = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcd%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcd = (TH1D*)fmcd->Get("hPROF_HC_mass_cnt_MC");
    TTree* tmcd = (TTree*)fmcd->Get("ana");
    tmcd->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcd->GetEntries(); ++it) { tmcd->GetEntry(it); cntd+=cntev; }
    double num_fact_d = static_cast<double>(cntd) / (std::log(1000.0) - std::log(0.2));
   
    double reg_p[2] = { 0.5417, 1.2510 };
    double sum_pas_p = 0;
    double sum_gen_p = 0;
    TH1D* haccp = new TH1D("haccp", "", hmcp->GetXaxis()->GetNbins(), hmcp->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= haccp->GetXaxis()->GetNbins(); ++ib) {
        double cnt = num_fact_p * (std::log(hmcp->GetXaxis()->GetBinUpEdge(ib)) - std::log(hmcp->GetXaxis()->GetBinLowEdge(ib)));
        double err = hmcp->GetBinError(ib) / hmcp->GetBinContent(ib);
        double acc = hmcp->GetBinContent(ib) / cnt;
        if (!std::isfinite(acc) || acc <= 0.0) continue;
        haccp->SetBinContent(ib, acc_geom_factor * acc);
        haccp->SetBinError  (ib, acc_geom_factor * acc * err);

        if (hmcp->GetXaxis()->GetBinLowEdge(ib) > reg_p[0] && hmcp->GetXaxis()->GetBinUpEdge(ib) < reg_p[1]) {
            sum_pas_p += hmcp->GetBinContent(ib);
            sum_gen_p += cnt;
        }
    }
    double acc_p = acc_geom_factor * sum_pas_p / sum_gen_p;
    
    double reg_d[2] = { 1.0832, 2.5015 };
    double sum_pas_d = 0;
    double sum_gen_d = 0;
    TH1D* haccd = new TH1D("haccd", "", hmcd->GetXaxis()->GetNbins(), hmcd->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= haccd->GetXaxis()->GetNbins(); ++ib) {
        double cnt = num_fact_d * (std::log(hmcd->GetXaxis()->GetBinUpEdge(ib)) - std::log(hmcd->GetXaxis()->GetBinLowEdge(ib)));
        double err = hmcd->GetBinError(ib) / hmcd->GetBinContent(ib);
        double acc = hmcd->GetBinContent(ib) / cnt;
        if (!std::isfinite(acc) || acc <= 0.0) continue;
        haccd->SetBinContent(ib, acc_geom_factor * acc);
        haccd->SetBinError  (ib, acc_geom_factor * acc * err);
        
        if (hmcd->GetXaxis()->GetBinLowEdge(ib) > reg_d[0] && hmcd->GetXaxis()->GetBinUpEdge(ib) < reg_d[1]) {
            sum_pas_d += hmcd->GetBinContent(ib);
            sum_gen_d += cnt;
        }
    }
    double acc_d = acc_geom_factor * sum_pas_d / sum_gen_d;

    std::cerr << Form("\n-------------------------------\n");
    std::cerr << Form("Acc P  %14.8f\n", acc_p);
    std::cerr << Form("Acc D  %14.8f\n", acc_d);
    std::cerr << Form("-------------------------------\n");

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/iss%s", subv.c_str()));

    Hist* hMassHC_pr = Hist::New("hPROF_HC_mass_pr", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcllpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_MC_FLUX10"));
    //Hist* hMassHC_pr = Hist::New("hPROF_HC_mass_pr", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_MC_FLUX10"));
    //Hist* hMassHC_ap = Hist::New("hPROF_HC_mass_ap", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcap%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_MC_FLUX10"));
    Hist* hMassHC_ap = Hist::New("hPROF_HC_mass_ap", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcap%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_MC_FLUX10"));
    Hist* hMassHC_d  = Hist::New("hPROF_HC_mass_d",  (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcd%s/YiMdst.root",  subv.c_str()))->Get("hPROF_HC_mass_MC_FLUX10"));
    
    Hist* hMassHC_pr_pure = Hist::New("hPROF_HC_mass_pr_pure", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcllpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_MC_PURE_FLUX10"));
    Hist* hMassHC_ap_pure = Hist::New("hPROF_HC_mass_ap_pure", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcap%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_MC_PURE_FLUX10"));
    Hist* hMassHC_d_pure  = Hist::New("hPROF_HC_mass_d_pure",  (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcd%s/YiMdst.root",  subv.c_str()))->Get("hPROF_HC_mass_MC_PURE_FLUX10"));
    
    Hist* hMassHC_pip = Hist::New("hPROF_HC_mass_pip", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcllpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_PI_POS_MC_FLUX10"));
    Hist* hMassHC_pin = Hist::New("hPROF_HC_mass_pin", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcllpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_PI_NEG_MC_FLUX10"));
    Hist* hMassHC_kp = Hist::New("hPROF_HC_mass_kp", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcllpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_K_POS_MC_FLUX10"));
    Hist* hMassHC_kn = Hist::New("hPROF_HC_mass_kn", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcllpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_K_NEG_MC_FLUX10"));

    //Hist* hMassHC_cnt_pr = Hist::New("hPROF_HC_mass_cnt_pr", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcllpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_cnt_MC_FLUX10"));
    Hist* hMassHC_cnt_pr = Hist::New("hPROF_HC_mass_cnt_pr", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcpr%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_cnt_MC_FLUX10"));
    Hist* hMassHC_cnt_ap = Hist::New("hPROF_HC_mass_cnt_ap", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcap%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_cnt_MC_FLUX10"));
    Hist* hMassHC_cnt_d  = Hist::New("hPROF_HC_mass_cnt_d",  (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcd%s/YiMdst.root",  subv.c_str()))->Get("hPROF_HC_mass_cnt_MC_FLUX10"));
    
    Hist* hMassHC_cnt_iss = Hist::New("hPROF_HC_mass_cnt_iss", (TH1D*) TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/iss%s/YiMdst.root", subv.c_str()))->Get("hPROF_HC_mass_cnt_ISS"));
    
    

    PdfEditor editor(Window(), "issprof", "out", PdfEditor::Type::kPDF);

    const Axis& AXrig = Hist::Head("hPROF_CK_mass_ISS")->xaxis();

    Hist* hMassCK = Hist::Head("hPROF_CK_mass_ISS");
    Hist* hMassKF = Hist::Head("hPROF_KF_mass_ISS");
    Hist* hMassHC  = Hist::Head("hPROF_HC_mass_ISS");
    Hist* hMassHC2 = Hist::Head("hPROF_HC_mass_ISS_CF");

    Hist* hRvarCK = Hist::Head("hPROF_CK_rvar_cut");
    Hist* hRvarHC = Hist::Head("hPROF_HC_rvar_cut");

    (*hRvarCK)()->Scale(1.0/(*hRvarCK)()->Integral());
    (*hRvarHC)()->Scale(1.0/(*hRvarHC)()->Integral());

    //THStack* stackM = Hist::Collect("hPROF_mass", HistList({ hMassCK, hMassKF, hMassHC }));
    THStack* stackM = Hist::Collect("hPROF_mass", HistList({ hMassCK, hMassHC }));
    stackM->SetTitle("|Z|=1 & (0.50 < #beta < 0.80)");
    
    THStack* stackM2 = Hist::Collect("hPROF_mass2", HistList({ hMassHC, hMassHC2 }));
    stackM2->SetTitle("|Z|=1 & (0.50 < #beta < 0.80)");

    hMassCK->style(Line(kBlue, 0, 1.0),  Marker(kBlue,  MarkerStyle(MarkerShape::kCircle), 0.75));
    hMassKF->style(Line(kGreen+1, 0, 1.0),  Marker(kGreen+1,  MarkerStyle(MarkerShape::kCircle), 0.75));
    hMassHC->style(Line(kRed, 0, 1.0),  Marker(kRed,  MarkerStyle(MarkerShape::kCircle), 0.75));
    hMassHC2->style(Line(kRed+2, 0, 1.0),  Marker(kRed+2,  MarkerStyle(MarkerShape::kCircle), 0.75));

    editor.create();
    editor.cd(0, PadAxis(0, 1));
    stackM->Draw("nostack hist");
    stackM->GetHistogram()->SetLineColor(0);
    stackM->GetHistogram()->SetMarkerColor(0);
    stackM->GetHistogram()->GetXaxis()->SetTitle("Mass/Z [GV/c^{2}]");
    stackM->GetHistogram()->GetYaxis()->SetTitle("Events / Bin");
    stackM->SetMinimum(0.1);
    stackM->Draw("nostack pe");
    Legend legM("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.70, 0.85));
    legM()->AddEntry((TObject*)nullptr, stackM->GetTitle(), "");
    //legM()->AddEntry((*hMassKF)(), Form("Official fitting %.2f B", (*hMassKF)()->Integral()*1.0e-9), "lp");
    //legM()->AddEntry((*hMassCK)(), Form("Official fitting %.2f B", (*hMassCK)()->Integral()*1.0e-9), "lp");
    //legM()->AddEntry((*hMassHC)(), Form("FPM fitting %.2f B", (*hMassHC)()->Integral()*1.0e-9), "lp");
    legM()->AddEntry((*hMassCK)(), Form("Official fitting"), "lp");
    legM()->AddEntry((*hMassHC)(), Form("FPM fitting"), "lp");
    legM()->SetTextFont(43);
    legM()->SetTextSize(15);
    legM()->SetFillColor(0);
    legM.draw();
    TextDraw("D^{+}"  , TextStyle(kBlack, 0.04), TextAlign(0.74, 0.77));
    TextDraw("p^{+}"  , TextStyle(kBlack, 0.04), TextAlign(0.62, 0.78));
    TextDraw("#pi^{+}", TextStyle(kBlack, 0.04), TextAlign(0.52, 0.75));
    TextDraw("#pi^{-}", TextStyle(kBlack, 0.04), TextAlign(0.48, 0.75));
    TextDraw("K^{-}"  , TextStyle(kBlack, 0.04), TextAlign(0.44, 0.52));
    TextDraw("p^{-}"  , TextStyle(kBlack, 0.04), TextAlign(0.38, 0.47));
    editor.save();

    editor.create();
    editor.cd(0, PadAxis(0, 1));
    stackM2->Draw("nostack hist");
    stackM2->GetHistogram()->SetLineColor(0);
    stackM2->GetHistogram()->SetMarkerColor(0);
    stackM2->GetHistogram()->GetXaxis()->SetTitle("Mass/Z [GV/c^{2}]");
    stackM2->GetHistogram()->GetYaxis()->SetTitle("Events / Bin");
    stackM2->SetMinimum(0.1);
    stackM2->Draw("nostack pe");
    Legend legM2("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.75, 0.85));
    legM2()->AddEntry((TObject*)nullptr, stackM2->GetTitle(), "");
    legM2()->AddEntry((*hMassHC)(), Form("Without Cutoff %.2f B", (*hMassHC)()->Integral()*1.0e-9), "lp");
    legM2()->AddEntry((*hMassHC2)(), Form("With Cutoff %.2f B", (*hMassHC2)()->Integral()*1.0e-9), "lp");
    legM2()->SetTextFont(43);
    legM2()->SetTextSize(15);
    legM2()->SetFillColor(0);
    legM2.draw();
    //editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    (*hMassHC)()->Draw("pe");
    (*hMassHC)()->SetMinimum(0.1);
    (*hMassHC)()->Draw("pe");
    Legend legM3("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.75, 0.85));
    legM3()->AddEntry((TObject*)nullptr, "|Z|=1 & (0.50 < #beta < 0.80)", "");
    legM3()->AddEntry((*hMassHC)(), Form("Events %.2f B", (*hMassHC)()->Integral()*1.0e-9), "lp");
    legM3()->SetTextFont(43);
    legM3()->SetTextSize(15);
    legM3()->SetFillColor(0);
    legM3.draw();
    //editor.save();
  
    TF1* gaus = new TF1("gaus", "gaus");
    
    ((TH1D*)((*hRvarCK)()->Clone()))->Fit(gaus, "q0", "", -0.5, 0.5);
    double sgmRvarCK = gaus->GetParameter(2);

    ((TH1D*)((*hRvarHC)()->Clone()))->Fit(gaus, "q0", "", -0.5, 0.5);
    double sgmRvarHC = gaus->GetParameter(2);

    THStack* stackR = Hist::Collect("hPROF_rvar", HistList({ hRvarCK, hRvarHC }));
    //stackR->SetTitle("|Z|=1 & |R_{L1L9}| > 100 [GV/c]");
    stackR->SetTitle("|Z|=1 & |R| > 100 [GV/c]");

    hRvarCK->style(Line(kBlue, 0, 1.5),  Marker(kBlue,  MarkerStyle(MarkerShape::kCircle  )));
    hRvarHC->style(Line(kRed, 0, 1.5),  Marker(kRed,  MarkerStyle(MarkerShape::kCircle  )));

    editor.create();
    editor.cd(0, PadAxis(0, 0));
    stackR->Draw("nostack hist");
    stackR->GetHistogram()->SetLineColor(0);
    stackR->GetHistogram()->SetMarkerColor(0);
    stackR->GetHistogram()->GetXaxis()->SetTitle("Scaled(1/R_{L9}-1/R_{L1})");
    stackR->GetHistogram()->GetYaxis()->SetTitle("Probability");
    stackR->Draw("nostack hist");
    Legend legR("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.35, 0.70, 0.85));
    legR()->AddEntry((TObject*)nullptr, stackR->GetTitle(), "");
    //legR()->AddEntry((TObject*)nullptr, "Gaussian fitting in [-0.5, 0.5]", "");
    legR()->AddEntry((*hRvarCK)(), Form("Official fitting (#sigma = %6.4f)", sgmRvarCK), "l");
    legR()->AddEntry((*hRvarHC)(), Form("FPM fitting (#sigma = %6.4f)", sgmRvarHC), "l");
    legR()->SetTextFont(43);
    legR()->SetTextSize(15);
    legR()->SetFillColor(0);
    legR.draw();
    editor.save();

    std::cerr << Form("CHECK %d   %d %d %d\n", hMassHC==0, hMassHC_pr==0, hMassHC_ap==0, hMassHC_d==0);

    // Fit
    Fit::RooVar roovar("sqrm", hMassHC, HistList({ hMassHC_pr, hMassHC_d, hMassHC_ap }));
    Fit::RooSysResult rlt(roovar, true, 100);
    Fit::RooVar var = rlt.var();
    Fit::RooPar stdpar = rlt.std_par();
    Fit::RooPar syspar = rlt.sys_par();

    std::cerr << Form("ROOFIT 1) %14.8f ERR %14.8f %14.8f\n", stdpar.val(0), stdpar.err(0), syspar.err(0));
    std::cerr << Form("ROOFIT 2) %14.8f ERR %14.8f %14.8f\n", stdpar.val(1), stdpar.err(1), syspar.err(1));
    std::cerr << Form("ROOFIT 3) %14.8f ERR %14.8f %14.8f\n", stdpar.val(2), stdpar.err(2), syspar.err(2));

    THStack* stackM_fit = Hist::Collect("hPROF_mass_fit", HistList({ var.temp(0), var.temp(1), var.temp(2) }));
    stackM_fit->SetTitle("|Z|=1 & (0.50 < #beta < 0.80)");

    hMassHC->style(Line(kBlack, 0, 1.0),  Marker(kBlack,  MarkerStyle(MarkerShape::kCircle ), 0.75));
    
    var.temp(0)->style(Line(kBlue, 0, 2.0),  Marker(kBlue,  MarkerStyle(MarkerShape::kCircle )), Fill(kBlue));
    var.temp(1)->style(Line(kGreen+1, 0, 2.0),  Marker(kGreen+1,  MarkerStyle(MarkerShape::kCircle )), Fill(kGreen+1));
    var.temp(2)->style(Line(kRed, 0, 2.0),  Marker(kRed,  MarkerStyle(MarkerShape::kCircle )), Fill(kRed));

    editor.create();
    editor.cd(0, PadAxis(0, 1));
    (*hMassHC)()->Draw("pe");
    (*hMassHC)()->SetMinimum(0.1);
    (*hMassHC)()->GetXaxis()->SetTitle("Mass/Z [GV/c^{2}]");
    (*hMassHC)()->GetYaxis()->SetTitle("Events / Bin");
    stackM_fit->Draw("nostack hist same");
    //stackM_fit->Draw("stack hist same"); // testcode
    (*hMassHC)()->Draw("pe same");
    Legend legM_fit("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.70, 0.85));
    legM_fit()->AddEntry((TObject*)nullptr, stackM_fit->GetTitle(), "");
    legM_fit()->AddEntry((*hMassHC)(), Form("Data (05/24/2011 ~ 12/14/2018) %.2f B", (*hMassHC)()->Integral()*1.0e-9), "lp");
    legM_fit()->AddEntry((*var.temp(0))(), Form("Proton Template %.2f B", (*var.temp(0))()->Integral()*1.0e-9), "l");
    legM_fit()->AddEntry((*var.temp(1))(), Form("Deuteron Template %.2f M", (*var.temp(1))()->Integral()*1.0e-6), "l");
    legM_fit()->AddEntry((*var.temp(2))(), Form("Antiproton Template %.2f K", (*var.temp(2))()->Integral()*1.0e-3), "l");
    legM_fit()->SetTextFont(43);
    legM_fit()->SetTextSize(15);
    legM_fit()->SetFillColor(0);
    legM_fit.draw();
    //editor.save();

    // Fit
    Fit::RooVar roovar2("sqrm2", hMassHC, HistList({ hMassHC_pr_pure, hMassHC_d_pure, hMassHC_ap_pure, hMassHC_kp, hMassHC_kn, hMassHC_pip, hMassHC_pin }));
    Fit::RooSysResult rlt2(roovar2, true, 100);
    Fit::RooVar var2 = rlt2.var();
    Fit::RooPar stdpar2 = rlt2.std_par();
    Fit::RooPar syspar2 = rlt2.sys_par();

    std::cerr << Form("ROOFIT 1) %14.8f ERR %14.8f %14.8f\n", stdpar2.val(0), stdpar2.err(0), syspar2.err(0));
    std::cerr << Form("ROOFIT 2) %14.8f ERR %14.8f %14.8f\n", stdpar2.val(1), stdpar2.err(1), syspar2.err(1));
    std::cerr << Form("ROOFIT 3) %14.8f ERR %14.8f %14.8f\n", stdpar2.val(2), stdpar2.err(2), syspar2.err(2));
    std::cerr << Form("ROOFIT 4) %14.8f ERR %14.8f %14.8f\n", stdpar2.val(3), stdpar2.err(3), syspar2.err(3));
    std::cerr << Form("ROOFIT 5) %14.8f ERR %14.8f %14.8f\n", stdpar2.val(4), stdpar2.err(4), syspar2.err(4));
    std::cerr << Form("ROOFIT 6) %14.8f ERR %14.8f %14.8f\n", stdpar2.val(5), stdpar2.err(5), syspar2.err(5));
    std::cerr << Form("ROOFIT 7) %14.8f ERR %14.8f %14.8f\n", stdpar2.val(6), stdpar2.err(6), syspar2.err(6));

    THStack* stackM_fit2 = Hist::Collect("hPROF_mass_fit2", HistList({ var2.temp(0), var2.temp(1), var2.temp(2), var2.temp(3), var2.temp(4), var2.temp(5), var2.temp(6) }));
    stackM_fit2->SetTitle("|Z|=1 & (0.50 < #beta < 0.80)");

    hMassHC->style(Line(kBlack, 0, 1.0),  Marker(kBlack,  MarkerStyle(MarkerShape::kCircle ), 0.75));
    
    var2.temp(0)->style(Line(kBlue, 0, 2.0),  Marker(kBlue,  MarkerStyle(MarkerShape::kCircle )), Fill(kBlue));
    var2.temp(1)->style(Line(kGreen+1, 0, 2.0),  Marker(kGreen+1,  MarkerStyle(MarkerShape::kCircle )), Fill(kGreen+1));
    var2.temp(2)->style(Line(kRed, 0, 2.0),  Marker(kRed,  MarkerStyle(MarkerShape::kCircle )), Fill(kRed));
    var2.temp(3)->style(Line(kYellow+1, 0, 2.0),  Marker(kYellow+1,  MarkerStyle(MarkerShape::kCircle )));
    var2.temp(4)->style(Line(kYellow+2, 0, 2.0),  Marker(kYellow+2,  MarkerStyle(MarkerShape::kCircle )));
    var2.temp(5)->style(Line(kMagenta+0, 0, 2.0),  Marker(kMagenta+0,  MarkerStyle(MarkerShape::kCircle )));
    var2.temp(6)->style(Line(kMagenta+1, 0, 2.0),  Marker(kMagenta+1,  MarkerStyle(MarkerShape::kCircle )));

    editor.create();
    editor.cd(0, PadAxis(0, 1));
    (*hMassHC)()->Draw("pe");
    (*hMassHC)()->SetMinimum(0.1);
    (*hMassHC)()->GetXaxis()->SetTitle("Mass/Z [GV/c^{2}]");
    (*hMassHC)()->GetYaxis()->SetTitle("Events / Bin");
    stackM_fit2->Draw("nostack hist same");
    (*hMassHC)()->Draw("pe same");
    Legend legM_fit2("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.40, 0.50, 0.85));
    //legM_fit2()->AddEntry((TObject*)nullptr, stackM_fit2->GetTitle(), "");
    //legM_fit2()->AddEntry((*hMassHC)(), Form("Data (05/24/2011 ~ 12/14/2018) %.2f B", (*hMassHC)()->Integral()*1.0e-9), "lp");
    //legM_fit2()->AddEntry((*var2.temp(0))(), Form("p^{+} Template %.2f B", (*var2.temp(0))()->Integral()*1.0e-9), "l");
    //legM_fit2()->AddEntry((*var2.temp(1))(), Form("D^{+} Template %.2f M", (*var2.temp(1))()->Integral()*1.0e-6), "l");
    //legM_fit2()->AddEntry((*var2.temp(2))(), Form("p^{-} Template %.2f K", (*var2.temp(2))()->Integral()*1.0e-3), "l");
    //legM_fit2()->AddEntry((*var2.temp(3))(), Form("K^{+} Template %.2f K", (*var2.temp(3))()->Integral()*1.0e-3), "l");
    //legM_fit2()->AddEntry((*var2.temp(4))(), Form("K^{-} Template %.2f K", (*var2.temp(4))()->Integral()*1.0e-3), "l");
    //legM_fit2()->AddEntry((*var2.temp(5))(), Form("#pi^{+} Template %.2f K", (*var2.temp(5))()->Integral()*1.0e-3), "l");
    //legM_fit2()->AddEntry((*var2.temp(6))(), Form("#pi^{-} Template %.2f K", (*var2.temp(6))()->Integral()*1.0e-3), "l");
    legM_fit2()->AddEntry((*hMassHC)(), Form("Data (05/24/2011 ~ 12/14/2018)"), "lp");
    legM_fit2()->AddEntry((*var2.temp(0))(), Form("p^{+} MC"), "l");
    legM_fit2()->AddEntry((*var2.temp(1))(), Form("D^{+} MC"), "l");
    legM_fit2()->AddEntry((*var2.temp(2))(), Form("p^{-} MC"), "l");
    legM_fit2()->AddEntry((*var2.temp(3))(), Form("K^{+} MC"), "l");
    legM_fit2()->AddEntry((*var2.temp(4))(), Form("K^{-} MC"), "l");
    legM_fit2()->AddEntry((*var2.temp(5))(), Form("#pi^{+} MC"), "l");
    legM_fit2()->AddEntry((*var2.temp(6))(), Form("#pi^{-} MC"), "l");
    legM_fit2()->SetTextFont(43);
    legM_fit2()->SetTextSize(15);
    legM_fit2()->SetFillColor(0);
    legM_fit2.draw();
    editor.save();

    editor.close();

    TFile * ofle = new TFile("out/issprof.root", "RECREATE");
    ofle->cd();

    haccp->Write();
    haccd->Write();

    (*hMassCK)()->Write();
    (*hMassKF)()->Write();
    (*hMassHC)()->Write();
    stackM->Write();

    stackM_fit->Write();
    
    (*hRvarCK)()->Write();
    (*hRvarHC)()->Write();
    stackR->Write();

    (*hMassHC_pr)()->Write();
    (*hMassHC_ap)()->Write();
    (*hMassHC_d )()->Write();
   
    (*var.temp(0))()->Write();
    (*var.temp(1))()->Write();
    (*var.temp(2))()->Write();

    (*hMassHC_cnt_pr)()->Write();
    (*hMassHC_cnt_ap)()->Write();
    (*hMassHC_cnt_d )()->Write();
    (*hMassHC_cnt_iss)()->Write();

    stackM_fit2->Write();
    (*var2.temp(0))()->Write();
    (*var2.temp(1))()->Write();
    (*var2.temp(2))()->Write();
    (*var2.temp(3))()->Write();
    (*var2.temp(4))()->Write();
    (*var2.temp(5))()->Write();
    (*var2.temp(6))()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
