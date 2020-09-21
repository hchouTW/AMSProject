#include <CPPLibs.h>
#include <ROOTLibs.h>
#include "Math/SpecFunc.h"

#include "HistFit1D.h"
#include "HistFit1D.C"

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

/*
std::array<long double, 2> PoissonProb(long double sig, long double bkg = 0.0, unsigned long obs = 0) {
    double normalized = 1.0;
    if (bkg > 0.0 && obs != 0) {
        long double den = 0.0;
        std::vector<long double> factorial_x;
        std::vector<long double> factorial_b;
        for (unsigned long iobs = 0; iobs <= obs; ++iobs) {
            if (iobs == 0) {
                factorial_x.push_back(1.0); 
                factorial_b.push_back(1.0); 
            }
            else {
                factorial_x.push_back(static_cast<long double>(iobs) * factorial_x.back());
                factorial_b.push_back(bkg * factorial_b.back());
            }
            den += (factorial_b.back() / factorial_x.back());
        }
        den *= factorial_x.back();
        if (std::isfinite(den)) normalized = (1.0 / den);
    }
    else if (bkg == 0.0 && obs != 0) {
        long double den = 1.0;
        for (unsigned long iobs = 1; iobs <= obs; ++iobs) {
            den *= static_cast<long double>(iobs);
        }
        if (std::isfinite(den)) normalized = (1.0 / den);
    }


    long double pdf = 0.0;
    if      (obs ==   0) pdf = std::exp(-sig);
    else if (bkg == 0.0) pdf = normalized * std::pow(sig, static_cast<long double>(obs)) * std::exp(-sig);
    else                 pdf = normalized * std::pow(sig + bkg, static_cast<long double>(obs)) * std::exp(-sig);
    if (!std::isfinite(pdf)) pdf = 0.0;
    
    long double cdf = 0.0;
    if      (obs ==   0) cdf = 1.0 - std::exp(-sig);
    else if (bkg == 0.0) cdf = ROOT::Math::inc_gamma(obs + 1, sig);
    else                 cdf = ROOT::Math::inc_gamma(obs + 1, sig + bkg);
    if (!std::isfinite(cdf)) cdf = 0.0;
    
    return std::array<long double, 2>({ pdf, cdf });
}
*/
long double PoissonCDF(long double sig, long double bkg = 0.0, unsigned long obs = 0) {
    long double cdf = 0.0;
    if      (obs ==   0) cdf = 1.0 - std::exp(-sig);
    else if (bkg == 0.0) cdf = ROOT::Math::inc_gamma(obs + 1, sig);
    else                 cdf = ROOT::Math::inc_gamma(obs + 1, sig + bkg);
    if (!std::isfinite(cdf)) cdf = 0.0;
    return cdf;
}

double FindPoissonCL(long double bkg = 0.0, unsigned long obs = 0, long double alpha = 0.95) {
    std::vector<double> seedx;
    for (int i = 0; i <= 10000; ++i)
        seedx.push_back(0.01*static_cast<double>(i));

    std::vector<double> seedv(seedx.size(), 0.0);
    for (int i = 0; i < seedx.size(); ++i) {
        seedv.at(i) = PoissonCDF(seedx.at(i), bkg, obs);
    }

    int initloc = 0;
    for (int i = 0; i < seedv.size() - 1; ++i) {
        if (alpha > seedv.at(i) && alpha < seedv.at(i+1)) { initloc = i; break; }
        if (i == 0 && alpha < seedv.at(i)) return 0.0;
        if (i == (seedv.size() - 2)) return -1.0;
    }
    std::array<double, 2> bdx({ seedx.at(initloc), seedx.at(initloc + 1) });
    std::array<double, 2> bdv({ seedv.at(initloc), seedv.at(initloc + 1) });
    double sig = seedx.at(initloc);

    bool status = false;
    const double tolerance = 1.0e-6;
    do {
        double wgtx = (alpha - bdv[0]) / (bdv[1] - bdv[0]);
        double newx = (1.0 - wgtx) * bdx[0] + wgtx * bdx[1];
        auto&& rlt = PoissonCDF(newx, bkg, obs);
        if (alpha < rlt) {
            bdx[1] = newx;
            bdv[1] = rlt;
        }
        else {
            bdx[0] = newx;
            bdv[0] = rlt;
        }

        sig = newx;
        status = (std::abs(rlt - alpha) < tolerance);
    } while(!status);

    return sig;
}


//static constexpr int BinMlw = 92;
//static constexpr int BinMup = 108;

static constexpr int NBinM = 240;
static constexpr int BinMlw = 1;
static constexpr int BinMup = 240;

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "66";
    
    TFile* file_acc = TFile::Open("out/adflux_acc.root");
    double TFAccFact = 1.0 / ((TH1D*)file_acc->Get("hTFcorr"))->Integral();
    double RHAccFact = 1.0 / ((TH1D*)file_acc->Get("hRHcorr"))->Integral();
    
    TFile* file_iss = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()));
    TTree* varsTF = (TTree*) file_iss->Get("varsTF");
    TTree* varsRH = (TTree*) file_iss->Get("varsRH");

    std::cerr << Form("\n============ TOF  ===============\n");
    varsTF->Scan("run:evt:sqrt(sqrm):rig:bta:llr:geom_lx:geom_ly:mutr_lx:mutr_ly:mutr_lb:a_phys_lx:a_phys_ly:a_phys_lb:b_phys_lx:b_phys_ly:b_phys_lb", "evt==133569||evt==67582||evt==409381||evt==574622||evt==808247||evt==527467||evt==391367", "col=10d:d:::::::::");
    std::cerr << Form("\n============ RICH ===============\n");
    varsRH->Scan("run:evt:sqrt(sqrm):rig:bta", "", "col=10d:d:::");
    std::cerr << Form("\n=================================\n");
    
    TH1D* hTFmutr = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hTF_mutr_lxyb"));
    TH1D* hTFphys = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hTF_phys_lxyb"));

    TH1D* hTFmassISS_PRE  = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hTF_mass_CF_PRE"));
    TH1D* hTFmassISS_MUTR = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hTF_mass_CF_SEL_MUTR"));
    TH1D* hTFmassISS_PHYS = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hTF_mass_CF_SEL_PHYS"));
    TH1D* hTFmassISS      = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hTF_mass_CF"));
    TH1D* hTFmassDE = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/mcd%s/YiMdst.root"  , subv.c_str()))->Get("hTF_mass_de"));
    TH1D* hTFmassPR = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/mcprL%s/YiMdst.root", subv.c_str()))->Get("hTF_mass_pr"));
    TH1D* hTFmassPP = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/mcprL%s/YiMdst.root", subv.c_str()))->Get("hTF_mass_pp"));
    TH1D* hTFmassAP = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/mcap%s/YiMdst.root" , subv.c_str()))->Get("hTF_mass_ap"));
    
    TH1D* hRHmutr0 = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hRH_tf_mutr_lxyb"));
    TH1D* hRHmutr1 = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hRH_rh_mutr_lxyb"));
    TH1D* hRHphys0 = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hRH_tf_phys_lxyb"));
    TH1D* hRHphys1 = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hRH_rh_phys_lxyb"));
    
    TH1D* hRHmassISS_PRE  = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hRH_mass_CF_PRE"));
    TH1D* hRHmassISS_MUTR = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hRH_mass_CF_SEL_MUTR"));
    TH1D* hRHmassISS_PHYS = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hRH_mass_CF_SEL_PHYS"));
    TH1D* hRHmassISS      = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/iss%s/YiMdst.root", subv.c_str()))->Get("hRH_mass_CF"));
    TH1D* hRHmassDE = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/mcd%s/YiMdst.root"  , subv.c_str()))->Get("hRH_mass_de"));
    TH1D* hRHmassPR = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/mcprL%s/YiMdst.root", subv.c_str()))->Get("hRH_mass_pr"));
    TH1D* hRHmassPP = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/mcprL%s/YiMdst.root", subv.c_str()))->Get("hRH_mass_pp"));
    TH1D* hRHmassAP = (TH1D*) (TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/adflux/20Jan15/mcap%s/YiMdst.root" , subv.c_str()))->Get("hRH_mass_ap"));

    PdfEditor editor(Window(WindowSize::kWideSliceLR), "adflux_flx", "out");

    // TOF
    HistFit::Axis1D AX1D_TFBKG(
        "Mass/Z [(GV/c^{2})]",
        "Events/Bin",
        hTFmassISS->GetXaxis(), 1, NBinM/2);
    
    TH1D*              hTFBKGlink_smp  = hTFmassISS;
    std::vector<TH1D*> hTFBKGlink_tmps({ hTFmassPP });
    
    std::string TFBKGprefix = "hTFBKG_";
    HistFit::HistFit1D TFBKGfit1D(hTFBKGlink_smp, hTFBKGlink_tmps, AX1D_TFBKG, TFBKGprefix, false, true);
    Hist* hTFpp = nullptr;

    HistFit::Axis1D AX1D_TF(
        "Mass/Z [(GV/c^{2})]",
        "Events/Bin",
        hTFmassISS->GetXaxis(), 1, NBinM);

    TH1D*              hTFlink_smp  = hTFmassISS;
    std::vector<TH1D*> hTFlink_tmps({ hTFmassAP, hTFmassDE, hTFmassPR });

    std::string TFprefix = "hTF_";
    HistFit::HistFit1D TFfit1D(hTFlink_smp, hTFlink_tmps, AX1D_TF, TFprefix, false, true);
    Hist* hTFsmp = nullptr; 
    Hist* hTFsum = nullptr; 
    Hist* hTFap  = nullptr; 
    Hist* hTFde  = nullptr; 
    Hist* hTFpr  = nullptr;
    THStack* hTFfit = nullptr;

    unsigned long TFobserved = 0;
    double TFlambda_sig = 0.0;
    double TFlambda_bkg = 0.0;
    double TFlambda_bkg_ap = 0.0;
    double TFlambda_bkg_pr = 0.0;
    double TFlambda_bkg_de = 0.0;
    double TFlambda_bkg_pp = 0.0;
    double TFlambda_sig_a90 = 0.0;
    double TFlambda_sig_a95 = 0.0;
    double TFlambda_sig_a99 = 0.0;

    Hist* hTFcanvas = nullptr;
    if (TFfit1D.status() && TFBKGfit1D.status()) {
        std::string prefix = TFprefix;
        HistFit::HistFit1D& fit1D = TFfit1D;
        HistFit::HistFit1D& BKGfit1D = TFBKGfit1D;

        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dap  = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dde  = fit1D.wgt_tmps(1);
        const HistFit::Hist1D& h1Dpr  = fit1D.wgt_tmps(2);
        const HistFit::Hist1D& h1Dpp  = BKGfit1D.wgt_tmps(0);
        
        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hap  = Hist::New(h1Dap.get());
        Hist* hde  = Hist::New(h1Dde.get());
        Hist* hpr  = Hist::New(h1Dpr.get());
        Hist* hpp  = Hist::New(h1Dpp.get());

		hsmp->style(Line(kBlack   , 0, 2), Marker(kBlack   , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kYellow+1, 0, 2), Marker(kYellow+1, MarkerStyle(MarkerShape::kCircle )));
		hap ->style(Line(kRed     , 0, 2), Marker(kRed     , MarkerStyle(MarkerShape::kCircle )));
		hde ->style(Line(kBlue    , 0, 2), Marker(kBlue    , MarkerStyle(MarkerShape::kCircle )));
		hpr ->style(Line(kGreen+2 , 0, 2), Marker(kGreen+2 , MarkerStyle(MarkerShape::kCircle )));
		hpp ->style(Line(kGreen+2 , 0, 2), Marker(kGreen+2 , MarkerStyle(MarkerShape::kCircle )));
       
        TFobserved   = (*hsmp)()->Integral(21, 56);
        TFlambda_bkg = (*hap)()->Integral(21, 56) + (*hpr)()->Integral(21, 56) + (*hde)()->Integral(21, 56);
        TFlambda_bkg_ap = (*hap)()->Integral();
        TFlambda_bkg_pr = (*hpr)()->Integral();
        TFlambda_bkg_de = (*hde)()->Integral();
        TFlambda_bkg_pp = (*hpp)()->Integral(21, 56);
        TFlambda_sig_a90 = FindPoissonCL(TFlambda_bkg, TFobserved, 0.90);
        TFlambda_sig_a95 = FindPoissonCL(TFlambda_bkg, TFobserved, 0.95);
        TFlambda_sig_a99 = FindPoissonCL(TFlambda_bkg, TFobserved, 0.99);

        THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hpr, hde, hap }));
        hTFfit = hfit;
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.55, 0.90));
        leg_table()->AddEntry((*hsmp)(), "Data"      , "lp");
        leg_table()->AddEntry((*hsum)(), "Sum"       , "l");
        leg_table()->AddEntry((*hap)() , "Antiproton", "l");
        leg_table()->AddEntry((*hde)() , "Deuterium" , "l");
        leg_table()->AddEntry((*hpr)() , "Proton"    , "l");
        leg_table()->AddEntry((TObject*)0, Form("#hat{N}=%ld", TFobserved), "");
        leg_table()->AddEntry((TObject*)0, Form("#lambda_{b}=%.2f", TFlambda_bkg), "");
        leg_table()->AddEntry((TObject*)0, Form("#lambda_{s,#alpha=0.95}=%.2f", TFlambda_sig_a95), "");
        leg_table()->AddEntry((TObject*)0, Form("#lambda_{s,#alpha=0.99}=%.2f", TFlambda_sig_a99), "");
        leg_table()->SetFillColor(0);

        Hist* hcanvas = Hist::New(
            Form("%scanvas", prefix.c_str()), 
            HistAxis(hsum->xaxis(), Axis("", 10000, 0.1, 5.0 * (*hsmp)()->GetMaximum(), AxisScale::kLog)));
        (*hcanvas)()->GetXaxis()->SetTitle(AX1D_TF.name_x().c_str());
        (*hcanvas)()->GetYaxis()->SetTitle(AX1D_TF.name_y().c_str());
        
        editor.create();
        editor.cd(0, PadAxis(0, 1));
        (*hcanvas)()->GetXaxis()->SetRange(1, NBinM/2);
        (*hcanvas)()->GetYaxis()->SetRangeUser((*hcanvas)()->GetYaxis()->GetXmin(), (*hcanvas)()->GetYaxis()->GetXmax() * TFlambda_bkg_ap/TFlambda_bkg_pr);
        (*hcanvas)()->Draw();
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        leg_table.draw();
        editor.save();

        editor.create();
        editor.cd(0, PadAxis(0, 1));
        (*hcanvas)()->GetXaxis()->UnZoom();
        (*hcanvas)()->GetYaxis()->UnZoom();
        (*hcanvas)()->Draw();
        hfit->Draw("nostack histi same");
        (*hsmp)()->Draw("pe same");
        leg_table.draw();
        editor.save();
        
        hTFsmp = hsmp;
        hTFsum = hsum;
        hTFap  = hap ;
        hTFde  = hde ;
        hTFpr  = hpr ;
        hTFpp  = hpp ;
        hTFcanvas = hcanvas;
    }
    TFlambda_sig = TFlambda_sig_a95;

    std::cerr << Form("\n\n");
    std::cerr << Form("<< TF >>   Num Observed Events : %ld\n"   , TFobserved);
    std::cerr << Form("<< TF >>        Antiproton Num : %14.2f\n", TFlambda_bkg_ap);
    std::cerr << Form("<< TF >>            Proton Num : %14.2f\n", TFlambda_bkg_pr);
    std::cerr << Form("<< TF >>          Deuteron Num : %14.2f\n", TFlambda_bkg_de);
    std::cerr << Form("<< TF >>         InvProton Num : %14.2f\n", TFlambda_bkg_pp);
    std::cerr << Form("<< TF >>  Background Parameter : %14.8f\n", TFlambda_bkg);
    std::cerr << Form("<< TF >>      Signal Parameter : %14.8f\n", TFlambda_sig);
    std::cerr << Form("<< TF >>           Flux Factor : %E\n"    , TFAccFact);
    std::cerr << Form("<< TF >> Signal Parameter (0.90) : %14.8f    FLUX(%E)\n", TFlambda_sig_a90, TFAccFact * TFlambda_sig_a90);
    std::cerr << Form("<< TF >> Signal Parameter (0.95) : %14.8f    FLUX(%E)\n", TFlambda_sig_a95, TFAccFact * TFlambda_sig_a95);
    std::cerr << Form("<< TF >> Signal Parameter (0.99) : %14.8f    FLUX(%E)\n", TFlambda_sig_a99, TFAccFact * TFlambda_sig_a99);
    std::cerr << Form("\n\n");
    
    // RICH
    HistFit::Axis1D AX1D_RHBKG(
        "Mass/Z [(GV/c^{2})]",
        "Events/Bin",
        hRHmassISS->GetXaxis(), 1, NBinM/2);
    
    TH1D*              hRHBKGlink_smp  = hRHmassISS;
    std::vector<TH1D*> hRHBKGlink_tmps({ hRHmassPP });
    
    std::string RHBKGprefix = "hRHBKG_";
    HistFit::HistFit1D RHBKGfit1D(hRHBKGlink_smp, hRHBKGlink_tmps, AX1D_RHBKG, RHBKGprefix, false, true);
    Hist* hRHpp = nullptr;

    HistFit::Axis1D AX1D_RH(
        "Mass/Z [(GV/c^{2})]",
        "Events/Bin",
        hRHmassISS->GetXaxis(), 1, NBinM);

    TH1D*              hRHlink_smp  = hRHmassISS;
    std::vector<TH1D*> hRHlink_tmps({ hRHmassAP, hRHmassDE, hRHmassPR });

    std::string RHprefix = "hRH_";
    HistFit::HistFit1D RHfit1D(hRHlink_smp, hRHlink_tmps, AX1D_RH, RHprefix, false, true);
    Hist* hRHsmp = nullptr; 
    Hist* hRHsum = nullptr; 
    Hist* hRHap  = nullptr; 
    Hist* hRHde  = nullptr; 
    Hist* hRHpr  = nullptr;
    THStack* hRHfit = nullptr;

    unsigned long RHobserved = 0;
    double RHlambda_sig = 0.0;
    double RHlambda_bkg = 0.0;
    double RHlambda_bkg_ap = 0.0;
    double RHlambda_bkg_pr = 0.0;
    double RHlambda_bkg_de = 0.0;
    double RHlambda_bkg_pp = 0.0;
    double RHlambda_sig_a90 = 0.0;
    double RHlambda_sig_a95 = 0.0;
    double RHlambda_sig_a99 = 0.0;
        
    Hist* hRHcanvas = nullptr;
    if (RHfit1D.status() && RHBKGfit1D.status()) {
        std::string prefix = RHprefix;
        HistFit::HistFit1D& fit1D = RHfit1D;
        HistFit::HistFit1D& BKGfit1D = RHBKGfit1D;

        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dap  = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dde  = fit1D.wgt_tmps(1);
        const HistFit::Hist1D& h1Dpr  = fit1D.wgt_tmps(2);
        const HistFit::Hist1D& h1Dpp  = BKGfit1D.wgt_tmps(0);
        
        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hap  = Hist::New(h1Dap.get());
        Hist* hde  = Hist::New(h1Dde.get());
        Hist* hpr  = Hist::New(h1Dpr.get());
        Hist* hpp  = Hist::New(h1Dpp.get());

		hsmp->style(Line(kBlack   , 0, 2), Marker(kBlack   , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kYellow+1, 0, 2), Marker(kYellow+1, MarkerStyle(MarkerShape::kCircle )));
		hap ->style(Line(kRed     , 0, 2), Marker(kRed     , MarkerStyle(MarkerShape::kCircle )));
		hde ->style(Line(kBlue    , 0, 2), Marker(kBlue    , MarkerStyle(MarkerShape::kCircle )));
		hpr ->style(Line(kGreen+2 , 0, 2), Marker(kGreen+2 , MarkerStyle(MarkerShape::kCircle )));
		hpp ->style(Line(kGreen+2 , 0, 2), Marker(kGreen+2 , MarkerStyle(MarkerShape::kCircle )));
        
        RHobserved   = (*hsmp)()->Integral(21, 56);
        RHlambda_bkg = (*hap)()->Integral(21, 56) + (*hpr)()->Integral(21, 56) + (*hde)()->Integral(21, 56);
        RHlambda_bkg_ap = (*hap)()->Integral();
        RHlambda_bkg_pr = (*hpr)()->Integral();
        RHlambda_bkg_de = (*hde)()->Integral();
        //RHlambda_bkg_pp = (*hpp)()->Integral(45, 54);
        RHlambda_bkg_pp = (*hpp)()->Integral(21, 56);
        RHlambda_sig_a90 = FindPoissonCL(RHlambda_bkg, RHobserved, 0.90);
        RHlambda_sig_a95 = FindPoissonCL(RHlambda_bkg, RHobserved, 0.95);
        RHlambda_sig_a99 = FindPoissonCL(RHlambda_bkg, RHobserved, 0.99);

        THStack* hfit = Hist::Collect(Form("%sFIT", prefix.c_str()), HistList({ hsum, hpr, hde, hap }));
        hRHfit = hfit;
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.55, 0.90));
        leg_table()->AddEntry((*hsmp)(), "Data"      , "lp");
        leg_table()->AddEntry((*hsum)(), "Sum"       , "l");
        leg_table()->AddEntry((*hap)() , "Antiproton", "l");
        leg_table()->AddEntry((*hde)() , "Deuterium" , "l");
        leg_table()->AddEntry((*hpr)() , "Proton"    , "l");
        leg_table()->AddEntry((TObject*)0, Form("#hat{N}=%ld", RHobserved), "");
        leg_table()->AddEntry((TObject*)0, Form("#lambda_{b}=%.2f", RHlambda_bkg), "");
        leg_table()->AddEntry((TObject*)0, Form("#lambda_{s,#alpha=0.95}=%.2f", RHlambda_sig_a95), "");
        leg_table()->AddEntry((TObject*)0, Form("#lambda_{s,#alpha=0.99}=%.2f", RHlambda_sig_a99), "");
        leg_table()->SetFillColor(0);

        Hist* hcanvas = Hist::New(
            Form("%scanvas", prefix.c_str()), 
            HistAxis(hsum->xaxis(), Axis("", 10000, 0.1, 5.0 * (*hsmp)()->GetMaximum(), AxisScale::kLog)));
        (*hcanvas)()->GetXaxis()->SetTitle(AX1D_RH.name_x().c_str());
        (*hcanvas)()->GetYaxis()->SetTitle(AX1D_RH.name_y().c_str());

        editor.create();
        editor.cd(0, PadAxis(0, 1));  
        (*hcanvas)()->GetXaxis()->SetRange(1, NBinM/2);
        (*hcanvas)()->GetYaxis()->SetRangeUser((*hcanvas)()->GetYaxis()->GetXmin(), (*hcanvas)()->GetYaxis()->GetXmax() * RHlambda_bkg_ap/RHlambda_bkg_pr);
        (*hcanvas)()->Draw();
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        leg_table()->Draw();
        editor.save();
        
        editor.create();
        editor.cd(0, PadAxis(0, 1));   
        (*hcanvas)()->GetXaxis()->UnZoom();
        (*hcanvas)()->GetYaxis()->UnZoom();
        (*hcanvas)()->Draw();
        hfit->Draw("nostack hist same");
        (*hsmp)()->Draw("pe same");
        leg_table()->Draw();
        editor.save();
        
        hRHsmp = hsmp;
        hRHsum = hsum;
        hRHap  = hap ;
        hRHde  = hde ;
        hRHpr  = hpr ;
        hRHpp  = hpp ;
        hRHcanvas = hcanvas;
    }
    RHlambda_sig = RHlambda_sig_a95;
    
    std::cerr << Form("\n\n");
    std::cerr << Form("<< RH >>   Num Observed Events : %ld\n"   , RHobserved);
    std::cerr << Form("<< RH >>        Antiproton Num : %14.2f\n", RHlambda_bkg_ap);
    std::cerr << Form("<< RH >>            Proton Num : %14.2f\n", RHlambda_bkg_pr);
    std::cerr << Form("<< RH >>          Deuteron Num : %14.2f\n", RHlambda_bkg_de);
    std::cerr << Form("<< RH >>         InvProton Num : %14.2f\n", RHlambda_bkg_pp);
    std::cerr << Form("<< RH >>  Background Parameter : %14.8f\n", RHlambda_bkg);
    std::cerr << Form("<< RH >>      Signal Parameter : %14.8f\n", RHlambda_sig);
    std::cerr << Form("<< RH >>           Flux Factor : %E\n"    , RHAccFact);
    std::cerr << Form("<< RH >> Signal Parameter (0.90) : %14.8f    FLUX(%E)\n", RHlambda_sig_a90, RHAccFact * RHlambda_sig_a90);
    std::cerr << Form("<< RH >> Signal Parameter (0.95) : %14.8f    FLUX(%E)\n", RHlambda_sig_a95, RHAccFact * RHlambda_sig_a95);
    std::cerr << Form("<< RH >> Signal Parameter (0.99) : %14.8f    FLUX(%E)\n", RHlambda_sig_a99, RHAccFact * RHlambda_sig_a99);
    std::cerr << Form("\n\n");
  

    // Set
    TF1* TFline = new TF1("TFline", "[0]", -5, 5);
    TFline->SetLineColor(kYellow+1);
    TFline->SetLineStyle(2);
    TFline->SetParameter(0, 0.0);

    hTFmassISS_PRE ->SetLineWidth(2);
    hTFmassISS_MUTR->SetLineWidth(2);
    hTFmassISS_PHYS->SetLineWidth(2);

    hTFmassISS_PRE ->SetLineColor(kBlack);
    hTFmassISS_MUTR->SetLineColor(kBlue);
    hTFmassISS_PHYS->SetLineColor(kRed);

    TF1* RHline = new TF1("RHline", "[0]", -5, 5);
    RHline->SetLineColor(kYellow+1);
    RHline->SetLineStyle(2);
    RHline->SetParameter(0, 0.0);
    
    hRHmassISS_PRE ->SetLineWidth(2);
    hRHmassISS_MUTR->SetLineWidth(2);
    hRHmassISS_PHYS->SetLineWidth(2);

    hRHmassISS_PRE ->SetLineColor(kBlack);
    hRHmassISS_MUTR->SetLineColor(kBlue);
    hRHmassISS_PHYS->SetLineColor(kRed);

    // Low energy region
    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    hTFmutr->Draw("colz");
    TFline->SetParameter(0, 1.75);
    TFline->Draw("l same");
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    hTFcanvas->draw();
    hTFmassISS_PRE->Draw("hist same");
    editor.save();
    
    double tf_sel_efft_mutr = 100.0 * hTFmassISS_MUTR->Integral(BinMlw, BinMup) / hTFmassISS_PRE->Integral(BinMlw, BinMup);
    Legend leg_tf_sel_mutr("", TextStyle(kBlack, 30, 43), PadWindow(0.15, 0.42, 0.70, 0.90));
    leg_tf_sel_mutr()->AddEntry(hTFmassISS_PRE , "After Selections" , "l");
    leg_tf_sel_mutr()->AddEntry(hTFmassISS_MUTR, Form("Free-Mass (#varepsilon = %5.1f%)", tf_sel_efft_mutr), "l");
    leg_tf_sel_mutr()->SetFillColor(0);
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    hTFcanvas->draw();
    hTFmassISS_PRE->Draw("hist same");
    hTFmassISS_MUTR->Draw("hist same");
    leg_tf_sel_mutr.draw();
    editor.save();
    
    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    hTFphys->Draw("colz");
    TFline->SetParameter(0, 1.75);
    TFline->Draw("l same");
    editor.save();
    
    double tf_sel_efft_phys = 100.0 * hTFmassISS_PHYS->Integral(BinMlw, BinMup) / hTFmassISS_MUTR->Integral(BinMlw, BinMup);
    Legend leg_tf_sel_phys("", TextStyle(kBlack, 30, 43), PadWindow(0.15, 0.42, 0.70, 0.90));
    leg_tf_sel_phys()->AddEntry(hTFmassISS_MUTR, "Free-Mass", "l");
    leg_tf_sel_phys()->AddEntry(hTFmassISS_PHYS, Form("Fixed-Mass (#varepsilon = %5.1f%)", tf_sel_efft_phys), "l");
    leg_tf_sel_phys()->SetFillColor(0);
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    hTFcanvas->draw();
    hTFmassISS_MUTR->Draw("hist same");
    hTFmassISS_PHYS->Draw("hist same");
    leg_tf_sel_phys.draw();
    editor.save();

    // High energy region
    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    hRHmutr0->Draw("colz");
    TFline->SetParameter(0, 2.00);
    TFline->Draw("l same");
    editor.save();

    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    hRHmutr1->Draw("colz");
    RHline->SetParameter(0, 1.75);
    RHline->Draw("l same");
    editor.save();
    
    Legend leg_rh_sel_mutr0("", TextStyle(kBlack, 20, 43), PadWindow(0.12, 0.42, 0.85, 0.90));
    leg_rh_sel_mutr0()->AddEntry(hRHmutr0, "Tracker-TOF", "");
    leg_rh_sel_mutr0()->SetFillColor(0);
    
    Legend leg_rh_sel_mutr1("", TextStyle(kBlack, 20, 43), PadWindow(0.12, 0.42, 0.85, 0.90));
    leg_rh_sel_mutr1()->AddEntry(hRHmutr1, "Tracker-TOF-RICH", "");
    leg_rh_sel_mutr1()->SetFillColor(0);
    
    editor.create("", 2, 1, PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(1, PadAxis(0, 0));
    hRHmutr0->Draw("colz");
    TFline->SetParameter(0, 2.00);
    TFline->Draw("l same");
    leg_rh_sel_mutr0.draw();
    editor.cd(2, PadAxis(0, 0));
    hRHmutr1->Draw("colz");
    RHline->SetParameter(0, 1.75);
    RHline->Draw("l same");
    leg_rh_sel_mutr1.draw();
    editor.save();
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    hRHcanvas->draw();
    hRHmassISS_PRE->Draw("hist same");
    editor.save();
    
    double rh_sel_efft_mutr = 100.0 * hRHmassISS_MUTR->Integral(BinMlw, BinMup) / hRHmassISS_PRE->Integral(BinMlw, BinMup);
    Legend leg_rh_sel_mutr("", TextStyle(kBlack, 30, 43), PadWindow(0.15, 0.42, 0.70, 0.90));
    leg_rh_sel_mutr()->AddEntry(hRHmassISS_PRE , "After Selections" , "l");
    leg_rh_sel_mutr()->AddEntry(hRHmassISS_MUTR, Form("Free-Mass (#varepsilon = %5.1f%)", tf_sel_efft_mutr), "l");
    leg_rh_sel_mutr()->SetFillColor(0);
    
    editor.create();
    editor.cd(0, PadAxis(0, 1));
    hRHcanvas->draw();
    hRHmassISS_PRE->Draw("hist same");
    hRHmassISS_MUTR->Draw("hist same");
    leg_rh_sel_mutr.draw();
    editor.save();
    
    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(1, PadAxis(0, 0));
    hRHphys0->Draw("colz");
    TFline->SetParameter(0, 2.00);
    TFline->Draw("l same");
    editor.save();
    
    editor.create("", PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(0, PadAxis(0, 0));
    hRHphys1->Draw("colz");
    RHline->SetParameter(0, 1.75);
    RHline->Draw("l same");
    editor.save();
    
    Legend leg_rh_sel_phys0("", TextStyle(kBlack, 20, 43), PadWindow(0.12, 0.42, 0.85, 0.9));
    leg_rh_sel_phys0()->AddEntry(hRHphys0, "Tracker-TOF", "");
    leg_rh_sel_phys0()->SetFillColor(0);
    
    Legend leg_rh_sel_phys1("", TextStyle(kBlack, 20, 43), PadWindow(0.12, 0.42, 0.85, 0.9));
    leg_rh_sel_phys1()->AddEntry(hRHphys1, "Tracker-TOF-RICH", "");
    leg_rh_sel_phys1()->SetFillColor(0);
    
    editor.create("", 2, 1, PadMargin(0.06, 0.1, 0.1, 0.12));
    editor.cd(1, PadAxis(0, 0));
    hRHphys0->Draw("colz");
    TFline->SetParameter(0, 2.00);
    TFline->Draw("l same");
    leg_rh_sel_phys0.draw();
    editor.cd(2, PadAxis(0, 0));
    hRHphys1->Draw("colz");
    RHline->SetParameter(0, 1.75);
    RHline->Draw("l same");
    leg_rh_sel_phys1.draw();
    editor.save();
    
    double rh_sel_efft_phys = 100.0 * hRHmassISS_PHYS->Integral(BinMlw, BinMup) / hRHmassISS_MUTR->Integral(BinMlw, BinMup);
    Legend leg_rh_sel_phys("", TextStyle(kBlack, 30, 43), PadWindow(0.15, 0.42, 0.70, 0.90));
    leg_rh_sel_phys()->AddEntry(hRHmassISS_MUTR, "Free-Mass", "l");
    leg_rh_sel_phys()->AddEntry(hRHmassISS_PHYS, Form("Fixed-Mass (#varepsilon = %5.1f%)", tf_sel_efft_phys), "l");
    leg_rh_sel_phys()->SetFillColor(0);

    editor.create();
    editor.cd(0, PadAxis(0, 1));
    hRHcanvas->draw();
    hRHmassISS_MUTR->Draw("hist same");
    hRHmassISS_PHYS->Draw("hist same");
    leg_rh_sel_phys.draw();
    editor.save();

    editor.close();

    TFile * ofle = new TFile("out/adflux_flx.root", "RECREATE");
    ofle->cd();
  
    (*hTFsmp)()->Write();
    (*hTFsum)()->Write();
    (*hTFap )()->Write();
    (*hTFde )()->Write();
    (*hTFpr )()->Write();
    (*hTFpp )()->Write();
    hTFfit->Write();

    (*hRHsmp)()->Write();
    (*hRHsum)()->Write();
    (*hRHap )()->Write();
    (*hRHde )()->Write();
    (*hRHpr )()->Write();
    (*hRHpp )()->Write();
    hRHfit->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
