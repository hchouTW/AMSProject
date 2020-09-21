#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "HistFit1D.h"
#include "HistFit1D.C"

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "45";
    
    const double acc_geom_factor = (3.9 * 3.9 * TMath::Pi()); // (~47.78 m^2 sr)

    std::string path = "/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15";

    UInt_t cntev = 0;
        
    UInt_t cntpr = 0;
    TFile* fmcpr = TFile::Open(Form("%s/mcpr%s/YiMdst.root", path.c_str(), subv.c_str()));
    TTree* tmcpr = (TTree*)fmcpr->Get("runlist");
    tmcpr->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
    double density_fact_pr = static_cast<double>(cntpr) / (std::log(1000.0) - std::log(0.2));

    TH1D* hmcpr_cnt = (TH1D*)fmcpr->Get("hFlx_cnt_MC");
    TH1D* hist_acc = new TH1D("hFlx_acc", "", hmcpr_cnt->GetXaxis()->GetNbins(), hmcpr_cnt->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= hist_acc->GetXaxis()->GetNbins(); ++ib) {
        double dlogR = std::log(hmcpr_cnt->GetXaxis()->GetBinUpEdge(ib)) - std::log(hmcpr_cnt->GetXaxis()->GetBinLowEdge(ib));
        double gcnt  = density_fact_pr * dlogR;
        double rerr  = hmcpr_cnt->GetBinError(ib) / hmcpr_cnt->GetBinContent(ib);
    
        double acc_val = acc_geom_factor * (hmcpr_cnt->GetBinContent(ib) / gcnt);
        double acc_err = acc_val * rerr;

        if (!std::isfinite(acc_val) || acc_val < 0.0) continue;
        if (!std::isfinite(acc_err) || acc_err < 0.0) continue;
        hist_acc->SetBinContent(ib, acc_val);
        hist_acc->SetBinError  (ib, acc_err);
    }
    hist_acc->GetXaxis()->SetTitle("|Rigidity| [GV/c]");
    hist_acc->GetYaxis()->SetTitle("Acceptance [m^2 sr]");

    TH1D* hist_dR = new TH1D("hFlx_dR", "", hist_acc->GetXaxis()->GetNbins(), hist_acc->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= hist_dR->GetXaxis()->GetNbins(); ++ib) {
        double dR = hist_dR->GetXaxis()->GetBinUpEdge(ib) - hist_dR->GetXaxis()->GetBinLowEdge(ib);
        if (!std::isfinite(dR) || dR <= 0.0) continue;
        hist_dR->SetBinContent(ib, dR);
        hist_dR->SetBinError  (ib, 0);
    }

    TFile* fiss = TFile::Open(Form("%s/iss%s/YiMdst.root", path.c_str(), subv.c_str()));
    TH1D* hiss_lv = (TH1D*)fiss->Get("hFlx_lv");
    TH1D* hist_lv = hiss_lv;

    TH1D*  hiss_alltrg = (TH1D*)fiss->Get("hFlxP_cnt_alltrg");
    TH1D*  hiss_phytrg = (TH1D*)fiss->Get("hFlxP_cnt");

    TH1D* hist_trg = new TH1D("hFlx_trg", "", hiss_alltrg->GetXaxis()->GetNbins(), hiss_alltrg->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= hist_trg->GetXaxis()->GetNbins(); ++ib) {
        double eff = hiss_phytrg->GetBinContent(ib) / hiss_alltrg->GetBinContent(ib);
        double err = (hiss_alltrg->GetBinError(ib) / hiss_alltrg->GetBinContent(ib)) * std::sqrt(eff * (1.0 - eff));
        
        if (!std::isfinite(eff) || eff < 0.0) continue;
        if (!std::isfinite(err) || err < 0.0) continue;
        hist_trg->SetBinContent(ib, eff);
        hist_trg->SetBinError  (ib, err);
    }
    
    TFile* ofle = new TFile("out/prflux_acc.root", "RECREATE");
    ofle->cd();

    hist_dR->Write();
    hist_acc->Write();

    hist_lv->Write();
    hist_trg->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
