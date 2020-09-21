#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "HistFit1D.h"
#include "HistFit1D.C"


int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "45";

    constexpr int NSET = 6;
    std::array<std::string, NSET> name_of_hcnt({ "hLtf_cnt_MC", "hLrh_cnt_MC", "hIin_cnt_MC", "hIex_cnt_MC", "hHex_cnt_MC", "hHfs_cnt_MC"});
    std::array<std::string, NSET> name_of_accp({ "hLtf_accp", "hLrh_accp", "hIin_accp", "hIex_accp", "hHex_accp", "hHfs_accp" });
    std::array<std::string, NSET> name_of_crss({ "hLtf_crss", "hLrh_crss", "hIin_crss", "hIex_crss", "hHex_crss", "hHfs_crss" });
    std::vector<std::array<double, 2>> range_of_accp({ { 0.5, 500.0 }, { 3.5, 500.0 }, { 8.0, 500.0 }, { 8.0, 500.0 }, { 8.0, 500.0 }, { 8.0, 500.0 } });
    std::array<TH1D*, NSET>       hist_accp_raw;
    std::array<TH1D*, NSET>       hist_crss_raw;
    std::array<TH1D*, NSET>       hist_accp;
    std::array<TH1D*, NSET>       hist_crss;

    TF1* accp_func = new TF1("accp_func", "[0] + TMath::Power((x / [1]), -[2])", 0.5, 10000);
    accp_func->SetParameters(1.0, 0.3, 0.5); 
    for (int is = 0; is < NSET; ++is) {
        UInt_t cntev = 0;
        
        UInt_t cntpr = 0;
        TFile* fmcpr = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr%s/YiMdst.root", subv.c_str()));
        TH1D*  hmcpr = (TH1D*)fmcpr->Get(name_of_hcnt[is].c_str());
        TTree* tmcpr = (TTree*)fmcpr->Get("runlist");
        tmcpr->SetBranchAddress("event", &cntev);
        for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
        
        UInt_t cntap = 0;
        TFile* fmcap = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap%s/YiMdst.root", subv.c_str()));
        TH1D*  hmcap = (TH1D*)fmcap->Get(name_of_hcnt[is].c_str());
        TTree* tmcap = (TTree*)fmcap->Get("runlist");
        tmcap->SetBranchAddress("event", &cntev);
        for (int it = 0; it < tmcap->GetEntries(); ++it) { tmcap->GetEntry(it); cntap+=cntev; }

        TH1D* haccp_raw = new TH1D(Form("%s_raw", name_of_accp[is].c_str()), "", hmcpr->GetXaxis()->GetNbins(), hmcpr->GetXaxis()->GetXbins()->GetArray());
        for (int ib = 1; ib <= haccp_raw->GetXaxis()->GetNbins(); ++ib) {
            double errpr = hmcpr->GetBinError(ib) / hmcpr->GetBinContent(ib);
            double errap = hmcap->GetBinError(ib) / hmcap->GetBinContent(ib);
            double error = std::sqrt(errpr * errpr + errap * errap);
            double accp = (hmcap->GetBinContent(ib) / static_cast<double>(cntap)) / (hmcpr->GetBinContent(ib) / static_cast<double>(cntpr));
            if (!std::isfinite(accp) || accp <= 0.0) continue;
            haccp_raw->SetBinContent(ib, 1.0/accp);
            haccp_raw->SetBinError  (ib, 1.0/accp * error);
        }
        haccp_raw->GetXaxis()->SetTitle("|Rigidity| [GV/c]");
        haccp_raw->GetYaxis()->SetTitle("Acceptance Correction");
        hist_accp_raw[is] = haccp_raw;
        
        UInt_t cntapp = 0;
        TFile* fmcapp = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap_plus%s/YiMdst.root", subv.c_str()));
        TH1D*  hmcapp = (TH1D*)fmcapp->Get(name_of_hcnt[is].c_str());
        TTree* tmcapp = (TTree*)fmcapp->Get("runlist");
        tmcapp->SetBranchAddress("event", &cntev);
        for (int it = 0; it < tmcapp->GetEntries(); ++it) { tmcapp->GetEntry(it); cntapp+=cntev; }
        
        UInt_t cntapm = 0;
        TFile* fmcapm = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap_minus%s/YiMdst.root", subv.c_str()));
        TH1D*  hmcapm = (TH1D*)fmcapm->Get(name_of_hcnt[is].c_str());
        TTree* tmcapm = (TTree*)fmcapm->Get("runlist");
        tmcapm->SetBranchAddress("event", &cntev);
        for (int it = 0; it < tmcapm->GetEntries(); ++it) { tmcapm->GetEntry(it); cntapm+=cntev; }
        
        TH1D* hcrss_raw = new TH1D(Form("%s_raw", name_of_crss[is].c_str()), "", hmcapp->GetXaxis()->GetNbins(), hmcapp->GetXaxis()->GetXbins()->GetArray());
        for (int ib = 1; ib <= hcrss_raw->GetXaxis()->GetNbins(); ++ib) {
            double errapp = hmcapp->GetBinError(ib) / hmcapp->GetBinContent(ib);
            double errapm = hmcapm->GetBinError(ib) / hmcapm->GetBinContent(ib);
            double error = std::sqrt(errapp * errapp + errapm * errapm);
            double cross = (hmcapp->GetBinContent(ib) / static_cast<double>(cntapp)) / (hmcapm->GetBinContent(ib) / static_cast<double>(cntapm));
            if (!std::isfinite(cross) || cross <= 0.0) continue;
            hcrss_raw->SetBinContent(ib, cross);
            hcrss_raw->SetBinError  (ib, cross * error);
        }
        hcrss_raw->GetXaxis()->SetTitle("|Rigidity| [GeV]");
        hcrss_raw->GetYaxis()->SetTitle("Acceptance Ratio (+10\%/-10\% Cross Section)");
        hist_crss_raw[is] = hcrss_raw;

        accp_func->SetParameters(1.0, 0.3, 0.5);
        TH1D* haccp_raw_data = (TH1D*) haccp_raw->Clone();
        haccp_raw_data->Fit(accp_func, "q0", "", range_of_accp[is][0], range_of_accp[is][1]);
        haccp_raw_data->Fit(accp_func, "q0", "", range_of_accp[is][0], range_of_accp[is][1]);
        haccp_raw_data->Fit(accp_func, "q0", "", range_of_accp[is][0], range_of_accp[is][1]);

        Axis AXrig("|Rigidity| [GV/c]", 10000, haccp_raw->GetXaxis()->GetXmin(), haccp_raw->GetXaxis()->GetXmax(), AxisScale::kLog);
        TH1D* haccp = new TH1D(name_of_accp[is].c_str(), "", AXrig.nbin(), &(AXrig(0)));
        for (int ib = 1; ib <= AXrig.nbin(); ++ib) {
            haccp->SetBinContent(ib, accp_func->Eval(AXrig.center(ib, AxisScale::kLog)));
        }
        haccp->GetXaxis()->SetTitle("|Rigidity| [GV/c]");
        haccp->GetYaxis()->SetTitle("Acceptance Correction");
        hist_accp[is] = haccp;
        
        TH1D* hcrss = new TH1D(name_of_crss[is].c_str(), "", hmcapp->GetXaxis()->GetNbins(), hmcapp->GetXaxis()->GetXbins()->GetArray());
        for (int ib = 1; ib <= hcrss->GetXaxis()->GetNbins(); ++ib) {
            hcrss->SetBinContent(ib, hcrss_raw->GetBinContent(ib) - 1.0);
            hcrss->SetBinError  (ib, hcrss_raw->GetBinError  (ib));
        }
        hcrss->GetXaxis()->SetTitle("|Rigidity| [GeV]");
        hcrss->GetYaxis()->SetTitle("Acceptance Difference (+10\%/-10\% Cross Section)");
        hist_crss[is] = hcrss;
    }

    TFile * ofle = new TFile("out/apflux_acc.root", "RECREATE");
    ofle->cd();

    for (int is = 0; is < NSET; ++is) hist_accp_raw[is]->Write();
    for (int is = 0; is < NSET; ++is) hist_crss_raw[is]->Write();
    for (int is = 0; is < NSET; ++is) hist_accp[is]->Write();
    for (int is = 0; is < NSET; ++is) hist_crss[is]->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
