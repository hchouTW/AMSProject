void draw() {
    TGraphAsymmErrors* hpbarp = (TGraphAsymmErrors*)(TFile::Open("20160406MIT.root")->Get("gpbarp"));
    TGraphAsymmErrors* hap  = (TGraphAsymmErrors*)(TFile::Open("ams02_ap_flux.root")->Get("gr_exp2"));
    TGraphAsymmErrors* hpr  = (TGraphAsymmErrors*)(TFile::Open("ams02_pr.root")->Get("gr_exp2"));
    
    TH1D* happ = (TH1D*)(TFile::Open("antip_all.root")->Get("hCrrR"));

    TH1D* hacc = (TH1D*)(TFile::Open("20160406MIT.root")->Get("hsyst_acc"));
    TH1D* hcount = (TH1D*)(TFile::Open("20160406MIT.root")->Get("hsyst_counts"));
    TH1D* hrigscale = (TH1D*)(TFile::Open("20160406MIT.root")->Get("hsyst_rigscale"));


    TH1D* hap_rat = new TH1D("hap_rat", "hap_rat", happ->GetXaxis()->GetNbins(), happ->GetXaxis()->GetXbins()->GetArray());
    TH1D* hap_flx = new TH1D("hap_flx", "hap_flx", happ->GetXaxis()->GetNbins(), happ->GetXaxis()->GetXbins()->GetArray());
    
    for (int ib = 2; ib <= happ->GetXaxis()->GetNbins()-4; ++ib) {
        double rig = happ->GetBinCenter(ib);
        double ap_rat = happ->GetBinContent(ib);
        double ap_flx = happ->GetBinContent(ib) * hpr->Eval(rig);

        double stat = (happ->GetBinError(ib) / happ->GetBinContent(ib)) * 100.;
        double acc  = hacc->Interpolate(rig);
        double cnt  = (rig < 8) ? 0 : hcount->Interpolate(rig);
        double rscl = hrigscale->Interpolate(rig);
        double err  = std::sqrt(stat * stat + acc * acc + cnt * cnt + rscl * rscl) * 0.01;
        //double err  = std::sqrt(stat * stat + cnt * cnt + rscl * rscl) * 0.01;
        std::cerr << Form("ID%d   VAL %14.8f %14.8f %14.8f %14.8f ERR %14.8f\n", ib, stat, acc, cnt, rscl, err);

        hap_rat->SetBinContent(ib, ap_rat);
        hap_rat->SetBinError  (ib, ap_rat * err);
        hap_flx->SetBinContent(ib, ap_flx);
        hap_flx->SetBinError  (ib, ap_flx * err);
    }

    hpbarp->SetLineColor(kBlue);
    hpbarp->SetMarkerColor(kBlue);
    
    hap->SetLineColor(kBlue);
    hap->SetMarkerColor(kBlue);

    hap_rat->SetLineColor(kRed);
    hap_rat->SetMarkerColor(kRed);
    hap_rat->SetMarkerStyle(20);
    hap_rat->GetXaxis()->SetTitle("|Rigidity| [GV]");
    hap_rat->GetYaxis()->SetTitle("#Phi_{#bar{p}}/#Phi_{p}");
    
    hap_flx->SetLineColor(kRed);
    hap_flx->SetMarkerColor(kRed);
    hap_flx->SetMarkerStyle(20);
    hap_flx->GetXaxis()->SetTitle("|Rigidity| [GV]");
    hap_flx->GetYaxis()->SetTitle("#Phi_{#bar{p}} #times R^{2.7} [m^{-2} sr^{-1} s^{-1} GV^{-1.7}]");
   
    TGraphAsymmErrors* gr_ap_rat = new TGraphAsymmErrors(hap_rat);
    TGraphAsymmErrors* gr_ap_flx = new TGraphAsymmErrors(hap_flx);
    for (int ib = gr_ap_rat->GetN()-1; ib >= 0; --ib) if (gr_ap_rat->GetY()[ib] < 1e-10) gr_ap_rat->RemovePoint(ib);
    for (int ib = gr_ap_flx->GetN()-1; ib >= 0; --ib) if (gr_ap_flx->GetY()[ib] < 1e-10) gr_ap_flx->RemovePoint(ib);

    for (int ib = 0; ib < gr_ap_rat->GetN(); ++ib) { gr_ap_rat->SetPointEXhigh(ib, 0); gr_ap_rat->SetPointEXlow(ib, 0); } 
    for (int ib = 0; ib < gr_ap_flx->GetN(); ++ib) { gr_ap_flx->SetPointEXhigh(ib, 0); gr_ap_flx->SetPointEXlow(ib, 0); } 
    
    for (int ib = 0; ib < hap->GetN(); ++ib) { hap->SetPointEXhigh(ib, 0); hap->SetPointEXlow(ib, 0); } 

    TMultiGraph* mg_rat = new TMultiGraph(); 
    mg_rat->SetName("rat");
    mg_rat->Add(hpbarp);
    mg_rat->Add(gr_ap_rat);
    
    TMultiGraph* mg_flx = new TMultiGraph(); 
    mg_flx->SetName("flx");
    mg_flx->Add(hap);
    mg_flx->Add(gr_ap_flx);


    TFile* fout = new TFile("fine.root", "RECREATE");
    hpbarp->Write();
    hap->Write();
    hap_rat->Write();
    hap_flx->Write();
    gr_ap_rat->Write();
    gr_ap_flx->Write();
    mg_rat->Write();
    mg_flx->Write();
    fout->Write();
    fout->Close();
}
