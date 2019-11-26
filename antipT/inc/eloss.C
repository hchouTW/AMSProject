void eloss() {
    TFile* file = TFile::Open("out/YiMdst.root");
    //TFile* file = TFile::Open("out/bk2_YiAna.root");
    TH2D* hMC      = (TH2D*) file->Get("hEa2MC");
    TH2D* hBethe   = (TH2D*) file->Get("hEa2Bethe");
    TH2D* hLandau  = (TH2D*) file->Get("hEa2Landau");
    TH2D* hMixture = (TH2D*) file->Get("hEa2Mixture");

    TH1D* fMC      = new TH1D("MC",      "MC",      hMC->GetXaxis()->GetNbins(), hMC->GetXaxis()->GetXbins()->GetArray());
    TH1D* fBethe   = new TH1D("Bethe",   "Bethe",   hMC->GetXaxis()->GetNbins(), hMC->GetXaxis()->GetXbins()->GetArray());
    TH1D* fLandau  = new TH1D("Landau",  "Landau",  hMC->GetXaxis()->GetNbins(), hMC->GetXaxis()->GetXbins()->GetArray());
    TH1D* fMixture = new TH1D("Mixture", "Mixture", hMC->GetXaxis()->GetNbins(), hMC->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= hMC->GetXaxis()->GetNbins(); ++ib) {
        double igb = 0.938272297*fMC->GetXaxis()->GetBinCenter(ib);
        double scl = (1.0 + igb * igb);
        TH1D* subMC = hMC     ->ProjectionY(Form("MC%d", ib), ib, ib);
        TH1D* subBB = hBethe  ->ProjectionY(Form("BB%d", ib), ib, ib);
        TH1D* subLL = hLandau ->ProjectionY(Form("LL%d", ib), ib, ib);
        TH1D* subMM = hMixture->ProjectionY(Form("MM%d", ib), ib, ib);
        fMC     ->SetBinContent(ib, scl * subMC->GetBinCenter(subMC->GetMaximumBin()));
        fBethe  ->SetBinContent(ib, scl * subBB->GetBinCenter(subBB->GetMaximumBin()));
        fLandau ->SetBinContent(ib, scl * subLL->GetBinCenter(subLL->GetMaximumBin()));
        fMixture->SetBinContent(ib, scl * subMM->GetBinCenter(subMM->GetMaximumBin()));
        std::cout << Form("BIN %d SCL %14.8f ELOSS    MC %14.8f BB %14.8f LL %14.8f MM %14.8f\n", ib, scl, fMC->GetBinContent(ib), fBethe->GetBinContent(ib) , fLandau->GetBinContent(ib), fMixture->GetBinContent(ib));
    }

    fMC->SetLineColor(kBlack);
    fBethe->SetLineColor(kGreen+2);
    fLandau->SetLineColor(kBlue);
    fMixture->SetLineColor(kRed);

    THStack* hstack = new THStack();
    hstack->Add(fBethe);
    hstack->Add(fLandau);
    hstack->Add(fMC);
    hstack->Add(fMixture);

    TCanvas* cvs = new TCanvas("cvs", "eloss", 1000, 600);
    gPad->SetLogx();
    hstack->Draw("nostack hist");

    hstack->GetHistogram()->GetXaxis()->SetTitle("1/Momentum [c/GeV]");
    hstack->GetHistogram()->GetYaxis()->SetTitle("Energy Loss [MeV]");
    hstack->GetHistogram()->GetXaxis()->SetMoreLogLabels();

    hstack->Draw("nostack hist");
   
    TLegend* leg = gPad->BuildLegend(0.15,0.65,0.35,0.85,"");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    cvs->Modified();

    cvs->SaveAs("out/eloss.eps");
}
