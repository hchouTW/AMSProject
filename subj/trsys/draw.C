void draw() {
    TChain*  data = new TChain("data");
    //data->Add("root://grid71.phy.ncu.edu.tw//ams02/user/hchou/BT_PR_400.B1082_18Oct11/*.root");
    data->Add("root://grid71.phy.ncu.edu.tw//ams02/user/hchou/BT_EL_180.B1082_18Oct11/*.root");

    Int_t patt = 0;
    //Double_t mom =  400;
    Double_t mom = -180.;
    //
    std::string ckCut = Form("ckTr.status[%d]", patt);
    std::string kfCut = Form("kfTr.status[%d]", patt);
    std::string hcCut = Form("hcTr.status[%d]", patt);
    //std::string ckCut = Form("ckTr.status[%d] && ckTr.nchi[%d][0]    < 8.0 && ckTr.nchi[%d][1]    < 8.0", patt, patt, patt);
    //std::string kfCut = Form("kfTr.status[%d] && kfTr.nchi[%d][0]    < 8.0 && kfTr.nchi[%d][1]    < 8.0", patt, patt, patt);
    //std::string hcCut = Form("hcTr.status[%d] && hcTr.quality[%d][0] < 2.0 && hcTr.quality[%d][1] < 2.0", patt, patt, patt);
    
    data->Draw(Form("1.0/ckTr.rig[%d]    - (1.0/(%f)) >> h1(400, -0.15, 0.15)", patt, mom), ckCut.c_str(), "");
    data->Draw(Form("1.0/kfTr.rig[%d][0] - (1.0/(%f)) >> h2(400, -0.15, 0.15)", patt, mom), kfCut.c_str(), "");
    data->Draw(Form("1.0/hcTr.rig[%d][0] - (1.0/(%f)) >> h3(400, -0.15, 0.15)", patt, mom), hcCut.c_str(), "");

    h1->SetLineColor(kGreen+1);
    h2->SetLineColor(kBlue);
    h3->SetLineColor(kRed);

    Double_t ckCC = static_cast<Double_t>(data->GetEntries(Form("%s && ckTr.rig[%d]    < 0.0", ckCut.c_str(), patt))) / static_cast<Double_t>(data->GetEntries(ckCut.c_str()));
    Double_t kfCC = static_cast<Double_t>(data->GetEntries(Form("%s && kfTr.rig[%d][0] < 0.0", kfCut.c_str(), patt))) / static_cast<Double_t>(data->GetEntries(kfCut.c_str()));
    Double_t hcCC = static_cast<Double_t>(data->GetEntries(Form("%s && hcTr.rig[%d][0] < 0.0", hcCut.c_str(), patt))) / static_cast<Double_t>(data->GetEntries(hcCut.c_str()));

    THStack* stack = new THStack("hist", "hist");
    stack->Add(h1);
    stack->Add(h2);
    stack->Add(h3);
    stack->Draw("nostack hist");

    TF1* fgaus = new TF1("fgaus", "gaus");

    fgaus->SetParameters(10., 0.0, 0.01);
    h1->Fit(fgaus, "q0", "", -0.005, 0.005);
    Double_t ckSgm = fgaus->GetParameter(2);

    fgaus->SetParameters(10., 0.0, 0.01);
    h2->Fit(fgaus, "q0", "", -0.005, 0.005);
    Double_t kfSgm = fgaus->GetParameter(2);
    
    fgaus->SetParameters(10., 0.0, 0.01);
    h3->Fit(fgaus, "q0", "", -0.005, 0.005);
    Double_t hcSgm = fgaus->GetParameter(2);

    std::cout << Form("\nSigma R   CK %14.8f KF %14.8f HC %14.8f\n\n", 1.0, kfSgm/ckSgm, hcSgm/ckSgm);
    std::cout << Form("\nCC Ratio  CK %14.8f KF %14.8f HC %14.8f\n\n", ckCC*100., kfCC*100., hcCC*100.);
}
