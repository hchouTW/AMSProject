void matmap() {
    TH1::AddDirectory(true);
    using namespace MGROOT;
    using namespace TrackSys;
    using namespace TrackSys::MatAms;
    
    std::string g4mat_file_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/g4mat_AMS02.root";
    TFile * root_file = TFile::Open(g4mat_file_path.c_str());
    if (root_file == nullptr || root_file->IsZombie()) return;
    
    Long64_t n[3];
    Float_t  min[3];
    Float_t  max[3];
    TTree * tree_inf = (TTree*) root_file->Get("g4mat_inf");
    tree_inf->SetBranchAddress("n",   &n);
    tree_inf->SetBranchAddress("min",  min);
    tree_inf->SetBranchAddress("max",  max);
    tree_inf->GetEntry(0);

    Bool_t  mat;
    Int_t   idx[3];
    Float_t coo[3];
    Bool_t  elm[9];
    Float_t den[9];
    TTree * tree_elm = (TTree*) root_file->Get("g4mat_elm");
    tree_elm->SetBranchAddress("mat", &mat);
    tree_elm->SetBranchAddress("idx",  idx);
    tree_elm->SetBranchAddress("coo",  coo);
    tree_elm->SetBranchAddress("elm",  elm);
    tree_elm->SetBranchAddress("den",  den);

    TFile * outf = new TFile("matmap.root", "RECREATE");
    
    TH1F * hMAT = new TH1F("hMAT", "hMAT", 800, -200., 200.);
    TH2F * hMATXZ = new TH2F("hMATXZ", "hMATXZ", 280, -140., 140., 800, -200., 200.);
    TH2F * hMATYZ = new TH2F("hMATYZ", "hMATYZ", 280, -140., 140., 800, -200., 200.);

    TH1F * hELM[9] { nullptr };
    hELM[0] = new TH1F("hELM1", "hELM1", 800, -200., 200.);
    hELM[1] = new TH1F("hELM2", "hELM2", 800, -200., 200.);
    hELM[2] = new TH1F("hELM3", "hELM3", 800, -200., 200.);
    hELM[3] = new TH1F("hELM4", "hELM4", 800, -200., 200.);
    hELM[4] = new TH1F("hELM5", "hELM5", 800, -200., 200.);
    hELM[5] = new TH1F("hELM6", "hELM6", 800, -200., 200.);
    hELM[6] = new TH1F("hELM7", "hELM7", 800, -200., 200.);
    hELM[7] = new TH1F("hELM8", "hELM8", 800, -200., 200.);
    hELM[8] = new TH1F("hELM9", "hELM9", 800, -200., 200.);
    
    TH2F * hRAD  = new TH2F("hRAD ", "hRAD ", RAD_N.at(0), RAD_MIN.at(0), RAD_MAX.at(0), RAD_N.at(1), RAD_MIN.at(1), RAD_MAX.at(1));
    TH2F * hTRL1 = new TH2F("hTRL1", "hTRL1", TRL1_N.at(0), TRL1_MIN.at(0), TRL1_MAX.at(0), TRL1_N.at(1), TRL1_MIN.at(1), TRL1_MAX.at(1));
    TH2F * hUTRD = new TH2F("hUTRD", "hUTRD", UTRD_N.at(0), UTRD_MIN.at(0), UTRD_MAX.at(0), UTRD_N.at(1), UTRD_MIN.at(1), UTRD_MAX.at(1));
    TH2F * hTRD  = new TH2F("hTRD ", "hTRD ", TRD_N.at(0), TRD_MIN.at(0), TRD_MAX.at(0), TRD_N.at(1), TRD_MIN.at(1), TRD_MAX.at(1));
    TH2F * hLTRD = new TH2F("hLTRD", "hLTRD", LTRD_N.at(0), LTRD_MIN.at(0), LTRD_MAX.at(0), LTRD_N.at(1), LTRD_MIN.at(1), LTRD_MAX.at(1));
    TH2F * hUTOF = new TH2F("hUTOF", "hUTOF", UTOF_N.at(0), UTOF_MIN.at(0), UTOF_MAX.at(0), UTOF_N.at(1), UTOF_MIN.at(1), UTOF_MAX.at(1));
    TH2F * hUITR = new TH2F("hUITR", "hUITR", UITR_N.at(0), UITR_MIN.at(0), UITR_MAX.at(0), UITR_N.at(1), UITR_MIN.at(1), UITR_MAX.at(1));
    TH2F * hTRS1 = new TH2F("hTRS1", "hTRS1", TRS1_N.at(0), TRS1_MIN.at(0), TRS1_MAX.at(0), TRS1_N.at(1), TRS1_MIN.at(1), TRS1_MAX.at(1));
    TH2F * hTRS2 = new TH2F("hTRS2", "hTRS2", TRS2_N.at(0), TRS2_MIN.at(0), TRS2_MAX.at(0), TRS2_N.at(1), TRS2_MIN.at(1), TRS2_MAX.at(1));
    TH2F * hTRS3 = new TH2F("hTRS3", "hTRS3", TRS3_N.at(0), TRS3_MIN.at(0), TRS3_MAX.at(0), TRS3_N.at(1), TRS3_MIN.at(1), TRS3_MAX.at(1));
    TH2F * hLITR = new TH2F("hLITR", "hLITR", LITR_N.at(0), LITR_MIN.at(0), LITR_MAX.at(0), LITR_N.at(1), LITR_MIN.at(1), LITR_MAX.at(1));
    TH2F * hLTOF = new TH2F("hLTOF", "hLTOF", LTOF_N.at(0), LTOF_MIN.at(0), LTOF_MAX.at(0), LTOF_N.at(1), LTOF_MIN.at(1), LTOF_MAX.at(1));
    TH2F * hNAF = new TH2F("hNAF", "hNAF", NAF_N.at(0), NAF_MIN.at(0), NAF_MAX.at(0), NAF_N.at(1), NAF_MIN.at(1), NAF_MAX.at(1));
    TH2F * hAGL = new TH2F("hAGL", "hAGL", AGL_N.at(0), AGL_MIN.at(0), AGL_MAX.at(0), AGL_N.at(1), AGL_MIN.at(1), AGL_MAX.at(1));
    TH2F * hPMT  = new TH2F("hPMT ", "hPMT ", PMT_N.at(0), PMT_MIN.at(0), PMT_MAX.at(0), PMT_N.at(1), PMT_MIN.at(1), PMT_MAX.at(1));
    TH2F * hTRL9 = new TH2F("hTRL9", "hTRL9", TRL9_N.at(0), TRL9_MIN.at(0), TRL9_MAX.at(0), TRL9_N.at(1), TRL9_MIN.at(1), TRL9_MAX.at(1));
    TH2F * hECAL = new TH2F("hECAL", "hECAL", ECAL_N.at(0), ECAL_MIN.at(0), ECAL_MAX.at(0), ECAL_N.at(1), ECAL_MIN.at(1), ECAL_MAX.at(1));

    for (Long64_t entry = 0; entry < tree_elm->GetEntries(); ++entry) {
        tree_elm->GetEntry(entry);
        if (!mat) continue;

        hMATXZ->Fill(coo[0], coo[2]);
        hMATYZ->Fill(coo[1], coo[2]);
        if (std::sqrt(coo[0]*coo[0] + coo[1]*coo[1]) < 35.) {
            hMAT->Fill(coo[2]);
            for (Int_t ie = 0; ie < 9; ++ie) {
                if (!elm[ie]) continue;
                hELM[ie]->Fill(coo[2], den[ie]);
            }
        }

        if (coo[2] > RAD_MIN.at(2) && coo[2] < RAD_MAX.at(2)) hRAD ->Fill(coo[0], coo[1]);
        if (coo[2] > TRL1_MIN.at(2) && coo[2] < TRL1_MAX.at(2)) hTRL1->Fill(coo[0], coo[1]);
        if (coo[2] > UTRD_MIN.at(2) && coo[2] < UTRD_MAX.at(2)) hUTRD->Fill(coo[0], coo[1]);
        if (coo[2] > TRD_MIN.at(2) && coo[2] < TRD_MAX.at(2)) hTRD ->Fill(coo[0], coo[1]);
        if (coo[2] > LTRD_MIN.at(2) && coo[2] < LTRD_MAX.at(2)) hLTRD->Fill(coo[0], coo[1]);
        if (coo[2] > UTOF_MIN.at(2) && coo[2] < UTOF_MAX.at(2)) hUTOF->Fill(coo[0], coo[1]);
        if (coo[2] > UITR_MIN.at(2) && coo[2] < UITR_MAX.at(2)) hUITR->Fill(coo[0], coo[1]);
        if (coo[2] > TRS1_MIN.at(2) && coo[2] < TRS1_MAX.at(2)) hTRS1->Fill(coo[0], coo[1]);
        if (coo[2] > TRS2_MIN.at(2) && coo[2] < TRS2_MAX.at(2)) hTRS2->Fill(coo[0], coo[1]);
        if (coo[2] > TRS3_MIN.at(2) && coo[2] < TRS3_MAX.at(2)) hTRS3->Fill(coo[0], coo[1]);
        if (coo[2] > LITR_MIN.at(2) && coo[2] < LITR_MAX.at(2)) hLITR->Fill(coo[0], coo[1]);
        if (coo[2] > LTOF_MIN.at(2) && coo[2] < LTOF_MAX.at(2)) hLTOF->Fill(coo[0], coo[1]);
        if (coo[2] > NAF_MIN.at(2) && coo[2] < NAF_MAX.at(2)) hNAF->Fill(coo[0], coo[1]);
        if (coo[2] > AGL_MIN.at(2) && coo[2] < AGL_MAX.at(2)) hAGL->Fill(coo[0], coo[1]);
        if (coo[2] > PMT_MIN.at(2) && coo[2] < PMT_MAX.at(2)) hPMT ->Fill(coo[0], coo[1]);
        if (coo[2] > TRL9_MIN.at(2) && coo[2] < TRL9_MAX.at(2)) hTRL9->Fill(coo[0], coo[1]);
        if (coo[2] > ECAL_MIN.at(2) && coo[2] < ECAL_MAX.at(2)) hECAL->Fill(coo[0], coo[1]);
    }

    outf->Write();
    outf->Close();
}
