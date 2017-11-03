#ifdef __HAS_AMS_OFFICE_LIBS__


#ifndef __TRACKLibs_MatEnvAms_C__
#define __TRACKLibs_MatEnvAms_C__


namespace TrackSys {


// Set to MatMgnt::Load()
Bool_t MatMgnt::Load() {
    if (is_load_ && reader_ != nullptr) return true;
    is_load_ = false;
    reader_ = nullptr;

    is_load_ = MatGeoBoxAms::Load();
    if (is_load_) reader_ = &MatGeoBoxAms::Reader();
    
    return is_load_;
};


Bool_t MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree() {
    //std::string g4mat_file_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/g4mat_AMS02.root";
    std::string g4mat_file_path = "/data1/hchou/material3/g4mscan.root";
    
    TFile * root_file = TFile::Open(g4mat_file_path.c_str());
    if (root_file == nullptr || root_file->IsZombie()) return false;

    //std::string dir_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/material";
    std::string dir_path = "/data1/hchou/material5";
   
    // TRACKER
    MatGeoBoxCreator creator_TRL1(MatAms::TRL1_N.at(0), MatAms::TRL1_MIN.at(0), MatAms::TRL1_MAX.at(0), MatAms::TRL1_N.at(1), MatAms::TRL1_MIN.at(1), MatAms::TRL1_MAX.at(1), MatAms::TRL1_N.at(2), MatAms::TRL1_MIN.at(2), MatAms::TRL1_MAX.at(2), CSTR_FMT("%s/AMS02TRL1.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRL2(MatAms::TRL2_N.at(0), MatAms::TRL2_MIN.at(0), MatAms::TRL2_MAX.at(0), MatAms::TRL2_N.at(1), MatAms::TRL2_MIN.at(1), MatAms::TRL2_MAX.at(1), MatAms::TRL2_N.at(2), MatAms::TRL2_MIN.at(2), MatAms::TRL2_MAX.at(2), CSTR_FMT("%s/AMS02TRL2.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRL3(MatAms::TRL3_N.at(0), MatAms::TRL3_MIN.at(0), MatAms::TRL3_MAX.at(0), MatAms::TRL3_N.at(1), MatAms::TRL3_MIN.at(1), MatAms::TRL3_MAX.at(1), MatAms::TRL3_N.at(2), MatAms::TRL3_MIN.at(2), MatAms::TRL3_MAX.at(2), CSTR_FMT("%s/AMS02TRL3.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRL4(MatAms::TRL4_N.at(0), MatAms::TRL4_MIN.at(0), MatAms::TRL4_MAX.at(0), MatAms::TRL4_N.at(1), MatAms::TRL4_MIN.at(1), MatAms::TRL4_MAX.at(1), MatAms::TRL4_N.at(2), MatAms::TRL4_MIN.at(2), MatAms::TRL4_MAX.at(2), CSTR_FMT("%s/AMS02TRL4.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRL5(MatAms::TRL5_N.at(0), MatAms::TRL5_MIN.at(0), MatAms::TRL5_MAX.at(0), MatAms::TRL5_N.at(1), MatAms::TRL5_MIN.at(1), MatAms::TRL5_MAX.at(1), MatAms::TRL5_N.at(2), MatAms::TRL5_MIN.at(2), MatAms::TRL5_MAX.at(2), CSTR_FMT("%s/AMS02TRL5.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRL6(MatAms::TRL6_N.at(0), MatAms::TRL6_MIN.at(0), MatAms::TRL6_MAX.at(0), MatAms::TRL6_N.at(1), MatAms::TRL6_MIN.at(1), MatAms::TRL6_MAX.at(1), MatAms::TRL6_N.at(2), MatAms::TRL6_MIN.at(2), MatAms::TRL6_MAX.at(2), CSTR_FMT("%s/AMS02TRL6.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRL7(MatAms::TRL7_N.at(0), MatAms::TRL7_MIN.at(0), MatAms::TRL7_MAX.at(0), MatAms::TRL7_N.at(1), MatAms::TRL7_MIN.at(1), MatAms::TRL7_MAX.at(1), MatAms::TRL7_N.at(2), MatAms::TRL7_MIN.at(2), MatAms::TRL7_MAX.at(2), CSTR_FMT("%s/AMS02TRL7.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRL8(MatAms::TRL8_N.at(0), MatAms::TRL8_MIN.at(0), MatAms::TRL8_MAX.at(0), MatAms::TRL8_N.at(1), MatAms::TRL8_MIN.at(1), MatAms::TRL8_MAX.at(1), MatAms::TRL8_N.at(2), MatAms::TRL8_MIN.at(2), MatAms::TRL8_MAX.at(2), CSTR_FMT("%s/AMS02TRL8.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRL9(MatAms::TRL9_N.at(0), MatAms::TRL9_MIN.at(0), MatAms::TRL9_MAX.at(0), MatAms::TRL9_N.at(1), MatAms::TRL9_MIN.at(1), MatAms::TRL9_MAX.at(1), MatAms::TRL9_N.at(2), MatAms::TRL9_MIN.at(2), MatAms::TRL9_MAX.at(2), CSTR_FMT("%s/AMS02TRL9.bin", dir_path.c_str()));
  
    // TRD
    MatGeoBoxCreator creator_TRDS(MatAms::TRDS_N.at(0), MatAms::TRDS_MIN.at(0), MatAms::TRDS_MAX.at(0), MatAms::TRDS_N.at(1), MatAms::TRDS_MIN.at(1), MatAms::TRDS_MAX.at(1), MatAms::TRDS_N.at(2), MatAms::TRDS_MIN.at(2), MatAms::TRDS_MAX.at(2), CSTR_FMT("%s/AMS02TRDS.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRDU(MatAms::TRDU_N.at(0), MatAms::TRDU_MIN.at(0), MatAms::TRDU_MAX.at(0), MatAms::TRDU_N.at(1), MatAms::TRDU_MIN.at(1), MatAms::TRDU_MAX.at(1), MatAms::TRDU_N.at(2), MatAms::TRDU_MIN.at(2), MatAms::TRDU_MAX.at(2), CSTR_FMT("%s/AMS02TRDU.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRDM(MatAms::TRDM_N.at(0), MatAms::TRDM_MIN.at(0), MatAms::TRDM_MAX.at(0), MatAms::TRDM_N.at(1), MatAms::TRDM_MIN.at(1), MatAms::TRDM_MAX.at(1), MatAms::TRDM_N.at(2), MatAms::TRDM_MIN.at(2), MatAms::TRDM_MAX.at(2), CSTR_FMT("%s/AMS02TRDM.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRDI(MatAms::TRDI_N.at(0), MatAms::TRDI_MIN.at(0), MatAms::TRDI_MAX.at(0), MatAms::TRDI_N.at(1), MatAms::TRDI_MIN.at(1), MatAms::TRDI_MAX.at(1), MatAms::TRDI_N.at(2), MatAms::TRDI_MIN.at(2), MatAms::TRDI_MAX.at(2), CSTR_FMT("%s/AMS02TRDI.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TRDL(MatAms::TRDL_N.at(0), MatAms::TRDL_MIN.at(0), MatAms::TRDL_MAX.at(0), MatAms::TRDL_N.at(1), MatAms::TRDL_MIN.at(1), MatAms::TRDL_MAX.at(1), MatAms::TRDL_N.at(2), MatAms::TRDL_MIN.at(2), MatAms::TRDL_MAX.at(2), CSTR_FMT("%s/AMS02TRDL.bin", dir_path.c_str()));

    // TOF
    MatGeoBoxCreator creator_TOFU(MatAms::TOFU_N.at(0), MatAms::TOFU_MIN.at(0), MatAms::TOFU_MAX.at(0), MatAms::TOFU_N.at(1), MatAms::TOFU_MIN.at(1), MatAms::TOFU_MAX.at(1), MatAms::TOFU_N.at(2), MatAms::TOFU_MIN.at(2), MatAms::TOFU_MAX.at(2), CSTR_FMT("%s/AMS02TOFU.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_TOFL(MatAms::TOFL_N.at(0), MatAms::TOFL_MIN.at(0), MatAms::TOFL_MAX.at(0), MatAms::TOFL_N.at(1), MatAms::TOFL_MIN.at(1), MatAms::TOFL_MAX.at(1), MatAms::TOFL_N.at(2), MatAms::TOFL_MIN.at(2), MatAms::TOFL_MAX.at(2), CSTR_FMT("%s/AMS02TOFL.bin", dir_path.c_str()));

    // RICH
    MatGeoBoxCreator creator_NAF(
            MatAms::NAF_N.at(0), MatAms::NAF_MIN.at(0), MatAms::NAF_MAX.at(0),
            MatAms::NAF_N.at(1), MatAms::NAF_MIN.at(1), MatAms::NAF_MAX.at(1),
            MatAms::NAF_N.at(2), MatAms::NAF_MIN.at(2), MatAms::NAF_MAX.at(2),
            CSTR_FMT("%s/AMS02NAF.bin", dir_path.c_str())
        );
    
    MatGeoBoxCreator creator_AGL(
            MatAms::AGL_N.at(0), MatAms::AGL_MIN.at(0), MatAms::AGL_MAX.at(0),
            MatAms::AGL_N.at(1), MatAms::AGL_MIN.at(1), MatAms::AGL_MAX.at(1),
            MatAms::AGL_N.at(2), MatAms::AGL_MIN.at(2), MatAms::AGL_MAX.at(2),
            CSTR_FMT("%s/AMS02AGL.bin", dir_path.c_str())
        );
    
    // ECAL
    MatGeoBoxCreator creator_ECAL(
            MatAms::ECAL_N.at(0), MatAms::ECAL_MIN.at(0), MatAms::ECAL_MAX.at(0),
            MatAms::ECAL_N.at(1), MatAms::ECAL_MIN.at(1), MatAms::ECAL_MAX.at(1),
            MatAms::ECAL_N.at(2), MatAms::ECAL_MIN.at(2), MatAms::ECAL_MAX.at(2),
            CSTR_FMT("%s/AMS02ECAL.bin", dir_path.c_str())
        );

    // RAD
    MatGeoBoxCreator creator_RAD(
            MatAms::RAD_N.at(0), MatAms::RAD_MIN.at(0), MatAms::RAD_MAX.at(0),
            MatAms::RAD_N.at(1), MatAms::RAD_MIN.at(1), MatAms::RAD_MAX.at(1),
            MatAms::RAD_N.at(2), MatAms::RAD_MIN.at(2), MatAms::RAD_MAX.at(2),
            CSTR_FMT("%s/AMS02RAD.bin", dir_path.c_str())
        );
    
    // PMT
    MatGeoBoxCreator creator_PMT(
            MatAms::PMT_N.at(0), MatAms::PMT_MIN.at(0), MatAms::PMT_MAX.at(0),
            MatAms::PMT_N.at(1), MatAms::PMT_MIN.at(1), MatAms::PMT_MAX.at(1),
            MatAms::PMT_N.at(2), MatAms::PMT_MIN.at(2), MatAms::PMT_MAX.at(2),
            CSTR_FMT("%s/AMS02PMT.bin", dir_path.c_str())
        );
    
    
    // SUPPORT U
    MatGeoBoxCreator creator_SUPU1(MatAms::SUPU1_N.at(0), MatAms::SUPU1_MIN.at(0), MatAms::SUPU1_MAX.at(0), MatAms::SUPU1_N.at(1), MatAms::SUPU1_MIN.at(1), MatAms::SUPU1_MAX.at(1), MatAms::SUPU1_N.at(2), MatAms::SUPU1_MIN.at(2), MatAms::SUPU1_MAX.at(2), CSTR_FMT("%s/AMS02SUPU1.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPU2(MatAms::SUPU2_N.at(0), MatAms::SUPU2_MIN.at(0), MatAms::SUPU2_MAX.at(0), MatAms::SUPU2_N.at(1), MatAms::SUPU2_MIN.at(1), MatAms::SUPU2_MAX.at(1), MatAms::SUPU2_N.at(2), MatAms::SUPU2_MIN.at(2), MatAms::SUPU2_MAX.at(2), CSTR_FMT("%s/AMS02SUPU2.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPU3(MatAms::SUPU3_N.at(0), MatAms::SUPU3_MIN.at(0), MatAms::SUPU3_MAX.at(0), MatAms::SUPU3_N.at(1), MatAms::SUPU3_MIN.at(1), MatAms::SUPU3_MAX.at(1), MatAms::SUPU3_N.at(2), MatAms::SUPU3_MIN.at(2), MatAms::SUPU3_MAX.at(2), CSTR_FMT("%s/AMS02SUPU3.bin", dir_path.c_str()));
    
    // SUPPORT M
    MatGeoBoxCreator creator_SUPM1(MatAms::SUPM1_N.at(0), MatAms::SUPM1_MIN.at(0), MatAms::SUPM1_MAX.at(0), MatAms::SUPM1_N.at(1), MatAms::SUPM1_MIN.at(1), MatAms::SUPM1_MAX.at(1), MatAms::SUPM1_N.at(2), MatAms::SUPM1_MIN.at(2), MatAms::SUPM1_MAX.at(2), CSTR_FMT("%s/AMS02SUPM1.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPM2(MatAms::SUPM2_N.at(0), MatAms::SUPM2_MIN.at(0), MatAms::SUPM2_MAX.at(0), MatAms::SUPM2_N.at(1), MatAms::SUPM2_MIN.at(1), MatAms::SUPM2_MAX.at(1), MatAms::SUPM2_N.at(2), MatAms::SUPM2_MIN.at(2), MatAms::SUPM2_MAX.at(2), CSTR_FMT("%s/AMS02SUPM2.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPM3(MatAms::SUPM3_N.at(0), MatAms::SUPM3_MIN.at(0), MatAms::SUPM3_MAX.at(0), MatAms::SUPM3_N.at(1), MatAms::SUPM3_MIN.at(1), MatAms::SUPM3_MAX.at(1), MatAms::SUPM3_N.at(2), MatAms::SUPM3_MIN.at(2), MatAms::SUPM3_MAX.at(2), CSTR_FMT("%s/AMS02SUPM3.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPM4(MatAms::SUPM4_N.at(0), MatAms::SUPM4_MIN.at(0), MatAms::SUPM4_MAX.at(0), MatAms::SUPM4_N.at(1), MatAms::SUPM4_MIN.at(1), MatAms::SUPM4_MAX.at(1), MatAms::SUPM4_N.at(2), MatAms::SUPM4_MIN.at(2), MatAms::SUPM4_MAX.at(2), CSTR_FMT("%s/AMS02SUPM4.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPM5(MatAms::SUPM5_N.at(0), MatAms::SUPM5_MIN.at(0), MatAms::SUPM5_MAX.at(0), MatAms::SUPM5_N.at(1), MatAms::SUPM5_MIN.at(1), MatAms::SUPM5_MAX.at(1), MatAms::SUPM5_N.at(2), MatAms::SUPM5_MIN.at(2), MatAms::SUPM5_MAX.at(2), CSTR_FMT("%s/AMS02SUPM5.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPM6(MatAms::SUPM6_N.at(0), MatAms::SUPM6_MIN.at(0), MatAms::SUPM6_MAX.at(0), MatAms::SUPM6_N.at(1), MatAms::SUPM6_MIN.at(1), MatAms::SUPM6_MAX.at(1), MatAms::SUPM6_N.at(2), MatAms::SUPM6_MIN.at(2), MatAms::SUPM6_MAX.at(2), CSTR_FMT("%s/AMS02SUPM6.bin", dir_path.c_str()));
    
    // SUPPORT T
    MatGeoBoxCreator creator_SUPT1(MatAms::SUPT1_N.at(0), MatAms::SUPT1_MIN.at(0), MatAms::SUPT1_MAX.at(0), MatAms::SUPT1_N.at(1), MatAms::SUPT1_MIN.at(1), MatAms::SUPT1_MAX.at(1), MatAms::SUPT1_N.at(2), MatAms::SUPT1_MIN.at(2), MatAms::SUPT1_MAX.at(2), CSTR_FMT("%s/AMS02SUPT1.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPT2(MatAms::SUPT2_N.at(0), MatAms::SUPT2_MIN.at(0), MatAms::SUPT2_MAX.at(0), MatAms::SUPT2_N.at(1), MatAms::SUPT2_MIN.at(1), MatAms::SUPT2_MAX.at(1), MatAms::SUPT2_N.at(2), MatAms::SUPT2_MIN.at(2), MatAms::SUPT2_MAX.at(2), CSTR_FMT("%s/AMS02SUPT2.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPT3(MatAms::SUPT3_N.at(0), MatAms::SUPT3_MIN.at(0), MatAms::SUPT3_MAX.at(0), MatAms::SUPT3_N.at(1), MatAms::SUPT3_MIN.at(1), MatAms::SUPT3_MAX.at(1), MatAms::SUPT3_N.at(2), MatAms::SUPT3_MIN.at(2), MatAms::SUPT3_MAX.at(2), CSTR_FMT("%s/AMS02SUPT3.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPT4(MatAms::SUPT4_N.at(0), MatAms::SUPT4_MIN.at(0), MatAms::SUPT4_MAX.at(0), MatAms::SUPT4_N.at(1), MatAms::SUPT4_MIN.at(1), MatAms::SUPT4_MAX.at(1), MatAms::SUPT4_N.at(2), MatAms::SUPT4_MIN.at(2), MatAms::SUPT4_MAX.at(2), CSTR_FMT("%s/AMS02SUPT4.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPT5(MatAms::SUPT5_N.at(0), MatAms::SUPT5_MIN.at(0), MatAms::SUPT5_MAX.at(0), MatAms::SUPT5_N.at(1), MatAms::SUPT5_MIN.at(1), MatAms::SUPT5_MAX.at(1), MatAms::SUPT5_N.at(2), MatAms::SUPT5_MIN.at(2), MatAms::SUPT5_MAX.at(2), CSTR_FMT("%s/AMS02SUPT5.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPT6(MatAms::SUPT6_N.at(0), MatAms::SUPT6_MIN.at(0), MatAms::SUPT6_MAX.at(0), MatAms::SUPT6_N.at(1), MatAms::SUPT6_MIN.at(1), MatAms::SUPT6_MAX.at(1), MatAms::SUPT6_N.at(2), MatAms::SUPT6_MIN.at(2), MatAms::SUPT6_MAX.at(2), CSTR_FMT("%s/AMS02SUPT6.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPT7(MatAms::SUPT7_N.at(0), MatAms::SUPT7_MIN.at(0), MatAms::SUPT7_MAX.at(0), MatAms::SUPT7_N.at(1), MatAms::SUPT7_MIN.at(1), MatAms::SUPT7_MAX.at(1), MatAms::SUPT7_N.at(2), MatAms::SUPT7_MIN.at(2), MatAms::SUPT7_MAX.at(2), CSTR_FMT("%s/AMS02SUPT7.bin", dir_path.c_str()));
    
    // SUPPORT L
    MatGeoBoxCreator creator_SUPL1(MatAms::SUPL1_N.at(0), MatAms::SUPL1_MIN.at(0), MatAms::SUPL1_MAX.at(0), MatAms::SUPL1_N.at(1), MatAms::SUPL1_MIN.at(1), MatAms::SUPL1_MAX.at(1), MatAms::SUPL1_N.at(2), MatAms::SUPL1_MIN.at(2), MatAms::SUPL1_MAX.at(2), CSTR_FMT("%s/AMS02SUPL1.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPL2(MatAms::SUPL2_N.at(0), MatAms::SUPL2_MIN.at(0), MatAms::SUPL2_MAX.at(0), MatAms::SUPL2_N.at(1), MatAms::SUPL2_MIN.at(1), MatAms::SUPL2_MAX.at(1), MatAms::SUPL2_N.at(2), MatAms::SUPL2_MIN.at(2), MatAms::SUPL2_MAX.at(2), CSTR_FMT("%s/AMS02SUPL2.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPL3(MatAms::SUPL3_N.at(0), MatAms::SUPL3_MIN.at(0), MatAms::SUPL3_MAX.at(0), MatAms::SUPL3_N.at(1), MatAms::SUPL3_MIN.at(1), MatAms::SUPL3_MAX.at(1), MatAms::SUPL3_N.at(2), MatAms::SUPL3_MIN.at(2), MatAms::SUPL3_MAX.at(2), CSTR_FMT("%s/AMS02SUPL3.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPL4(MatAms::SUPL4_N.at(0), MatAms::SUPL4_MIN.at(0), MatAms::SUPL4_MAX.at(0), MatAms::SUPL4_N.at(1), MatAms::SUPL4_MIN.at(1), MatAms::SUPL4_MAX.at(1), MatAms::SUPL4_N.at(2), MatAms::SUPL4_MIN.at(2), MatAms::SUPL4_MAX.at(2), CSTR_FMT("%s/AMS02SUPL4.bin", dir_path.c_str()));
    MatGeoBoxCreator creator_SUPL5(MatAms::SUPL5_N.at(0), MatAms::SUPL5_MIN.at(0), MatAms::SUPL5_MAX.at(0), MatAms::SUPL5_N.at(1), MatAms::SUPL5_MIN.at(1), MatAms::SUPL5_MAX.at(1), MatAms::SUPL5_N.at(2), MatAms::SUPL5_MIN.at(2), MatAms::SUPL5_MAX.at(2), CSTR_FMT("%s/AMS02SUPL5.bin", dir_path.c_str()));
    
    
    root_file->cd();

    Bool_t   mat;
    Double_t coo[3];
    Bool_t   elm[9];
    Double_t mol[9];
    TTree * tree_elm = (TTree*) root_file->Get("g4mat_elm");
    tree_elm->SetBranchAddress("mat", &mat);
    tree_elm->SetBranchAddress("coo",  coo);
    tree_elm->SetBranchAddress("elm",  elm);
    tree_elm->SetBranchAddress("mol",  mol);

    COUT("====================================================================\n");
    COUT("===========  MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree  ===========\n");
    COUT("====================================================================\n");
    Long64_t nentries = tree_elm->GetEntries();
    Long64_t printRat = (nentries / 20);
    for (Long64_t entry = 0; entry < nentries; ++entry) {
        if ((entry%printRat) == 0 || (entry == nentries-1)) {
            COUT("Entry %lld/%lld  (%6.2f %)\n", entry+1, nentries, 100. * static_cast<Double_t>(entry+1)/static_cast<Double_t>(nentries));
        }
        tree_elm->GetEntry(entry);
        if (!mat) continue;
        Float_t radius = std::sqrt(coo[0]*coo[0] + coo[1]*coo[1]);
        Float_t absx = std::fabs(coo[0]);
        Float_t absy = std::fabs(coo[1]);
        Float_t minxy = std::min(absx, absy);
        Float_t maxxy = std::max(absx, absy);
    
        // TRACKER
        Bool_t is_inner_tracker = (radius < MatAms::INNER_TRACKER_RADIUS);
        creator_TRL1.fill(coo, elm, mol);
        creator_TRL2.fill(coo, elm, mol);
        creator_TRL3.fill(coo, elm, mol, is_inner_tracker);
        creator_TRL4.fill(coo, elm, mol, is_inner_tracker);
        creator_TRL5.fill(coo, elm, mol, is_inner_tracker);
        creator_TRL6.fill(coo, elm, mol, is_inner_tracker);
        creator_TRL7.fill(coo, elm, mol, is_inner_tracker);
        creator_TRL8.fill(coo, elm, mol, is_inner_tracker);
        creator_TRL9.fill(coo, elm, mol);
       
        // TRD
        Float_t trd_z = coo[2];
        if      (coo[2] > MatAms::TRDU_MAX[2]) trd_z = MatAms::TRDU_MAX[2];
        else if (coo[2] < MatAms::TRDL_MIN[2]) trd_z = MatAms::TRDL_MIN[2];
        Float_t trd_bound = (MatAms::TRDL_RADIUS + (trd_z - MatAms::TRDL_Z) * MatAms::TRD_SLOPE);
        Bool_t is_trd = ( (maxxy < trd_bound) && (((absx + absy) / trd_bound) < MatAms::TRD_FACTOR) );
        
        creator_TRDS.fill(coo, elm, mol, is_trd);
        creator_TRDU.fill(coo, elm, mol, is_trd);
        creator_TRDM.fill(coo, elm, mol, is_trd);
        creator_TRDI.fill(coo, elm, mol, is_trd);
        creator_TRDL.fill(coo, elm, mol, is_trd);
        
        // TOF
        creator_TOFU.fill(coo, elm, mol);
        creator_TOFL.fill(coo, elm, mol);
        
        // RICH
        Bool_t is_naf = (maxxy < MatAms::RICH_BOUND_INNER);
        Bool_t is_agl = (!is_naf && (radius < MatAms::RICH_BOUND_OUTER));
        
        creator_NAF.fill(coo, elm, mol, is_naf);
        creator_AGL.fill(coo, elm, mol, is_agl);
        
        // ECAL
        creator_ECAL.fill(coo, elm, mol);

        // RAD
        creator_RAD.fill(coo, elm, mol);

        // PMT
        creator_PMT.fill(coo, elm, mol);
        
        // SUPPORT U
        creator_SUPU1.fill(coo, elm, mol);
        creator_SUPU2.fill(coo, elm, mol);
        creator_SUPU3.fill(coo, elm, mol);
        
        // SUPPORT M
        creator_SUPM1.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPM2.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPM3.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPM4.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPM5.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPM6.fill(coo, elm, mol, is_inner_tracker);
        
        creator_SUPT1.fill(coo, elm, mol);
        creator_SUPT2.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPT3.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPT4.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPT5.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPT6.fill(coo, elm, mol, is_inner_tracker);
        creator_SUPT7.fill(coo, elm, mol, is_inner_tracker);
       
        // SUPPORT L
        creator_SUPL1.fill(coo, elm, mol);
        creator_SUPL2.fill(coo, elm, mol);
        creator_SUPL3.fill(coo, elm, mol);
        creator_SUPL4.fill(coo, elm, mol);
        creator_SUPL5.fill(coo, elm, mol);
    }
    COUT("====================================================================\n");
    
    creator_TRL1.save_and_close();
    creator_TRL2.save_and_close();
    creator_TRL3.save_and_close();
    creator_TRL4.save_and_close();
    creator_TRL5.save_and_close();
    creator_TRL6.save_and_close();
    creator_TRL7.save_and_close();
    creator_TRL8.save_and_close();
    creator_TRL9.save_and_close();
    
    creator_TRDS.save_and_close();
    creator_TRDU.save_and_close();
    creator_TRDM.save_and_close();
    creator_TRDI.save_and_close();
    creator_TRDL.save_and_close();
    
    creator_TOFU.save_and_close();
    creator_TOFL.save_and_close();
    
    creator_NAF.save_and_close();
    creator_AGL.save_and_close();
    
    creator_ECAL.save_and_close();
    
    creator_RAD.save_and_close();
    
    creator_PMT.save_and_close();
    
    creator_SUPU1.save_and_close();
    creator_SUPU2.save_and_close();
    creator_SUPU3.save_and_close();
    
    creator_SUPM1.save_and_close();
    creator_SUPM2.save_and_close();
    creator_SUPM3.save_and_close();
    creator_SUPM4.save_and_close();
    creator_SUPM5.save_and_close();
    creator_SUPM6.save_and_close();
    
    creator_SUPT1.save_and_close();
    creator_SUPT2.save_and_close();
    creator_SUPT3.save_and_close();
    creator_SUPT4.save_and_close();
    creator_SUPT5.save_and_close();
    creator_SUPT6.save_and_close();
    creator_SUPT7.save_and_close();
    
    creator_SUPL1.save_and_close();
    creator_SUPL2.save_and_close();
    creator_SUPL3.save_and_close();
    creator_SUPL4.save_and_close();
    creator_SUPL5.save_and_close();

    root_file->Close();

    return true;
}


Bool_t MatGeoBoxAms::Load() {
    if (is_load_) return is_load_;
    //std::string g4mat_dir_path = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/material"; // at CERN
    std::string g4mat_dir_path = "/data1/hchou/material5"; // at NCU

    reader_TRL1_.load(STR_FMT("%s/AMS02TRL1.bin" , g4mat_dir_path.c_str()));
    reader_TRL2_.load(STR_FMT("%s/AMS02TRL2.bin" , g4mat_dir_path.c_str()));
    reader_TRL3_.load(STR_FMT("%s/AMS02TRL3.bin" , g4mat_dir_path.c_str()));
    reader_TRL4_.load(STR_FMT("%s/AMS02TRL4.bin" , g4mat_dir_path.c_str()));
    reader_TRL5_.load(STR_FMT("%s/AMS02TRL5.bin" , g4mat_dir_path.c_str()));
    reader_TRL6_.load(STR_FMT("%s/AMS02TRL6.bin" , g4mat_dir_path.c_str()));
    reader_TRL7_.load(STR_FMT("%s/AMS02TRL7.bin" , g4mat_dir_path.c_str()));
    reader_TRL8_.load(STR_FMT("%s/AMS02TRL8.bin" , g4mat_dir_path.c_str()));
    reader_TRL9_.load(STR_FMT("%s/AMS02TRL9.bin" , g4mat_dir_path.c_str()));
    
    reader_TRDS_.load(STR_FMT("%s/AMS02TRDS.bin" , g4mat_dir_path.c_str()));
    reader_TRDU_.load(STR_FMT("%s/AMS02TRDU.bin" , g4mat_dir_path.c_str()));
    reader_TRDM_.load(STR_FMT("%s/AMS02TRDM.bin" , g4mat_dir_path.c_str()));
    reader_TRDI_.load(STR_FMT("%s/AMS02TRDI.bin" , g4mat_dir_path.c_str()));
    reader_TRDL_.load(STR_FMT("%s/AMS02TRDL.bin" , g4mat_dir_path.c_str()));
    
    reader_TOFU_.load(STR_FMT("%s/AMS02TOFU.bin" , g4mat_dir_path.c_str()));
    reader_TOFL_.load(STR_FMT("%s/AMS02TOFL.bin" , g4mat_dir_path.c_str()));
    
    reader_NAF_.load(STR_FMT("%s/AMS02NAF.bin" , g4mat_dir_path.c_str()));
    reader_AGL_.load(STR_FMT("%s/AMS02AGL.bin" , g4mat_dir_path.c_str()));
    
    reader_ECAL_.load(STR_FMT("%s/AMS02ECAL.bin" , g4mat_dir_path.c_str()));
    
    reader_RAD_.load(STR_FMT("%s/AMS02RAD.bin" , g4mat_dir_path.c_str()));
    
    reader_PMT_.load(STR_FMT("%s/AMS02PMT.bin" , g4mat_dir_path.c_str()));
    
    reader_SUPU1_.load(STR_FMT("%s/AMS02SUPU1.bin" , g4mat_dir_path.c_str()));
    reader_SUPU2_.load(STR_FMT("%s/AMS02SUPU2.bin" , g4mat_dir_path.c_str()));
    reader_SUPU3_.load(STR_FMT("%s/AMS02SUPU3.bin" , g4mat_dir_path.c_str()));
    
    reader_SUPM1_.load(STR_FMT("%s/AMS02SUPM1.bin" , g4mat_dir_path.c_str()));
    reader_SUPM2_.load(STR_FMT("%s/AMS02SUPM2.bin" , g4mat_dir_path.c_str()));
    reader_SUPM3_.load(STR_FMT("%s/AMS02SUPM3.bin" , g4mat_dir_path.c_str()));
    reader_SUPM4_.load(STR_FMT("%s/AMS02SUPM4.bin" , g4mat_dir_path.c_str()));
    reader_SUPM5_.load(STR_FMT("%s/AMS02SUPM5.bin" , g4mat_dir_path.c_str()));
    reader_SUPM6_.load(STR_FMT("%s/AMS02SUPM6.bin" , g4mat_dir_path.c_str()));
    
    reader_SUPT1_.load(STR_FMT("%s/AMS02SUPT1.bin" , g4mat_dir_path.c_str()));
    reader_SUPT2_.load(STR_FMT("%s/AMS02SUPT2.bin" , g4mat_dir_path.c_str()));
    reader_SUPT3_.load(STR_FMT("%s/AMS02SUPT3.bin" , g4mat_dir_path.c_str()));
    reader_SUPT4_.load(STR_FMT("%s/AMS02SUPT4.bin" , g4mat_dir_path.c_str()));
    reader_SUPT5_.load(STR_FMT("%s/AMS02SUPT5.bin" , g4mat_dir_path.c_str()));
    reader_SUPT6_.load(STR_FMT("%s/AMS02SUPT6.bin" , g4mat_dir_path.c_str()));
    reader_SUPT7_.load(STR_FMT("%s/AMS02SUPT7.bin" , g4mat_dir_path.c_str()));
    
    reader_SUPL1_.load(STR_FMT("%s/AMS02SUPL1.bin" , g4mat_dir_path.c_str()));
    reader_SUPL2_.load(STR_FMT("%s/AMS02SUPL2.bin" , g4mat_dir_path.c_str()));
    reader_SUPL3_.load(STR_FMT("%s/AMS02SUPL3.bin" , g4mat_dir_path.c_str()));
    reader_SUPL4_.load(STR_FMT("%s/AMS02SUPL4.bin" , g4mat_dir_path.c_str()));
    reader_SUPL5_.load(STR_FMT("%s/AMS02SUPL5.bin" , g4mat_dir_path.c_str()));
   
    reader_.clear();
    
    reader_.push_back(&reader_TRL1_);
    reader_.push_back(&reader_TRL2_);
    reader_.push_back(&reader_TRL3_);
    reader_.push_back(&reader_TRL4_);
    reader_.push_back(&reader_TRL5_);
    reader_.push_back(&reader_TRL6_);
    reader_.push_back(&reader_TRL7_);
    reader_.push_back(&reader_TRL8_);
    reader_.push_back(&reader_TRL9_);
    
    reader_.push_back(&reader_TRDS_);
    reader_.push_back(&reader_TRDU_);
    reader_.push_back(&reader_TRDM_);
    reader_.push_back(&reader_TRDI_);
    reader_.push_back(&reader_TRDL_);
    
    reader_.push_back(&reader_TOFU_);
    reader_.push_back(&reader_TOFL_);
    
    reader_.push_back(&reader_NAF_);
    reader_.push_back(&reader_AGL_);
    
    reader_.push_back(&reader_ECAL_);
    
    reader_.push_back(&reader_RAD_);
    
    reader_.push_back(&reader_PMT_);
    
    reader_.push_back(&reader_SUPU1_);
    reader_.push_back(&reader_SUPU2_);
    reader_.push_back(&reader_SUPU3_);
    
    reader_.push_back(&reader_SUPM1_);
    reader_.push_back(&reader_SUPM2_);
    reader_.push_back(&reader_SUPM3_);
    reader_.push_back(&reader_SUPM4_);
    reader_.push_back(&reader_SUPM5_);
    reader_.push_back(&reader_SUPM6_);
    
    reader_.push_back(&reader_SUPT1_);
    reader_.push_back(&reader_SUPT2_);
    reader_.push_back(&reader_SUPT3_);
    reader_.push_back(&reader_SUPT4_);
    reader_.push_back(&reader_SUPT5_);
    reader_.push_back(&reader_SUPT6_);
    reader_.push_back(&reader_SUPT7_);
    
    reader_.push_back(&reader_SUPL1_);
    reader_.push_back(&reader_SUPL2_);
    reader_.push_back(&reader_SUPL3_);
    reader_.push_back(&reader_SUPL4_);
    reader_.push_back(&reader_SUPL5_);
   
    is_load_ = true;
    for (auto&& reader : reader_) {
        if (!reader->exist()) {
            is_load_ = false;
            break;
        }
    }

    return is_load_;
}


} // namespace TrackSys


#endif // __TRACKLibs_MatEnvAms_C__


#endif // __HAS_AMS_OFFICE_LIBS__
