#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

//#include "/ams_home/hchou/AMSCore/prod/19Jan21/src/ClassDef.h"
#include "/afs/cern.ch/work/h/hchou/AMSCore/prod/19Jan31/src/ClassDef.h"

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

    FLAGS_logtostderr = true;
    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(google::GLOG_FATAL);

    //TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");
    
    //TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    
    //TrackSys::MagMgnt::Load();
    //TrackSys::MatMgnt::Load();

    MGConfig::JobOpt opt(argc, argv);

    TChain* dst = new TChain("data");
    for (auto&& file : opt.flist()) dst->Add(file.c_str());

    LIST * fList = new LIST;
    G4MC * fG4mc = (opt.mode() == MGConfig::JobOpt::MODE::MC ) ? new G4MC : nullptr;
    RTI  * fRti  = (opt.mode() == MGConfig::JobOpt::MODE::ISS) ? new RTI  : nullptr;
    TRG  * fTrg  = new TRG ;
    TOF  * fTof  = new TOF ;
    ACC  * fAcc  = new ACC ;
    TRK  * fTrk  = new TRK ;
    TRD  * fTrd  = new TRD ;
    RICH * fRich = new RICH;
    ECAL * fEcal = new ECAL;
    HYC  * fHyc  = new HYC ;

    dst->SetBranchAddress("list", &fList);
    if (opt.mode() == MGConfig::JobOpt::MODE::MC)
        dst->SetBranchAddress("g4mc", &fG4mc);
    if (opt.mode() == MGConfig::JobOpt::MODE::ISS)
        dst->SetBranchAddress("rti",  &fRti);
    dst->SetBranchAddress("trg",  &fTrg);
    dst->SetBranchAddress("tof",  &fTof);
    dst->SetBranchAddress("acc",  &fAcc);
    dst->SetBranchAddress("trk",  &fTrk);
    dst->SetBranchAddress("trd",  &fTrd);
    dst->SetBranchAddress("rich", &fRich);
    dst->SetBranchAddress("ecal", &fEcal);
    dst->SetBranchAddress("hyc",  &fHyc);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    TFile* ofle = new TFile(Form("%s/fill%05ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    TTree* tree = dst->CloneTree(0);

    PartInfo info(PartType::Proton);
    //PartInfo info(PartType::Helium4);
    
    PartInfo::SetDefault(info.type());
    PhyArg::SetOpt(true, true);
   
    std::vector<Double_t> vmom;
    if (info.type() == PartType::Proton) {
        vmom = std::vector<Double_t>( {
               0.50,   0.80,
               1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
               3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
               8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
               19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
               41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
               93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00,
               800.00 } );
    }
    if (info.type() == PartType::Helium4) {
        vmom = std::vector<Double_t>( {
               0.50,   0.80,
               1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
               3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
               8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
               19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
               41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
               93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00,
               800.00 } );
    }
    
    std::vector<Double_t> vrig;
    for (auto&& val : vmom) vrig.push_back(std::fabs(val / info.chrg()));

    Axis AXmom("Momentum [GeV]", vmom);
    Axis AXrig("Rigidity [GV]", vrig);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    
    // AntiD
    Axis AXmu("Mass", 500, -5.0, 5.0);
    Axis AXmass("Mass", 500, 0.0, 5.0);
    Axis AXnh("TRD NHit", 200, 0, 200);
    Axis AXns("TRD NSeg", 20, 0, 20);
    Axis AXtrd("TRD Estimator", 200, 0.0, 1.8);
    Axis AXqlt("Quality", 200, -2.0, 10.0);
    
    Hist* hAD_TF_mass    = Hist::New("hAD_TF_mass",    HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_TF_massQ   = Hist::New("hAD_TF_massQ",   HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_TF_CKmass  = Hist::New("hAD_TF_CKmass",  HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_TF_CKmassQ = Hist::New("hAD_TF_CKmassQ", HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_TF_TDnh    = Hist::New("hAD_TF_TDnh",    HistAxis(AXmu, AXnh, "Events/Bin"));
    Hist* hAD_TF_TDns    = Hist::New("hAD_TF_TDns",    HistAxis(AXmu, AXns, "Events/Bin"));
    Hist* hAD_TF_TDllr   = Hist::New("hAD_TF_TDllr",   HistAxis(AXmu, AXtrd, "Events/Bin"));
    Hist* hAD_TF_M1qltx  = Hist::New("hAD_TF_M1qltx",  HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_TF_M1qlty  = Hist::New("hAD_TF_M1qlty",  HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_TF_M1qlt   = Hist::New("hAD_TF_M1qlt",   HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_TF_M2qltx  = Hist::New("hAD_TF_M2qltx",  HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_TF_M2qlty  = Hist::New("hAD_TF_M2qlty",  HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_TF_M2qlt   = Hist::New("hAD_TF_M2qlt",   HistAxis(AXmu, AXqlt, "Events/Bin"));
    
    Hist* hAD_RH_mass    = Hist::New("hAD_RH_mass",    HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_massQ   = Hist::New("hAD_RH_massQ",   HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_massQ1   = Hist::New("hAD_RH_massQ1",   HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_massQ2   = Hist::New("hAD_RH_massQ2",   HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_massQ3   = Hist::New("hAD_RH_massQ3",   HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_massC   = Hist::New("hAD_RH_massC",   HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_CKmass  = Hist::New("hAD_RH_CKmass",  HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_CKmassQ = Hist::New("hAD_RH_CKmassQ", HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_CKmassQ1 = Hist::New("hAD_RH_CKmassQ1", HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_CKmassQ2 = Hist::New("hAD_RH_CKmassQ2", HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_CKmassQ3 = Hist::New("hAD_RH_CKmassQ3", HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_CKmassC = Hist::New("hAD_RH_CKmassC", HistAxis(AXmu, "Events/Bin"));
    Hist* hAD_RH_TDnh    = Hist::New("hAD_RH_TDnh",    HistAxis(AXmu, AXnh, "Events/Bin"));
    Hist* hAD_RH_TDns    = Hist::New("hAD_RH_TDns",    HistAxis(AXmu, AXns, "Events/Bin"));
    Hist* hAD_RH_TDllr   = Hist::New("hAD_RH_TDllr",   HistAxis(AXmu, AXtrd, "Events/Bin"));
    Hist* hAD_RH_M1qltx  = Hist::New("hAD_RH_M1qltx",  HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_RH_M1qlty  = Hist::New("hAD_RH_M1qlty",  HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_RH_M1qlt   = Hist::New("hAD_RH_M1qlt",   HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_RH_M2qltx  = Hist::New("hAD_RH_M2qltx",  HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_RH_M2qlty  = Hist::New("hAD_RH_M2qlty",  HistAxis(AXmu, AXqlt, "Events/Bin"));
    Hist* hAD_RH_M2qlt   = Hist::New("hAD_RH_M2qlt",   HistAxis(AXmu, AXqlt, "Events/Bin"));

    // Antiproton
    Hist* hAP_LPR_flux  = Hist::New("hAP_LPR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_LNR_flux  = Hist::New("hAP_LNR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_LPR_mass  = Hist::New("hAP_LPR_mass",  HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hAP_LNR_mass  = Hist::New("hAP_LNR_mass",  HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hAP_LPR_mass2 = Hist::New("hAP_LPR_mass2", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hAP_LNR_mass2 = Hist::New("hAP_LNR_mass2", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hAP_LPR_CKmass = Hist::New("hAP_LPR_CKmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hAP_LNR_CKmass = Hist::New("hAP_LNR_CKmass", HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hAP_LPR_TDllr = Hist::New("hAD_LPR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_LNR_TDllr = Hist::New("hAD_LNR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_LPR_qltx  = Hist::New("hAD_LPR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_LNR_qltx  = Hist::New("hAD_LNR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_LPR_qlty  = Hist::New("hAD_LPR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_LNR_qlty  = Hist::New("hAD_LNR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));

    Hist* hAP_MPR_flux  = Hist::New("hAP_MPR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_MNR_flux  = Hist::New("hAP_MNR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_MPR_mass  = Hist::New("hAP_MPR_mass",  HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hAP_MNR_mass  = Hist::New("hAP_MNR_mass",  HistAxis(AXrig, AXmass, "Events/Bin"));
    Hist* hAP_MPR_TDllr = Hist::New("hAD_MPR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_MNR_TDllr = Hist::New("hAD_MNR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_MPR_qltx  = Hist::New("hAD_MPR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_MNR_qltx  = Hist::New("hAD_MNR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_MPR_qlty  = Hist::New("hAD_MPR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_MNR_qlty  = Hist::New("hAD_MNR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));

    Hist* hAP_INPR_flux  = Hist::New("hAP_INPR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_INNR_flux  = Hist::New("hAP_INNR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_INPR_TDllr = Hist::New("hAD_INPR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_INNR_TDllr = Hist::New("hAD_INNR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_INPR_qltx  = Hist::New("hAD_INPR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_INNR_qltx  = Hist::New("hAD_INNR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_INPR_qlty  = Hist::New("hAD_INPR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_INNR_qlty  = Hist::New("hAD_INNR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    
    Hist* hAP_L1PR_flux  = Hist::New("hAP_L1PR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_L1NR_flux  = Hist::New("hAP_L1NR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_L1PR_TDllr = Hist::New("hAD_L1PR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_L1NR_TDllr = Hist::New("hAD_L1NR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_L1PR_qltx  = Hist::New("hAD_L1PR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_L1NR_qltx  = Hist::New("hAD_L1NR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_L1PR_qlty  = Hist::New("hAD_L1PR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_L1NR_qlty  = Hist::New("hAD_L1NR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));

    Hist* hAP_L9PR_flux  = Hist::New("hAP_L9PR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_L9NR_flux  = Hist::New("hAP_L9NR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_L9PR_TDllr = Hist::New("hAD_L9PR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_L9NR_TDllr = Hist::New("hAD_L9NR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_L9PR_qltx  = Hist::New("hAD_L9PR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_L9NR_qltx  = Hist::New("hAD_L9NR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_L9PR_qlty  = Hist::New("hAD_L9PR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_L9NR_qlty  = Hist::New("hAD_L9NR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    
    Hist* hAP_FSPR_flux  = Hist::New("hAP_FSPR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_FSNR_flux  = Hist::New("hAP_FSNR_flux",  HistAxis(AXrig, "Events/Bin"));
    Hist* hAP_FSPR_TDllr = Hist::New("hAD_FSPR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_FSNR_TDllr = Hist::New("hAD_FSNR_TDllr", HistAxis(AXrig, AXtrd, "Events/Bin"));
    Hist* hAP_FSPR_qltx  = Hist::New("hAD_FSPR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_FSNR_qltx  = Hist::New("hAD_FSNR_qltx",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_FSPR_qlty  = Hist::New("hAD_FSPR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));
    Hist* hAP_FSNR_qlty  = Hist::New("hAD_FSNR_qlty",  HistAxis(AXrig, AXqlt, "Events/Bin"));


    TrackSys::Sys::HrsStopwatch sw; sw.start();
    Long64_t printRate = dst->GetEntries()/20;
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) {
            sw.stop();
            COUT("Entry %lld/%lld  TIME %10.2f\n", entry, dst->GetEntries(), sw.time());
        }
        dst->GetEntry(entry);
        
        // IGRF RTI
        double maxCF = (opt.mode() == MGConfig::JobOpt::MODE::ISS) ? fRti->maxCfIGRF : 0;
        
        // Geometry (TRK)
        if (fTrk->numOfTrack != 1) continue;
     
        // Geometry (TOF)
        if (fTof->numOfBetaH != 1) continue;
        if (!fTof->statusBetaH) continue;
        
        // Geometry (TRD)
        if (!fTrd->LLRstatus[0]) continue;
        if (fTrd->LLRnh[0] < 10) continue;
        if (fTrd->ITnh[0]  < 6) continue;
        
        // Geometry (ACC)
        if (fAcc->clusters.size() != 0) continue;
        
        // Down-going
        if (fTof->betaH < 0.) continue;

        // Charge
        if (std::abs(info.chrg()) == 1) {
            if (fTof->Qall < 0.8 || fTof->Qall > 1.3) continue;
            if (fTrk->QIn < 0.8 || fTrk->QIn > 1.3) continue;
            if (fTrk->QInMin < 0.7) continue;
            
            if (fTrk->QL2 > 0 && (fTrk->QL2 < 0.7 || fTrk->QL2 > 1.8)) continue;
            if (fTrk->QL1 > 0 && (fTrk->QL1 < 0.7 || fTrk->QL1 > 1.8)) continue;
            if (fTrk->QL9 > 0 && (fTrk->QL9 < 0.7 || fTrk->QL9 > 1.8)) continue;
        }
        if (std::abs(info.chrg()) == 2) {
            if (fTof->Qall < 1.7 || fTof->Qall > 2.4) continue;
            if (fTrk->QIn < 1.7 || fTrk->QIn > 2.4) continue;
            if (fTrk->QInMin < 1.7) continue;
            
            if (fTrk->QL2 > 0 && (fTrk->QL2 < 1.7 || fTrk->QL2 > 2.8)) continue;
            if (fTrk->QL1 > 0 && (fTrk->QL1 < 1.7 || fTrk->QL1 > 2.8)) continue;
            if (fTrk->QL9 > 0 && (fTrk->QL9 < 1.7 || fTrk->QL9 > 2.8)) continue;
        }

        // TOF
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        if (fTof->noiseExtCls != 0) continue;
        if (fTof->numOfInTimeCls > 4) continue;
        
        // TRD
        //if (fTrd.VTXstatus && fTrd.VTXncls >= 15 && fTrd.VTXnseg >= 2) continue;
       
        // ECAL
        if (fEcal->shower.energyE > 0 && fEcal->shower.PisaBDT > -0.2) continue;

        // Track In
        CKTrackInfo& cktrIn = fTrk->cktr.at(0);
        HCTrInfo&    hctrIn = fHyc->trM1.at(0);
        HCTrInfo&    hctrM1In = fHyc->trM1All.at(0);
        HCTrInfo&    hctrM2In = fHyc->trM2All.at(0);
        
        // Track L1
        CKTrackInfo& cktrL1 = fTrk->cktr.at(1);
        HCTrInfo&    hctrL1 = fHyc->trM1.at(1);
        HCTrInfo&    hctrM1L1 = fHyc->trM1All.at(1);
        HCTrInfo&    hctrM2L1 = fHyc->trM2All.at(1);
        
        // Track L9
        CKTrackInfo& cktrL9 = fTrk->cktr.at(2);
        HCTrInfo&    hctrL9 = fHyc->trM1.at(2);
        HCTrInfo&    hctrM1L9 = fHyc->trM1All.at(2);
        HCTrInfo&    hctrM2L9 = fHyc->trM2All.at(2);
        
        // Track Fs
        CKTrackInfo& cktrFs = fTrk->cktr.at(3);
        HCTrInfo&    hctrFs = fHyc->trM1.at(3);
        HCTrInfo&    hctrM1Fs = fHyc->trM1All.at(3);
        HCTrInfo&    hctrM2Fs = fHyc->trM2All.at(3);
       
        // AntiD
        if (cktrIn.status && fTof->betaH > 0.5 && fTof->betaH < 0.8) {
            Double_t ckmass = std::fabs(info.chrg() * cktrIn.rig * std::sqrt(1.0/fTof->betaH/fTof->betaH-1.0));
            Short_t  cksign = TrackSys::Numc::Compare(cktrIn.rig);
            Double_t ckqltx = std::log(cktrIn.nchi[0]);
            Double_t ckqlty = std::log(cktrIn.nchi[1]);
            Bool_t   overcf = (std::fabs(cktrIn.rig) > 1.0 * maxCF);
            Bool_t   llrc   = (fTrd->LLRep[0] > (0.2 + 0.6 * fTof->betaH * fTof->betaH));
            hAD_TF_CKmass->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.01:1.0));
            if (llrc && ckqltx < 2.0 && ckqlty < 2.0) hAD_TF_CKmassQ->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.01:1.0));
        }
        
        if (cktrIn.status && fRich->status && fRich->kind == 0 && fRich->isGood && fRich->beta > 0.95 && fRich->beta < 0.98 && fRich->nhit > 2 && fRich->npmt > 2) {
            Bool_t rhcut[3] = { fRich->prob > 0.01, fRich->numOfExpPE > 2.0, fRich->eftOfColPE > 0.4 };
            
            Short_t  cksign = TrackSys::Numc::Compare(cktrIn.rig);
            Double_t ckmass = std::fabs(info.chrg() * cktrIn.rig * std::sqrt(1.0/fRich->beta/fRich->beta-1.0));
            Double_t ckqltx = std::log(cktrIn.nchi[0]);
            Double_t ckqlty = std::log(cktrIn.nchi[1]);
            Bool_t   overcf = (std::fabs(cktrIn.rig) > 1.0 * maxCF);
            Bool_t   llrc   = (fTrd->LLRep[0] > 0.8);
            if (llrc) hAD_RH_CKmass->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.01:1.0));
            if (llrc && ckqltx < 2.0 && ckqlty < 2.0) hAD_RH_CKmassQ->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.01:1.0));
            if (llrc && ckqltx < 2.0 && ckqlty < 2.0 && rhcut[0]) hAD_RH_CKmassQ1->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.01:1.0));
            if (llrc && ckqltx < 2.0 && ckqlty < 2.0 && rhcut[0] && rhcut[1]) hAD_RH_CKmassQ2->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.01:1.0));
            if (llrc && ckqltx < 2.0 && ckqlty < 2.0 && rhcut[0] && rhcut[1] && rhcut[2]) hAD_RH_CKmassQ3->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.01:1.0));
            if (llrc && overcf && ckqltx < 2.0 && ckqlty < 2.0) hAD_RH_CKmassC->fillH1D(cksign * ckmass, fList->weight * (cksign>0?0.01:1.0));
        }
        
        HCMuInfo& mutr = fHyc->mutr;
        Short_t musign = (mutr.status ? (mutr.rig[0] > 0 ? 1 : -1) : 0);
        Float_t muarig = (mutr.status ? std::fabs(mutr.rig[0]) : 0);
        Bool_t  mugood = (mutr.status && mutr.muQlt[0] < 2.0 && mutr.muQlt[1] < 2.0 && mutr.muQlt[2] < 2.0 && mutr.qlt[0] < 2.0 && mutr.qlt[1] < 2.0 && mutr.qlt[2] < 2.0);
        Bool_t  mullrc = (mutr.status ? (fTrd->LLRep[0] > (0.2 + 0.6 * mutr.bta[0] * mutr.bta[0])) : true);

        if (mutr.status && mutr.bta[0] > 0.5 && mutr.bta[0] < 0.8) {
            Bool_t overcf = (muarig > 1.0 * maxCF);
            hAD_TF_mass->fillH1D(musign * mutr.mass, fList->weight * (musign>0?0.01:1.0));
            if (mugood && mullrc) hAD_TF_massQ->fillH1D(musign * mutr.mass, fList->weight * (musign>0?0.01:1.0));
            if (mugood) hAD_TF_TDnh->fillH2D(musign * mutr.mass, fTrd->numOfCls, fList->weight * (musign>0?0.01:1.0));
            if (mugood) hAD_TF_TDns->fillH2D(musign * mutr.mass, fTrd->numOfSegment, fList->weight * (musign>0?0.01:1.0));
            if (mugood) hAD_TF_TDllr->fillH2D(musign * mutr.mass, fTrd->LLRep[0], fList->weight * (musign>0?0.01:1.0));
            if (mugood && mullrc && hctrM1In.status) hAD_TF_M1qltx->fillH2D(musign * mutr.mass, hctrM1In.qlt[0], fList->weight * (musign>0?0.01:1.0));
            if (mugood && mullrc && hctrM1In.status) hAD_TF_M1qlty->fillH2D(musign * mutr.mass, hctrM1In.qlt[1], fList->weight * (musign>0?0.01:1.0));
            if (mugood && mullrc && hctrM1In.status) hAD_TF_M1qlt ->fillH2D(musign * mutr.mass, hctrM1In.qlt[2], fList->weight * (musign>0?0.01:1.0));
            if (mugood && mullrc && hctrM2In.status) hAD_TF_M2qltx->fillH2D(musign * mutr.mass, hctrM2In.qlt[0], fList->weight * (musign>0?0.01:1.0));
            if (mugood && mullrc && hctrM2In.status) hAD_TF_M2qlty->fillH2D(musign * mutr.mass, hctrM2In.qlt[1], fList->weight * (musign>0?0.01:1.0));
            if (mugood && mullrc && hctrM2In.status) hAD_TF_M2qlt ->fillH2D(musign * mutr.mass, hctrM2In.qlt[2], fList->weight * (musign>0?0.01:1.0));

            // new tree
            if (mugood && (musign * mutr.mass) < -1.5) tree->Fill();
        }
        
        if (hctrIn.status && hctrIn.statusRh && fRich->status && fRich->kind == 0 && fRich->isGood && fRich->beta > 0.95 && fRich->beta < 0.98 && fRich->nhit > 2 && fRich->npmt > 2) {
            Bool_t rhcut[3] = { fRich->prob > 0.01, fRich->numOfExpPE > 2.0, fRich->eftOfColPE > 0.4 };

            Bool_t overcf   = (std::fabs(hctrIn.rig[0]) > 1.0 * maxCF);
            Short_t  rhsign = TrackSys::Numc::Compare(hctrIn.rig[4]);
            Double_t rhmass = std::fabs(info.chrg() * hctrIn.rig[4] * std::sqrt(1.0/fRich->beta/fRich->beta-1.0));
            Double_t qltx   = hctrIn.qlt[0];
            Double_t qlty   = hctrIn.qlt[1];
            Double_t qlt    = hctrIn.qlt[2];
            Bool_t   rhgood = (qltx < 2.0 && qlty < 2.0 && qlt < 2.0);
            Bool_t   llrc   = (fTrd->LLRep[0] > 0.8);
            if (llrc) hAD_RH_mass->fillH1D(rhsign * rhmass, fList->weight * (rhsign>0?0.01:1.0));
            if (rhgood && llrc) hAD_RH_massQ->fillH1D(rhsign * rhmass, fList->weight * (rhsign>0?0.01:1.0));
            if (rhgood && llrc && rhcut[0]) hAD_RH_massQ1->fillH1D(rhsign * rhmass, fList->weight * (rhsign>0?0.01:1.0));
            if (rhgood && llrc && rhcut[0] && rhcut[1]) hAD_RH_massQ2->fillH1D(rhsign * rhmass, fList->weight * (rhsign>0?0.01:1.0));
            if (rhgood && llrc && rhcut[0] && rhcut[1] && rhcut[2]) hAD_RH_massQ3->fillH1D(rhsign * rhmass, fList->weight * (rhsign>0?0.01:1.0));
            if (rhgood && llrc && overcf) hAD_RH_massC->fillH1D(rhsign * rhmass, fList->weight * (rhsign>0?0.01:1.0));
            if (rhgood) hAD_RH_TDnh->fillH2D(rhsign * rhmass, fTrd->numOfCls, fList->weight * (rhsign>0?0.01:1.0));
            if (rhgood) hAD_RH_TDns->fillH2D(rhsign * rhmass, fTrd->numOfSegment, fList->weight * (rhsign>0?0.01:1.0));
            if (rhgood) hAD_RH_TDllr->fillH2D(rhsign * rhmass, fTrd->LLRep[0], fList->weight * (rhsign>0?0.01:1.0));
            if (hctrM1In.status && llrc) hAD_RH_M1qltx->fillH2D(rhsign * rhmass, hctrM1In.qlt[0], fList->weight * (rhsign>0?0.01:1.0));
            if (hctrM1In.status && llrc) hAD_RH_M1qlty->fillH2D(rhsign * rhmass, hctrM1In.qlt[1], fList->weight * (rhsign>0?0.01:1.0));
            if (hctrM1In.status && llrc) hAD_RH_M1qlt ->fillH2D(rhsign * rhmass, hctrM1In.qlt[2], fList->weight * (rhsign>0?0.01:1.0));
            if (hctrM2In.status && llrc) hAD_RH_M2qltx->fillH2D(rhsign * rhmass, hctrM2In.qlt[0], fList->weight * (rhsign>0?0.01:1.0));
            if (hctrM2In.status && llrc) hAD_RH_M2qlty->fillH2D(rhsign * rhmass, hctrM2In.qlt[1], fList->weight * (rhsign>0?0.01:1.0));
            if (hctrM2In.status && llrc) hAD_RH_M2qlt ->fillH2D(rhsign * rhmass, hctrM2In.qlt[2], fList->weight * (rhsign>0?0.01:1.0));
        }
        
        // Antiproton
        if (mutr.status && hctrIn.status) {
            Short_t  sign = TrackSys::Numc::Compare(hctrIn.rig[4]);
            Double_t arig = std::fabs(hctrIn.rig[0]);
            Bool_t   good = (hctrIn.qlt[0] < 2.0 && hctrIn.qlt[1] < 2.0 && hctrIn.qlt[2] < 2.0);
            Bool_t   llrc = (fTrd->LLRep[0] > 0.8);
            Bool_t overcf = (arig > 1.0 * maxCF);
            Double_t ckmass = (fTof->betaH < 1.0 ? std::fabs(info.chrg() * cktrIn.rig * std::sqrt(1.0/fTof->betaH/fTof->betaH-1.0)) : -1);
            if (overcf && sign>0 && mugood && good && mullrc) hAP_LPR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign<0 && mugood && good && mullrc) hAP_LNR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign>0 && mugood && good && mullrc) hAP_LPR_mass->fillH2D(arig, mutr.mass, fList->weight);
            if (overcf && sign<0 && mugood && good && mullrc) hAP_LNR_mass->fillH2D(arig, mutr.mass, fList->weight);
            if (overcf && sign>0 && mugood && good && mullrc && mutr.bta[0] < 0.99) hAP_LPR_mass2->fillH2D(arig, mutr.mass, fList->weight);
            if (overcf && sign<0 && mugood && good && mullrc && mutr.bta[0] < 0.99) hAP_LNR_mass2->fillH2D(arig, mutr.mass, fList->weight);
            if (overcf && sign>0 && good && mullrc && ckmass > 0) hAP_LPR_CKmass->fillH2D(arig, ckmass, fList->weight);
            if (overcf && sign<0 && good && mullrc && ckmass > 0) hAP_LNR_CKmass->fillH2D(arig, ckmass, fList->weight);
            if (overcf && sign>0 && mugood && good) hAP_LPR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign<0 && mugood && good) hAP_LNR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign>0 && mugood && mullrc) hAP_LPR_qltx->fillH2D(arig, hctrIn.qlt[0], fList->weight);
            if (overcf && sign<0 && mugood && mullrc) hAP_LNR_qltx->fillH2D(arig, hctrIn.qlt[0], fList->weight);
            if (overcf && sign>0 && mugood && mullrc) hAP_LPR_qlty->fillH2D(arig, hctrIn.qlt[1], fList->weight);
            if (overcf && sign<0 && mugood && mullrc) hAP_LNR_qlty->fillH2D(arig, hctrIn.qlt[1], fList->weight);
        } 
        
        if (hctrIn.status && hctrIn.statusRh && fRich->status && fRich->kind == 0 && fRich->isGood) {
            Short_t  rhsign = TrackSys::Numc::Compare(hctrIn.rig[4]);
            Double_t rharig = std::fabs(hctrIn.rig[0]);
            Double_t rhmass = (fRich->beta >= 1.0) ? 0.0 : std::fabs(info.chrg() * hctrIn.rig[4] * std::sqrt(1.0/fRich->beta/fRich->beta-1.0));
            Bool_t   good   = (hctrIn.qlt[0] < 2.0 && hctrIn.qlt[1] < 2.0 && hctrIn.qlt[2] < 2.0);
            Bool_t   llrc   = (fTrd->LLRep[0] > 0.8);
            Bool_t   overcf = (rharig > 1.0 * maxCF);
            if (overcf && rhsign>0 && good && llrc) hAP_MPR_flux->fillH1D(rharig, fList->weight);
            if (overcf && rhsign<0 && good && llrc) hAP_MNR_flux->fillH1D(rharig, fList->weight);
            if (overcf && rhsign>0 && good && llrc) hAP_MPR_mass->fillH2D(rharig, rhmass, fList->weight);
            if (overcf && rhsign<0 && good && llrc) hAP_MNR_mass->fillH2D(rharig, rhmass, fList->weight);
            if (overcf && rhsign>0 && good) hAP_MPR_TDllr->fillH2D(rharig, fTrd->LLRep[0], fList->weight);
            if (overcf && rhsign<0 && good) hAP_MNR_TDllr->fillH2D(rharig, fTrd->LLRep[0], fList->weight);
            if (overcf && rhsign>0 && llrc) hAP_MPR_qltx->fillH2D(rharig, hctrIn.qlt[0], fList->weight);
            if (overcf && rhsign<0 && llrc) hAP_MNR_qltx->fillH2D(rharig, hctrIn.qlt[0], fList->weight);
            if (overcf && rhsign>0 && llrc) hAP_MPR_qlty->fillH2D(rharig, hctrIn.qlt[1], fList->weight);
            if (overcf && rhsign<0 && llrc) hAP_MNR_qlty->fillH2D(rharig, hctrIn.qlt[1], fList->weight);
        }

        if (hctrIn.status) {
            HCTrInfo* hctr = &hctrIn;
            Short_t   sign = TrackSys::Numc::Compare(hctr->rig[0]);
            Double_t  arig = std::fabs(hctr->rig[0]);
            Double_t  qltx = hctr->qlt[0];
            Double_t  qlty = hctr->qlt[1];
            Bool_t    good = (hctr->qlt[0] < 2.0 && hctr->qlt[1] < 2.0 && hctr->qlt[2] < 2.0);
            Bool_t    llrc = (fTrd->LLRep[0] > 0.75);
            Bool_t   overcf = (arig > 1.0 * maxCF);
            if (overcf && sign>0 && good && llrc) hAP_INPR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign<0 && good && llrc) hAP_INNR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign>0 && good) hAP_INPR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign<0 && good) hAP_INNR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign>0 && llrc) hAP_INPR_qltx->fillH2D(arig, qltx, fList->weight);
            if (overcf && sign<0 && llrc) hAP_INNR_qltx->fillH2D(arig, qltx, fList->weight);
            if (overcf && sign>0 && llrc) hAP_INPR_qlty->fillH2D(arig, qlty, fList->weight);
            if (overcf && sign<0 && llrc) hAP_INNR_qlty->fillH2D(arig, qlty, fList->weight);
        }
        
        if (hctrL1.status) {
            HCTrInfo* hctr = &hctrL1;
            Short_t   sign = TrackSys::Numc::Compare(hctr->rig[0]);
            Double_t  arig = std::fabs(hctr->rig[0]);
            Double_t  qltx = hctr->qlt[0];
            Double_t  qlty = hctr->qlt[1];
            Bool_t    good = (hctr->qlt[0] < 2.0 && hctr->qlt[1] < 2.0 && hctr->qlt[2] < 2.0);
            Bool_t    llrc = (fTrd->LLRep[0] > 0.75);
            Bool_t overcf = (arig > 1.0 * maxCF);
            if (overcf && sign>0 && good && llrc) hAP_L1PR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign<0 && good && llrc) hAP_L1NR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign>0 && good) hAP_L1PR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign<0 && good) hAP_L1NR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign>0 && llrc) hAP_L1PR_qltx->fillH2D(arig, qltx, fList->weight);
            if (overcf && sign<0 && llrc) hAP_L1NR_qltx->fillH2D(arig, qltx, fList->weight);
            if (overcf && sign>0 && llrc) hAP_L1PR_qlty->fillH2D(arig, qlty, fList->weight);
            if (overcf && sign<0 && llrc) hAP_L1NR_qlty->fillH2D(arig, qlty, fList->weight);
        }
        
        if (hctrL9.status) {
            HCTrInfo* hctr = &hctrL9;
            Short_t   sign = TrackSys::Numc::Compare(hctr->rig[0]);
            Double_t  arig = std::fabs(hctr->rig[0]);
            Double_t  qltx = hctr->qlt[0];
            Double_t  qlty = hctr->qlt[1];
            Bool_t    good = (hctr->qlt[0] < 2.0 && hctr->qlt[1] < 2.0 && hctr->qlt[2] < 2.0);
            Bool_t    llrc = (fTrd->LLRep[0] > 0.75);
            Bool_t overcf = (arig > 1.0 * maxCF);
            if (overcf && sign>0 && good && llrc) hAP_L9PR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign<0 && good && llrc) hAP_L9NR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign>0 && good) hAP_L9PR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign<0 && good) hAP_L9NR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign>0 && llrc) hAP_L9PR_qltx->fillH2D(arig, qltx, fList->weight);
            if (overcf && sign<0 && llrc) hAP_L9NR_qltx->fillH2D(arig, qltx, fList->weight);
            if (overcf && sign>0 && llrc) hAP_L9PR_qlty->fillH2D(arig, qlty, fList->weight);
            if (overcf && sign<0 && llrc) hAP_L9NR_qlty->fillH2D(arig, qlty, fList->weight);
        }
        
        if (hctrFs.status) {
            HCTrInfo* hctr = &hctrFs;
            Short_t   sign = TrackSys::Numc::Compare(hctr->rig[0]);
            Double_t  arig = std::fabs(hctr->rig[0]);
            Double_t  qltx = hctr->qlt[0];
            Double_t  qlty = hctr->qlt[1];
            Bool_t    good = (hctr->qlt[0] < 2.0 && hctr->qlt[1] < 2.0 && hctr->qlt[2] < 2.0);
            Bool_t    llrc = (fTrd->LLRep[0] > 0.75);
            Bool_t overcf = (arig > 1.0 * maxCF);
            if (overcf && sign>0 && good && llrc) hAP_FSPR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign<0 && good && llrc) hAP_FSNR_flux->fillH1D(arig, fList->weight);
            if (overcf && sign>0 && good) hAP_FSPR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign<0 && good) hAP_FSNR_TDllr->fillH2D(arig, fTrd->LLRep[0], fList->weight);
            if (overcf && sign>0 && llrc) hAP_FSPR_qltx->fillH2D(arig, qltx, fList->weight);
            if (overcf && sign<0 && llrc) hAP_FSNR_qltx->fillH2D(arig, qltx, fList->weight);
            if (overcf && sign>0 && llrc) hAP_FSPR_qlty->fillH2D(arig, qlty, fList->weight);
            if (overcf && sign<0 && llrc) hAP_FSNR_qlty->fillH2D(arig, qlty, fList->weight);
        }
    }

    ofle->Write();
    ofle->Close();

    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    if (fList) { delete fList; fList = nullptr; }
    if (fG4mc) { delete fG4mc; fG4mc = nullptr; }
    if (fRti ) { delete fRti ; fRti  = nullptr; }
    if (fTrg ) { delete fTrg ; fTrg  = nullptr; }
    if (fTof ) { delete fTof ; fTof  = nullptr; }
    if (fAcc ) { delete fAcc ; fAcc  = nullptr; }
    if (fTrk ) { delete fTrk ; fTrk  = nullptr; }
    if (fTrd ) { delete fTrd ; fTrd  = nullptr; }
    if (fRich) { delete fRich; fRich = nullptr; }
    if (fEcal) { delete fEcal; fEcal = nullptr; }
    if (fHyc ) { delete fHyc;  fHyc  = nullptr; }

    google::ShutdownGoogleLogging();
    return 0;
}
