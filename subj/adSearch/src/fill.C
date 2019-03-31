#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKSys.h>

#include "/ams_home/hchou/AMSCore/prod/19Mar29/src/ClassDef.h"

#include "TMultiGraph.h"


int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();

    FLAGS_logtostderr = true;
    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(google::GLOG_FATAL);

    TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/ams_home/hchou/AMSData/magnetic/AMS02Mag.bin");
    TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/ams_home/hchou/AMSData/material");
    
    //TrackSys::Sys::SetEnv("TRACKSys_MagBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/magnetic/AMS02Mag.bin");
    //TrackSys::Sys::SetEnv("TRACKSys_MatBox", "/eos/ams/user/h/hchou/ExternalLibs/DB/material");
    
    //TrackSys::Sys::ShowMsg( TrackSys::Sys::GetEnv("TRACKSys_MagBox") );
    //TrackSys::Sys::ShowMsg( TrackSys::Sys::GetEnv("TRACKSys_MatBox") );
    
    //TrackSys::MagMgnt::Load();
    //TrackSys::MatMgnt::Load();

    MGConfig::JobOpt opt(argc, argv);

    TChain * dst = new TChain("data");
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
    //HYC  * fHyc  = new HYC ;

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
    //dst->SetBranchAddress("hyc",  &fHyc);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    TFile * ofle = new TFile(Form("%s/fill%04ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    PdfEditor editor(WindowSize::kSliceLR, Form("events%04ld", opt.gi()), opt.opath().c_str());

    PartInfo info(PartType::Proton);
    //PartInfo info(PartType::Helium4);
    
    PartInfo::SetDefault(info.type());
    PhyArg::SetOpt(true, true);
    
    std::vector<Double_t> vmom;
    if (info.type() == PartType::Proton) {
        vmom = std::vector<Double_t>( {
               0.50,
               1.00,    1.33,   1.71,   2.15,   2.67, 
               3.29,    4.02,   4.88,   5.90,   7.09,
               8.48,   10.10,  12.00,  14.10,  16.60, 
               19.50,  22.80,  26.70,  31.10,  36.10, 
               41.90,  48.50,  56.10,  64.80,  74.90, 
               93.00, 125.00, 175.00, 259.00, 
               800.00 } );
    }
    if (info.type() == PartType::Helium4) {
        vmom = std::vector<Double_t>( {
               0.50,
               1.00,    1.33,   1.71,   2.15,   2.67, 
               3.29,    4.02,   4.88,   5.90,   7.09,
               8.48,   10.10,  12.00,  14.10,  16.60, 
               19.50,  22.80,  26.70,  31.10,  36.10, 
               41.90,  48.50,  56.10,  64.80,  74.90, 
               93.00, 125.00, 175.00, 259.00, 
               800.00 } );
    }
    
    std::vector<Double_t> vrig;
    for (auto&& val : vmom) vrig.push_back(std::fabs(val / info.chrg()));

    Axis AXmom("Momentum [GeV]", vmom);
    Axis AXrig("Rigidity [GV]", vrig);
    Axis AXirig("1/Rigidity [1/GV]", AXrig, 1, true);
    Axis AXbtas("Beta", 20, 0.95, 0.995);
   
    Axis AXbta("Beta", 400, 0.95, 1.02);
    Hist* hOFFbeta = Hist::New("hOFFbeta", HistAxis(AXrig, AXbta, "Events/Bin"));
    Hist* hNEWbeta = Hist::New("hNEWbeta", HistAxis(AXrig, AXbta, "Events/Bin"));
    
    Hist* hOFFbetaH = Hist::New("hOFFbetaH", HistAxis(AXbta, "Events/Bin"));
    Hist* hNEWbetaH = Hist::New("hNEWbetaH", HistAxis(AXbta, "Events/Bin"));
    
    Hist* hOFFbetaCut = Hist::New("hOFFbetaCut", HistAxis(AXrig, AXbta, "Events/Bin"));
    Hist* hNEWbetaCut = Hist::New("hNEWbetaCut", HistAxis(AXrig, AXbta, "Events/Bin"));
    Hist* hNEWbetaCut2 = Hist::New("hNEWbetaCut2", HistAxis(AXrig, AXbta, "Events/Bin"));
    
    Hist* hOFFbetaHCut = Hist::New("hOFFbetaHCut", HistAxis(AXbta, "Events/Bin"));
    Hist* hNEWbetaHCut = Hist::New("hNEWbetaHCut", HistAxis(AXbta, "Events/Bin"));
    Hist* hNEWbetaHCut2 = Hist::New("hNEWbetaHCut2", HistAxis(AXbta, "Events/Bin"));
    
    Axis AXnh("Nhit", 60, -30, 30);
    Axis AXu("Uniformity", 100, 0, 1.00001);
    Hist* hNEWbetaHCutOthNH = Hist::New("hNEWbetaHCutOthNH", HistAxis(AXbta, AXnh, "Events/Bin"));
    Hist* hNEWbetaHCutCldNH = Hist::New("hNEWbetaHCutCldNH", HistAxis(AXbta, AXnh, "Events/Bin"));
    Hist* hNEWbetaHCutCldu = Hist::New("hNEWbetaHCutCldu", HistAxis(AXbta, AXu, "Events/Bin"));
   
    Axis AXchi("NChi", 100, 0, 20.0);
    Hist* hNEWbetaHCutStnChi = Hist::New("hNEWbetaHCutStnChi", HistAxis(AXbta, AXchi, "Events/Bin"));

    Axis AXc("Coo [cm]", 30, 0, 102);
    Hist* hNEWbetaHCutGxy = Hist::New("hNEWbetaHCutGxy", HistAxis(AXc, "Events/Bin"));
    Hist* hNEWbetaHCutBxy = Hist::New("hNEWbetaHCutBxy", HistAxis(AXc, "Events/Bin"));
   
    Hist* hNEWbetaHCutGu = Hist::New("hNEWbetaHCutGu", HistAxis(AXu, "Events/Bin"));
    Hist* hNEWbetaHCutBu = Hist::New("hNEWbetaHCutBu", HistAxis(AXu, "Events/Bin"));
    
    Hist* hNEWbetaHCutGu3 = Hist::New("hNEWbetaHCutGu3", HistAxis(AXu, "Events/Bin"));
    Hist* hNEWbetaHCutBu3 = Hist::New("hNEWbetaHCutBu3", HistAxis(AXu, "Events/Bin"));
    
    Hist* hNEWbetaHCutGu4 = Hist::New("hNEWbetaHCutGu4", HistAxis(AXu, "Events/Bin"));
    Hist* hNEWbetaHCutBu4 = Hist::New("hNEWbetaHCutBu4", HistAxis(AXu, "Events/Bin"));
    
    Hist* hNEWbetaHCutGchi = Hist::New("hNEWbetaHCutGchi", HistAxis(AXchi, "Events/Bin"));
    Hist* hNEWbetaHCutBchi = Hist::New("hNEWbetaHCutBchi", HistAxis(AXchi, "Events/Bin"));
    
    Axis AXmass("mass [GeV]", 400, 0.0, 5.0);
    Hist* hOFFmass = Hist::New("hOFFmass", HistAxis(AXbtas, AXmass, "Events/Bin"));
    Hist* hNEWmass = Hist::New("hNEWmass", HistAxis(AXbtas, AXmass, "Events/Bin"));

    Hist* hOFFmassCut = Hist::New("hOFFmassCut", HistAxis(AXbtas, AXmass, "Events/Bin"));
    Hist* hNEWmassCut = Hist::New("hNEWmassCut", HistAxis(AXbtas, AXmass, "Events/Bin"));
    Hist* hNEWmassCut2 = Hist::New("hNEWmassCut2", HistAxis(AXbtas, AXmass, "Events/Bin"));

    TGraph grmir;
    for (int it = 0; it <= 360; ++it) {
        double tha = 2.0 * TMath::Pi() * (it / 360.0);        
        grmir.SetPoint(it, 67.00 * std::cos(tha), 67.00 * std::sin(tha));
    }
    grmir.SetLineColor(kBlack);
    TGraph grhole;
    for (int it = -10; it <= 10; ++it) grhole.SetPoint(grhole.GetN(), 31.9*0.1*it, -32.15);
    for (int it = -10; it <= 10; ++it) grhole.SetPoint(grhole.GetN(), 31.9, 32.15*0.1*it);
    for (int it = -10; it <= 10; ++it) grhole.SetPoint(grhole.GetN(), -31.9*0.1*it, 32.15);
    for (int it = -10; it <= 10; ++it) grhole.SetPoint(grhole.GetN(), -31.9, -32.15*0.1*it);

    MGClock::HrsStopwatch hrssw; hrssw.start();
    Long64_t printRate = static_cast<Long64_t>(0.1 * dst->GetEntries());
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) {
            hrssw.stop();
            COUT("Entry %lld/%lld Time %14.8f\n", entry, dst->GetEntries(), hrssw.time());
        }
        dst->GetEntry(entry);
    
        // Reweight (MC)
        Double_t wgt = fList->weight * ((opt.mode() != MGConfig::JobOpt::MODE::MC) ? 1.0 : std::pow(fG4mc->primPart.mom/AXmom.min(), -1.7));
        
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
        if (!cktrIn.status) continue;
        if (cktrIn.rig < 0) continue;
        if (cktrIn.nchi[0] > 8) continue;
        if (cktrIn.nchi[1] > 8) continue;

        ChFitInfo& chfit = fRich->chfit;
        if (!chfit.status || chfit.kind != 1 || chfit.is_bad_tile || !chfit.is_good_geom) continue;

        if (fRich->status && fRich->kind == 1 && fRich->npmt >= 3) {
            double mass = (fRich->beta < 1) ? std::sqrt((cktrIn.rig * cktrIn.rig) * (1.0 / fRich->beta / fRich->beta - 1.0)) : -1.0;
            hOFFbeta->fillH2D(cktrIn.rig, fRich->beta, wgt);
            if (cktrIn.rig > 20) hOFFbetaH->fillH1D(fRich->beta, wgt);
            hOFFmass->fillH2D(fRich->beta, mass, wgt);
            
            bool cut = (fRich->numOfCls <= 3 && fRich->numOfCls >= 1 && fRich->prob > 0.01 && fRich->cstcq < 5 && fRich->eftOfColPE > 0.4);
            if (cut) {
                hOFFbetaCut->fillH2D(cktrIn.rig, fRich->beta, wgt);
                if (cktrIn.rig > 20) hOFFbetaHCut->fillH1D(fRich->beta, wgt);
                hOFFmassCut->fillH2D(fRich->beta, mass, wgt);
            }
        }

        if (chfit.cloud.status) {
            double mass = (chfit.cloud.cbta < 1) ? std::sqrt((cktrIn.rig * cktrIn.rig) * (1.0 / chfit.cloud.cbta / chfit.cloud.cbta - 1.0)) : -1.0;
            hNEWbeta->fillH2D(cktrIn.rig, chfit.cloud.cbta, wgt);
            if (cktrIn.rig > 20) hNEWbetaH->fillH1D(chfit.cloud.cbta, wgt);
            hNEWmass->fillH2D(chfit.cloud.cbta, mass, wgt);

            bool cut_cloud = (chfit.ncld == 1 && chfit.cloud.status && chfit.cloud.nchi < 3.5 && chfit.cloud.misjudge < 3.5 && chfit.cloud.border > 0.35 && chfit.cloud.trace > 0.1 && chfit.cloud.accuracy > 0.99);
            bool cut_stone = (chfit.nstn <= 1 && chfit.stone.nchi < 3.5 && chfit.stone.dist < 3.4);
            bool cut = (cut_cloud && cut_stone && chfit.ntmr == 0 && chfit.ngst == 0);
            if (cut) {
                hNEWbetaCut->fillH2D(cktrIn.rig, chfit.cloud.cbta, wgt);
                if (cktrIn.rig > 20) hNEWbetaHCut->fillH1D(chfit.cloud.cbta, wgt);
                hNEWmassCut->fillH2D(chfit.cloud.cbta, mass, wgt);
                
                if (chfit.nhit_oth - chfit.cloud.nhit <= 0) {
                    hNEWbetaCut2->fillH2D(cktrIn.rig, chfit.cloud.cbta, wgt);
                    if (cktrIn.rig > 20) hNEWbetaHCut2->fillH1D(chfit.cloud.cbta, wgt);
                    hNEWmassCut2->fillH2D(chfit.cloud.cbta, mass, wgt);
                }
                
                if (cktrIn.rig > 20) hNEWbetaHCutOthNH->fillH2D(chfit.cloud.cbta, chfit.nhit_oth - chfit.cloud.nhit, wgt);
                if (cktrIn.rig > 20) hNEWbetaHCutCldNH->fillH2D(chfit.cloud.cbta, chfit.cloud.nhit_dir - chfit.cloud.nhit_rfl, wgt);
                if (cktrIn.rig > 20) hNEWbetaHCutCldu->fillH2D(chfit.cloud.cbta, chfit.cloud.uniform, wgt);
                
                if (cktrIn.rig > 20 && chfit.nstn == 1) hNEWbetaHCutStnChi->fillH2D(chfit.cloud.cbta, chfit.stone.chic, wgt);
            }

            bool is_bad  = chfit.cloud.cbta < 0.975;
            bool is_good = std::fabs(chfit.cloud.cbta-1.0) < 0.005;
            if (cut && cktrIn.rig > 20 && (is_good || is_bad)) {
                Hist* hNEWxy = (is_good ? hNEWbetaHCutGxy : hNEWbetaHCutBxy);
                for (auto&& hit : chfit.hits) {
                    if (hit.cls != 4) continue;
                    hNEWxy->fillH1D(std::hypot(hit.cx, hit.cy), wgt);
                }
                Hist* hNEWu = (is_good ? hNEWbetaHCutGu : hNEWbetaHCutBu);
                hNEWu->fillH1D(chfit.cloud.uniform, wgt);
                
                Hist* hNEWu3 = (is_good ? hNEWbetaHCutGu3 : hNEWbetaHCutBu3);
                if (chfit.cloud.nhit == 3) hNEWu3->fillH1D(chfit.cloud.uniform, wgt);
                
                Hist* hNEWu4 = (is_good ? hNEWbetaHCutGu4 : hNEWbetaHCutBu4);
                if (chfit.cloud.nhit == 4) hNEWu4->fillH1D(chfit.cloud.uniform, wgt);
                
                Hist* hNEWchi = (is_good ? hNEWbetaHCutGchi : hNEWbetaHCutBchi);
                if (chfit.nstn == 1) hNEWchi->fillH1D(chfit.stone.chic, wgt);
            }

            if (cut && cktrIn.rig > 20 && (is_bad || (is_good && TrackSys::Rndm::DecimalUniform() < 0.0002))) {
                TGraph grstn;
                grstn.SetMarkerColor(kBlue);
                grstn.SetMarkerStyle(29);
                grstn.SetMarkerSize(0.4);

                TGraph grcldd;
                TGraph grcldr;
                grcldd.SetMarkerColor(kRed);
                grcldr.SetMarkerColor(kRed);
                grcldd.SetMarkerStyle(20);
                grcldr.SetMarkerStyle(21);
                grcldd.SetMarkerSize(0.8);
                grcldr.SetMarkerSize(0.8);
                
                TGraph grinn;
                TGraph grout;
                grinn.SetMarkerColor(kGreen+1);
                grout.SetMarkerColor(kGreen+1);
                grinn.SetMarkerStyle(33);
                grout.SetMarkerStyle(27);
                grinn.SetMarkerSize(0.8);
                grout.SetMarkerSize(0.8);
                
                TGraph grpart;
                grpart.SetMarkerColor(kBlack);
                grpart.SetMarkerStyle(30);
                grpart.SetMarkerSize(1.0);
                
                TGraph grray;
                grray.SetMarkerColor(kCyan);
                grray.SetMarkerStyle(20);
                grray.SetMarkerSize(0.05);

                grpart.SetPoint(0, chfit.pmtp[0], chfit.pmtp[1]);

                auto&& rays = RichAms::RayTrace(
                              std::array<double, 6>({chfit.radp[0], chfit.radp[1], chfit.radp[2], chfit.radd[0], chfit.radd[1], chfit.radd[2]}), 
                              chfit.cloud.cbta, 
                              chfit.index, 
                              2.5, 
                              chfit.tile);
                for (auto&& ray : rays) {
                    grray.SetPoint(grray.GetN(), ray[0], ray[1]);
                }
                
                for (auto&& hit : chfit.hits) {
                    if (hit.cls == 0) grstn.SetPoint(grstn.GetN(), hit.cx, hit.cy);
                    if (hit.cls == 1 && hit.mode == 0) grcldd.SetPoint(grcldd.GetN(), hit.cx, hit.cy);
                    if (hit.cls == 1 && hit.mode != 0) grcldr.SetPoint(grcldr.GetN(), hit.cx, hit.cy);
                    if (hit.cls == 4 && hit.mode != -1) grinn.SetPoint(grinn.GetN(), hit.cx, hit.cy);
                    if (hit.cls == 4 && hit.mode == -1) grout.SetPoint(grout.GetN(), hit.cx, hit.cy);
                }

                TMultiGraph mg("mg", "mg");
                if (grpart.GetN() > 0) mg.Add(&grpart);
                if (grray.GetN() > 0) mg.Add(&grray);
                if (grstn.GetN() > 0) mg.Add(&grstn);
                if (grcldd.GetN() > 0) mg.Add(&grcldd);
                if (grcldr.GetN() > 0) mg.Add(&grcldr);
                if (grinn.GetN() > 0) mg.Add(&grinn);
                if (grout.GetN() > 0) mg.Add(&grout);

                editor.create();
                mg.Draw("ap");
                mg.GetHistogram()->SetLineColor(0);
                mg.GetHistogram()->GetXaxis()->SetLimits(-80., 80.);
                mg.GetHistogram()->SetMinimum(-80.);
                mg.GetHistogram()->SetMaximum(80.);
                mg.GetHistogram()->GetXaxis()->SetTitle("X [cm]");
                mg.GetHistogram()->GetYaxis()->SetTitle("Y [cm]");
                mg.Draw("ap");
                grmir.Draw("l");
                grhole.Draw("l");
                TitleDraw(Form("(%s) Beta %6.3f (Cloud) NHIT %2d NPMT %2d (Other) INN %2d OUT %2d", is_good ? "GOOD" : "BAD", chfit.cloud.cbta, chfit.cloud.nhit, chfit.cloud.npmt, chfit.nhit_oth_inn, chfit.nhit_oth_out));
                editor.save();
            }
        }
    }
    

    editor.close();

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
    //if (fHyc ) { delete fHyc;  fHyc  = nullptr; }

    google::ShutdownGoogleLogging();
    return 0;
}
