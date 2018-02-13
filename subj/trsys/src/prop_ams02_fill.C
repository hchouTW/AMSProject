//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/ams_home/hchou/AMSCore/prod/18Feb05/src/ClassDef.h"
#include "/ams_home/hchou/AMSCore/prod/18Feb05/src/ClassDef.C"

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();
   
    //MatGeoBoxAms::CreateMatGeoBoxFromG4MatTree();
    //
    //MatFld&& mf1 = MatMgnt::Get(SVecD<3>(0, 0, 175), SVecD<3>(0, 0, 50));
    //MatFld&& mf2 = MatMgnt::Get(SVecD<3>(0, 0, -30), SVecD<3>(0, 0, -100));
    //MatFld&& mf3 = MatMgnt::Get(SVecD<3>(0, 0, 50), SVecD<3>(0, 0, -50));
    //mf1.print();
    //mf2.print();
    //mf3.print();

    //return 0;

    MGConfig::JobOpt opt(argc, argv);

    TChain * dst = new TChain("data");
    for (auto&& file : opt.flist()) dst->Add(file.c_str());

    LIST * fList = new LIST;
    G4MC * fG4mc = (opt.type() == "MC" ) ? new G4MC : nullptr;
    RTI  * fRti  = (opt.type() == "ISS") ? new RTI  : nullptr;
    TRG  * fTrg  = new TRG ;
    TOF  * fTof  = new TOF ;
    ACC  * fAcc  = new ACC ;
    TRK  * fTrk  = new TRK ;
    TRD  * fTrd  = new TRD ;
    RICH * fRich = new RICH;
    ECAL * fEcal = new ECAL;

    dst->SetBranchAddress("list", &fList);
    if (opt.type() == "MC")
        dst->SetBranchAddress("g4mc", &fG4mc);
    //if (opt.type() == "ISS")
    //    dst->SetBranchAddress("rti",  &fRti);
    //dst->SetBranchAddress("trg",  &fTrg);
    dst->SetBranchAddress("tof",  &fTof);
    //dst->SetBranchAddress("acc",  &fAcc);
    dst->SetBranchAddress("trk",  &fTrk);
    //dst->SetBranchAddress("trd",  &fTrd);
    //dst->SetBranchAddress("rich", &fRich);
    //dst->SetBranchAddress("ecal", &fEcal);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //PhyArg::SetOpt(true, false);
    PhyArg::SetOpt(true, true);
    Int_t laySat = 5;
    Int_t layEnd = 6;
    
    TFile * ofle = new TFile(Form("%s/prop_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 50, 0.5, 4000., AxisScale::kLog);
    
    Axis AXcxy("Coo [cm]", 260, -65., 65.);
    Hist * hEvt = Hist::New("hEvt", "hEvt", HistAxis(AXcxy, AXcxy));
    Hist * hVtx = Hist::New("hVtx", "hVtx", HistAxis(AXcxy, AXcxy));

    // Coo
    Axis AXres("res [#mum]", 1600, -300., 300.);
    Hist * hMrx = Hist::New("hMrx", "hMrx", HistAxis(AXmom, AXres));
    Hist * hMry = Hist::New("hMry", "hMry", HistAxis(AXmom, AXres));
    
    Hist * hMrxNN = Hist::New("hMrxNN", "hMrxNN", HistAxis(AXres));
    Hist * hMryNN = Hist::New("hMryNN", "hMryNN", HistAxis(AXres));
    Hist * hMrxN1 = Hist::New("hMrxN1", "hMrxN1", HistAxis(AXres));
    Hist * hMryN1 = Hist::New("hMryN1", "hMryN1", HistAxis(AXres));
    Hist * hMrxN2 = Hist::New("hMrxN2", "hMrxN2", HistAxis(AXres));
    Hist * hMryN2 = Hist::New("hMryN2", "hMryN2", HistAxis(AXres));
    Hist * hMrxN3 = Hist::New("hMrxN3", "hMrxN3", HistAxis(AXres));
    Hist * hMryN3 = Hist::New("hMryN3", "hMryN3", HistAxis(AXres));
    Hist * hMryN4 = Hist::New("hMryN4", "hMryN4", HistAxis(AXres));

    // Prop
    Axis AXcoo("Residual [cm * p#beta/Q * L^-1]", 800, -0.2, 0.2);
    Axis AXagl("Residual [p#beta/Q]", 800, -0.2, 0.2);
    Axis AXels("Eloss [MeV * #beta^{2}/Q^{2}]", 1200, 0.2, 18);
    
    Hist * hMcx = Hist::New("hMcx", "hMcx", HistAxis(AXmom, AXcoo)); // (TH2D) MC: residual x
    Hist * hMcy = Hist::New("hMcy", "hMcy", HistAxis(AXmom, AXcoo)); // (TH2D) MC: residual y
    Hist * hTcx = Hist::New("hTcx", "hTcx", HistAxis(AXmom, AXcoo)); // (TH2D) ToyMC: residual x
    Hist * hTcy = Hist::New("hTcy", "hTcy", HistAxis(AXmom, AXcoo)); // (TH2D) ToyMC: residual y

    Hist * hMux = Hist::New("hMux", "hMux", HistAxis(AXmom, AXagl)); // (TH2D) MC: cosine angle x
    Hist * hMuy = Hist::New("hMuy", "hMuy", HistAxis(AXmom, AXagl)); // (TH2D) MC: cosine angle y
    Hist * hTux = Hist::New("hTux", "hTux", HistAxis(AXmom, AXagl)); // (TH2D) ToyMC: cosine angle x
    Hist * hTuy = Hist::New("hTuy", "hTuy", HistAxis(AXmom, AXagl)); // (TH2D) ToyMC: cosine angle y
    
    Hist * hMcux = Hist::New("hMcux", "hMcux", HistAxis(AXcoo, AXagl)); // (TH2D) MC: residual x vs. cosine angle x
    Hist * hMcuy = Hist::New("hMcuy", "hMcuy", HistAxis(AXcoo, AXagl)); // (TH2D) MC: residual y vs. cosine angle y
    Hist * hTcux = Hist::New("hTcux", "hTcux", HistAxis(AXcoo, AXagl)); // (TH2D) ToyMC: residual x vs. cosine angle x
    Hist * hTcuy = Hist::New("hTcuy", "hTcuy", HistAxis(AXcoo, AXagl)); // (TH2D) ToyMC: residual y vs. cosine angle y

    Hist * hMee = Hist::New("hMee", "hMee", HistAxis(AXmom, AXels)); // (TH2D) MC: kinetic energy difference
    Hist * hTee = Hist::New("hTee", "hTee", HistAxis(AXmom, AXels)); // (TH2D) ToyMC: kinetic energy difference
    
    Axis AXcoo2("Residual [cm * p#beta/Q * L^-1]", 2000, -0.5, 0.5);
    Axis AXagl2("Residual [p#beta/Q]", 2000, -0.5, 0.5);
    Hist * hMcx2 = Hist::New("hMcx2", "hMcx2", HistAxis(AXcoo2)); // (TH1D) MC: residual x
    Hist * hMcy2 = Hist::New("hMcy2", "hMcy2", HistAxis(AXcoo2)); // (TH1D) MC: residual y
    Hist * hTcx2 = Hist::New("hTcx2", "hTcx2", HistAxis(AXcoo2)); // (TH1D) ToyMC: residual x
    Hist * hTcy2 = Hist::New("hTcy2", "hTcy2", HistAxis(AXcoo2)); // (TH1D) ToyMC: residual y
    
    Hist * hMux2 = Hist::New("hMux2", "hMux2", HistAxis(AXagl2)); // (TH1D) MC: cosine angle x
    Hist * hMuy2 = Hist::New("hMuy2", "hMuy2", HistAxis(AXagl2)); // (TH1D) MC: cosine angle y
    Hist * hTux2 = Hist::New("hTux2", "hTux2", HistAxis(AXagl2)); // (TH1D) ToyMC: cosine angle x
    Hist * hTuy2 = Hist::New("hTuy2", "hTuy2", HistAxis(AXagl2)); // (TH1D) ToyMC: cosine angle y

    Long64_t printRate = dst->GetEntries();
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        dst->GetEntry(entry);
        
        if (fTof->numOfBetaH != 1) continue;
        if (!fTof->statusBetaH) continue;
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        if (fTof->betaHPatt != 15) continue;
        //if (fTof->numOfInTimeCls > 4) continue;
        //if ((fTof->extClsN[0]+fTof->extClsN[1]) > 0 || 
        //    (fTof->extClsN[2]+fTof->extClsN[3]) > 1) continue; 
        if (fTof->Zall != 1) continue;
        if (fTof->betaH < 0.) continue;
        
        if (fAcc->clusters.size() != 0) continue;
        
        if ((fTrg->bit&8) != 8) continue;
        
        if (fTrd->numOfTrack != 1 && fTrd->numOfHTrack != 1) continue;
        if (!fTrd->statusKCls[0]) continue;
        if (fTrd->LLReh[0] > 0.3) continue;
        
        TrackInfo& track = fTrk->track;

        // No Interaction
        //if (fG4mc->primVtx.status) {
        //    if (std::fabs(fG4mc->primVtx.coo[2]) < 40.) hVtx->fill(fG4mc->primVtx.coo[0], fG4mc->primVtx.coo[1]);
        //}
        //else continue;

        // REC hit
        HitTRKInfo * rec[9]; std::fill_n(rec, 9, nullptr);
        HitTRKInfo * recU = nullptr;
        HitTRKInfo * recL = nullptr;

        for (auto&& hit : track.hits) {
            if (hit.layJ == laySat) recU = &hit;
            if (hit.layJ == layEnd) recL = &hit;
            rec[hit.layJ-1] = &hit;
        }

        // MC seg from primary particle
        SegPARTMCInfo * mcseg[9]; std::fill_n(mcseg, 9, nullptr);
        SegPARTMCInfo * mcsegU = nullptr;
        SegPARTMCInfo * mcsegL = nullptr;
        for (auto&& seg : fG4mc->primPart.segs) {
            if (seg.dec != 0) continue;
            if (seg.lay == laySat) mcsegU = &seg;
            if (seg.lay == layEnd) mcsegL = &seg;
            mcseg[seg.lay-1] = &seg;
        }
        
        HitTRKMCInfo * mchit[9]; std::fill_n(mchit, 9, nullptr);
        for (auto&& hit : fG4mc->primPart.hits) {
            mchit[hit.layJ-1] = &hit;
        }

        for (Int_t it = 2; it < 8; ++it) {
            if (!rec[it] || !mchit[it]) continue;
            Short_t  ntp[2] = { rec[it]->nsr[0], rec[it]->nsr[1] };
            Double_t res[2] = { rec[it]->coo[0] - mchit[it]->coo[0], rec[it]->coo[1] - mchit[it]->coo[1] };
            
            constexpr Double_t CM2UM = 1.0e4;
            if (ntp[0]!=0) hMrx->fill(mchit[it]->mom, CM2UM * res[0]);
            if (ntp[1]!=0) hMry->fill(mchit[it]->mom, CM2UM * res[1]);
           
            if (mcseg[it]->mom > 5.0) {
                if (ntp[0]>=1) hMrxNN->fill(CM2UM * res[0]);
                if (ntp[1]>=1) hMryNN->fill(CM2UM * res[1]);
                if (ntp[0]==1) hMrxN1->fill(CM2UM * res[0]);
                if (ntp[1]==1) hMryN1->fill(CM2UM * res[1]);
                if (ntp[0]==2) hMrxN2->fill(CM2UM * res[0]);
                if (ntp[1]==2) hMryN2->fill(CM2UM * res[1]);
                if (ntp[0]>=3) hMrxN3->fill(CM2UM * res[0]);
                if (ntp[1]==3) hMryN3->fill(CM2UM * res[1]);
                if (ntp[1]>=4) hMryN4->fill(CM2UM * res[1]);
            }
        }

        if (mcsegU && mcsegL) {
            PhySt part(PartType::Proton);
            part.set_state_with_cos(
                mcsegU->coo[0], mcsegU->coo[1], mcsegU->coo[2],
                mcsegU->dir[0], mcsegU->dir[1], mcsegU->dir[2]
            );
            part.set_mom(mcsegU->mom);
            Double_t mc_mom = part.mom();
            
            MatFld mfld;       // Material information
            PhySt ppst(part);  // Particle Status
            PropMgnt::PropToZ(mcsegL->coo[2], ppst, &mfld); // Propagate to Z with magnetic field
            Double_t len = std::fabs(mcsegL->coo[2]-mcsegU->coo[2]); // Delta Z
            Double_t nrl = mfld.nrl();  // Number of Radiator Length [1]
            Double_t ela = mfld.ela();  // Electron Abundance [mol/cm^2]
            SVecD<3> refc = ppst.c();   // coord
            SVecD<3> refu = ppst.u();   // cosine angle
    
            const Double_t GeV2MeV = 1.0e3;
            Double_t scl_eloss = (part.bta() * part.bta()) / (part.chrg() * part.chrg()) / ela;       // normalized factor (energy loss)
            Double_t scl_mscat = (part.mom() * part.bta()) / std::fabs(part.chrg()) / std::sqrt(nrl); // normalized factor (multiple-scattering)
            
            ppst = part;
            PropMgnt::PropToZWithMC(mcsegL->coo[2], ppst);
            Double_t mc_resc[2] = { mcsegL->coo[0] - refc(0), mcsegL->coo[1] - refc(1) }; // MC: residual xy [cm]
            Double_t mc_resu[2] = { mcsegL->dir[0] - refu(0), mcsegL->dir[1] - refu(1) }; // MC: cosine angle xy [1]
            Double_t mc_elsm    = GeV2MeV * (mcsegU->ke - mcsegL->ke);                    // MC: kinetic energy difference [GeV]
            Double_t tm_resc[2] = { ppst.cx() - refc(0), ppst.cy() - refc(1) };           // ToyMC: residual xy [cm]
            Double_t tm_resu[2] = { ppst.ux() - refu(0), ppst.uy() - refu(1) };           // ToyMC: cosine angle xy [1]
            Double_t tm_elsm    = GeV2MeV * (part.ke() - ppst.ke());                      // ToyMC: kinetic energy difference [GeV]

            hMcx->fill(mc_mom, scl_mscat * mc_resc[0] / len);
            hMcy->fill(mc_mom, scl_mscat * mc_resc[1] / len);
            hMux->fill(mc_mom, scl_mscat * mc_resu[0]);
            hMuy->fill(mc_mom, scl_mscat * mc_resu[1]);
            hMee->fill(mc_mom, scl_eloss * mc_elsm);
            hTcx->fill(mc_mom, scl_mscat * tm_resc[0] / len);
            hTcy->fill(mc_mom, scl_mscat * tm_resc[1] / len);
            hTux->fill(mc_mom, scl_mscat * tm_resu[0]);
            hTuy->fill(mc_mom, scl_mscat * tm_resu[1]);
            hTee->fill(mc_mom, scl_eloss * tm_elsm);
            
            if (mc_mom > 5.0) hMcux->fill(scl_mscat * mc_resc[0] / len, scl_mscat * mc_resu[0]);
            if (mc_mom > 5.0) hMcuy->fill(scl_mscat * mc_resc[1] / len, scl_mscat * mc_resu[1]);
            if (mc_mom > 5.0) hTcux->fill(scl_mscat * tm_resc[0] / len, scl_mscat * tm_resu[0]);
            if (mc_mom > 5.0) hTcuy->fill(scl_mscat * tm_resc[1] / len, scl_mscat * tm_resu[1]);
            
            if (mc_mom > 5.0) hMcx2->fill(scl_mscat * mc_resc[0] / len);
            if (mc_mom > 5.0) hMcy2->fill(scl_mscat * mc_resc[1] / len);
            if (mc_mom > 5.0) hMux2->fill(scl_mscat * mc_resu[0]);
            if (mc_mom > 5.0) hMuy2->fill(scl_mscat * mc_resu[1]);
            if (mc_mom > 5.0) hTcx2->fill(scl_mscat * tm_resc[0] / len);
            if (mc_mom > 5.0) hTcy2->fill(scl_mscat * tm_resc[1] / len);
            if (mc_mom > 5.0) hTux2->fill(scl_mscat * tm_resu[0]);
            if (mc_mom > 5.0) hTuy2->fill(scl_mscat * tm_resu[1]);
            
            Double_t cx = 0.5 * (mcsegU->coo[0] + mcsegL->coo[0]);
            Double_t cy = 0.5 * (mcsegU->coo[1] + mcsegL->coo[1]);
            hEvt->fill(cx, cy);
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

    return 0;
}
