//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/data1/hchou/17Sep08/src/ClassDef.h"
#include "/data1/hchou/17Sep08/src/ClassDef.C"

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    TH1::AddDirectory(true);
    
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
    //RICH * fRich = new RICH;
    //ECAL * fEcal = new ECAL;

    dst->SetBranchAddress("list", &fList);
    if (opt.type() == "MC")
        dst->SetBranchAddress("g4mc", &fG4mc);
    if (opt.type() == "ISS")
        dst->SetBranchAddress("rti",  &fRti);
    dst->SetBranchAddress("trg",  &fTrg);
    dst->SetBranchAddress("tof",  &fTof);
    dst->SetBranchAddress("acc",  &fAcc);
    dst->SetBranchAddress("trk",  &fTrk);
    dst->SetBranchAddress("trd",  &fTrd);
    //dst->SetBranchAddress("rich", &fRich);
    //dst->SetBranchAddress("ecal", &fEcal);









    TFile * ofle = new TFile("fit.root", "RECREATE");
    TH1D * hmom = new TH1D("hmom", "hmom", 100, 20, 2000);
    TH1D * hrx = new TH1D("hrx", "hrx", 1600, -300, 300);
    TH1D * hry = new TH1D("hry", "hry", 1600, -200, 200);
    
    TH2D * hrig = new TH2D("hrig", "hrig", 10, 20, 2000, 400, -0.05, 0.05);

    std::cout << Form("Entries %lld\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        //if (entry > 100) break; // testcode
        dst->GetEntry(entry);
        
        if (fTrk->tracks.size() != 1) continue;
        TrackInfo& track = fTrk->tracks.at(0);

        hmom->Fill(fG4mc->primPart.mom);
        for (auto&& hit : fG4mc->primPart.hits) {
            if (hit.layJ < 3 || hit.layJ > 8) continue;
            for (auto&& trhit : track.hits) {
                if (trhit.layJ != hit.layJ) continue;
                if (trhit.side != 3) continue;
                Double_t resx = 1e4 * (trhit.coo[0] - hit.coo[0]);
                Double_t resy = 1e4 * (trhit.coo[1] - hit.coo[1]);
                hrx->Fill(resx);
                hry->Fill(resy);
            }
        }
       
        std::vector<HitSt> hits;
        for (auto&& trhit : track.hits) {
            if (trhit.layJ < 3 || trhit.layJ > 8) continue;
            if (trhit.side != 3) continue;
            HitSt hit(trhit.coo[0], trhit.coo[1], trhit.coo[2]);
            hits.push_back(hit);
        }
        if (hits.size() < 4) continue;

        PhyTr tr(hits);
        //tr.print();
        tr.fit();

        Double_t diff = (tr.part().irig() - 1./fG4mc->primPart.mom);
        hrig->Fill(fG4mc->primPart.mom, diff);
    }

    ofle->Write();
    ofle->Close();















    if (fList) { delete fList; fList = nullptr; }
    if (fG4mc) { delete fG4mc; fG4mc = nullptr; }
    if (fRti ) { delete fRti ; fRti  = nullptr; }
    if (fTrg ) { delete fTrg ; fTrg  = nullptr; }
    if (fTof ) { delete fTof ; fTof  = nullptr; }
    if (fAcc ) { delete fAcc ; fAcc  = nullptr; }
    if (fTrk ) { delete fTrk ; fTrk  = nullptr; }
    if (fTrd ) { delete fTrd ; fTrd  = nullptr; }
    //if (fRich) { delete fRich; fRich = nullptr; }
    //if (fEcal) { delete fEcal; fEcal = nullptr; }

    return 0;
}
