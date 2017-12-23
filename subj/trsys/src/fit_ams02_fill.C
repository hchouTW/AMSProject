//#define __HAS_TESTPROP__
//#define __HAS_TESTFIT__
#define __HAS_AMS_OFFICE_LIBS__
#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

#include "/data3/hchou/AMSCore/prod/17Dec12/src/ClassDef.h"
#include "/data3/hchou/AMSCore/prod/17Dec12/src/ClassDef.C"

using namespace std;

int main(int argc, char * argv[]) {
    using namespace MGROOT;
    using namespace TrackSys;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory();
    
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
    dst->SetBranchAddress("rich", &fRich);
    //dst->SetBranchAddress("ecal", &fEcal);
    
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    //---------------------------------------------------------------//
    PhyArg::SetOpt(false, true);
    Bool_t optL1 = false;
    Bool_t optL9 = false;
    
    TFile * ofle = new TFile(Form("%s/fit_ams02_fill%03ld.root", opt.opath().c_str(), opt.gi()), "RECREATE");
    
    Axis AXmom("Momentum [GeV]", 100, 0.5, 4000., AxisScale::kLog);
    //Axis AXmom("Momentum [GeV]", 10, 20., 100., AxisScale::kLog);
    Axis AXimm("1/Momentum [1/GeV]", AXmom, 1, true);

    // Fit
    Axis AXRrso("(1/Rm - 1/Rt) [1/GV]", 2000, -10.0, 10.0);
    Hist * hCKRrso = Hist::New("hCKRrso", "hCKRrso", HistAxis(AXmom, AXRrso));
    Hist * hCNRrso = Hist::New("hCNRrso", "hCNRrso", HistAxis(AXmom, AXRrso));
    Hist * hHCRrso = Hist::New("hHCRrso", "hHCRrso", HistAxis(AXmom, AXRrso));
    
    Hist * hCKRrsoI0 = Hist::New("hCKRrsoI0", "hCKRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI0 = Hist::New("hCNRrsoI0", "hCNRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI0 = Hist::New("hHCRrsoI0", "hHCRrsoI0", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI1 = Hist::New("hCKRrsoI1", "hCKRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI1 = Hist::New("hCNRrsoI1", "hCNRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI1 = Hist::New("hHCRrsoI1", "hHCRrsoI1", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI2 = Hist::New("hCKRrsoI2", "hCKRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI2 = Hist::New("hCNRrsoI2", "hCNRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI2 = Hist::New("hHCRrsoI2", "hHCRrsoI2", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoI3 = Hist::New("hCKRrsoI3", "hCKRrsoI3", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoI3 = Hist::New("hCNRrsoI3", "hCNRrsoI3", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoI3 = Hist::New("hHCRrsoI3", "hHCRrsoI3", HistAxis(AXmom, AXRrso));
    
    Axis AXRchi("Log-Chi-square [1]", 800, -8.0, 8.0);
    Hist * hCKRchi = Hist::New("hCKRchi", "hCKRchi", HistAxis(AXmom, AXRchi));
    Hist * hCNRchi = Hist::New("hCNRchi", "hCNRchi", HistAxis(AXmom, AXRchi));
    Hist * hHCRchi = Hist::New("hHCRchi", "hHCRchi", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRchiI0 = Hist::New("hCKRchiI0", "hCKRchiI0", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiI0 = Hist::New("hCNRchiI0", "hCNRchiI0", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiI0 = Hist::New("hHCRchiI0", "hHCRchiI0", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiI1 = Hist::New("hCKRchiI1", "hCKRchiI1", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiI1 = Hist::New("hCNRchiI1", "hCNRchiI1", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiI1 = Hist::New("hHCRchiI1", "hHCRchiI1", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiI2 = Hist::New("hCKRchiI2", "hCKRchiI2", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiI2 = Hist::New("hCNRchiI2", "hCNRchiI2", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiI2 = Hist::New("hHCRchiI2", "hHCRchiI2", HistAxis(AXmom, AXRchi));
    Hist * hCKRchiI3 = Hist::New("hCKRchiI3", "hCKRchiI3", HistAxis(AXmom, AXRchi));
    Hist * hCNRchiI3 = Hist::New("hCNRchiI3", "hCNRchiI3", HistAxis(AXmom, AXRchi));
    Hist * hHCRchiI3 = Hist::New("hHCRchiI3", "hHCRchiI3", HistAxis(AXmom, AXRchi));
    
    Hist * hCKRrsoCut = Hist::New("hCKRrsoCut", "hCKRrsoCut", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCut = Hist::New("hCNRrsoCut", "hCNRrsoCut", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCut = Hist::New("hHCRrsoCut", "hHCRrsoCut", HistAxis(AXmom, AXRrso));
    
    Hist * hCKRrsoCutI0 = Hist::New("hCKRrsoCutI0", "hCKRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI0 = Hist::New("hCNRrsoCutI0", "hCNRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI0 = Hist::New("hHCRrsoCutI0", "hHCRrsoCutI0", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI1 = Hist::New("hCKRrsoCutI1", "hCKRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI1 = Hist::New("hCNRrsoCutI1", "hCNRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI1 = Hist::New("hHCRrsoCutI1", "hHCRrsoCutI1", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI2 = Hist::New("hCKRrsoCutI2", "hCKRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI2 = Hist::New("hCNRrsoCutI2", "hCNRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI2 = Hist::New("hHCRrsoCutI2", "hHCRrsoCutI2", HistAxis(AXmom, AXRrso));
    Hist * hCKRrsoCutI3 = Hist::New("hCKRrsoCutI3", "hCKRrsoCutI3", HistAxis(AXmom, AXRrso));
    Hist * hCNRrsoCutI3 = Hist::New("hCNRrsoCutI3", "hCNRrsoCutI3", HistAxis(AXmom, AXRrso));
    Hist * hHCRrsoCutI3 = Hist::New("hHCRrsoCutI3", "hHCRrsoCutI3", HistAxis(AXmom, AXRrso));
    
/*    
    Axis AXTOFrso("#betat/#betam-1 [1]", 2000, -0.5, 1.5);
    Hist * hTOFrso = Hist::New("hTOFrso", "hTOFrso", HistAxis(AXmom, AXTOFrso));
    
    Hist * hTOFrsoI0 = Hist::New("hTOFrsoI0", "hTOFrsoI0", HistAxis(AXmom, AXTOFrso));
    Hist * hTOFrsoI1 = Hist::New("hTOFrsoI1", "hTOFrsoI1", HistAxis(AXmom, AXTOFrso));
    Hist * hTOFrsoI2 = Hist::New("hTOFrsoI2", "hTOFrsoI2", HistAxis(AXmom, AXTOFrso));
    Hist * hTOFrsoI3 = Hist::New("hTOFrsoI3", "hTOFrsoI3", HistAxis(AXmom, AXTOFrso));
    
    Axis AXRICHrso("#betat/#betam-1 [1]", 4000, -0.7, 0.4);
    Hist * hRICHrso = Hist::New("hRICHrso", "hRICHrso", HistAxis(AXmom, AXRICHrso));
    
    Hist * hRICHrsoI0 = Hist::New("hRICHrsoI0", "hRICHrsoI0", HistAxis(AXmom, AXRICHrso));
    Hist * hRICHrsoI1 = Hist::New("hRICHrsoI1", "hRICHrsoI1", HistAxis(AXmom, AXRICHrso));
    Hist * hRICHrsoI2 = Hist::New("hRICHrsoI2", "hRICHrsoI2", HistAxis(AXmom, AXRICHrso));
    Hist * hRICHrsoI3 = Hist::New("hRICHrsoI3", "hRICHrsoI3", HistAxis(AXmom, AXRICHrso));
   
    TGraph* fAntiDMass = new TGraph();
    TGraph* fKaonMass = new TGraph();
    TGraph* fPionMass = new TGraph();
    TGraph* fElMass = new TGraph();
    fAntiDMass->SetName("fAntiDMass");
    fKaonMass->SetName("fKaonMass");
    fPionMass->SetName("fPionMass");
    fElMass->SetName("fElMass");
    for (Int_t ibin = 1; ibin <= AXmom.nbin(); ++ibin) {
        PartInfo H(PartType::Proton);
        PartInfo AntiD(PartType::AntiDeuterium);
        PartInfo Kaon(PartType::KaonPlus);
        PartInfo Pion(PartType::PionPlus);
        PartInfo El(PartType::Electron);
        Double_t bincen = AXmom.center(ibin, AxisScale::kLog);
        Double_t sclmass = std::pow((1.0 + (H.mass()/bincen)*(H.mass()/bincen)), -2.0);
        Double_t ad_mass = sclmass * (1.0/bincen/bincen) * (1.0 - (AntiD.mass()/H.mass())*(AntiD.mass()/H.mass()));
        Double_t kn_mass = sclmass * (1.0/bincen/bincen) * (1.0 - (Kaon.mass()/H.mass())*(Kaon.mass()/H.mass()));
        Double_t pi_mass = sclmass * (1.0/bincen/bincen) * (1.0 - (Pion.mass()/H.mass())*(Pion.mass()/H.mass()));
        Double_t el_mass = sclmass * (1.0/bincen/bincen) * (1.0 - (El.mass()/H.mass())*(El.mass()/H.mass()));
        fAntiDMass->SetPoint(ibin-1, bincen, ad_mass);
        fKaonMass->SetPoint(ibin-1, bincen, kn_mass);
        fPionMass->SetPoint(ibin-1, bincen, pi_mass);
        fElMass->SetPoint(ibin-1, bincen, el_mass);
    }
    fAntiDMass->Write();
    fKaonMass->Write();
    fPionMass->Write();
    fElMass->Write();

    Axis AXMrso("1/M^2 [1/GeV^2]", 2000, -2.0, 2.0);
    Hist * hCKMrso = Hist::New("hCKMrso", "hCKMrso", HistAxis(AXmom, AXMrso));
    Hist * hCNMrso = Hist::New("hCNMrso", "hCNMrso", HistAxis(AXmom, AXMrso));
    Hist * hHCMrso = Hist::New("hHCMrso", "hHCMrso", HistAxis(AXmom, AXMrso));
    
    Hist * hCKMrsoI0 = Hist::New("hCKMrsoI0", "hCKMrsoI0", HistAxis(AXmom, AXMrso));
    Hist * hCNMrsoI0 = Hist::New("hCNMrsoI0", "hCNMrsoI0", HistAxis(AXmom, AXMrso));
    Hist * hHCMrsoI0 = Hist::New("hHCMrsoI0", "hHCMrsoI0", HistAxis(AXmom, AXMrso));
    Hist * hCKMrsoI1 = Hist::New("hCKMrsoI1", "hCKMrsoI1", HistAxis(AXmom, AXMrso));
    Hist * hCNMrsoI1 = Hist::New("hCNMrsoI1", "hCNMrsoI1", HistAxis(AXmom, AXMrso));
    Hist * hHCMrsoI1 = Hist::New("hHCMrsoI1", "hHCMrsoI1", HistAxis(AXmom, AXMrso));
    Hist * hCKMrsoI2 = Hist::New("hCKMrsoI2", "hCKMrsoI2", HistAxis(AXmom, AXMrso));
    Hist * hCNMrsoI2 = Hist::New("hCNMrsoI2", "hCNMrsoI2", HistAxis(AXmom, AXMrso));
    Hist * hHCMrsoI2 = Hist::New("hHCMrsoI2", "hHCMrsoI2", HistAxis(AXmom, AXMrso));
    Hist * hCKMrsoI3 = Hist::New("hCKMrsoI3", "hCKMrsoI3", HistAxis(AXmom, AXMrso));
    Hist * hCNMrsoI3 = Hist::New("hCNMrsoI3", "hCNMrsoI3", HistAxis(AXmom, AXMrso));
    Hist * hHCMrsoI3 = Hist::New("hHCMrsoI3", "hHCMrsoI3", HistAxis(AXmom, AXMrso));
    
    Axis AXMRrso("1/M^2 [1/GeV^2]", 2000, -2.0, 2.0);
    Hist * hCKMRrso = Hist::New("hCKMRrso", "hCKMRrso", HistAxis(AXmom, AXMRrso));
    Hist * hCNMRrso = Hist::New("hCNMRrso", "hCNMRrso", HistAxis(AXmom, AXMRrso));
    Hist * hHCMRrso = Hist::New("hHCMRrso", "hHCMRrso", HistAxis(AXmom, AXMRrso));
    
    Hist * hCKMRrsoI0 = Hist::New("hCKMRrsoI0", "hCKMRrsoI0", HistAxis(AXmom, AXMRrso));
    Hist * hCNMRrsoI0 = Hist::New("hCNMRrsoI0", "hCNMRrsoI0", HistAxis(AXmom, AXMRrso));
    Hist * hHCMRrsoI0 = Hist::New("hHCMRrsoI0", "hHCMRrsoI0", HistAxis(AXmom, AXMRrso));
    Hist * hCKMRrsoI1 = Hist::New("hCKMRrsoI1", "hCKMRrsoI1", HistAxis(AXmom, AXMRrso));
    Hist * hCNMRrsoI1 = Hist::New("hCNMRrsoI1", "hCNMRrsoI1", HistAxis(AXmom, AXMRrso));
    Hist * hHCMRrsoI1 = Hist::New("hHCMRrsoI1", "hHCMRrsoI1", HistAxis(AXmom, AXMRrso));
    Hist * hCKMRrsoI2 = Hist::New("hCKMRrsoI2", "hCKMRrsoI2", HistAxis(AXmom, AXMRrso));
    Hist * hCNMRrsoI2 = Hist::New("hCNMRrsoI2", "hCNMRrsoI2", HistAxis(AXmom, AXMRrso));
    Hist * hHCMRrsoI2 = Hist::New("hHCMRrsoI2", "hHCMRrsoI2", HistAxis(AXmom, AXMRrso));
    Hist * hCKMRrsoI3 = Hist::New("hCKMRrsoI3", "hCKMRrsoI3", HistAxis(AXmom, AXMRrso));
    Hist * hCNMRrsoI3 = Hist::New("hCNMRrsoI3", "hCNMRrsoI3", HistAxis(AXmom, AXMRrso));
    Hist * hHCMRrsoI3 = Hist::New("hHCMRrsoI3", "hHCMRrsoI3", HistAxis(AXmom, AXMRrso));
*/    
    Hist * hCKRrso2 = Hist::New("hCKRrso2", "hCKRrso2", HistAxis(AXRrso));
    Hist * hCNRrso2 = Hist::New("hCNRrso2", "hCNRrso2", HistAxis(AXRrso));
    Hist * hHCRrso2 = Hist::New("hHCRrso2", "hHCRrso2", HistAxis(AXRrso));
    
    Hist * hCKflux = Hist::New("hCKflux", "hCKflux", HistAxis(AXimm));
    Hist * hCNflux = Hist::New("hCNflux", "hCNflux", HistAxis(AXimm));
    Hist * hHCflux = Hist::New("hHCflux", "hHCflux", HistAxis(AXimm));
    Hist * hCKflux2 = Hist::New("hCKflux2", "hCKflux2", HistAxis(AXimm));
    Hist * hCNflux2 = Hist::New("hCNflux2", "hCNflux2", HistAxis(AXimm));
    Hist * hHCflux2 = Hist::New("hHCflux2", "hHCflux2", HistAxis(AXimm));

    Long64_t printRate = dst->GetEntries();
    std::cout << Form("\n==== Totally Entries %lld ====\n", dst->GetEntries());
    for (Long64_t entry = 0; entry < dst->GetEntries(); ++entry) {
        if (entry%printRate==0) COUT("Entry %lld/%lld\n", entry, dst->GetEntries());
        dst->GetEntry(entry);
        
        // No Interaction
        Int_t IntType = 0;
        if (fG4mc->primVtx.status) {
            if      (std::fabs(fG4mc->primVtx.vtx[2]) < 55.) IntType = 1;
            else if (fG4mc->primVtx.vtx[2] > 55.)            IntType = 2;
            else if (fG4mc->primVtx.vtx[2] > -80.)           IntType = 3;
        }
       
        if (fTof->betaH < 0.3 || fTof->betaH > 1.3) continue;
        if (fTof->normChisqT > 10.) continue;
        if (fTof->normChisqC > 10.) continue;
        if (fTof->Qall < 0.8 || fTof->Qall > 1.8) continue;
        if (fTof->betaHPatt != 15) continue;
        
        if (fTrk->tracks.size() != 1) continue;
        TrackInfo& track = fTrk->tracks.at(0);
        
        Bool_t hasL1 = false;
        Bool_t hasL9 = false;
        std::vector<HitSt> mhits;
        for (auto&& hit : track.hits) {
            //if (hit.side != 3) continue; // testcode
            if (hit.layJ == 1) hasL1 = true;
            if (hit.layJ == 9) hasL9 = true;
            Bool_t isInn = (hit.layJ >= 2 && hit.layJ <= 8);
            HitSt mhit(hit.side%2==1, hit.side/2==1);
            mhit.set_coo(hit.coo[0], hit.coo[1], hit.coo[2]);
            mhit.set_err(hit.stripSigX.size(), hit.stripSigY.size());
          
            if (isInn) mhits.push_back(mhit);
            else {
                if (optL1 && hit.layJ == 1) mhits.push_back(mhit);
                if (optL9 && hit.layJ == 9) mhits.push_back(mhit);
            }
            
            Bool_t hasMCL1 = false;
            Bool_t hasMCL9 = false;
            for (auto&& mchit : fG4mc->primPart.hits) {
                if (mchit.layJ == 1) hasMCL1 = true;
                if (mchit.layJ == 9) hasMCL9 = true;
            }
            if (hasL1 && !hasMCL1) hasL1 = false;
            if (hasL9 && !hasMCL9) hasL9 = false;
        }
        Short_t cutNHit = 4 + optL1 + optL9;
        if (mhits.size() <= cutNHit) continue;

        if (optL1 && !hasL1) continue;
        if (optL9 && !hasL9) continue;
        Short_t patt = (optL1 + optL9 * 2);

        if (fG4mc->primPart.hits.size() == 0) continue;
        HitTRKMCInfo* topmc = nullptr;
        for (auto&& hit : fG4mc->primPart.hits) {
            //if (patt == 0 && (hit.layJ == 1 || hit.layJ == 9)) continue;
            //if (patt == 1 && (hit.layJ == 9)) continue;
            //if (patt == 2 && (hit.layJ == 1)) continue;
            topmc = &hit;
            break;
        }
        if (topmc == nullptr) continue;

        Double_t mc_mom = topmc->mom;
        Double_t bincen = AXmom.center(AXmom.find(mc_mom), AxisScale::kLog);
       
        //-------------------------------------//
        PhyTr tr(mhits);
        Bool_t succ = tr.fit();
        //Bool_t succ = false;
        Double_t hc_irig = tr.part().irig();
        //-------------------------------------//

        Double_t mc_irig = (fG4mc->primPart.chrg / mc_mom);
        Double_t ck_irig = ((track.status[0][patt]) ? MGMath::ONE/track.rigidity[0][patt] : 0.);
        Double_t cn_irig = ((track.status[1][patt]) ? MGMath::ONE/track.rigidity[1][patt] : 0.);
        
        Double_t ck_lchi = ((track.status[0][patt]) ? std::log(track.chisq[0][patt][1]) : 0.); 
        Double_t cn_lchi = ((track.status[1][patt]) ? std::log(track.chisq[1][patt][1]) : 0.); 
        Double_t hc_lchi = (                 (succ) ? std::log(tr.nchi())               : 0.); 

        Double_t ck_cut = ((track.status[0][patt]) ? (ck_lchi < 2.0) : false); // 96%
        Double_t cn_cut = ((track.status[1][patt]) ? (cn_lchi < 2.0) : false); // 96%
        Double_t hc_cut = (                 (succ) ? (hc_lchi < 0.9) : false); // 96%

/*
        Double_t mc_mass = fG4mc->primPart.mass;
        Double_t mc_ibta = std::sqrt(1.0 + (mc_mass/mc_mom)*(mc_mass/mc_mom));
        Double_t tofibta = (1.0 / fTof->betaH);
        Double_t richibta = (fRich->status) ? (1.0/fRich->beta) : 1.0;

        Double_t sclmass = std::pow((1.0 + (mc_mass/bincen)*(mc_mass/bincen)), -2.0); 
        Double_t tofmass = (tofibta*tofibta-1.0)*(1.0/mc_mass/mc_mass);
        Double_t ck_mass = ((track.status[0][patt]) ? sclmass * ((ck_irig*ck_irig) - tofmass) : 0.); 
        Double_t cn_mass = ((track.status[1][patt]) ? sclmass * ((cn_irig*cn_irig) - tofmass) : 0.); 
        Double_t hc_mass = (                 (succ) ? sclmass * ((hc_irig*hc_irig) - tofmass) : 0.); 
        
        Double_t richmass = (richibta*richibta-1.0)*(1.0/mc_mass/mc_mass);
        Double_t ck_massr = ((track.status[0][patt] && fRich->status) ? sclmass * ((ck_irig*ck_irig) - richmass) : 0.); 
        Double_t cn_massr = ((track.status[1][patt] && fRich->status) ? sclmass * ((cn_irig*cn_irig) - richmass) : 0.); 
        Double_t hc_massr = (                 (succ && fRich->status) ? sclmass * ((hc_irig*hc_irig) - richmass) : 0.); 

        hTOFrso->fill(mc_mom, (tofibta/mc_ibta - 1.0));
        if (IntType == 0) hTOFrsoI0->fill(mc_mom, (tofibta/mc_ibta - 1.0));
        if (IntType == 1) hTOFrsoI1->fill(mc_mom, (tofibta/mc_ibta - 1.0));
        if (IntType == 2) hTOFrsoI2->fill(mc_mom, (tofibta/mc_ibta - 1.0));
        if (IntType == 3) hTOFrsoI3->fill(mc_mom, (tofibta/mc_ibta - 1.0));

        if (fRich->status && fRich->kindOfRad == 0) hRICHrso->fill(mc_mom, (richibta/mc_ibta - 1.0));
        if (fRich->status && fRich->kindOfRad == 0 && IntType == 0) hRICHrsoI0->fill(mc_mom, (richibta/mc_ibta - 1.0));
        if (fRich->status && fRich->kindOfRad == 0 && IntType == 1) hRICHrsoI1->fill(mc_mom, (richibta/mc_ibta - 1.0));
        if (fRich->status && fRich->kindOfRad == 0 && IntType == 2) hRICHrsoI2->fill(mc_mom, (richibta/mc_ibta - 1.0));
        if (fRich->status && fRich->kindOfRad == 0 && IntType == 3) hRICHrsoI3->fill(mc_mom, (richibta/mc_ibta - 1.0));
*/
        if (track.status[0][patt]) hCKRrso->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (track.status[1][patt]) hCNRrso->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (                 succ) hHCRrso->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (track.status[0][patt] && IntType == 0) hCKRrsoI0->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (track.status[1][patt] && IntType == 0) hCNRrsoI0->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (                 succ && IntType == 0) hHCRrsoI0->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (track.status[0][patt] && IntType == 1) hCKRrsoI1->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (track.status[1][patt] && IntType == 1) hCNRrsoI1->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (                 succ && IntType == 1) hHCRrsoI1->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (track.status[0][patt] && IntType == 2) hCKRrsoI2->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (track.status[1][patt] && IntType == 2) hCNRrsoI2->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (                 succ && IntType == 2) hHCRrsoI2->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (track.status[0][patt] && IntType == 3) hCKRrsoI3->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (track.status[1][patt] && IntType == 3) hCNRrsoI3->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (                 succ && IntType == 3) hHCRrsoI3->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (track.status[0][patt]) hCKRchi->fillH2D(mc_mom, ck_lchi);
        if (track.status[1][patt]) hCNRchi->fillH2D(mc_mom, cn_lchi);
        if (                 succ) hHCRchi->fillH2D(mc_mom, hc_lchi);
        
        if (track.status[0][patt] && IntType == 0) hCKRchiI0->fillH2D(mc_mom, ck_lchi);
        if (track.status[1][patt] && IntType == 0) hCNRchiI0->fillH2D(mc_mom, cn_lchi);
        if (                 succ && IntType == 0) hHCRchiI0->fillH2D(mc_mom, hc_lchi);
        if (track.status[0][patt] && IntType == 1) hCKRchiI1->fillH2D(mc_mom, ck_lchi);
        if (track.status[1][patt] && IntType == 1) hCNRchiI1->fillH2D(mc_mom, cn_lchi);
        if (                 succ && IntType == 1) hHCRchiI1->fillH2D(mc_mom, hc_lchi);
        if (track.status[0][patt] && IntType == 2) hCKRchiI2->fillH2D(mc_mom, ck_lchi);
        if (track.status[1][patt] && IntType == 2) hCNRchiI2->fillH2D(mc_mom, cn_lchi);
        if (                 succ && IntType == 2) hHCRchiI2->fillH2D(mc_mom, hc_lchi);
        if (track.status[0][patt] && IntType == 3) hCKRchiI3->fillH2D(mc_mom, ck_lchi);
        if (track.status[1][patt] && IntType == 3) hCNRchiI3->fillH2D(mc_mom, cn_lchi);
        if (                 succ && IntType == 3) hHCRchiI3->fillH2D(mc_mom, hc_lchi);
        
        if (ck_cut) hCKRrsoCut->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut) hCNRrsoCut->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut) hHCRrsoCut->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
        if (ck_cut && IntType == 0) hCKRrsoCutI0->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut && IntType == 0) hCNRrsoCutI0->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut && IntType == 0) hHCRrsoCutI0->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_cut && IntType == 1) hCKRrsoCutI1->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut && IntType == 1) hCNRrsoCutI1->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut && IntType == 1) hHCRrsoCutI1->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_cut && IntType == 2) hCKRrsoCutI2->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut && IntType == 2) hCNRrsoCutI2->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut && IntType == 2) hHCRrsoCutI2->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        if (ck_cut && IntType == 3) hCKRrsoCutI3->fillH2D(mc_mom, bincen * (ck_irig - mc_irig));
        if (cn_cut && IntType == 3) hCNRrsoCutI3->fillH2D(mc_mom, bincen * (cn_irig - mc_irig));
        if (hc_cut && IntType == 3) hHCRrsoCutI3->fillH2D(mc_mom, bincen * (hc_irig - mc_irig));
        
/*        
        if (track.status[0][patt]) hCKMrso->fill(mc_mom, ck_mass);
        if (track.status[1][patt]) hCNMrso->fill(mc_mom, cn_mass);
        if (                 succ) hHCMrso->fill(mc_mom, hc_mass);
        
        if (track.status[0][patt] && IntType == 0) hCKMrsoI0->fill(mc_mom, ck_mass);
        if (track.status[1][patt] && IntType == 0) hCNMrsoI0->fill(mc_mom, cn_mass);
        if (                 succ && IntType == 0) hHCMrsoI0->fill(mc_mom, hc_mass);
        if (track.status[0][patt] && IntType == 1) hCKMrsoI1->fill(mc_mom, ck_mass);
        if (track.status[1][patt] && IntType == 1) hCNMrsoI1->fill(mc_mom, cn_mass);
        if (                 succ && IntType == 1) hHCMrsoI1->fill(mc_mom, hc_mass);
        if (track.status[0][patt] && IntType == 2) hCKMrsoI2->fill(mc_mom, ck_mass);
        if (track.status[1][patt] && IntType == 2) hCNMrsoI2->fill(mc_mom, cn_mass);
        if (                 succ && IntType == 2) hHCMrsoI2->fill(mc_mom, hc_mass);
        if (track.status[0][patt] && IntType == 3) hCKMrsoI3->fill(mc_mom, ck_mass);
        if (track.status[1][patt] && IntType == 3) hCNMrsoI3->fill(mc_mom, cn_mass);
        if (                 succ && IntType == 3) hHCMrsoI3->fill(mc_mom, hc_mass);
        
        if (track.status[0][patt] && fRich->status && fRich->kindOfRad == 0) hCKMRrso->fill(mc_mom, ck_massr);
        if (track.status[1][patt] && fRich->status && fRich->kindOfRad == 0) hCNMRrso->fill(mc_mom, cn_massr);
        if (                 succ && fRich->status && fRich->kindOfRad == 0) hHCMRrso->fill(mc_mom, hc_massr);
        
        if (track.status[0][patt] && fRich->status && fRich->kindOfRad == 0 && IntType == 0) hCKMRrsoI0->fill(mc_mom, ck_massr);
        if (track.status[1][patt] && fRich->status && fRich->kindOfRad == 0 && IntType == 0) hCNMRrsoI0->fill(mc_mom, cn_massr);
        if (                 succ && fRich->status && fRich->kindOfRad == 0 && IntType == 0) hHCMRrsoI0->fill(mc_mom, hc_massr);
        if (track.status[0][patt] && fRich->status && fRich->kindOfRad == 0 && IntType == 1) hCKMRrsoI1->fill(mc_mom, ck_massr);
        if (track.status[1][patt] && fRich->status && fRich->kindOfRad == 0 && IntType == 1) hCNMRrsoI1->fill(mc_mom, cn_massr);
        if (                 succ && fRich->status && fRich->kindOfRad == 0 && IntType == 1) hHCMRrsoI1->fill(mc_mom, hc_massr);
        if (track.status[0][patt] && fRich->status && fRich->kindOfRad == 0 && IntType == 2) hCKMRrsoI2->fill(mc_mom, ck_massr);
        if (track.status[1][patt] && fRich->status && fRich->kindOfRad == 0 && IntType == 2) hCNMRrsoI2->fill(mc_mom, cn_massr);
        if (                 succ && fRich->status && fRich->kindOfRad == 0 && IntType == 2) hHCMRrsoI2->fill(mc_mom, hc_massr);
        if (track.status[0][patt] && fRich->status && fRich->kindOfRad == 0 && IntType == 3) hCKMRrsoI3->fill(mc_mom, ck_massr);
        if (track.status[1][patt] && fRich->status && fRich->kindOfRad == 0 && IntType == 3) hCNMRrsoI3->fill(mc_mom, cn_massr);
        if (                 succ && fRich->status && fRich->kindOfRad == 0 && IntType == 3) hHCMRrsoI3->fill(mc_mom, hc_massr);
*/
        Double_t pow27 = std::pow(10., 1.7) * std::pow(mc_mom, -1.7);
        if (mc_mom > 10. && track.status[0][patt]) hCKflux->fillH1D(ck_irig, pow27);
        if (mc_mom > 10. && track.status[1][patt]) hCNflux->fillH1D(cn_irig, pow27);
        if (mc_mom > 10. &&                  succ) hHCflux->fillH1D(hc_irig, pow27);
        
        if (mc_mom > 10. && ck_cut) hCKflux2->fillH1D(ck_irig, pow27);
        if (mc_mom > 10. && cn_cut) hCNflux2->fillH1D(cn_irig, pow27);
        if (mc_mom > 10. && hc_cut) hHCflux2->fillH1D(hc_irig, pow27);
        
        if (mc_mom > 10. && mc_mom < 20.) {
            if (track.status[0][patt]) hCKRrso2->fillH1D(10. * (ck_irig - mc_irig));
            if (track.status[1][patt]) hCNRrso2->fillH1D(10. * (cn_irig - mc_irig));
            if (                 succ) hHCRrso2->fillH1D(10. * (hc_irig - mc_irig));
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
