#include <CPPLibs/CPPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>
int main(int argc, char * argv[]) {
    //TH1::AddDirectory(true);
    using namespace MGROOT;
    using namespace TrackSys;
    
    TChain * chain = new TChain("Ntuple");
    
    //chain->Add("/data3/hchou/AMSData/Lorentz/proton_0.3GeV_990k.root");
    //chain->Add("/data3/hchou/AMSData/Lorentz/proton_0.5GeV_1000k.root");
    //chain->Add("/data3/hchou/AMSData/Lorentz/proton_1GeV_1000k.root");
    //chain->Add("/data3/hchou/AMSData/Lorentz/proton_10GeV_1000k.root");
    //chain->Add("/data3/hchou/AMSData/Lorentz/proton_100GeV_1000k.root");
    chain->Add("/data3/hchou/AMSData/Lorentz/e-_0.3GeV_1000k.root");
    //chain->Add("/data3/hchou/AMSData/Lorentz/e-_100GeV_1000k.root");

    //MGConfig::JobOpt opt(argc, argv);
    //for (auto&& str : opt.flist()) chain->Add(str.c_str());
  
    std::vector<double> * pos_x = nullptr;
    std::vector<double> * pos_y = nullptr;
    std::vector<double> * pos_z = nullptr;
    std::vector<double> * mom_x = nullptr;
    std::vector<double> * mom_y = nullptr;
    std::vector<double> * mom_z = nullptr;

    chain->SetBranchAddress("pos_x", &pos_x);
    chain->SetBranchAddress("pos_y", &pos_y);
    chain->SetBranchAddress("pos_z", &pos_z);
    chain->SetBranchAddress("mom_x", &mom_x);
    chain->SetBranchAddress("mom_y", &mom_y);
    chain->SetBranchAddress("mom_z", &mom_z);

    TH1D * hx = new TH1D("hx", "hx", 400, -0.1, 0.1);
    TH1D * hy = new TH1D("hy", "hy", 400, -0.1, 0.1);
    TH1D * dx = new TH1D("dx", "dx", 400, -0.1, 0.1);
    TH1D * dy = new TH1D("dy", "dy", 400, -0.1, 0.1);
    TH1D * mm = new TH1D("mm", "mm", 400,  0.0, 0.001);
    TH1D * mr = new TH1D("mr", "mr", 400,  0.0, 1.0);


    for (Long64_t it = 0; it < chain->GetEntries(); ++it) {
        //if (it > 1) break;
        chain->GetEntry(it++);
        std::vector<double> px = *pos_x;
        std::vector<double> py = *pos_y;
        std::vector<double> pz = *pos_z;
        chain->GetEntry(it);
        std::vector<double> mx = *mom_x;
        std::vector<double> my = *mom_y;
        std::vector<double> mz = *mom_z;

        //COUT("Entry %8ld Vector-Size POS %d %d %d MOM %d %d %d (POSZ, MOMZ) (%8.2f, %8.2f)\n",
        //    it,
        //    pos_x->size(), pos_y->size(), pos_z->size(),
        //    mom_x->size(), mom_y->size(), mom_z->size(),
        //    (pos_y->size() ? pos_y->at(0) : 0), (mom_y->size() ? mom_y->at(0) : 0)
        //);
        if (px.size() < 2) continue;
        hx->Fill(0.1*px.at(1));
        hy->Fill(0.1*py.at(1));
       
        //PhySt part(PartType::Proton);
        PhySt part(PartType::Electron);
        part.set_state(0.1*px.at(0), 0.1*py.at(0), 0.1*pz.at(0), 0.001*mx.at(0), 0.001*my.at(0), 0.001*mz.at(0));
        PropMgnt::PropToZ(0.1*pz.at(1), part);

        dx->Fill(part.cx()-0.1*px.at(1));
        dy->Fill(part.cy()-0.1*py.at(1));
        
        double mm0 = 0.001 * std::sqrt(mx.at(0) * mx.at(0) + my.at(0) * my.at(0) + mz.at(0) * mz.at(0));
        double mm1 = 0.001 * std::sqrt(mx.at(1) * mx.at(1) + my.at(1) * my.at(1) + mz.at(1) * mz.at(1));
        mm->Fill(mm0 - mm1);
        mr->Fill((mm0 - mm1)/mm0);

    }
    

    

    //for (Long64_t it = 0; it < 100000; ++it) {
    //    part.set_mom(10.);
    //    part.set_state_with_cos(0., 0., 50.);
    //    PropMgnt::PropToZ(0., part);
    //    tx->Fill(part.cx());
    //    ty->Fill(part.cy());
    //}

    TFile * file = new TFile("out.root", "RECREATE");
    hx->Write();
    hy->Write();
    dx->Write();
    dy->Write();
    mm->Write();
    mr->Write();
    file->Write();
    file->Close();

    return 0;
}
