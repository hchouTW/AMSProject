// STL c++
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdarg>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>

// ROOT
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

// User defination library
#include "ClassDefRun.h"

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <CPPLibs.h>
#include <ROOTLibs.h>
#include <TrSys.h>

static constexpr double Mproton = 0.938272297;
static constexpr double Mdeuterium = 1.876123915;

int main(int argc, char** argv) {
    if (argc != 2) return 0;

    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    
    Axis AXrig("Rigidity [GV]", 400*10, 0.5, 50.0, AxisScale::kLog);
    Axis AXbta("Velocity", 701*10, 0.3, 1.001);
    
    Hist* hRIG_expt = Hist::New("hRIG_expt", HistAxis(AXrig));
    Hist* hBTA_expt = Hist::New("hBTA_expt", HistAxis(AXbta));
    
    // TOF
    Axis AXTFbta("Velocity", 30, 0.50, 0.80);
    std::vector<double> AXTFrig_bins;
    std::vector<double> AXTFken_bins;
    for (int ii = 0; ii <= AXTFbta.nbin(); ++ii) {
        double ibta = (1.0 / AXTFbta(ii));
        double rig_fact = 1.0 / std::sqrt((ibta - 1.0) * (ibta + 1.0));
        double rig = rig_fact * Mdeuterium;
        AXTFrig_bins.push_back(rig);
        
        double ken_fact = (std::hypot(1.0, rig_fact) - 1.0);
        double ken = ken_fact * Mproton;
        AXTFken_bins.push_back(ken);
    }
    Axis AXTFrig("Rigidity [GV]", AXTFrig_bins);
    Axis AXTFken("Kinetic Energy per Nucleon [GeV/n]", AXTFken_bins);
    
    Hist* hTFrig_expt = Hist::New("hTFrig_expt", HistAxis(AXTFrig));
    Hist* hTFken_expt = Hist::New("hTFken_expt", HistAxis(AXTFken));
   
    // RICH
    Axis AXRHbta("Velocity", 20, 0.96, 0.98);
    std::vector<double> AXRHrig_bins;
    std::vector<double> AXRHken_bins;
    for (int ii = 0; ii <= AXRHbta.nbin(); ++ii) {
        double ibta = (1.0 / AXRHbta(ii));
        double rig_fact = 1.0 / std::sqrt((ibta - 1.0) * (ibta + 1.0));
        double rig = rig_fact * Mdeuterium;
        AXRHrig_bins.push_back(rig);
        
        double ken_fact = (std::hypot(1.0, rig_fact) - 1.0);
        double ken = ken_fact * Mproton;
        AXRHken_bins.push_back(ken);
    }
    Axis AXRHrig("Rigidity [GV]", AXRHrig_bins);
    Axis AXRHken("Kinetic Energy per Nucleon [GeV/n]", AXRHken_bins);
    
    Hist* hRHrig_expt = Hist::New("hRHrig_expt", HistAxis(AXRHrig));
    Hist* hRHken_expt = Hist::New("hRHken_expt", HistAxis(AXRHken));

    // Read RUN 
    std::vector<std::string>&& filelist = MGIO::ReadFileContent("lst/flist.cern.mrun.iss.pass7");
    int fidx = std::atoi(argv[1]);
    int fsat = fidx * 100;
    int fend = (fidx+1) * 100;
    if (fsat > filelist.size()) fsat = filelist.size();
    if (fend > filelist.size()) fend = filelist.size();
    if (fsat == filelist.size()) return -1;

    TChain* mdst = new TChain("mdst");
    for (int idx = fsat; idx < fend; ++idx) mdst->Add(filelist.at(idx).c_str());

    LIST* list = new LIST;
    RTI*  rti  = new RTI;

    //mdst->SetBranchAddress("list", &list);
    mdst->SetBranchAddress("rti" , &rti);

    constexpr double stable_factor = 1.2;
    unsigned long nentries = mdst->GetEntries();
    for (unsigned long entry = 0; entry < nentries; ++entry) {
        mdst->GetEntry(entry);

        if (entry%864000==0) std::cerr << Form("PROCESS (%5.1f)  %ld/%ld\n", (100.0*entry)/nentries, entry, nentries);

        if (!rti->good) continue;
        if (rti->is_in_SAA) continue;
        if (rti->zenith > 40.0) continue;
        if (rti->livetime < 0.5) continue;

        if (rti->tk_align[0][0] > 35.0) continue;
        if (rti->tk_align[0][1] > 35.0) continue;
        if (rti->tk_align[1][0] > 45.0) continue;
        if (rti->tk_align[1][1] > 45.0) continue;
        
        double cfrig = stable_factor * rti->max_IGRF;
        
        int RIGcfbin = AXrig.find(cfrig) + 1;
        for (int ir = RIGcfbin; ir <= AXrig.nbin(); ++ir) {
            double rig = AXrig.center(ir, AxisScale::kLog);
            hRIG_expt->fillH1D(rig, rti->livetime);
        }
        
        double cfbta = 1.0 / std::hypot(1.0, (2.0*Mproton) / (stable_factor * rti->max_IGRF));
        
        int BTAcfbin = AXbta.find(cfbta) + 1;
        for (int ib = BTAcfbin; ib <= AXbta.nbin(); ++ib) {
            double bta = AXbta.center(ib);
            hBTA_expt->fillH1D(bta, rti->livetime);
        }

        // Rigidity
        {
            int TFcfbin = AXTFrig.find(cfrig) + 1;
            for (int ik = TFcfbin; ik <= AXTFrig.nbin(); ++ik) {
                double rig = AXTFrig.center(ik, AxisScale::kLog);
                hTFrig_expt->fillH1D(rig, rti->livetime);
            }
            
            int RHcfbin = AXRHrig.find(cfrig) + 1;
            for (int ik = RHcfbin; ik <= AXRHrig.nbin(); ++ik) {
                double rig = AXRHrig.center(ik, AxisScale::kLog);
                hRHrig_expt->fillH1D(rig, rti->livetime);
            }
        }

        // Ken
        double cfken = Mproton * (std::hypot(1.0, (stable_factor * rti->max_IGRF) / (2.0*Mproton)) - 1.0);
        {
            int TFcfbin = AXTFken.find(cfken) + 1;
            for (int ik = TFcfbin; ik <= AXTFken.nbin(); ++ik) {
                double ken = AXTFken.center(ik, AxisScale::kLog);
                hTFken_expt->fillH1D(ken, rti->livetime);
            }
            
            int RHcfbin = AXRHken.find(cfken) + 1;
            for (int ik = RHcfbin; ik <= AXRHken.nbin(); ++ik) {
                double ken = AXRHken.center(ik, AxisScale::kLog);
                hRHken_expt->fillH1D(ken, rti->livetime);
            }
        }
    }

    TFile* file = new TFile(Form("out/lvtme/lvtme%02d.root", fidx), "RECREATE");
    file->cd();

    (*hRIG_expt)()->Write();
    (*hBTA_expt)()->Write();
    (*hTFrig_expt)()->Write();
    (*hTFken_expt)()->Write();
    (*hRHrig_expt)()->Write();
    (*hRHken_expt)()->Write();

    file->Write();
    file->Close();

    return 0;
}
