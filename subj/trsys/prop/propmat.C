void propmat() {
    TH1::AddDirectory(true);
    using namespace MGROOT;
    using namespace TrackSys;

    Double_t zcoo = 0.;
    Double_t step = 200.;
    Double_t mom = 1.0;

    //MatArg marg;
    MatArg marg(true, true);

    PhySt part;
    PhyJb phyJb;

    COUT("NORM\n");
    part.set_mom(mom);
    part.set_state_with_cos(0., 0., 200.);
    part.print();
    PropMgnt::PropToZ(0., part, marg, &phyJb);
    part.print();
    
    COUT("MC\n");
    part.set_mom(mom);
    part.set_state_with_cos(0., 0., 200.);
    PropMgnt::PropToZWithMC(0., part, marg);
    part.print();

    Long64_t ntimes = 10000000 / 1000;
    
    MGClock::HrsStopwatch sw1;
    for (Long64_t it = 0; it < ntimes; ++it) {
        part.set_mom(mom);
        part.set_state_with_cos(0., 0., 200.);
        PropMgnt::PropToZ(zcoo, part, marg);
    }
    sw1.stop();
   
    MGClock::HrsStopwatch sw2;
    for (Long64_t it = 0; it < ntimes; ++it) {
        part.set_mom(mom);
        part.set_state_with_cos(0., 0., 200.);
        PropMgnt::PropToZ(zcoo, part, marg, &phyJb);
    }
    sw2.stop();

    sw1.print();
    sw2.print();

    TFile * outf = new TFile("propmat.root", "RECREATE");
    
    std::function<double()> unfp = MGRndm::Uniform<double>(-30, 30);
    std::function<double()> unfd = MGRndm::Uniform<double>(-1e-4, 1e-4);
    
    TGraph * graph = new TGraph();
    graph->SetNameTitle("eloss", "eloss");
    for (Long64_t it = 1; it < 100; ++it) {
        part.set_mom(it*0.1);
        part.set_state_with_cos(0., 0., 200.);
        PropMgnt::PropToZ(zcoo, part, marg);
        Double_t eloss = 1000. * (0.1*it - part.mom());
        graph->SetPoint(it-1, 0.1*it, eloss);
    }
    
    graph->Write();
    outf->Write();
    outf->Close();
}
