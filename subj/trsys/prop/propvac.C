void propvac() {
    TH1::AddDirectory(true);
    using namespace MGROOT;
    using namespace TrackSys;

    Double_t zcoo = 0.;
    Double_t step = 200.;

    PhySt part;
    part.set_mom(10.);
    part.set_state_with_cos(0., 0., 200.);

    COUT("AMS Libs\n");
    part.print();
    PropMgnt::PropToZ_AMSLibs(zcoo, part);
    part.print();
    
    COUT("Self PropToZ\n");
    part.set_mom(10.);
    part.set_state_with_cos(0., 0., 200.);
    part.print();
    PropMgnt::PropToZ(zcoo, part);
    part.print();

    Long64_t ntimes = 10000000 / 100;
   
    MGClock::HrsStopwatch sw1;
    for (Long64_t it = 0; it < ntimes; ++it) {
        part.set_mom(10.);
        part.set_state_with_cos(0., 0., 200.);
        PropMgnt::PropToZ_AMSLibs(zcoo, part);
    }
    sw1.stop();
    
    MGClock::HrsStopwatch sw2;
    for (Long64_t it = 0; it < ntimes; ++it) {
        part.set_mom(10.);
        part.set_state_with_cos(0., 0., 200.);
        PropMgnt::PropToZ(zcoo, part);
    }
    sw2.stop();
    
    MGClock::HrsStopwatch sw3;
    for (Long64_t it = 0; it < ntimes; ++it) {
        part.set_mom(10.);
        part.set_state_with_cos(0., 0., 200.);
        PropMgnt::Prop(step, part);
    }
    sw3.stop();

    sw1.print();
    sw2.print();
    sw3.print();
    
    TFile * outf = new TFile("propvac.root", "RECREATE");
    
    std::function<double()> unfp = MGRndm::Uniform<double>(-30, 30);
    std::function<double()> unfd = MGRndm::Uniform<double>(-1e-4, 1e-4);
    TH1D * dif = new TH1D("dif", "dif", 100, -0.0004, 0.0004);

    for (Long64_t it = 0; it < ntimes; ++it) {
        Double_t x = unfp();
        Double_t y = unfp();
        Double_t dx = unfd();
        Double_t dy = unfd();
        
        part.set_mom(10.);
        part.set_state_with_cos(x, y, 200., dx, dy);
        PropMgnt::PropToZ_AMSLibs(zcoo, part);
        Double_t y1 = part.cy();
        
        part.set_mom(10.);
        part.set_state_with_cos(x, y, 200., dx, dy);
        PropMgnt::PropToZ(zcoo, part);
        Double_t y2 = part.cy();

        dif->Fill(y2 - y1);
    }
    
    outf->Write();
    outf->Close();
}
