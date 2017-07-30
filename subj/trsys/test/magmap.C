void magmap() {
    TH1::AddDirectory(true);
    using namespace MGROOT;
    using namespace TrackSys;

    std::function<double()> unf = MGRndm::Uniform<double>(-50, 50);
    TH1D * dif = new TH1D("dif", "dif", 100, -0.001, 0.001);
    for (Int_t i = 0; i < 10000; ++i) {
        SVecD<3> coo(unf(), unf(), unf());
        MagFld&& ams = MagGeoBoxAms::Get(coo);
        MagFld&& slf = MagMgnt::Get(coo);
        dif->Fill(ams.x() - slf.x());
    }

    Long64_t ntimes = 10000000;

    SVecD<3> coo(unf(), unf(), unf());
    MGClock::HrsStopwatch sw1;
    for (Long64_t it = 0; it < ntimes; ++it) {
        MagFld&& ams = MagGeoBoxAms::Get(coo);
    }
    sw1.stop();
    
    MGClock::HrsStopwatch sw2;
    for (Long64_t it = 0; it < ntimes; ++it) {
        MagFld&& slf = MagMgnt::Get(coo);
    }
    sw2.stop();
    
    sw1.print();
    sw2.print();

    dif->Draw("hist");
}
