void matget() {
    TH1::AddDirectory(true);
    using namespace MGROOT;
    using namespace TrackSys;

    MatGeoBoxAms::Load();
    
    Long64_t ntimes = 10000000 / 10;
    
    std::function<double()> unf = MGRndm::Uniform<double>(-50, 50);

    SVecD<3> vcoo(12, 0, 140);
    SVecD<3> wcoo(0, 12,  86);
    
    MatFld&& mattest1 = MatGeoBoxAms::Get(vcoo, wcoo, false);
    MatFld&& mattest2 = MatGeoBoxAms::Get(vcoo, wcoo);
    mattest1.print();
    mattest2.print();
    
    MGClock::HrsStopwatch sw1;
    for (Long64_t it = 0; it < ntimes; ++it) {
        MatFld&& mat = MatGeoBoxAms::Get(vcoo, wcoo, false);
    }
    sw1.stop();
    
    MGClock::HrsStopwatch sw2;
    for (Long64_t it = 0; it < ntimes; ++it) {
        MatFld&& mat = MatGeoBoxAms::Get(vcoo, wcoo);
    }
    sw2.stop();
    
    sw1.print();
    sw2.print();
    
    TFile * outf = new TFile("matget.root", "RECREATE");
   
    TH1D * dif = new TH1D("dif", "dif", 100, -0.08, 0.08);
    for (Long64_t it = 0; it < ntimes; ++it) {
        SVecD<3> ivcoo(unf(), unf(), 140);
        SVecD<3> iwcoo(unf(), unf(),  86);
        
        MatFld&& matfst = MatGeoBoxAms::Get(ivcoo, iwcoo, false);
        MatFld&& matstd = MatGeoBoxAms::Get(ivcoo, iwcoo);

        Double_t sum = (matfst.num_rad_len() + matstd.num_rad_len());
        Double_t sub = (matfst.num_rad_len() - matstd.num_rad_len());
        dif->Fill(sub/sum);
    }

    outf->Write();
    outf->Close();
}
