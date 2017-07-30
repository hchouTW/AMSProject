void fitvac() {
    TH1::AddDirectory(true);
    using namespace MGROOT;
    using namespace TrackSys;

    Double_t mom = 1.0;

    const Int_t nlay = 9;
    const Double_t zcoo[nlay] = { 150., 50., 30., 25., 2., -3., -25., -30., -130. };

    //std::function<double()> unfx = MGRndm::Uniform<double>(20.0e-4, 40.0e-4);
    //std::function<double()> unfy = MGRndm::Uniform<double>(5.0e-4, 20.0e-4);
  
    const Int_t ntimes = 100000*20;

    TFile * outf = new TFile("fitvac.root", "RECREATE");
    
    TH1D * resx = new TH1D("resx", "resx", 800, -0.04, 0.04);
    TH1D * resy = new TH1D("resy", "resy", 800, -0.02, 0.02);
    TH1D * dif1 = new TH1D("dif1", "dif1", 800, -5, 5);
    TH1D * dif2 = new TH1D("dif2", "dif2", 800, -5, 5);
    TH1D * dif3 = new TH1D("dif3", "dif3", 800, -5, 5);

    const Double_t wgt1 = 0.9;
    const Double_t sgm1 = 1.0;
    const Double_t wgt2 = 0.1;
    const Double_t sgm2 = 3.0;
    TF1 * func = new TF1("func", "[0]/[1]*exp(-0.5*x*x/[1]/[1])+[2]/[3]*exp(-0.5*x*x/[3]/[3])", -30, 30);
    func->SetNpx(100000);
    func->SetParameters(wgt1, sgm1, wgt2, sgm2);


    for (Int_t it = 0; it < ntimes; ++it) {
        if (it%10000==0) COUT("TIME %d/%d\n", it, ntimes);

        PhySt part;
        part.set_mom(mom);
        part.set_state_with_cos(0, 0, 150.);
        
        //TrFit trFit;
        std::vector<HitSt> hits1;
        std::vector<HitSt> hits2;
        for (Int_t ih = 0; ih < nlay; ++ih) {
            //PropMgnt::PropToZ_AMSLibs(zcoo[ih], part);
            PropMgnt::PropToZWithMC(zcoo[ih], part, MatArg(true, false));
            //PropMgnt::PropToZWithMC(zcoo[ih], part, MatArg(false, true));
            
            //Double_t res[2] = { hit1.ex() * MGRndm::NormalGaussian(), hit1.ey() * MGRndm::NormalGaussian() };
            Double_t res[2] = { HitSt::GetDefaultErrX(0.0) * func->GetRandom(), HitSt::GetDefaultErrY(0.0) * func->GetRandom() };
            
            Double_t x = part.cx() + res[0];
            Double_t y = part.cy() + res[1];
            Double_t z = part.cz();
            
            HitSt hit1;
            hit1.set_coo(x, y, z);
            hit1.set_err(MultiGauss(24.0e-4), MultiGauss(10.0e-4));
            hits1.push_back(hit1);
           
            HitSt hit2;
            hit2.set_coo(x, y, z);
            hit2.set_err(MultiGauss(wgt1, sgm1*24.0e-4, wgt2, sgm2*24.0e-4), MultiGauss(wgt1, sgm1*10.0e-4, wgt2, sgm2*10.0e-4));
            hits2.push_back(hit2);
            
            resx->Fill(res[0]);
            resy->Fill(res[1]);
            //trFit.Add(x, y, z, hit.ex(), hit.ey(), hit.ez());
            //
            COUT("=== ID %d Z %8.2f ===\n", ih, zcoo[ih]);
            part.print();
        }
        //trFit.SimpleFit();

        PhyTr tr1(hits1);
        if (tr1.exist()) dif1->Fill( (tr1.part().irig() - 1.0/mom) * mom);
        
        PhyTr tr2(hits2);
        if (tr2.exist()) dif2->Fill( (tr2.part().irig() - 1.0/mom) * mom);

        if (tr1.exist() && tr2.exist()) dif3->Fill( (tr1.part().irig() - tr2.part().irig()) * mom);
        if (it == 0) break;
    }

    TH1D * rat = new TH1D("rat", "rat", 800, -8, 8);
    for (Int_t it = 1; it <= rat->GetNbinsX(); ++it) {
        if (dif1->GetBinContent(it) > 0.1 && dif2->GetBinContent(it) > 0.1) { 
            rat->SetBinContent(it, dif2->GetBinContent(it)/dif1->GetBinContent(it)); 
            rat->SetBinError(it, rat->GetBinContent(it)*std::sqrt(1./dif1->GetBinContent(it) + 1./dif2->GetBinContent(it))); 
        }
    }

    outf->Write();
    outf->Close();
}
