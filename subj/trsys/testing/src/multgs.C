void multgs() {
    TFile * file = new TFile("multgs.root", "RECREATE");
    TH1D * tgs = new TH1D("tgs", "tgs", 400, -5, 5);
    TH1D * sgs = new TH1D("sgs", "sgs", 400, -5, 5);
    TH1D * mgs = new TH1D("mgs", "mgs", 400, -5, 5);
    TH1D * cgs = new TH1D("cgs", "cgs", 400, -5, 5);
   
    const int ns = 6;
    const double es[ns] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.5 };
    double sums = 0, wgts = 0;
    for (Int_t s = 0; s < ns; ++s) { sums += es[s]; wgts += es[s]*es[s]; }
    wgts = std::sqrt(wgts);

    TRandom3 rndm(0);
    for (Int_t it = 0; it < 1000000; ++it) {
        double tp = wgts * rndm.Gaus();
        double sp = 0;
        double mp = 0;
        double cp = 0;
        double mr = rndm.Gaus();
        for (Int_t s = 0; s < ns; ++s) {
            double rs = es[s] * rndm.Gaus();
            sp += rs;
            mp += es[s] * mr;
            cp += ((es[s] * mr) * (wgts / sums));
        }
        tgs->Fill(tp);
        sgs->Fill(sp);
        mgs->Fill(mp);
        cgs->Fill(cp);
    }

    file->Write();
    file->Close();
}
