double landau(double norm) {
	double kLandau0   =  1.80655634e-01;
	double kLandau0_X = -2.22782980e-01;
	double kWidthScl  =  1.17741002e+00; // sqrt(2*log(2))
	double kLandauArgs[4] = { 4.90120e-01, 1.16366e+00, -3.99384e+00, 1.27500e-01 }; // Fit landau with "nrm" from -3 to 37
	
	double neg_ln = kLandauArgs[0] * (kLandauArgs[1] * norm + std::exp(-kLandauArgs[1] * norm) - 1.0) +
	                kLandauArgs[2] * (kLandauArgs[3] * norm + std::exp(-kLandauArgs[3] * norm) - 1.0);
	return neg_ln;
}

double LandauGausFunc(double *xx, double *param) {
	const int kN = 41;                                                                                                                                                                              
	const double kX[kN] = {
	        -3.00, -2.85, -2.70, -2.55, -2.40, -2.25, -2.10, -1.95, -1.80, -1.65, 
	        -1.50, -1.35, -1.20, -1.05, -0.90, -0.75, -0.60, -0.45, -0.30, -0.15, 
	         0.00,  0.15,  0.30,  0.45,  0.60,  0.75,  0.90,  1.05,  1.20,  1.35, 
	         1.50,  1.65,  1.80,  1.95,  2.10,  2.25,  2.40,  2.55,  2.70,  2.85, 3.00 };
	const double kP[kN] = { 
	        0.000666, 0.001033, 0.001566, 0.002322, 0.003366, 0.004771, 0.006611, 0.008958, 0.011867, 0.015372, 
	        0.019468, 0.024108, 0.029189, 0.034554, 0.039996, 0.045265, 0.050088, 0.054192, 0.057328, 0.059296, 
	        0.059966, 0.059296, 0.057328, 0.054192, 0.050088, 0.045265, 0.039996, 0.034554, 0.029189, 0.024108, 
	        0.019468, 0.015372, 0.011867, 0.008958, 0.006611, 0.004771, 0.003366, 0.002322, 0.001566, 0.001033, 0.000666 };
	
	double x    = xx[0];
	double wgt  = param[0];
	double kpa  = param[1];
	double mpv  = param[2];
	double sgm  = param[3];
	double fluc = param[4];

	if(sgm <= 0.0 || kpa < 0 || kpa > 1.0 || fluc < 0.0) return 0.0;
   	bool is_fluc = (fluc > 0.0);
    
	double value = 0.0;
	double lg_norm = (x - mpv) / sgm;
	if (!is_fluc) {
	    double gs = -0.5 * lg_norm * lg_norm;
	    double ld = -landau(lg_norm);
	    double lg = kpa * gs + (1.0 - kpa) * ld; 
	    value = std::exp(lg);
	}   
	else {
	    for (int it = 0; it < kN; ++it) {
	        double xx = lg_norm + (fluc / sgm) * kX[it];
	        double gs = -0.5 * xx * xx; 
	        double ld = -landau(xx);
	        double lg = kpa * gs + (1.0 - kpa) * ld; 
	        value += std::exp(lg) * kP[it];
	    }
	}   
	if (!std::isfinite(value)) value = 0.0;
	return wgt * value;
}

void fit_edep() {
	TF1* func = new TF1("func", LandauGausFunc, 0, 100, 5);
	TF1* ffit = new TF1("ffit", "[0] + [1] * (1+x*x) - [2]*log([3] + (x*x))", 0.0001, 10);
	
	TF1* kfit = new TF1("kfit", "0.5*(1+TMath::Erf([0]*log(x*x)-[1]))+[2]*TMath::Erfc([3]*log(x*x)+[4])", 0.0001, 10);
	kfit->SetParameters(5.88012e-01, 2.39019e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00); // kpv
	
	TFile* file = TFile::Open("YiMdst_HE.root");
	TH2D* hist = (TH2D*)file->Get("hTd");

	//TH1D* htime = new TH1D("htime", "", hist->GetXaxis()->GetNbins(), hist->GetXaxis()->GetXbins()->GetArray());
	TGraphErrors* gkpa  = new TGraphErrors();
	TGraphErrors* gmpv  = new TGraphErrors();
	TGraphErrors* gsgm  = new TGraphErrors();
	TGraphErrors* gfluc = new TGraphErrors();
	for (int ti = 30; ti <= hist->GetXaxis()->GetNbins()-1; ++ti) {
		double igb = hist->GetXaxis()->GetBinCenter(ti);
		//kfit->SetParameters(5.88012e-01, 2.39019e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00); // pr kpv
		kfit->SetParameters(5.37505e-01, 2.10937e+00, 2.80664e-03, 0.00000e+00, 0.00000e+00); // he kpv
		double kpa = kfit->Eval(igb);
		//ffit->SetParameters(-3.52442e-01, 9.67636e-01, 8.30227e-02, 1.15207e-05); // pr mpv
		ffit->SetParameters(-3.13533e-01, 4.78310e+00, 4.10680e-01, 1.15207e-05); // he mpv
		double mpv = ffit->Eval(igb);
		//ffit->SetParameters(-4.06982e-02, 3.69372e-01, 3.06650e-02, 1.15207e-05); // pr sgm
		ffit->SetParameters(4.21523e-01, 1.03661e+00, 5.44214e-02, 1.15207e-05); // he sgm
		double sgm = ffit->Eval(igb);
		func->SetParameters(1000.0, 0.05, mpv, sgm, 0.001);
		func->FixParameter(1, kpa);
		func->FixParameter(2, mpv);
		func->FixParameter(3, sgm);
		func->FixParameter(4, 0.000);
		TH1D* proj = (TH1D*)(hist->ProjectionY(Form("Proj%d", ti), ti, ti));
		proj->Fit(func, "", "q0", 0.7 * mpv, 100);
		//htime->SetBinContent(ti, func->GetParameter(1));
		//if (ti == 1 || ti == 45 || ti == 46) continue;
		gkpa->SetPoint(gkpa->GetN(), hist->GetXaxis()->GetBinCenter(ti), func->GetParameter(1));
		gkpa->SetPointError(gkpa->GetN()-1, 0, 0.001);
		gmpv->SetPoint(gmpv->GetN(), hist->GetXaxis()->GetBinCenter(ti), func->GetParameter(2));
		gmpv->SetPointError(gmpv->GetN()-1, 0, 0.001);
		gsgm->SetPoint(gsgm->GetN(), hist->GetXaxis()->GetBinCenter(ti), func->GetParameter(3));
		gsgm->SetPointError(gsgm->GetN()-1, 0, 0.001);
		gfluc->SetPoint(gfluc->GetN(), hist->GetXaxis()->GetBinCenter(ti), func->GetParameter(4));
		gfluc->SetPointError(gfluc->GetN()-1, 0, 0.001);
		TCanvas* cv = new TCanvas;
		TH1D* pp = (TH1D*)proj->Clone();
		TF1* ff = (TF1*)func->Clone();
		pp->Draw("hist");
		ff->Draw("same");
		if (ti == 30) cv->Print("td.pdf(");
		else cv->Print("td.pdf");
	}
	//htime->Draw("hist");
	TCanvas* cvs = new TCanvas;
	gkpa->Draw("ape");
	gkpa->SetLineColor(kRed);
	cvs->SaveAs("td.pdf");
	cvs = new TCanvas;
	gmpv->Draw("ape");
	gmpv->SetLineColor(kRed);
	cvs->SaveAs("td.pdf");
	cvs = new TCanvas;
	gsgm->Draw("ape");
	gsgm->SetLineColor(kRed);
	cvs->SaveAs("td.pdf");
	cvs = new TCanvas;
	gfluc->Draw("ape");
	gfluc->SetLineColor(kRed);
	cvs->SaveAs("td.pdf)");
}
