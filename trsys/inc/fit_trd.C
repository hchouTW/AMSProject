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
	TF1* func = new TF1("func", LandauGausFunc, 0, 30, 5);
	
	TFile* file = TFile::Open("YiMdst56.root");
	TH2D* hist = (TH2D*)file->Get("hTd");

	//TH1D* htime = new TH1D("htime", "", hist->GetXaxis()->GetNbins(), hist->GetXaxis()->GetXbins()->GetArray());
	TGraphErrors* gkpa = new TGraphErrors();
	TGraphErrors* gmpv = new TGraphErrors();
	TGraphErrors* gsgm = new TGraphErrors();
	for (int ti = 3; ti <= hist->GetXaxis()->GetNbins(); ++ti) {
		func->SetParameters(0.5, 1000.0, 0.1, 1.0, 0.05);
		func->FixParameter(1, 0.1);
		//func->FixParameter(3, 0.3);
		func->FixParameter(4, 0.0);
		TH1D* proj = (TH1D*)((hist->ProjectionY(Form("Proj%d", ti), ti, ti))->Clone());
		proj->Fit(func, "", "q0");
		//htime->SetBinContent(ti, func->GetParameter(1));
		//if (ti == 1 || ti == 45 || ti == 46) continue;
		gkpa->SetPoint(gkpa->GetN(), hist->GetXaxis()->GetBinCenter(ti), func->GetParameter(1));
		gkpa->SetPointError(gkpa->GetN()-1, 0, 0.01);
		gmpv->SetPoint(gmpv->GetN(), hist->GetXaxis()->GetBinCenter(ti), func->GetParameter(2));
		gmpv->SetPointError(gmpv->GetN()-1, 0, 0.01);
		gsgm->SetPoint(gsgm->GetN(), hist->GetXaxis()->GetBinCenter(ti), func->GetParameter(3));
		gsgm->SetPointError(gsgm->GetN()-1, 0, 0.01);
		//new TCanvas;
		//proj->Draw("hist");
		//func->Draw("same");
	}
	//htime->Draw("hist");
	new TCanvas;
	gkpa->Draw("ape");
	gkpa->SetLineColor(kRed);
	new TCanvas;
	gmpv->Draw("ape");
	gmpv->SetLineColor(kRed);
	new TCanvas;
	gsgm->Draw("ape");
	gsgm->SetLineColor(kRed);
}
