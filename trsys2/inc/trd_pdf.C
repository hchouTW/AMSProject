double GausLogGausFunc(double *xx, double *param) {
	double logx = xx[0];
	double x    = std::exp(logx);
	
	double wgt1 = param[0];
	double men  = param[1];
	double sgm  = param[2];
	
	double wgt2 = param[3];
	double lnm  = param[4];
	double lns  = param[5];

	if(sgm <= 0.0 || lns <= 0.0) return 0.0;

	double gs_norm = (x - men) / sgm;
	double func1   = std::exp(-0.5 * gs_norm * gs_norm) / (sgm * std::sqrt(2.0 * TMath::Pi())); 

	double ln_norm = (logx - lnm) / lns;
	double func2   = std::exp(-0.5 * ln_norm * ln_norm) / (x * lns * std::sqrt(2.0 * TMath::Pi()));

	double func = wgt1 * func1 + wgt2 * func2;
	return func;
}

TF1* load_func() {
	static TF1* gaus_ln_gaus_func = new TF1("gaus_ln_gaus_func", GausLogGausFunc, -3., 5., 6);
	gaus_ln_gaus_func->SetParameters(1.42869e+05, 2.78050e+00, 1.23789e+00, 5.55351e+05, 2.53627e+00, 5.95098e-01);
	return gaus_ln_gaus_func;
}
