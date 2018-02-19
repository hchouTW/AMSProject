void Unfold(Hist * hrat, Hist * horg, RootFile * fPr, RootFile * fAp) {
	std::cerr << "\n========  Unfold  ========\n";
	std::cerr << "== Load Tree\n";
	fPr->LoadTree();
	fAp->LoadTree();

	// Proton Flux
	TGraphAsymmErrors * PrFlux = (TGraphAsymmErrors*)(TFile::Open("/afs/cern.ch/work/h/hchou/public/DATABASE/physics/database_antipp/database_pr.root")->Get("gr_exp2"));

	// Antiproton-to-Proton Flux Ratio
	TGraphAsymmErrors * AppRatio = new TGraphAsymmErrors((*hrat)());

	// Unfolding Ratio
	Hist * unfoldRat = Hist::New("unfoldRat", "", AppRlt::AXnr);
	for (Int_t it = 1; it <= AppRlt::AXnr.nbin(); ++it) {
		(*unfoldRat)()->SetBinContent(it, 1);
	}
		
	//---- Proton ----//
	Hist * hPrTue = Hist::New("hPrTue", "", AppRlt::AXnr); 
	Hist * hPrRec = Hist::New("hPrRec", "", AppRlt::AXnr);

	std::cerr << "== Fill PR Event\n";
	for (Long64_t ev = 0; ev < fPr->fTree->GetEntries(); ++ev) {
		fPr->fTree->GetEntry(ev);
		Double_t mcrig = fPr->fMCRig;
		Double_t trrig = fPr->fTRRig;
		if (!fPr->fOptL) continue;
		if (trrig < 0) continue;
		Double_t prwgt = std::pow(mcrig, -1.7) * PrFlux->Eval(mcrig);
		hPrTue->fill(mcrig, prwgt);
		hPrRec->fill(trrig, prwgt);
	}

	//---- Iteration ----//
	const Int_t maxIter = 3;
	Int_t iter = 1;
	while (iter <= maxIter) {
		std::cerr << CStrFmt("== Iter (%02d/%02d)\n", iter, maxIter);
		
		//---- Antiproton ----//
		Hist * hApTue = Hist::New(CStrFmt("hApTue%02d", iter), "", AppRlt::AXnr); 
		Hist * hApRec = Hist::New(CStrFmt("hApRec%02d", iter), "", AppRlt::AXnr); 
		std::cerr << "== Fill Ap Event\n";
		for (Long64_t ev = 0; ev < fAp->fTree->GetEntries(); ++ev) {
			fAp->fTree->GetEntry(ev);
			Double_t mcrig = -fAp->fMCRig;
			Double_t trrig = -fAp->fTRRig;
			if (!fAp->fOptL) continue;
			if (trrig < 0) continue;
			Double_t prwgt = std::pow(mcrig, -1.7) * PrFlux->Eval(mcrig);
			Double_t apwgt = prwgt * AppRatio->Eval(mcrig) * (*unfoldRat)()->Interpolate(mcrig);
			hApTue->fill(mcrig, apwgt);
			hApRec->fill(trrig, apwgt);
		}

		//---- Flux Ratio ----//
		Hist * hAppTue = Hist::New(StrFmt("hAppTue%02d", iter), "", AppRlt::AXnr); 
		Hist * hAppRec = Hist::New(StrFmt("hAppRec%02d", iter), "", AppRlt::AXnr);
		Hist * hAppWgt = Hist::New(StrFmt("hAppWgt%02d", iter), "", AppRlt::AXnr);
		for (Int_t it = 1; it <= AppRlt::AXnr.nbin(); ++it) {
			(*hAppTue)()->SetBinContent(it, 
					(*hApTue)()->GetBinContent(it) / (*hPrTue)()->GetBinContent(it));
			(*hAppRec)()->SetBinContent(it, 
					(*hApRec)()->GetBinContent(it) / (*hPrRec)()->GetBinContent(it));
			(*hAppWgt)()->SetBinContent(it, 
					(*horg)()->GetBinContent(it) / (*hAppRec)()->GetBinContent(it));
			(*hAppWgt)()->SetBinError(it, 0.05*(*hAppWgt)()->GetBinContent(it));
		}
		//TF1 * fAppWgt = new TF1(CStrFmt("fAppWgt%02d", iter), "1.0+[0]*TMath::Power([1]*(x+[2]), [3])", 1, 800);
		//fAppWgt->SetParameters(3.56281e+09, 7.92588e-03, 1.43510e+02, -1.77527e+02);
		//(*hAppWgt)()->Fit(fAppWgt, "q0", "", 1, 30);
		//(*hAppWgt)()->Fit(fAppWgt, "q0", "", 1, 30);

		for (Int_t it = 1; it <= AppRlt::AXnr.nbin(); ++it) {
			//(*unfoldRat)()->SetBinContent(it, 
			//		(*unfoldRat)()->GetBinContent(it) * fAppWgt->Eval((*hAppTue)()->GetBinCenter(it)));
			(*unfoldRat)()->SetBinContent(it, 
					(*unfoldRat)()->GetBinContent(it) * (*hAppWgt)()->GetBinContent(it));
		}
		
		iter++;
	}
}
