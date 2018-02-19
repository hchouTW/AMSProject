#ifndef __Unfold_C__
#define __Unfold_C__

#include <omp.h>
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/ROOT/selflib.h" 
#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/Unfold/Include.h" 

int main() {
	// Pamela Result
	TString pamelaFilePath = "/afs/cern.ch/user/h/hchou/public/DATASET/physics/antipp3/others/antipp_pamela.root";
	TFile * inFile_Pamela = TFile::Open(pamelaFilePath.Data());
	TGraphErrors * hRat_Pamela = (TGraphErrors *) inFile_Pamela->Get("antipp_pamela");
	TGraph * RatPamela = MgntGraph::New("Ratio_Pamela", "Antiproton / Proton Ratio;Absolute Rigidity [GV];Antiproton / Proton Ratio", hRat_Pamela);
	inFile_Pamela->Close();
	
	// AMS Day Result
	TString amsdayFilePath = "/afs/cern.ch/user/h/hchou/public/DATASET/physics/antipp3/others/ams02_acceptance_corrected_20150330.root";
	TFile * inFile_Amsday = TFile::Open(amsdayFilePath.Data());
	TGraphErrors * hRat_AMSDay = (TGraphErrors *) inFile_Amsday->Get("gratio");
	TGraph * RatAMSDay = MgntGraph::New("Ratio_AMSDay", "Antiproton / Proton Ratio;Absolute Rigidity [GV];Antiproton / Proton Ratio", hRat_AMSDay);
	inFile_Amsday->Close();
	
	// AMS Pub Result
	TString amspubFilePath = "/afs/cern.ch/user/h/hchou/public/DATASET/physics/antipp3/others/20160406MIT.root";
	TFile * inFile_Amspub = TFile::Open(amspubFilePath.Data());
	TGraphErrors * hRat_AMSPub = (TGraphErrors *) inFile_Amspub->Get("gpbarp");
	TGraph * RioAMSPub = MgntGraph::New("Ratio_AMSPub", "Antiproton / Proton Ratio;Absolute Rigidity [GV];Antiproton / Proton Ratio", hRat_AMSPub);
	inFile_Amspub->Close();

	// Antipp Binning
	TString sbn = "/afs/cern.ch/user/h/hchou/public/DATABASE/physics/binning/antipp2_rig.root";
	TFile * fBin = TFile::Open(sbn);
	TH1D  * hBn0 = (TH1D*)fBin->Get("hbin0");
	TH1D  * hBn1 = (TH1D*)fBin->Get("hbin1");
	AxisBin AXBnr = AxisBin::AxisBin(hBn0, AxisBin::kX);
  AxisBin AXBir = AxisBin::AxisBin(hBn1, AxisBin::kX);
	fBin->Close();

	Int_t mgFB = 3;
	Int_t mgLB = 42;
	AxisBin axis = AXBnr;
	axis.merge(mgFB, mgLB);

	// Input File	
	std::cout << "**** LOAD : AMS Pr Flux\n";
	TString inAMSFlx = "/afs/cern.ch/user/h/hchou/public/DATABASE/physics/pflux/AMS_proton_flux.root";
	TFile * fAMS = TFile::Open(inAMSFlx.Data());
	TGraph * grAMSFlx = (TGraph *) fAMS->Get("Proton_flux")->Clone("Proton_flux_Copy");
	fAMS->Close();
	ISSFlx * amsFlx = new ISSFlx(axis);
	amsFlx->buildFinFlx(grAMSFlx);
	
	// ISS Event Ap / Pr Ratio
	std::cout << "**** LOAD : ISS Ap/Pr Ratio\n";
	TString inISSDir = "/afs/cern.ch/user/h/hchou/public/DATASET/physics/antipp3/final";
	TFile * fISS = TFile::Open(Form("%s/YiAna.projHist.root", inISSDir.Data()));
	fISS->cd();
	
	TH1D * hRio = (TH1D *) fISS->Get("hALL_Ratio");
	TH1D * haprio = new TH1D("haprio", "haprio", axis.Nbin(), axis.Bins());
	for (Int_t ib = mgFB; ib <= mgLB; ++ib) {
		Int_t ibin = ib - mgFB + 1;
		haprio->SetBinContent(ibin, hRio->GetBinContent(ib));
		haprio->SetBinError  (ibin, hRio->GetBinError  (ib));
	}
	
	ISSRio * issRio = new ISSRio(axis);
	issRio->buildRecRio(haprio);
	issRio->buildFinRio(haprio);

	// Antiproton flux
	TH1D * hapflx = new TH1D("hapflx", "hapflx", axis.Nbin(), axis.Bins());
	for (Int_t ib = 1; ib <= axis.Nbin(); ++ib) {
		Double_t cen = axis.BinCenter(ib, AxisBin::kLog);
		Double_t rer = haprio->GetBinError(ib) / haprio->GetBinContent(ib);
		Double_t val = haprio->GetBinContent(ib) * amsFlx->fFinFlx.eval(cen);
		Double_t err = val * rer;
		hapflx->SetBinContent(ib, val);
		hapflx->SetBinError  (ib, err);
	}

	// MC Trigger
	std::cout << "**** LOAD : MC\n";
	Int_t ChrgPr = 1;
	Int_t ChrgAp = -1;
	Double_t MomLU[2][2] = { { 0.5, 10 }, { 1, 800 } };
	TString inMCDir = "/nas09/hchou/V.2016Apr13v1/ntuple";
	
	TFile * fPr[2] = { 0, 0 };
	fPr[0] = TFile::Open(Form("%s/YiAnalytics_Pr0510.root", inMCDir.Data()));
	fPr[1] = TFile::Open(Form("%s/YiAnalytics_Pr1800.root", inMCDir.Data()));
	
	TFile * fAp[2] = { 0, 0 };
	fAp[0] = TFile::Open(Form("%s/YiAnalytics_Ap0510.root", inMCDir.Data()));
	fAp[1] = TFile::Open(Form("%s/YiAnalytics_Ap1800.root", inMCDir.Data()));

	MCFlx * mcPrFlx[2] = { new MCFlx(axis), new MCFlx(axis) };
	for (Int_t v = 0; v < 2; ++v) {
		fPr[v]->cd();
		UInt_t numTrg = 0;
		UInt_t varTrg = 0;
		TTree * tPr = (TTree *) fPr[v]->Get("DSTRun");
		tPr->SetBranchAddress("numOfTrgEvent", &varTrg);
		for (Long64_t it = 0; it < tPr->GetEntries(); ++it) { tPr->GetEntry(it), numTrg += varTrg; }
		mcPrFlx[v]->buildIntFlx(numTrg, MomLU[v][0], MomLU[v][1], ChrgPr);
		mcPrFlx[v]->buildFinFlx(&amsFlx->fFinFlx);
	}

	MCFlx * mcApFlx[2] = { new MCFlx(axis), new MCFlx(axis) };
	for (Int_t v = 0; v < 2; ++v) {
		fAp[v]->cd();
		UInt_t numTrg = 0;
		UInt_t varTrg = 0;
		TTree * tAp = (TTree *) fAp[v]->Get("DSTRun");
		tAp->SetBranchAddress("numOfTrgEvent", &varTrg);
		for (Long64_t it = 0; it < tAp->GetEntries(); ++it) { tAp->GetEntry(it), numTrg += varTrg; }
		mcApFlx[v]->buildIntFlx(numTrg, MomLU[v][0], MomLU[v][1], ChrgAp);
		mcApFlx[v]->buildFinFlx(hapflx);
	}	

	TTree * mcPrTre[2] = { 0, 0 };
	Float_t mcPrGen[2] = { 0, 0 };
	Float_t mcPrRec[2] = { 0, 0 };

	// For Unfold
	const Int_t NPt = 6;
	const TString PtNm[6] = { "dstSL", "dstLM", "dstMH", "dstHA", "dstHB", "dstHC" };
	const Int_t PtBn[6][2] = { {2, 13}, {14, 31}, {14, 31}, {32, 51}, {32, 51}, {32, 56} };
	const TString GenNm = "McRig";
	const TString RecNm = "Rig";
	Float_t GenVar = 0;
	Float_t RecVar = 0;

	const UInt_t NIt = 8;

	// Level 01 loop -- pattern
	// Level 02 loop -- event
	std::cout << Form("*****************\n");
	std::cout << Form("**  MC Proton  **\n");
	std::cout << Form("*****************\n");
	mcPrFlx[1]->buildMigMtx();
	for (Int_t ipt = 0; ipt < NPt; ++ipt) { // Level 02 loop -- pattern
		if (ipt == 3 || ipt == 4) continue;
		TTree * tree = (TTree *) fPr[1]->Get(PtNm[ipt].Data());
		tree->SetBranchAddress(GenNm.Data(), &GenVar);
		tree->SetBranchAddress(RecNm.Data(), &RecVar);
		Long64_t Nevt = tree->GetEntries();
		for (Long64_t ievt = 0; ievt < Nevt; ++ievt) { // Level 03 loop -- event
			tree->GetEntry(ievt);
			Int_t iorg = AXBnr.findBin(mcPrFlx[1]->signChrg() * RecVar);
			if (iorg < mgFB || iorg > mgLB) continue;
			if (iorg < PtBn[ipt][0] || iorg > PtBn[ipt][1]) continue;
			mcPrFlx[1]->fillEvent(GenVar, RecVar);
		} // Level 02 loop -- event
	} // Level 01 loop -- pattern
	mcPrFlx[1]->finBuild();
	
	TFile * outPr = new TFile(Form("%s_%03d.root", "OutPr", 0), "RECREATE");
	outPr->cd();
	amsFlx->write();
	issRio->write();
	mcPrFlx[1]->write();
	mcApFlx[1]->write();
	outPr->Close();
	
	// Antiproton flux
	TH1D * haprecflx = new TH1D("haprecflx", "haprecflx", axis.Nbin(), axis.Bins());
	for (Int_t ib = 1; ib <= axis.Nbin(); ++ib) {
		Double_t cen = axis.BinCenter(ib, AxisBin::kLog);
		Double_t rer = haprio->GetBinError(ib) / haprio->GetBinContent(ib);
		Double_t val = haprio->GetBinContent(ib) * mcPrFlx[1]->hRecRat->GetBinContent(ib);
		Double_t err = val * rer;
		haprecflx->SetBinContent(ib, val);
		haprecflx->SetBinError  (ib, err);
	}

	ISSFlx * issFlx = new ISSFlx(axis);
	issFlx->buildRecRat(haprecflx);
	issFlx->buildFinFlx(haprecflx);

	// Level 01 loop -- iterator
	// Level 02 loop -- pattern
	// Level 03 loop -- event
	for (Int_t it = 0; it <= NIt; ++it) { // Level 01 loop -- iterator
		std::cout << Form("*****************\n");
		std::cout << Form("**  Iter %04d  **\n", it);
		std::cout << Form("*****************\n");
		
		if (it != 0) UnfoldFlx::Update(issFlx, mcApFlx[1]);
		mcApFlx[1]->buildMigMtx();
		for (Int_t ipt = 0; ipt < NPt; ++ipt) { // Level 02 loop -- pattern
			if (ipt == 3 || ipt == 4) continue;
			TTree * tree = (TTree *) fAp[1]->Get(PtNm[ipt].Data());
			tree->SetBranchAddress(GenNm.Data(), &GenVar);
			tree->SetBranchAddress(RecNm.Data(), &RecVar);
			Long64_t Nevt = tree->GetEntries();
			for (Long64_t ievt = 0; ievt < Nevt; ++ievt) { // Level 03 loop -- event
				tree->GetEntry(ievt);
				Int_t ibin = axis.findBin(mcApFlx[1]->signChrg() * RecVar);
				Int_t iorg = AXBnr.findBin(mcApFlx[1]->signChrg() * RecVar);
				if (iorg < mgFB || iorg > mgLB) continue;
				if (iorg < PtBn[ipt][0] || iorg > PtBn[ipt][1]) continue;
				mcApFlx[1]->fillEvent(GenVar, RecVar);
			} // Level 03 loop -- event
		} // Level 02 loop -- pattern
		mcApFlx[1]->finBuild();

		// Acc
		TH1D * haccr = (TH1D *) mcPrFlx[1]->hEffAcc->Clone(Form("haccr%04d", it));
		haccr->SetTitle(Form("haccr%04d", it));
		if (it != 0) haccr->Divide(mcApFlx[1]->hEffAcc);
		if (it != 0) {
			SplFitPar crrWgt;
			crrWgt.reset();
			crrWgt.fName = "Spl_CrrWgt";
			crrWgt.fRegL = axis.Min(); crrWgt.fRegU = axis.Max();
			crrWgt.fLogX = 1; crrWgt.fLogY = 0;
			crrWgt.fBlxL = 0; crrWgt.fBlxU = 0;
			crrWgt.fN = 0; crrWgt.fX.clear();
			for (Int_t it = 0; it < 5 * UnfoldFlx::NSet; ++it) {
				Double_t val = std::pow(10., UnfoldFlx::XSet[it%UnfoldFlx::NSet] + Double_t(it/UnfoldFlx::NSet));
				if (val < crrWgt.fRegL || val > crrWgt.fRegU) continue;
				crrWgt.fX.push_back(val);
				crrWgt.fN++;
			}
			crrWgt.fit(haccr);
			for (Int_t ib = 1; ib <= axis.Nbin(); ++ib) {
				Double_t cen = axis.BinCenter(ib, AxisBin::kLog);
				haccr->SetBinContent(ib, crrWgt.eval(cen));
				haccr->SetBinError  (ib, 0);
			}
		}


		TFile * outAp = new TFile(Form("%s_%03d.root", "OutAp", it), "RECREATE");
		outAp->cd();
		amsFlx->write();
		issRio->write();
		issFlx->write();
		mcPrFlx[1]->write();
		mcApFlx[1]->write();
		UnfoldFlx::Write();
		RatPamela->Write();
		RatAMSDay->Write();
		haccr->Write();
		outAp->Close();
	} // Level 01 loop -- iterator 

	return 1;
}

#endif // __Unfold_C__
