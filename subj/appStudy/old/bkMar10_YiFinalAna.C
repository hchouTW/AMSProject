#ifndef __YiFinalAna_C__
#define __YiFinalAna_C__

// User defination library
#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/STL/stl.h"
#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/ROOT/selflib.h"

using namespace MgntROOT;

static UInt_t  Acc090_NTrg = 0;
static UInt_t  Acc110_NTrg = 0;

static UInt_t  Pr0510_NTrg = 0;
static UInt_t  Pr1800_NTrg = 0;

static UInt_t  Ap0510_NTrg = 0;
static UInt_t  Ap1800_NTrg = 0;

static TF1 * PdfFlux = new TF1("PdfFlux", "[0]/(TMath::Log([2])-TMath::Log([1]))*TMath::Power(x,-1)");
static TF1 * CdfFlux = new TF1("CdfFlux", "[0]/(TMath::Log([2])-TMath::Log([1]))*TMath::Log(x)");
static TF1 * PdfPowFlux = new TF1("PdfPowFlux", "[0]*([3]+1)/(TMath::Power([2],[3]+1)-TMath::Power([1],[3]+1))*TMath::Power(x,[3])");
static TF1 * CdfPowFlux = new TF1("CdfPowFlux", "[0]/(TMath::Power([2],[3]+1)-TMath::Power([1],[3]+1))*TMath::Power(x,[3]+1)");

int main(int argc, const char ** argv) {
	Style::LoadDefaultEnvironment();

	std::string vtme = (argc==1) ? "" : StrFmt("%03d", std::atoi(argv[1]));

	// Pamela Result
	TString pamelaFilePath = "/afs/cern.ch/user/h/hchou/public/DATASET/physics/antipp3/others/antipp_pamela.root";
	TFile * inFile_Pamela = TFile::Open(pamelaFilePath.Data());
	TGraphErrors * Apprlt_Pamela = (TGraphErrors *) inFile_Pamela->Get("antipp_pamela");
	Apprlt_Pamela->SetMarkerColor(kBlack);
	Apprlt_Pamela->SetLineColor(kBlack);
	//inFile_Pamela->Close();

	// AMS Pub Result
	TString amspubFilePath = "/afs/cern.ch/user/h/hchou/public/DATASET/physics/antipp3/others/20160406MIT.root";
	TFile * inFile_Amspub = TFile::Open(amspubFilePath.Data());
	TGraphErrors * Apprlt_AMSPub = (TGraphErrors *) inFile_Amspub->Get("gpbarp");
	Apprlt_AMSPub->SetMarkerColor(kYellow+1);
	Apprlt_AMSPub->SetLineColor(kYellow+1);
	//inFile_Amspub->Close();

	std::string tagdir = "/afs/cern.ch/work/h/hchou/BSUB_ANALYSIS/antipp5/NOW4__";
	//std::string tagdir = "/afs/cern.ch/work/h/hchou/BSUB_ANALYSIS/antipp5/NOW3Mar02_time/NOW3__";
	std::string goddir = "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/rlt7";

	// ISS
	TFile * fIss = TFile::Open(CStrFmt("%sISS/YiAnalytics.root", tagdir.c_str()));
	if (fIss == nullptr) MgntSys::Error("ISS Files is not exist.");

	// Charge Confusion
	TFile * fPrL1a9 = TFile::Open(CStrFmt("%sL1a9/YiAnalytics.root", tagdir.c_str()));
	TFile * fPrL1o9 = TFile::Open(CStrFmt("%sL1o9/YiAnalytics.root", tagdir.c_str()));
	if (fPrL1a9 == nullptr) MgntSys::Error("L1a9 Files is not exist.");
	if (fPrL1o9 == nullptr) MgntSys::Error("L1o9 Files is not exist.");

	// Cross Section
	TFile * fPrC090 = TFile::Open(CStrFmt("%sC090/YiAnalytics.root", tagdir.c_str()));
	if (fPrC090 == nullptr) MgntSys::Error("C090 Files is not exist.");
	TTree * tPrC090 = (TTree *) fPrC090->Get("DSTRun");
	UInt_t  trgC090 = 0;
	{ UInt_t trg = 0; tPrC090->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tPrC090->GetEntries(); ++it)
		{ tPrC090->GetEntry(it); trgC090 += trg; } }

	TFile * fPrC110 = TFile::Open(CStrFmt("%sC110/YiAnalytics.root", tagdir.c_str()));
	if (fPrC110 == nullptr) MgntSys::Error("C110 Files is not exist.");
	TTree * tPrC110 = (TTree *) fPrC110->Get("DSTRun");
	UInt_t  trgC110 = 0;
	{ UInt_t trg = 0; tPrC110->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tPrC110->GetEntries(); ++it)
		{ tPrC110->GetEntry(it); trgC110 += trg; } }

	// Unfolding
	TFile * fPr0510 = TFile::Open(CStrFmt("%sPr0_510/YiAnalytics.root", tagdir.c_str()));
	if (fPr0510 == nullptr) MgntSys::Error("Pr0510 Files is not exist.");
	TTree * tPr0510 = (TTree *) fPr0510->Get("DSTRun");
	UInt_t  trgPr0510 = 0;
	{ UInt_t trg = 0; tPr0510->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tPr0510->GetEntries(); ++it)
		{ tPr0510->GetEntry(it); trgPr0510 += trg; } }
	
	TFile * fPr1800 = TFile::Open(CStrFmt("%sPr1800/YiAnalytics.root", tagdir.c_str()));
	if (fPr1800 == nullptr) MgntSys::Error("Pr1800 Files is not exist.");
	TTree * tPr1800 = (TTree *) fPr1800->Get("DSTRun");
	UInt_t  trgPr1800 = 0;
	{ UInt_t trg = 0; tPr1800->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tPr1800->GetEntries(); ++it)
		{ tPr1800->GetEntry(it); trgPr1800 += trg; } }
	
	TFile * fAp0510 = TFile::Open(CStrFmt("%sAp0_510/YiAnalytics.root", tagdir.c_str()));
	if (fAp0510 == nullptr) MgntSys::Error("Ap0510 Files is not exist.");
	TTree * tAp0510 = (TTree *) fAp0510->Get("DSTRun");
	UInt_t  trgAp0510 = 0;
	{ UInt_t trg = 0; tAp0510->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tAp0510->GetEntries(); ++it)
		{ tAp0510->GetEntry(it); trgAp0510 += trg; } }
	
	TFile * fAp1800 = TFile::Open(CStrFmt("%sAp1800/YiAnalytics.root", tagdir.c_str()));
	if (fAp1800 == nullptr) MgntSys::Error("Ap1800 Files is not exist.");
	TTree * tAp1800 = (TTree *) fAp1800->Get("DSTRun");
	UInt_t  trgAp1800 = 0;
	{ UInt_t trg = 0; tAp1800->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tAp1800->GetEntries(); ++it)
		{ tAp1800->GetEntry(it); trgAp1800 += trg; } }


		std::cout << "TRG " << trgAp1800 << std::endl;

	// Binning
	Int_t lBin[2] = {  1+2, 19+2};
	Int_t iBin[2] = {  9+2, 42+2};
	Int_t hBin[2] = { 32+2, 57+2}; // L1 52(53), L9 52(54), Fs 57
	Int_t cbBin[4] = { 0+2, 13+2, 35+2, 57+2 };
	
	// Binning
	//Int_t lBin[2] = {  1,  6 };
	//Int_t iBin[2] = {  5, 12 };
	//Int_t hBin[2] = { 10, 16 };
	//Int_t cbBin[4] = { 0, 5, 10, 16 };

	//TFile * fAppBin = TFile::Open("/afs/cern.ch/user/h/hchou/public/DATABASE/physics/binning/antipp3_rig.root");
	//Axis AXnr((TH1*)fAppBin->Get("hbin0"), Axis::kX);
  //Axis AXir((TH1*)fAppBin->Get("hbin1"), Axis::kX);
	//fAppBin->Close();
			
	// rigidity binning
	Axis AXnr = Axis("Rigidity [GV]",
		{   0.50,   0.80, // extern bins
		    1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
		    3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
			  8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
			 19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
			 41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
			 93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00, 
			800.00 } ); // extern bin
	Axis AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
	
	//Axis AXnr = Axis("Rigidity [GV]",
	//	{   1.00,   1.51,   2.00,   3.00,   4.02, 
	//	    5.37,   7.09,   9.26,  12.00,  15.30, 
	//		 19.50,  24.70,  31.10,  38.90,  48.50, 
	//		 60.30,  74.90 } );
	//Axis AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
	
	//Axis AXnr = Axis("Rigidity [GV]",
	//	{   1.00,   1.46,   2.00,   3.00,   4.12, 
	//	    5.00,   6.00,   7.10,   8.30,   9.62, 
	//	   11.04,  12.59,  14.25,  16.05,  17.98, 
	//	   20.04,  22.25,  24.62,  27.25,  30.21, 
	//	   35.36,  40.00 } );
	//Axis AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);

	PdfEditor editor(PdfEditor::kSlice, StrFmt("YiAna%s", vtme.c_str()), goddir);
	
	// Antiproton-to-Proton Flux Ratio
	Graph * Appnum_pr   = Graph::New("apprlt_pr"  , "", "|Rigidity| [GV]", "N_{p}");
	Graph * Appnum_ap   = Graph::New("apprlt_ap"  , "", "|Rigidity| [GV]", "N_{#bar{p}}");
	Graph * Apprlt_stat = Graph::New("apprlt_stat", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * Apprlt_totl = Graph::New("apprlt_totl", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * Appcrr_accp = Graph::New("appcrr_accp", "", "|Rigidity| [GV]", "Acceptance Correction");
	Graph * Apperr_totl = Graph::New("apperr_totl", "", "|Rigidity| [GV]", "Total Error [%]");
	Graph * Apperr_stat = Graph::New("apperr_stat", "", "|Rigidity| [GV]", "Statistics Error [%]");
	Graph * Apperr_sysm = Graph::New("apperr_sysm", "", "|Rigidity| [GV]", "Systematic Error [%]");
	Graph * Apperr_temp = Graph::New("apperr_temp", "", "|Rigidity| [GV]", "Template Error [%]");
	Graph * Apperr_accp = Graph::New("apperr_accp", "", "|Rigidity| [GV]", "Acceptance Error [%]");
	Graph * Apperr_xsec = Graph::New("apperr_xsec", "", "|Rigidity| [GV]", "Cross-Section Error [%]");
	Graph * Appfit_nchi = Graph::New("appfit_nchi", "", "|Rigidity| [GV]", "Chisquare/NDF");


	//----  Low Energy  ----//
	std::cout << "\n=======  Low Energy  =======\n";
	Graph * LAppnum_pr   = Graph::New("lapprlt_pr"  , "", "|Rigidity| [GV]", "N_{p}");
	Graph * LAppnum_ap   = Graph::New("lapprlt_ap"  , "", "|Rigidity| [GV]", "N_{#bar{p}}");
	Graph * LApprlt_stat = Graph::New("lapprlt_stat", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * LApprlt_totl = Graph::New("lapprlt_totl", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * LAppcrr_accp = Graph::New("lappcrr_accp", "", "|Rigidity| [GV]", "Acceptance Correction");
	Graph * LApperr_totl = Graph::New("lapperr_totl", "", "|Rigidity| [GV]", "Total Error [%]");
	Graph * LApperr_stat = Graph::New("lapperr_stat", "", "|Rigidity| [GV]", "Statistics Error [%]");
	Graph * LApperr_sysm = Graph::New("lapperr_sysm", "", "|Rigidity| [GV]", "Systematic Error [%]");
	Graph * LApperr_temp = Graph::New("lapperr_temp", "", "|Rigidity| [GV]", "Template Error [%]");
	Graph * LApperr_accp = Graph::New("lapperr_accp", "", "|Rigidity| [GV]", "Acceptance Error [%]");
	Graph * LApperr_xsec = Graph::New("lapperr_xsec", "", "|Rigidity| [GV]", "Cross-Section Error [%]");
	Graph * LAppfit_nchi = Graph::New("lappfit_nchi", "", "|Rigidity| [GV]", "Chisquare/NDF");
	
	Hist * hL2D_pos = Hist::New((TH1*)fIss->Get(CStrFmt("hLp_TofM%s", vtme.c_str())));
	Hist * hL2D_neg = Hist::New((TH1*)fIss->Get(CStrFmt("hLn_TofM%s", vtme.c_str())));
	Hist * hL2D_sig = Hist::New((TH1*)fIss->Get(CStrFmt("hLs_TofM%s", vtme.c_str())));
	Hist * hL2D_bkg = Hist::New((TH1*)fIss->Get(CStrFmt("hLb_TofM%s", vtme.c_str())));
		
	std::vector<Hist *>&& hLVec_pos = Hist::Project(Hist::kProjY, hL2D_pos);
	std::vector<Hist *>&& hLVec_neg = Hist::Project(Hist::kProjY, hL2D_neg);
	std::vector<Hist *>&& hLVec_sig = Hist::Project(Hist::kProjY, hL2D_sig);
	std::vector<Hist *>&& hLVec_bkg = Hist::Project(Hist::kProjY, hL2D_bkg);

	Hist * hL_C090 = Hist::New("hLp_C090_MCEvt", (TH1*)fPrC090->Get("hLp_MCEvt"));
	Hist * hL_C110 = Hist::New("hLn_C110_MCEvt", (TH1*)fPrC110->Get("hLp_MCEvt"));
	
	Hist * hL_ACCPR = Hist::New("hLp_ACCPR_MCEvt", (TH1*)fPr1800->Get("hLp_MCEvt"));
	Hist * hL_ACCAP = Hist::New("hLn_ACCAP_MCEvt", (TH1*)fAp1800->Get("hLn_MCEvt"));

	Hist * hLp_Cutflow = Hist::New("hLp_Cutflow", (TH1*)fIss->Get(CStrFmt("hLp_Cutflow%s", vtme.c_str())));
	Hist * hLn_Cutflow = Hist::New("hLn_Cutflow", (TH1*)fIss->Get(CStrFmt("hLn_Cutflow%s", vtme.c_str())));

	for (Int_t ib = lBin[0]; ib <= lBin[1]; ++ib) {
		std::cout << CStrFmt("L %02d FIT\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		Hist * hL_pos = hLVec_pos.at(ib);
		Hist * hL_neg = hLVec_neg.at(ib);
		Hist * hL_sig = hLVec_sig.at(ib);
		Hist * hL_bkg = hLVec_bkg.at(ib);

		Double_t numOfPr     = (*hL_pos)()->Integral();
		Double_t errOfPrStat = std::sqrt(numOfPr);
		
		Double_t numOfSmp = (*hL_neg)()->Integral();
		if (numOfSmp < 20.) continue;

		Fit::RooParam param("bta", hL_neg, Fit::RooParam::VLIST({ hL_bkg, hL_sig }));
		Fit::RooSysResult rlt(param, true, 100);
		if (!rlt.valid()) continue;
		if (!rlt.sysValid()) continue;
		Double_t numOfAp     = rlt.val(1);
		Double_t errOfApStat = rlt.err(1);
		Double_t errOfApTemp = rlt.sysErr(1);

		//Double_t crrOfAppAccp = 3.38751e-01 * TMath::Power(cen, -7.73513e-01) + 1.03987e+00;
		Double_t crrOfAppAccp = 1.03640e+00 + 2.68627e-01 * TMath::Power(cen, -5.83324e-01) + 1.13567e-01 * TMath::Power(cen, -5.60999e+00);
		Double_t valOfApp     = crrOfAppAccp * numOfAp / numOfPr; 
		Double_t errOfAppStat = valOfApp * std::sqrt(
		                        errOfPrStat*errOfPrStat/numOfPr/numOfPr + 
														errOfApStat*errOfApStat/numOfAp/numOfAp
														);
		Double_t errOfAppTemp = valOfApp * std::sqrt(
														errOfApTemp*errOfApTemp/numOfAp/numOfAp
														);
		Double_t errOfAppAccp = valOfApp * 0.015;
		Double_t errOfAppXsec = valOfApp * 0.03;
		Double_t errOfAppSysm = std::sqrt(
														errOfAppTemp*errOfAppTemp+
														errOfAppAccp*errOfAppAccp+
														errOfAppXsec*errOfAppXsec
														);
		Double_t errOfAppTotl = std::sqrt(
		                        errOfAppStat*errOfAppStat+
														errOfAppSysm*errOfAppSysm
														);

		LAppnum_pr->pushPoint(Graph::Point_t(cen, numOfPr));
		LAppnum_ap->pushPoint(Graph::Point_t(cen, numOfAp));
		LApprlt_stat->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppStat));
		LApprlt_totl->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
		//LAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		LApperr_totl->pushPoint(Graph::Point_t(cen, 100. * errOfAppTotl/valOfApp));
		LApperr_stat->pushPoint(Graph::Point_t(cen, 100. * errOfAppStat/valOfApp));
		LApperr_sysm->pushPoint(Graph::Point_t(cen, 100. * errOfAppSysm/valOfApp));
		LApperr_temp->pushPoint(Graph::Point_t(cen, 100. * errOfAppTemp/valOfApp));
		//LApperr_accp->pushPoint(Graph::Point_t(cen, 100. * errOfAppAccp/valOfApp));
		//LApperr_xsec->pushPoint(Graph::Point_t(cen, 100. * errOfAppXsec/valOfApp));
		LAppfit_nchi->pushPoint(Graph::Point_t(cen, rlt.chisq()/rlt.ndf()));

		editor.create(StrFmt("L%02d", ib));
		editor.cd(1, Canvas::AxisScl_t(0, 1));
		rlt.samp() ->setStyle(Area(), Line(kBlack), Marker(kBlack));
		rlt.sumt() ->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
		rlt.temp(0)->setStyle(Area(), Line(kBlue), Marker(kBlue));
		rlt.temp(1)->setStyle(Area(), Line(kRed), Marker(kRed));
		rlt.samp()->draw("pe");
		Hist::Collect(StrFmt("fit%03d", ib), "", Hist::VLIST({ rlt.sumt(), rlt.temp(0), rlt.temp(1) }))->Draw("same nostack hist");
		Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.85, 0.92, 32));
		Text::Draw(Text::Txt_t(StrFmt("N_{Background} ( %.2f +- stat %.2f sys %.2f )", rlt.val(0), rlt.err(0), rlt.sysErr(0)), kBlue, 0.03), Text::Att_t(0.15, 0.85, 12));
		Text::Draw(Text::Txt_t(StrFmt("N_{Antiproton} ( %.2f +- stat %.2f sys %.2f )", rlt.val(1), rlt.err(1), rlt.sysErr(1)), kRed, 0.03), Text::Att_t(0.15, 0.81, 12));
		Text::Draw(Text::Txt_t(StrFmt("#bar{p}/p Ratio ( %E +- stat %E )", valOfApp, errOfAppStat), kBlack, 0.03), Text::Att_t(0.15, 0.77, 12));
		Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.15, 0.73, 12));
		//Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.80, 0.92, 32));
		//Text::Draw(Text::Txt_t(StrFmt("Background"), kBlue, 0.03), Text::Att_t(0.20, 0.85, 12));
		//Text::Draw(Text::Txt_t(StrFmt("Antiproton"), kRed, 0.03), Text::Att_t(0.20, 0.81, 12));
		//Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.20, 0.77, 12));
		editor.save();

		if (ib > cbBin[0] && ib <= cbBin[1]) {
			Apprlt_stat->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppStat));
			Apprlt_totl->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
			Appfit_nchi->pushPoint(Graph::Point_t(cen, rlt.chisq()/rlt.ndf()));
		}
	}
	
	for (Int_t ib = 1; ib <= AXnr.nbin(); ++ib) {
		std::cout << CStrFmt("L %02d XSEC\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		if (AXnr.bins(ib-1) >= 1.0 && AXnr.bins(ib) <= 200.) {
			CdfFlux->SetParameters(trgC090, 1.0, 200.0);
			Double_t numOfSel_C090 = (*hL_C090)()->GetBinContent(ib);
			Double_t numOfTrg_C090 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratC090 = numOfSel_C090 / numOfTrg_C090;
			CdfFlux->SetParameters(trgC110, 1.0, 200.0);
			Double_t numOfSel_C110 = (*hL_C110)()->GetBinContent(ib);
			Double_t numOfTrg_C110 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratC110 = numOfSel_C110 / numOfTrg_C110;
		
			if (numOfSel_C090 < 100 || numOfSel_C110 < 100) continue;

			Double_t raterr  = std::sqrt(ratC090/numOfTrg_C090 + ratC110/numOfTrg_C110);
			Double_t ratmius = std::fabs(ratC090 - ratC110);
			Double_t ratplus = std::fabs(ratC090 + ratC110);

			Double_t errOfAppXsec = ratmius / ratplus;
			Double_t errOfAppXsecErr = errOfAppXsec * raterr * std::sqrt((1./ratmius/ratmius + 1./ratplus/ratplus));
			if ((errOfAppXsecErr / errOfAppXsec) > 0.5 || !MgntNum::Valid(errOfAppXsec)) continue;
			LApperr_xsec->pushPointWithError(Graph::Point_t(cen, 100. * errOfAppXsec), Graph::Error_t(0.0, 100. * errOfAppXsecErr));
		}
	}

	for (Int_t ib = 1; ib <= AXnr.nbin(); ++ib) {
		std::cout << CStrFmt("L %02d ACCP\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		if (AXnr.bins(ib-1) >= 1.0 && AXnr.bins(ib) <= 800.) {
			CdfFlux->SetParameters(trgPr1800, 1.0, 800.0);
			Double_t numOfSel_ACCPR = (*hL_ACCPR)()->GetBinContent(ib);
			Double_t numOfTrg_ACCPR = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratACCPR = numOfSel_ACCPR / numOfTrg_ACCPR;
			CdfFlux->SetParameters(trgAp1800, 1.0, 800.0);
			Double_t numOfSel_ACCAP = (*hL_ACCAP)()->GetBinContent(ib);
			Double_t numOfTrg_ACCAP = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratACCAP = numOfSel_ACCAP / numOfTrg_ACCAP;
	
			std::cout << "XXXXXXXXXXX\n";
			std::cout << numOfSel_ACCPR << " " << numOfSel_ACCAP << std::endl;
			if (numOfSel_ACCPR < 100 || numOfSel_ACCAP < 100) continue;
			std::cout << "YYYYYYYYYY\n";

			Double_t crrOfAppAccp = ratACCPR / ratACCAP;
			Double_t crrOfAppAccpErr = crrOfAppAccp * std::sqrt(1./numOfSel_ACCPR+1./numOfSel_ACCAP);
			if ((crrOfAppAccpErr / crrOfAppAccp) > 0.02 || !MgntNum::Valid(crrOfAppAccp)) continue;
			LAppcrr_accp->pushPointWithError(Graph::Point_t(cen, crrOfAppAccp), Graph::Error_t(0.0,  crrOfAppAccpErr));
			
			//Double_t crrOfAppAccpPred = 3.38751e-01 * TMath::Power(cen, -7.73513e-01) + 1.03987e+00;
			Double_t crrOfAppAccpPred = 1.03640e+00 + 2.68627e-01 * TMath::Power(cen, -5.83324e-01) + 1.13567e-01 * TMath::Power(cen, -5.60999e+00);
			Double_t errOfAppAccp = (crrOfAppAccp - crrOfAppAccpPred) / crrOfAppAccpPred;
			LApperr_accp->pushPointWithError(Graph::Point_t(cen, 100. * errOfAppAccp), Graph::Error_t(0.0, 100. * crrOfAppAccpErr / crrOfAppAccpPred));
		}
	}

	Int_t  nLCut = (*hLp_Cutflow)()->GetNbinsY()-1;
	Axis AXLcut("cut", nLCut, 0.0, Double_t(nLCut));
	Hist * hLp_Eff = Hist::New("hLp_Eff", "COMp Efficiency Rig", AXnr, AXLcut);
	for (Int_t irig = 1; irig <= (*hLp_Eff)()->GetNbinsX(); ++irig) {
		for (Int_t icut = 1; icut <= (*hLp_Eff)()->GetNbinsY(); ++icut) {
			Double_t pw  = (*hLp_Cutflow)()->GetBinContent(irig, icut + 1);
			Double_t tw  = (*hLp_Cutflow)()->GetBinContent(irig, icut);
			if (MgntNum::EqualToZero(pw) || MgntNum::EqualToZero(tw)) continue; 
			Double_t eff = pw / tw;
			Double_t sgm = std::sqrt(eff * (1. - eff) / tw);
			Double_t eru = ((eff + sgm) > 1.0) ? (1.0 - eff) : sgm;
			Double_t erl = ((eff - sgm) > 0.0) ? eff         : sgm;
			(*hLp_Eff)()->SetBinContent(irig, icut, eff);
			(*hLp_Eff)()->SetBinError  (irig, icut, sgm);
		}
	}
	Hist * hLn_Eff = Hist::New("hLn_Eff", "COMn Efficiency Rig", AXnr, AXLcut);
	for (Int_t irig = 1; irig <= (*hLn_Eff)()->GetNbinsX(); ++irig) {
		for (Int_t icut = 1; icut <= (*hLn_Eff)()->GetNbinsY(); ++icut) {
			Double_t pw  = (*hLn_Cutflow)()->GetBinContent(irig, icut + 1);
			Double_t tw  = (*hLn_Cutflow)()->GetBinContent(irig, icut);
			if (MgntNum::EqualToZero(pw) || MgntNum::EqualToZero(tw)) continue; 
			Double_t eff = pw / tw;
			Double_t sgm = std::sqrt(eff * (1. - eff) / tw);
			Double_t eru = ((eff + sgm) > 1.0) ? (1.0 - eff) : sgm;
			Double_t erl = ((eff - sgm) > 0.0) ? eff         : sgm;
			(*hLn_Eff)()->SetBinContent(irig, icut, eff);
			(*hLn_Eff)()->SetBinError  (irig, icut, sgm);
		}
	}

	//----  Intermedia Energy  ----//
	std::cout << "\n=======  Intermedia Energy  =======\n";
	Graph * IAppnum_pr   = Graph::New("iapprlt_ap"  , "", "|Rigidity| [GV]", "N_{p}");
	Graph * IAppnum_ap   = Graph::New("iapprlt_pr"  , "", "|Rigidity| [GV]", "N_{#bar{p}}");
	Graph * IApprlt_stat = Graph::New("iapprlt_stat", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * IApprlt_totl = Graph::New("iapprlt_totl", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * IAppcrr_accp = Graph::New("iappcrr_accp", "", "|Rigidity| [GV]", "Acceptance Correction");
	Graph * IApperr_totl = Graph::New("iapperr_totl", "", "|Rigidity| [GV]", "Total Error [%]");
	Graph * IApperr_stat = Graph::New("iapperr_stat", "", "|Rigidity| [GV]", "Statistics Error [%]");
	Graph * IApperr_sysm = Graph::New("iapperr_sysm", "", "|Rigidity| [GV]", "Systematic Error [%]");
	Graph * IApperr_temp = Graph::New("iapperr_temp", "", "|Rigidity| [GV]", "Template Error [%]");
	Graph * IApperr_accp = Graph::New("iapperr_accp", "", "|Rigidity| [GV]", "Acceptance Error [%]");
	Graph * IApperr_xsec = Graph::New("iapperr_xsec", "", "|Rigidity| [GV]", "Cross-Section Error [%]");
	Graph * IAppfit_nchi = Graph::New("iappfit_nchi", "", "|Rigidity| [GV]", "Chisquare/NDF");
	
	Hist * hI2D_pos = Hist::New((TH1*)fIss->Get(CStrFmt("hIp_TrdEst%s", vtme.c_str())));
	Hist * hI2D_neg = Hist::New((TH1*)fIss->Get(CStrFmt("hIn_TrdEst%s", vtme.c_str())));
	Hist * hI2D_sig = Hist::New((TH1*)fIss->Get(CStrFmt("hIs_TrdEst%s", vtme.c_str())));
	Hist * hI2D_bkg = Hist::New((TH1*)fIss->Get(CStrFmt("hIb_TrdEst%s", vtme.c_str())));
		
	std::vector<Hist *>&& hIVec_pos = Hist::Project(Hist::kProjY, hI2D_pos);
	std::vector<Hist *>&& hIVec_neg = Hist::Project(Hist::kProjY, hI2D_neg);
	std::vector<Hist *>&& hIVec_sig = Hist::Project(Hist::kProjY, hI2D_sig);
	std::vector<Hist *>&& hIVec_bkg = Hist::Project(Hist::kProjY, hI2D_bkg);
	
	Hist * hI_C090 = Hist::New("hIp_C090_MCEvt", (TH1*)fPrC090->Get("hIp_MCEvt"));
	Hist * hI_C110 = Hist::New("hIn_C110_MCEvt", (TH1*)fPrC110->Get("hIp_MCEvt"));
	
	Hist * hI_ACCPR = Hist::New("hIp_ACCPR_MCEvt", (TH1*)fPr1800->Get("hIp_MCEvt"));
	Hist * hI_ACCAP = Hist::New("hIn_ACCAP_MCEvt", (TH1*)fAp1800->Get("hIn_MCEvt"));
	
	Hist * hIp_Cutflow = Hist::New("hIp_Cutflow", (TH1*)fIss->Get(CStrFmt("hIp_Cutflow%s", vtme.c_str())));
	Hist * hIn_Cutflow = Hist::New("hIn_Cutflow", (TH1*)fIss->Get(CStrFmt("hIn_Cutflow%s", vtme.c_str())));

	for (Int_t ib = iBin[0]; ib <= iBin[1]; ++ib) {
		std::cout << CStrFmt("I %02d FIT\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		Hist * hI_pos = hIVec_pos.at(ib);
		Hist * hI_neg = hIVec_neg.at(ib);
		Hist * hI_sig = hIVec_sig.at(ib);
		Hist * hI_bkg = hIVec_bkg.at(ib);

		Double_t numOfPr     = (*hI_pos)()->Integral();
		Double_t errOfPrStat = std::sqrt(numOfPr);

		Double_t numOfSmp = (*hI_neg)()->Integral();
		if (numOfSmp < 20.) continue;

		Fit::RooParam param("trd", hI_neg, Fit::RooParam::VLIST({ hI_bkg, hI_sig }));
		Fit::RooSysResult rlt(param, true, 100);
		if (!rlt.valid()) continue;
		if (!rlt.sysValid()) continue;
		Double_t numOfAp     = rlt.val(1);
		Double_t errOfApStat = rlt.err(1);
		Double_t errOfApTemp = rlt.sysErr(1);

		//Double_t crrOfAppAccp = 2.34627e-01 * TMath::Power(cen, -6.23656e-01) + 1.03143e+00;
		Double_t crrOfAppAccp = 1.03092e+00 + 2.40482e-01 * TMath::Power(cen, -6.23688e-01);
		Double_t valOfApp     = crrOfAppAccp * numOfAp / numOfPr; 
		Double_t errOfAppStat = valOfApp * std::sqrt(
		                        errOfPrStat*errOfPrStat/numOfPr/numOfPr + 
														errOfApStat*errOfApStat/numOfAp/numOfAp
														);
		Double_t errOfAppTemp = valOfApp * std::sqrt(
														errOfApTemp*errOfApTemp/numOfAp/numOfAp
														);
		Double_t errOfAppAccp = valOfApp * 0.015;
		Double_t errOfAppXsec = valOfApp * 0.03;
		Double_t errOfAppSysm = std::sqrt(
														errOfAppTemp*errOfAppTemp+
														errOfAppAccp*errOfAppAccp+
														errOfAppXsec*errOfAppXsec
														);
		Double_t errOfAppTotl = std::sqrt(
		                        errOfAppStat*errOfAppStat+
														errOfAppSysm*errOfAppSysm
														);

		IAppnum_pr->pushPoint(Graph::Point_t(cen, numOfPr));
		IAppnum_ap->pushPoint(Graph::Point_t(cen, numOfAp));
		IApprlt_stat->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppStat));
		IApprlt_totl->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
		//IAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		IApperr_totl->pushPoint(Graph::Point_t(cen, 100. * errOfAppTotl/valOfApp));
		IApperr_stat->pushPoint(Graph::Point_t(cen, 100. * errOfAppStat/valOfApp));
		IApperr_sysm->pushPoint(Graph::Point_t(cen, 100. * errOfAppSysm/valOfApp));
		IApperr_temp->pushPoint(Graph::Point_t(cen, 100. * errOfAppTemp/valOfApp));
		//IApperr_accp->pushPoint(Graph::Point_t(cen, 100. * errOfAppAccp/valOfApp));
		//IApperr_xsec->pushPoint(Graph::Point_t(cen, 100. * errOfAppXsec/valOfApp));
		IAppfit_nchi->pushPoint(Graph::Point_t(cen, rlt.chisq()/rlt.ndf()));

		editor.create(StrFmt("I%02d", ib));
		editor.cd(1, Canvas::AxisScl_t(0, 1));
		rlt.samp() ->setStyle(Area(), Line(kBlack), Marker(kBlack));
		rlt.sumt() ->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
		rlt.temp(0)->setStyle(Area(), Line(kBlue), Marker(kBlue));
		rlt.temp(1)->setStyle(Area(), Line(kRed), Marker(kRed));
		rlt.samp()->draw("pe");
		Hist::Collect(StrFmt("fit%03d", ib), "", Hist::VLIST({ rlt.sumt(), rlt.temp(0), rlt.temp(1) }))->Draw("same nostack hist");
		Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.85, 0.92, 32));
		Text::Draw(Text::Txt_t(StrFmt("N_{Background} ( %.2f +- stat %.2f sys %.2f )", rlt.val(0), rlt.err(0), rlt.sysErr(0)), kBlue, 0.03), Text::Att_t(0.15, 0.85, 12));
		Text::Draw(Text::Txt_t(StrFmt("N_{Antiproton} ( %.2f +- stat %.2f sys %.2f )", rlt.val(1), rlt.err(1), rlt.sysErr(1)), kRed, 0.03), Text::Att_t(0.15, 0.81, 12));
		Text::Draw(Text::Txt_t(StrFmt("#bar{p}/p Ratio ( %E +- stat %E )", valOfApp, errOfAppStat), kBlack, 0.03), Text::Att_t(0.15, 0.77, 12));
		Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.15, 0.73, 12));
		//Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.80, 0.92, 32));
		//Text::Draw(Text::Txt_t(StrFmt("Background"), kBlue, 0.03), Text::Att_t(0.20, 0.85, 12));
		//Text::Draw(Text::Txt_t(StrFmt("Antiproton"), kRed, 0.03), Text::Att_t(0.20, 0.81, 12));
		//Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.20, 0.77, 12));
		editor.save();
		
		if (ib > cbBin[1] && ib <= cbBin[2]) {
			Apprlt_stat->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppStat));
			Apprlt_totl->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
			Appfit_nchi->pushPoint(Graph::Point_t(cen, rlt.chisq()/rlt.ndf()));
		}
	}
		
	for (Int_t ib = 1; ib <= AXnr.nbin(); ++ib) {
		std::cout << CStrFmt("I %02d XSEC\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		if (AXnr.bins(ib-1) >= 1.0 && AXnr.bins(ib) <= 200.) {
			CdfFlux->SetParameters(trgC090, 1.0, 200.0);
			Double_t numOfSel_C090 = (*hI_C090)()->GetBinContent(ib);
			Double_t numOfTrg_C090 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratC090 = numOfSel_C090 / numOfTrg_C090;
			CdfFlux->SetParameters(trgC110, 1.0, 200.0);
			Double_t numOfSel_C110 = (*hI_C110)()->GetBinContent(ib);
			Double_t numOfTrg_C110 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratC110 = numOfSel_C110 / numOfTrg_C110;

			if (numOfSel_C090 < 100 || numOfSel_C110 < 100) continue;

			Double_t raterr  = std::sqrt(ratC090/numOfTrg_C090 + ratC110/numOfTrg_C110);
			Double_t ratmius = std::fabs(ratC090 - ratC110);
			Double_t ratplus = std::fabs(ratC090 + ratC110);

			Double_t errOfAppXsec = ratmius / ratplus;
			Double_t errOfAppXsecErr = errOfAppXsec * raterr * std::sqrt((1./ratmius/ratmius + 1./ratplus/ratplus));
			if ((errOfAppXsecErr / errOfAppXsec) > 0.5 || !MgntNum::Valid(errOfAppXsec)) continue;
			IApperr_xsec->pushPointWithError(Graph::Point_t(cen, 100. * errOfAppXsec), Graph::Error_t(0.0, 100. * errOfAppXsecErr));
		}
	}

	for (Int_t ib = 1; ib <= AXnr.nbin(); ++ib) {
		std::cout << CStrFmt("I %02d ACCP\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		if (AXnr.bins(ib-1) >= 1.0 && AXnr.bins(ib) <= 800.) {
			CdfFlux->SetParameters(trgPr1800, 1.0, 800.0);
			Double_t numOfSel_ACCPR = (*hI_ACCPR)()->GetBinContent(ib);
			Double_t numOfTrg_ACCPR = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratACCPR = numOfSel_ACCPR / numOfTrg_ACCPR;
			CdfFlux->SetParameters(trgAp1800, 1.0, 800.0);
			Double_t numOfSel_ACCAP = (*hI_ACCAP)()->GetBinContent(ib);
			Double_t numOfTrg_ACCAP = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratACCAP = numOfSel_ACCAP / numOfTrg_ACCAP;
			
			if (numOfSel_ACCPR < 100 || numOfSel_ACCAP < 100) continue;
			
			Double_t crrOfAppAccp = ratACCPR / ratACCAP;
			Double_t crrOfAppAccpErr = crrOfAppAccp * std::sqrt(1./numOfSel_ACCPR+1./numOfSel_ACCAP);
			if ((crrOfAppAccpErr / crrOfAppAccp) > 0.02 || !MgntNum::Valid(crrOfAppAccp)) continue;
			IAppcrr_accp->pushPointWithError(Graph::Point_t(cen, crrOfAppAccp), Graph::Error_t(0.0,  crrOfAppAccpErr));
			
			//Double_t crrOfAppAccpPred = 2.34627e-01 * TMath::Power(cen, -6.23656e-01) + 1.03143e+00;
			Double_t crrOfAppAccpPred = 1.03092e+00 + 2.40482e-01 * TMath::Power(cen, -6.23688e-01);
			Double_t errOfAppAccp = (crrOfAppAccp - crrOfAppAccpPred) / crrOfAppAccpPred;
			IApperr_accp->pushPointWithError(Graph::Point_t(cen, 100. * errOfAppAccp), Graph::Error_t(0.0, 100. * crrOfAppAccpErr / crrOfAppAccpPred));
		}
	}
	
	Int_t  nICut = (*hIp_Cutflow)()->GetNbinsY()-1;
	Axis AXIcut("cut", nICut, 0.0, Double_t(nICut));
	Hist * hIp_Eff = Hist::New("hIp_Eff", "COMp Efficiency Rig", AXnr, AXIcut);
	for (Int_t irig = 1; irig <= (*hIp_Eff)()->GetNbinsX(); ++irig) {
		for (Int_t icut = 1; icut <= (*hIp_Eff)()->GetNbinsY(); ++icut) {
			Double_t pw  = (*hIp_Cutflow)()->GetBinContent(irig, icut + 1);
			Double_t tw  = (*hIp_Cutflow)()->GetBinContent(irig, icut);
			if (MgntNum::EqualToZero(pw) || MgntNum::EqualToZero(tw)) continue; 
			Double_t eff = pw / tw;
			Double_t sgm = std::sqrt(eff * (1. - eff) / tw);
			Double_t eru = ((eff + sgm) > 1.0) ? (1.0 - eff) : sgm;
			Double_t erl = ((eff - sgm) > 0.0) ? eff         : sgm;
			(*hIp_Eff)()->SetBinContent(irig, icut, eff);
			(*hIp_Eff)()->SetBinError  (irig, icut, sgm);
		}
	}
	Hist * hIn_Eff = Hist::New("hIn_Eff", "COMn Efficiency Rig", AXnr, AXIcut);
	for (Int_t irig = 1; irig <= (*hIn_Eff)()->GetNbinsX(); ++irig) {
		for (Int_t icut = 1; icut <= (*hIn_Eff)()->GetNbinsY(); ++icut) {
			Double_t pw  = (*hIn_Cutflow)()->GetBinContent(irig, icut + 1);
			Double_t tw  = (*hIn_Cutflow)()->GetBinContent(irig, icut);
			if (MgntNum::EqualToZero(pw) || MgntNum::EqualToZero(tw)) continue; 
			Double_t eff = pw / tw;
			Double_t sgm = std::sqrt(eff * (1. - eff) / tw);
			Double_t eru = ((eff + sgm) > 1.0) ? (1.0 - eff) : sgm;
			Double_t erl = ((eff - sgm) > 0.0) ? eff         : sgm;
			(*hIn_Eff)()->SetBinContent(irig, icut, eff);
			(*hIn_Eff)()->SetBinError  (irig, icut, sgm);
		}
	}
	
	//----  High Energy  ----//
	std::cout << "\n=======  High Energy  =======\n";
	Graph * HAppnum_pr   = Graph::New("happrlt_ap"  , "", "|Rigidity| [GV]", "N_{p}");
	Graph * HAppnum_ap   = Graph::New("happrlt_pr"  , "", "|Rigidity| [GV]", "N_{#bar{p}}");
	Graph * HApprlt_stat = Graph::New("happrlt_stat", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * HApprlt_totl = Graph::New("happrlt_totl", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * HAppcrr_accp = Graph::New("happcrr_accp", "", "|Rigidity| [GV]", "Acceptance Correction");
	Graph * HApperr_totl = Graph::New("happerr_totl", "", "|Rigidity| [GV]", "Total Error [%]");
	Graph * HApperr_stat = Graph::New("happerr_stat", "", "|Rigidity| [GV]", "Statistics Error [%]");
	Graph * HApperr_sysm = Graph::New("happerr_sysm", "", "|Rigidity| [GV]", "Systematic Error [%]");
	Graph * HApperr_temp = Graph::New("happerr_temp", "", "|Rigidity| [GV]", "Template Error [%]");
	Graph * HApperr_accp = Graph::New("happerr_accp", "", "|Rigidity| [GV]", "Acceptance Error [%]");
	Graph * HApperr_xsec = Graph::New("happerr_xsec", "", "|Rigidity| [GV]", "Cross-Section Error [%]");

	Hist * hH2D_L1pos = Hist::New((TH1*)fIss->Get(CStrFmt("hHp_CCestL1%s", vtme.c_str())));
	Hist * hH2D_L1neg = Hist::New((TH1*)fIss->Get(CStrFmt("hHn_CCestL1%s", vtme.c_str())));
	Hist * hH2D_L1sig = Hist::New("hHs_CCestL1", (TH1*)fIss->Get(CStrFmt("hHp_CCestL1%s", vtme.c_str())));
	Hist * hH2D_L1bkg = Hist::New("hHb_CCestL1", (TH1*)fPrL1o9->Get("hHn_CCestL1"));
	
	std::vector<Hist *>&& hHVec_L1pos = Hist::Project(Hist::kProjY, hH2D_L1pos);
	std::vector<Hist *>&& hHVec_L1neg = Hist::Project(Hist::kProjY, hH2D_L1neg);
	std::vector<Hist *>&& hHVec_L1sig = Hist::Project(Hist::kProjY, hH2D_L1sig);
	std::vector<Hist *>&& hHVec_L1bkg = Hist::Project(Hist::kProjY, hH2D_L1bkg);
	
	Hist * hH2D_L9pos = Hist::New((TH1*)fIss->Get(CStrFmt("hHp_CCestL9%s", vtme.c_str())));
	Hist * hH2D_L9neg = Hist::New((TH1*)fIss->Get(CStrFmt("hHn_CCestL9%s", vtme.c_str())));
	Hist * hH2D_L9sig = Hist::New("hHs_CCestL9", (TH1*)fIss->Get(CStrFmt("hHp_CCestL9%s", vtme.c_str())));
	Hist * hH2D_L9bkg = Hist::New("hHb_CCestL9", (TH1*)fPrL1o9->Get("hHn_CCestL9"));
	
	std::vector<Hist *>&& hHVec_L9pos = Hist::Project(Hist::kProjY, hH2D_L9pos);
	std::vector<Hist *>&& hHVec_L9neg = Hist::Project(Hist::kProjY, hH2D_L9neg);
	std::vector<Hist *>&& hHVec_L9sig = Hist::Project(Hist::kProjY, hH2D_L9sig);
	std::vector<Hist *>&& hHVec_L9bkg = Hist::Project(Hist::kProjY, hH2D_L9bkg);
	
	Hist * hH2D_FSpos = Hist::New((TH1*)fIss->Get(CStrFmt("hHp_CCestFs%s", vtme.c_str())));
	Hist * hH2D_FSneg = Hist::New((TH1*)fIss->Get(CStrFmt("hHn_CCestFs%s", vtme.c_str())));
	Hist * hH2D_FSsig = Hist::New("hHs_CCestFs", (TH1*)fIss->Get(CStrFmt("hHp_CCestFs%s", vtme.c_str())));
	Hist * hH2D_FSbkg = Hist::New("hHb_CCestFs", (TH1*)fPrL1a9->Get("hHn_CCestFs"));
		
	std::vector<Hist *>&& hHVec_FSpos = Hist::Project(Hist::kProjY, hH2D_FSpos);
	std::vector<Hist *>&& hHVec_FSneg = Hist::Project(Hist::kProjY, hH2D_FSneg);
	std::vector<Hist *>&& hHVec_FSsig = Hist::Project(Hist::kProjY, hH2D_FSsig);
	std::vector<Hist *>&& hHVec_FSbkg = Hist::Project(Hist::kProjY, hH2D_FSbkg);
	
	Hist * hH_C090 = Hist::New("hHp_C090_MCEvt", (TH1*)fPrC090->Get("hHp_MCEvt"));
	Hist * hH_C110 = Hist::New("hHp_C110_MCEvt", (TH1*)fPrC110->Get("hHp_MCEvt"));
	
	Hist * hH_ACCPR = Hist::New("hHp_ACCPR_MCEvt", (TH1*)fPr1800->Get("hHp_MCEvt"));
	Hist * hH_ACCAP = Hist::New("hHn_ACCAP_MCEvt", (TH1*)fAp1800->Get("hHn_MCEvt"));
	
	Hist * hHp_Cutflow = Hist::New("hHp_Cutflow", (TH1*)fIss->Get(CStrFmt("hHp_Cutflow%s", vtme.c_str())));
	Hist * hHn_Cutflow = Hist::New("hHn_Cutflow", (TH1*)fIss->Get(CStrFmt("hHn_Cutflow%s", vtme.c_str())));

	Axis AXcc(hH2D_FSsig, Axis::kY);

	for (Int_t ib = hBin[0]; ib <= hBin[1]; ++ib) {
		std::cout << CStrFmt("H %02d FIT\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		
		Double_t numOfPr     = 0.0;
		Double_t errOfPrStat = 0.0;
		Double_t numOfAp     = 0.0;
		Double_t errOfApStat = 0.0;
		Double_t errOfApTemp = 0.0;
		
		Hist * hH_L1pos = hHVec_L1pos.at(ib);
		Hist * hH_L1neg = hHVec_L1neg.at(ib);
		Hist * hH_L1sig = hHVec_L1sig.at(ib);
		Hist * hH_L1bkg = hHVec_L1bkg.at(ib);

		Double_t numOfL1Pr     = (*hH_L1pos)()->Integral();
		Double_t errOfL1PrStat = std::sqrt(numOfL1Pr);
		Double_t numOfL1Smp    = (*hH_L1neg)()->Integral();
		//if (false && numOfL1Smp > 20. && cen < 147.) {
		if (numOfL1Smp > 20. && cen < 147.) {
			Fit::RooParam param("ccest", hH_L1neg, Fit::RooParam::VLIST({ hH_L1bkg, hH_L1sig }));
			Fit::RooSysResult rlt(param, true, 100);
			if (rlt.valid() && rlt.sysValid()) {
				Double_t numOfL1Ap     = rlt.val(1);
				Double_t errOfL1ApStat = rlt.err(1);
				Double_t errOfL1ApTemp = rlt.sysErr(1);
			
				numOfPr += numOfL1Pr;
				errOfPrStat = std::sqrt(errOfPrStat * errOfPrStat + errOfL1PrStat * errOfL1PrStat);
				numOfAp += numOfL1Ap;
				errOfApStat = std::sqrt(errOfApStat * errOfApStat + errOfL1ApStat * errOfL1ApStat);
				errOfApTemp = std::sqrt(errOfApTemp * errOfApTemp + errOfL1ApTemp * errOfL1ApTemp);

				editor.create(StrFmt("H%02dL1", ib));
				editor.cd(1, Canvas::AxisScl_t(0, 1));
				rlt.samp() ->setStyle(Area(), Line(kBlack), Marker(kBlack));
				rlt.sumt() ->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
				rlt.temp(0)->setStyle(Area(), Line(kBlue), Marker(kBlue));
				rlt.temp(1)->setStyle(Area(), Line(kRed), Marker(kRed));
				rlt.samp()->draw("pe");
				Hist::Collect(StrFmt("fit%03d", ib), "", Hist::VLIST({ rlt.sumt(), rlt.temp(0), rlt.temp(1) }))->Draw("same nostack hist");
				Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.85, 0.92, 32));
				Text::Draw(Text::Txt_t(StrFmt("N_{Background} ( %.2f +- stat %.2f sys %.2f )", rlt.val(0), rlt.err(0), rlt.sysErr(0)), kBlue, 0.03), Text::Att_t(0.15, 0.85, 12));
				Text::Draw(Text::Txt_t(StrFmt("N_{Antiproton} ( %.2f +- stat %.2f sys %.2f )", rlt.val(1), rlt.err(1), rlt.sysErr(1)), kRed, 0.03), Text::Att_t(0.15, 0.81, 12));
				Text::Draw(Text::Txt_t(StrFmt("L1 #bar{p}/p Ratio ( %E )", (numOfL1Ap/numOfL1Pr)), kBlack, 0.03), Text::Att_t(0.15, 0.77, 12));
				Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.15, 0.73, 12));
				editor.save();
			}
		}

		Hist * hH_L9pos = hHVec_L9pos.at(ib);
		Hist * hH_L9neg = hHVec_L9neg.at(ib);
		Hist * hH_L9sig = hHVec_L9sig.at(ib);
		Hist * hH_L9bkg = hHVec_L9bkg.at(ib);

		Double_t numOfL9Pr     = (*hH_L9pos)()->Integral();
		Double_t errOfL9PrStat = std::sqrt(numOfL9Pr);
		Double_t numOfL9Smp    = (*hH_L9neg)()->Integral();
		//if (false && numOfL9Smp > 20. && cen < 175.) {
		if (numOfL9Smp > 20. && cen < 175.) {
			Fit::RooParam param("ccest", hH_L9neg, Fit::RooParam::VLIST({ hH_L9bkg, hH_L9sig }));
			Fit::RooSysResult rlt(param, true, 100);
			if (rlt.valid() && rlt.sysValid()) {
				Double_t numOfL9Ap     = rlt.val(1);
				Double_t errOfL9ApStat = rlt.err(1);
				Double_t errOfL9ApTemp = rlt.sysErr(1);
				
				numOfPr += numOfL9Pr;
				errOfPrStat = std::sqrt(errOfPrStat * errOfPrStat + errOfL9PrStat * errOfL9PrStat);
				numOfAp += numOfL9Ap;
				errOfApStat = std::sqrt(errOfApStat * errOfApStat + errOfL9ApStat * errOfL9ApStat);
				errOfApTemp = std::sqrt(errOfApTemp * errOfApTemp + errOfL9ApTemp * errOfL9ApTemp);
				
				editor.create(StrFmt("H%02dL9", ib));
				editor.cd(1, Canvas::AxisScl_t(0, 1));
				rlt.samp() ->setStyle(Area(), Line(kBlack), Marker(kBlack));
				rlt.sumt() ->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
				rlt.temp(0)->setStyle(Area(), Line(kBlue), Marker(kBlue));
				rlt.temp(1)->setStyle(Area(), Line(kRed), Marker(kRed));
				rlt.samp()->draw("pe");
				Hist::Collect(StrFmt("fit%03d", ib), "", Hist::VLIST({ rlt.sumt(), rlt.temp(0), rlt.temp(1) }))->Draw("same nostack hist");
				Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.85, 0.92, 32));
				Text::Draw(Text::Txt_t(StrFmt("N_{Background} ( %.2f +- stat %.2f sys %.2f )", rlt.val(0), rlt.err(0), rlt.sysErr(0)), kBlue, 0.03), Text::Att_t(0.15, 0.85, 12));
				Text::Draw(Text::Txt_t(StrFmt("N_{Antiproton} ( %.2f +- stat %.2f sys %.2f )", rlt.val(1), rlt.err(1), rlt.sysErr(1)), kRed, 0.03), Text::Att_t(0.15, 0.81, 12));
				Text::Draw(Text::Txt_t(StrFmt("L9 #bar{p}/p Ratio ( %E )", (numOfL9Ap/numOfL9Pr)), kBlack, 0.03), Text::Att_t(0.15, 0.77, 12));
				Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.15, 0.73, 12));
				editor.save();
			}
		}

		Hist * hH_FSpos = hHVec_FSpos.at(ib);
		Hist * hH_FSneg = hHVec_FSneg.at(ib);
		Hist * hH_FSsig = hHVec_FSsig.at(ib);
		Hist * hH_FSbkg = hHVec_FSbkg.at(ib);

		Double_t numOfFSPr     = (*hH_FSpos)()->Integral();
		Double_t errOfFSPrStat = std::sqrt(numOfFSPr);
		Double_t numOfFSSmp    = (*hH_FSneg)()->Integral();
		//if (false && numOfFSSmp > 20.) {
		if (numOfFSSmp > 20.) {
			Fit::RooParam param("ccest", hH_FSneg, Fit::RooParam::VLIST({ hH_FSbkg, hH_FSsig }));
			Fit::RooSysResult rlt(param, true, 100);
			if (rlt.valid() && rlt.sysValid()) {
				Double_t numOfFSAp     = rlt.val(1);
				Double_t errOfFSApStat = rlt.err(1);
				Double_t errOfFSApTemp = rlt.sysErr(1);
				
				numOfPr += numOfFSPr;
				errOfPrStat = std::sqrt(errOfPrStat * errOfPrStat + errOfFSPrStat * errOfFSPrStat);
				numOfAp += numOfFSAp;
				errOfApStat = std::sqrt(errOfApStat * errOfApStat + errOfFSApStat * errOfFSApStat);
				errOfApTemp = std::sqrt(errOfApTemp * errOfApTemp + errOfFSApTemp * errOfFSApTemp);
				
				editor.create(StrFmt("H%02dFS", ib));
				editor.cd(1, Canvas::AxisScl_t(0, 1));
				rlt.samp() ->setStyle(Area(), Line(kBlack), Marker(kBlack));
				rlt.sumt() ->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
				rlt.temp(0)->setStyle(Area(), Line(kBlue), Marker(kBlue));
				rlt.temp(1)->setStyle(Area(), Line(kRed), Marker(kRed));
				rlt.samp()->draw("pe");
				Hist::Collect(StrFmt("fit%03d", ib), "", Hist::VLIST({ rlt.sumt(), rlt.temp(0), rlt.temp(1) }))->Draw("same nostack hist");
				Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.85, 0.92, 32));
				Text::Draw(Text::Txt_t(StrFmt("N_{Background} ( %.2f +- stat %.2f sys %.2f )", rlt.val(0), rlt.err(0), rlt.sysErr(0)), kBlue, 0.03), Text::Att_t(0.15, 0.85, 12));
				Text::Draw(Text::Txt_t(StrFmt("N_{Antiproton} ( %.2f +- stat %.2f sys %.2f )", rlt.val(1), rlt.err(1), rlt.sysErr(1)), kRed, 0.03), Text::Att_t(0.15, 0.81, 12));
				Text::Draw(Text::Txt_t(StrFmt("FS #bar{p}/p Ratio ( %E )", (numOfFSAp/numOfFSPr)), kBlack, 0.03), Text::Att_t(0.15, 0.77, 12));
				Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.15, 0.73, 12));
				editor.save();
			}
		}

		if (MgntNum::EqualToZero(numOfPr) || MgntNum::EqualToZero(numOfAp)) continue;

		//Double_t crrOfAppAccp = 3.16258e-01 * TMath::Power(cen, -6.30666e-01) + 1.03140e+00;
		Double_t crrOfAppAccp = 1.02394e+00 + 2.77708e-01 * TMath::Power(cen, -5.21422e-01) + 1.01791e-01 * TMath::Power(cen, -3.07641e+00);
		Double_t valOfApp     = crrOfAppAccp * numOfAp / numOfPr; 
		Double_t errOfAppStat = valOfApp * std::sqrt(
		                        errOfPrStat*errOfPrStat/numOfPr/numOfPr + 
														errOfApStat*errOfApStat/numOfAp/numOfAp
														);
		Double_t errOfAppTemp = valOfApp * std::sqrt(
														errOfApTemp*errOfApTemp/numOfAp/numOfAp
														);
		Double_t errOfAppAccp = valOfApp * 0.015;
		Double_t errOfAppXsec = valOfApp * 0.03;
		Double_t errOfAppSysm = std::sqrt(
														errOfAppTemp*errOfAppTemp+
														errOfAppAccp*errOfAppAccp+
														errOfAppXsec*errOfAppXsec
														);
		Double_t errOfAppTotl = std::sqrt(
		                        errOfAppStat*errOfAppStat+
														errOfAppSysm*errOfAppSysm
														);

		HAppnum_pr->pushPoint(Graph::Point_t(cen, numOfPr));
		HAppnum_ap->pushPoint(Graph::Point_t(cen, numOfAp));
		HApprlt_stat->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppStat));
		HApprlt_totl->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
		//HAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		HApperr_totl->pushPoint(Graph::Point_t(cen, 100. * errOfAppTotl/valOfApp));
		HApperr_stat->pushPoint(Graph::Point_t(cen, 100. * errOfAppStat/valOfApp));
		HApperr_sysm->pushPoint(Graph::Point_t(cen, 100. * errOfAppSysm/valOfApp));
		HApperr_temp->pushPoint(Graph::Point_t(cen, 100. * errOfAppTemp/valOfApp));
		//HApperr_accp->pushPoint(Graph::Point_t(cen, 100. * errOfAppAccp/valOfApp));
		//HApperr_xsec->pushPoint(Graph::Point_t(cen, 100. * errOfAppXsec/valOfApp));

		if (ib > cbBin[2] && ib <= cbBin[3]) {
			Apprlt_stat->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppStat));
			Apprlt_totl->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
		}
	}
	
	for (Int_t ib = 1; ib <= AXnr.nbin(); ++ib) {
		std::cout << CStrFmt("H %02d XSEC\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		if (AXnr.bins(ib-1) >= 1.0 && AXnr.bins(ib) <= 200.) {
			CdfFlux->SetParameters(trgC090, 1.0, 200.0);
			Double_t numOfSel_C090 = (*hH_C090)()->GetBinContent(ib);
			Double_t numOfTrg_C090 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratC090 = numOfSel_C090 / numOfTrg_C090;
			CdfFlux->SetParameters(trgC110, 1.0, 200.0);
			Double_t numOfSel_C110 = (*hH_C110)()->GetBinContent(ib);
			Double_t numOfTrg_C110 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratC110 = numOfSel_C110 / numOfTrg_C110;
			
			if (numOfSel_C090 < 100 || numOfSel_C110 < 100) continue;

			Double_t raterr  = std::sqrt(ratC090/numOfTrg_C090 + ratC110/numOfTrg_C110);
			Double_t ratmius = std::fabs(ratC090 - ratC110);
			Double_t ratplus = std::fabs(ratC090 + ratC110);

			Double_t errOfAppXsec = ratmius / ratplus;
			Double_t errOfAppXsecErr = errOfAppXsec * raterr * std::sqrt((1./ratmius/ratmius + 1./ratplus/ratplus));
			if ((errOfAppXsecErr / errOfAppXsec) > 0.5 || !MgntNum::Valid(errOfAppXsec)) continue;
			HApperr_xsec->pushPointWithError(Graph::Point_t(cen, 100. * errOfAppXsec), Graph::Error_t(0.0, 100. * errOfAppXsecErr));
		}
	}

	for (Int_t ib = 1; ib <= AXnr.nbin(); ++ib) {
		std::cout << CStrFmt("H %02d ACCP\n", ib);
		Double_t cen = AXnr.center(ib, Axis::kLog);
		if (AXnr.bins(ib-1) >= 1.0 && AXnr.bins(ib) <= 800.) {
			CdfFlux->SetParameters(trgPr1800, 1.0, 800.0);
			Double_t numOfSel_ACCPR = (*hH_ACCPR)()->GetBinContent(ib);
			Double_t numOfTrg_ACCPR = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratACCPR = numOfSel_ACCPR / numOfTrg_ACCPR;
			CdfFlux->SetParameters(trgAp1800, 1.0, 800.0);
			Double_t numOfSel_ACCAP = (*hH_ACCAP)()->GetBinContent(ib);
			Double_t numOfTrg_ACCAP = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
			Double_t ratACCAP = numOfSel_ACCAP / numOfTrg_ACCAP;
			
			if (numOfSel_ACCPR < 100 || numOfSel_ACCAP < 100) continue;
			
			Double_t crrOfAppAccp = ratACCPR / ratACCAP;
			Double_t crrOfAppAccpErr = crrOfAppAccp * std::sqrt(1./numOfSel_ACCPR+1./numOfSel_ACCAP);
			if ((crrOfAppAccpErr / crrOfAppAccp) > 0.02 || !MgntNum::Valid(crrOfAppAccp)) continue;
			HAppcrr_accp->pushPointWithError(Graph::Point_t(cen, crrOfAppAccp), Graph::Error_t(0.0,  crrOfAppAccpErr));
			
			//Double_t crrOfAppAccpPred = 3.16258e-01 * TMath::Power(cen, -6.30666e-01) + 1.03140e+00;
			Double_t crrOfAppAccpPred = 1.02394e+00 + 2.77708e-01 * TMath::Power(cen, -5.21422e-01) + 1.01791e-01 * TMath::Power(cen, -3.07641e+00);
			Double_t errOfAppAccp = (crrOfAppAccp - crrOfAppAccpPred) / crrOfAppAccpPred;
			HApperr_accp->pushPointWithError(Graph::Point_t(cen, 100. * errOfAppAccp), Graph::Error_t(0.0, 100. * crrOfAppAccpErr / crrOfAppAccpPred));
		}	
	}

	Int_t  nHCut = (*hHp_Cutflow)()->GetNbinsY()-1;
	Axis AXHcut("cut", nHCut, 0.0, Double_t(nHCut));
	Hist * hHp_Eff = Hist::New("hHp_Eff", "COMp Efficiency Rig", AXnr, AXHcut);
	for (Int_t irig = 1; irig <= (*hHp_Eff)()->GetNbinsX(); ++irig) {
		for (Int_t icut = 1; icut <= (*hHp_Eff)()->GetNbinsY(); ++icut) {
			Double_t pw  = (*hHp_Cutflow)()->GetBinContent(irig, icut + 1);
			Double_t tw  = (*hHp_Cutflow)()->GetBinContent(irig, icut);
			if (MgntNum::EqualToZero(pw) || MgntNum::EqualToZero(tw)) continue; 
			Double_t eff = pw / tw;
			Double_t sgm = std::sqrt(eff * (1. - eff) / tw);
			Double_t eru = ((eff + sgm) > 1.0) ? (1.0 - eff) : sgm;
			Double_t erl = ((eff - sgm) > 0.0) ? eff         : sgm;
			(*hHp_Eff)()->SetBinContent(irig, icut, eff);
			(*hHp_Eff)()->SetBinError  (irig, icut, sgm);
		}
	}
	Hist * hHn_Eff = Hist::New("hHn_Eff", "COMn Efficiency Rig", AXnr, AXHcut);
	for (Int_t irig = 1; irig <= (*hHn_Eff)()->GetNbinsX(); ++irig) {
		for (Int_t icut = 1; icut <= (*hHn_Eff)()->GetNbinsY(); ++icut) {
			Double_t pw  = (*hHn_Cutflow)()->GetBinContent(irig, icut + 1);
			Double_t tw  = (*hHn_Cutflow)()->GetBinContent(irig, icut);
			if (MgntNum::EqualToZero(pw) || MgntNum::EqualToZero(tw)) continue; 
			Double_t eff = pw / tw;
			Double_t sgm = std::sqrt(eff * (1. - eff) / tw);
			Double_t eru = ((eff + sgm) > 1.0) ? (1.0 - eff) : sgm;
			Double_t erl = ((eff - sgm) > 0.0) ? eff         : sgm;
			(*hHn_Eff)()->SetBinContent(irig, icut, eff);
			(*hHn_Eff)()->SetBinError  (irig, icut, sgm);
		}
	}


	// Result
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	LAppnum_pr->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IAppnum_pr->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HAppnum_pr->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CAppnum_pr", ";|Rigidity| [GV];N_{p}", Graph::VLIST({ LAppnum_pr, IAppnum_pr, HAppnum_pr }))->Draw("ap");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	LAppnum_ap->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IAppnum_ap->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HAppnum_ap->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CAppnum_ap", ";|Rigidity| [GV];N_{#bar{p}}", Graph::VLIST({ LAppnum_ap, IAppnum_ap, HAppnum_ap }))->Draw("ap");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	LApprlt_stat->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApprlt_stat->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApprlt_stat->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CApprlt_stat", ";|Rigidity| [GV];#bar{p}/p Ratio", Graph::VLIST({ LApprlt_stat, IApprlt_stat, HApprlt_stat }))->Draw("ap");
	//Apprlt_Pamela->Draw("same p");
	Apprlt_AMSPub->Draw("same p");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	LApprlt_totl->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApprlt_totl->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApprlt_totl->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CApprlt_totl", ";|Rigidity| [GV];#bar{p}/p Ratio", Graph::VLIST({ LApprlt_totl, IApprlt_totl, HApprlt_totl }))->Draw("ap");
	//Apprlt_Pamela->Draw("same p");
	Apprlt_AMSPub->Draw("same p");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	LAppcrr_accp->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IAppcrr_accp->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HAppcrr_accp->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CAppcrr_accp", ";|Rigidity| [GV];Acceptance Correction", Graph::VLIST({ LAppcrr_accp, IAppcrr_accp, HAppcrr_accp }))->Draw("ap");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	LApperr_stat->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApperr_stat->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApperr_stat->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CApperr_stat", ";|Rigidity| [GV];Statistics Error [%]", Graph::VLIST({ LApperr_stat, IApperr_stat, HApperr_stat }))->Draw("ap");
	editor.save();
	
	editor.create();
	LApperr_temp->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApperr_temp->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApperr_temp->setStyle(Area(), Line(kRed), Marker(kRed));
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	Graph::Collect("CApperr_temp", ";|Rigidity| [GV];Template Error [%]", Graph::VLIST({ LApperr_temp, IApperr_temp, HApperr_temp }))->Draw("ap");
	editor.save();
	
	editor.create();
	LApperr_accp->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApperr_accp->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApperr_accp->setStyle(Area(), Line(kRed), Marker(kRed));
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	Graph::Collect("CApperr_accp", ";|Rigidity| [GV];Acceptance Error [%]", Graph::VLIST({ LApperr_accp, IApperr_accp, HApperr_accp }))->Draw("ap");
	editor.save();
	
	editor.create();
	LApperr_xsec->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApperr_xsec->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApperr_xsec->setStyle(Area(), Line(kRed), Marker(kRed));
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	Graph::Collect("CApperr_xsec", ";|Rigidity| [GV];Cross Section Error [%]", Graph::VLIST({ LApperr_xsec, IApperr_xsec, HApperr_xsec }))->Draw("ap");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	LApperr_totl->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApperr_totl->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApperr_totl->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CApperr_totl", ";|Rigidity| [GV];Total Error [%]", Graph::VLIST({ LApperr_totl, IApperr_totl, HApperr_totl }))->Draw("ap");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	Apprlt_stat->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("FApprlt_stat", ";|Rigidity| [GV];#bar{p}/p Ratio", Graph::VLIST({ Apprlt_stat }))->Draw("ap");
	//Apprlt_Pamela->Draw("same p");
	Apprlt_AMSPub->Draw("same p");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(1, 0));
	Apprlt_totl->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("FApprlt_totl", ";|Rigidity| [GV];#bar{p}/p Ratio", Graph::VLIST({ Apprlt_totl }))->Draw("ap");
	//Apprlt_Pamela->Draw("same p");
	Apprlt_AMSPub->Draw("same p");
	editor.save();
	
	editor.close();

	TFile * file = new TFile(CStrFmt("%s/YiAna%s.root", goddir.c_str(), vtme.c_str()), "RECREATE");
	hLp_Cutflow->write();
	hLn_Cutflow->write();
	hLp_Eff->write();
	hLn_Eff->write();
	
	hIp_Cutflow->write();
	hIn_Cutflow->write();
	hIp_Eff->write();
	hIn_Eff->write();
	
	hHp_Cutflow->write();
	hHn_Cutflow->write();
	hHp_Eff->write();
	hHn_Eff->write();

	Graph::Write();
	file->Write();
	file->Close();

	return 1;
}

#endif // __YiFinalAna_C__
