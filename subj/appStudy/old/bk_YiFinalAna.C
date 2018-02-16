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

	std::string tagdir = "/afs/cern.ch/work/h/hchou/BSUB_ANALYSIS/antipp5/NOW__";

	// ISS
	TFile * fIss = TFile::Open(CStrFmt("%sISS/YiAnalytics.root", tagdir.c_str()));

	// Charge Confusion
	TFile * fPrL1a9 = TFile::Open(CStrFmt("%sL1a9/YiAnalytics.root", tagdir.c_str()));
	TFile * fPrL1o9 = TFile::Open(CStrFmt("%sL1o9/YiAnalytics.root", tagdir.c_str()));

	// Cross Section
	TFile * fPrC090 = TFile::Open(CStrFmt("%sC090/YiAnalytics.root", tagdir.c_str()));
	TTree * tPrC090 = (TTree *) fPrC090->Get("DSTRun");
	UInt_t  trgC090 = 0;
	{ UInt_t trg = 0; tPrC090->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tPrC090->GetEntries(); ++it)
		{ tPrC090->GetEntry(it); trgC090 += trg; } }

	TFile * fPrC110 = TFile::Open(CStrFmt("%sC110/YiAnalytics.root", tagdir.c_str()));
	TTree * tPrC110 = (TTree *) fPrC110->Get("DSTRun");
	UInt_t  trgC110 = 0;
	{ UInt_t trg = 0; tPrC110->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tPrC110->GetEntries(); ++it)
		{ tPrC110->GetEntry(it); trgC110 += trg; } }

	// Unfolding
	TFile * fPr0510 = TFile::Open(CStrFmt("%sPr0_510/YiAnalytics.root", tagdir.c_str()));
	TTree * tPr0510 = (TTree *) fPr0510->Get("DSTRun");
	UInt_t  trgPr0510 = 0;
	{ UInt_t trg = 0; tPr0510->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tPr0510->GetEntries(); ++it)
		{ tPr0510->GetEntry(it); trgPr0510 += trg; } }
	
	TFile * fPr1800 = TFile::Open(CStrFmt("%sPr1800/YiAnalytics.root", tagdir.c_str()));
	TTree * tPr1800 = (TTree *) fPr1800->Get("DSTRun");
	UInt_t  trgPr1800 = 0;
	{ UInt_t trg = 0; tPr1800->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tPr1800->GetEntries(); ++it)
		{ tPr1800->GetEntry(it); trgPr1800 += trg; } }
	
	TFile * fAp0510 = TFile::Open(CStrFmt("%sAp0_510/YiAnalytics.root", tagdir.c_str()));
	TTree * tAp0510 = (TTree *) fAp0510->Get("DSTRun");
	UInt_t  trgAp0510 = 0;
	{ UInt_t trg = 0; tAp0510->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tAp0510->GetEntries(); ++it)
		{ tAp0510->GetEntry(it); trgAp0510 += trg; } }
	
	TFile * fAp1800 = TFile::Open(CStrFmt("%sAp1800/YiAnalytics.root", tagdir.c_str()));
	TTree * tAp1800 = (TTree *) fAp1800->Get("DSTRun");
	UInt_t  trgAp1800 = 0;
	{ UInt_t trg = 0; tAp1800->SetBranchAddress("trgEV", &trg);
		for (Long64_t it = 0; it < tAp1800->GetEntries(); ++it)
		{ tAp1800->GetEntry(it); trgAp1800 += trg; } }

	// Binning
	Int_t lBin[2] = {  1, 16 };
	Int_t iBin[2] = { 10, 30 };
	Int_t hBin[2] = { 35, 57 };
	Int_t combineBin[3] = { 1, 16, 57 };

	TFile * fAppBin = TFile::Open("/afs/cern.ch/user/h/hchou/public/DATABASE/physics/binning/antipp3_rig.root");
	Axis AXnr((TH1*)fAppBin->Get("hbin0"), Axis::kX);
  Axis AXir((TH1*)fAppBin->Get("hbin1"), Axis::kX);
	fAppBin->Close();
	
	PdfEditor editor(PdfEditor::kSlice, StrFmt("YiAna%s", vtme.c_str()));
	
	// Antiproton-to-Proton Flux Ratio
	Graph * Appnum_pr   = Graph::New("apprlt_pr"  , "", "|Rigidity| [GV]", "N_{p}");
	Graph * Appnum_ap   = Graph::New("apprlt_ap"  , "", "|Rigidity| [GV]", "N_{#bar{p}}");
	Graph * Apprlt_stat = Graph::New("apprlt_stat", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * Apprlt_full = Graph::New("apprlt_full", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * Appcrr_accp = Graph::New("appcrr_accp", "", "|Rigidity| [GV]", "Acceptance Correction");
	Graph * Apperr_totl = Graph::New("apperr_totl", "", "|Rigidity| [GV]", "Total Error [%]");
	Graph * Apperr_stat = Graph::New("apperr_stat", "", "|Rigidity| [GV]", "Statistics Error [%]");
	Graph * Apperr_sysm = Graph::New("apperr_sysm", "", "|Rigidity| [GV]", "Systematic Error [%]");
	Graph * Apperr_temp = Graph::New("apperr_temp", "", "|Rigidity| [GV]", "Template Error [%]");
	Graph * Apperr_accp = Graph::New("apperr_accp", "", "|Rigidity| [GV]", "Acceptance Error [%]");
	Graph * Apperr_xsec = Graph::New("apperr_xsec", "", "|Rigidity| [GV]", "Cross-Section Error [%]");
	
	//----  Low Energy  ----//
	Graph * LAppnum_pr   = Graph::New("lapprlt_pr"  , "", "|Rigidity| [GV]", "N_{p}");
	Graph * LAppnum_ap   = Graph::New("lapprlt_ap"  , "", "|Rigidity| [GV]", "N_{#bar{p}}");
	Graph * LApprlt_stat = Graph::New("lapprlt_stat", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * LApprlt_full = Graph::New("lapprlt_full", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * LAppcrr_accp = Graph::New("lappcrr_accp", "", "|Rigidity| [GV]", "Acceptance Correction");
	Graph * LApperr_totl = Graph::New("lapperr_totl", "", "|Rigidity| [GV]", "Total Error [%]");
	Graph * LApperr_stat = Graph::New("lapperr_stat", "", "|Rigidity| [GV]", "Statistics Error [%]");
	Graph * LApperr_sysm = Graph::New("lapperr_sysm", "", "|Rigidity| [GV]", "Systematic Error [%]");
	Graph * LApperr_temp = Graph::New("lapperr_temp", "", "|Rigidity| [GV]", "Template Error [%]");
	Graph * LApperr_accp = Graph::New("lapperr_accp", "", "|Rigidity| [GV]", "Acceptance Error [%]");
	Graph * LApperr_xsec = Graph::New("lapperr_xsec", "", "|Rigidity| [GV]", "Cross-Section Error [%]");
	
	Hist * hL2D_pos = Hist::New((TH1*)fIss->Get(CStrFmt("hLp_Tofb%s", vtme.c_str())));
	Hist * hL2D_neg = Hist::New((TH1*)fIss->Get(CStrFmt("hLn_Tofb%s", vtme.c_str())));
	Hist * hL2D_sig = Hist::New((TH1*)fIss->Get(CStrFmt("hLs_Tofb%s", vtme.c_str())));
	Hist * hL2D_bkg = Hist::New((TH1*)fIss->Get(CStrFmt("hLb_Tofb%s", vtme.c_str())));
		
	std::vector<Hist *>&& hLVec_pos = Hist::Project(Hist::kProjY, hL2D_pos);
	std::vector<Hist *>&& hLVec_neg = Hist::Project(Hist::kProjY, hL2D_neg);
	std::vector<Hist *>&& hLVec_sig = Hist::Project(Hist::kProjY, hL2D_sig);
	std::vector<Hist *>&& hLVec_bkg = Hist::Project(Hist::kProjY, hL2D_bkg);

	Hist * hL_C090 = Hist::New("hLp_C090_Evt", (TH1*)fPrC090->Get("hLp_Evt"));
	Hist * hL_C110 = Hist::New("hLn_C110_Evt", (TH1*)fPrC110->Get("hLp_Evt"));
	
	Hist * hL_ACCPR = Hist::New("hLp_ACCPR_MCEvt", (TH1*)fPr1800->Get("hLp_MCEvt"));
	Hist * hL_ACCAP = Hist::New("hLn_ACCAP_MCEvt", (TH1*)fAp1800->Get("hLn_MCEvt"));

	for (Int_t ib = lBin[0]; ib <= lBin[1]; ++ib) {
		Double_t cen = AXnr.center(ib, Axis::kLog);
		Hist * hL_pos = hLVec_pos.at(ib);
		Hist * hL_neg = hLVec_neg.at(ib);
		Hist * hL_sig = hLVec_sig.at(ib);
		Hist * hL_bkg = hLVec_bkg.at(ib);

		Double_t numOfPr     = (*hL_pos)()->Integral();
		Double_t errOfPrStat = std::sqrt(numOfPr);
		Double_t errOfPr     = errOfPrStat;
		
		Double_t numOfSmp = (*hL_neg)()->Integral();
		if (numOfSmp < 20.) continue;

		Fit::RooParam param("bta", hL_neg, Fit::RooParam::VLIST({ hL_bkg, hL_sig }));
		Fit::RooSysResult rlt(param);
		if (!rlt.valid()) continue;
		if (!rlt.sysValid()) continue;
		Double_t numOfAp     = rlt.val(1);
		Double_t errOfApStat = rlt.err(1);
		Double_t errOfApTemp = rlt.sysErr(1);
		Double_t errOfAp     = std::sqrt(errOfApStat * errOfApStat + errOfApTemp * errOfApTemp);

		Double_t crrOfAppAccp = 0.362262 * std::exp(-0.857079 * std::log(cen)) + 1.06241;
		Double_t valOfApp     = crrOfAppAccp * numOfAp / numOfPr; 
		Double_t errOfAppStat = valOfApp * std::sqrt(
		                        errOfPrStat*errOfPrStat/numOfPr/numOfPr + 
														errOfApStat*errOfApStat/numOfAp/numOfAp
														);
		Double_t errOfAppTemp = valOfApp * std::sqrt(
														errOfApTemp*errOfApTemp/numOfAp/numOfAp
														);
		Double_t errOfAppAccp = valOfApp * 0.015;
		Double_t errOfAppXsec = valOfApp * 0.04;
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
		LApprlt_full->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
		LAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		LApperr_totl->pushPoint(Graph::Point_t(cen, 100. * errOfAppTotl/valOfApp));
		LApperr_stat->pushPoint(Graph::Point_t(cen, 100. * errOfAppStat/valOfApp));
		LApperr_sysm->pushPoint(Graph::Point_t(cen, 100. * errOfAppSysm/valOfApp));
		LApperr_temp->pushPoint(Graph::Point_t(cen, 100. * errOfAppTemp/valOfApp));
		LApperr_accp->pushPoint(Graph::Point_t(cen, 100. * errOfAppAccp/valOfApp));
		LApperr_xsec->pushPoint(Graph::Point_t(cen, 100. * errOfAppXsec/valOfApp));

		editor.create(StrFmt("L%02d", ib));
		editor.cd(1, Canvas::AxisScl_t(0, 0));
		rlt.samp() ->setStyle(Area(), Line(kBlack), Marker(kBlack));
		rlt.sumt() ->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
		rlt.temp(0)->setStyle(Area(), Line(kBlue), Marker(kBlue));
		rlt.temp(1)->setStyle(Area(), Line(kRed), Marker(kRed));
		Hist::Collect(StrFmt("fit%03d", ib), "", Hist::VLIST({ rlt.samp(), rlt.sumt(), rlt.temp(0), rlt.temp(1) }))->Draw("nostack hist");
		Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GeV ~ %.2f GeV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.85, 0.92, 32));
		Text::Draw(Text::Txt_t(StrFmt("N_{Background} ( %.2f +- stat %.2f sys %.2f )", rlt.val(0), rlt.err(0), rlt.sysErr(0)), kBlue, 0.03), Text::Att_t(0.15, 0.85, 12));
		Text::Draw(Text::Txt_t(StrFmt("N_{Antiproton} ( %.2f +- stat %.2f sys %.2f )", rlt.val(1), rlt.err(1), rlt.sysErr(1)), kRed, 0.03), Text::Att_t(0.15, 0.81, 12));
		Text::Draw(Text::Txt_t(StrFmt("#bar{p}/p Ratio ( %E +- stat %E )", valOfApp, errOfAppStat), kBlack, 0.03), Text::Att_t(0.15, 0.77, 12));
		Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.15, 0.73, 12));
		editor.save();

		CdfFlux->SetParameters(trgC090, 1.0, 200.0);
		Double_t numOfSel_C090 = (*hL_C090)()->GetBinContent(ib);
		Double_t numOfTrg_C090 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratC090 = numOfSel_C090 / numOfTrg_C090;
		CdfFlux->SetParameters(trgC110, 1.0, 200.0);
		Double_t numOfSel_C110 = (*hL_C110)()->GetBinContent(ib);
		Double_t numOfTrg_C110 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratC110 = numOfSel_C110 / numOfTrg_C110;
		//Double_t errOfAppXsec = std::sqrt(2.0) * (ratC090 - ratC110) / (ratC110 + ratC090);
		//LApperr_xsec->pushPoint(Graph::Point_t(cen, 100. * errOfAppXsec));
		
		CdfFlux->SetParameters(trgPr1800, 1.0, 800.0);
		Double_t numOfSel_ACCPR = (*hL_ACCPR)()->GetBinContent(ib);
		Double_t numOfTrg_ACCPR = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratACCPR = numOfSel_ACCPR / numOfTrg_ACCPR;
		CdfFlux->SetParameters(trgAp1800, 1.0, 800.0);
		Double_t numOfSel_ACCAP = (*hL_ACCAP)()->GetBinContent(ib);
		Double_t numOfTrg_ACCAP = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratACCAP = numOfSel_ACCAP / numOfTrg_ACCAP;
		//Double_t crrOfAppAccp = ratACCPR / ratACCAP;
		//LAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		//Double_t errOfAppAccp = std::fabs((ratACCPR / ratACCAP) - crrOfAppAccp);
		//LApperr_accp->pushPoint(Graph::Point_t(cen, errOfAppAccp));
	}
	
	//----  Intermedia Energy  ----//
	Graph * IAppnum_pr   = Graph::New("iapprlt_ap"  , "", "|Rigidity| [GV]", "N_{p}");
	Graph * IAppnum_ap   = Graph::New("iapprlt_pr"  , "", "|Rigidity| [GV]", "N_{#bar{p}}");
	Graph * IApprlt_stat = Graph::New("iapprlt_stat", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * IApprlt_full = Graph::New("iapprlt_full", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * IAppcrr_accp = Graph::New("iappcrr_accp", "", "|Rigidity| [GV]", "Acceptance Correction");
	Graph * IApperr_totl = Graph::New("iapperr_totl", "", "|Rigidity| [GV]", "Total Error [%]");
	Graph * IApperr_stat = Graph::New("iapperr_stat", "", "|Rigidity| [GV]", "Statistics Error [%]");
	Graph * IApperr_sysm = Graph::New("iapperr_sysm", "", "|Rigidity| [GV]", "Systematic Error [%]");
	Graph * IApperr_temp = Graph::New("iapperr_temp", "", "|Rigidity| [GV]", "Template Error [%]");
	Graph * IApperr_accp = Graph::New("iapperr_accp", "", "|Rigidity| [GV]", "Acceptance Error [%]");
	Graph * IApperr_xsec = Graph::New("iapperr_xsec", "", "|Rigidity| [GV]", "Cross-Section Error [%]");
	
	Hist * hI2D_pos = Hist::New((TH1*)fIss->Get(CStrFmt("hIp_Trdl%s", vtme.c_str())));
	Hist * hI2D_neg = Hist::New((TH1*)fIss->Get(CStrFmt("hIn_Trdl%s", vtme.c_str())));
	Hist * hI2D_sig = Hist::New((TH1*)fIss->Get(CStrFmt("hIs_Trdl%s", vtme.c_str())));
	Hist * hI2D_bkg = Hist::New((TH1*)fIss->Get(CStrFmt("hIb_Trdl%s", vtme.c_str())));
		
	std::vector<Hist *>&& hIVec_pos = Hist::Project(Hist::kProjY, hI2D_pos);
	std::vector<Hist *>&& hIVec_neg = Hist::Project(Hist::kProjY, hI2D_neg);
	std::vector<Hist *>&& hIVec_sig = Hist::Project(Hist::kProjY, hI2D_sig);
	std::vector<Hist *>&& hIVec_bkg = Hist::Project(Hist::kProjY, hI2D_bkg);
	
	Hist * hI_C090 = Hist::New("hIp_C090_Evt", (TH1*)fPrC090->Get("hIp_Evt"));
	Hist * hI_C110 = Hist::New("hIn_C110_Evt", (TH1*)fPrC110->Get("hIp_Evt"));
	
	Hist * hI_ACCPR = Hist::New("hIp_ACCPR_MCEvt", (TH1*)fPr1800->Get("hIp_MCEvt"));
	Hist * hI_ACCAP = Hist::New("hIn_ACCAP_MCEvt", (TH1*)fAp1800->Get("hIn_MCEvt"));

	for (Int_t ib = iBin[0]; ib <= iBin[1]; ++ib) {
		Double_t cen = AXnr.center(ib, Axis::kLog);
		Hist * hI_pos = hIVec_pos.at(ib);
		Hist * hI_neg = hIVec_neg.at(ib);
		Hist * hI_sig = hIVec_sig.at(ib);
		Hist * hI_bkg = hIVec_bkg.at(ib);

		Double_t numOfPr     = (*hI_pos)()->Integral();
		Double_t errOfPrStat = std::sqrt(numOfPr);
		Double_t errOfPr     = errOfPrStat;

		Double_t numOfSmp = (*hI_neg)()->Integral();
		if (numOfSmp < 20.) continue;

		Fit::RooParam param("bta", hI_neg, Fit::RooParam::VLIST({ hI_bkg, hI_sig }));
		Fit::RooSysResult rlt(param);
		if (!rlt.valid()) continue;
		if (!rlt.sysValid()) continue;
		Double_t numOfAp     = rlt.val(1);
		Double_t errOfApStat = rlt.err(1);
		Double_t errOfApTemp = rlt.sysErr(1);
		Double_t errOfAp     = std::sqrt(errOfApStat * errOfApStat + errOfApTemp * errOfApTemp);

		Double_t crrOfAppAccp = 0.362262 * std::exp(-0.857079 * std::log(cen)) + 1.06241;
		Double_t valOfApp     = crrOfAppAccp * numOfAp / numOfPr; 
		Double_t errOfAppStat = valOfApp * std::sqrt(
		                        errOfPrStat*errOfPrStat/numOfPr/numOfPr + 
														errOfApStat*errOfApStat/numOfAp/numOfAp
														);
		Double_t errOfAppTemp = valOfApp * std::sqrt(
														errOfApTemp*errOfApTemp/numOfAp/numOfAp
														);
		Double_t errOfAppAccp = valOfApp * 0.015;
		Double_t errOfAppXsec = valOfApp * 0.04;
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
		IApprlt_full->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
		IAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		IApperr_totl->pushPoint(Graph::Point_t(cen, 100. * errOfAppTotl/valOfApp));
		IApperr_stat->pushPoint(Graph::Point_t(cen, 100. * errOfAppStat/valOfApp));
		IApperr_sysm->pushPoint(Graph::Point_t(cen, 100. * errOfAppSysm/valOfApp));
		IApperr_temp->pushPoint(Graph::Point_t(cen, 100. * errOfAppTemp/valOfApp));
		IApperr_accp->pushPoint(Graph::Point_t(cen, 100. * errOfAppAccp/valOfApp));
		IApperr_xsec->pushPoint(Graph::Point_t(cen, 100. * errOfAppXsec/valOfApp));

		editor.create(StrFmt("I%02d", ib));
		editor.cd(1, Canvas::AxisScl_t(0, 0));
		rlt.samp() ->setStyle(Area(), Line(kBlack), Marker(kBlack));
		rlt.sumt() ->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
		rlt.temp(0)->setStyle(Area(), Line(kBlue), Marker(kBlue));
		rlt.temp(1)->setStyle(Area(), Line(kRed), Marker(kRed));
		Hist::Collect(StrFmt("fit%03d", ib), "", Hist::VLIST({ rlt.samp(), rlt.sumt(), rlt.temp(0), rlt.temp(1) }))->Draw("nostack hist");
		Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GeV ~ %.2f GeV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.85, 0.92, 32));
		Text::Draw(Text::Txt_t(StrFmt("N_{Background} ( %.2f +- stat %.2f sys %.2f )", rlt.val(0), rlt.err(0), rlt.sysErr(0)), kBlue, 0.03), Text::Att_t(0.15, 0.85, 12));
		Text::Draw(Text::Txt_t(StrFmt("N_{Antiproton} ( %.2f +- stat %.2f sys %.2f )", rlt.val(1), rlt.err(1), rlt.sysErr(1)), kRed, 0.03), Text::Att_t(0.15, 0.81, 12));
		Text::Draw(Text::Txt_t(StrFmt("#bar{p}/p Ratio ( %E +- stat %E )", valOfApp, errOfAppStat), kBlack, 0.03), Text::Att_t(0.15, 0.77, 12));
		Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.15, 0.73, 12));
		editor.save();
		
		CdfFlux->SetParameters(trgC090, 1.0, 200.0);
		Double_t numOfSel_C090 = (*hI_C090)()->GetBinContent(ib);
		Double_t numOfTrg_C090 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratC090 = numOfSel_C090 / numOfTrg_C090;
		CdfFlux->SetParameters(trgC110, 1.0, 200.0);
		Double_t numOfSel_C110 = (*hI_C110)()->GetBinContent(ib);
		Double_t numOfTrg_C110 = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratC110 = numOfSel_C110 / numOfTrg_C110;
		//Double_t errOfAppXsec = std::sqrt(2.0) * (ratC090 - ratC110) / (ratC110 + ratC090);
		//IApperr_xsec->pushPoint(Graph::Point_t(cen, 100. * errOfAppXsec));
		
		CdfFlux->SetParameters(trgPr1800, 1.0, 800.0);
		Double_t numOfSel_ACCPR = (*hI_ACCPR)()->GetBinContent(ib);
		Double_t numOfTrg_ACCPR = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratACCPR = numOfSel_ACCPR / numOfTrg_ACCPR;
		CdfFlux->SetParameters(trgAp1800, 1.0, 800.0);
		Double_t numOfSel_ACCAP = (*hI_ACCAP)()->GetBinContent(ib);
		Double_t numOfTrg_ACCAP = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratACCAP = numOfSel_ACCAP / numOfTrg_ACCAP;
		//Double_t crrOfAppAccp = ratACCPR / ratACCAP;
		//IAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		//Double_t errOfAppAccp = std::fabs((ratACCPR / ratACCAP) - crrOfAppAccp);
		//IApperr_accp->pushPoint(Graph::Point_t(cen, errOfAppAccp));
	}
	
	//----  High Energy  ----//
	Graph * HAppnum_pr   = Graph::New("happrlt_ap"  , "", "|Rigidity| [GV]", "N_{p}");
	Graph * HAppnum_ap   = Graph::New("happrlt_pr"  , "", "|Rigidity| [GV]", "N_{#bar{p}}");
	Graph * HApprlt_stat = Graph::New("happrlt_stat", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * HApprlt_full = Graph::New("happrlt_full", "", "|Rigidity| [GV]", "#bar{p}/p Ratio");
	Graph * HAppcrr_accp = Graph::New("happcrr_accp", "", "|Rigidity| [GV]", "Acceptance Correction");
	Graph * HApperr_totl = Graph::New("happerr_totl", "", "|Rigidity| [GV]", "Total Error [%]");
	Graph * HApperr_stat = Graph::New("happerr_stat", "", "|Rigidity| [GV]", "Statistics Error [%]");
	Graph * HApperr_sysm = Graph::New("happerr_sysm", "", "|Rigidity| [GV]", "Systematic Error [%]");
	Graph * HApperr_temp = Graph::New("happerr_temp", "", "|Rigidity| [GV]", "Template Error [%]");
	Graph * HApperr_accp = Graph::New("happerr_accp", "", "|Rigidity| [GV]", "Acceptance Error [%]");
	Graph * HApperr_xsec = Graph::New("happerr_xsec", "", "|Rigidity| [GV]", "Cross-Section Error [%]");
/*	
	Hist * hH2D_pos = Hist::New((TH1*)fIss->Get(CStrFmt("hHp_CCest%s", vtme.c_str())));
	Hist * hH2D_neg = Hist::New((TH1*)fIss->Get(CStrFmt("hHn_CCest%s", vtme.c_str())));
	Hist * hH2D_sig = Hist::New("hHs_CCest", (TH1*)fIss->Get(CStrFmt("hHp_CCest%s", vtme.c_str())));
	Hist * hH2D_bkg = Hist::New("hHb_CCest", (TH1*)fPrL1a9->Get("hHn_CCest"));
		
	std::vector<Hist *>&& hHVec_pos = Hist::Project(Hist::kProjY, hH2D_pos);
	std::vector<Hist *>&& hHVec_neg = Hist::Project(Hist::kProjY, hH2D_neg);
	std::vector<Hist *>&& hHVec_sig = Hist::Project(Hist::kProjY, hH2D_sig);
	std::vector<Hist *>&& hHVec_bkg = Hist::Project(Hist::kProjY, hH2D_bkg);
	
	Hist * hH_ACCPR = Hist::New("hHp_ACCPR_MCEvt", (TH1*)fPr1800->Get("hHp_MCEvt"));
	Hist * hH_ACCAP = Hist::New("hHn_ACCAP_MCEvt", (TH1*)fAp1800->Get("hHn_MCEvt"));

	for (Int_t ib = hBin[0]; ib <= hBin[1]; ++ib) {
		Double_t cen = AXnr.center(ib, Axis::kLog);
		Hist * hH_pos = hHVec_pos.at(ib);
		Hist * hH_neg = hHVec_neg.at(ib);
		Hist * hH_sig = hHVec_sig.at(ib);
		Hist * hH_bkg = hHVec_bkg.at(ib);

		Double_t numOfPr     = (*hH_pos)()->Integral();
		Double_t errOfPrStat = std::sqrt(numOfPr);
		Double_t errOfPr     = errOfPrStat;

		Double_t numOfSmp = (*hH_neg)()->Integral();
		if (numOfSmp < 20.) continue;

		Fit::RooParam param("bta", hH_neg, Fit::RooParam::VLIST({ hH_bkg, hH_sig }));
		Fit::RooSysResult rlt(param);
		if (!rlt.valid()) continue;
		if (!rlt.sysValid()) continue;
		Double_t numOfAp     = rlt.val(1);
		Double_t errOfApStat = rlt.err(1);
		Double_t errOfApTemp = rlt.sysErr(1);
		Double_t errOfAp     = std::sqrt(errOfApStat * errOfApStat + errOfApTemp * errOfApTemp);

		Double_t crrOfAppAccp = 0.362262 * std::exp(-0.857079 * std::log(cen)) + 1.06241;
		Double_t valOfApp     = crrOfAppAccp * numOfAp / numOfPr; 
		Double_t errOfAppStat = valOfApp * std::sqrt(
		                        errOfPrStat*errOfPrStat/numOfPr/numOfPr + 
														errOfApStat*errOfApStat/numOfAp/numOfAp
														);
		Double_t errOfAppTemp = valOfApp * std::sqrt(
														errOfApTemp*errOfApTemp/numOfAp/numOfAp
														);
		Double_t errOfAppAccp = valOfApp * 0.015;
		Double_t errOfAppXsec = valOfApp * 0.04;
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
		HApprlt_full->pushPointWithError(Graph::Point_t(cen, valOfApp), Graph::Error_t(0.0, errOfAppTotl));
		HAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		HApperr_totl->pushPoint(Graph::Point_t(cen, 100. * errOfAppTotl/valOfApp));
		HApperr_stat->pushPoint(Graph::Point_t(cen, 100. * errOfAppStat/valOfApp));
		HApperr_sysm->pushPoint(Graph::Point_t(cen, 100. * errOfAppSysm/valOfApp));
		HApperr_temp->pushPoint(Graph::Point_t(cen, 100. * errOfAppTemp/valOfApp));
		HApperr_accp->pushPoint(Graph::Point_t(cen, 100. * errOfAppAccp/valOfApp));
		HApperr_xsec->pushPoint(Graph::Point_t(cen, 100. * errOfAppXsec/valOfApp));

		editor.create(StrFmt("H%02d", ib));
		editor.cd(1, Canvas::AxisScl_t(0, 0));
		rlt.samp() ->setStyle(Area(), Line(kBlack), Marker(kBlack));
		rlt.sumt() ->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
		rlt.temp(0)->setStyle(Area(), Line(kBlue), Marker(kBlue));
		rlt.temp(1)->setStyle(Area(), Line(kRed), Marker(kRed));
		Hist::Collect(StrFmt("fit%03d", ib), "", Hist::VLIST({ rlt.samp(), rlt.sumt(), rlt.temp(0), rlt.temp(1) }))->Draw("nostack hist");
		Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GeV ~ %.2f GeV )", AXnr.bins(ib-1), AXnr.bins(ib)), kBlack, 0.03), Text::Att_t(0.85, 0.92, 32));
		Text::Draw(Text::Txt_t(StrFmt("N_{Background} ( %.2f +- stat %.2f sys %.2f )", rlt.val(0), rlt.err(0), rlt.sysErr(0)), kBlue, 0.03), Text::Att_t(0.15, 0.85, 12));
		Text::Draw(Text::Txt_t(StrFmt("N_{Antiproton} ( %.2f +- stat %.2f sys %.2f )", rlt.val(1), rlt.err(1), rlt.sysErr(1)), kRed, 0.03), Text::Att_t(0.15, 0.81, 12));
		Text::Draw(Text::Txt_t(StrFmt("#bar{p}/p Ratio ( %E +- stat %E )", valOfApp, errOfAppStat), kBlack, 0.03), Text::Att_t(0.15, 0.77, 12));
		Text::Draw(Text::Txt_t(StrFmt("#chi^{2}/NDF ( %.2f )", rlt.chisq()/rlt.ndf()), kBlack, 0.03), Text::Att_t(0.15, 0.73, 12));
		editor.save();
		
		CdfFlux->SetParameters(trgPr1800, 1.0, 800.0);
		Double_t numOfSel_ACCPR = (*hH_ACCPR)()->GetBinContent(ib);
		Double_t numOfTrg_ACCPR = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratACCPR = numOfSel_ACCPR / numOfTrg_ACCPR;
		CdfFlux->SetParameters(trgAp1800, 1.0, 800.0);
		Double_t numOfSel_ACCAP = (*hH_ACCAP)()->GetBinContent(ib);
		Double_t numOfTrg_ACCAP = CdfFlux->Eval(AXnr.bins(ib)) - CdfFlux->Eval(AXnr.bins(ib-1));
		Double_t ratACCAP = numOfSel_ACCAP / numOfTrg_ACCAP;
		//Double_t crrOfAppAccp = ratACCPR / ratACCAP;
		//HAppcrr_accp->pushPoint(Graph::Point_t(cen, crrOfAppAccp));
		//Double_t errOfAppAccp = std::fabs((ratACCPR / ratACCAP) - crrOfAppAccp);
		//HApperr_accp->pushPoint(Graph::Point_t(cen, errOfAppAccp));
	}
*/

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
	editor.cd(1, Canvas::AxisScl_t(0, 0));
	LApprlt_stat->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApprlt_stat->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApprlt_stat->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CApprlt_stat", ";|Rigidity| [GV];#bar{p}/p Ratio", Graph::VLIST({ LApprlt_stat, IApprlt_stat, HApprlt_stat }))->Draw("ap");
	//Apprlt_Pamela->Draw("same p");
	Apprlt_AMSPub->Draw("same p");
	editor.save();
	
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(0, 0));
	LApprlt_full->setStyle(Area(), Line(kGreen+2), Marker(kGreen+2));
	IApprlt_full->setStyle(Area(), Line(kBlue), Marker(kBlue));
	HApprlt_full->setStyle(Area(), Line(kRed), Marker(kRed));
	Graph::Collect("CApprlt_full", ";|Rigidity| [GV];#bar{p}/p Ratio", Graph::VLIST({ LApprlt_full, IApprlt_full, HApprlt_full }))->Draw("ap");
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
	
	editor.close();

	TFile * file = new TFile(CStrFmt("YiAna%s.root", vtme.c_str()), "RECREATE");
	Graph::Write();
	file->Write();
	file->Close();

	return 1;
}

#endif // __YiFinalAna_C__
