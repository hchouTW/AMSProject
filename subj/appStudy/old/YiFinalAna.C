#ifndef __YiFinalAna_C__
#define __YiFinalAna_C__

// User defination library
//#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/STL/stl.h"
//#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/ROOT/selflib.h"
#include "/afs/cern.ch/user/h/hchou/private/AMSProject/libraries/CPPLibs/CPPLibs.h"
#include "/afs/cern.ch/user/h/hchou/private/AMSProject/libraries/ROOTLibs/ROOTLibs.h"

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/src/AppRlt.h"
//#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/src/Unfold.h"

#include "TGaxis.h"

using namespace MGROOT;

int main(int argc, const char ** argv) {
    Style::LoadDefaultEnvironment();

	std::string cvs = (argc<=1) ? "" : STR_FMT("_%s", argv[1]);

	std::string goddir = "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/rlt17";
	
	// AMS Pub Result
	TString amspubFilePath = "/afs/cern.ch/user/h/hchou/public/DATASET/physics/antipp3/others/20160406MIT.root";
	TFile * inFile_Amspub = TFile::Open(amspubFilePath.Data());
	TGraphErrors * Apprlt_AMSPub = (TGraphErrors *) inFile_Amspub->Get("gpbarp");
	Apprlt_AMSPub->SetMarkerColor(kCyan+1);
	Apprlt_AMSPub->SetLineColor(kCyan+1);
	Apprlt_AMSPub->SetMarkerSize(3.0);
	Apprlt_AMSPub->SetLineWidth(2.0);
	Apprlt_AMSPub->GetXaxis()->SetTitle("|Rigidity| [GV]");
	Apprlt_AMSPub->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");

	for (Int_t it = Apprlt_AMSPub->GetN() -1; it >= 0; --it) {
		//if      (Apprlt_AMSPub->GetX()[it] > 30.) Apprlt_AMSPub->RemovePoint(it);
		//else if (Apprlt_AMSPub->GetX()[it] <  1.) Apprlt_AMSPub->RemovePoint(it);
		//else Apprlt_AMSPub->GetX()[it] *= 0.99;
	}

	std::cout << "\n=======  Load ROOT File  =======\n";
	//RootFile::SourceDir = "/afs/cern.ch/work/h/hchou/BSUB_ANALYSIS/antipp5/NEW07rat/NEW07__";
	//RootFile::SourceDir = "/afs/cern.ch/work/h/hchou/BSUB_ANALYSIS/antipp5/NEW06__"; // time
	RootFile::SourceDir = "/afs/cern.ch/work/h/hchou/BSUB_ANALYSIS/antipp5/NEW08__";
	RootFile fIss("ISS");
	RootFile fPrL1a9("L1a9", 20, 16000, &AppRlt::AXnr);
	RootFile fPrL1o9("L1o9", 20, 16000, &AppRlt::AXnr);
	RootFile fPrC090("C090", 1, 200, &AppRlt::AXnr);
	RootFile fPrC110("C110", 1, 200, &AppRlt::AXnr);
	RootFile fPr0510("Pr0_510", 0.5,  10, &AppRlt::AXnr);
	RootFile fPr1800("Pr1800" ,   1, 800, &AppRlt::AXnr);
	RootFile fAp0510("Ap0_510", 0.5,  10, &AppRlt::AXnr);
	RootFile fAp1800("Ap1800",    1, 800, &AppRlt::AXnr);
	
	PdfEditor editor(Window(WindowSize::kMac), STR_FMT("YiAna%s", cvs.c_str()), goddir);
	AppFitRlt::editor = &editor;

	//----  Low Energy  ----//
	std::cout << "\n=======  Low Energy  =======\n";
	Hist * hLXsecL = Hist::New("hLXsecL_MCEvt", fPrC090.GetHist("hLp_MCEvt"));
	Hist * hLXsecU = Hist::New("hLXsecU_MCEvt", fPrC110.GetHist("hLp_MCEvt"));
	Hist * hLAccPR = Hist::New("hLAccPR_MCEvt", fPr1800.GetHist("hLp_MCEvt"));
	Hist * hLAccAP = Hist::New("hLAccAP_MCEvt", fAp1800.GetHist("hLn_MCEvt"));
	
	Hist * hLXsec = AppRlt::CalXsecErr("hLXsec", hLXsecL, fPrC090.fHGen, hLXsecU, fPrC110.fHGen);
	Hist * hLAccp = AppRlt::CalAccpCrr("hLAccp", hLAccPR, fPr1800.fHGen, hLAccAP, fAp1800.fHGen);

	Hist * hL2Dp = Hist::New(fIss.GetHist(STR_FMT("hLp_TofM%s", cvs.c_str())));
	Hist * hL2Dn = Hist::New(fIss.GetHist(STR_FMT("hLn_TofM%s", cvs.c_str())));
	Hist * hL2Ds = Hist::New(fIss.GetHist(STR_FMT("hLs_TofM%s", cvs.c_str())));
	Hist * hL2Db = Hist::New(fIss.GetHist(STR_FMT("hLb_TofM%s", cvs.c_str())));
	
	AppRlt appl(1, 19, "L", hLAccp);
	//AppRlt appl(1, 6, "L", hLAccp);
	
	AppFitRlt::savePdf = true;
	appl.fit("Mass Estimator", hL2Dp, hL2Dn, hL2Ds, hL2Db);
	AppFitRlt::savePdf = false;
	
	Hist * hTL3Dp = Hist::New(fIss.GetHist("hTLp_TofM"));
	Hist * hTL3Dn = Hist::New(fIss.GetHist("hTLn_TofM"));
	Hist * hTL3Ds = Hist::New(fIss.GetHist("hTLs_TofM"));
	Hist * hTL3Db = Hist::New(fIss.GetHist("hTLb_TofM"));
	//AppFitRlt::savePdf = true;
	//appl.tmefit("Mass Estimator", hTL3Dp, hTL3Dn, hL2Ds, hL2Db);
	//AppFitRlt::savePdf = false;
	

	//Hist * hTL3DpCut = Hist::New(fIss.GetHist(STR_FMT("hTL%sp_Cutflow", vL.c_str())));
	//Hist * hTL3DnCut = Hist::New(fIss.GetHist(STR_FMT("hTL%sn_Cutflow", vL.c_str())));
	//Hist * hLEff = AppRlt::CalCutEff("hLEff", hTL3DpCut, hTL3DnCut);

	//----  Intermedia Energy  ----//
	std::cout << "\n=======  Intermedia Energy  =======\n";
	Hist * hIXsecL = Hist::New("hIXsecL_MCEvt", fPrC090.GetHist("hIp_MCEvt"));
	Hist * hIXsecU = Hist::New("hIXsecU_MCEvt", fPrC110.GetHist("hIp_MCEvt"));
	Hist * hIAccPR = Hist::New("hIAccPR_MCEvt", fPr1800.GetHist("hIp_MCEvt"));
	Hist * hIAccAP = Hist::New("hIAccAP_MCEvt", fAp1800.GetHist("hIn_MCEvt"));
	
	Hist * hIXsec = AppRlt::CalXsecErr("hIXsec", hIXsecL, fPrC090.fHGen, hIXsecU, fPrC110.fHGen);
	Hist * hIAccp = AppRlt::CalAccpCrr("hIAccp", hIAccPR, fPr1800.fHGen, hIAccAP, fAp1800.fHGen);

	Hist * hI2Dp = Hist::New(fIss.GetHist(STR_FMT("hIp_TrdEst%s", cvs.c_str())));
	Hist * hI2Dn = Hist::New(fIss.GetHist(STR_FMT("hIn_TrdEst%s", cvs.c_str())));
	Hist * hI2Ds = Hist::New(fIss.GetHist(STR_FMT("hIs_TrdEst%s", cvs.c_str())));
	Hist * hI2Db = Hist::New(fIss.GetHist(STR_FMT("hIb_TrdEst%s", cvs.c_str())));
	
	AppRlt appi(12, 41, "I", hIAccp);
	//AppRlt appi(5, 18, "I", hIAccp);
	
	AppFitRlt::savePdf = true;
	appi.fit("TRD Estimator", hI2Dp, hI2Dn, hI2Ds, hI2Db);
	AppFitRlt::savePdf = false;
	
	Hist * hTI3Dp = Hist::New(fIss.GetHist("hTIp_TrdEst"));
	Hist * hTI3Dn = Hist::New(fIss.GetHist("hTIn_TrdEst"));
	Hist * hTI3Ds = Hist::New(fIss.GetHist("hTIs_TrdEst"));
	Hist * hTI3Db = Hist::New(fIss.GetHist("hTIb_TrdEst"));
	//AppFitRlt::savePdf = true;
	//appi.tmefit("TRD Estimator", hTI3Dp, hTI3Dn, hI2Ds, hI2Db);
	//AppFitRlt::savePdf = false;

	//Hist * hTI3DpCut = Hist::New(fIss.GetHist(STR_FMT("hTI%sp_Cutflow", vI.c_str())));
	//Hist * hTI3DnCut = Hist::New(fIss.GetHist(STR_FMT("hTI%sn_Cutflow", vI.c_str())));
	//Hist * hIEff = AppRlt::CalCutEff("hIEff", hTI3DpCut, hTI3DnCut);
	
	//----  Intermedia Energy  ----//
	std::cout << "\n=======  Intermedia Energy  =======\n";
	Hist * hMXsecL = Hist::New("hMXsecL_MCEvt", fPrC090.GetHist("hMp_MCEvt"));
	Hist * hMXsecU = Hist::New("hMXsecU_MCEvt", fPrC110.GetHist("hMp_MCEvt"));
	Hist * hMAccPR = Hist::New("hMAccPR_MCEvt", fPr1800.GetHist("hMp_MCEvt"));
	Hist * hMAccAP = Hist::New("hMAccAP_MCEvt", fAp1800.GetHist("hMn_MCEvt"));
	
	Hist * hMXsec = AppRlt::CalXsecErr("hMXsec", hMXsecL, fPrC090.fHGen, hMXsecU, fPrC110.fHGen);
	Hist * hMAccp = AppRlt::CalAccpCrr("hMAccp", hMAccPR, fPr1800.fHGen, hMAccAP, fAp1800.fHGen);

	Hist * hM2Dp = Hist::New(fIss.GetHist(STR_FMT("hMp_TrdEst%s", cvs.c_str())));
	Hist * hM2Dn = Hist::New(fIss.GetHist(STR_FMT("hMn_TrdEst%s", cvs.c_str())));
	Hist * hM2Ds = Hist::New(fIss.GetHist(STR_FMT("hMs_TrdEst%s", cvs.c_str())));
	Hist * hM2Db = Hist::New(fIss.GetHist(STR_FMT("hMb_TrdEst%s", cvs.c_str())));
	
	AppRlt appm(18, 41, "M", hMAccp);
	//AppRlt appm(6, 18, "M", hMAccp);
	
	AppFitRlt::savePdf = true;
	appm.fit("TRD Estimator", hM2Dp, hM2Dn, hM2Ds, hM2Db);
	AppFitRlt::savePdf = false;

	Hist * hTM3Dp = Hist::New(fIss.GetHist("hTMp_TrdEst"));
	Hist * hTM3Dn = Hist::New(fIss.GetHist("hTMn_TrdEst"));
	Hist * hTM3Ds = Hist::New(fIss.GetHist("hTMs_TrdEst"));
	Hist * hTM3Db = Hist::New(fIss.GetHist("hTMb_TrdEst"));
	//AppFitRlt::savePdf = true;
	//appm.tmefit("TRD Estimator", hTM3Dp, hTM3Dn, hM2Ds, hM2Db);
	//AppFitRlt::savePdf = false;


	//Hist * hTM3DpCut = Hist::New(fIss.GetHist(STR_FMT("hTM%sp_Cutflow", vM.c_str())));
	//Hist * hTM3DnCut = Hist::New(fIss.GetHist(STR_FMT("hTM%sn_Cutflow", vM.c_str())));
	//Hist * hMEff = AppRlt::CalCutEff("hMEff", hTM3DpCut, hTM3DnCut);











	//----  High Energy  ----//
	
	std::cout << "\n=======  High Energy  =======\n";
	Hist * hHXsecLL1 = Hist::New("hHXsecL_MCEvtL1", fPrC090.GetHist("hHp_MCEvtL1"));
	Hist * hHXsecUL1 = Hist::New("hHXsecU_MCEvtL1", fPrC110.GetHist("hHp_MCEvtL1"));
	Hist * hHAccPRL1 = Hist::New("hHAccPR_MCEvtL1", fPr1800.GetHist("hHp_MCEvtL1"));
	Hist * hHAccAPL1 = Hist::New("hHAccAP_MCEvtL1", fAp1800.GetHist("hHn_MCEvtL1"));
	
	Hist * hHXsecL1 = AppRlt::CalXsecErr("hHL1Xsec", hHXsecLL1, fPrC090.fHGen, hHXsecUL1, fPrC110.fHGen);
	Hist * hHAccpL1 = AppRlt::CalAccpCrr("hHL1Accp", hHAccPRL1, fPr1800.fHGen, hHAccAPL1, fAp1800.fHGen);

	Hist * hH2DpL1 = Hist::New(fIss   .GetHist(STR_FMT("hHp_TrEstL1%s", cvs.c_str())));
	Hist * hH2DnL1 = Hist::New(fIss   .GetHist(STR_FMT("hHn_TrEstL1%s", cvs.c_str())));
	Hist * hH2DsL1 = Hist::New(STR_FMT("hHs_TrEstL1%s", cvs.c_str()), fIss   .GetHist(STR_FMT("hHp_TrEstL1%s", cvs.c_str())));
	Hist * hH2DbL1 = Hist::New(STR_FMT("hHb_TrEstL1%s", cvs.c_str()), fPrL1o9.GetHist(STR_FMT("hHn_TrEstL1%s", cvs.c_str())));
	
	AppRlt apphL1(35, 57, "HL1", hHAccpL1);
	AppFitRlt::savePdf = true;
	apphL1.fit("CC Estimator", hH2DpL1, hH2DnL1, hH2DsL1, hH2DbL1);
	AppFitRlt::savePdf = false;



	Hist * hHXsecLL9 = Hist::New("hHXsecL_MCEvtL9", fPrC090.GetHist("hHp_MCEvtL9"));
	Hist * hHXsecUL9 = Hist::New("hHXsecU_MCEvtL9", fPrC110.GetHist("hHp_MCEvtL9"));
	Hist * hHAccPRL9 = Hist::New("hHAccPR_MCEvtL9", fPr1800.GetHist("hHp_MCEvtL9"));
	Hist * hHAccAPL9 = Hist::New("hHAccAP_MCEvtL9", fAp1800.GetHist("hHn_MCEvtL9"));
	
	Hist * hHXsecL9 = AppRlt::CalXsecErr("hHL9Xsec", hHXsecLL9, fPrC090.fHGen, hHXsecUL9, fPrC110.fHGen);
	Hist * hHAccpL9 = AppRlt::CalAccpCrr("hHL9Accp", hHAccPRL9, fPr1800.fHGen, hHAccAPL9, fAp1800.fHGen);

	Hist * hH2DpL9 = Hist::New(fIss   .GetHist(STR_FMT("hHp_TrEstL9%s", cvs.c_str())));
	Hist * hH2DnL9 = Hist::New(fIss   .GetHist(STR_FMT("hHn_TrEstL9%s", cvs.c_str())));
	Hist * hH2DsL9 = Hist::New(STR_FMT("hHs_TrEstL9%s", cvs.c_str()), fIss   .GetHist(STR_FMT("hHp_TrEstL9%s", cvs.c_str())));
	Hist * hH2DbL9 = Hist::New(STR_FMT("hHb_TrEstL9%s", cvs.c_str()), fPrL1o9.GetHist(STR_FMT("hHn_TrEstL9%s", cvs.c_str())));
	
	AppRlt apphL9(35, 57, "HL9", hHAccpL9);
	AppFitRlt::savePdf = true;
	apphL9.fit("CC Estimator", hH2DpL9, hH2DnL9, hH2DsL9, hH2DbL9);
	AppFitRlt::savePdf = false;

	Hist * hHXsecLFs = Hist::New("hHXsecL_MCEvtFs", fPrC090.GetHist("hHp_MCEvtFs"));
	Hist * hHXsecUFs = Hist::New("hHXsecU_MCEvtFs", fPrC110.GetHist("hHp_MCEvtFs"));
	Hist * hHAccPRFs = Hist::New("hHAccPR_MCEvtFs", fPr1800.GetHist("hHp_MCEvtFs"));
	Hist * hHAccAPFs = Hist::New("hHAccAP_MCEvtFs", fAp1800.GetHist("hHn_MCEvtFs"));
	
	Hist * hHXsecFs = AppRlt::CalXsecErr("hHFsXsec", hHXsecLFs, fPrC090.fHGen, hHXsecUFs, fPrC110.fHGen);
	Hist * hHAccpFs = AppRlt::CalAccpCrr("hHFsAccp", hHAccPRFs, fPr1800.fHGen, hHAccAPFs, fAp1800.fHGen);

	Hist * hH2DpFs = Hist::New(fIss   .GetHist(STR_FMT("hHp_TrEstFs%s", cvs.c_str())));
	Hist * hH2DnFs = Hist::New(fIss   .GetHist(STR_FMT("hHn_TrEstFs%s", cvs.c_str())));
	Hist * hH2DsFs = Hist::New(STR_FMT("hHs_TrEstFs%s", cvs.c_str()), fIss   .GetHist(STR_FMT("hHp_TrEstFs%s", cvs.c_str())));
	Hist * hH2DbFs = Hist::New(STR_FMT("hHb_TrEstFs%s", cvs.c_str()), fPrL1a9.GetHist(STR_FMT("hHn_TrEstFs%s", cvs.c_str())));
	
	AppRlt apphfs(35, 57, "HFs", hHAccpFs);
	AppFitRlt::savePdf = true;
	apphfs.fit("CC Estimator", hH2DpFs, hH2DnFs, hH2DsFs, hH2DbFs);
	AppFitRlt::savePdf = false;


	
	(*apphfs.hErrStat)()->SetBinContent(57, (*apphfs.hErrStat)()->GetBinContent(57) / sqrt(2.));
	(*apphfs.hErrCntn)()->SetBinContent(57, (*apphfs.hErrCntn)()->GetBinContent(57) / sqrt(2.));
	(*apphfs.hErrAccp)()->SetBinContent(57, (*apphfs.hErrAccp)()->GetBinContent(57) / 2.);
	(*apphfs.hErrRscl)()->SetBinContent(57, (*apphfs.hErrRscl)()->GetBinContent(57) / 2.);
	(*apphfs.hErrSysm)()->SetBinContent(57, 
			std::sqrt((*apphfs.hErrCntn)()->GetBinContent(57)*(*apphfs.hErrCntn)()->GetBinContent(57) + 
				      (*apphfs.hErrAccp)()->GetBinContent(57)*(*apphfs.hErrAccp)()->GetBinContent(57) + 
					  (*apphfs.hErrRscl)()->GetBinContent(57)*(*apphfs.hErrRscl)()->GetBinContent(57)));
	(*apphfs.hErrTotl)()->SetBinContent(57, 
			std::sqrt((*apphfs.hErrStat)()->GetBinContent(57)*(*apphfs.hErrStat)()->GetBinContent(57) + 
					  (*apphfs.hErrSysm)()->GetBinContent(57)*(*apphfs.hErrSysm)()->GetBinContent(57)));
	
	(*apphfs.hRatStat)()->SetBinContent(57, (*apphfs.hRatStat)()->GetBinContent(57) / 2.);
	(*apphfs.hRatStat)()->SetBinError  (57, (*apphfs.hErrStat)()->GetBinContent(57));
	(*apphfs.hRatTotl)()->SetBinContent(57, (*apphfs.hRatTotl)()->GetBinContent(57) / 2.);
	(*apphfs.hRatTotl)()->SetBinError  (57, (*apphfs.hErrTotl)()->GetBinContent(57));
	













	AppRlt appc(-1, -1, "C");
	appc.addRlt(appl,  1, 16);
	appc.addRlt(appi, 17, 23);
	appc.addRlt(appm, 24, 36);
	appc.addRlt(apphL1, apphL9, apphfs, 37, 54);
	appc.addRlt(apphfs, 55, 57);



	//AppRlt appc(-1, -1, "C");
	//appc.addRlt(appl, 1, 5);
	//appc.addRlt(appi, 6, 9);
	//appc.addRlt(appm, 10, 18);





/*
	// Sp Index
	const Float_t bdConst[2] = { 1.0, 2.2 };
	const Float_t bdSlope[2] = { -1.0, 1.6 };
	const Float_t transSC = ((bdConst[1]-bdConst[0]) / (bdSlope[1]-bdSlope[0]));
	TGraphErrors * spConst = new TGraphErrors();
	TGraphErrors * spSlope = new TGraphErrors();
	spConst->SetNameTitle("SpConst", "");
	spSlope->SetNameTitle("SpSlope", "");
	TF1 * fspidx = new TF1("fspidx", "[0]*TMath::Power(x, [1])");
	for (Int_t irig = 3; irig <= AppRlt::AXnr.nbin()-4; irig+=5) {
		TGraphErrors * slice = new TGraphErrors();
		Double_t wx = 0;
		Double_t ws = 0;
		//for (Int_t it = irig; it <= AppRlt::AXnr.nbin(); ++it) {
		for (Int_t it = irig; it <= irig+4; ++it) {
			Double_t err = std::sqrt((*appc.hErrStat)()->GetBinContent(it)*(*appc.hErrStat)()->GetBinContent(it) + 
					                 (*appc.hErrCntn)()->GetBinContent(it)*(*appc.hErrCntn)()->GetBinContent(it));
			Double_t wgt = 1. / err / err;
			wx += wgt * (*appc.hRatStat)()->GetBinCenter(it);
			ws += wgt;
			slice->SetPoint(it,
					(*appc.hRatStat)()->GetBinCenter(it),
					(*appc.hRatStat)()->GetBinContent(it));
			slice->SetPointError(it, 0, err);
		}
		wx /= ws;
		fspidx->SetParameters(1e-6, 1);
		//fspidx->SetParameters(2e-4, 1e-6, wx);
		//fspidx->FixParameter(2, wx);
		slice->Fit(fspidx, "q0S");
		spConst->SetPoint     ((irig-3)/5, wx, fspidx->GetParameter(1));
		spConst->SetPointError((irig-3)/5, 0,  fspidx->GetParError (1));
		//spConst->SetPoint     ((irig-22), wx, fspidx->GetParameter(0) * 1e4);
		//spConst->SetPointError((irig-22), 0,  fspidx->GetParError (0) * 1e4);
		//spSlope->SetPoint     ((irig-22), wx, (fspidx->GetParameter(1) * 1e6 - bdSlope[0]) * transSC + bdConst[0]);
		//spSlope->SetPointError((irig-22), 0,  (fspidx->GetParError (1) * 1e6 * transSC));
		delete slice;
	}
	delete fspidx;
	
	editor.create();
	editor.cd(1, PadAxis(1, 0));
	spConst->Draw("ap");
	spConst->GetHistogram()->SetLineColor(0);
	spConst->SetMarkerSize(3.0);
	spConst->SetLineWidth(2.0);
	spConst->SetLineColor(kRed);
	spConst->SetMarkerColor(kRed);
	spConst->GetXaxis()->SetTitle("|Rigidity| [GV]");
	spConst->GetYaxis()->SetTitle("Spectrum Index #gamma");
	spConst->GetXaxis()->SetMoreLogLabels();
	spConst->GetXaxis()->SetLimits(1, 450);
	//spConst->SetMinimum(bdConst[0]);
	//spConst->SetMaximum(bdConst[1]);
	spConst->Draw("ap");
	editor()().SetTicks(0, 0);
	editor()().Modified();
	editor()().Update();
	//spSlope->Draw("p");

	//TGaxis gaxis(spConst->GetXaxis()->GetXmax(), bdConst[0], spConst->GetXaxis()->GetXmax(), bdConst[1], bdSlope[0], bdSlope[1], 510, "+");
	//gaxis.SetTitle("k [10^{-6}]");
	//gaxis.Draw();
	editor.save();
	


*/





	editor.create();
	editor.cd(1, PadAxis(1, 0));
	(*appc.hErrAccp)()->Divide((*appc.hRatStat)()); (*appc.hErrAccp)()->Scale(100.0);
	(*appc.hErrCntn)()->Divide((*appc.hRatStat)()); (*appc.hErrCntn)()->Scale(100.0);
	(*appc.hErrRscl)()->Divide((*appc.hRatStat)()); (*appc.hErrRscl)()->Scale(100.0);
	(*appc.hErrSysm)()->Divide((*appc.hRatStat)()); (*appc.hErrSysm)()->Scale(100.0);
	(*appc.hErrStat)()->Divide((*appc.hRatStat)()); (*appc.hErrStat)()->Scale(100.0);
	(*appc.hErrTotl)()->Divide((*appc.hRatStat)()); (*appc.hErrTotl)()->Scale(100.0);
	appc.hErrRscl->style(Fill(), Line(kYellow+1), Marker(kYellow+1));
	appc.hErrAccp->style(Fill(), Line(kAzure-2) , Marker(kAzure-2) );
	appc.hErrCntn->style(Fill(), Line(kMagenta) , Marker(kMagenta) );
	appc.hErrSysm->style(Fill(), Line(kGreen+1) , Marker(kGreen+1) );
	THStack * hErrSys = Hist::Collect("CApperr_sysm", ";|Rigidity| [GV];Relative Error [%]", HistList({ appc.hErrRscl, appc.hErrAccp, appc.hErrCntn, appc.hErrSysm }));
	hErrSys->Draw("nostack hist");
	hErrSys->GetXaxis()->SetMoreLogLabels();
	hErrSys->Draw("nostack hist");
	TextDraw("Total Systematic Errors", TextStyle(kGreen+1  , 0.04), TextAlign(0.60, 0.82, 32));
	TextDraw(    "Antiproton Counting", TextStyle(kMagenta  , 0.04), TextAlign(0.60, 0.77, 32));
	TextDraw(  "Acceptance Correction", TextStyle(kAzure-2  , 0.04), TextAlign(0.60, 0.72, 32));
	TextDraw("Absolute Rigidity Scale", TextStyle(kYellow+1 , 0.04), TextAlign(0.60, 0.67, 32));

	editor.save();
	
	editor.create();
	editor.cd(1, PadAxis(1, 0));
	appc.hErrSysm->style(Fill(), Line(kGreen+1), Marker(kGreen+1));
	appc.hErrStat->style(Fill(), Line(kBlue)   , Marker(kBlue));
	appc.hErrTotl->style(Fill(), Line(kRed)    , Marker(kRed));
	THStack * hErrTot = Hist::Collect("CApperr_totl", ";|Rigidity| [GV];Relative Error [%]", HistList({ appc.hErrSysm, appc.hErrStat, appc.hErrTotl }));
	hErrTot->Draw("nostack hist");
	hErrTot->GetXaxis()->SetMoreLogLabels();
	hErrTot->Draw("nostack hist");
	TextDraw(      "Total Errors", TextStyle(kRed    , 0.04), TextAlign(0.60, 0.82, 32));
	TextDraw("Statistical Errors", TextStyle(kBlue   , 0.04), TextAlign(0.60, 0.77, 32));
	TextDraw( "Systematic Errors", TextStyle(kGreen+1, 0.04), TextAlign(0.60, 0.72, 32));
	editor.save();

	editor.create();
	editor.cd(1, PadAxis(1, 0));
	appc.hRatStat->style(Fill(), Line(kRed), Marker(kRed));
	appc.hRatStat->draw("pe");
	(*appc.hRatStat)()->GetXaxis()->SetMoreLogLabels();
	(*appc.hRatStat)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
	(*appc.hRatStat)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
	appc.hRatStat->draw("pe");
	editor.save();
	
	editor.create();
	editor.cd(1, PadAxis(1, 0));
	appc.hRatTotl->style(Fill(), Line(kRed), Marker(kRed));
	appc.hRatTotl->draw("pe");
	(*appc.hRatTotl)()->GetXaxis()->SetMoreLogLabels();
	(*appc.hRatTotl)()->GetXaxis()->SetTitle("|Rigidity| [GV]");
	(*appc.hRatTotl)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
	appc.hRatTotl->draw("pe");
	editor.save();
	
	editor.create();
	editor.cd(1, PadAxis(1, 0));
	Apprlt_AMSPub->SetMarkerColor(kBlue);
	Apprlt_AMSPub->SetLineColor(kBlue);
	Apprlt_AMSPub->Draw("ap");
	Apprlt_AMSPub->GetXaxis()->SetMoreLogLabels();
	//Apprlt_AMSPub->GetXaxis()->SetTitle("|Rigidity| [GV]");
	//Apprlt_AMSPub->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
	Apprlt_AMSPub->Draw("ap");
	appc.hRatTotl->style(Fill(), Line(kRed), Marker(kRed));
	appc.hRatTotl->draw("same");
	TextDraw("This Analysis", TextStyle(kRed  , 0.05), TextAlign(0.40, 0.80, 32));
	TextDraw("AMS Published", TextStyle(kBlue , 0.05), TextAlign(0.40, 0.72, 32));
	editor.save();

	editor.create();
	editor.cd(1, PadAxis(1, 0));
	appl.hRatStat->style(Fill(), Line(kGreen+2), Marker(kGreen+2));
	appi.hRatStat->style(Fill(), Line(kBlue), Marker(kBlue));
	appm.hRatStat->style(Fill(), Line(kRed), Marker(kRed));
	apphL1.hRatStat->style(Fill(), Line(kYellow+1), Marker(kYellow+1));
	apphL9.hRatStat->style(Fill(), Line(kYellow+2), Marker(kYellow+2));
	apphfs.hRatStat->style(Fill(), Line(kYellow+3), Marker(kYellow+3));
	//Apprlt_AMSPub->GetXaxis()->SetRangeUser(1, 40);
	Apprlt_AMSPub->GetXaxis()->SetMoreLogLabels();
	Apprlt_AMSPub->Draw("ap");
	Hist::Collect("CApprlt_stat", ";|Rigidity| [GV];#bar{p}/p Flux Ratio", HistList({ appl.hRatStat, appi.hRatStat, appm.hRatStat, apphL1.hRatStat, apphL9.hRatStat, apphfs.hRatStat }))->Draw("nostack pe same");
	//Hist::Collect("CApprlt_stat", ";|Rigidity| [GV];#bar{p}/p Ratio", HistList({ appl.hRatStat, appi.hRatStat, appm.hRatStat }))->Draw("same nostack pe");
	editor.save();

	

	//Unfold(appc.hRatStat, appc.hRatOrg, &fPr1800, &fAp1800);
	//Unfold(1, appc.hRatStat, appl.hNumPr, &fPr1800, &fAp1800);












/*
	std::cerr << "\n==== COMPARISION ====\n";

	const Int_t TTCS_OFF = 16;
	std::vector<Hist*>&& hAppVec = Hist::ProjectAll(HistProj::kX, appc.hTRatStat);
	for (Int_t irig = 1; irig <= AppRlt::AXnr.nbin(); ++irig) {
		Float_t avgRat = (*appc.hRatStat)()->GetBinContent(irig);
		Float_t avgErr = (*appc.hRatStat)()->GetBinError(irig);

		Hist * happ = hAppVec.at(irig);
		(*happ)()->SetBinContent(TTCS_OFF, 0.0);
		(*happ)()->SetBinError(TTCS_OFF, 0.0);

		SetTimeAxis((*happ)()->GetXaxis());
		happ->style(Fill(), Line(kRed), Marker(kRed));

		editor.create();
		editor.cd(1, PadAxis(0, 0));
		(*happ)()->GetXaxis()->SetTitle("Date");
		(*happ)()->GetYaxis()->SetTitle("#bar{p}/p Flux Ratio");
		happ->draw("pe");
		(*happ)()->GetXaxis()->SetMoreLogLabels();
		(*happ)()->SetMaximum(avgRat * 1.5);
		(*happ)()->SetMinimum(avgRat * 0.6);
		happ->draw("pe");
		TextDraw(STR_FMT("|Rigidity| ( %.2f GV ~ %.2f GV )", AppRlt::AXnr(irig-1), AppRlt::AXnr(irig)), TextStyle(kBlack, 0.03), TextAlign(0.85, 0.92, 32));
		editor.save();
	}
	
	
	std::string eppath = "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/data/ElectronPositron_AntiPBin_3Months_PGKIT.root";
	TFile * fElPs = TFile::Open(eppath.c_str());
	Hist * hElPs = Hist::New((TH1*)fElPs->Get("hh_Ratio_StatErr"));
	std::vector<Hist*>&& hElPsVec = Hist::ProjectAll(HistProj::kX, hElPs);
	
	std::vector<Hist *> hAppVec2;
	std::vector<Hist *> hElpsVec;
	std::cout << "TEST O1\n";	
	for (Int_t irig = 1; irig <= AppRlt::AXnr.nbin(); ++irig) {
		if (!(irig >= 1 && irig <= 9)) continue;
		Float_t avgRat = (*appc.hRatStat)()->GetBinContent(irig);
		Float_t avgErr = (*appc.hRatStat)()->GetBinError(irig);
		Hist * happ = hAppVec.at(irig);
		(*happ)()->Scale(1./avgRat);
		hAppVec2.push_back(happ);
		//happ->style(Fill(), Line(Color(irig-1, 10)), Marker(Color(irig-1, 10)));
		
		Float_t avgRat2v = 0;
		Float_t avgRat2s = 0;
		Hist * helps = hElPsVec.at(irig+3);
		SetTimeAxis((*helps)()->GetXaxis());
		helps->style(Fill(), Line(kGreen+2), Marker(kGreen+2));
		for (Int_t itme = 1; itme <= AppRlt::AXtme.nbin(); ++itme) {
			Double_t val = (*helps)()->GetBinContent(itme);
			Double_t err = (*helps)()->GetBinError  (itme);
			Double_t wgt = (1. / err /err);
			if (MGNumc::EqualToZero(val)) continue;
			avgRat2v += wgt * val;
			avgRat2s += wgt;
		}
		avgRat2v /= avgRat2s;
		(*helps)()->Scale(1./avgRat2v);

		hElpsVec.push_back(helps);
	}
		
	std::cout << "TEST O2\n";	
	for (Int_t ibin = 0; ibin < hAppVec2.size(); ++ibin) {
		Int_t irig = 1 + ibin;
		Hist * happ = hAppVec2.at(ibin);
		Hist * heps = hElpsVec.at(ibin);

		editor.create();
		editor.cd(1, PadAxis(0, 0));
		(*heps)()->GetXaxis()->SetTitle("Date");
		(*heps)()->GetYaxis()->SetTitle("Normalized Flux Ratio by Average Flux Ratio");
		heps->draw("pe");
		(*heps)()->GetXaxis()->SetMoreLogLabels();
		(*heps)()->SetMaximum(1.5);
		(*heps)()->SetMinimum(0.6);
		heps->draw("pe");
		happ->draw("pe same");
		TextDraw(STR_FMT("|Rigidity| ( %.2f GV ~ %.2f GV )", AppRlt::AXnr(irig-1), AppRlt::AXnr(irig)), TextStyle(kBlack, 0.03), TextAlign(0.85, 0.92, 32));
		TextDraw("#bar{p}/p Flux Ratio", TextStyle(kRed  , 0.05), TextAlign(0.80, 0.80, 32));
		TextDraw("e^{-}/e^{+} Flux Ratio", TextStyle(kGreen+2 , 0.05), TextAlign(0.80, 0.72, 32));
		editor.save();
	}
	

	// testcode
	///////////////////////
	Hist * hEl = Hist::New((TH1*)fElPs->Get("hh_EleFlux_StatErr"));
	Hist * hPs = Hist::New((TH1*)fElPs->Get("hh_PosFlux_StatErr"));
	
	std::string prpath = "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/data/ProtonFlux_AntiPBin_3Months_MiB.root";
	TFile * fPr = TFile::Open(prpath.c_str());
	Hist * hPr = Hist::New((TH1*)fPr->Get("hh_ProtonFlux_TotErr"));

	Hist * hAp = Hist::New("hTC_ApFlxStat", "", HistAxis(AppRlt::AXtme, AppRlt::AXnr));

	std::cout << "TEST1\n";	
	for (Int_t itme = 1; itme <= AppRlt::AXtme.nbin(); ++itme) {
		for (Int_t irig = 1; irig <= AppRlt::AXnr.nbin(); ++irig) {
			double appv = (*appc.hTRatStat)()->GetBinContent(itme, irig);
			double apps = (*appc.hTRatStat)()->GetBinError  (itme, irig);
			double prv = (*hPr)()->GetBinContent(itme, irig+3);
			double prs = (*hPr)()->GetBinError  (itme, irig+3);
			double apv = prv * appv;
			double aps = prv * apps;
			(*hAp)()->SetBinContent(itme, irig, apv);
			(*hAp)()->SetBinError  (itme, irig, aps);
		}
	}



	std::vector<Hist *>&& hflxPsVec = Hist::ProjectAll(HistProj::kX, hPs);
	std::vector<Hist *>&& hflxElVec = Hist::ProjectAll(HistProj::kX, hEl);
	std::vector<Hist *>&& hflxPrVec = Hist::ProjectAll(HistProj::kX, hPr);
	std::vector<Hist *>&& hflxApVec = Hist::ProjectAll(HistProj::kX, hAp);

	std::vector< std::vector<Hist *> > hflx;

	std::cout << "TEST2\n";	
	for (Int_t irig = 1; irig <= AppRlt::AXnr.nbin(); ++irig) {
		if (!(irig >= 1 && irig <= 9)) continue;
		Hist * hparts[4] = { hflxPsVec.at(irig+3), hflxElVec.at(irig+3), hflxPrVec.at(irig+3), hflxApVec.at(irig) };
		Color_t cols[4] = { kGreen+2, kYellow+2, kBlue, kRed };
		
		std::vector<Hist *> elm;
		for (Int_t it = 0; it < 4; ++it) {
			Hist * hpart = hparts[it];
			SetTimeAxis((*hpart)()->GetXaxis());
			Float_t avgFlxv = 0;
			Float_t avgFlxs = 0;
			hpart->style(Fill(), Line(cols[it]), Marker(cols[it]));
			for (Int_t itme = 1; itme <= AppRlt::AXtme.nbin(); ++itme) {
				if (itme==TTCS_OFF) {
					(*hpart)()->SetBinContent(TTCS_OFF, 0.0);
					(*hpart)()->SetBinError  (TTCS_OFF, 0.0);
					continue;
				}
				Double_t val = (*hpart)()->GetBinContent(itme);
				Double_t err = (*hpart)()->GetBinError  (itme);
				Double_t wgt = (1. / err /err);
				if (MGNumc::EqualToZero(val)) continue;
				avgFlxv += wgt * val;
				avgFlxs += wgt;
			}
			avgFlxv /= avgFlxs;
			(*hpart)()->Scale(1./avgFlxv);
			elm.push_back(hpart);
		}


		hflx.push_back(elm);
	}
		
	std::cout << "TEST2\n";	
	for (Int_t ibin = 0; ibin < hflx.size(); ++ibin) {
		Int_t irig = 1 + ibin;
		std::cout << "TEST FIN 00\n";	
		std::vector<Hist *> hparts = hflx.at(ibin);

		editor.create();
		editor.cd(1, PadAxis(0, 0));
		std::cout << "TEST FIN 01\n";	
		std::cout << "TEST FIN 02\n";	
		hparts.at(0)->draw("pe");
		(*hparts.at(0))()->GetXaxis()->SetTitle("Date");
		(*hparts.at(0))()->GetYaxis()->SetTitle("Normalized Flux by Average Flux");
		(*hparts.at(0))()->GetXaxis()->SetMoreLogLabels();
		(*hparts.at(0))()->SetMaximum(1.5);
		(*hparts.at(0))()->SetMinimum(0.6);
		hparts.at(0)->draw("pe");
		hparts.at(1)->draw("pe same");
		hparts.at(2)->draw("pe same");
		hparts.at(3)->draw("pe same");
		std::cout << "TEST FIN 03\n";	
		TextDraw(STR_FMT("|Rigidity| ( %.2f GV ~ %.2f GV )", AppRlt::AXnr(irig-1), AppRlt::AXnr(irig)), TextStyle(kBlack, 0.03), TextAlign(0.85, 0.92, 32));
		TextDraw( "#bar{p} Flux", TextStyle(kRed      , 0.05), TextAlign(0.55, 0.80, 32));
		TextDraw(       "p Flux", TextStyle(kBlue     , 0.05), TextAlign(0.55, 0.74, 32));
		TextDraw(   "e^{-} Flux", TextStyle(kYellow+2 , 0.05), TextAlign(0.70, 0.80, 32));
		TextDraw(   "e^{+} Flux", TextStyle(kGreen+2  , 0.05), TextAlign(0.70, 0.74, 32));
		editor.save();
	}
	std::cout << "TEST FIN 04\n";	

	////////////////

*/


















	//editor.create();
	//editor.cd(1, PadAxis(0, 0));
	//THStack * hTapp = Hist::Collect("TApprlt", ";Date;#bar{p}/p Normalized Flux Ratio", hAppVec2);
	//hTapp->Draw("nostack pe");
	//SetTimeAxis(hTapp->GetXaxis());
	//hTapp->SetMaximum(1.4);
	//hTapp->SetMinimum(0.6);
	//hTapp->Draw("nostack pe");
	//editor.save();



	editor.close();

	TFile * file = new TFile(CSTR_FMT("%s/YiAna%s.root", goddir.c_str(), cvs.c_str()), "RECREATE");
	Hist::Write();
    Apprlt_AMSPub->Write();
	file->Write();
	file->Close();

	return 1;
}

#endif // __YiFinalAna_C__
