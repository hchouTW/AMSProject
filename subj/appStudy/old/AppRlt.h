#ifndef __AppRlt_H__
#define __AppRlt_H__

// User defination library
//#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/STL/stl.h"
//#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/ROOT/selflib.h"
#include "/afs/cern.ch/user/h/hchou/private/AMSProject/libraries/CPPLibs/CPPLibs.h"
#include "/afs/cern.ch/user/h/hchou/private/AMSProject/libraries/ROOTLibs/ROOTLibs.h"

using namespace MGROOT;

void SetTimeAxis(TAxis * axis) {
	axis->SetTitle("Time");
	axis->SetTimeDisplay(1);
	axis->SetTimeFormat("#splitline{%b}{%Y}%F1970-01-01 00:00:00s0");
	axis->SetNdivisions(1010);
	axis->SetLabelOffset(0.02);
	axis->SetTitleOffset(1.5);
	axis->SetTitleFont(42);
}

class RootFile {
	public :
		RootFile(const std::string& patt, Double_t moml = 0, Double_t momu = 0, Axis * axis = nullptr);
		~RootFile() {}

		Double_t GetNTrg  (Double_t bdl, Double_t bdu);
		Double_t GetNTrg27(Double_t bdl, Double_t bdu);
		
		TH1 *    GetHist(const std::string& name);

		void LoadTree();

	protected :
		TFile *  fFile;
		UInt_t   fNTrg;
		Double_t fMoml;
		Double_t fMomu;

	public :
		Hist * fHGen;
		Hist * fHGen27;

	public :
		TTree *  fTree;
		Float_t  fMCRig;
		Float_t  fTRRig;
		Bool_t   fOptL;
		Bool_t   fOptI;
		Bool_t   fOptM;
		Bool_t   fOptHL1;
		Bool_t   fOptHL9;
		Bool_t   fOptHFs;

	public :
		static std::string SourceDir;
		static TF1 *       CdfFlux; 
		static TF1 *       CdfFlux27; 
};

std::string RootFile::SourceDir = "/afs/cern.ch/work/h/hchou/BSUB_ANALYSIS/antipp5/NEW07__";
TF1 *       RootFile::CdfFlux   = new TF1("CdfFlux"  , "[0]/(TMath::Log([2])-TMath::Log([1]))*TMath::Log(x)");
TF1 *       RootFile::CdfFlux27 = new TF1("CdfFlux27", "[0]/(TMath::Log([2])-TMath::Log([1]))*(-TMath::Power(x, -1.7)/1.7)");

RootFile::RootFile(const std::string& patt, Double_t moml, Double_t momu, Axis * axis) : fFile(0), fNTrg(0), fMoml(0), fMomu(0), fHGen(0), fHGen27(0), fTree(0), fMCRig(0), fTRRig(0), fOptL(0), fOptI(0), fOptM(0), fOptHL1(0), fOptHL9(0), fOptHFs(0) {
	fFile = TFile::Open(CSTR_FMT("%s%s/YiAnalytics.root", SourceDir.c_str(), patt.c_str()));
	if (fFile == nullptr) MGSys::ShowErrorAndExit(STR_FMT("%s Files is not exist.", patt.c_str()));

	TTree * tree = (TTree *) fFile->Get("DSTRun");
	UInt_t  ntrg = 0, trg = 0;
	tree->SetBranchAddress("trgEV", &trg);
	for (Long64_t it = 0; it < tree->GetEntries(); ++it) {
		tree->GetEntry(it); ntrg += trg;
	}
	fNTrg = ntrg;
	fMoml = moml;
	fMomu = momu;

	if (axis && MGNumc::Compare(moml, momu) != 0) {
		UInt_t nTrgSum   = 0;
		UInt_t nTrgSum27 = 0;
		fHGen   = Hist::New(STR_FMT("hGen_%s"  , patt.c_str()), "", HistAxis((*axis)));
		fHGen27 = Hist::New(STR_FMT("hGen27_%s", patt.c_str()), "", HistAxis((*axis)));
		for (Int_t irig = 1; irig <= axis->nbin(); ++irig) {
			if (momu < (*axis)(irig-1) || moml > (*axis)(irig)) continue;
			(*fHGen  )()->SetBinContent(irig, GetNTrg  ((*axis)(irig-1), (*axis)(irig)));
			(*fHGen27)()->SetBinContent(irig, GetNTrg27((*axis)(irig-1), (*axis)(irig)));
			nTrgSum   += (*fHGen  )()->GetBinContent(irig);
			nTrgSum27 += (*fHGen27)()->GetBinContent(irig);
		}
		std::cout << CSTR_FMT("Info (%12s) : NTrg %14ld %14ld\n", patt.c_str(), fNTrg, nTrgSum);
	}
}

Double_t RootFile::GetNTrg(Double_t bdl, Double_t bdu) {
	Double_t vl = (bdl > fMoml) ? bdl : fMoml;
	Double_t vu = (bdu < fMomu) ? bdu : fMomu;
	CdfFlux->SetParameters(fNTrg, fMoml, fMomu);
	Double_t num = (CdfFlux->Eval(vu) - CdfFlux->Eval(vl));
	return num;
}

Double_t RootFile::GetNTrg27(Double_t bdl, Double_t bdu) {
	Double_t vl = (bdl > fMoml) ? bdl : fMoml;
	Double_t vu = (bdu < fMomu) ? bdu : fMomu;
	CdfFlux27->SetParameters(fNTrg, fMoml, fMomu);
	Double_t num = (CdfFlux27->Eval(vu) - CdfFlux27->Eval(vl));
	return num;
}
		
TH1 * RootFile::GetHist(const std::string& name) {
	if (fFile == nullptr) return nullptr;
	return (TH1*)fFile->Get(name.c_str());
}

void RootFile::LoadTree() {
	fTree = (TTree*)fFile->Get("data");
	if (fTree) {
		fTree->SetBranchAddress("MCRig" , &fMCRig);
		fTree->SetBranchAddress("TRRig" , &fTRRig);
		fTree->SetBranchAddress("OptL"  , &fOptL);
		fTree->SetBranchAddress("OptI"  , &fOptI);
		fTree->SetBranchAddress("OptM"  , &fOptM);
		fTree->SetBranchAddress("OptHL1", &fOptHL1);
		fTree->SetBranchAddress("OptHL9", &fOptHL9);
		fTree->SetBranchAddress("OptHFs", &fOptHFs);
	}
}

class AppFitRlt {
	public :
		AppFitRlt(const std::string& name, Hist * hPos, Hist * hNeg, Hist * hSig, Hist * hBkg, Double_t rigl = 0, Double_t rigu = 0, Double_t refcrr = 1, Int_t time = -1);
		~AppFitRlt() {}

	public :
		Bool_t   isSucc;

		Double_t numPr;
		Double_t numAp;
		Double_t orgAppr;
		Double_t orgAppe;
		Double_t ratAppc;
		Double_t crrAccp;
		Double_t errStat;
		Double_t errCntn;
		Double_t errAccp;
		Double_t errRscl;
		Double_t errSysm;
		Double_t errTotl;

		Double_t fitSN;
		Double_t fitNChi;
		Double_t fitChi;
		Int_t    fitNdf;

	public :
		static Int_t       ntimes;
		static PdfEditor * editor;
		static Bool_t      savePdf;
};

PdfEditor * AppFitRlt::editor = nullptr;
Int_t       AppFitRlt::ntimes = 400;
Bool_t      AppFitRlt::savePdf = false;

AppFitRlt::AppFitRlt(const std::string& name, Hist * hPos, Hist * hNeg, Hist * hSig, Hist * hBkg, Double_t rigl, Double_t rigu, Double_t refcrr, Int_t time) : isSucc(false), numPr(0), numAp(0), orgAppr(0), orgAppe(0), ratAppc(0), crrAccp(1), errStat(0), errCntn(0), errAccp(0), errRscl(0), errSysm(0), errTotl(0), fitSN(0), fitNChi(0), fitChi(0), fitNdf(0) {
	Double_t numOfPr     = (*hPos)()->Integral();
	Double_t errOfPrStat = std::sqrt(numOfPr);
	Double_t numOfSmp    = (*hNeg)()->Integral();
	if (numOfPr < 20 || numOfSmp < 20.) return;

	Fit::RooVar roovar(name, hNeg, HistList({ hSig, hBkg }));
	Fit::RooSysResult rlt(roovar, true, ntimes);
	if (!rlt.exist()) return;
    Fit::RooVar var = rlt.var();
    Fit::RooPar stdpar = rlt.std_par();
    Fit::RooPar syspar = rlt.sys_par();
	Double_t numOfAp     = stdpar.val(0);
	Double_t errOfApStat = stdpar.err(0);
	Double_t errOfApCntn = syspar.sys(0);

	crrAccp = refcrr;
	
	TF1 UnfoldError("UnfoldError", "[0]*TMath::Power([1]*(x+[2]), [3])", 1, 1000);
	UnfoldError.SetParameters(8.036300e+06, 2.94097e-03, 4.54156e+02, -6.74121e+01);

	numPr   = numOfPr;
	numAp   = numOfAp;
	orgAppr = (numAp / numPr);
	orgAppe = std::sqrt(errOfPrStat*errOfPrStat/numOfPr/numOfPr +
			            errOfApStat*errOfApStat/numOfAp/numOfAp);
	ratAppc = crrAccp * orgAppr;
	errStat = ratAppc * orgAppe;
	errCntn = ratAppc * std::sqrt(errOfApCntn*errOfApCntn/numOfAp/numOfAp);
	errAccp = ratAppc * (0.02 + UnfoldError.Eval(std::sqrt(rigl*rigu)));
	errRscl = ratAppc * (std::sqrt(2.) * 5.16454e-05 * std::pow(std::sqrt(rigl*rigu), 9.23852e-01));
	errSysm = std::sqrt(errCntn*errCntn + errAccp*errAccp + errRscl*errRscl);
	errTotl = std::sqrt(errStat*errStat + errSysm*errSysm);
	fitSN   = stdpar.val(0) / (stdpar.val(0)+stdpar.val(1)); 
	fitNChi = (stdpar.chi() / stdpar.ndf());
	fitChi  = stdpar.chi();
	fitNdf  = stdpar.ndf();
	isSucc  = true;

	if (editor && savePdf) {
		var.samp() ->style(Fill(), Line(kBlack), Marker(kBlack));
		var.sumt() ->style(Fill(), Line(kGreen+2), Marker(kGreen+2));
		var.temp(0)->style(Fill(), Line(kRed), Marker(kRed));
		var.temp(1)->style(Fill(), Line(kBlue), Marker(kBlue));

		((TH1D*)(*var.samp())())->GetYaxis()->SetTitle("Events/Bin");
		((TH1D*)(*var.samp())())->GetXaxis()->SetTitle(name.c_str());

		editor->create();
		editor->cd(1, PadAxis(0, 0));

		var.samp()->draw("pe");
		var.sumt()->draw("same hist");
		var.temp(1)->draw("same hist");
		var.temp(0)->draw("same hist");
		
		if (time != -1)
			TextDraw(STR_FMT("Date ( %02d )", time), TextStyle(kBlack, 0.03), TextAlign(0.85, 0.96, 32));
		if (!MGNumc::EqualToZero(rigl) && !MGNumc::EqualToZero(rigu))
			TextDraw(STR_FMT("Rigidity ( %.2f GV ~ %.2f GV )", rigl, rigu), TextStyle(kBlack, 0.03), TextAlign(0.85, 0.92, 32));
	
		TextDraw(      "Data", TextStyle(kBlack  , 0.03), TextAlign(0.82, 0.83, 32));
		TextDraw("Antiproton", TextStyle(kRed    , 0.03), TextAlign(0.82, 0.79, 32));
		TextDraw("Background", TextStyle(kBlue   , 0.03), TextAlign(0.82, 0.75, 32));
		TextDraw(STR_FMT("#chi^{2}/NDF ( %.2f )", fitNChi), TextStyle(kGreen+2, 0.03), TextAlign(0.82, 0.71, 32));

		TextDraw(STR_FMT("N_{Ap} ( %.0f #pm %.0f )", stdpar.val(0), stdpar.err(0)), TextStyle(kRed , 0.03), TextAlign(0.82, 0.65, 32));
		TextDraw(STR_FMT("N_{Bk} ( %.0f #pm %.0f )", stdpar.val(1), stdpar.err(1)), TextStyle(kBlue, 0.03), TextAlign(0.82, 0.61, 32));
		TextDraw(STR_FMT("#bar{p}/p ( %E )", ratAppc), TextStyle(kGreen+2, 0.03), TextAlign(0.82, 0.57, 32));
		
		editor->save();
	}
    
    COUT("RIG(%8.2f %8.2f) SAMP %14.2f AP %14.2f (%14.2f) BK %14.2f (%14.2f)\n", rigl, rigu, numOfSmp, stdpar.val(0), syspar.val(0), stdpar.val(1), syspar.val(1));
}

class AppRlt {
	public :
		AppRlt(Int_t sbin = -1, Int_t ebin = -1, const std::string& mark = "", Hist * href = 0);
		~AppRlt() { if (fRefCrrAccp) { delete fRefCrrAccp; fRefCrrAccp = 0; } }

		void setPoint(AppFitRlt& fitRlt, Int_t irig, Int_t itme = -1);
		void addRlt(AppRlt& rlt, Int_t sbin = -1, Int_t ebin = -1);
		void addRlt(AppRlt& rltA, AppRlt& rltB, AppRlt& rltC, Int_t sbin = -1, Int_t ebin = -1);

		void fit(const std::string& name, Hist * hPos, Hist * hNeg, Hist * hSig, Hist * hBkg, Int_t itme = -1);
		void tmefit(const std::string& name, Hist * hPos, Hist * hNeg, Hist * hSig, Hist * hBkg);

	protected :
		Int_t fSbin;
		Int_t fEbin;

	public :
		Hist * hNumPr;
		Hist * hNumAp;

		Hist * hErrStat; 
		Hist * hErrCntn; 
		Hist * hErrAccp; 
		Hist * hErrRscl; 
		Hist * hErrSysm; 
		Hist * hErrTotl;
		
		Hist * hFitSN;
		Hist * hFitNChi;
		
		Hist * hRatOrg;
		Hist * hRatStat;
		Hist * hRatTotl;

		Hist * hCrrAccp;

		// Time
		Hist * hTNumPr;
		Hist * hTNumAp;

		Hist * hTErrStat; 
		Hist * hTErrCntn; 
		Hist * hTErrAccp; 
		Hist * hTErrRscl; 
		Hist * hTErrSysm; 
		Hist * hTErrTotl;
		
		Hist * hTFitSN;
		Hist * hTFitNChi;
		
		Hist * hTRatOrg;
		Hist * hTRatStat;
		Hist * hTRatTotl;
		
		Hist * hTCrrAccp;

	public :
		Hist * hRefCrrAccp;
		TF1  * fRefCrrAccp;

	public :
		static Hist * CalXsecErr(const std::string& name, Hist * selC090, Hist * trgC090, Hist * selC110, Hist * trgC110);
		static Hist * CalAccpCrr(const std::string& name, Hist * selPr, Hist * trgPr, Hist * selAp, Hist * trgAp);
		static Hist * CalCutEff(const std::string& name, Hist * cutPr, Hist * cutAp);

	public :
		static Axis AXnr;
		static Axis AXtme;
};

Axis AppRlt::AXnr("|Rigidity| [GV]", BinList(
		{   1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
		    3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
			  8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
			 19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
			 41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
			 93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } ) );
			
//Axis AppRlt::AXnr("|Rigidity| [GV]", BinList(
//	{   1.00,   2.00,   3.00,   4.12,   5.00,   
//	    6.00,   7.10,   8.30,   9.62,  11.04,  
//		 12.59,  14.25,  16.05,  17.98,  20.04,  
//		 22.25,  24.62,  27.25,  30.21  } ) );
			
// 3 month
Axis AppRlt::AXtme("Date", BinList(
		{ 1305417600, 1312416000, 1319414400, 1326412800, 1333411200,
		  1340409600, 1347408000, 1354406400, 1361404800, 1368403200,
			1375401600, 1382400000, 1389398400, 1396396800, 1403395200,
			1410393600, 1417392000, 1424390400, 1431388800, 1438387200,
			1445385600, 1452384000, 1459382400, 1466380800, 1473379200, 
			1480377600 } ) );

AppRlt::AppRlt(Int_t sbin, Int_t ebin, const std::string& mark, Hist * href) : fSbin(0), fEbin(0), hNumPr(0), hNumAp(0), hErrStat(0), hErrCntn(0), hErrAccp(0), hErrRscl(0), hErrSysm(0), hErrTotl(0), hFitNChi(0), hRatOrg(0), hRatStat(0), hRatTotl(0), hCrrAccp(0), hRefCrrAccp(href), fRefCrrAccp(0) {
	fSbin = (sbin == -1 || sbin <           1) ?           1 : sbin;
	fEbin = (ebin == -1 || ebin > AXnr.nbin()) ? AXnr.nbin() : ebin;

	hNumPr   = Hist::New(STR_FMT("h%s_NumPr"  , mark.c_str()), "", HistAxis(AXnr));
	hNumAp   = Hist::New(STR_FMT("h%s_NumAp"  , mark.c_str()), "", HistAxis(AXnr));

	hErrStat = Hist::New(STR_FMT("h%s_ErrStat", mark.c_str()), "", HistAxis(AXnr));
	hErrCntn = Hist::New(STR_FMT("h%s_ErrCntn", mark.c_str()), "", HistAxis(AXnr));
	hErrAccp = Hist::New(STR_FMT("h%s_ErrAccp", mark.c_str()), "", HistAxis(AXnr));
	hErrRscl = Hist::New(STR_FMT("h%s_ErrRscl", mark.c_str()), "", HistAxis(AXnr));
	hErrSysm = Hist::New(STR_FMT("h%s_ErrSysm", mark.c_str()), "", HistAxis(AXnr));
	hErrTotl = Hist::New(STR_FMT("h%s_ErrTotl", mark.c_str()), "", HistAxis(AXnr));
	
	hFitSN   = Hist::New(STR_FMT("h%s_FitSN"  , mark.c_str()), "", HistAxis(AXnr));
	hFitNChi = Hist::New(STR_FMT("h%s_FitNChi", mark.c_str()), "", HistAxis(AXnr));
	
	hRatOrg  = Hist::New(STR_FMT("h%s_RatOrg" , mark.c_str()), "", HistAxis(AXnr));
	hRatStat = Hist::New(STR_FMT("h%s_RatStat", mark.c_str()), "", HistAxis(AXnr));
	hRatTotl = Hist::New(STR_FMT("h%s_RatTotl", mark.c_str()), "", HistAxis(AXnr));
	
	hCrrAccp = Hist::New(STR_FMT("h%s_CrrAccp", mark.c_str()), "", HistAxis(AXnr));
	
	// Time
	hTNumPr   = Hist::New(STR_FMT("hT%s_NumPr"  , mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTNumAp   = Hist::New(STR_FMT("hT%s_NumAp"  , mark.c_str()), "", HistAxis(AXtme, AXnr));

	hTErrStat = Hist::New(STR_FMT("hT%s_ErrStat", mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTErrCntn = Hist::New(STR_FMT("hT%s_ErrCntn", mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTErrAccp = Hist::New(STR_FMT("hT%s_ErrAccp", mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTErrRscl = Hist::New(STR_FMT("hT%s_ErrRscl", mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTErrSysm = Hist::New(STR_FMT("hT%s_ErrSysm", mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTErrTotl = Hist::New(STR_FMT("hT%s_ErrTotl", mark.c_str()), "", HistAxis(AXtme, AXnr));
	
	hTFitSN   = Hist::New(STR_FMT("hT%s_FitSN"  , mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTFitNChi = Hist::New(STR_FMT("hT%s_FitNChi", mark.c_str()), "", HistAxis(AXtme, AXnr));
	
	hTRatOrg  = Hist::New(STR_FMT("hT%s_RatOrg" , mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTRatStat = Hist::New(STR_FMT("hT%s_RatStat", mark.c_str()), "", HistAxis(AXtme, AXnr));
	hTRatTotl = Hist::New(STR_FMT("hT%s_RatTotl", mark.c_str()), "", HistAxis(AXtme, AXnr));
	
	hTCrrAccp = Hist::New(STR_FMT("hT%s_CrrAccp", mark.c_str()), "", HistAxis(AXtme, AXnr));

	SetTimeAxis((*hTNumPr)()->GetXaxis());
	
	SetTimeAxis((*hTErrStat)()->GetXaxis());
	SetTimeAxis((*hTErrCntn)()->GetXaxis());
	SetTimeAxis((*hTErrAccp)()->GetXaxis());
	SetTimeAxis((*hTErrRscl)()->GetXaxis());
	SetTimeAxis((*hTErrSysm)()->GetXaxis());
	SetTimeAxis((*hTErrTotl)()->GetXaxis());
	
	SetTimeAxis((*hTFitSN)()->GetXaxis());
	SetTimeAxis((*hTFitNChi)()->GetXaxis());
	
	SetTimeAxis((*hTRatOrg)() ->GetXaxis());
	SetTimeAxis((*hTRatStat)()->GetXaxis());
	SetTimeAxis((*hTRatTotl)()->GetXaxis());
	
	SetTimeAxis((*hTCrrAccp)()->GetXaxis());

	// Acceptance correction
	if (hRefCrrAccp) {
		fRefCrrAccp = new TF1(CSTR_FMT("fRefCrrAccp_%s", mark.c_str()), "[0] + [1]*TMath::Power(x, [2]) + [3]*TMath::Power(x, [4])", 1, 1800);
		fRefCrrAccp->SetParameters(1.05, 0.3, -0.5, 0.1, -0.5);
		(*hRefCrrAccp)()->Fit(fRefCrrAccp, "qn");
	}
}

void AppRlt::setPoint(AppFitRlt& fitRlt, Int_t irig, Int_t itme) {
	if (itme == -1) {
		(*hNumPr)()->SetBinContent(irig, fitRlt.numPr);
		(*hNumAp)()->SetBinContent(irig, fitRlt.numAp);
		
		(*hErrStat)()->SetBinContent(irig, fitRlt.errStat);
		(*hErrCntn)()->SetBinContent(irig, fitRlt.errCntn);
		(*hErrAccp)()->SetBinContent(irig, fitRlt.errAccp);
		(*hErrRscl)()->SetBinContent(irig, fitRlt.errRscl);
		(*hErrSysm)()->SetBinContent(irig, fitRlt.errSysm);
		(*hErrTotl)()->SetBinContent(irig, fitRlt.errTotl);
		
		(*hFitSN  )()->SetBinContent(irig, fitRlt.fitSN);
		(*hFitNChi)()->SetBinContent(irig, fitRlt.fitNChi);
		
		(*hRatOrg)()->SetBinContent(irig, fitRlt.orgAppr);
		(*hRatOrg)()->SetBinError  (irig, fitRlt.orgAppe);
		
		(*hRatOrg)()->SetBinContent(irig, fitRlt.orgAppr);
		(*hRatOrg)()->SetBinError  (irig, fitRlt.orgAppe);

		(*hRatStat)()->SetBinContent(irig, fitRlt.ratAppc);
		(*hRatStat)()->SetBinError  (irig, fitRlt.errStat);
		
		(*hRatTotl)()->SetBinContent(irig, fitRlt.ratAppc);
		(*hRatTotl)()->SetBinError  (irig, fitRlt.errTotl);
		
		(*hCrrAccp)()->SetBinContent(irig, fitRlt.crrAccp);
	}
	else {
		(*hTNumPr)()->SetBinContent(itme, irig, fitRlt.numPr);
		(*hTNumAp)()->SetBinContent(itme, irig, fitRlt.numAp);
		
		(*hTErrStat)()->SetBinContent(itme, irig, fitRlt.errStat);
		(*hTErrCntn)()->SetBinContent(itme, irig, fitRlt.errCntn);
		(*hTErrAccp)()->SetBinContent(itme, irig, fitRlt.errAccp);
		(*hTErrRscl)()->SetBinContent(itme, irig, fitRlt.errRscl);
		(*hTErrSysm)()->SetBinContent(itme, irig, fitRlt.errSysm);
		(*hTErrTotl)()->SetBinContent(itme, irig, fitRlt.errTotl);
		
		(*hTFitSN  )()->SetBinContent(itme, irig, fitRlt.fitSN);
		(*hTFitNChi)()->SetBinContent(itme, irig, fitRlt.fitNChi);
		
		(*hTRatOrg)()->SetBinContent(itme, irig, fitRlt.orgAppr);
		(*hTRatOrg)()->SetBinError  (itme, irig, fitRlt.orgAppe);
		
		(*hTRatStat)()->SetBinContent(itme, irig, fitRlt.ratAppc);
		(*hTRatStat)()->SetBinError  (itme, irig, fitRlt.errStat);
		
		(*hTRatTotl)()->SetBinContent(itme, irig, fitRlt.ratAppc);
		(*hTRatTotl)()->SetBinError  (itme, irig, fitRlt.errTotl);
		
		(*hTCrrAccp)()->SetBinContent(itme, irig, fitRlt.crrAccp);
	}
}
		
void AppRlt::addRlt(AppRlt& rlt, Int_t sbin, Int_t ebin) {
	Int_t _sbin = (sbin == -1 || sbin <           1) ?           1 : sbin;
	Int_t _ebin = (ebin == -1 || ebin > AXnr.nbin()) ? AXnr.nbin() : ebin;

	for (Int_t irig = _sbin; irig <= _ebin; ++irig) {
		(*hNumPr)()->SetBinContent(irig, (*rlt.hNumPr)()->GetBinContent(irig));
		(*hNumAp)()->SetBinContent(irig, (*rlt.hNumAp)()->GetBinContent(irig));
		
		(*hErrStat)()->SetBinContent(irig, (*rlt.hErrStat)()->GetBinContent(irig));
		(*hErrCntn)()->SetBinContent(irig, (*rlt.hErrCntn)()->GetBinContent(irig));
		(*hErrAccp)()->SetBinContent(irig, (*rlt.hErrAccp)()->GetBinContent(irig));
		(*hErrRscl)()->SetBinContent(irig, (*rlt.hErrRscl)()->GetBinContent(irig));
		(*hErrSysm)()->SetBinContent(irig, (*rlt.hErrSysm)()->GetBinContent(irig));
		(*hErrTotl)()->SetBinContent(irig, (*rlt.hErrTotl)()->GetBinContent(irig));
		                                                 
		(*hFitSN  )()->SetBinContent(irig, (*rlt.hFitSN  )()->GetBinContent(irig));
		(*hFitNChi)()->SetBinContent(irig, (*rlt.hFitNChi)()->GetBinContent(irig));
		                                                 
		(*hRatOrg)()->SetBinContent(irig, (*rlt.hRatOrg)()->GetBinContent(irig));
		(*hRatOrg)()->SetBinError  (irig, (*rlt.hRatOrg)()->GetBinError  (irig));
		
		(*hRatStat)()->SetBinContent(irig, (*rlt.hRatStat)()->GetBinContent(irig));
		(*hRatStat)()->SetBinError  (irig, (*rlt.hRatStat)()->GetBinError  (irig));
		                                                 
		(*hRatTotl)()->SetBinContent(irig, (*rlt.hRatTotl)()->GetBinContent(irig));
		(*hRatTotl)()->SetBinError  (irig, (*rlt.hRatTotl)()->GetBinError  (irig));
		                                                 
		(*hCrrAccp)()->SetBinContent(irig, (*rlt.hCrrAccp)()->GetBinContent(irig));
		
		for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {	
			(*hTNumPr)()->SetBinContent(irig, (*rlt.hTNumPr)()->GetBinContent(irig));
			(*hTNumAp)()->SetBinContent(irig, (*rlt.hTNumAp)()->GetBinContent(irig));
			
			(*hTErrStat)()->SetBinContent(itme, irig, (*rlt.hTErrStat)()->GetBinContent(itme, irig));
			(*hTErrCntn)()->SetBinContent(itme, irig, (*rlt.hTErrCntn)()->GetBinContent(itme, irig));
			(*hTErrAccp)()->SetBinContent(itme, irig, (*rlt.hTErrAccp)()->GetBinContent(itme, irig));
			(*hTErrRscl)()->SetBinContent(itme, irig, (*rlt.hTErrRscl)()->GetBinContent(itme, irig));
			(*hTErrSysm)()->SetBinContent(itme, irig, (*rlt.hTErrSysm)()->GetBinContent(itme, irig));
			(*hTErrTotl)()->SetBinContent(itme, irig, (*rlt.hTErrTotl)()->GetBinContent(itme, irig));
			
			(*hTFitSN  )()->SetBinContent(itme, irig, (*rlt.hTFitSN  )()->GetBinContent(itme, irig));
			(*hTFitNChi)()->SetBinContent(itme, irig, (*rlt.hTFitNChi)()->GetBinContent(itme, irig));
			
			(*hTRatOrg)()->SetBinContent(itme, irig, (*rlt.hTRatOrg)()->GetBinContent(itme, irig));
			(*hTRatOrg)()->SetBinError  (itme, irig, (*rlt.hTRatOrg)()->GetBinError  (itme, irig));
			
			(*hTRatStat)()->SetBinContent(itme, irig, (*rlt.hTRatStat)()->GetBinContent(itme, irig));
			(*hTRatStat)()->SetBinError  (itme, irig, (*rlt.hTRatStat)()->GetBinError  (itme, irig));
			
			(*hTRatTotl)()->SetBinContent(itme, irig, (*rlt.hTRatTotl)()->GetBinContent(itme, irig));
			(*hTRatTotl)()->SetBinError  (itme, irig, (*rlt.hTRatTotl)()->GetBinError  (itme, irig));
			
			(*hTCrrAccp)()->SetBinContent(itme, irig, (*rlt.hTCrrAccp)()->GetBinContent(itme, irig));
		}
	}
}


void AppRlt::addRlt(AppRlt& rltA, AppRlt& rltB, AppRlt& rltC, Int_t sbin, Int_t ebin) {
	Int_t _sbin = (sbin == -1 || sbin <           1) ?           1 : sbin;
	Int_t _ebin = (ebin == -1 || ebin > AXnr.nbin()) ? AXnr.nbin() : ebin;
	
	for (Int_t irig = _sbin; irig <= _ebin; ++irig) {
		Double_t numPrA = (*rltA.hNumPr)()->GetBinContent(irig);
		Double_t numPrB = (*rltB.hNumPr)()->GetBinContent(irig);
		Double_t numPrC = (*rltC.hNumPr)()->GetBinContent(irig);
		Double_t numPrT = (numPrA + numPrB + numPrC);
		Double_t numApA = (*rltA.hNumAp)()->GetBinContent(irig);
		Double_t numApB = (*rltB.hNumAp)()->GetBinContent(irig);
		Double_t numApC = (*rltC.hNumAp)()->GetBinContent(irig);
		Double_t numApT = (numApA + numApB + numApC);
		Double_t orgAppr = (numApT / numPrT);
		Double_t orgAppe = std::sqrt((*rltA.hErrStat)()->GetBinContent(irig)*(*rltA.hErrStat)()->GetBinContent(irig) +
				                     (*rltB.hErrStat)()->GetBinContent(irig)*(*rltB.hErrStat)()->GetBinContent(irig) +
				                     (*rltC.hErrStat)()->GetBinContent(irig)*(*rltC.hErrStat)()->GetBinContent(irig));
		Double_t crrAccp = ((*rltA.hCrrAccp)()->GetBinContent(irig) * numPrA +
				            (*rltB.hCrrAccp)()->GetBinContent(irig) * numPrB +
				            (*rltC.hCrrAccp)()->GetBinContent(irig) * numPrC) / numPrT;
		Double_t ratStat = ((*rltA.hRatStat)()->GetBinContent(irig) * numPrA +
				            (*rltB.hRatStat)()->GetBinContent(irig) * numPrB +
				            (*rltC.hRatStat)()->GetBinContent(irig) * numPrC) / numPrT;
		Double_t ratTotl = ((*rltA.hRatTotl)()->GetBinContent(irig) * numPrA +
				            (*rltB.hRatTotl)()->GetBinContent(irig) * numPrB +
				            (*rltC.hRatTotl)()->GetBinContent(irig) * numPrC) / numPrT;
		Double_t errStat = std::sqrt((*rltA.hErrStat)()->GetBinContent(irig)*(*rltA.hErrStat)()->GetBinContent(irig) * numPrB*numPrB +
				                     (*rltB.hErrStat)()->GetBinContent(irig)*(*rltB.hErrStat)()->GetBinContent(irig) * numPrB*numPrB +
				                     (*rltC.hErrStat)()->GetBinContent(irig)*(*rltC.hErrStat)()->GetBinContent(irig) * numPrC*numPrC) / numPrT;
		Double_t errCntn = std::sqrt((*rltA.hErrCntn)()->GetBinContent(irig)*(*rltA.hErrCntn)()->GetBinContent(irig) * numPrB*numPrB +
				                     (*rltB.hErrCntn)()->GetBinContent(irig)*(*rltB.hErrCntn)()->GetBinContent(irig) * numPrB*numPrB +
				                     (*rltC.hErrCntn)()->GetBinContent(irig)*(*rltC.hErrCntn)()->GetBinContent(irig) * numPrC*numPrC) / numPrT;
		//Double_t errAccp = std::sqrt((*rltA.hErrAccp)()->GetBinContent(irig)*(*rltA.hErrAccp)()->GetBinContent(irig) * numPrA*numPrA +
		//		                     (*rltB.hErrAccp)()->GetBinContent(irig)*(*rltB.hErrAccp)()->GetBinContent(irig) * numPrB*numPrB +
		//		                     (*rltC.hErrAccp)()->GetBinContent(irig)*(*rltC.hErrAccp)()->GetBinContent(irig) * numPrC*numPrC) / numPrT;
		Double_t errAccp = ((*rltA.hErrAccp)()->GetBinContent(irig) * numPrA +
				            (*rltB.hErrAccp)()->GetBinContent(irig) * numPrB +
				            (*rltC.hErrAccp)()->GetBinContent(irig) * numPrC) / numPrT;
		//Double_t errRscl = std::sqrt((*rltA.hErrRscl)()->GetBinContent(irig)*(*rltA.hErrRscl)()->GetBinContent(irig) +
		//		                     (*rltB.hErrRscl)()->GetBinContent(irig)*(*rltB.hErrRscl)()->GetBinContent(irig) +
		//		                     (*rltC.hErrRscl)()->GetBinContent(irig)*(*rltC.hErrRscl)()->GetBinContent(irig));
		Double_t errRscl = ((*rltA.hErrRscl)()->GetBinContent(irig) * numPrA +
				            (*rltB.hErrRscl)()->GetBinContent(irig) * numPrB +
				            (*rltC.hErrRscl)()->GetBinContent(irig) * numPrC) / numPrT;
		Double_t errSysm = std::sqrt(errCntn*errCntn + errAccp*errAccp + errRscl*errRscl);
		Double_t errTotl = std::sqrt(errStat*errStat + errSysm*errSysm);
		//Double_t errSysm = std::sqrt((*rltA.hErrSysm)()->GetBinContent(irig)*(*rltA.hErrSysm)()->GetBinContent(irig) +
		//		                     (*rltB.hErrSysm)()->GetBinContent(irig)*(*rltB.hErrSysm)()->GetBinContent(irig) +
		//		                     (*rltC.hErrSysm)()->GetBinContent(irig)*(*rltC.hErrSysm)()->GetBinContent(irig));
		//Double_t errTotl = std::sqrt((*rltA.hErrTotl)()->GetBinContent(irig)*(*rltA.hErrTotl)()->GetBinContent(irig) +
		//		                     (*rltB.hErrTotl)()->GetBinContent(irig)*(*rltB.hErrTotl)()->GetBinContent(irig) +
		//		                     (*rltC.hErrTotl)()->GetBinContent(irig)*(*rltC.hErrTotl)()->GetBinContent(irig));

		(*hNumPr)()->SetBinContent(irig, numPrT);
		(*hNumAp)()->SetBinContent(irig, numApT);

		(*hErrStat)()->SetBinContent(irig, errStat);
		(*hErrCntn)()->SetBinContent(irig, errCntn);
		(*hErrAccp)()->SetBinContent(irig, errAccp);
		(*hErrRscl)()->SetBinContent(irig, errRscl);
		(*hErrSysm)()->SetBinContent(irig, errSysm);
		(*hErrTotl)()->SetBinContent(irig, errTotl);
		                                                 
		(*hRatOrg)()->SetBinContent(irig, orgAppr);
		(*hRatOrg)()->SetBinError  (irig, orgAppe);
		
		(*hRatStat)()->SetBinContent(irig, ratStat);
		(*hRatStat)()->SetBinError  (irig, errStat);
		                                                 
		(*hRatTotl)()->SetBinContent(irig, ratTotl);
		(*hRatTotl)()->SetBinError  (irig, errTotl);
		                                                 
		(*hCrrAccp)()->SetBinContent(irig, crrAccp);
	}
}


void AppRlt::fit(const std::string& name, Hist * hPos, Hist * hNeg, Hist * hSig, Hist * hBkg, Int_t itme) {
	std::vector<Hist *>&& hPosVec = Hist::ProjectAll(HistProj::kY, hPos);
	std::vector<Hist *>&& hNegVec = Hist::ProjectAll(HistProj::kY, hNeg);
	std::vector<Hist *>&& hSigVec = Hist::ProjectAll(HistProj::kY, hSig);
	std::vector<Hist *>&& hBkgVec = Hist::ProjectAll(HistProj::kY, hBkg);
	
	for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
		Double_t prnum = (*hPosVec.at(irig))()->Integral();
		if (itme == -1) (*hNumPr)()->SetBinContent(irig, prnum);
		else            (*hTNumPr)()->SetBinContent(itme, irig, prnum);
	}

	for (Int_t irig = fSbin; irig <= fEbin; ++irig) {
		AppFitRlt fitRlt
			(name, hPosVec.at(irig), hNegVec.at(irig), 
			       hSigVec.at(irig), hBkgVec.at(irig),
						 AXnr(irig-1), AXnr(irig),
						 (fRefCrrAccp?fRefCrrAccp->Eval(AXnr.center(irig, AxisScale::kLog)):1.0),
						 //((hRefCrrAccp)?(*hRefCrrAccp)()->GetBinContent(irig):1.0),
						 itme );
		if (!fitRlt.isSucc) continue;
		setPoint(fitRlt, irig, itme);
	}
	
	Hist::Delete(hPosVec);
	Hist::Delete(hNegVec);
	Hist::Delete(hSigVec);
	Hist::Delete(hBkgVec);
}
		
void AppRlt::tmefit(const std::string& name, Hist * hPos, Hist * hNeg, Hist * hSig, Hist * hBkg) {
	std::vector<Hist *>&& hPosVec = Hist::ProjectAll(HistProj::kZY, hPos);
	std::vector<Hist *>&& hNegVec = Hist::ProjectAll(HistProj::kZY, hNeg);
	std::vector<Hist *>&& hSigVec = Hist::ProjectAll(HistProj::kZY, hSig);
	std::vector<Hist *>&& hBkgVec = Hist::ProjectAll(HistProj::kZY, hBkg);
	
	for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
		Hist * _sig = (hSig->dim() == HistDim::k3D) ? hSigVec.at(itme) : hSig;
		Hist * _bkg = (hBkg->dim() == HistDim::k3D) ? hBkgVec.at(itme) : hBkg;
		fit(name, hPosVec.at(itme), hNegVec.at(itme), _sig, _bkg, itme);
	}

	Hist::Delete(hPosVec);
	Hist::Delete(hNegVec);
	if (hSig->dim() == HistDim::k3D) Hist::Delete(hSigVec);
	if (hBkg->dim() == HistDim::k3D) Hist::Delete(hBkgVec);
}


Hist * AppRlt::CalXsecErr(const std::string& name, Hist * selC090, Hist * trgC090, Hist * selC110, Hist * trgC110) {
	Hist * hXsec = Hist::New(name, "", HistAxis(AXnr));
	for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
		Double_t sel[2] = { (*selC090)()->GetBinContent(irig),  (*selC110)()->GetBinContent(irig) };
		Double_t trg[2] = { (*trgC090)()->GetBinContent(irig),  (*trgC110)()->GetBinContent(irig) };
		if (MGNumc::EqualToZero(sel[0]) || MGNumc::EqualToZero(sel[1])) continue;
		if (MGNumc::EqualToZero(trg[0]) || MGNumc::EqualToZero(trg[1])) continue;
		if (sel[0] < 1000 || sel[1] < 1000) continue;
		Double_t rat[2] = { (sel[0]/trg[0]), (sel[1]/trg[1]) };
		Double_t raterr = std::sqrt(rat[0]/trg[0] + rat[1]/trg[1]);
		Double_t ratmp[2] = { std::fabs(rat[0]-rat[1]), std::fabs(rat[0]+rat[1]) };
		Double_t crr = (ratmp[0] / ratmp[1]);
		Double_t err = crr * raterr * std::sqrt(1./ratmp[0]/ratmp[0] + 1./ratmp[1]/ratmp[1]);
		if ((err/crr) > 0.5 || !MGNumc::Valid(crr)) continue;
		(*hXsec)()->SetBinContent(irig, crr);
		(*hXsec)()->SetBinError  (irig, err);
	}
	return hXsec;
}


Hist * AppRlt::CalAccpCrr(const std::string& name, Hist * selPr, Hist * trgPr, Hist * selAp, Hist * trgAp) {
	TF1 * UnfoldRatio = new TF1("UnfoldRatio", "1.0+[0]*TMath::Power([1]*(x+[2]), [3])", 1, 1000);
	UnfoldRatio->SetParameters(-6.29727e+02, 3.49266e-02, 3.26913e+01, -5.17749e+01);

	Hist * hAccp = Hist::New(name, "", HistAxis(AXnr));
	for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
		Double_t sel[2] = { (*selPr)()->GetBinContent(irig),  (*selAp)()->GetBinContent(irig) };
		Double_t trg[2] = { (*trgPr)()->GetBinContent(irig),  (*trgAp)()->GetBinContent(irig) };
		if (MGNumc::EqualToZero(sel[0]) || MGNumc::EqualToZero(sel[1])) continue;
		if (MGNumc::EqualToZero(trg[0]) || MGNumc::EqualToZero(trg[1])) continue;
		if (sel[0] < 5000 || sel[1] < 5000) continue;
		Double_t ratPr = (sel[0]/trg[0]);
		Double_t ratAp = (sel[1]/trg[1]);
		Double_t crr   = ratPr / ratAp;
		Double_t err   = crr * std::sqrt(1.0/sel[0] + 1.0/sel[1]);
		if ((err/crr) > 0.1 || !MGNumc::Valid(crr)) continue;
		Double_t unfRat = UnfoldRatio->Eval(AXnr.center(irig, AxisScale::kLog));
		(*hAccp)()->SetBinContent(irig, unfRat * crr);
		(*hAccp)()->SetBinError  (irig, unfRat * err);
	}
	
	delete UnfoldRatio;
	return hAccp;
}


Hist * AppRlt::CalCutEff(const std::string& name, Hist * cutPr, Hist * cutAp) {
	Axis AXeff("Eff", cutPr->axis().z().nbin()-1, 1., (cutPr->axis().z().nbin()));
	Hist * hPrEff = Hist::New(STR_FMT("%s_PR", name.c_str()), "", HistAxis(AXtme, AXnr, AXeff));
	Hist * hApEff = Hist::New(STR_FMT("%s_AP", name.c_str()), "", HistAxis(AXtme, AXnr, AXeff));
	Hist * hEff = Hist::New(name.c_str(), "", HistAxis(AXtme, AXnr, AXeff));

	(*hPrEff)()->GetXaxis()->SetTimeDisplay(1);
	(*hPrEff)()->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}%F1970-01-01 00:00:00s0");
	(*hApEff)()->GetXaxis()->SetTimeDisplay(1);
	(*hApEff)()->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}%F1970-01-01 00:00:00s0");
	(*hEff)()->GetXaxis()->SetTimeDisplay(1);
	(*hEff)()->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}%F1970-01-01 00:00:00s0");

	for (Int_t ieff = 1; ieff <= AXeff.nbin(); ++ieff) {
		for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
			std::vector<Float_t> prefvec;
			std::vector<Float_t> prsgvec;
			std::vector<Float_t> apefvec;
			std::vector<Float_t> apsgvec;
			std::vector<Float_t> appefvec;
			std::vector<Float_t> appsgvec;
			for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
				Float_t prtw = (*cutPr)()->GetBinContent(itme, irig, ieff);
				Float_t prpw = (*cutPr)()->GetBinContent(itme, irig, ieff+1);
				Float_t pref = prpw/ prtw;
				Float_t prsg = std::sqrt(pref * (1.0 - pref) / prtw);
				
				Float_t aptw = (*cutAp)()->GetBinContent(itme, irig, ieff);
				Float_t appw = (*cutAp)()->GetBinContent(itme, irig, ieff+1);
				Float_t apef = appw/ aptw;
				Float_t apsg = std::sqrt(apef * (1.0 - apef) / aptw);

				Float_t appef = apef / pref;
				Float_t appsg = appef * std::sqrt((prsg/pref)*(prsg/pref) + (apsg/apef)*(apsg/apef));

				prefvec.push_back(pref);
				prsgvec.push_back(prsg);
				apefvec.push_back(apef);
				apsgvec.push_back(apsg);
				appefvec.push_back(appef);
				appsgvec.push_back(appsg);
			}
			Float_t prmn = 0.;
			Float_t apmn = 0.;
			Float_t appmn = 0.;
			for (Int_t it = 0; it < appefvec.size(); ++it) {
				prmn += prefvec.at(it);
				apmn += apefvec.at(it);
				appmn += appefvec.at(it);
			}
			prmn /= prefvec.size();
			apmn /= apefvec.size();
			appmn /= appefvec.size();

			for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
				//(*hPrEff)()->SetBinContent(itme, irig, ieff, prefvec.at(itme-1)/prmn);
				//(*hPrEff)()->SetBinError  (itme, irig, ieff, prsgvec.at(itme-1)/prmn);
				//(*hApEff)()->SetBinContent(itme, irig, ieff, apefvec.at(itme-1)/apmn);
				//(*hApEff)()->SetBinError  (itme, irig, ieff, apsgvec.at(itme-1)/apmn);
				(*hPrEff)()->SetBinContent(itme, irig, ieff, prefvec.at(itme-1));
				(*hPrEff)()->SetBinError  (itme, irig, ieff, prsgvec.at(itme-1));
				(*hApEff)()->SetBinContent(itme, irig, ieff, apefvec.at(itme-1));
				(*hApEff)()->SetBinError  (itme, irig, ieff, apsgvec.at(itme-1));
				(*hEff)()->SetBinContent(itme, irig, ieff, appefvec.at(itme-1)/appmn);
				(*hEff)()->SetBinError  (itme, irig, ieff, appsgvec.at(itme-1)/appmn);
			}
		}
	}
	
	Axis AXcut("Cut", cutPr->axis().z().nbin(), 0., (cutPr->axis().z().nbin()));
	Hist * hRat = Hist::New(CSTR_FMT("%s_PreRatio", name.c_str()), "", HistAxis(AXtme, AXnr, AXcut));
	(*hRat)()->GetXaxis()->SetTimeDisplay(1);
	(*hRat)()->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}%F1970-01-01 00:00:00s0");
	for (Int_t icut = 1; icut <= AXcut.nbin(); ++icut) {
		for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
			for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
				Float_t pr = (*cutPr)()->GetBinContent(itme, irig, icut);
				Float_t ap = (*cutAp)()->GetBinContent(itme, irig, icut);
				Float_t val = (ap / pr);
				Float_t err = val * std::sqrt(1./ap+1./pr);
				(*hRat)()->SetBinContent(itme, irig, icut, val);
				(*hRat)()->SetBinError  (itme, irig, icut, err);
			}
		}
	}

	return hEff;
}

#endif // __AppRlt_H__
