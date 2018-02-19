#ifndef __YiMergeAna_C__
#define __YiMergeAna_C__

// User defination library
#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/STL/stl.h"
#include "/afs/cern.ch/user/h/hchou/libraries/CPlusPlus/ROOT/selflib.h"

using namespace MgntROOT;
	
void SetTimeAxis(TAxis * axis) {
	axis->SetTitle("Time");
	axis->SetTimeDisplay(1);
	axis->SetTimeFormat("#splitline{%b}{%Y}%F1970-01-01 00:00:00s0");
	//axis->SetNdivisions(-216);
	axis->SetLabelOffset(0.02);
	axis->SetTitleOffset(1.5);
	axis->SetTitleFont(42);
}

int main(int argc, const char ** argv) {
	Style::LoadDefaultEnvironment();
	//TString sbn = "/afs/cern.ch/user/h/hchou/public/DATABASE/physics/binning/antipp3_rig.root";
	//TFile * fBin = TFile::Open(sbn);
	//TH1D  * hBn0 = (TH1D*)fBin->Get("hbin0");
	//TH1D  * hBn1 = (TH1D*)fBin->Get("hbin1");
	//Axis AXnr = Axis(hBn0, Axis::kX);
  //Axis AXir = Axis(hBn1, Axis::kX);
	//fBin->Close();
	
	// rigidity binning
	//Axis AXnr = Axis("Rigidity [GV]",
	//	{   0.50,   0.80, // extern bins
	//	    1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
	//	    3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
	//		  8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
	//		 19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
	//		 41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
	//		 93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00, 
	//		800.00 } ); // extern bin
	//Axis AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
	
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
			
	Axis AXnr = Axis("Rigidity [GV]",
		{   1.00,   1.46,   2.00,   3.00,   4.12, 
		    5.00,   6.00,   7.10,   8.30,   9.62, 
		   11.04,  12.59,  14.25,  16.05,  17.98, 
		   20.04,  22.25,  24.62,  27.25,  30.21, 
		   35.36,  40.00 } );
	Axis AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
	
	// 6 month
	//Axis AXtme = Axis("Time", 
	//	{ 1312416000, 1326412800, 1340409600, 1354406400, 1368403200, 
	//	  1382400000, 1396396800, 1410393600 ,1424390400, 1438387200, 
	//		1452384000, 1464299612, 1480377600 } );
			
	// 3 month
	Axis AXtme = Axis("Time", 
		{ 1305417600, 1312416000, 1319414400, 1326412800, 1333411200,
		  1340409600, 1347408000, 1354406400, 1361404800, 1368403200,
			1375401600, 1382400000, 1389398400, 1396396800, 1403395200,
			1410393600, 1417392000, 1424390400, 1431388800, 1438387200,
			1445385600, 1452384000, 1459382400, 1466380800, 1473379200, 
			1480377600 } );
	
	std::string tagdir = "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/rlt7";
	std::string goddir = "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/rlt7";
	
	PdfEditor editor(PdfEditor::kSlice, "YiFin", goddir);
	
	const Int_t MaxTme = 100;
	TFile * fTmeVec[MaxTme] = { nullptr };
	for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
		fTmeVec[itme] = TFile::Open(CStrFmt("%s/YiAna%03d.root", tagdir.c_str(), itme));
	}

	// Efficiency
	const Int_t MaxEff = 31;
	Int_t  nCOMEff = (((TH2D*)(fTmeVec[1]->Get("hCOMp_Eff")))->GetNbinsY());
	Hist * hCOMp_Eff[MaxEff] = { nullptr };
	Hist * hCOMn_Eff[MaxEff] = { nullptr };
	for (Int_t ieff = 1; ieff <= nCOMEff; ++ieff) {
		hCOMp_Eff[ieff] = Hist::New(StrFmt("hCOMp_Eff%02d", ieff), StrFmt("COMp(%02d) Efficiency Rig vs Time", ieff), AXnr, AXtme);
		hCOMn_Eff[ieff] = Hist::New(StrFmt("hCOMn_Eff%02d", ieff), StrFmt("COMn(%02d) Efficiency Rig vs Time", ieff), AXnr, AXtme);
		for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
			TH2D * histp = ((TH2D*)(fTmeVec[itme]->Get("hCOMp_Eff")));
			TH2D * histn = ((TH2D*)(fTmeVec[itme]->Get("hCOMn_Eff")));
			for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
				Float_t conp = histp->GetBinContent(irig, ieff);
				Float_t errp = histp->GetBinError  (irig, ieff);
				hCOMp_Eff[ieff]->setContentWithError(conp, errp, irig, itme);
				Float_t conn = histn->GetBinContent(irig, ieff);
				Float_t errn = histn->GetBinError  (irig, ieff);
				hCOMn_Eff[ieff]->setContentWithError(conn, errn, irig, itme);
			}
		}
		//(*hCOMp_Eff[ieff])()->SetMaximum(1.0001);
		//(*hCOMp_Eff[ieff])()->SetMinimum(1.0e-2);
		SetTimeAxis((*hCOMp_Eff[ieff])()->GetYaxis());
		//(*hCOMn_Eff[ieff])()->SetMaximum(1.0001);
		//(*hCOMn_Eff[ieff])()->SetMinimum(1.0e-2);
		SetTimeAxis((*hCOMn_Eff[ieff])()->GetYaxis());
	}



	Hist::VLIST hCOMColl;
	Double_t max = 0.0, min = 1.0;
	for (Int_t ieff = 1; ieff <= nCOMEff; ++ieff) {
		std::vector<Hist *>&& hCOMp_ProjY = Hist::Project(Hist::kProjY, hCOMp_Eff[ieff]);
		std::vector<Hist *>&& hCOMn_ProjY = Hist::Project(Hist::kProjY, hCOMn_Eff[ieff]);
		for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
			Float_t cen = AXnr.center(irig, Axis::kLog);
			SetTimeAxis((*hCOMp_ProjY.at(irig))()->GetXaxis());
			SetTimeAxis((*hCOMn_ProjY.at(irig))()->GetXaxis());
			(*hCOMp_ProjY.at(irig))()->SetMaximum(1.01);
			(*hCOMp_ProjY.at(irig))()->SetMinimum(0.80);
			hCOMp_ProjY.at(irig)->setStyle(Area(), Line(kRed), Marker(kRed));
			hCOMn_ProjY.at(irig)->setStyle(Area(), Line(kBlue), Marker(kBlue));
			if (ieff == 1) continue;
			if (cen < 2 || cen > 3) continue;
			editor.create();
			editor.cd(1, Canvas::AxisScl_t(0, 0));
			hCOMp_ProjY.at(irig)->draw("hist");
			hCOMn_ProjY.at(irig)->draw("same hist");
			Text::Draw(Text::Txt_t(StrFmt("COM(%02d) Time vs Efficiency ( %6.2f ~ %6.2f GV )", ieff, AXnr.bins(irig-1), AXnr.bins(irig)), kBlack, 0.03), Text::Att_t(0.50, 0.92, 22));
			editor.save();
			hCOMColl.push_back(hCOMp_ProjY.at(irig));
			if ((*(hCOMp_ProjY.at(irig)))()->GetMaximum() > max) max = (*(hCOMp_ProjY.at(irig)))()->GetMaximum();
			if ((*(hCOMp_ProjY.at(irig)))()->GetMinimum() < min) min = (*(hCOMp_ProjY.at(irig)))()->GetMinimum();
			hCOMp_ProjY.at(irig)->setStyle(Area(), Line(Color(ieff, nCOMEff)), Marker(Color(ieff, nCOMEff)));
		}
	}
	Double_t width = ((max-min) < 1e-3) ? 1e-3 : (max-min)*5e-2;
	editor.create();
	editor.cd(1, Canvas::AxisScl_t(0, 0));
	THStack * hColl = Hist::Collect(StrFmt("hCOMp_EffColl"), StrFmt("COM Time vs Efficiency"), hCOMColl);
	hColl->Draw("nostack hist");
	SetTimeAxis(hColl->GetXaxis());
	(hColl->GetHistogram())->SetLineColor(0);
	//hColl->SetMaximum(max*(1.0+width));
	//hColl->SetMinimum(min*(1.0-width));
	hColl->SetMaximum(1.01);
	hColl->SetMinimum(0.80);
	Text::Draw(Text::Txt_t(StrFmt("COM Time vs Efficiency"), kBlack, 0.03), Text::Att_t(0.50, 0.92, 22));
	editor.save();


	/*
	for (Int_t icut = 1; icut <= nCOMCut; ++icut) {
		Hist::VLIST hCOMColl;
		std::vector<Hist *>&& hCOM_Proj = Hist::Project(Hist::kProjX, hCOMp_Eff[icut]);
		Double_t max = 0.0, min = 1.0;
		for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
			Hist * hCOM = hCOM_Proj.at(itme);
			(*hCOM)()->SetNameTitle(CStrFmt("hCOMp_Eff%02d_TME%03d", icut, itme), CStrFmt("COM(%02d) Rig vs Efficiency", icut));
			hCOM->setStyle(Area(), Line(Color(itme, AXtme.nbin())), Marker(Color(itme, AXtme.nbin())));
			SetTimeAxis((*hCOM)()->GetXaxis());
			hCOMColl.push_back(hCOM);
			if ((*hCOM)()->GetMaximum() > max) max = (*hCOM)()->GetMaximum();
			if ((*hCOM)()->GetMinimum() < min) min = (*hCOM)()->GetMinimum();
		}
		Double_t width = ((max-min) < 1e-3) ? 1e-3 : (max-min)*5e-2;
		editor.create();
		editor.cd(1, Canvas::AxisScl_t(0, 0));
		THStack * hColl = Hist::Collect(StrFmt("hCOMp_Eff%02dColl", icut), StrFmt("COM(%02d) Rig vs Efficiency", icut), hCOMColl);
		hColl->Draw("nostack hist");
		SetTimeAxis(hColl->GetXaxis());
		(hColl->GetHistogram())->SetLineColor(0);
		hColl->SetMaximum(max*(1.0+width));
		hColl->SetMinimum(min*(1.0-width));
		Text::Draw(Text::Txt_t(StrFmt("COM(%02d) Rig vs Efficiency", icut), kBlack, 0.03), Text::Att_t(0.50, 0.92, 22));
		editor.save();
	}
	*/

	// Two thistograms : total error and systematic only errots  
	Hist * hStatRat = Hist::New("hStatRat", "#bar{p}/p Ratio (StatErr) Rig vs Time", AXtme, AXnr);
	Hist * hTotlRat = Hist::New("hTotlRat", "#bar{p}/p Ratio (TotlErr) Rig vs Time",  AXtme, AXnr);
	Hist * hFitNChi = Hist::New("hFitNChi", "Fit Chisquare/NDF Rig vs Time",  AXtme, AXnr);
	
	for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
		//if (itme == 16) continue; // TTCS off
		TFile * fTme = fTmeVec[itme]; 
		TGraphAsymmErrors * appstat = (TGraphAsymmErrors*)fTme->Get("apprlt_stat");
		TGraphAsymmErrors * appfull = (TGraphAsymmErrors*)fTme->Get("apprlt_totl");
		TGraphAsymmErrors * fitnchi = (TGraphAsymmErrors*)fTme->Get("appfit_nchi");
	
		for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
			Float_t cen = AXnr.center(irig, Axis::kLog);
			if (cen > 6) continue;
			Float_t val  = appstat->GetY()[irig-1];
			Float_t stat = appstat->GetEYlow()[irig-1];
			Float_t totl = appfull->GetEYlow()[irig-1];
		
			hStatRat->setContentWithError(val, stat, itme, irig);
			hTotlRat->setContentWithError(val, totl, itme, irig);
			
			hFitNChi->setContentWithError(fitnchi->GetY()[irig-1], 0.0, itme, irig);
		}
	}
	SetTimeAxis((*hStatRat)()->GetXaxis());
	SetTimeAxis((*hTotlRat)()->GetXaxis());
	SetTimeAxis((*hFitNChi)()->GetXaxis());

	std::vector<Hist *>&& hStatRatVec = Hist::Project(Hist::kProjX, hStatRat);
	std::vector<Hist *>&& hTotlRatVec = Hist::Project(Hist::kProjX, hTotlRat);
	std::vector<Hist *>&& hFitNChiVec = Hist::Project(Hist::kProjX, hFitNChi);
	for (Int_t irig = 1; irig <= AXnr.nbin(); ++irig) {
		SetTimeAxis((*hStatRatVec.at(irig))()->GetXaxis());
		//(*hStatRatVec.at(irig))()->SetMaximum(->GetMaximum());
		//(*hStatRatVec.at(irig))()->SetMinimum(->GetMinimum());
		hStatRatVec.at(irig)->setStyle(Area(), Line(kRed), Marker(kRed));
		hTotlRatVec.at(irig)->setStyle(Area(), Line(kRed), Marker(kRed));
		hFitNChiVec.at(irig)->setStyle(Area(), Line(kRed), Marker(kRed));
		
		//(*hStatRatVec.at(irig))()->SetMaximum(60.*1e-6);
		//(*hStatRatVec.at(irig))()->SetMinimum(25.*1e-6);
		//
		//(*hTotlRatVec.at(irig))()->SetMaximum(60.*1e-6);
		//(*hTotlRatVec.at(irig))()->SetMinimum(25.*1e-6);
		
		editor.create();
		editor.cd(1, Canvas::AxisScl_t(0, 0));
		hStatRatVec.at(irig)->draw("pe");
		Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(irig-1), AXnr.bins(irig)), kBlack, 0.03), Text::Att_t(0.80, 0.92, 32));
		//Text::Draw(Text::Txt_t(StrFmt("Stat"), kBlue, 0.03), Text::Att_t(0.20, 0.85, 12));
		editor.save();
		
		SetTimeAxis((*hTotlRatVec.at(irig))()->GetXaxis());
		editor.create();
		editor.cd(1, Canvas::AxisScl_t(0, 0));
		hTotlRatVec.at(irig)->draw("pe");
		Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(irig-1), AXnr.bins(irig)), kBlack, 0.03), Text::Att_t(0.80, 0.92, 32));
		//Text::Draw(Text::Txt_t(StrFmt("Totl"), kBlue, 0.03), Text::Att_t(0.20, 0.85, 12));
		editor.save();
		
		SetTimeAxis((*hFitNChiVec.at(irig))()->GetXaxis());
		editor.create();
		editor.cd(1, Canvas::AxisScl_t(0, 0));
		hFitNChiVec.at(irig)->draw("pe");
		Text::Draw(Text::Txt_t(StrFmt("Rigidity ( %.2f GV ~ %.2f GV )", AXnr.bins(irig-1), AXnr.bins(irig)), kBlack, 0.03), Text::Att_t(0.80, 0.92, 32));
		//Text::Draw(Text::Txt_t(StrFmt("Totl"), kBlue, 0.03), Text::Att_t(0.20, 0.85, 12));
		editor.save();
	}

	for (Int_t itme = 1; itme <= AXtme.nbin(); ++itme) {
		if (fTmeVec[itme] == nullptr) continue;
		fTmeVec[itme]->Close();
		fTmeVec[itme] = nullptr;
	}

	TFile * file = new TFile(CStrFmt("%s/YiFin.root", goddir.c_str()), "RECREATE");
	//Graph::Write();
	Hist::Write();
	file->Write();
	file->Close();

	editor.close();

	return 1;
}

#endif // __YiMergeAna_C__
