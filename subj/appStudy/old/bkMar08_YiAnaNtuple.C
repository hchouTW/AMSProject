#ifndef __YiAnaNtuple_C__
#define __YiAnaNtuple_C__

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.h"
#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.tcc"

#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jan22/src/ClassDef.h"
#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jan22/src/ClassDef.C"

using namespace MgntROOT;

// User defination macro
#ifdef Debug
	#undef Debug
#endif
#define Debug false

/*-------*/
/*  DST  */
/*-------*/
class DST {
	public :
		DST() : fDST(0) {}
		~DST() {}

		void initDST(TFile * file, TChain * runChain = 0, TChain * dataChain = 0) {
			if (file == 0 || !gSaveDST) return;
			resetDST();
			file->cd();

			std::cout << "\n<<  init DST Info  >>\n";

			if (runChain != 0) {
				std::string files = "";
				UInt_t runID = 0;
				UInt_t trgEV = 0;
				UInt_t selEV = 0;
				UInt_t year = 0;
				UInt_t yday = 0;
				UInt_t month = 0;
				UInt_t mday = 0;
				UInt_t hour = 0;
				UInt_t min = 0;
				TTree * DSTRun = new TTree("DSTRun", "DST Run Info");
				DSTRun->Branch("files", &files);
				DSTRun->Branch("runID", &runID, "runID/i");
				DSTRun->Branch("trgEV", &trgEV, "trgEV/i");
				DSTRun->Branch("selEV", &selEV, "selEV/i");
				if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
					DSTRun->Branch("year",  &year,  "year/i");
					DSTRun->Branch("yday",  &yday,  "yday/i");
					DSTRun->Branch("month", &month, "month/i");
					DSTRun->Branch("mday",  &mday,  "mday/i");
					DSTRun->Branch("hour",  &hour,  "hour/i");
					DSTRun->Branch("min",   &min,   "min/i");
				}
				
				RunTagInfo * runTag = new RunTagInfo;
				runChain->SetBranchAddress("runTag", &runTag);
				for (Int_t i = 0; i < runChain->GetEntries(); ++i) {
					runChain->GetEntry(i);
					files = runChain->GetCurrentFile()->GetName();
					runID = runTag->runID;
					trgEV  = runTag->numOfTrgEvent;
					selEV  = runTag->numOfSelEvent;
					if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
						MgntClock::TTime * ttime = MgntClock::ConvertFromUTimeToTTime(runID, MgntClock::ClockType::UTC);
						year  = ttime->tm_year+1900;
						yday  = ttime->tm_yday;
						month = ttime->tm_mon+1;
						mday  = ttime->tm_mday;
						hour  = ttime->tm_hour;
						min   = ttime->tm_min;
					}
					DSTRun->Fill();
				}
				delete runTag;
			}

			if (dataChain != 0 && gSaveDSTClone) fDST = dataChain->CloneTree(0);
			if (fDST == 0      && gSaveDSTTree)  fDST = new TTree("data", "DST data");
			file->cd();
			
			Axis AXrso("Rigidity Resolution", 400, -1.5, 1.5);

			// rigidity binning
			//AXnr = Axis("Rigidity [GV]",
			//	{   0.50,   0.80, // extern bins
			//	    1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
			//	    3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
			//		  8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
			//		 19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
			//		 41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
			//		 93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00, 
			//		800.00 } ); // extern bin
			//AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
		
			// merge 3 bins
			//AXnr = Axis("Rigidity [GV]",
			//	{   1.00,   1.51,   2.15,   2.97,   4.02, 
			//	    5.37,   7.09,   9.26,  12.00,  15.30, 
			//		 19.50,  24.70,  31.10,  38.90,  48.50 } );
			//AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
			
			//AXtnr = Axis("Rigidity [GV]",
			//	{   1.00,   1.51,   2.00,   3.00,   4.02, 
			//	    5.37,   7.09,   9.26,  12.00,  15.30, 
			//		 19.50,  24.70,  31.10,  38.90,  48.50 } );
			//AXtir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
			
			AXnr = Axis("Rigidity [GV]",
				{   1.00,   1.46,   2.00,   3.00,   4.12, 
				    5.00,   6.00,   7.10,   8.30,   9.62, 
				   11.04,  12.59,  14.25,  16.05,  17.98, 
				   20.04,  22.25,  24.62,  27.25,  30.21, 
				   35.36,  40.00 } );
			AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);

			// 6 month
			//AXtme = Axis("Time", 
			//	{ 1312416000, 1326412800, 1340409600, 1354406400, 1368403200, 
			//	  1382400000, 1396396800, 1410393600 ,1424390400, 1438387200, 
			//		1452384000, 1464299612, 1480377600 } );
			
			// 3 month
			AXtme = Axis("Time", 
				{ 1305417600, 1312416000, 1319414400, 1326412800, 1333411200,
				  1340409600, 1347408000, 1354406400, 1361404800, 1368403200,
					1375401600, 1382400000, 1389398400, 1396396800, 1403395200,
					1410393600, 1417392000, 1424390400, 1431388800, 1438387200,
					1445385600, 1452384000, 1459382400, 1466380800, 1473379200, 
					1480377600 } );

			//TString sbnPH = "/afs/cern.ch/user/h/hchou/public/DATABASE/physics/binning/phe_bin2.root";
			//TFile * fPHBin = TFile::Open(sbnPH);
			//TH1D  * hPHBn0 = (TH1D*)fPHBin->Get("hist2");
			//Axis AXPHnr("Rigidity [GV]", hPHBn0, Axis::kX);
      //Axis AXPHir = MgntROOT::Axis::Invert("Inverse Rigidity [1/GV]", AXPHnr);
			//Axis AXnr = AXPHnr;
			//Axis AXir = AXPHir;
			//fPHBin->Close();
			//file->cd();

			
			//---- Common Cut ----//
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Long64_t ntme = Long64_t(((AXtme.max() - AXtme.min()) / 86400.) / 7.) + 1;
				Axis AXRTItme = Axis("Time", ntme, AXtme.min(), AXtme.max());
				
				Axis AXRTIcut("Cutflow", 4, 0., 4.);
				Hist::New("hRTI_Cutflow", "", AXRTItme, AXRTIcut);
				
				Axis AXRTIlivetime("LiveTime", 100, 0.5, 1.0);
				Hist::New("hRTI_LiveTime", "", AXRTItme, AXRTIlivetime);
				
				Hist::New("hRTI_ExpT25Deg", "", AXRTItme, AXnr);
				Hist::New("hRTI_ExpT30Deg", "", AXRTItme, AXnr);
				Hist::New("hRTI_ExpT35Deg", "", AXRTItme, AXnr);
				Hist::New("hRTI_ExpT40Deg", "", AXRTItme, AXnr);

				Hist::New("hRTI_AllExpT25Deg", "", AXnr);
				Hist::New("hRTI_AllExpT30Deg", "", AXnr);
				Hist::New("hRTI_AllExpT35Deg", "", AXnr);
				Hist::New("hRTI_AllExpT40Deg", "", AXnr);
				for (Int_t it = 1; it <= AXtme.nbin() && gTimeStudy; ++it) {
					Hist::New(StrFmt("hRTI_AllExpT25Deg%03d", it), "", AXnr);
					Hist::New(StrFmt("hRTI_AllExpT30Deg%03d", it), "", AXnr);
					Hist::New(StrFmt("hRTI_AllExpT35Deg%03d", it), "", AXnr);
					Hist::New(StrFmt("hRTI_AllExpT40Deg%03d", it), "", AXnr);
				}
			}
			
			Axis AXCOMcut("Cutflow", 8, 0., 8.);
			Hist::New("hCOM_Cutflow", "", AXCOMcut);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOM_Cutflow%03d", it), "", AXCOMcut);
			}
			
			Axis AXCOMcutflow("Cutflow", 19, 0., 19.);
			Hist::New("hCOMi_Cutflow", "", AXir, AXCOMcutflow);
			Hist::New("hCOMp_Cutflow", "", AXnr, AXCOMcutflow);
			Hist::New("hCOMn_Cutflow", "", AXnr, AXCOMcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_Cutflow%03d", it), "", AXir, AXCOMcutflow);
				Hist::New(StrFmt("hCOMp_Cutflow%03d", it), "", AXnr, AXCOMcutflow);
				Hist::New(StrFmt("hCOMn_Cutflow%03d", it), "", AXnr, AXCOMcutflow);
			}
			
			Axis AXCOMcutflowPt("Cutflow", 7, 0., 7.);
			Hist::New("hCOMi_CutflowIn", "", AXir, AXCOMcutflowPt);
			Hist::New("hCOMp_CutflowIn", "", AXnr, AXCOMcutflowPt);
			Hist::New("hCOMn_CutflowIn", "", AXnr, AXCOMcutflowPt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_CutflowIn%03d", it), "", AXir, AXCOMcutflowPt);
				Hist::New(StrFmt("hCOMp_CutflowIn%03d", it), "", AXnr, AXCOMcutflowPt);
				Hist::New(StrFmt("hCOMn_CutflowIn%03d", it), "", AXnr, AXCOMcutflowPt);
			}

			Hist::New("hCOMi_CutflowL1", "", AXir, AXCOMcutflowPt);
			Hist::New("hCOMp_CutflowL1", "", AXnr, AXCOMcutflowPt);
			Hist::New("hCOMn_CutflowL1", "", AXnr, AXCOMcutflowPt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_CutflowL1%03d", it), "", AXir, AXCOMcutflowPt);
				Hist::New(StrFmt("hCOMp_CutflowL1%03d", it), "", AXnr, AXCOMcutflowPt);
				Hist::New(StrFmt("hCOMn_CutflowL1%03d", it), "", AXnr, AXCOMcutflowPt);
			}

			Hist::New("hCOMi_CutflowL9", "", AXir, AXCOMcutflowPt);
			Hist::New("hCOMp_CutflowL9", "", AXnr, AXCOMcutflowPt);
			Hist::New("hCOMn_CutflowL9", "", AXnr, AXCOMcutflowPt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_CutflowL9%03d", it), "", AXir, AXCOMcutflowPt);
				Hist::New(StrFmt("hCOMp_CutflowL9%03d", it), "", AXnr, AXCOMcutflowPt);
				Hist::New(StrFmt("hCOMn_CutflowL9%03d", it), "", AXnr, AXCOMcutflowPt);
			}

			Hist::New("hCOMi_CutflowFs", "", AXir, AXCOMcutflowPt);
			Hist::New("hCOMp_CutflowFs", "", AXnr, AXCOMcutflowPt);
			Hist::New("hCOMn_CutflowFs", "", AXnr, AXCOMcutflowPt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_CutflowFs%03d", it), "", AXir, AXCOMcutflowPt);
				Hist::New(StrFmt("hCOMp_CutflowFs%03d", it), "", AXnr, AXCOMcutflowPt);
				Hist::New(StrFmt("hCOMn_CutflowFs%03d", it), "", AXnr, AXCOMcutflowPt);
			}

			Axis AXCOMtrInQ("Chrg", 400, 0.5, 2.6);
			Hist::New("hCOMi_TrInQ", "", AXir, AXCOMtrInQ);
			Hist::New("hCOMp_TrInQ", "", AXnr, AXCOMtrInQ);
			Hist::New("hCOMn_TrInQ", "", AXnr, AXCOMtrInQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrInQ%03d", it), "", AXir, AXCOMtrInQ);
				Hist::New(StrFmt("hCOMp_TrInQ%03d", it), "", AXnr, AXCOMtrInQ);
				Hist::New(StrFmt("hCOMn_TrInQ%03d", it), "", AXnr, AXCOMtrInQ);
			}
			
			Axis AXCOMtrExQ("Chrg", 400, 0.5, 3.2);
			Hist::New("hCOMi_TrL2Q", "", AXir, AXCOMtrExQ);
			Hist::New("hCOMp_TrL2Q", "", AXnr, AXCOMtrExQ);
			Hist::New("hCOMn_TrL2Q", "", AXnr, AXCOMtrExQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrL2Q%03d", it), "", AXir, AXCOMtrExQ);
				Hist::New(StrFmt("hCOMp_TrL2Q%03d", it), "", AXnr, AXCOMtrExQ);
				Hist::New(StrFmt("hCOMn_TrL2Q%03d", it), "", AXnr, AXCOMtrExQ);
			}

			Hist::New("hCOMi_TrL1Q", "", AXir, AXCOMtrExQ);
			Hist::New("hCOMp_TrL1Q", "", AXnr, AXCOMtrExQ);
			Hist::New("hCOMn_TrL1Q", "", AXnr, AXCOMtrExQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrL1Q%03d", it), "", AXir, AXCOMtrExQ);
				Hist::New(StrFmt("hCOMp_TrL1Q%03d", it), "", AXnr, AXCOMtrExQ);
				Hist::New(StrFmt("hCOMn_TrL1Q%03d", it), "", AXnr, AXCOMtrExQ);
			}
			
			Hist::New("hCOMi_TrL9Q", "", AXir, AXCOMtrExQ);
			Hist::New("hCOMp_TrL9Q", "", AXnr, AXCOMtrExQ);
			Hist::New("hCOMn_TrL9Q", "", AXnr, AXCOMtrExQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrL9Q%03d", it), "", AXir, AXCOMtrExQ);
				Hist::New(StrFmt("hCOMp_TrL9Q%03d", it), "", AXnr, AXCOMtrExQ);
				Hist::New(StrFmt("hCOMn_TrL9Q%03d", it), "", AXnr, AXCOMtrExQ);
			}
			
			Hist::New("hCOMi_TrMxQ", "", AXir, AXCOMtrExQ);
			Hist::New("hCOMp_TrMxQ", "", AXnr, AXCOMtrExQ);
			Hist::New("hCOMn_TrMxQ", "", AXnr, AXCOMtrExQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrMxQ%03d", it), "", AXir, AXCOMtrExQ);
				Hist::New(StrFmt("hCOMp_TrMxQ%03d", it), "", AXnr, AXCOMtrExQ);
				Hist::New(StrFmt("hCOMn_TrMxQ%03d", it), "", AXnr, AXCOMtrExQ);
			}
	
			const Double_t AreaFT = 180.;
			std::vector<Double_t> radius; radius.push_back(0.);
			while (radius.at(radius.size()-1) < 52.) {
				Double_t r0 = radius.at(radius.size()-1);
				Double_t dr = 0.5 * (std::sqrt(r0*r0 + 2.*TMath::InvPi()*AreaFT) - r0);
				radius.push_back( (r0+dr) );
			}
			Axis AXCOMradius("Radius", radius);

			Hist::New("hCOMi_TrRadius", "", AXir, AXCOMradius);
			Hist::New("hCOMp_TrRadius", "", AXnr, AXCOMradius);
			Hist::New("hCOMn_TrRadius", "", AXnr, AXCOMradius);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrRadius%03d", it), "", AXir, AXCOMradius);
				Hist::New(StrFmt("hCOMp_TrRadius%03d", it), "", AXnr, AXCOMradius);
				Hist::New(StrFmt("hCOMn_TrRadius%03d", it), "", AXnr, AXCOMradius);
			}
			
			const Double_t CosFT = 0.03;
			std::vector<Double_t> aglcos; aglcos.push_back(1.);
			std::vector<Double_t> agldeg; agldeg.push_back(0.);
			while (agldeg.at(agldeg.size()-1) < 40.) {
				Double_t cos0 = aglcos.at(aglcos.size()-1);
				Double_t dcos = 0.5 * (std::sqrt(cos0*cos0 + 2.*TMath::InvPi()*CosFT) - cos0);
				Double_t cosv = (cos0 - dcos);
				Double_t degv = std::acos(cosv) * TMath::RadToDeg();
				aglcos.push_back(cosv);
				agldeg.push_back(degv);
			}
			Axis AXCOMangle("Angle", agldeg);
			
			Hist::New("hCOMi_TrAngle", "", AXir, AXCOMangle);
			Hist::New("hCOMp_TrAngle", "", AXnr, AXCOMangle);
			Hist::New("hCOMn_TrAngle", "", AXnr, AXCOMangle);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrAngle%03d", it), "", AXir, AXCOMangle);
				Hist::New(StrFmt("hCOMp_TrAngle%03d", it), "", AXnr, AXCOMangle);
				Hist::New(StrFmt("hCOMn_TrAngle%03d", it), "", AXnr, AXCOMangle);
			}

			Axis AXCOMtofQ("Chrg", 400, 0.0, 2.6);
			Hist::New("hCOMi_TofQ", "", AXir, AXCOMtofQ);
			Hist::New("hCOMp_TofQ", "", AXnr, AXCOMtofQ);
			Hist::New("hCOMn_TofQ", "", AXnr, AXCOMtofQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TofQ%03d", it), "", AXir, AXCOMtofQ);
				Hist::New(StrFmt("hCOMp_TofQ%03d", it), "", AXnr, AXCOMtofQ);
				Hist::New(StrFmt("hCOMn_TofQ%03d", it), "", AXnr, AXCOMtofQ);
			}
			
			Hist::New("hCOMi_TofQu", "", AXir, AXCOMtofQ);
			Hist::New("hCOMp_TofQu", "", AXnr, AXCOMtofQ);
			Hist::New("hCOMn_TofQu", "", AXnr, AXCOMtofQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TofQu%03d", it), "", AXir, AXCOMtofQ);
				Hist::New(StrFmt("hCOMp_TofQu%03d", it), "", AXnr, AXCOMtofQ);
				Hist::New(StrFmt("hCOMn_TofQu%03d", it), "", AXnr, AXCOMtofQ);
			}
			
			Hist::New("hCOMi_TofQl", "", AXir, AXCOMtofQ);
			Hist::New("hCOMp_TofQl", "", AXnr, AXCOMtofQ);
			Hist::New("hCOMn_TofQl", "", AXnr, AXCOMtofQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TofQl%03d", it), "", AXir, AXCOMtofQ);
				Hist::New(StrFmt("hCOMp_TofQl%03d", it), "", AXnr, AXCOMtofQ);
				Hist::New(StrFmt("hCOMn_TofQl%03d", it), "", AXnr, AXCOMtofQ);
			}
			
			Axis AXCOMdltTofQ("Delta Chrg", 400, -1.4, 1.4);
			Hist::New("hCOMi_TofQul", "", AXir, AXCOMdltTofQ);
			Hist::New("hCOMp_TofQul", "", AXnr, AXCOMdltTofQ);
			Hist::New("hCOMn_TofQul", "", AXnr, AXCOMdltTofQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TofQul%03d", it), "", AXir, AXCOMdltTofQ);
				Hist::New(StrFmt("hCOMp_TofQul%03d", it), "", AXnr, AXCOMdltTofQ);
				Hist::New(StrFmt("hCOMn_TofQul%03d", it), "", AXnr, AXCOMdltTofQ);
			}
			
			Axis AXCOMtofM("Mass Estimator", 400, -1.3, 1.3);
			Hist::New("hCOMi_TofM", "", AXir, AXCOMtofM);
			Hist::New("hCOMp_TofM", "", AXnr, AXCOMtofM);
			Hist::New("hCOMn_TofM", "", AXnr, AXCOMtofM);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TofM%03d", it), "", AXir, AXCOMtofM);
				Hist::New(StrFmt("hCOMp_TofM%03d", it), "", AXnr, AXCOMtofM);
				Hist::New(StrFmt("hCOMn_TofM%03d", it), "", AXnr, AXCOMtofM);
			}
			
			Axis AXCOMnacc("NAcc", 6, 0., 6.);
			Hist::New("hCOMi_NAcc", "", AXir, AXCOMnacc);
			Hist::New("hCOMp_NAcc", "", AXnr, AXCOMnacc);
			Hist::New("hCOMn_NAcc", "", AXnr, AXCOMnacc);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_NAcc%03d", it), "", AXir, AXCOMnacc);
				Hist::New(StrFmt("hCOMp_NAcc%03d", it), "", AXnr, AXCOMnacc);
				Hist::New(StrFmt("hCOMn_NAcc%03d", it), "", AXnr, AXCOMnacc);
			}
			
			Axis AXCOMntrtdVtx("NTrkTrdVtx", 25, 0., 25.);
			Hist::New("hCOMi_NTrkTrdVtx", "", AXir, AXCOMntrtdVtx);
			Hist::New("hCOMp_NTrkTrdVtx", "", AXnr, AXCOMntrtdVtx);
			Hist::New("hCOMn_NTrkTrdVtx", "", AXnr, AXCOMntrtdVtx);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_NTrkTrdVtx%03d", it), "", AXir, AXCOMntrtdVtx);
				Hist::New(StrFmt("hCOMp_NTrkTrdVtx%03d", it), "", AXnr, AXCOMntrtdVtx);
				Hist::New(StrFmt("hCOMn_NTrkTrdVtx%03d", it), "", AXnr, AXCOMntrtdVtx);
			}
			
			Axis AXCOMtrdCls("NTrdCls", 300, 0., 300.);
			Hist::New("hCOMi_NTrdCls", "", AXir, AXCOMtrdCls);
			Hist::New("hCOMp_NTrdCls", "", AXnr, AXCOMtrdCls);
			Hist::New("hCOMn_NTrdCls", "", AXnr, AXCOMtrdCls);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_NTrdCls%03d", it), "", AXir, AXCOMtrdCls);
				Hist::New(StrFmt("hCOMp_NTrdCls%03d", it), "", AXnr, AXCOMtrdCls);
				Hist::New(StrFmt("hCOMn_NTrdCls%03d", it), "", AXnr, AXCOMtrdCls);
			}
			
			Axis AXCOMtrdl("Trdl", 200, 0.2, 1.6);
			Hist::New("hCOMi_Trdl", "", AXir, AXCOMtrdl);
			Hist::New("hCOMp_Trdl", "", AXnr, AXCOMtrdl);
			Hist::New("hCOMn_Trdl", "", AXnr, AXCOMtrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_Trdl%03d", it), "", AXir, AXCOMtrdl);
				Hist::New(StrFmt("hCOMp_Trdl%03d", it), "", AXnr, AXCOMtrdl);
				Hist::New(StrFmt("hCOMn_Trdl%03d", it), "", AXnr, AXCOMtrdl);
			}
			
			Axis AXCOMtrdQ("TrdQ", 200, 0.0, 2.2);
			Hist::New("hCOMi_TrdQ", "", AXir, AXCOMtrdQ);
			Hist::New("hCOMp_TrdQ", "", AXnr, AXCOMtrdQ);
			Hist::New("hCOMn_TrdQ", "", AXnr, AXCOMtrdQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrdQ%03d", it), "", AXir, AXCOMtrdQ);
				Hist::New(StrFmt("hCOMp_TrdQ%03d", it), "", AXnr, AXCOMtrdQ);
				Hist::New(StrFmt("hCOMn_TrdQ%03d", it), "", AXnr, AXCOMtrdQ);
			}
			
			Axis AXCOMbdt("EcalBDT", 400, -1.0, 1.0);
			Hist::New("hCOMi_EcalBDT", "", AXir, AXCOMbdt);
			Hist::New("hCOMp_EcalBDT", "", AXnr, AXCOMbdt);
			Hist::New("hCOMn_EcalBDT", "", AXnr, AXCOMbdt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_EcalBDT%03d", it), "", AXir, AXCOMbdt);
				Hist::New(StrFmt("hCOMp_EcalBDT%03d", it), "", AXnr, AXCOMbdt);
				Hist::New(StrFmt("hCOMn_EcalBDT%03d", it), "", AXnr, AXCOMbdt);
			}
			
			Axis AXCOMaglbPr("AglbPr", 400, -0.04, 0.02);
			Hist::New("hCOMi_AglbPr", "", AXir, AXCOMaglbPr);
			Hist::New("hCOMp_AglbPr", "", AXnr, AXCOMaglbPr);
			Hist::New("hCOMn_AglbPr", "", AXnr, AXCOMaglbPr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_AglbPr%03d", it), "", AXir, AXCOMaglbPr);
				Hist::New(StrFmt("hCOMp_AglbPr%03d", it), "", AXnr, AXCOMaglbPr);
				Hist::New(StrFmt("hCOMn_AglbPr%03d", it), "", AXnr, AXCOMaglbPr);
			}
			
			Axis AXCOMnafbPr("NafbPr", 400, -0.06, 0.03);
			Hist::New("hCOMi_NafbPr", "", AXir, AXCOMnafbPr);
			Hist::New("hCOMp_NafbPr", "", AXnr, AXCOMnafbPr);
			Hist::New("hCOMn_NafbPr", "", AXnr, AXCOMnafbPr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_NafbPr%03d", it), "", AXir, AXCOMnafbPr);
				Hist::New(StrFmt("hCOMp_NafbPr%03d", it), "", AXnr, AXCOMnafbPr);
				Hist::New(StrFmt("hCOMn_NafbPr%03d", it), "", AXnr, AXCOMnafbPr);
			}
			
			Axis AXCOMaglbEl("AglbEl", 400, -0.02, 0.04);
			Hist::New("hCOMi_AglbEl", "", AXir, AXCOMaglbEl);
			Hist::New("hCOMp_AglbEl", "", AXnr, AXCOMaglbEl);
			Hist::New("hCOMn_AglbEl", "", AXnr, AXCOMaglbEl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_AglbEl%03d", it), "", AXir, AXCOMaglbEl);
				Hist::New(StrFmt("hCOMp_AglbEl%03d", it), "", AXnr, AXCOMaglbEl);
				Hist::New(StrFmt("hCOMn_AglbEl%03d", it), "", AXnr, AXCOMaglbEl);
			}
			
			Axis AXCOMnafbEl("NafbEl", 400, -0.03, 0.06);
			Hist::New("hCOMi_NafbEl", "", AXir, AXCOMnafbEl);
			Hist::New("hCOMp_NafbEl", "", AXnr, AXCOMnafbEl);
			Hist::New("hCOMn_NafbEl", "", AXnr, AXCOMnafbEl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_NafbEl%03d", it), "", AXir, AXCOMnafbEl);
				Hist::New(StrFmt("hCOMp_NafbEl%03d", it), "", AXnr, AXCOMnafbEl);
				Hist::New(StrFmt("hCOMn_NafbEl%03d", it), "", AXnr, AXCOMnafbEl);
			}
		
			Axis AXCOMaglMPr("AglMPr", 200, -0.03, 0.05);
			Hist::New("hCOMi_AglMPr", "", AXir, AXCOMaglMPr);
			Hist::New("hCOMp_AglMPr", "", AXnr, AXCOMaglMPr);
			Hist::New("hCOMn_AglMPr", "", AXnr, AXCOMaglMPr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_AglMPr%03d", it), "", AXir, AXCOMaglMPr);
				Hist::New(StrFmt("hCOMp_AglMPr%03d", it), "", AXnr, AXCOMaglMPr);
				Hist::New(StrFmt("hCOMn_AglMPr%03d", it), "", AXnr, AXCOMaglMPr);
			}
			
			Axis AXCOMnafMPr("NafMPr", 200, -0.06, 0.10);
			Hist::New("hCOMi_NafMPr", "", AXir, AXCOMnafMPr);
			Hist::New("hCOMp_NafMPr", "", AXnr, AXCOMnafMPr);
			Hist::New("hCOMn_NafMPr", "", AXnr, AXCOMnafMPr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_NafMPr%03d", it), "", AXir, AXCOMnafMPr);
				Hist::New(StrFmt("hCOMp_NafMPr%03d", it), "", AXnr, AXCOMnafMPr);
				Hist::New(StrFmt("hCOMn_NafMPr%03d", it), "", AXnr, AXCOMnafMPr);
			}
			
			Axis AXCOMaglMEl("AglMEl", 200, -0.03, 0.05);
			Hist::New("hCOMi_AglMEl", "", AXir, AXCOMaglMEl);
			Hist::New("hCOMp_AglMEl", "", AXnr, AXCOMaglMEl);
			Hist::New("hCOMn_AglMEl", "", AXnr, AXCOMaglMEl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_AglMEl%03d", it), "", AXir, AXCOMaglMEl);
				Hist::New(StrFmt("hCOMp_AglMEl%03d", it), "", AXnr, AXCOMaglMEl);
				Hist::New(StrFmt("hCOMn_AglMEl%03d", it), "", AXnr, AXCOMaglMEl);
			}
			
			Axis AXCOMnafMEl("NafMEl", 200, -0.06, 0.10);
			Hist::New("hCOMi_NafMEl", "", AXir, AXCOMnafMEl);
			Hist::New("hCOMp_NafMEl", "", AXnr, AXCOMnafMEl);
			Hist::New("hCOMn_NafMEl", "", AXnr, AXCOMnafMEl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_NafMEl%03d", it), "", AXir, AXCOMnafMEl);
				Hist::New(StrFmt("hCOMp_NafMEl%03d", it), "", AXnr, AXCOMnafMEl);
				Hist::New(StrFmt("hCOMn_NafMEl%03d", it), "", AXnr, AXCOMnafMEl);
			}
			
			Axis AXCOMlchi("lchi", 400, -7., 7.);
			
			Hist::New("hCOMi_LchixIn", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchixIn", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchixIn", "", AXnr, AXCOMlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LchixIn%03d", it), "", AXir, AXCOMlchi);
				Hist::New(StrFmt("hCOMp_LchixIn%03d", it), "", AXnr, AXCOMlchi);
				Hist::New(StrFmt("hCOMn_LchixIn%03d", it), "", AXnr, AXCOMlchi);
			}
			
			Hist::New("hCOMi_LchixL1", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchixL1", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchixL1", "", AXnr, AXCOMlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LchixL1%03d", it), "", AXir, AXCOMlchi);
				Hist::New(StrFmt("hCOMp_LchixL1%03d", it), "", AXnr, AXCOMlchi);
				Hist::New(StrFmt("hCOMn_LchixL1%03d", it), "", AXnr, AXCOMlchi);
			}
			
			Hist::New("hCOMi_LchixL9", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchixL9", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchixL9", "", AXnr, AXCOMlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LchixL9%03d", it), "", AXir, AXCOMlchi);
				Hist::New(StrFmt("hCOMp_LchixL9%03d", it), "", AXnr, AXCOMlchi);
				Hist::New(StrFmt("hCOMn_LchixL9%03d", it), "", AXnr, AXCOMlchi);
			}
			
			Hist::New("hCOMi_LchixFs", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchixFs", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchixFs", "", AXnr, AXCOMlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LchixFs%03d", it), "", AXir, AXCOMlchi);
				Hist::New(StrFmt("hCOMp_LchixFs%03d", it), "", AXnr, AXCOMlchi);
				Hist::New(StrFmt("hCOMn_LchixFs%03d", it), "", AXnr, AXCOMlchi);
			}
		
			Hist::New("hCOMi_LchiyIn", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchiyIn", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyIn", "", AXnr, AXCOMlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LchiyIn%03d", it), "", AXir, AXCOMlchi);
				Hist::New(StrFmt("hCOMp_LchiyIn%03d", it), "", AXnr, AXCOMlchi);
				Hist::New(StrFmt("hCOMn_LchiyIn%03d", it), "", AXnr, AXCOMlchi);
			}

			Hist::New("hCOMi_LchiyL1", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchiyL1", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyL1", "", AXnr, AXCOMlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LchiyL1%03d", it), "", AXir, AXCOMlchi);
				Hist::New(StrFmt("hCOMp_LchiyL1%03d", it), "", AXnr, AXCOMlchi);
				Hist::New(StrFmt("hCOMn_LchiyL1%03d", it), "", AXnr, AXCOMlchi);
			}

			Hist::New("hCOMi_LchiyL9", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchiyL9", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyL9", "", AXnr, AXCOMlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LchiyL9%03d", it), "", AXir, AXCOMlchi);
				Hist::New(StrFmt("hCOMp_LchiyL9%03d", it), "", AXnr, AXCOMlchi);
				Hist::New(StrFmt("hCOMn_LchiyL9%03d", it), "", AXnr, AXCOMlchi);
			}

			Hist::New("hCOMi_LchiyFs", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchiyFs", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyFs", "", AXnr, AXCOMlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LchiyFs%03d", it), "", AXir, AXCOMlchi);
				Hist::New(StrFmt("hCOMp_LchiyFs%03d", it), "", AXnr, AXCOMlchi);
				Hist::New(StrFmt("hCOMn_LchiyFs%03d", it), "", AXnr, AXCOMlchi);
			}

			Axis AXCOMlasym("Lasym", 400, -6., 4.);
			
			Hist::New("hCOMi_LasymUL", "", AXir, AXCOMlasym);
			Hist::New("hCOMp_LasymUL", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_LasymUL", "", AXnr, AXCOMlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_LasymUL%03d", it), "", AXir, AXCOMlasym);
				Hist::New(StrFmt("hCOMp_LasymUL%03d", it), "", AXnr, AXCOMlasym);
				Hist::New(StrFmt("hCOMn_LasymUL%03d", it), "", AXnr, AXCOMlasym);
			}
			
			Hist::New("hCOMi_Lasym1I", "", AXir, AXCOMlasym);
			Hist::New("hCOMp_Lasym1I", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym1I", "", AXnr, AXCOMlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_Lasym1I%03d", it), "", AXir, AXCOMlasym);
				Hist::New(StrFmt("hCOMp_Lasym1I%03d", it), "", AXnr, AXCOMlasym);
				Hist::New(StrFmt("hCOMn_Lasym1I%03d", it), "", AXnr, AXCOMlasym);
			}
			
			Hist::New("hCOMi_Lasym9I", "", AXir, AXCOMlasym);
			Hist::New("hCOMp_Lasym9I", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym9I", "", AXnr, AXCOMlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_Lasym9I%03d", it), "", AXir, AXCOMlasym);
				Hist::New(StrFmt("hCOMp_Lasym9I%03d", it), "", AXnr, AXCOMlasym);
				Hist::New(StrFmt("hCOMn_Lasym9I%03d", it), "", AXnr, AXCOMlasym);
			}
			
			Hist::New("hCOMi_Lasym91", "", AXir, AXCOMlasym);
			Hist::New("hCOMp_Lasym91", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym91", "", AXnr, AXCOMlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_Lasym91%03d", it), "", AXir, AXCOMlasym);
				Hist::New(StrFmt("hCOMp_Lasym91%03d", it), "", AXnr, AXCOMlasym);
				Hist::New(StrFmt("hCOMn_Lasym91%03d", it), "", AXnr, AXCOMlasym);
			}

			Axis AXCOMccest("CCest", 200, -5.0, 3.0);
			
			Hist::New("hCOMi_CCestIn", "", AXir, AXCOMccest);
			Hist::New("hCOMp_CCestIn", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestIn", "", AXnr, AXCOMccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_CCestIn%03d", it), "", AXir, AXCOMccest);
				Hist::New(StrFmt("hCOMp_CCestIn%03d", it), "", AXnr, AXCOMccest);
				Hist::New(StrFmt("hCOMn_CCestIn%03d", it), "", AXnr, AXCOMccest);
			}
			
			Hist::New("hCOMi_CCestL1", "", AXir, AXCOMccest);
			Hist::New("hCOMp_CCestL1", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestL1", "", AXnr, AXCOMccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_CCestL1%03d", it), "", AXir, AXCOMccest);
				Hist::New(StrFmt("hCOMp_CCestL1%03d", it), "", AXnr, AXCOMccest);
				Hist::New(StrFmt("hCOMn_CCestL1%03d", it), "", AXnr, AXCOMccest);
			}
			
			Hist::New("hCOMi_CCestL9", "", AXir, AXCOMccest);
			Hist::New("hCOMp_CCestL9", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestL9", "", AXnr, AXCOMccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_CCestL9%03d", it), "", AXir, AXCOMccest);
				Hist::New(StrFmt("hCOMp_CCestL9%03d", it), "", AXnr, AXCOMccest);
				Hist::New(StrFmt("hCOMn_CCestL9%03d", it), "", AXnr, AXCOMccest);
			}
			
			Hist::New("hCOMi_CCestFs", "", AXir, AXCOMccest);
			Hist::New("hCOMp_CCestFs", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestFs", "", AXnr, AXCOMccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_CCestFs%03d", it), "", AXir, AXCOMccest);
				Hist::New(StrFmt("hCOMp_CCestFs%03d", it), "", AXnr, AXCOMccest);
				Hist::New(StrFmt("hCOMn_CCestFs%03d", it), "", AXnr, AXCOMccest);
			}
		
			Hist::New("hCOMi_EvtIn", "", AXir);
			Hist::New("hCOMp_EvtIn", "", AXnr);
			Hist::New("hCOMn_EvtIn", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_EvtIn%03d", it), "", AXir);
				Hist::New(StrFmt("hCOMp_EvtIn%03d", it), "", AXnr);
				Hist::New(StrFmt("hCOMn_EvtIn%03d", it), "", AXnr);
			}
			
			Hist::New("hCOMi_EvtL1", "", AXir);
			Hist::New("hCOMp_EvtL1", "", AXnr);
			Hist::New("hCOMn_EvtL1", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_EvtL1%03d", it), "", AXir);
				Hist::New(StrFmt("hCOMp_EvtL1%03d", it), "", AXnr);
				Hist::New(StrFmt("hCOMn_EvtL1%03d", it), "", AXnr);
			}
			
			Hist::New("hCOMi_EvtL9", "", AXir);
			Hist::New("hCOMp_EvtL9", "", AXnr);
			Hist::New("hCOMn_EvtL9", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_EvtL9%03d", it), "", AXir);
				Hist::New(StrFmt("hCOMp_EvtL9%03d", it), "", AXnr);
				Hist::New(StrFmt("hCOMn_EvtL9%03d", it), "", AXnr);
			}
			
			Hist::New("hCOMi_EvtFs", "", AXir);
			Hist::New("hCOMp_EvtFs", "", AXnr);
			Hist::New("hCOMn_EvtFs", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_EvtFs%03d", it), "", AXir);
				Hist::New(StrFmt("hCOMp_EvtFs%03d", it), "", AXnr);
				Hist::New(StrFmt("hCOMn_EvtFs%03d", it), "", AXnr);
			}

			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hCOMi_MCTrRsoIn", "", AXir, AXir);
				Hist::New("hCOMp_MCTrRsoIn", "", AXnr, AXrso);
				Hist::New("hCOMn_MCTrRsoIn", "", AXnr, AXrso);
				
				Hist::New("hCOMi_MCTrRsoL1", "", AXir, AXir);
				Hist::New("hCOMp_MCTrRsoL1", "", AXnr, AXrso);
				Hist::New("hCOMn_MCTrRsoL1", "", AXnr, AXrso);
				
				Hist::New("hCOMi_MCTrRsoL9", "", AXir, AXir);
				Hist::New("hCOMp_MCTrRsoL9", "", AXnr, AXrso);
				Hist::New("hCOMn_MCTrRsoL9", "", AXnr, AXrso);
				
				Hist::New("hCOMi_MCTrRsoFs", "", AXir, AXir);
				Hist::New("hCOMp_MCTrRsoFs", "", AXnr, AXrso);
				Hist::New("hCOMn_MCTrRsoFs", "", AXnr, AXrso);
				
				Hist::New("hCOMi_MCEvtIn", "", AXir);
				Hist::New("hCOMp_MCEvtIn", "", AXnr);
				Hist::New("hCOMn_MCEvtIn", "", AXnr);
				
				Hist::New("hCOMi_MCEvtL1", "", AXir);
				Hist::New("hCOMp_MCEvtL1", "", AXnr);
				Hist::New("hCOMn_MCEvtL1", "", AXnr);
				
				Hist::New("hCOMi_MCEvtL9", "", AXir);
				Hist::New("hCOMp_MCEvtL9", "", AXnr);
				Hist::New("hCOMn_MCEvtL9", "", AXnr);
				
				Hist::New("hCOMi_MCEvtFs", "", AXir);
				Hist::New("hCOMp_MCEvtFs", "", AXnr);
				Hist::New("hCOMn_MCEvtFs", "", AXnr);
			}

			//----  Low Energy  ----//
			Axis AXLcutflow("Cutflow", 7, 0., 7.);
			Hist::New("hLi_Cutflow", "", AXir, AXLcutflow);
			Hist::New("hLp_Cutflow", "", AXnr, AXLcutflow);
			Hist::New("hLn_Cutflow", "", AXnr, AXLcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_Cutflow%03d", it), "", AXir, AXLcutflow);
				Hist::New(StrFmt("hLp_Cutflow%03d", it), "", AXnr, AXLcutflow);
				Hist::New(StrFmt("hLn_Cutflow%03d", it), "", AXnr, AXLcutflow);
			}
			
			Axis AXLccest("CCest", 100, -3.5, 2.5);
			Hist::New("hLi_CCest", "", AXir, AXLccest);
			Hist::New("hLp_CCest", "", AXnr, AXLccest);
			Hist::New("hLn_CCest", "", AXnr, AXLccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_CCest%03d", it), "", AXir, AXLccest);
				Hist::New(StrFmt("hLp_CCest%03d", it), "", AXnr, AXLccest);
				Hist::New(StrFmt("hLn_CCest%03d", it), "", AXnr, AXLccest);
			}
			
			Axis AXLprph("Prph", 400, 0., 10.);
			Hist::New("hLi_AglprphHV", "", AXir, AXLprph);
			Hist::New("hLp_AglprphHV", "", AXnr, AXLprph);
			Hist::New("hLn_AglprphHV", "", AXnr, AXLprph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_AglprphHV%03d", it), "", AXir, AXLprph);
				Hist::New(StrFmt("hLp_AglprphHV%03d", it), "", AXnr, AXLprph);
				Hist::New(StrFmt("hLn_AglprphHV%03d", it), "", AXnr, AXLprph);
			}
			
			Hist::New("hLi_AglprphNO", "", AXir, AXLprph);
			Hist::New("hLp_AglprphNO", "", AXnr, AXLprph);
			Hist::New("hLn_AglprphNO", "", AXnr, AXLprph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_AglprphNO%03d", it), "", AXir, AXLprph);
				Hist::New(StrFmt("hLp_AglprphNO%03d", it), "", AXnr, AXLprph);
				Hist::New(StrFmt("hLn_AglprphNO%03d", it), "", AXnr, AXLprph);
			}
			
			Hist::New("hLi_NafprphHV", "", AXir, AXLprph);
			Hist::New("hLp_NafprphHV", "", AXnr, AXLprph);
			Hist::New("hLn_NafprphHV", "", AXnr, AXLprph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_NafprphHV%03d", it), "", AXir, AXLprph);
				Hist::New(StrFmt("hLp_NafprphHV%03d", it), "", AXnr, AXLprph);
				Hist::New(StrFmt("hLn_NafprphHV%03d", it), "", AXnr, AXLprph);
			}
			
			Hist::New("hLi_NafprphNO", "", AXir, AXLprph);
			Hist::New("hLp_NafprphNO", "", AXnr, AXLprph);
			Hist::New("hLn_NafprphNO", "", AXnr, AXLprph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_NafprphNO%03d", it), "", AXir, AXLprph);
				Hist::New(StrFmt("hLp_NafprphNO%03d", it), "", AXnr, AXLprph);
				Hist::New(StrFmt("hLn_NafprphNO%03d", it), "", AXnr, AXLprph);
			}
			
			Axis AXLelph("Elph", 400, 0., 10.);
			Hist::New("hLi_AglelphHV", "", AXir, AXLelph);
			Hist::New("hLp_AglelphHV", "", AXnr, AXLelph);
			Hist::New("hLn_AglelphHV", "", AXnr, AXLelph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_AglelphHV%03d", it), "", AXir, AXLelph);
				Hist::New(StrFmt("hLp_AglelphHV%03d", it), "", AXnr, AXLelph);
				Hist::New(StrFmt("hLn_AglelphHV%03d", it), "", AXnr, AXLelph);
			}
			
			Hist::New("hLi_AglelphNO", "", AXir, AXLelph);
			Hist::New("hLp_AglelphNO", "", AXnr, AXLelph);
			Hist::New("hLn_AglelphNO", "", AXnr, AXLelph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_AglelphNO%03d", it), "", AXir, AXLelph);
				Hist::New(StrFmt("hLp_AglelphNO%03d", it), "", AXnr, AXLelph);
				Hist::New(StrFmt("hLn_AglelphNO%03d", it), "", AXnr, AXLelph);
			}
			
			Hist::New("hLi_NafelphHV", "", AXir, AXLelph);
			Hist::New("hLp_NafelphHV", "", AXnr, AXLelph);
			Hist::New("hLn_NafelphHV", "", AXnr, AXLelph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_NafelphHV%03d", it), "", AXir, AXLelph);
				Hist::New(StrFmt("hLp_NafelphHV%03d", it), "", AXnr, AXLelph);
				Hist::New(StrFmt("hLn_NafelphHV%03d", it), "", AXnr, AXLelph);
			}
			
			Hist::New("hLi_NafelphNO", "", AXir, AXLelph);
			Hist::New("hLp_NafelphNO", "", AXnr, AXLelph);
			Hist::New("hLn_NafelphNO", "", AXnr, AXLelph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_NafelphNO%03d", it), "", AXir, AXLelph);
				Hist::New(StrFmt("hLp_NafelphNO%03d", it), "", AXnr, AXLelph);
				Hist::New(StrFmt("hLn_NafelphNO%03d", it), "", AXnr, AXLelph);
			}
			
			Axis AXLph("Ph", 50, 0., 50.);
			Hist::New("hLi_AglphHV", "", AXir, AXLph);
			Hist::New("hLp_AglphHV", "", AXnr, AXLph);
			Hist::New("hLn_AglphHV", "", AXnr, AXLph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_AglphHV%03d", it), "", AXir, AXLph);
				Hist::New(StrFmt("hLp_AglphHV%03d", it), "", AXnr, AXLph);
				Hist::New(StrFmt("hLn_AglphHV%03d", it), "", AXnr, AXLph);
			}
			
			Hist::New("hLi_AglphNO", "", AXir, AXLph);
			Hist::New("hLp_AglphNO", "", AXnr, AXLph);
			Hist::New("hLn_AglphNO", "", AXnr, AXLph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_AglphNO%03d", it), "", AXir, AXLph);
				Hist::New(StrFmt("hLp_AglphNO%03d", it), "", AXnr, AXLph);
				Hist::New(StrFmt("hLn_AglphNO%03d", it), "", AXnr, AXLph);
			}
			
			Hist::New("hLi_NafphHV", "", AXir, AXLph);
			Hist::New("hLp_NafphHV", "", AXnr, AXLph);
			Hist::New("hLn_NafphHV", "", AXnr, AXLph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_NafphHV%03d", it), "", AXir, AXLph);
				Hist::New(StrFmt("hLp_NafphHV%03d", it), "", AXnr, AXLph);
				Hist::New(StrFmt("hLn_NafphHV%03d", it), "", AXnr, AXLph);
			}
			
			Hist::New("hLi_NafphNO", "", AXir, AXLph);
			Hist::New("hLp_NafphNO", "", AXnr, AXLph);
			Hist::New("hLn_NafphNO", "", AXnr, AXLph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_NafphNO%03d", it), "", AXir, AXLph);
				Hist::New(StrFmt("hLp_NafphNO%03d", it), "", AXnr, AXLph);
				Hist::New(StrFmt("hLn_NafphNO%03d", it), "", AXnr, AXLph);
			}
			
			Axis AXLtofb("Tofb", 50, -0.16, 0.4);
			Hist::New("hLs_Tofb", "", AXnr, AXLtofb);
			Hist::New("hLb_Tofb", "", AXnr, AXLtofb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLs_Tofb%03d", it), "", AXnr, AXLtofb);
				Hist::New(StrFmt("hLb_Tofb%03d", it), "", AXnr, AXLtofb);
			}

			Hist::New("hLi_Tofb", "", AXir, AXLtofb);
			Hist::New("hLp_Tofb", "", AXnr, AXLtofb);
			Hist::New("hLn_Tofb", "", AXnr, AXLtofb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_Tofb%03d", it), "", AXir, AXLtofb);
				Hist::New(StrFmt("hLp_Tofb%03d", it), "", AXnr, AXLtofb);
				Hist::New(StrFmt("hLn_Tofb%03d", it), "", AXnr, AXLtofb);
			}
			
			Axis AXLtofM("Mass Estimator", 50, -1.1, 0.35);
			Hist::New("hLs_TofM", "", AXnr, AXLtofM);
			Hist::New("hLb_TofM", "", AXnr, AXLtofM);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLs_TofM%03d", it), "", AXnr, AXLtofM);
				Hist::New(StrFmt("hLb_TofM%03d", it), "", AXnr, AXLtofM);
			}
			
			Hist::New("hLi_TofM", "", AXir, AXLtofM);
			Hist::New("hLp_TofM", "", AXnr, AXLtofM);
			Hist::New("hLn_TofM", "", AXnr, AXLtofM);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_TofM%03d", it), "", AXir, AXLtofM);
				Hist::New(StrFmt("hLp_TofM%03d", it), "", AXnr, AXLtofM);
				Hist::New(StrFmt("hLn_TofM%03d", it), "", AXnr, AXLtofM);
			}

			Hist::New("hLi_Evt", "", AXir);
			Hist::New("hLp_Evt", "", AXnr);
			Hist::New("hLn_Evt", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_Evt%03d", it), "", AXir);
				Hist::New(StrFmt("hLp_Evt%03d", it), "", AXnr);
				Hist::New(StrFmt("hLn_Evt%03d", it), "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hLi_MCEvt", "", AXir);
				Hist::New("hLp_MCEvt", "", AXnr);
				Hist::New("hLn_MCEvt", "", AXnr);
			}

			//----  Intermedia Energy  ----//
			Axis AXIcutflow("Cutflow", 4, 0., 4.);
			Hist::New("hIi_Cutflow", "", AXir, AXIcutflow);
			Hist::New("hIp_Cutflow", "", AXnr, AXIcutflow);
			Hist::New("hIn_Cutflow", "", AXnr, AXIcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Cutflow%03d", it), "", AXir, AXIcutflow);
				Hist::New(StrFmt("hIp_Cutflow%03d", it), "", AXnr, AXIcutflow);
				Hist::New(StrFmt("hIn_Cutflow%03d", it), "", AXnr, AXIcutflow);
			}
			
			Axis AXIccest("CCest", 100, -3.5, 2.5);
			Hist::New("hIi_CCest", "", AXir, AXIccest);
			Hist::New("hIp_CCest", "", AXnr, AXIccest);
			Hist::New("hIn_CCest", "", AXnr, AXIccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_CCest%03d", it), "", AXir, AXIccest);
				Hist::New(StrFmt("hIp_CCest%03d", it), "", AXnr, AXIccest);
				Hist::New(StrFmt("hIn_CCest%03d", it), "", AXnr, AXIccest);
			}
			
			Axis AXItrdl("Trdl", 50, 0.2, 1.6);
			Hist::New("hIs_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIb_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_Trdl", "", AXir, AXItrdl);
			Hist::New("hIp_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIn_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Trdl%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_TrdlI", "", AXnr, AXItrdl);
			Hist::New("hIb_TrdlI", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_TrdlI%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_TrdlI%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_TrdlI", "", AXir, AXItrdl);
			Hist::New("hIp_TrdlI", "", AXnr, AXItrdl);
			Hist::New("hIn_TrdlI", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_TrdlI%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_TrdlI%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_TrdlI%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_TrdlL", "", AXnr, AXItrdl);
			Hist::New("hIb_TrdlL", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_TrdlL%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_TrdlL%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_TrdlL", "", AXir, AXItrdl);
			Hist::New("hIp_TrdlL", "", AXnr, AXItrdl);
			Hist::New("hIn_TrdlL", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_TrdlL%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_TrdlL%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_TrdlL%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_TrdlH", "", AXnr, AXItrdl);
			Hist::New("hIb_TrdlH", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_TrdlH%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_TrdlH%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_TrdlH", "", AXir, AXItrdl);
			Hist::New("hIp_TrdlH", "", AXnr, AXItrdl);
			Hist::New("hIn_TrdlH", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_TrdlH%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_TrdlH%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_TrdlH%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_Evt", "", AXir);
			Hist::New("hIp_Evt", "", AXnr);
			Hist::New("hIn_Evt", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Evt%03d", it), "", AXir);
				Hist::New(StrFmt("hIp_Evt%03d", it), "", AXnr);
				Hist::New(StrFmt("hIn_Evt%03d", it), "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hIi_MCEvt", "", AXir);
				Hist::New("hIp_MCEvt", "", AXnr);
				Hist::New("hIn_MCEvt", "", AXnr);
			}
			
			//----  High Energy  ----//
			Axis AXHcutflow("Cutflow", 4, 0., 4.);
			Hist::New("hHi_Cutflow", "", AXir, AXHcutflow);
			Hist::New("hHp_Cutflow", "", AXnr, AXHcutflow);
			Hist::New("hHn_Cutflow", "", AXnr, AXHcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_Cutflow%03d", it), "", AXir, AXHcutflow);
				Hist::New(StrFmt("hHp_Cutflow%03d", it), "", AXnr, AXHcutflow);
				Hist::New(StrFmt("hHn_Cutflow%03d", it), "", AXnr, AXHcutflow);
			}

			Axis AXHccest("CCest", 50, -3.5, 2.5);
			Hist::New("hHi_CCestL1", "", AXir, AXHccest);
			Hist::New("hHp_CCestL1", "", AXnr, AXHccest);
			Hist::New("hHn_CCestL1", "", AXnr, AXHccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_CCestL1%03d", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_CCestL1%03d", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_CCestL1%03d", it), "", AXnr, AXHccest);
			}
			
			Hist::New("hHi_CCestL9", "", AXir, AXHccest);
			Hist::New("hHp_CCestL9", "", AXnr, AXHccest);
			Hist::New("hHn_CCestL9", "", AXnr, AXHccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_CCestL9%03d", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_CCestL9%03d", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_CCestL9%03d", it), "", AXnr, AXHccest);
			}
			
			Hist::New("hHi_CCestFs", "", AXir, AXHccest);
			Hist::New("hHp_CCestFs", "", AXnr, AXHccest);
			Hist::New("hHn_CCestFs", "", AXnr, AXHccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_CCestFs%03d", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_CCestFs%03d", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_CCestFs%03d", it), "", AXnr, AXHccest);
			}

			for (Int_t it = 0; it <= 20; ++it) {
				Hist::New(StrFmt("hHi_CCest%03dL1", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_CCest%03dL1", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_CCest%03dL1", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHi_CCest%03dL9", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_CCest%03dL9", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_CCest%03dL9", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHi_CCest%03dFs", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_CCest%03dFs", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_CCest%03dFs", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHi_Evt%03dL1", it), "", AXir);
				Hist::New(StrFmt("hHp_Evt%03dL1", it), "", AXnr);
				Hist::New(StrFmt("hHn_Evt%03dL1", it), "", AXnr);
				Hist::New(StrFmt("hHi_Evt%03dL9", it), "", AXir);
				Hist::New(StrFmt("hHp_Evt%03dL9", it), "", AXnr);
				Hist::New(StrFmt("hHn_Evt%03dL9", it), "", AXnr);
				Hist::New(StrFmt("hHi_Evt%03dFs", it), "", AXir);
				Hist::New(StrFmt("hHp_Evt%03dFs", it), "", AXnr);
				Hist::New(StrFmt("hHn_Evt%03dFs", it), "", AXnr);
				Hist::New(StrFmt("hHi_MCEvt%03dL1", it), "", AXir);
				Hist::New(StrFmt("hHp_MCEvt%03dL1", it), "", AXnr);
				Hist::New(StrFmt("hHn_MCEvt%03dL1", it), "", AXnr);
				Hist::New(StrFmt("hHi_MCEvt%03dL9", it), "", AXir);
				Hist::New(StrFmt("hHp_MCEvt%03dL9", it), "", AXnr);
				Hist::New(StrFmt("hHn_MCEvt%03dL9", it), "", AXnr);
				Hist::New(StrFmt("hHi_MCEvt%03dFs", it), "", AXir);
				Hist::New(StrFmt("hHp_MCEvt%03dFs", it), "", AXnr);
				Hist::New(StrFmt("hHn_MCEvt%03dFs", it), "", AXnr);
			}

			Hist::New("hHi_Evt", "", AXir);
			Hist::New("hHp_Evt", "", AXnr);
			Hist::New("hHn_Evt", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_Evt%03d", it), "", AXir);
				Hist::New(StrFmt("hHp_Evt%03d", it), "", AXnr);
				Hist::New(StrFmt("hHn_Evt%03d", it), "", AXnr);
			}
			
			Hist::New("hHi_EvtL1", "", AXir);
			Hist::New("hHp_EvtL1", "", AXnr);
			Hist::New("hHn_EvtL1", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_EvtL1%03d", it), "", AXir);
				Hist::New(StrFmt("hHp_EvtL1%03d", it), "", AXnr);
				Hist::New(StrFmt("hHn_EvtL1%03d", it), "", AXnr);
			}
			
			Hist::New("hHi_EvtL9", "", AXir);
			Hist::New("hHp_EvtL9", "", AXnr);
			Hist::New("hHn_EvtL9", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_EvtL9%03d", it), "", AXir);
				Hist::New(StrFmt("hHp_EvtL9%03d", it), "", AXnr);
				Hist::New(StrFmt("hHn_EvtL9%03d", it), "", AXnr);
			}
			
			Hist::New("hHi_EvtFs", "", AXir);
			Hist::New("hHp_EvtFs", "", AXnr);
			Hist::New("hHn_EvtFs", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_EvtFs%03d", it), "", AXir);
				Hist::New(StrFmt("hHp_EvtFs%03d", it), "", AXnr);
				Hist::New(StrFmt("hHn_EvtFs%03d", it), "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hHi_MCEvt", "", AXir);
				Hist::New("hHp_MCEvt", "", AXnr);
				Hist::New("hHn_MCEvt", "", AXnr);
				
				Hist::New("hHi_MCEvtL1", "", AXir);
				Hist::New("hHp_MCEvtL1", "", AXnr);
				Hist::New("hHn_MCEvtL1", "", AXnr);
				
				Hist::New("hHi_MCEvtL9", "", AXir);
				Hist::New("hHp_MCEvtL9", "", AXnr);
				Hist::New("hHn_MCEvtL9", "", AXnr);
				
				Hist::New("hHi_MCEvtFs", "", AXir);
				Hist::New("hHp_MCEvtFs", "", AXnr);
				Hist::New("hHn_MCEvtFs", "", AXnr);
			}

			std::cout << "\n<<  init DST End  >>\n";
			file->cd();
			return;
		}

		void resetDST() {
#if Debug == true
	std::cerr << "DST::resetDST()\n";
#endif
		}

		void fillDST() {
#if Debug == true
	std::cerr << "DST::fillDST()\n";
#endif
			if (!gSaveDST) return;
			if (fDST == 0 || !gSaveDSTTree) return;
			fDST->Fill();
		}

		void finishDST() {
			Hist::Write();
		}

	public :
		TTree * fDST;
	
	public :
		static Bool_t gSaveDST;
		static Bool_t gSaveDSTTree;
		static Bool_t gSaveDSTClone;
		static UInt_t gUTime[2]; // (pre, cur)
	
	public :
		static Axis AXnr;
		static Axis AXir;
		static Axis AXtnr;
		static Axis AXtir;
		static Axis AXtme;
		static const Float_t CfStableFT = 1.2;

		static Bool_t gTimeStudy;

	public :
		static Float_t GetStdMDR(Int_t ipatt);
		static Float_t GetStdSqrMDR(Int_t ipatt);
		static Float_t GetMDR(Int_t ipatt);
		static Float_t GetSqrMDR(Int_t ipatt);
		static Float_t GetRigSgm(Int_t ipatt, Float_t rig, Float_t mass);
};

Bool_t DST::gSaveDST = true;
Bool_t DST::gSaveDSTTree = false;
Bool_t DST::gSaveDSTClone = false;
UInt_t DST::gUTime[2] = {0, 0};

Axis DST::AXnr;
Axis DST::AXir;
Axis DST::AXtnr;
Axis DST::AXtir;
Axis DST::AXtme;

// TmStudy
Bool_t DST::gTimeStudy = true;
//Bool_t DST::gTimeStudy = false;

Float_t DST::GetStdMDR(Int_t ipatt) {
	Float_t MDR = 1800.; // (GeV)
	if      (ipatt == 2) MDR =  240.;
	if      (ipatt == 3) MDR =  540.;
	if      (ipatt == 4) MDR =  750.;
	if      (ipatt == 5) MDR = 1800.;
	return MDR;
}

Float_t DST::GetStdSqrMDR(Int_t ipatt) {
	Float_t SqrMDR = 3.24; // (TV)
	if      (ipatt == 2) SqrMDR = 5.760000e-02;
	else if (ipatt == 3) SqrMDR = 2.916000e-01;
	else if (ipatt == 4) SqrMDR = 5.625000e-01;
	else if (ipatt == 5) SqrMDR = 3.240000e+00;
	return SqrMDR;
}

Float_t DST::GetMDR(Int_t ipatt) {
	Float_t MDR = 2.0; // (TV)
	if      (ipatt == 0) MDR = 1.398205e-01;
	else if (ipatt == 1) MDR = 1.619846e-01;
	else if (ipatt == 2) MDR = 2.885081e-01;
	else if (ipatt == 3) MDR = 8.571684e-01;
	else if (ipatt == 4) MDR = 1.079867e+00;
	else if (ipatt == 5) MDR = 2.000000e+00;
	return MDR;
}

Float_t DST::GetSqrMDR(Int_t ipatt) {
	Float_t SqrMDR = 4.0; // (TV)
	if      (ipatt == 0) SqrMDR = 1.398205e-01 * 1.398205e-01;
	else if (ipatt == 1) SqrMDR = 1.619846e-01 * 1.619846e-01;
	else if (ipatt == 2) SqrMDR = 2.885081e-01 * 2.885081e-01;
	else if (ipatt == 3) SqrMDR = 8.571684e-01 * 8.571684e-01;
	else if (ipatt == 4) SqrMDR = 1.079867e+00 * 1.079867e+00;
	else if (ipatt == 5) SqrMDR = 2.000000e+00 * 2.000000e+00;
	return SqrMDR;
}

Float_t DST::GetRigSgm(Int_t ipatt, Float_t nrig, Float_t mass) {
	Float_t trSgm = 1.0;
	if      (ipatt == 0) {
		Float_t trParIU[3] = { 1.44615e-02, 1.58178e-03, 5.11515e-05 };
		Float_t trSgmIU = std::sqrt(trParIU[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIU[1] * nrig + trParIU[2] * nrig * nrig) / nrig;
		trSgm = trSgmIU;
	}
	else if (ipatt == 1) {
		Float_t trParIL[3] = { 9.60391e-03, 7.54520e-04, 3.81112e-05 };
		Float_t trSgmIL = std::sqrt(trParIL[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIL[1] * nrig + trParIL[2] * nrig * nrig) / nrig;
		trSgm = trSgmIL;
	}
	else if (ipatt == 2) {
		Float_t trParIn[3] = { 8.39123e-03, 3.36288e-04, 1.20139e-05 };
		Float_t trSgmIn = std::sqrt(trParIn[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIn[1] * nrig + trParIn[2] * nrig * nrig) / nrig;
		trSgm = trSgmIn;
	}
	else if (ipatt == 3) {
		Float_t trParL1[3] = { 7.64252e-03, 5.34561e-04, 1.36103e-06 };
		Float_t trSgmL1 = std::sqrt(trParL1[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParL1[1] * nrig + trParL1[2] * nrig * nrig) / nrig;
		trSgm = trSgmL1;
	}
	else if (ipatt == 4) {
		Float_t trParL9[3] = { 8.23647e-03, 2.97335e-04, 8.57550e-07 };
		Float_t trSgmL9 = std::sqrt(trParL9[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParL9[1] * nrig + trParL9[2] * nrig * nrig) / nrig;
		trSgm = trSgmL9;
	}
	else if (ipatt == 5) {
		Float_t trParFs[3] = { 9.29537e-03, 1.14609e-04, 2.50000e-07 };
		Float_t trSgmFs = std::sqrt(trParFs[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParFs[1] * nrig + trParFs[2] * nrig * nrig) / nrig;
		trSgm = trSgmFs;
	}
	return trSgm;
}


/*---------*/
/*  YiAna  */
/*---------*/
class YiAna : public YiNtuple, public DST {
	public :
		YiAna();
		~YiAna();

		void setBranchAddress();
		void freeBranchAddress();
		void analyzeEvent();
		bool analyzeRunInfo();
		bool analyzeCore();

	public :
		LIST * fList;
		G4MC * fG4mc;
		RTI  * fRti;
		TRG  * fTrg;
		TOF  * fTof;
		ACC  * fAcc;
		TRK  * fTrk;
		TRD  * fTrd;
		RICH * fRich;
		ECAL * fEcal;
};

YiAna::YiAna() {
	DST::fDST = 0;
	YiNtuple::fRunChain = 0;
	YiNtuple::fDataChain = 0;
	YiNtuple::init();
	YiNtuple::fTimer.start();
	fList = 0;
	fG4mc = 0;
	fRti  = 0;
	fTrg  = 0;
	fTof  = 0;
	fAcc  = 0;
	fTrk  = 0;
	fTrd  = 0;
	fRich = 0;
	fEcal = 0;
}

YiAna::~YiAna() {
	freeBranchAddress();
}


void YiAna::setBranchAddress() {
	std::cout << "\n<<  Set Branch Address Info  >>\n";

	/*******************/
	/**  USER DEFINE  **/
	/*******************/
	freeBranchAddress();
	fList = new LIST;
	fG4mc = (YiNtuple::CheckEventMode(YiNtuple::MC)) ? (new G4MC) : 0;
	fRti  = (YiNtuple::CheckEventMode(YiNtuple::ISS)) ? (new RTI) : 0;
	fTrg  = new TRG;
	fTof  = new TOF;  
	fAcc  = new ACC;  
	fTrk  = new TRK;  
	fTrd  = new TRD;  
	fRich = new RICH;
	fEcal = new ECAL;

	fDataChain->SetBranchAddress("list", &fList);
	if (YiNtuple::CheckEventMode(YiNtuple::MC))
		fDataChain->SetBranchAddress("g4mc", &fG4mc);
	if (YiNtuple::CheckEventMode(YiNtuple::ISS))
		fDataChain->SetBranchAddress("rti",  &fRti);
	fDataChain->SetBranchAddress("trg",  &fTrg);
	fDataChain->SetBranchAddress("tof",  &fTof);
	fDataChain->SetBranchAddress("acc",  &fAcc);
	fDataChain->SetBranchAddress("trk",  &fTrk);
	fDataChain->SetBranchAddress("trd",  &fTrd);
	fDataChain->SetBranchAddress("rich", &fRich);
	fDataChain->SetBranchAddress("ecal", &fEcal);
	/*******************/
	/**               **/
	/*******************/

	std::cout << "\n<<  Set Branch Address End  >>\n";
}

void YiAna::freeBranchAddress() {
	if (fList != 0) { delete fList; fList = 0; }
	if (fG4mc != 0) { delete fG4mc; fG4mc = 0; }
	if (fRti  != 0) { delete fRti ; fRti  = 0; }
	if (fTrg  != 0) { delete fTrg ; fTrg  = 0; }
	if (fTof  != 0) { delete fTof ; fTof  = 0; }
	if (fAcc  != 0) { delete fAcc ; fAcc  = 0; }
	if (fTrk  != 0) { delete fTrk ; fTrk  = 0; }
	if (fTrd  != 0) { delete fTrd ; fTrd  = 0; }
	if (fRich != 0) { delete fRich; fRich = 0; }
	if (fEcal != 0) { delete fEcal; fEcal = 0; }
}

void YiAna::analyzeEvent() {
	std::cout << "\n<<  Analyze Event Info  >>\n";

	/*******************/
	/**  USER DEFINE  **/
	/*******************/
	initDST(fFile, fRunChain, fDataChain);
	fFile->cd();
	/*******************/
	/**               **/
	/*******************/

	Long64_t nentries = fDataChain->GetEntries();
	Long64_t npassed = 0;
	Long64_t nprocessed = 0;
	Long64_t printRate = nentries / 200;
	if (printRate < 250000) printRate = 250000;
	if (printRate > 2500000) printRate = 2500000;

	for (Long64_t ientry = 0; ientry < nentries; ientry++) {
		if (nprocessed%printRate == 0) {
			const UInt_t MemSize = 1024;
			ProcInfo_t procinfo;
			gSystem->GetProcInfo(&procinfo);
			Long64_t memRes = procinfo.fMemResident / MemSize;
			Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
			fTimer.stop();

			std::cout << Form("Info :: %lf %\n", 100. * float(nprocessed)/float(nentries));
			std::cout << Form("        Processed       : %ld / %ld\n", nprocessed, nentries);
			std::cout << Form("        Passed          : %ld / %ld\n", npassed, nprocessed);
			std::cout << Form("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
			std::cout << Form("        Real Time       : %9.2f (second)\n", fTimer.time());
			std::cout << Form("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fTimer.time());
			std::cout << Form("        Cpu    System   : %4.1f %\n", procinfo.fCpuSys);
			std::cout << Form("               User     : %4.1f %\n", procinfo.fCpuUser);
			std::cout << Form("        Memory Resident : %2ld GB %4ld MB\n", memRes / MemSize, memRes % MemSize);
			std::cout << Form("               Virtual  : %2ld GB %4ld MB\n", memVrl / MemSize, memVrl % MemSize);
		}
		nprocessed++;

		fDataChain->GetEntry(ientry);

		/*******************/
		/**  USER DEFINE  **/
		/*******************/
		resetDST();
		//if (!analyzeRunInfo()) continue;
		if (!analyzeCore()) continue;
		fillDST();
		/*******************/
		/**               **/
		/*******************/

		npassed++;
	}

	/*******************/
	/**  USER DEFINE  **/
	/*******************/
	finishDST();
	/*******************/
	/**               **/
	/*******************/

	if (nprocessed == nentries) {
		const UInt_t MemSize = 1024;
		ProcInfo_t procinfo;
		gSystem->GetProcInfo(&procinfo);
		Long64_t memRes = procinfo.fMemResident / MemSize;
		Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
		fTimer.stop();
		
		std::cout << Form("Info :: %lf %\n", 100. * float(nprocessed)/float(nentries));
		std::cout << Form("        Processed       : %ld / %ld\n", nprocessed, nentries);
		std::cout << Form("        Passed          : %ld / %ld\n", npassed, nprocessed);
		std::cout << Form("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
		std::cout << Form("        Real Time       : %9.2f (second)\n", fTimer.time());
		std::cout << Form("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fTimer.time());
		std::cout << Form("        Cpu    System   : %4.1f %\n", procinfo.fCpuSys);
		std::cout << Form("               User     : %4.1f %\n", procinfo.fCpuUser);
		std::cout << Form("        Memory Resident : %2ld GB %4ld MB\n", memRes / MemSize, memRes % MemSize);
		std::cout << Form("               Virtual  : %2ld GB %4ld MB\n", memVrl / MemSize, memVrl % MemSize);
		std::cout << Form("Info :: Root Files Processed Successfully Finished.\n");
	}
	else {
		std::cout << Form("Info :: Root Files Processed Seems Failed.\n");
		std::cout << Form("        Processed %ld in %ld\n", nprocessed, nentries);
	}

	std::cout << "\n<<  Analyze Event End  >>\n";
}


/*******************/
/**  USER DEFINE  **/
/*******************/
bool YiAna::analyzeRunInfo() {
#if Debug == true
	std::cerr << "YiAna::analyzeRunInfo()\n";
#endif

	if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
		// preselection (Mini-RTI Requirement)
		Bool_t UpdateTime = false; 
		if (gUTime[0] != fRti->uTime) {
			gUTime[0] = (gUTime[1] == 0) ? gUTime[1] : fRti->uTime;
			gUTime[1] = fRti->uTime;
			UpdateTime = true;
		}

		// Time Binning
		Bool_t tvStudy = gTimeStudy && (fRti->uTime >= AXtme.min() && fRti->uTime < AXtme.max());
		Int_t tvBin = (!tvStudy) ? 0 : AXtme.find(fRti->uTime);

		// Independ on Events
		if (UpdateTime) Hist::Head("hRTI_Cutflow")->fill(fRti->uTime, 0);
		
		if (!fRti->flagRun) return false;
		if (UpdateTime) Hist::Head("hRTI_Cutflow")->fill(fRti->uTime, 1);
		
		if (!fRti->isGoodSecond) return false;
		if (UpdateTime) Hist::Head("hRTI_Cutflow")->fill(fRti->uTime, 2);
		
		if (fRti->isInSAA) return false;
		if (UpdateTime) Hist::Head("hRTI_Cutflow")->fill(fRti->uTime, 3);
		

		if (UpdateTime) Hist::Head("hRTI_LiveTime")->fill(fRti->uTime, fRti->liveTime);

		// Exposure Time
		if (UpdateTime) {
			const Int_t nDeg = 4;
			const Int_t sDeg[nDeg] = { 25, 30, 35, 40 };
			for (Int_t iDeg = 0; iDeg < nDeg; ++iDeg) {
				Float_t     cfRig = DST::CfStableFT * fRti->cutoffIGRF[iDeg];
				Int_t       nrBin = AXnr.find(cfRig) + 1;
				if (nrBin > AXnr.nbin()) continue;
				
				std::string nmSecDeg = StrFmt("hRTI_ExpT%02dDeg", sDeg[iDeg]);
				std::string nmAllDeg = StrFmt("hRTI_AllExpT%02dDeg", sDeg[iDeg]);
				for (Int_t ibin = nrBin; ibin <= AXnr.nbin(); ++ibin) {
					Double_t cen = AXnr.center(ibin, Axis::kLog);
					
					Hist::Head(nmSecDeg)->fill(fRti->uTime, cen, fRti->liveTime);
					Hist::Head(nmAllDeg)->fill(cen, fRti->liveTime);
					if (tvBin!=0) Hist::Head(StrFmt("%s%03d", nmAllDeg.c_str(), tvBin))->fill(cen, fRti->liveTime);
				}
			}
		}

		// Depend on Events
		if (fRti->isInShadow) return false;
	}

	return true;
}
/*******************/
/**               **/
/*******************/


/*******************/
/**  USER DEFINE  **/
/*******************/
bool YiAna::analyzeCore() {
#if Debug == true
	std::cerr << "YiAna::analyzeCore()\n";
#endif
	
	Long64_t runID   = fList->runID;
	Long64_t eventID = fList->eventID;
	Long64_t entryID = fList->entryID;
	Double_t weight  = fList->weight;

	// Monte Carlo
	Float_t MCTRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
	                 0.0 : (fG4mc->primPart.mom / fG4mc->primPart.chrg);
	Float_t MCNRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
	                 0.0 : std::fabs(MCTRig);
	Float_t MCIRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ?
	                 0.0 : (1. / MCTRig);
	Int_t   MCSign = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ?
	                 0 : MgntNum::Compare(MCTRig);
	Float_t MCrwgt = TMath::Power(MCNRig, -1.7);

	Int_t   MCNRigBin = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
	                    0   : AXnr.find(MCNRig);
	Float_t MCNRigCen = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
	                    1.0 : AXnr.center(MCNRigBin, Axis::kLog);

	// Time Binning
	Bool_t tvStudy = gTimeStudy && YiNtuple::CheckEventMode(YiNtuple::ISS) && (fRti->uTime >= AXtme.min() && fRti->uTime < AXtme.max());
	Int_t tvBin = (!tvStudy) ? 0 : AXtme.find(fRti->uTime);

#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] COM Selection [][]\n";
#endif
	
	// Sample
	Hist::Head("hCOM_Cutflow")->fill(0, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(0, weight);

	// preselection (Mini-RTI Requirement)
	if (!YiAna::analyzeRunInfo()) return false;
	
	Hist::Head("hCOM_Cutflow")->fill(1, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(1, weight);
	
	// preselection (Mini-Track Requirement)
	TrackInfo * track = (fTrk->tracks.size() == 1) ? (&fTrk->tracks.at(0)) : nullptr;
	if (track == nullptr || (track->bitPatt&2)!=2 || !track->status[0][2]) return false; // Inner XY
	if (MgntNum::EqualToZero(track->Qinner)) return false;
	Int_t nTrInHit = 0;
	for (Int_t ih = 0; ih < track->hits.size(); ++ih) {
		HitTRKInfo &hit = track->hits.at(ih);
		if (hit.layJ == 1 || hit.layJ == 2 || hit.layJ == 9) continue;
		nTrInHit++;
	}
	
	Hist::Head("hCOM_Cutflow")->fill(2, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(2, weight);
	
	// preselection (Min-Track Pattern Requirement) 
	Bool_t isTrIn = ((track->bitPatt& 10)== 10 && 
	                 track->status[0][2] && 
									 track->status[0][0] && track->status[0][1]);  // Inner XY + L2 XY          ( 10 = 2 +  8)
	Bool_t isTrL1 = ((track->bitPatt& 42)== 42 && 
	                 track->status[0][3]);                         // Inner XY + L1 XY + L2 XY  ( 42 = 2 +  8 +  32)
	Bool_t isTrL9 = ((track->bitPatt&138)==138 && 
	                 track->status[0][4]);                         // Inner XY + L2 XY + L9 XY  (138 = 2 +  8 + 128)
	Bool_t isTrFs = ((track->bitPatt&162)==162 && 
	                 track->status[0][5] && 
									 track->status[0][3] && track->status[0][4]);  // Inner XY + L1 XY + L9 XY  (162 = 2 + 32 + 128)
	Bool_t isTrPt = (isTrIn || isTrL1 || isTrL9 || isTrFs);

	if (!isTrPt) return false;
	
	Hist::Head("hCOM_Cutflow")->fill(3, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(3, weight);

	if (isTrIn) Hist::Head("hCOM_Cutflow")->fill(4, weight);
	if (isTrIn && tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(4, weight);

	if (isTrL1) Hist::Head("hCOM_Cutflow")->fill(5, weight);
	if (isTrL1 && tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(5, weight);

	if (isTrL9) Hist::Head("hCOM_Cutflow")->fill(6, weight);
	if (isTrL9 && tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(6, weight);

	if (isTrFs) Hist::Head("hCOM_Cutflow")->fill(7, weight);
	if (isTrFs && tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(7, weight);

	Int_t       trPt = 2;
	std::string trNm = "";
	if      (isTrFs) { trPt = 5; trNm = "Fs"; }
	else if (isTrL9) { trPt = 4; trNm = "L9"; }
	else if (isTrL1) { trPt = 3; trNm = "L1"; }
	else if (isTrIn) { trPt = 2; trNm = "In"; }

	Float_t trXY[2] = { track->state[0][trPt][0], track->state[0][trPt][1] };
	Float_t trRad   = std::sqrt(trXY[0] * trXY[0] + trXY[1] * trXY[1]); 
	Float_t trAgl   = std::acos( std::fabs(track->state[0][trPt][5]) ) * TMath::RadToDeg();
	Float_t trig    = track->rigidity[0][trPt];
	Int_t   sign    = MgntNum::Compare(trig);
	Float_t nrig    = std::fabs(trig);
	Float_t irig    = 1. / trig;

	Float_t trigIn = track->rigidity[0][2];
	Float_t nrigIn = std::fabs(trigIn);
	Float_t irigIn = 1. / trig;

	const Float_t massPr = 0.938272297; 
	const Float_t massPi = 0.139570180; 
	const Float_t massEl = 0.000510999; 
	Float_t       tbtaPr  = 1. / std::sqrt(1.+(massPr/nrig)*(massPr/nrig));
	Float_t       tbtaPi  = 1. / std::sqrt(1.+(massPi/nrig)*(massPi/nrig));
	Float_t       tbtaEl  = 1. / std::sqrt(1.+(massEl/nrig)*(massEl/nrig));

	// preselection (IGRF Cutoff Requirement) 
	Bool_t isOverCf = true;
	if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
		Float_t trIRig[4] = { (isTrIn?std::fabs(1./track->rigidity[0][2]):0.), 
	                        (isTrL1?std::fabs(1./track->rigidity[0][3]):0.), 
											  	(isTrL9?std::fabs(1./track->rigidity[0][4]):0.),
											  	(isTrFs?std::fabs(1./track->rigidity[0][5]):0.) };
		Float_t minRig = 1.0 / (*std::max_element(trIRig, trIRig+4));
		Float_t cosAgl = std::acos( std::fabs(track->stateL1[0][trPt][5]) ) * TMath::RadToDeg();
		Int_t   patAgl = -1;
		if      (cosAgl < 25.0) patAgl = 0; 
		else if (cosAgl < 30.0) patAgl = 1; 
		else if (cosAgl < 35.0) patAgl = 2; 
		else if (cosAgl < 40.0) patAgl = 3;
		if (patAgl != -1) {	
			Float_t cfRig = DST::CfStableFT * fRti->cutoffIGRF[patAgl];
			if (MgntNum::Compare(minRig, cfRig) >= 0) isOverCf = true;
		}
		else isOverCf = false;
	}
	if (!isOverCf) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 0, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 0, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 0, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
	
	// preselection (Min-Track Charge Requirement) 
	if (track->Qinner < 0.8 || track->Qinner > 1.4) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 1, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 1, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 1, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
	
	// preselection (Mini-FriducialVolume Requirement)
	Hist::Head("hCOMi_TrRadius")->fill(irig, trRad, weight);
	if (sign>0) Hist::Head("hCOMp_TrRadius")->fill(nrig, trRad, weight);
	if (sign<0) Hist::Head("hCOMn_TrRadius")->fill(nrig, trRad, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TrRadius%03d", tvBin))->fill(irig, trRad, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TrRadius%03d", tvBin))->fill(nrig, trRad, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TrRadius%03d", tvBin))->fill(nrig, trRad, weight);
	
	Hist::Head("hCOMi_TrAngle")->fill(irig, trAgl, weight);
	if (sign>0) Hist::Head("hCOMp_TrAngle")->fill(nrig, trAgl, weight);
	if (sign<0) Hist::Head("hCOMn_TrAngle")->fill(nrig, trAgl, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TrAngle%03d", tvBin))->fill(irig, trAgl, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TrAngle%03d", tvBin))->fill(nrig, trAgl, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TrAngle%03d", tvBin))->fill(nrig, trAgl, weight);

	const Float_t trFiducialVolume = 47.0;
	if (trRad > trFiducialVolume) return false;

	Hist::Head("hCOMi_Cutflow")->fill(irig, 2, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 2, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 2, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
	
	// preselection (Mini-Trigger Requirement)
	if ((fTrg->bit&8) != 8) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 3, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 3, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 3, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
	
	// preselection (Mini-TOF Requirement)
	if (!fTof->statusBetaH) return false;
	if (fTof->betaHPatt != 15) return false;
	if (fTof->betaH < 0.3 || fTof->betaH > 1.3) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 4, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 4, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 4, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 4, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
	
	// preselection (Min-Track Charge Requirement) 
	HitTRKInfo * trHitL1 = nullptr; Float_t trL1Q = 0.0; Float_t trL1Qcrr = std::sqrt(1.0+0.58*(massPr/nrig)*(massPr/nrig));
	HitTRKInfo * trHitL2 = nullptr; Float_t trL2Q = 0.0; Float_t trL2Qcrr = std::sqrt(1.0+0.71*(massPr/nrig)*(massPr/nrig));
	HitTRKInfo * trHitL9 = nullptr; Float_t trL9Q = 0.0; Float_t trL9Qcrr = std::sqrt(1.0+0.67*(massPr/nrig)*(massPr/nrig));
	HitTRKInfo * trHitMx = nullptr; Float_t trMxQ = 0.0; Float_t trMxQcrr = std::sqrt(1.0+0.71*(massPr/nrig)*(massPr/nrig));
	for (Int_t ih = 0; ih < track->hits.size(); ++ih) {
		HitTRKInfo &hit = track->hits.at(ih);
		if (hit.clsId[0] < 0 || hit.clsId[1] < 0) continue;
		Float_t chrg = 0.5 * (hit.chrg[0] + hit.chrg[1]);
		if      (hit.layJ == 1) { trHitL1 = &hit; trL1Q = (chrg / trL1Qcrr); } 
		else if (hit.layJ == 2) { trHitL2 = &hit; trL2Q = (chrg / trL2Qcrr); } 
		else if (hit.layJ == 9) { trHitL9 = &hit; trL9Q = (chrg / trL9Qcrr); }
		else if (chrg > trMxQ)  { trHitMx = &hit; trMxQ = chrg; }
  }
	if (trHitMx) trMxQ /= trMxQcrr;

	Hist::Head("hCOMi_TrInQ")->fill(irig, track->Qinner, weight);
	if (sign>0) Hist::Head("hCOMp_TrInQ")->fill(nrig, track->Qinner, weight);
	if (sign<0) Hist::Head("hCOMn_TrInQ")->fill(nrig, track->Qinner, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TrInQ%03d", tvBin))->fill(irig, track->Qinner, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TrInQ%03d", tvBin))->fill(nrig, track->Qinner, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TrInQ%03d", tvBin))->fill(nrig, track->Qinner, weight);
	
	if (trHitMx) Hist::Head("hCOMi_TrMxQ")->fill(irig, trMxQ, weight);
	if (trHitMx && sign>0) Hist::Head("hCOMp_TrMxQ")->fill(nrig, trMxQ, weight);
	if (trHitMx && sign<0) Hist::Head("hCOMn_TrMxQ")->fill(nrig, trMxQ, weight);
	if (tvBin!=0 && trHitMx) Hist::Head(StrFmt("hCOMi_TrMxQ%03d", tvBin))->fill(irig, trMxQ, weight);
	if (tvBin!=0 && trHitMx && sign>0) Hist::Head(StrFmt("hCOMp_TrMxQ%03d", tvBin))->fill(nrig, trMxQ, weight);
	if (tvBin!=0 && trHitMx && sign<0) Hist::Head(StrFmt("hCOMn_TrMxQ%03d", tvBin))->fill(nrig, trMxQ, weight);

	if (trHitMx && (trMxQ < 0.90 || trMxQ > 2.3)) return false;
	
	if (trHitL2) Hist::Head("hCOMi_TrL2Q")->fill(irig, trL2Q, weight);
	if (trHitL2 && sign>0) Hist::Head("hCOMp_TrL2Q")->fill(nrig, trL2Q, weight);
	if (trHitL2 && sign<0) Hist::Head("hCOMn_TrL2Q")->fill(nrig, trL2Q, weight);
	if (tvBin!=0 && trHitL2) Hist::Head(StrFmt("hCOMi_TrL2Q%03d", tvBin))->fill(irig, trL2Q, weight);
	if (tvBin!=0 && trHitL2 && sign>0) Hist::Head(StrFmt("hCOMp_TrL2Q%03d", tvBin))->fill(nrig, trL2Q, weight);
	if (tvBin!=0 && trHitL2 && sign<0) Hist::Head(StrFmt("hCOMn_TrL2Q%03d", tvBin))->fill(nrig, trL2Q, weight);
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 5, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 5, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 5, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 5, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
	
	if (trHitL2 && (trL2Q < 0.75 || trL2Q > 2.1)) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 6, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 6, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 6, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 6, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
	
	if (trHitL1) Hist::Head("hCOMi_TrL1Q")->fill(irig, trL1Q, weight);
	if (trHitL1 && sign>0) Hist::Head("hCOMp_TrL1Q")->fill(nrig, trL1Q, weight);
	if (trHitL1 && sign<0) Hist::Head("hCOMn_TrL1Q")->fill(nrig, trL1Q, weight);
	if (tvBin!=0 && trHitL1) Hist::Head(StrFmt("hCOMi_TrL1Q%03d", tvBin))->fill(irig, trL1Q, weight);
	if (tvBin!=0 && trHitL1 && sign>0) Hist::Head(StrFmt("hCOMp_TrL1Q%03d", tvBin))->fill(nrig, trL1Q, weight);
	if (tvBin!=0 && trHitL1 && sign<0) Hist::Head(StrFmt("hCOMn_TrL1Q%03d", tvBin))->fill(nrig, trL1Q, weight);

	if (trHitL9) Hist::Head("hCOMi_TrL9Q")->fill(irig, trL9Q, weight);
	if (trHitL9 && sign>0) Hist::Head("hCOMp_TrL9Q")->fill(nrig, trL9Q, weight);
	if (trHitL9 && sign<0) Hist::Head("hCOMn_TrL9Q")->fill(nrig, trL9Q, weight);
	if (tvBin!=0 && trHitL9) Hist::Head(StrFmt("hCOMi_TrL9Q%03d", tvBin))->fill(irig, trL9Q, weight);
	if (tvBin!=0 && trHitL9 && sign>0) Hist::Head(StrFmt("hCOMp_TrL9Q%03d", tvBin))->fill(nrig, trL9Q, weight);
	if (tvBin!=0 && trHitL9 && sign<0) Hist::Head(StrFmt("hCOMn_TrL9Q%03d", tvBin))->fill(nrig, trL9Q, weight);
	
	if (trHitL1 && (trL1Q < 0.8 || trL1Q > 2.1)) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 7, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 7, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 7, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 7, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 7, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 7, weight);
	
	if (trHitL9 && (trL9Q < 0.8 || trL9Q > 2.1)) return false;

	Hist::Head("hCOMi_Cutflow")->fill(irig, 8, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 8, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 8, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 8, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 8, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 8, weight);
	
	// preselection (Mini-TOF Requirement)
	Float_t tofb   = (fTof->betaH - tbtaPr);
	Float_t tofM   = ((1.0/fTof->betaH/fTof->betaH - 1.0) - (massPr/nrig)*(massPr/nrig));
	Float_t tofQ   = fTof->Qall;
	Float_t tofQu  = 0.50 * (fTof->Q[0] + fTof->Q[1]);
	Float_t tofQl  = 0.50 * (fTof->Q[2] + fTof->Q[3]);
	Float_t tofQul = tofQu - tofQl;

	if (fTof->normChisqT > 10) return false;

	Hist::Head("hCOMi_Cutflow")->fill(irig, 9, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 9, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 9, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 9, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 9, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 9, weight);
	
	if (fTof->normChisqC > 10) return false;

	Hist::Head("hCOMi_Cutflow")->fill(irig, 10, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 10, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 10, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 10, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 10, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 10, weight);
	
	if (fTof->numOfInTimeCluster > 4) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 10, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 11, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 11, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 11, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 11, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 11, weight);
	
	Hist::Head("hCOMi_TofQ")->fill(irig, tofQ, weight);
	if (sign>0) Hist::Head("hCOMp_TofQ")->fill(nrig, tofQ, weight);
	if (sign<0) Hist::Head("hCOMn_TofQ")->fill(nrig, tofQ, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TofQ%03d", tvBin))->fill(irig, tofQ, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TofQ%03d", tvBin))->fill(nrig, tofQ, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TofQ%03d", tvBin))->fill(nrig, tofQ, weight);

	if (tofQ  < 0.80 || tofQ  > 1.60) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 12, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 12, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 12, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 12, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 12, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 12, weight);
	
	Hist::Head("hCOMi_TofQu")->fill(irig, tofQu, weight);
	if (sign>0) Hist::Head("hCOMp_TofQu")->fill(nrig, tofQu, weight);
	if (sign<0) Hist::Head("hCOMn_TofQu")->fill(nrig, tofQu, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TofQu%03d", tvBin))->fill(irig, tofQu, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TofQu%03d", tvBin))->fill(nrig, tofQu, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TofQu%03d", tvBin))->fill(nrig, tofQu, weight);
	
	Hist::Head("hCOMi_TofQl")->fill(irig, tofQl, weight);
	if (sign>0) Hist::Head("hCOMp_TofQl")->fill(nrig, tofQl, weight);
	if (sign<0) Hist::Head("hCOMn_TofQl")->fill(nrig, tofQl, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TofQl%03d", tvBin))->fill(irig, tofQl, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TofQl%03d", tvBin))->fill(nrig, tofQl, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TofQl%03d", tvBin))->fill(nrig, tofQl, weight);
	
	Hist::Head("hCOMi_TofM")->fill(irig, tofM, weight);
	if (sign>0) Hist::Head("hCOMp_TofM")->fill(nrig, tofM, weight);
	if (sign<0) Hist::Head("hCOMn_TofM")->fill(nrig, tofM, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TofM%03d", tvBin))->fill(irig, tofM, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TofM%03d", tvBin))->fill(nrig, tofM, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TofM%03d", tvBin))->fill(nrig, tofM, weight);
	
	if (tofQu < 0.75 || tofQu > 1.8) return false;
	if (tofQl < 0.75 || tofQl > 1.8) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 13, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 13, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 13, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 13, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 13, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 13, weight);
	
	//Float_t CutTofQul = 0.250 * (TMath::Erf(1.5 * (std::sqrt(nrig) - 1.8)) + 1.0) + 0.15;
	Float_t CutTofQul = 0.30 * (TMath::Erf(1.5 * (std::sqrt(nrig) - 1.8)) + 1.0) + 0.10;
	
	Hist::Head("hCOMi_TofQul")->fill(irig, tofQul, weight);
	if (sign>0) Hist::Head("hCOMp_TofQul")->fill(nrig, tofQul, weight);
	if (sign<0) Hist::Head("hCOMn_TofQul")->fill(nrig, tofQul, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TofQul%03d", tvBin))->fill(irig, tofQul, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TofQul%03d", tvBin))->fill(nrig, tofQul, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TofQul%03d", tvBin))->fill(nrig, tofQul, weight);

	if (std::fabs(tofQul) > CutTofQul) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 14, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 14, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 14, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 14, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 14, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 14, weight);
	
	// preselection (Mini-Acc Requirement)
	Hist::Head("hCOMi_NAcc")->fill(irig, fAcc->numOfCluster, weight);
	if (sign>0) Hist::Head("hCOMp_NAcc")->fill(nrig, fAcc->numOfCluster, weight);
	if (sign<0) Hist::Head("hCOMn_NAcc")->fill(nrig, fAcc->numOfCluster, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_NAcc%03d", tvBin))->fill(irig, fAcc->numOfCluster, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_NAcc%03d", tvBin))->fill(nrig, fAcc->numOfCluster, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_NAcc%03d", tvBin))->fill(nrig, fAcc->numOfCluster, weight);
	
	Int_t nAccCut = 4;
	if (fAcc->numOfCluster > nAccCut) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 15, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 15, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 15, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 15, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 15, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 15, weight);
	
	// preselection (Mini-TRD Requirement)
	Bool_t isCleanTrdTrack = 
		(fTrd->numOfTrack <= 2 && fTrd->numOfHTrack <= 2) && 
		(fTrd->numOfTrack == 1 || fTrd->numOfHTrack == 1);
	if (!isCleanTrdTrack) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 16, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 16, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 16, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 16, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 16, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 16, weight);
	
	Hist::Head("hCOMi_NTrdCls")->fill(irig, fTrd->numOfCluster, weight);
	if (sign>0) Hist::Head("hCOMp_NTrdCls")->fill(nrig, fTrd->numOfCluster, weight);
	if (sign<0) Hist::Head("hCOMn_NTrdCls")->fill(nrig, fTrd->numOfCluster, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_NTrdCls%03d", tvBin))->fill(irig, fTrd->numOfCluster, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_NTrdCls%03d", tvBin))->fill(nrig, fTrd->numOfCluster, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_NTrdCls%03d", tvBin))->fill(nrig, fTrd->numOfCluster, weight);
	
	Hist::Head("hCOMi_NTrkTrdVtx")->fill(irig, fTrd->numOfVertexWithTrTrack, weight);
	if (sign>0) Hist::Head("hCOMp_NTrkTrdVtx")->fill(nrig, fTrd->numOfVertexWithTrTrack, weight);
	if (sign<0) Hist::Head("hCOMn_NTrkTrdVtx")->fill(nrig, fTrd->numOfVertexWithTrTrack, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_NTrkTrdVtx%03d", tvBin))->fill(irig, fTrd->numOfVertexWithTrTrack, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_NTrkTrdVtx%03d", tvBin))->fill(nrig, fTrd->numOfVertexWithTrTrack, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_NTrkTrdVtx%03d", tvBin))->fill(nrig, fTrd->numOfVertexWithTrTrack, weight);
	
	Int_t trdPt = -1; // trd 0, trk 1
	if      (fTrd->statusKCls[0] && fTrd->LLR_nhit[0] >= 8) trdPt = 0;
	else if (fTrd->statusKCls[1] && fTrd->LLR_nhit[1] >= 8) trdPt = 1;
	Bool_t  hasTrd  = (trdPt != -1); 
	Float_t trdl    = (hasTrd) ? (fTrd->LLR[trdPt][0]) : 0.0;
	Float_t trdQ    = (hasTrd) ? (fTrd->Q[trdPt]) : 0.0;
	Bool_t  istrdHe = (hasTrd) ? (fTrd->LLR[trdPt][2] > 0.3) : false;
	Bool_t  istrdPr = (hasTrd) ? (trdl > (0.50 + 0.25 * tbtaPr)) : false;
	Bool_t  istrdSQ = (hasTrd) ? (trdQ < 0.4 || trdQ > (1.45 - 0.30 / std::log(2.72+nrig))) : false;
	if (!hasTrd) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 17, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 17, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 17, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 17, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 17, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 19, weight);
	
	Hist::Head("hCOMi_TrdQ")->fill(irig, trdQ, weight);
	if (sign>0) Hist::Head("hCOMp_TrdQ")->fill(nrig, trdQ, weight);
	if (sign<0) Hist::Head("hCOMn_TrdQ")->fill(nrig, trdQ, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TrdQ%03d", tvBin))->fill(irig, trdQ, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_TrdQ%03d", tvBin))->fill(nrig, trdQ, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_TrdQ%03d", tvBin))->fill(nrig, trdQ, weight);
	
	if (istrdHe || istrdSQ) return false;

	Hist::Head("hCOMi_Trdl")->fill(irig, trdl, weight);
	if (sign>0) Hist::Head("hCOMp_Trdl")->fill(nrig, trdl, weight);
	if (sign<0) Hist::Head("hCOMn_Trdl")->fill(nrig, trdl, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 18, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 18, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 18, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%03d", tvBin))->fill(irig, 18, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%03d", tvBin))->fill(nrig, 18, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%03d", tvBin))->fill(nrig, 18, weight);

	// preselection (Mini-Shower Requirement)
	ShowerInfo * shower = (fEcal->showers.size() >= 1) ? (&fEcal->showers.at(0)) : nullptr;
	Bool_t  hasEcal  = (shower != nullptr);
	Float_t ecalbdt  = (!hasEcal) ? -2.0  : shower->PisaBDT;
	Bool_t  isecalPr = (!hasEcal) ? false : (ecalbdt < -0.8);
	Bool_t  isecalEl = (!hasEcal) ? false : (ecalbdt >  0.6);
	
	if (hasEcal) Hist::Head("hCOMi_EcalBDT")->fill(irig, ecalbdt, weight);
	if (hasEcal && sign>0) Hist::Head("hCOMp_EcalBDT")->fill(nrig, ecalbdt, weight);
	if (hasEcal && sign<0) Hist::Head("hCOMn_EcalBDT")->fill(nrig, ecalbdt, weight);
	if (tvBin!=0 && hasEcal) Hist::Head(StrFmt("hCOMi_EcalBDT%03d", tvBin))->fill(irig, ecalbdt, weight);
	if (tvBin!=0 && hasEcal && sign>0) Hist::Head(StrFmt("hCOMp_EcalBDT%03d", tvBin))->fill(nrig, ecalbdt, weight);
	if (tvBin!=0 && hasEcal && sign<0) Hist::Head(StrFmt("hCOMn_EcalBDT%03d", tvBin))->fill(nrig, ecalbdt, weight);

	// preselection (Mini-Rich Requirement)
	Int_t   richph    = fRich->numOfHit;
	Float_t richelph  = (fRich->kindOfRad == -1) ? 0.0 : fRich->numOfExpPE[0];
	Float_t richprph  = (fRich->kindOfRad == -1) ? 0.0 : fRich->numOfExpPE[3];
	
	Bool_t  hasRich   = (fRich->status && fRich->kindOfRad != -1);
	Float_t richbPr   = (!hasRich) ? 0.0 : (1./fRich->beta - 1./tbtaPr);
	Float_t richbPi   = (!hasRich) ? 0.0 : (1./fRich->beta - 1./tbtaPi);
	Float_t richbEl   = (!hasRich) ? 0.0 : (1./fRich->beta - 1./tbtaEl);
	Float_t richMPr   = (!hasRich) ? 0.0 : ((1.0/fRich->beta/fRich->beta - 1.0) - (massPr/nrig)*(massPr/nrig));
	Float_t richMPi   = (!hasRich) ? 0.0 : ((1.0/fRich->beta/fRich->beta - 1.0) - (massPi/nrig)*(massPi/nrig));
	Float_t richMEl   = (!hasRich) ? 0.0 : ((1.0/fRich->beta/fRich->beta - 1.0) - (massEl/nrig)*(massEl/nrig));

	Float_t richbTh   = (fRich->kindOfRad == 0) ? 0.003 : 0.007;
	Float_t richbSh   = (std::sqrt(1.0 + 0.075 * (massPr/nrig)*(massPr/nrig)) - 1.0);
	//Bool_t  isrichNo  = (fRich->kindOfRad == -1) ? false : (!fRich->status && fRich->kindOfRad == 0 && MgntNum::Compare(richelph, 2.0f) > 0 && MgntNum::Compare(richprph, 0.5f) < 0 && richph <= 2);
	//Bool_t  isrichPr  = (fRich->kindOfRad == -1) ? false : ( fRich->status && richph <= 25 && std::fabs(richbPr) < (richbTh + richbSh));
	//Bool_t  isrichPi  = (fRich->kindOfRad == -1) ? false : ( fRich->status && richph <= 25 && std::fabs(richbPi) < (richbTh));
	//Bool_t  isrichEl  = (fRich->kindOfRad == -1) ? false : ( fRich->status && richph <= 25 && std::fabs(richbEl) < (richbTh));
	
	Bool_t  isrichNo  = (fRich->kindOfRad == -1) ? false : (!fRich->status && fRich->kindOfRad == 0 && MgntNum::Compare(richelph, 2.0f) > 0 && MgntNum::Compare(richprph, 0.5f) < 0 && richph == 0);
	Bool_t  isrichPr  = (fRich->kindOfRad == -1) ? false : ( fRich->status && std::fabs(richbPr) < (richbTh + richbSh));
	Bool_t  isrichPi  = (fRich->kindOfRad == -1) ? false : ( fRich->status && std::fabs(richbPi) < (richbTh));
	Bool_t  isrichEl  = (fRich->kindOfRad == -1) ? false : ( fRich->status && std::fabs(richbEl) < (richbTh));

	if (hasRich && fRich->kindOfRad == 0) Hist::Head("hCOMi_AglbPr")->fill(irig, richbPr, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head("hCOMp_AglbPr")->fill(nrig, richbPr, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head("hCOMn_AglbPr")->fill(nrig, richbPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0) Hist::Head(StrFmt("hCOMi_AglbPr%03d", tvBin))->fill(irig, richbPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head(StrFmt("hCOMp_AglbPr%03d", tvBin))->fill(nrig, richbPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head(StrFmt("hCOMn_AglbPr%03d", tvBin))->fill(nrig, richbPr, weight);
	
	if (hasRich && fRich->kindOfRad == 1) Hist::Head("hCOMi_NafbPr")->fill(irig, richbPr, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head("hCOMp_NafbPr")->fill(nrig, richbPr, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head("hCOMn_NafbPr")->fill(nrig, richbPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1) Hist::Head(StrFmt("hCOMi_NafbPr%03d", tvBin))->fill(irig, richbPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head(StrFmt("hCOMp_NafbPr%03d", tvBin))->fill(nrig, richbPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head(StrFmt("hCOMn_NafbPr%03d", tvBin))->fill(nrig, richbPr, weight);
	
	if (hasRich && fRich->kindOfRad == 0) Hist::Head("hCOMi_AglbEl")->fill(irig, richbEl, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head("hCOMp_AglbEl")->fill(nrig, richbEl, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head("hCOMn_AglbEl")->fill(nrig, richbEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0) Hist::Head(StrFmt("hCOMi_AglbEl%03d", tvBin))->fill(irig, richbEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head(StrFmt("hCOMp_AglbEl%03d", tvBin))->fill(nrig, richbEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head(StrFmt("hCOMn_AglbEl%03d", tvBin))->fill(nrig, richbEl, weight);

	if (hasRich && fRich->kindOfRad == 1) Hist::Head("hCOMi_NafbEl")->fill(irig, richbEl, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head("hCOMp_NafbEl")->fill(nrig, richbEl, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head("hCOMn_NafbEl")->fill(nrig, richbEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1) Hist::Head(StrFmt("hCOMi_NafbEl%03d", tvBin))->fill(irig, richbEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head(StrFmt("hCOMp_NafbEl%03d", tvBin))->fill(nrig, richbEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head(StrFmt("hCOMn_NafbEl%03d", tvBin))->fill(nrig, richbEl, weight);
	
	if (hasRich && fRich->kindOfRad == 0) Hist::Head("hCOMi_AglMPr")->fill(irig, richMPr, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head("hCOMp_AglMPr")->fill(nrig, richMPr, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head("hCOMn_AglMPr")->fill(nrig, richMPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0) Hist::Head(StrFmt("hCOMi_AglMPr%03d", tvBin))->fill(irig, richMPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head(StrFmt("hCOMp_AglMPr%03d", tvBin))->fill(nrig, richMPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head(StrFmt("hCOMn_AglMPr%03d", tvBin))->fill(nrig, richMPr, weight);
	
	if (hasRich && fRich->kindOfRad == 1) Hist::Head("hCOMi_NafMPr")->fill(irig, richMPr, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head("hCOMp_NafMPr")->fill(nrig, richMPr, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head("hCOMn_NafMPr")->fill(nrig, richMPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1) Hist::Head(StrFmt("hCOMi_NafMPr%03d", tvBin))->fill(irig, richMPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head(StrFmt("hCOMp_NafMPr%03d", tvBin))->fill(nrig, richMPr, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head(StrFmt("hCOMn_NafMPr%03d", tvBin))->fill(nrig, richMPr, weight);
	
	if (hasRich && fRich->kindOfRad == 0) Hist::Head("hCOMi_AglMEl")->fill(irig, richMEl, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head("hCOMp_AglMEl")->fill(nrig, richMEl, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head("hCOMn_AglMEl")->fill(nrig, richMEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0) Hist::Head(StrFmt("hCOMi_AglMEl%03d", tvBin))->fill(irig, richMEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head(StrFmt("hCOMp_AglMEl%03d", tvBin))->fill(nrig, richMEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head(StrFmt("hCOMn_AglMEl%03d", tvBin))->fill(nrig, richMEl, weight);

	if (hasRich && fRich->kindOfRad == 1) Hist::Head("hCOMi_NafMEl")->fill(irig, richMEl, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head("hCOMp_NafMEl")->fill(nrig, richMEl, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head("hCOMn_NafMEl")->fill(nrig, richMEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1) Hist::Head(StrFmt("hCOMi_NafMEl%03d", tvBin))->fill(irig, richMEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head(StrFmt("hCOMp_NafMEl%03d", tvBin))->fill(nrig, richMEl, weight);
	if (tvBin!=0 && hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head(StrFmt("hCOMn_NafMEl%03d", tvBin))->fill(nrig, richMEl, weight);
	
	// Charge Confusion
	Float_t tuneCCPnt = 0.0;
	if      (trPt == 2) tuneCCPnt = 40.0;
	else if (trPt == 3) tuneCCPnt = 40.0;
	else if (trPt == 4) tuneCCPnt = 60.0;
	else if (trPt == 5) tuneCCPnt = 75.0;
	
	Hist::Head(StrFmt("hCOMi_Cutflow%s", trNm.c_str()))->fill(irig, 0, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s", trNm.c_str()))->fill(nrig, 0, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s", trNm.c_str()))->fill(nrig, 0, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(irig, 0, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 0, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 0, weight);
		
	// Charge Confusion ---- Chisq
	Float_t lchix    = std::log(track->chisq[0][trPt][0]);
	Float_t lchiy    = std::log(track->chisq[0][trPt][1]);
	Float_t lchiWgt  = GetSqrMDR(trPt);
	Float_t lchixCut   = 3.0; // ~98%
	Float_t lchiyCutTh[4] = { 3.6, 3.4, 2.7, 2.6 }; // ~99%
	//Float_t lchiyCut98 = 0.0; // ~98%
	//Float_t lchiyCut95 = 0.0; // ~95%
	//if      (trPt == 2) { lchiyCut98 = 3.0; lchiyCut95 = 2.3; } 
	//else if (trPt == 3) { lchiyCut98 = 2.8; lchiyCut95 = 2.1; }
	//else if (trPt == 4) { lchiyCut98 = 2.3; lchiyCut95 = 1.9; }
	//else if (trPt == 5) { lchiyCut98 = 2.3; lchiyCut95 = 1.9; }
	Float_t lchiyCut98 = 0.0; // ~98%
	Float_t lchiyCut95 = 0.0; // ~95%
	if      (trPt == 2) { lchiyCut98 = 2.35; lchiyCut95 = 1.95; } 
	else if (trPt == 3) { lchiyCut98 = 2.20; lchiyCut95 = 1.80; }
	else if (trPt == 4) { lchiyCut98 = 2.00; lchiyCut95 = 1.65; }
	else if (trPt == 5) { lchiyCut98 = 2.20; lchiyCut95 = 1.75; }
	//Float_t lchiyCut = std::sqrt(lchiyCut98 * lchiyCut95) * (0.5 * (std::erf((std::log(GetStdMDR(trPt)) - std::log(nrig)) / 1.5) + 1.0));
	Float_t lchiyCut = lchiyCut98;
	
	Float_t lchiyIU = (track->status[0]) ? std::log(track->chisq[0][0][1]) : 0.0;
	Float_t lchiyIL = (track->status[1]) ? std::log(track->chisq[0][1][1]) : 0.0;
	Float_t lchiyIn = (track->status[2]) ? std::log(track->chisq[0][2][1]) : 0.0;
	Float_t lchiyL1 = (track->status[3]) ? std::log(track->chisq[0][3][1]) : 0.0;
	Float_t lchiyL9 = (track->status[4]) ? std::log(track->chisq[0][4][1]) : 0.0;
	Float_t lchiyFs = (track->status[5]) ? std::log(track->chisq[0][5][1]) : 0.0;

	// Charge Confusion ---- AsymRig
	std::string lasymNm    = "";
	Float_t     lasym      = 0.0;
	Float_t     lasymWgt   = 0.0;
	Float_t     lasymCutTh[4] = { 1.75, 1.90, 1.60, 1.90 }; // ~99%
	//Float_t     lasymCut98 = 0.0; // ~98%
	//Float_t     lasymCut95 = 0.0; // ~95%
	//if      (trPt == 2) { lasymCut98 = 1.45; lasymCut95 = 1.05; }
	//else if (trPt == 3) { lasymCut98 = 1.70; lasymCut95 = 1.35; }
	//else if (trPt == 4) { lasymCut98 = 1.40; lasymCut95 = 1.10; }
	//else if (trPt == 5) { lasymCut98 = 1.50; lasymCut95 = 1.30; }
	Float_t     lasymCut98 = 0.0; // ~98%
	Float_t     lasymCut95 = 0.0; // ~95%
	if      (trPt == 2) { lasymCut98 = 1.25; lasymCut95 = 0.95; }
	else if (trPt == 3) { lasymCut98 = 1.55; lasymCut95 = 1.20; }
	else if (trPt == 4) { lasymCut98 = 1.35; lasymCut95 = 1.05; }
	else if (trPt == 5) { lasymCut98 = 1.70; lasymCut95 = 1.30; }
	//Float_t lasymCut = std::sqrt(lasymCut98 * lasymCut95) * (0.5 * (std::erf((std::log(GetStdMDR(trPt)) - std::log(nrig)) / 1.5) + 1.0));
	Float_t lasymCut = lasymCut98;

	const Float_t lasymLMT = 1.0e-8;
	const Float_t lasymSGM = 1.8;
	Float_t trSgmIU = GetRigSgm(0, nrig, massPr);
	Float_t trSgmIL = GetRigSgm(1, nrig, massPr);
	Float_t trSgmIn = GetRigSgm(2, nrig, massPr);
	Float_t trSgmL1 = GetRigSgm(3, nrig, massPr);
	Float_t trSgmL9 = GetRigSgm(4, nrig, massPr);
	Float_t trSgmFs = GetRigSgm(5, nrig, massPr);

	Float_t crrParUL[6] = { 7.63668e-01, 2.13220e-01, -1.40775e+00, 9.60481e-01, 3.50012e+00, 1.49281e+00 };
	Float_t crrSgmUL    = crrParUL[0]*(std::erf(crrParUL[1]*(std::log(nrig)+crrParUL[2]))+1.0) + TMath::Landau(crrParUL[3], crrParUL[4], crrParUL[5]);
	Float_t trSgmUL     = std::sqrt(trSgmIU * trSgmIU + trSgmIL * trSgmIL) * crrSgmUL;
	Float_t asymUL      = (1.0/track->rigidity[0][0] - 1.0/track->rigidity[0][1]) / trSgmUL;
	Float_t lasymUL     = std::log(asymUL * asymUL + lasymLMT) / lasymSGM;
	if (trPt == 2) {
		lasymNm  = "UL";
		lasym    = lasymUL;
		lasymWgt = GetStdSqrMDR(2);
	}

	Float_t crrPar1I[5] = { 4.93251e-01, 1.90417e+00, 9.18502e+01, 4.86389e+00, 3.29747e-01 };
	Float_t crrSgm1I    = crrPar1I[0]*std::erf(crrPar1I[1]*(std::log(nrig+crrPar1I[2])-crrPar1I[3]))+crrPar1I[4];
	Float_t trSgm1I     = std::sqrt(trSgmIn * trSgmIn + trSgmL1 * trSgmL1) * crrSgm1I;
	Float_t asym1I      = (1.0/track->rigidity[0][2] - 1.0/track->rigidity[0][3]) / trSgm1I; 
	Float_t lasym1I     = std::log(asym1I * asym1I + lasymLMT) / lasymSGM;
	if (trPt == 3) {
		lasymNm  = "1I";
		lasym    = lasym1I;
		lasymWgt = GetStdSqrMDR(3);
	}
	
	Float_t crrPar9I[5] = { 4.17294e-01, 1.15603e+00, 1.62194e+01, 3.81780e+00, 3.92404e-01 };
	Float_t crrSgm9I    = crrPar9I[0]*std::erf(crrPar9I[1]*(std::log(nrig+crrPar9I[2])-crrPar9I[3]))+crrPar9I[4];
	Float_t trSgm9I     = std::sqrt(trSgmIn * trSgmIn + trSgmL9 * trSgmL9) * crrSgm9I;
	Float_t asym9I      = (1.0/track->rigidity[0][2] - 1.0/track->rigidity[0][4]) / trSgm9I;
	Float_t lasym9I     = std::log(asym9I * asym9I + lasymLMT) / lasymSGM;
	if (trPt == 4) {
		lasymNm  = "9I";
		lasym    = lasym9I;
		lasymWgt = GetStdSqrMDR(4);
	}

	Float_t crrPar91[5] = { 5.25352e-01, 9.88059e-01, 1.21200e+01, 3.72027e+00, 4.99424e-01 };
	Float_t crrSgm91    = crrPar91[0]*std::erf(crrPar91[1]*(std::log(nrig+crrPar91[2])-crrPar91[3]))+crrPar91[4];
	Float_t trSgm91     = std::sqrt(trSgmL1 * trSgmL1 + trSgmL9 * trSgmL9) * crrSgm91;
	Float_t asym91      = (1.0/track->rigidity[0][3] - 1.0/track->rigidity[0][4]) / trSgm91;
	Float_t lasym91     = std::log(asym91 * asym91 + lasymLMT) / lasymSGM;
	if (trPt == 5) {
		lasymNm  = "91";
		lasym    = lasym91;
		lasymWgt = GetStdSqrMDR(5);
	}
	
	// Charge Confusion ---- CC Estimator
	Float_t ccest    = ((lchiy * lchiWgt + lasym * lasymWgt) / (lchiWgt + lasymWgt));
	Float_t ccestCut = 2.0 * std::erf((std::log(tuneCCPnt) - std::log(nrig)) / 1.5);

	Bool_t ccUL = (lchiyIU > lchiyCutTh[0] || lchiyIL > lchiyCutTh[0] || 
	               (lchiyIU > 0 && lchiyIL > 0 && ((lchiyIU*lchiyIU + lchiyIL*lchiyIL)/lchiyCutTh[0]/lchiyCutTh[0]) > 1) );
	Bool_t ccIn = (lchiyIn > lchiyCutTh[0] || lasymUL > lasymCutTh[0] ||
	               (lchiyIn > 0 && lasymUL > 0 && (lchiyIn*lchiyIn/lchiyCutTh[0]/lchiyCutTh[0] + lasymUL*lasymUL/lasymCutTh[0]/lasymCutTh[0]) > 1) );
	Bool_t ccL1 = (lchiyL1 > lchiyCutTh[1] || lasym1I > lasymCutTh[1] ||
	               (lchiyL1 > 0 && lasym1I > 0 && (lchiyL1*lchiyL1/lchiyCutTh[1]/lchiyCutTh[1] + lasym1I*lasym1I/lasymCutTh[1]/lasymCutTh[1]) > 1) );
	Bool_t ccL9 = (lchiyL9 > lchiyCutTh[2] || lasym9I > lasymCutTh[2] ||
	               (lchiyL9 > 0 && lasym9I > 0 && (lchiyL9*lchiyL9/lchiyCutTh[2]/lchiyCutTh[2] + lasym9I*lasym9I/lasymCutTh[2]/lasymCutTh[2]) > 1) );
	Bool_t ccFs = (lchiyFs > lchiyCutTh[3] || lasym91 > lasymCutTh[3] ||
	               (lchiyFs > 0 && lasym91 > 0 && (lchiyFs*lchiyFs/lchiyCutTh[3]/lchiyCutTh[3] + lasym91*lasym91/lasymCutTh[3]/lasymCutTh[3]) > 1) );

	Bool_t isPrLike = (istrdPr && (hasEcal ? isecalPr : true) && (hasRich ? isrichPr : true));

	if (isPrLike) {
		Hist::Head(StrFmt("hCOMi_Lchix%s", trNm.c_str()))->fill(irig, lchix, weight);
		if (sign>0) Hist::Head(StrFmt("hCOMp_Lchix%s", trNm.c_str()))->fill(nrig, lchix, weight);
		if (sign<0) Hist::Head(StrFmt("hCOMn_Lchix%s", trNm.c_str()))->fill(nrig, lchix, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Lchix%s%03d", trNm.c_str(), tvBin))->fill(irig, lchix, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Lchix%s%03d", trNm.c_str(), tvBin))->fill(nrig, lchix, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Lchix%s%03d", trNm.c_str(), tvBin))->fill(nrig, lchix, weight);
	}

	if (lchix > lchixCut) return false;
	
	Hist::Head(StrFmt("hCOMi_Cutflow%s", trNm.c_str()))->fill(irig, 1, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s", trNm.c_str()))->fill(nrig, 1, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s", trNm.c_str()))->fill(nrig, 1, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(irig, 1, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 1, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 1, weight);
	
	if (trPt == 2 && ccUL) return false;
	if (trPt == 3 && ccIn) return false;
	if (trPt == 4 && ccIn) return false;
	if (trPt == 5 && (ccL1 || ccL9)) return false;

	Hist::Head(StrFmt("hCOMi_Cutflow%s", trNm.c_str()))->fill(irig, 2, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s", trNm.c_str()))->fill(nrig, 2, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s", trNm.c_str()))->fill(nrig, 2, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(irig, 2, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 2, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 2, weight);
	
	if (trPt == 2 && ccIn) return false;
	if (trPt == 3 && ccL1) return false;
	if (trPt == 4 && ccL9) return false;
	if (trPt == 5 && ccFs) return false;
	
	Hist::Head(StrFmt("hCOMi_Cutflow%s", trNm.c_str()))->fill(irig, 3, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s", trNm.c_str()))->fill(nrig, 3, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s", trNm.c_str()))->fill(nrig, 3, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(irig, 3, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 3, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 3, weight);
	
	if (isPrLike) {
		Hist::Head(StrFmt("hCOMi_Lchiy%s", trNm.c_str()))->fill(irig, lchiy, weight);
		if (sign>0) Hist::Head(StrFmt("hCOMp_Lchiy%s", trNm.c_str()))->fill(nrig, lchiy, weight);
		if (sign<0) Hist::Head(StrFmt("hCOMn_Lchiy%s", trNm.c_str()))->fill(nrig, lchiy, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Lchiy%s%03d", trNm.c_str(), tvBin))->fill(irig, lchiy, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Lchiy%s%03d", trNm.c_str(), tvBin))->fill(nrig, lchiy, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Lchiy%s%03d", trNm.c_str(), tvBin))->fill(nrig, lchiy, weight);
	}	
	
	if (isPrLike) {
		Hist::Head(StrFmt("hCOMi_Lasym%s", lasymNm.c_str()))->fill(irig, lasym, weight);
		if (sign>0) Hist::Head(StrFmt("hCOMp_Lasym%s", lasymNm.c_str()))->fill(nrig, lasym, weight);
		if (sign<0) Hist::Head(StrFmt("hCOMn_Lasym%s", lasymNm.c_str()))->fill(nrig, lasym, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Lasym%s%03d", lasymNm.c_str(), tvBin))->fill(irig, lasym, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Lasym%s%03d", lasymNm.c_str(), tvBin))->fill(nrig, lasym, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Lasym%s%03d", lasymNm.c_str(), tvBin))->fill(nrig, lasym, weight);
	}
	
	if (lchiy > lchiyCut) return false;
	
	Hist::Head(StrFmt("hCOMi_Cutflow%s", trNm.c_str()))->fill(irig, 4, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s", trNm.c_str()))->fill(nrig, 4, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s", trNm.c_str()))->fill(nrig, 4, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(irig, 4, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 4, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 4, weight);
	
	if (lasym > lasymCut) return false;
	
	Hist::Head(StrFmt("hCOMi_Cutflow%s", trNm.c_str()))->fill(irig, 5, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s", trNm.c_str()))->fill(nrig, 5, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s", trNm.c_str()))->fill(nrig, 5, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(irig, 5, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 5, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 5, weight);
	
	//if (lchiy > 0 && lasym > 0 && ((lchiy/lchiyCut)*(lchiy/lchiyCut) + (lasym/lasymCut)*(lasym/lasymCut)) > 1) return false;
	
	if (isPrLike) {
		Hist::Head(StrFmt("hCOMi_CCest%s", trNm.c_str()))->fill(irig, ccest, weight);
		if (sign>0) Hist::Head(StrFmt("hCOMp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hCOMn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hCOMi_CCest%s%03d", trNm.c_str(), tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_CCest%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_CCest%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
	}
	
	Hist::Head(StrFmt("hCOMi_Cutflow%s", trNm.c_str()))->fill(irig, 6, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s", trNm.c_str()))->fill(nrig, 6, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s", trNm.c_str()))->fill(nrig, 6, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(irig, 6, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 6, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Cutflow%s%03d", trNm.c_str(), tvBin))->fill(nrig, 6, weight);
		
	Hist::Head(StrFmt("hCOMi_Evt%s", trNm.c_str()))->fill(irig, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Evt%s", trNm.c_str()))->fill(nrig, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Evt%s", trNm.c_str()))->fill(nrig, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_Evt%s%03d", trNm.c_str(), tvBin))->fill(irig, weight);
	if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hCOMp_Evt%s%03d", trNm.c_str(), tvBin))->fill(nrig, weight);
	if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hCOMn_Evt%s%03d", trNm.c_str(), tvBin))->fill(nrig, weight);

	if (YiNtuple::CheckEventMode(YiNtuple::MC) && isPrLike) {
		Hist::Head(StrFmt("hCOMi_MCTrRso%s", trNm.c_str()))->fill(MCIRig, irig, weight);
		if (MCSign >0) Hist::Head(StrFmt("hCOMp_MCTrRso%s", trNm.c_str()))->fill(MCNRig, MCNRigCen * (irig - MCIRig), weight);
		if (MCSign <0) Hist::Head(StrFmt("hCOMn_MCTrRso%s", trNm.c_str()))->fill(MCNRig, MCNRigCen * (irig - MCIRig), weight);
	
		Hist::Head(StrFmt("hCOMi_MCEvt%s", trNm.c_str()))->fill(MCIRig, weight);
		if (MCSign >0) Hist::Head(StrFmt("hCOMp_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
		if (MCSign <0) Hist::Head(StrFmt("hCOMn_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
	}

#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] L Selection [][]\n";
#endif
	/**** Low Energy Region ****/
	while (true) {
		// TRK
		Hist::Head("hLi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		
		// ECAL
		if (hasEcal && !isecalPr) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		// TRD
		if (!istrdPr) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);

		// Charge Confusion
		Hist::Head("hLi_CCest")->fill(irig, ccest, weight);
		if (sign>0) Hist::Head("hLp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hLn_CCest")->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_CCest%03d", tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		
		//if (ccest > ccestCut) break;

		Hist::Head("hLi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		
		// RICH
		if (fRich->kindOfRad == -1) break;
		Float_t aglPhTh   = std::min(Int_t(1.0+2.5*nrig), 16);
		Float_t nafPhTh   = std::min(Int_t(1.0+2.5*nrig),  6);
		
		std::string richPt  = (fRich->status) ? "HV" : "NO";
		std::string richNm  = (fRich->kindOfRad == 0) ? "Agl" : "Naf";
		Bool_t      isElph  = (fRich->kindOfRad == 0) ? (MgntNum::Compare(richelph, 2.0f) > 0) : (MgntNum::Compare(richelph, 0.5f) > 0);
		Bool_t      isPhTh  = (fRich->kindOfRad == 0) ? (richph <= aglPhTh) : (richph <= nafPhTh); 
		
		Hist::Head(StrFmt("hLi_%sprph%s", richNm.c_str(), richPt.c_str()))->fill(irig, richprph, weight);
		if (sign>0) Hist::Head(StrFmt("hLp_%sprph%s", richNm.c_str(), richPt.c_str()))->fill(nrig, richprph, weight);
		if (sign<0) Hist::Head(StrFmt("hLn_%sprph%s", richNm.c_str(), richPt.c_str()))->fill(nrig, richprph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_%sprph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(irig, richprph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_%sprph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(nrig, richprph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_%sprph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(nrig, richprph, weight);
		
		Hist::Head(StrFmt("hLi_%selph%s", richNm.c_str(), richPt.c_str()))->fill(irig, richelph, weight);
		if (sign>0) Hist::Head(StrFmt("hLp_%selph%s", richNm.c_str(), richPt.c_str()))->fill(nrig, richelph, weight);
		if (sign<0) Hist::Head(StrFmt("hLn_%selph%s", richNm.c_str(), richPt.c_str()))->fill(nrig, richelph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_%selph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(irig, richelph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_%selph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(nrig, richelph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_%selph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(nrig, richelph, weight);
		
		Hist::Head("hLi_Cutflow")->fill(irig, 4, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 4, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 4, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 4, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		
		Hist::Head(StrFmt("hLi_%sph%s", richNm.c_str(), richPt.c_str()))->fill(irig, richph, weight);
		if (sign>0) Hist::Head(StrFmt("hLp_%sph%s", richNm.c_str(), richPt.c_str()))->fill(nrig, richph, weight);
		if (sign<0) Hist::Head(StrFmt("hLn_%sph%s", richNm.c_str(), richPt.c_str()))->fill(nrig, richph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_%sph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(irig, richph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_%sph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(nrig, richph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_%sph%s%03d", richNm.c_str(), richPt.c_str(), tvBin))->fill(nrig, richph, weight);

		if ((!fRich->status) && !(isElph && isPhTh)) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 5, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 5, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 5, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 5, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 5, weight);

		//Bool_t      richCut = (fRich->status) ? (std::fabs(richbPr) > 1.5 * (richbTh + richbSh) || richph > 25) :
		//                                        ((fRich->kindOfRad == 0) ? !(MgntNum::Compare(richelph, 2.0f) > 0 && richph <= aglPhTh) : 
		//																				                           !(MgntNum::Compare(richelph, 0.5f) > 0 && richph <= nafPhTh));

		// TOF
		Bool_t isSignal     = (sign > 0 && (isrichNo || (fRich->kindOfRad == 0 && MgntNum::Compare(richprph, 1.0f) > 0 && isrichPr)));
		Bool_t isBackground = (sign < 0 && (fRich->kindOfRad == 0 && (isrichPi || isrichEl)));
		if (isSignal)     Hist::Head("hLs_Tofb")->fill(nrig, tofb, weight);
		if (isBackground) Hist::Head("hLb_Tofb")->fill(nrig, tofb, weight);
		if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hLs_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hLb_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		
		if (isSignal)     Hist::Head("hLs_TofM")->fill(nrig, tofM, weight);
		if (isBackground) Hist::Head("hLb_TofM")->fill(nrig, tofM, weight);
		if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hLs_TofM%03d", tvBin))->fill(nrig, tofM, weight);
		if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hLb_TofM%03d", tvBin))->fill(nrig, tofM, weight);
		
		Bool_t richCut = (fRich->status && std::fabs(richbPr) > 1.5 * (richbTh + richbSh));
		if (richCut) break;

		Hist::Head("hLi_Tofb")->fill(irig, tofb, weight);
		if (sign>0) Hist::Head("hLp_Tofb")->fill(nrig, tofb, weight);
		if (sign<0) Hist::Head("hLn_Tofb")->fill(nrig, tofb, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Tofb%03d", tvBin))->fill(irig, tofb, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		
		Hist::Head("hLi_TofM")->fill(irig, tofM, weight);
		if (sign>0) Hist::Head("hLp_TofM")->fill(nrig, tofM, weight);
		if (sign<0) Hist::Head("hLn_TofM")->fill(nrig, tofM, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_TofM%03d", tvBin))->fill(irig, tofM, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_TofM%03d", tvBin))->fill(nrig, tofM, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_TofM%03d", tvBin))->fill(nrig, tofM, weight);

		Hist::Head("hLi_Evt")->fill(irig, weight);
		if (sign>0) Hist::Head("hLp_Evt")->fill(nrig, weight);
		if (sign<0) Hist::Head("hLn_Evt")->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Evt%03d", tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Evt%03d", tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Evt%03d", tvBin))->fill(nrig, weight);
		
		Hist::Head("hLi_Cutflow")->fill(irig, 6, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 6, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 6, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 6, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 6, weight);

		if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
			Hist::Head("hLi_MCEvt")->fill(MCIRig, weight);
			if (MCSign>0) Hist::Head("hLp_MCEvt")->fill(MCNRig, weight);
			if (MCSign<0) Hist::Head("hLn_MCEvt")->fill(MCNRig, weight);
		}

		break;
	}


#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] I Selection [][]\n";
#endif
	/**** Intermedia Energy Region ****/
	while (true) {
		// TRK
		Hist::Head("hIi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);

		// Charge Confusion
		Hist::Head("hIi_CCest")->fill(irig, ccest, weight);
		if (sign>0) Hist::Head("hIp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hIn_CCest")->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_CCest%03d", tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		
		if (ccest > ccestCut) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		// RICH
		if (fRich->kindOfRad == -1) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);

		// Template Fitting (3-Method Merge)
		const Float_t gapReg[2] = { 5.90, 10.10 };
		//const Float_t gapReg[2] = { 6.00, 9.62 };
		Int_t richPt = 0;
		if      (nrig  < gapReg[0]) richPt = 1; // (Low        Part)
		else if (nrig >= gapReg[1]) richPt = 3; // (Intermedia Part)
		else                        richPt = 2; // (High       Part)

		Bool_t isPassedCut = false;

		while (richPt == 1) {
			if (fRich->kindOfRad != 0) break;
			if (hasEcal && !isecalPr) break;

			Bool_t isSignal     = (sign > 0 && (isrichNo || (fRich->kindOfRad == 0 && MgntNum::Compare(richprph, 1.0f) > 0 && isrichPr)));
			Bool_t isBackground = (sign < 0 && (fRich->kindOfRad == 0 && (isrichEl || isrichPi)));
			if (isSignal)     Hist::Head("hIs_Trdl")->fill(nrig, trdl, weight);
			if (isBackground) Hist::Head("hIb_Trdl")->fill(nrig, trdl, weight);
			if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			
			if (hasRich && !isrichPr) break;

			Hist::Head("hIi_Trdl")->fill(irig, trdl, weight);
			if (sign>0) Hist::Head("hIp_Trdl")->fill(nrig, trdl, weight);
			if (sign<0) Hist::Head("hIn_Trdl")->fill(nrig, trdl, weight);
			if (tvBin!=0) Hist::Head(StrFmt("hIi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			
			isPassedCut = true;
			break;
		}
		
		while (richPt == 2) {
			if (!hasRich || fRich->kindOfRad != 0) break;

			Bool_t isSignal     = (sign > 0 && isecalPr);
			Bool_t isBackground = (sign < 0 && isecalEl);
			if (isSignal)     Hist::Head("hIs_Trdl")->fill(nrig, trdl, weight);
			if (isBackground) Hist::Head("hIb_Trdl")->fill(nrig, trdl, weight);
			if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			
			if (!isrichPr) break;
			if (hasEcal && !isecalPr) break;
		
			Hist::Head("hIi_Trdl")->fill(irig, trdl, weight);
			if (sign>0) Hist::Head("hIp_Trdl")->fill(nrig, trdl, weight);
			if (sign<0) Hist::Head("hIn_Trdl")->fill(nrig, trdl, weight);
			if (tvBin!=0) Hist::Head(StrFmt("hIi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);

			isPassedCut = true;
			break;
		}	
		
		while (richPt == 3) {
			if (hasRich && !isrichPr) break;
			
			Bool_t isSignal     = (sign > 0 && isecalPr);
			Bool_t isBackground = (sign < 0 && isecalEl);
			if (isSignal)     Hist::Head("hIs_Trdl")->fill(nrig, trdl, weight);
			if (isBackground) Hist::Head("hIb_Trdl")->fill(nrig, trdl, weight);
			if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_Trdl%03d", tvBin))->fill(nrig, trdl, weight);

			if (hasEcal && !isecalPr) break;
		
			Hist::Head("hIi_Trdl")->fill(irig, trdl, weight);
			if (sign>0) Hist::Head("hIp_Trdl")->fill(nrig, trdl, weight);
			if (sign<0) Hist::Head("hIn_Trdl")->fill(nrig, trdl, weight);
			if (tvBin!=0) Hist::Head(StrFmt("hIi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);

			isPassedCut = true;
			break;
		}

		// Version L  // low    ( < 5.37)
		while (fRich->kindOfRad == 0) {
			if (hasEcal && !isecalPr) break;

			Bool_t isSignal     = (sign > 0 && (fRich->kindOfRad == 0 && MgntNum::Compare(richprph, 1.0f) > 0 && isrichPr));
			Bool_t isBackground = (sign < 0 && (fRich->kindOfRad == 0 && (isrichEl || isrichPi)));
			if (isSignal)     Hist::Head("hIs_TrdlL")->fill(nrig, trdl, weight);
			if (isBackground) Hist::Head("hIb_TrdlL")->fill(nrig, trdl, weight);
			if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_TrdlL%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_TrdlL%03d", tvBin))->fill(nrig, trdl, weight);
			
			if (hasRich && !isrichPr) break;

			Hist::Head("hIi_TrdlL")->fill(irig, trdl, weight);
			if (sign>0) Hist::Head("hIp_TrdlL")->fill(nrig, trdl, weight);
			if (sign<0) Hist::Head("hIn_TrdlL")->fill(nrig, trdl, weight);
			if (tvBin!=0) Hist::Head(StrFmt("hIi_TrdlL%03d", tvBin))->fill(irig, trdl, weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_TrdlL%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_TrdlL%03d", tvBin))->fill(nrig, trdl, weight);

			break;
		}
		
		// Version I   // all ( > 5.37    < 9.26)
		while (hasRich && fRich->kindOfRad == 0) {
			Bool_t isSignal     = (sign > 0 && isecalPr);
			Bool_t isBackground = (sign < 0 && isecalEl);
			if (isSignal)     Hist::Head("hIs_TrdlI")->fill(nrig, trdl, weight);
			if (isBackground) Hist::Head("hIb_TrdlI")->fill(nrig, trdl, weight);
			if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_TrdlI%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_TrdlI%03d", tvBin))->fill(nrig, trdl, weight);
			
			if (!isrichPr) break;
			if (hasEcal && !isecalPr) break;
		
			Hist::Head("hIi_TrdlI")->fill(irig, trdl, weight);
			if (sign>0) Hist::Head("hIp_TrdlI")->fill(nrig, trdl, weight);
			if (sign<0) Hist::Head("hIn_TrdlI")->fill(nrig, trdl, weight);
			if (tvBin!=0) Hist::Head(StrFmt("hIi_TrdlI%03d", tvBin))->fill(irig, trdl, weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_TrdlI%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_TrdlI%03d", tvBin))->fill(nrig, trdl, weight);

			break;
		}
		
		// Version H // high  ( > 9.26)
		while (true) {
			if (hasRich && !isrichPr) break;
			
			Bool_t isSignal     = (sign > 0 && isecalPr);
			Bool_t isBackground = (sign < 0 && isecalEl);
			if (isSignal)     Hist::Head("hIs_TrdlH")->fill(nrig, trdl, weight);
			if (isBackground) Hist::Head("hIb_TrdlH")->fill(nrig, trdl, weight);
			if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_TrdlH%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_TrdlH%03d", tvBin))->fill(nrig, trdl, weight);

			if (hasEcal && !isecalPr) break;
		
			Hist::Head("hIi_TrdlH")->fill(irig, trdl, weight);
			if (sign>0) Hist::Head("hIp_TrdlH")->fill(nrig, trdl, weight);
			if (sign<0) Hist::Head("hIn_TrdlH")->fill(nrig, trdl, weight);
			if (tvBin!=0) Hist::Head(StrFmt("hIi_TrdlH%03d", tvBin))->fill(irig, trdl, weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_TrdlH%03d", tvBin))->fill(nrig, trdl, weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_TrdlH%03d", tvBin))->fill(nrig, trdl, weight);

			break;
		}
			
		if (!isPassedCut) break;
		
		Hist::Head("hIi_Evt")->fill(irig, weight);
		if (sign>0) Hist::Head("hIp_Evt")->fill(nrig, weight);
		if (sign<0) Hist::Head("hIn_Evt")->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Evt%03d", tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Evt%03d", tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Evt%03d", tvBin))->fill(nrig, weight);
		
		Hist::Head("hIi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);

		if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
			Hist::Head("hIi_MCEvt")->fill(MCIRig, weight);
			if (MCSign>0) Hist::Head("hIp_MCEvt")->fill(MCNRig, weight);
			if (MCSign<0) Hist::Head("hIn_MCEvt")->fill(MCNRig, weight);
		}

		break;
	}

#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] H Selection [][]\n";
#endif
	/**** High Energy Region ****/
	while (true) {
		// TRK
		if (trPt == 2) break; 
		
		Hist::Head("hHi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		
		// ECAL	
		if (hasEcal && !isecalPr) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		// TRD
		if (!istrdPr) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		
		// Charge Confusion ---- CC Estimator
		Hist::Head(StrFmt("hHi_CCest%s", trNm.c_str()))->fill(irig, ccest, weight);
		if (sign>0) Hist::Head(StrFmt("hHp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_CCest%s%03d", trNm.c_str(), tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_CCest%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_CCest%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		
		// testcode
		for (Int_t it = 0; it <= 20; ++it) {
			Float_t rat = (1.0 - it * 0.02); 
			//Float_t lchiyCutv = lchiyCut98 * ((0.5 * (std::erf((std::log(GetStdMDR(trPt)) - std::log(nrig)) / 1.5) + 1.0)) * (1.-rat) + rat);
			//Float_t lasymCutv = lasymCut98 * ((0.5 * (std::erf((std::log(GetStdMDR(trPt)) - std::log(nrig)) / 1.5) + 1.0)) * (1.-rat) + rat);
			Float_t lchiyCutv = lchiyCut98 * rat;
			Float_t lasymCutv = lasymCut98 * rat;
			if (lchiy > lchiyCutv || lasym > lasymCutv) continue;
			//if (lchiy > 0 && lasym > 0 && ((lchiy/lchiyCutv)*(lchiy/lchiyCutv) + (lasym/lasymCutv)*(lasym/lasymCutv)) > 1) continue;
			Hist::Head(StrFmt("hHi_CCest%03d%s", it, trNm.c_str()))->fill(irig, ccest, weight);
			if (sign>0) Hist::Head(StrFmt("hHp_CCest%03d%s", it, trNm.c_str()))->fill(nrig, ccest, weight);
			if (sign<0) Hist::Head(StrFmt("hHn_CCest%03d%s", it, trNm.c_str()))->fill(nrig, ccest, weight);
			Hist::Head(StrFmt("hHi_Evt%03d%s", it, trNm.c_str()))->fill(irig, weight);
			if (sign>0) Hist::Head(StrFmt("hHp_Evt%03d%s", it, trNm.c_str()))->fill(nrig, weight);
			if (sign<0) Hist::Head(StrFmt("hHn_Evt%03d%s", it, trNm.c_str()))->fill(nrig, weight);
			Hist::Head(StrFmt("hHi_MCEvt%03d%s", it, trNm.c_str()))->fill(MCIRig, weight);
			if (MCSign>0) Hist::Head(StrFmt("hHp_MCEvt%03d%s", it, trNm.c_str()))->fill(MCNRig, weight);
			if (MCSign<0) Hist::Head(StrFmt("hHn_MCEvt%03d%s", it, trNm.c_str()))->fill(MCNRig, weight);
		}

		Hist::Head("hHi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		
		Hist::Head("hHi_Evt")->fill(irig, weight);
		if (sign>0) Hist::Head("hHp_Evt")->fill(nrig, weight);
		if (sign<0) Hist::Head("hHn_Evt")->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Evt%03d", tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Evt%03d", tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Evt%03d", tvBin))->fill(nrig, weight);
		
		Hist::Head(StrFmt("hHi_Evt%s", trNm.c_str()))->fill(irig, weight);
		if (sign>0) Hist::Head(StrFmt("hHp_Evt%s", trNm.c_str()))->fill(nrig, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_Evt%s", trNm.c_str()))->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Evt%s%03d", trNm.c_str(), tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Evt%s%03d", trNm.c_str(), tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Evt%s%03d", trNm.c_str(), tvBin))->fill(nrig, weight);
		
		if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
			Hist::Head("hHi_MCEvt")->fill(MCIRig, weight);
			if (MCSign>0) Hist::Head("hHp_MCEvt")->fill(MCNRig, weight);
			if (MCSign<0) Hist::Head("hHn_MCEvt")->fill(MCNRig, weight);
			
			Hist::Head(StrFmt("hHi_MCEvt%s", trNm.c_str()))->fill(MCIRig, weight);
			if (MCSign>0) Hist::Head(StrFmt("hHp_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
			if (MCSign<0) Hist::Head(StrFmt("hHn_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
		}
	
		break;
	}

	return true;
}
/*******************/
/**               **/
/*******************/


/*-----------------*/
/*  Main Function  */
/*-----------------*/
int main(int argc, const char ** argv) {
	std::cout << "\n**--------------------------**\n";
	std::cout << "\n**    YiAnaNtuple START    **\n";
	std::cout << "\n**--------------------------**\n";
	MgntROOT::Style::LoadDefaultEnvironment();

	std::cout << std::endl << std::endl;
	std::cout << "Usage : YiAnaNtuple event_mode file_list group_th group_size (path)\n";
	std::cout << "    Parameters : \n";
	std::cout << "    event_mode [ISS BT MC]\n";
	std::cout << "    file_list\n";
	std::cout << "    group_th\n";
	std::cout << "    group_size\n";
	std::cout << "    (path)\n";
	std::cout << std::endl << std::endl;

	if (argc != 5 && argc != 6) {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("Number of argument is not conform! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}

	std::string event_mode = argv[1];
	std::string file_list = argv[2];
	Long64_t group_th = atol(argv[3]);
	Long64_t group_size = atol(argv[4]);

	std::use_facet<std::ctype<char> >(std::locale()).toupper(&event_mode[0], &event_mode[0] + event_mode.size());
	if (event_mode == "ISS") YiNtuple::SetEventMode(YiNtuple::ISS);
	else if (event_mode == "BT") YiNtuple::SetEventMode(YiNtuple::BT);
	else if (event_mode == "MC") YiNtuple::SetEventMode(YiNtuple::MC);
	else {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("Can't find event mode (ISS, BT, MC)! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}

	std::string outputFile = StrFmt("YiAnalytics_%s.%07ld.root", event_mode.c_str(), group_th);
	std::string path = ".";
	if (argc == 6) path = argv[5];

	YiAna * ntuple = new YiAna();
	ntuple->setOutputFile(outputFile.c_str(), path.c_str());
	ntuple->readDataFrom(file_list.c_str(), group_th, group_size);
	ntuple->loopEventChain();
	if (ntuple != 0) delete ntuple;
	ntuple = 0;

	std::cout << "\n**------------------------**\n";
	std::cout << "\n**    YiAnaNtuple END    **\n";
	std::cout << "\n**------------------------**\n";
	return 0;
}
#endif // __YiAnaNtuple_C__
