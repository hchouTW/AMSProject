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
			
			// rigidity binning
			AXnr = Axis("Rigidity [GV]",
				{   1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
				    3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
					  8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
					 19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
					 41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
					 93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } );
			AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
		
			// merge 3 bins
			//AXnr = Axis("Rigidity [GV]",
			//	{   1.00,   1.51,   2.15,   2.97,   4.02, 
			//	    5.37,   7.09,   9.26,  12.00,  15.30, 
			//		 19.50,  24.70,  31.10,  38.90,  48.50, 
			//		 60.30,  74.90 } );
			//AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
			
			//AXtnr = Axis("Rigidity [GV]",
			//	{   1.00,   1.51,   2.00,   3.00,   4.02, 
			//	    5.37,   7.09,   9.26,  12.00,  15.30, 
			//		 19.50,  24.70,  31.10,  38.90,  48.50, 
			//		 60.30,  74.90 } );
			//AXtir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
			
			//AXnr = Axis("Rigidity [GV]",
			//	{   1.00,   1.46,   2.00,   3.00,   4.12, 
			//	    5.00,   6.00,   7.10,   8.30,   9.62, 
			//	   11.04,  12.59,  14.25,  16.05,  17.98, 
			//	   20.04,  22.25,  24.62,  27.25,  30.21, 
			//	   35.36,  40.00 } );
			//AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);

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

			//---- Exposure Time ----//
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Hist::New("hExpT25Deg", "", AXnr);
				Hist::New("hExpT30Deg", "", AXnr);
				Hist::New("hExpT35Deg", "", AXnr);
				Hist::New("hExpT40Deg", "", AXnr);
				for (Int_t it = 1; it <= AXtme.nbin() && gTimeStudy; ++it) {
					Hist::New(StrFmt("hExpT25Deg%03d", it), "", AXnr);
					Hist::New(StrFmt("hExpT30Deg%03d", it), "", AXnr);
					Hist::New(StrFmt("hExpT35Deg%03d", it), "", AXnr);
					Hist::New(StrFmt("hExpT40Deg%03d", it), "", AXnr);
				}
			}
		
			//---- Common Cut ----//
			Axis AXCOMcutflow("Cutflow", 13, 0., 13.);
			Hist::New("hCOM_Cutflow", "", AXCOMcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOM_Cutflow%03d", it), "", AXCOMcutflow);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Axis AXCOMliveTime("LiveTime", 100, 0.5, 1.0);
				Hist::New("hCOM_LiveTime", "", AXCOMliveTime);
				for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
					Hist::New(StrFmt("hCOM_LiveTime%03d", it), "", AXCOMliveTime);
				}
			}
			
			Axis AXCOMtrInQ("Chrg", 400, 0., 8.5);
			Hist::New("hCOMi_TrInQ", "", AXir, AXCOMtrInQ);
			Hist::New("hCOMp_TrInQ", "", AXnr, AXCOMtrInQ);
			Hist::New("hCOMn_TrInQ", "", AXnr, AXCOMtrInQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TrInQ%03d", it), "", AXir, AXCOMtrInQ);
				Hist::New(StrFmt("hCOMp_TrInQ%03d", it), "", AXnr, AXCOMtrInQ);
				Hist::New(StrFmt("hCOMn_TrInQ%03d", it), "", AXnr, AXCOMtrInQ);
			}
			
			Axis AXCOMtrExQ("Chrg", 200, 0., 3.4);
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
			
			Axis AXCOMtofQ("Chrg", 200, 0., 3.4);
			Hist::New("hCOMi_TofQ", "", AXir, AXCOMtofQ);
			Hist::New("hCOMp_TofQ", "", AXnr, AXCOMtofQ);
			Hist::New("hCOMn_TofQ", "", AXnr, AXCOMtofQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TofQ%03d", it), "", AXir, AXCOMtofQ);
				Hist::New(StrFmt("hCOMp_TofQ%03d", it), "", AXnr, AXCOMtofQ);
				Hist::New(StrFmt("hCOMn_TofQ%03d", it), "", AXnr, AXCOMtofQ);
			}
			
			Axis AXCOMdltTofQ("Delta Chrg", 200, -1.2, 1.2);
			Hist::New("hCOMi_TofQul", "", AXir, AXCOMdltTofQ);
			Hist::New("hCOMp_TofQul", "", AXnr, AXCOMdltTofQ);
			Hist::New("hCOMn_TofQul", "", AXnr, AXCOMdltTofQ);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hCOMi_TofQul%03d", it), "", AXir, AXCOMdltTofQ);
				Hist::New(StrFmt("hCOMp_TofQul%03d", it), "", AXnr, AXCOMdltTofQ);
				Hist::New(StrFmt("hCOMn_TofQul%03d", it), "", AXnr, AXCOMdltTofQ);
			}
			
			//----  Low Energy  ----//
			Axis AXLcutflow("Cutflow", 9, 0., 9.);
			Hist::New("hLi_Cutflow", "", AXir, AXLcutflow);
			Hist::New("hLp_Cutflow", "", AXnr, AXLcutflow);
			Hist::New("hLn_Cutflow", "", AXnr, AXLcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_Cutflow%03d", it), "", AXir, AXLcutflow);
				Hist::New(StrFmt("hLp_Cutflow%03d", it), "", AXnr, AXLcutflow);
				Hist::New(StrFmt("hLn_Cutflow%03d", it), "", AXnr, AXLcutflow);
			}
			
			Axis AXLlchi("Lchi", 400, -8., 8.);
			Hist::New("hLi_LchixIn", "", AXir, AXLlchi);
			Hist::New("hLp_LchixIn", "", AXnr, AXLlchi);
			Hist::New("hLn_LchixIn", "", AXnr, AXLlchi);
			Hist::New("hLi_LchiyIn", "", AXir, AXLlchi);
			Hist::New("hLp_LchiyIn", "", AXnr, AXLlchi);
			Hist::New("hLn_LchiyIn", "", AXnr, AXLlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_LchixIn%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_LchixIn%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_LchixIn%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLi_LchiyIn%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_LchiyIn%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_LchiyIn%03d", it), "", AXnr, AXLlchi);
			}
			
			Axis AXLlasym("Lasym", 400, -8., 8.);
			Hist::New("hLi_LasymUL", "", AXir, AXLlasym);
			Hist::New("hLp_LasymUL", "", AXnr, AXLlasym);
			Hist::New("hLn_LasymUL", "", AXnr, AXLlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_LasymUL%03d", it), "", AXir, AXLlasym);
				Hist::New(StrFmt("hLp_LasymUL%03d", it), "", AXnr, AXLlasym);
				Hist::New(StrFmt("hLn_LasymUL%03d", it), "", AXnr, AXLlasym);
			}
			
			Axis AXLccest("CCest", 100, -4.0, 2.5);
			Hist::New("hLi_CCestIn", "", AXir, AXLccest);
			Hist::New("hLp_CCestIn", "", AXnr, AXLccest);
			Hist::New("hLn_CCestIn", "", AXnr, AXLccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_CCestIn%03d", it), "", AXir, AXLccest);
				Hist::New(StrFmt("hLp_CCestIn%03d", it), "", AXnr, AXLccest);
				Hist::New(StrFmt("hLn_CCestIn%03d", it), "", AXnr, AXLccest);
			}
			
			Axis AXLdchrg("DChrg", 400, -0.8, 0.8);
			Hist::New("hLi_TofulChrg", "", AXir, AXLdchrg);
			Hist::New("hLp_TofulChrg", "", AXnr, AXLdchrg);
			Hist::New("hLn_TofulChrg", "", AXnr, AXLdchrg);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_TofulChrg%03d", it), "", AXir, AXLdchrg);
				Hist::New(StrFmt("hLp_TofulChrg%03d", it), "", AXnr, AXLdchrg);
				Hist::New(StrFmt("hLn_TofulChrg%03d", it), "", AXnr, AXLdchrg);
			}
			
			Axis AXLtrdl("Trdl", 400, 0.2, 1.6);
			Hist::New("hLi_Trdl", "", AXir, AXLtrdl);
			Hist::New("hLp_Trdl", "", AXnr, AXLtrdl);
			Hist::New("hLn_Trdl", "", AXnr, AXLtrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_Trdl%03d", it), "", AXir, AXLtrdl);
				Hist::New(StrFmt("hLp_Trdl%03d", it), "", AXnr, AXLtrdl);
				Hist::New(StrFmt("hLn_Trdl%03d", it), "", AXnr, AXLtrdl);
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
			
			Axis AXLrichb("Richb", 400, -0.04, 0.02);
			Hist::New("hLi_Aglb", "", AXir, AXLrichb);
			Hist::New("hLp_Aglb", "", AXnr, AXLrichb);
			Hist::New("hLn_Aglb", "", AXnr, AXLrichb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_Aglb%03d", it), "", AXir, AXLrichb);
				Hist::New(StrFmt("hLp_Aglb%03d", it), "", AXnr, AXLrichb);
				Hist::New(StrFmt("hLn_Aglb%03d", it), "", AXnr, AXLrichb);
			}
			
			Hist::New("hLi_Nafb", "", AXir, AXLrichb);
			Hist::New("hLp_Nafb", "", AXnr, AXLrichb);
			Hist::New("hLn_Nafb", "", AXnr, AXLrichb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_Nafb%03d", it), "", AXir, AXLrichb);
				Hist::New(StrFmt("hLp_Nafb%03d", it), "", AXnr, AXLrichb);
				Hist::New(StrFmt("hLn_Nafb%03d", it), "", AXnr, AXLrichb);
			}
			
			Axis AXLtofrichb("TofRichb", 400, -0.16, 0.16);
			Hist::New("hLi_TofAglb", "", AXir, AXLtofrichb);
			Hist::New("hLp_TofAglb", "", AXnr, AXLtofrichb);
			Hist::New("hLn_TofAglb", "", AXnr, AXLtofrichb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_TofAglb%03d", it), "", AXir, AXLtofrichb);
				Hist::New(StrFmt("hLp_TofAglb%03d", it), "", AXnr, AXLtofrichb);
				Hist::New(StrFmt("hLn_TofAglb%03d", it), "", AXnr, AXLtofrichb);
			}
			
			Hist::New("hLi_TofNafb", "", AXir, AXLtofrichb);
			Hist::New("hLp_TofNafb", "", AXnr, AXLtofrichb);
			Hist::New("hLn_TofNafb", "", AXnr, AXLtofrichb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_TofNafb%03d", it), "", AXir, AXLtofrichb);
				Hist::New(StrFmt("hLp_TofNafb%03d", it), "", AXnr, AXLtofrichb);
				Hist::New(StrFmt("hLn_TofNafb%03d", it), "", AXnr, AXLtofrichb);
			}
			
			Axis AXLtofb("Tofb", 50, -0.16, 0.4);
			Hist::New("hLi_Tofb", "", AXir, AXLtofb);
			Hist::New("hLp_Tofb", "", AXnr, AXLtofb);
			Hist::New("hLn_Tofb", "", AXnr, AXLtofb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLi_Tofb%03d", it), "", AXir, AXLtofb);
				Hist::New(StrFmt("hLp_Tofb%03d", it), "", AXnr, AXLtofb);
				Hist::New(StrFmt("hLn_Tofb%03d", it), "", AXnr, AXLtofb);
			}
			
			Hist::New("hLs_Tofb", "", AXnr, AXLtofb);
			Hist::New("hLb_Tofb", "", AXnr, AXLtofb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hLs_Tofb%03d", it), "", AXnr, AXLtofb);
				Hist::New(StrFmt("hLb_Tofb%03d", it), "", AXnr, AXLtofb);
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
				Hist::New("hLi_MCTrRso", "", AXir, AXir);
				Hist::New("hLp_MCTrRso", "", AXnr, AXir);
				Hist::New("hLn_MCTrRso", "", AXnr, AXir);
				
				Hist::New("hLi_MCEvt", "", AXir);
				Hist::New("hLp_MCEvt", "", AXnr);
				Hist::New("hLn_MCEvt", "", AXnr);
			}
			
			//----  Intermedia Energy  ----//
			Axis AXIcutflow("Cutflow", 9, 0., 9.);
			Hist::New("hIi_Cutflow", "", AXir, AXIcutflow);
			Hist::New("hIp_Cutflow", "", AXnr, AXIcutflow);
			Hist::New("hIn_Cutflow", "", AXnr, AXIcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Cutflow%03d", it), "", AXir, AXIcutflow);
				Hist::New(StrFmt("hIp_Cutflow%03d", it), "", AXnr, AXIcutflow);
				Hist::New(StrFmt("hIn_Cutflow%03d", it), "", AXnr, AXIcutflow);
			}
			
			Axis AXIaglb("Aglb", 400, -0.04, 0.02);
			Hist::New("hIi_Aglb", "", AXir, AXIaglb);
			Hist::New("hIp_Aglb", "", AXnr, AXIaglb);
			Hist::New("hIn_Aglb", "", AXnr, AXIaglb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Aglb%03d", it), "", AXir, AXIaglb);
				Hist::New(StrFmt("hIp_Aglb%03d", it), "", AXnr, AXIaglb);
				Hist::New(StrFmt("hIn_Aglb%03d", it), "", AXnr, AXIaglb);
			}
			
			Axis AXIecalbdt("Ecalbdt", 400, -1.0, 1.0);
			Hist::New("hIi_Ecalbdt", "", AXir, AXIecalbdt);
			Hist::New("hIp_Ecalbdt", "", AXnr, AXIecalbdt);
			Hist::New("hIn_Ecalbdt", "", AXnr, AXIecalbdt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Ecalbdt%03d", it), "", AXir, AXIecalbdt);
				Hist::New(StrFmt("hIp_Ecalbdt%03d", it), "", AXnr, AXIecalbdt);
				Hist::New(StrFmt("hIn_Ecalbdt%03d", it), "", AXnr, AXIecalbdt);
			}

			Axis AXIlasym("Lasym", 400, -8., 8.);
			Hist::New("hIi_LasymUL", "", AXir, AXIlasym);
			Hist::New("hIp_LasymUL", "", AXnr, AXIlasym);
			Hist::New("hIn_LasymUL", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_LasymUL%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_LasymUL%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_LasymUL%03d", it), "", AXnr, AXIlasym);
			}
			
			Hist::New("hIi_Lasym1I", "", AXir, AXIlasym);
			Hist::New("hIp_Lasym1I", "", AXnr, AXIlasym);
			Hist::New("hIn_Lasym1I", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Lasym1I%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_Lasym1I%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_Lasym1I%03d", it), "", AXnr, AXIlasym);
			}
			
			Hist::New("hIi_Lasym9I", "", AXir, AXIlasym);
			Hist::New("hIp_Lasym9I", "", AXnr, AXIlasym);
			Hist::New("hIn_Lasym9I", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Lasym9I%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_Lasym9I%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_Lasym9I%03d", it), "", AXnr, AXIlasym);
			}
			
			Hist::New("hIi_Lasym91", "", AXir, AXIlasym);
			Hist::New("hIp_Lasym91", "", AXnr, AXIlasym);
			Hist::New("hIn_Lasym91", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Lasym91%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_Lasym91%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_Lasym91%03d", it), "", AXnr, AXIlasym);
			}
			
			Axis AXIlchi("Lchi", 400, -8., 8.);
			Hist::New("hIi_LchixIn", "", AXir, AXIlchi);
			Hist::New("hIp_LchixIn", "", AXnr, AXIlchi);
			Hist::New("hIn_LchixIn", "", AXnr, AXIlchi);
			Hist::New("hIi_LchiyIn", "", AXir, AXIlchi);
			Hist::New("hIp_LchiyIn", "", AXnr, AXIlchi);
			Hist::New("hIn_LchiyIn", "", AXnr, AXIlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_LchixIn%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchixIn%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchixIn%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIi_LchiyIn%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchiyIn%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchiyIn%03d", it), "", AXnr, AXIlchi);
			}
			
			Hist::New("hIi_LchixL1", "", AXir, AXIlchi);
			Hist::New("hIp_LchixL1", "", AXnr, AXIlchi);
			Hist::New("hIn_LchixL1", "", AXnr, AXIlchi);
			Hist::New("hIi_LchiyL1", "", AXir, AXIlchi);
			Hist::New("hIp_LchiyL1", "", AXnr, AXIlchi);
			Hist::New("hIn_LchiyL1", "", AXnr, AXIlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_LchixL1%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchixL1%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchixL1%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIi_LchiyL1%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchiyL1%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchiyL1%03d", it), "", AXnr, AXIlchi);
			}
			
			Hist::New("hIi_LchixL9", "", AXir, AXIlchi);
			Hist::New("hIp_LchixL9", "", AXnr, AXIlchi);
			Hist::New("hIn_LchixL9", "", AXnr, AXIlchi);
			Hist::New("hIi_LchiyL9", "", AXir, AXIlchi);
			Hist::New("hIp_LchiyL9", "", AXnr, AXIlchi);
			Hist::New("hIn_LchiyL9", "", AXnr, AXIlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_LchixL9%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchixL9%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchixL9%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIi_LchiyL9%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchiyL9%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchiyL9%03d", it), "", AXnr, AXIlchi);
			}
			
			Hist::New("hIi_LchixFs", "", AXir, AXIlchi);
			Hist::New("hIp_LchixFs", "", AXnr, AXIlchi);
			Hist::New("hIn_LchixFs", "", AXnr, AXIlchi);
			Hist::New("hIi_LchiyFs", "", AXir, AXIlchi);
			Hist::New("hIp_LchiyFs", "", AXnr, AXIlchi);
			Hist::New("hIn_LchiyFs", "", AXnr, AXIlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_LchixFs%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchixFs%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchixFs%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIi_LchiyFs%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchiyFs%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchiyFs%03d", it), "", AXnr, AXIlchi);
			}
			
			Axis AXIccest("CCest", 100, -4.0, 2.0);
			Hist::New("hIi_CCestIn", "", AXir, AXIccest);
			Hist::New("hIp_CCestIn", "", AXnr, AXIccest);
			Hist::New("hIn_CCestIn", "", AXnr, AXIccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_CCestIn%03d", it), "", AXir, AXIccest);
				Hist::New(StrFmt("hIp_CCestIn%03d", it), "", AXnr, AXIccest);
				Hist::New(StrFmt("hIn_CCestIn%03d", it), "", AXnr, AXIccest);
			}
			
			Hist::New("hIi_CCestL1", "", AXir, AXIccest);
			Hist::New("hIp_CCestL1", "", AXnr, AXIccest);
			Hist::New("hIn_CCestL1", "", AXnr, AXIccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_CCestL1%03d", it), "", AXir, AXIccest);
				Hist::New(StrFmt("hIp_CCestL1%03d", it), "", AXnr, AXIccest);
				Hist::New(StrFmt("hIn_CCestL1%03d", it), "", AXnr, AXIccest);
			}
			
			Hist::New("hIi_CCestL9", "", AXir, AXIccest);
			Hist::New("hIp_CCestL9", "", AXnr, AXIccest);
			Hist::New("hIn_CCestL9", "", AXnr, AXIccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_CCestL9%03d", it), "", AXir, AXIccest);
				Hist::New(StrFmt("hIp_CCestL9%03d", it), "", AXnr, AXIccest);
				Hist::New(StrFmt("hIn_CCestL9%03d", it), "", AXnr, AXIccest);
			}
			
			Hist::New("hIi_CCestFs", "", AXir, AXIccest);
			Hist::New("hIp_CCestFs", "", AXnr, AXIccest);
			Hist::New("hIn_CCestFs", "", AXnr, AXIccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_CCestFs%03d", it), "", AXir, AXIccest);
				Hist::New(StrFmt("hIp_CCestFs%03d", it), "", AXnr, AXIccest);
				Hist::New(StrFmt("hIn_CCestFs%03d", it), "", AXnr, AXIccest);
			}
			
			Axis AXItrdl("Trdl", 50, 0.2, 1.6);
			Hist::New("hIs_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIb_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_TrdlIn", "", AXnr, AXItrdl);
			Hist::New("hIb_TrdlIn", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_TrdlIn%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_TrdlIn%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_TrdlL1", "", AXnr, AXItrdl);
			Hist::New("hIb_TrdlL1", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_TrdlL1%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_TrdlL1%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_TrdlL9", "", AXnr, AXItrdl);
			Hist::New("hIb_TrdlL9", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_TrdlL9%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_TrdlL9%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_TrdlFs", "", AXnr, AXItrdl);
			Hist::New("hIb_TrdlFs", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_TrdlFs%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_TrdlFs%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_TrdlSp", "", AXnr, AXItrdl);
			Hist::New("hIb_TrdlSp", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIs_TrdlSp%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_TrdlSp%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_Trdl", "", AXir, AXItrdl);
			Hist::New("hIp_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIn_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Trdl%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_TrdlIn", "", AXir, AXItrdl);
			Hist::New("hIp_TrdlIn", "", AXnr, AXItrdl);
			Hist::New("hIn_TrdlIn", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_TrdlIn%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_TrdlIn%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_TrdlIn%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_TrdlL1", "", AXir, AXItrdl);
			Hist::New("hIp_TrdlL1", "", AXnr, AXItrdl);
			Hist::New("hIn_TrdlL1", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_TrdlL1%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_TrdlL1%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_TrdlL1%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_TrdlL9", "", AXir, AXItrdl);
			Hist::New("hIp_TrdlL9", "", AXnr, AXItrdl);
			Hist::New("hIn_TrdlL9", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_TrdlL9%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_TrdlL9%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_TrdlL9%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_TrdlFs", "", AXir, AXItrdl);
			Hist::New("hIp_TrdlFs", "", AXnr, AXItrdl);
			Hist::New("hIn_TrdlFs", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_TrdlFs%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_TrdlFs%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_TrdlFs%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_TrdlSp", "", AXir, AXItrdl);
			Hist::New("hIp_TrdlSp", "", AXnr, AXItrdl);
			Hist::New("hIn_TrdlSp", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_TrdlSp%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_TrdlSp%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_TrdlSp%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_Evt", "", AXir);
			Hist::New("hIp_Evt", "", AXnr);
			Hist::New("hIn_Evt", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_Evt%03d", it), "", AXir);
				Hist::New(StrFmt("hIp_Evt%03d", it), "", AXnr);
				Hist::New(StrFmt("hIn_Evt%03d", it), "", AXnr);
			}
			
			Hist::New("hIi_EvtIn", "", AXir);
			Hist::New("hIp_EvtIn", "", AXnr);
			Hist::New("hIn_EvtIn", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_EvtIn%03d", it), "", AXir);
				Hist::New(StrFmt("hIp_EvtIn%03d", it), "", AXnr);
				Hist::New(StrFmt("hIn_EvtIn%03d", it), "", AXnr);
			}
			
			Hist::New("hIi_EvtL1", "", AXir);
			Hist::New("hIp_EvtL1", "", AXnr);
			Hist::New("hIn_EvtL1", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_EvtL1%03d", it), "", AXir);
				Hist::New(StrFmt("hIp_EvtL1%03d", it), "", AXnr);
				Hist::New(StrFmt("hIn_EvtL1%03d", it), "", AXnr);
			}
			
			Hist::New("hIi_EvtL9", "", AXir);
			Hist::New("hIp_EvtL9", "", AXnr);
			Hist::New("hIn_EvtL9", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_EvtL9%03d", it), "", AXir);
				Hist::New(StrFmt("hIp_EvtL9%03d", it), "", AXnr);
				Hist::New(StrFmt("hIn_EvtL9%03d", it), "", AXnr);
			}
			
			Hist::New("hIi_EvtFs", "", AXir);
			Hist::New("hIp_EvtFs", "", AXnr);
			Hist::New("hIn_EvtFs", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_EvtFs%03d", it), "", AXir);
				Hist::New(StrFmt("hIp_EvtFs%03d", it), "", AXnr);
				Hist::New(StrFmt("hIn_EvtFs%03d", it), "", AXnr);
			}
			
			Hist::New("hIi_EvtSp", "", AXir);
			Hist::New("hIp_EvtSp", "", AXnr);
			Hist::New("hIn_EvtSp", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hIi_EvtSp%03d", it), "", AXir);
				Hist::New(StrFmt("hIp_EvtSp%03d", it), "", AXnr);
				Hist::New(StrFmt("hIn_EvtSp%03d", it), "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hIi_MCTrRso", "", AXir, AXir);
				Hist::New("hIp_MCTrRso", "", AXnr, AXir);
				Hist::New("hIn_MCTrRso", "", AXnr, AXir);
				
				Hist::New("hIi_MCTrRsoIn", "", AXir, AXir);
				Hist::New("hIp_MCTrRsoIn", "", AXnr, AXir);
				Hist::New("hIn_MCTrRsoIn", "", AXnr, AXir);
				
				Hist::New("hIi_MCTrRsoL1", "", AXir, AXir);
				Hist::New("hIp_MCTrRsoL1", "", AXnr, AXir);
				Hist::New("hIn_MCTrRsoL1", "", AXnr, AXir);
				
				Hist::New("hIi_MCTrRsoL9", "", AXir, AXir);
				Hist::New("hIp_MCTrRsoL9", "", AXnr, AXir);
				Hist::New("hIn_MCTrRsoL9", "", AXnr, AXir);
				
				Hist::New("hIi_MCTrRsoFs", "", AXir, AXir);
				Hist::New("hIp_MCTrRsoFs", "", AXnr, AXir);
				Hist::New("hIn_MCTrRsoFs", "", AXnr, AXir);
				
				Hist::New("hIi_MCTrRsoSp", "", AXir, AXir);
				Hist::New("hIp_MCTrRsoSp", "", AXnr, AXir);
				Hist::New("hIn_MCTrRsoSp", "", AXnr, AXir);
				
				Hist::New("hIi_MCEvt", "", AXir);
				Hist::New("hIp_MCEvt", "", AXnr);
				Hist::New("hIn_MCEvt", "", AXnr);
				
				Hist::New("hIi_MCEvtIn", "", AXir);
				Hist::New("hIp_MCEvtIn", "", AXnr);
				Hist::New("hIn_MCEvtIn", "", AXnr);
				
				Hist::New("hIi_MCEvtL1", "", AXir);
				Hist::New("hIp_MCEvtL1", "", AXnr);
				Hist::New("hIn_MCEvtL1", "", AXnr);
				
				Hist::New("hIi_MCEvtL9", "", AXir);
				Hist::New("hIp_MCEvtL9", "", AXnr);
				Hist::New("hIn_MCEvtL9", "", AXnr);
				
				Hist::New("hIi_MCEvtFs", "", AXir);
				Hist::New("hIp_MCEvtFs", "", AXnr);
				Hist::New("hIn_MCEvtFs", "", AXnr);
				
				Hist::New("hIi_MCEvtSp", "", AXir);
				Hist::New("hIp_MCEvtSp", "", AXnr);
				Hist::New("hIn_MCEvtSp", "", AXnr);
			}
			
			//----  High Energy  ----//
			Axis AXHcutflow("Cutflow", 7, 0., 7.);
			Hist::New("hHi_Cutflow", "", AXir, AXHcutflow);
			Hist::New("hHp_Cutflow", "", AXnr, AXHcutflow);
			Hist::New("hHn_Cutflow", "", AXnr, AXHcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_Cutflow%03d", it), "", AXir, AXHcutflow);
				Hist::New(StrFmt("hHp_Cutflow%03d", it), "", AXnr, AXHcutflow);
				Hist::New(StrFmt("hHn_Cutflow%03d", it), "", AXnr, AXHcutflow);
			}
			
			Axis AXHtrdl("Trdl", 400, 0.2, 1.6);
			Hist::New("hHi_Trdl", "", AXir, AXHtrdl);
			Hist::New("hHp_Trdl", "", AXnr, AXHtrdl);
			Hist::New("hHn_Trdl", "", AXnr, AXHtrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_Trdl%03d", it), "", AXir, AXHtrdl);
				Hist::New(StrFmt("hHp_Trdl%03d", it), "", AXnr, AXHtrdl);
				Hist::New(StrFmt("hHn_Trdl%03d", it), "", AXnr, AXHtrdl);
			}
			
			Axis AXHecalbdt("Ecalbdt", 400, -1.0, 1.0);
			Hist::New("hHi_Ecalbdt", "", AXir, AXHecalbdt);
			Hist::New("hHp_Ecalbdt", "", AXnr, AXHecalbdt);
			Hist::New("hHn_Ecalbdt", "", AXnr, AXHecalbdt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_Ecalbdt%03d", it), "", AXir, AXHecalbdt);
				Hist::New(StrFmt("hHp_Ecalbdt%03d", it), "", AXnr, AXHecalbdt);
				Hist::New(StrFmt("hHn_Ecalbdt%03d", it), "", AXnr, AXHecalbdt);
			}
			
			Axis AXHlasym("Lasym", 400, -6., 6.);
			Hist::New("hHi_Lasym1I", "", AXir, AXHlasym);
			Hist::New("hHp_Lasym1I", "", AXnr, AXHlasym);
			Hist::New("hHn_Lasym1I", "", AXnr, AXHlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_Lasym1I%03d", it), "", AXir, AXHlasym);
				Hist::New(StrFmt("hHp_Lasym1I%03d", it), "", AXnr, AXHlasym);
				Hist::New(StrFmt("hHn_Lasym1I%03d", it), "", AXnr, AXHlasym);
			}
			
			Hist::New("hHi_Lasym9I", "", AXir, AXHlasym);
			Hist::New("hHp_Lasym9I", "", AXnr, AXHlasym);
			Hist::New("hHn_Lasym9I", "", AXnr, AXHlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_Lasym9I%03d", it), "", AXir, AXHlasym);
				Hist::New(StrFmt("hHp_Lasym9I%03d", it), "", AXnr, AXHlasym);
				Hist::New(StrFmt("hHn_Lasym9I%03d", it), "", AXnr, AXHlasym);
			}
			
			Hist::New("hHi_Lasym91", "", AXir, AXHlasym);
			Hist::New("hHp_Lasym91", "", AXnr, AXHlasym);
			Hist::New("hHn_Lasym91", "", AXnr, AXHlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_Lasym91%03d", it), "", AXir, AXHlasym);
				Hist::New(StrFmt("hHp_Lasym91%03d", it), "", AXnr, AXHlasym);
				Hist::New(StrFmt("hHn_Lasym91%03d", it), "", AXnr, AXHlasym);
			}
			
			Axis AXHlchi("Lchi", 400, -6., 6.);
			Hist::New("hHi_LchixL1", "", AXir, AXHlchi);
			Hist::New("hHp_LchixL1", "", AXnr, AXHlchi);
			Hist::New("hHn_LchixL1", "", AXnr, AXHlchi);
			Hist::New("hHi_LchiyL1", "", AXir, AXHlchi);
			Hist::New("hHp_LchiyL1", "", AXnr, AXHlchi);
			Hist::New("hHn_LchiyL1", "", AXnr, AXHlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_LchixL1%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_LchixL1%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_LchixL1%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHi_LchiyL1%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_LchiyL1%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_LchiyL1%03d", it), "", AXnr, AXHlchi);
			}
			
			Hist::New("hHi_LchixL9", "", AXir, AXHlchi);
			Hist::New("hHp_LchixL9", "", AXnr, AXHlchi);
			Hist::New("hHn_LchixL9", "", AXnr, AXHlchi);
			Hist::New("hHi_LchiyL9", "", AXir, AXHlchi);
			Hist::New("hHp_LchiyL9", "", AXnr, AXHlchi);
			Hist::New("hHn_LchiyL9", "", AXnr, AXHlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_LchixL9%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_LchixL9%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_LchixL9%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHi_LchiyL9%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_LchiyL9%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_LchiyL9%03d", it), "", AXnr, AXHlchi);
			}
			
			Hist::New("hHi_LchixFs", "", AXir, AXHlchi);
			Hist::New("hHp_LchixFs", "", AXnr, AXHlchi);
			Hist::New("hHn_LchixFs", "", AXnr, AXHlchi);
			Hist::New("hHi_LchiyFs", "", AXir, AXHlchi);
			Hist::New("hHp_LchiyFs", "", AXnr, AXHlchi);
			Hist::New("hHn_LchiyFs", "", AXnr, AXHlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_LchixFs%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_LchixFs%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_LchixFs%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHi_LchiyFs%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_LchiyFs%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_LchiyFs%03d", it), "", AXnr, AXHlchi);
			}
			
			Axis AXHccest("CCest", 40, -4.0, 2.5);
			Hist::New("hHi_CCest", "", AXir, AXHccest);
			Hist::New("hHp_CCest", "", AXnr, AXHccest);
			Hist::New("hHn_CCest", "", AXnr, AXHccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_CCest%03d", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_CCest%03d", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_CCest%03d", it), "", AXnr, AXHccest);
			}
			
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
			
			Axis AXHccestq("CCestQ", 100, -5.0, 4.0);
			Hist::New("hHi_CCestQ", "", AXir, AXHccestq);
			Hist::New("hHp_CCestQ", "", AXnr, AXHccestq);
			Hist::New("hHn_CCestQ", "", AXnr, AXHccestq);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_CCestQ%03d", it), "", AXir, AXHccestq);
				Hist::New(StrFmt("hHp_CCestQ%03d", it), "", AXnr, AXHccestq);
				Hist::New(StrFmt("hHn_CCestQ%03d", it), "", AXnr, AXHccestq);
			}
			
			Hist::New("hHi_CCestQL1", "", AXir, AXHccestq);
			Hist::New("hHp_CCestQL1", "", AXnr, AXHccestq);
			Hist::New("hHn_CCestQL1", "", AXnr, AXHccestq);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_CCestQL1%03d", it), "", AXir, AXHccestq);
				Hist::New(StrFmt("hHp_CCestQL1%03d", it), "", AXnr, AXHccestq);
				Hist::New(StrFmt("hHn_CCestQL1%03d", it), "", AXnr, AXHccestq);
			}
			
			Hist::New("hHi_CCestQL9", "", AXir, AXHccestq);
			Hist::New("hHp_CCestQL9", "", AXnr, AXHccestq);
			Hist::New("hHn_CCestQL9", "", AXnr, AXHccestq);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_CCestQL9%03d", it), "", AXir, AXHccestq);
				Hist::New(StrFmt("hHp_CCestQL9%03d", it), "", AXnr, AXHccestq);
				Hist::New(StrFmt("hHn_CCestQL9%03d", it), "", AXnr, AXHccestq);
			}
			
			Hist::New("hHi_CCestQFs", "", AXir, AXHccestq);
			Hist::New("hHp_CCestQFs", "", AXnr, AXHccestq);
			Hist::New("hHn_CCestQFs", "", AXnr, AXHccestq);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_CCestQFs%03d", it), "", AXir, AXHccestq);
				Hist::New(StrFmt("hHp_CCestQFs%03d", it), "", AXnr, AXHccestq);
				Hist::New(StrFmt("hHn_CCestQFs%03d", it), "", AXnr, AXHccestq);
			}
			
			Axis AXHlresx("Lres", 100., -6., 2.5);
			Axis AXHlresy("Lres", 100., -6., 2.5);

			Hist::New("hHi_L1LresxL1", "", AXir, AXHlresx);
			Hist::New("hHp_L1LresxL1", "", AXnr, AXHlresx);
			Hist::New("hHn_L1LresxL1", "", AXnr, AXHlresx);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_L1LresxL1%03d", it), "", AXir, AXHlresx);
				Hist::New(StrFmt("hHp_L1LresxL1%03d", it), "", AXnr, AXHlresx);
				Hist::New(StrFmt("hHn_L1LresxL1%03d", it), "", AXnr, AXHlresx);
			}
			
			Hist::New("hHi_L1LresyL1", "", AXir, AXHlresy);
			Hist::New("hHp_L1LresyL1", "", AXnr, AXHlresy);
			Hist::New("hHn_L1LresyL1", "", AXnr, AXHlresy);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_L1LresyL1%03d", it), "", AXir, AXHlresy);
				Hist::New(StrFmt("hHp_L1LresyL1%03d", it), "", AXnr, AXHlresy);
				Hist::New(StrFmt("hHn_L1LresyL1%03d", it), "", AXnr, AXHlresy);
			}

			Hist::New("hHi_L9LresxL9", "", AXir, AXHlresx);
			Hist::New("hHp_L9LresxL9", "", AXnr, AXHlresx);
			Hist::New("hHn_L9LresxL9", "", AXnr, AXHlresx);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_L9LresxL9%03d", it), "", AXir, AXHlresx);
				Hist::New(StrFmt("hHp_L9LresxL9%03d", it), "", AXnr, AXHlresx);
				Hist::New(StrFmt("hHn_L9LresxL9%03d", it), "", AXnr, AXHlresx);
			}
			
			Hist::New("hHi_L9LresyL9", "", AXir, AXHlresy);
			Hist::New("hHp_L9LresyL9", "", AXnr, AXHlresy);
			Hist::New("hHn_L9LresyL9", "", AXnr, AXHlresy);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_L9LresyL9%03d", it), "", AXir, AXHlresy);
				Hist::New(StrFmt("hHp_L9LresyL9%03d", it), "", AXnr, AXHlresy);
				Hist::New(StrFmt("hHn_L9LresyL9%03d", it), "", AXnr, AXHlresy);
			}
			
			Hist::New("hHi_FsLresxL1", "", AXir, AXHlresx);
			Hist::New("hHp_FsLresxL1", "", AXnr, AXHlresx);
			Hist::New("hHn_FsLresxL1", "", AXnr, AXHlresx);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_FsLresxL1%03d", it), "", AXir, AXHlresx);
				Hist::New(StrFmt("hHp_FsLresxL1%03d", it), "", AXnr, AXHlresx);
				Hist::New(StrFmt("hHn_FsLresxL1%03d", it), "", AXnr, AXHlresx);
			}
			
			Hist::New("hHi_FsLresyL1", "", AXir, AXHlresy);
			Hist::New("hHp_FsLresyL1", "", AXnr, AXHlresy);
			Hist::New("hHn_FsLresyL1", "", AXnr, AXHlresy);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_FsLresyL1%03d", it), "", AXir, AXHlresy);
				Hist::New(StrFmt("hHp_FsLresyL1%03d", it), "", AXnr, AXHlresy);
				Hist::New(StrFmt("hHn_FsLresyL1%03d", it), "", AXnr, AXHlresy);
			}
			
			Hist::New("hHi_FsLresxL9", "", AXir, AXHlresx);
			Hist::New("hHp_FsLresxL9", "", AXnr, AXHlresx);
			Hist::New("hHn_FsLresxL9", "", AXnr, AXHlresx);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_FsLresxL9%03d", it), "", AXir, AXHlresx);
				Hist::New(StrFmt("hHp_FsLresxL9%03d", it), "", AXnr, AXHlresx);
				Hist::New(StrFmt("hHn_FsLresxL9%03d", it), "", AXnr, AXHlresx);
			}
			
			Hist::New("hHi_FsLresyL9", "", AXir, AXHlresy);
			Hist::New("hHp_FsLresyL9", "", AXnr, AXHlresy);
			Hist::New("hHn_FsLresyL9", "", AXnr, AXHlresy);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS) && gTimeStudy; ++it) {
				Hist::New(StrFmt("hHi_FsLresyL9%03d", it), "", AXir, AXHlresy);
				Hist::New(StrFmt("hHp_FsLresyL9%03d", it), "", AXnr, AXHlresy);
				Hist::New(StrFmt("hHn_FsLresyL9%03d", it), "", AXnr, AXHlresy);
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
				Hist::New("hHi_MCTrRso", "", AXir, AXir);
				Hist::New("hHp_MCTrRso", "", AXnr, AXir);
				Hist::New("hHn_MCTrRso", "", AXnr, AXir);
				
				Hist::New("hHi_MCTrRsoL1", "", AXir, AXir);
				Hist::New("hHp_MCTrRsoL1", "", AXnr, AXir);
				Hist::New("hHn_MCTrRsoL1", "", AXnr, AXir);
				
				Hist::New("hHi_MCTrRsoL9", "", AXir, AXir);
				Hist::New("hHp_MCTrRsoL9", "", AXnr, AXir);
				Hist::New("hHn_MCTrRsoL9", "", AXnr, AXir);
				
				Hist::New("hHi_MCTrRsoFs", "", AXir, AXir);
				Hist::New("hHp_MCTrRsoFs", "", AXnr, AXir);
				Hist::New("hHn_MCTrRsoFs", "", AXnr, AXir);
				
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

Bool_t DST::gTimeStudy = true;


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
	Float_t SqrMDR = 2.0; // (TV)
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
		bool analyzeRunInfo(Long64_t ientry);
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
		if (!analyzeRunInfo(ientry)) continue;
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
bool YiAna::analyzeRunInfo(Long64_t ientry) {
#if Debug == true
	std::cerr << "YiAna::analyzeRunInfo()\n";
#endif

	if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
		// RTI
		if (gUTime[0] != fRti->uTime) {
			gUTime[0] = (gUTime[1] == 0) ? gUTime[1] : fRti->uTime;
			gUTime[1] = fRti->uTime;
	
			// Time Binning
			Bool_t tvStudy = gTimeStudy && (fRti->uTime > AXtme.min() && fRti->uTime < AXtme.max());
			Int_t tvBin = (!tvStudy) ? 0 : AXtme.find(fRti->uTime);

			// Exposure Time
			const Int_t nDeg = 4;
			const Int_t sDeg[nDeg] = { 25, 30, 35, 40 };
			for (Int_t iDeg = 0; iDeg < nDeg; ++iDeg) {
				std::string nmDeg = StrFmt("hExpT%02dDeg", sDeg[iDeg]);
				Float_t     cfRig = DST::CfStableFT * fRti->cutoffIGRF[iDeg];
				Int_t       nrBin = AXnr.find(cfRig) + 1;
				if (cfRig > AXnr.max()) continue;	
				if (cfRig < AXnr.min()) nrBin = 1;
				for (Int_t ibin = nrBin; ibin <= AXnr.nbin(); ++ibin) {
					Double_t cen = AXnr.center(ibin, Axis::kLog);
					Hist::Head(nmDeg)->fill(cen, fRti->liveTime);
					if (tvBin!=0) Hist::Head(StrFmt("%s%03d", nmDeg.c_str(), tvBin))->fill(cen, fRti->liveTime);
				}
			}
		} // if algo --- utime
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

	// Constant
	const Float_t massPr = 0.938272297; 
	const Float_t massPi = 0.139570180; 
	const Float_t massEl = 0.000510999; 

	// Monte Carlo
	Float_t MCTRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
	                 0.0 : (fG4mc->primPart.mom / fG4mc->primPart.chrg);
	Float_t MCNRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
	                 0.0 : std::fabs(MCTRig);
	Float_t MCIRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ?
	                 0.0 : (1. / MCTRig);
	Float_t MCrwgt = TMath::Power(MCNRig, -1.7);

	// Time Binning
	Bool_t tvStudy = gTimeStudy && YiNtuple::CheckEventMode(YiNtuple::ISS) && (fRti->uTime > AXtme.min() && fRti->uTime < AXtme.max());
	Int_t tvBin = (!tvStudy) ? 0 : AXtme.find(fRti->uTime);

#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] COM Selection [][]\n";
#endif
	
	// Sample
	Hist::Head("hCOM_Cutflow")->fill(0, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(0, weight);

	// preselection (Mini-RTI Requirement)
	const Float_t CfStableFT = 1.2;
	if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
		if (!fRti->flagRun) return false;
		if (!fRti->isGoodSecond) return false;
		if (fRti->isInSAA) return false;
		if (fRti->isInShadow) return false;
		
		Hist::Head("hCOM_LiveTime")->fill(fRti->liveTime, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hCOM_LiveTime%03d", tvBin))->fill(fRti->liveTime, weight);
		if (fRti->liveTime > 0.7) {
			Hist::Head("hCOM_Cutflow")->fill(12, weight);
			if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(12, weight);
		}
	}
	
	Hist::Head("hCOM_Cutflow")->fill(1, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(1, weight);
	
	// preselection (Mini-Track Requirement)
	TrackInfo * track = (fTrk->tracks.size() == 1) ? (&fTrk->tracks.at(0)) : nullptr;
	Int_t       refPt = -1;
	for (Int_t ipt = 5; track != nullptr && ipt >= 3; --ipt) {
		if (!track->status[0][ipt]) continue;
		refPt = ipt; break;
	}
	Float_t refTRig = (refPt == -1) ? 0.0 : track->rigidity[0][refPt];
	Float_t refNRig = (refPt == -1) ? 0.0 : std::fabs(refTRig);
	Float_t refIRig = (refPt == -1) ? 0.0 : 1. / refTRig;
	Int_t   refSign = (refPt == -1) ? 0   : MgntNum::Compare(refTRig);

	if (track == nullptr || refPt == -1) return false;
	if (MgntNum::EqualToZero(track->Qinner)) return false;
	Hist::Head("hCOM_Cutflow")->fill(2, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(2, weight);
		
	Hist::Head("hCOMi_TrInQ")->fill(refIRig, track->Qinner, weight);
	if (refSign>0) Hist::Head("hCOMp_TrInQ")->fill(refNRig, track->Qinner, weight);
	if (refSign<0) Hist::Head("hCOMn_TrInQ")->fill(refNRig, track->Qinner, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TrInQ%03d", tvBin))->fill(refIRig, track->Qinner, weight);
	if (tvBin!=0 && refSign>0) Hist::Head(StrFmt("hCOMp_TrInQ%03d", tvBin))->fill(refNRig, track->Qinner, weight);
	if (tvBin!=0 && refSign<0) Hist::Head(StrFmt("hCOMn_TrInQ%03d", tvBin))->fill(refNRig, track->Qinner, weight);

	if (track->Qinner < 0.8 || track->Qinner > 1.4) return false;
	Hist::Head("hCOM_Cutflow")->fill(3, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(3, weight);
	
	HitTRKInfo * trHitL1 = nullptr; Float_t trL1Q = 0.0; 
	HitTRKInfo * trHitL2 = nullptr; Float_t trL2Q = 0.0;
	HitTRKInfo * trHitL9 = nullptr; Float_t trL9Q = 0.0;
	for (Int_t ih = 0; track != nullptr && ih < track->hits.size(); ++ih) {
		HitTRKInfo &hit = track->hits.at(ih);
		//if (hit.clsId[0] >= 0 && (hit.chrg[0] < 0.6 || hit.chrg[0] > 3.4)) return false;
		//if (hit.clsId[1] >= 0 && (hit.chrg[1] < 0.6 || hit.chrg[1] > 3.4)) return false;
		if (hit.clsId[0] < 0 || hit.clsId[1] < 0) continue;
		Float_t chrg = 0.5 * (hit.chrg[0] + hit.chrg[1]);
		if (hit.layJ == 1) { trHitL1 = &hit; trL1Q = chrg; } 
		if (hit.layJ == 2) { trHitL2 = &hit; trL2Q = chrg; } 
		if (hit.layJ == 9) { trHitL9 = &hit; trL9Q = chrg; }
  }
	Hist::Head("hCOM_Cutflow")->fill(4, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(4, weight);
	
	// preselection (Mini-Trigger Requirement)
	if ((fTrg->bit&8) != 8) return false;
	Hist::Head("hCOM_Cutflow")->fill(5, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(5, weight);
  
	// preselection (Mini-TOF Requirement)
	if (!fTof->statusBetaH) return false;
	if (fTof->betaHPatt != 15) return false;
	if (fTof->betaH < 0.5 || fTof->betaH > 1.3) return false;
	Hist::Head("hCOM_Cutflow")->fill(6, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(6, weight);

  if (fTof->normChisqT > 10) return false;
  if (fTof->normChisqC > 10) return false;
	Hist::Head("hCOM_Cutflow")->fill(7, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(7, weight);
	
	Float_t tofQ   = 0.25 * (fTof->Q[0] + fTof->Q[1] + fTof->Q[2] + fTof->Q[3]);
	Float_t tofQu  = 0.50 * (fTof->Q[0] + fTof->Q[1]);
	Float_t tofQl  = 0.50 * (fTof->Q[2] + fTof->Q[3]);
	Float_t tofQdu = (fTof->Q[0] - fTof->Q[1]);
	Float_t tofQdl = (fTof->Q[2] - fTof->Q[3]);
	Float_t tofQul = tofQu - tofQl;
	
	Hist::Head("hCOMi_TofQ")->fill(refIRig, tofQ, weight);
	if (refSign>0) Hist::Head("hCOMp_TofQ")->fill(refNRig, tofQ, weight);
	if (refSign<0) Hist::Head("hCOMn_TofQ")->fill(refNRig, tofQ, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TofQ%03d", tvBin))->fill(refIRig, tofQ, weight);
	if (tvBin!=0 && refSign>0) Hist::Head(StrFmt("hCOMp_TofQ%03d", tvBin))->fill(refNRig, tofQ, weight);
	if (tvBin!=0 && refSign<0) Hist::Head(StrFmt("hCOMn_TofQ%03d", tvBin))->fill(refNRig, tofQ, weight);
	
	Hist::Head("hCOMi_TofQul")->fill(refIRig, tofQul, weight);
	if (refSign>0) Hist::Head("hCOMp_TofQul")->fill(refNRig, tofQul, weight);
	if (refSign<0) Hist::Head("hCOMn_TofQul")->fill(refNRig, tofQul, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOMi_TofQul%03d", tvBin))->fill(refIRig, tofQul, weight);
	if (tvBin!=0 && refSign>0) Hist::Head(StrFmt("hCOMp_TofQul%03d", tvBin))->fill(refNRig, tofQul, weight);
	if (tvBin!=0 && refSign<0) Hist::Head(StrFmt("hCOMn_TofQul%03d", tvBin))->fill(refNRig, tofQul, weight);
	
	if (tofQu < 0.8 || tofQu > 1.6) return false;
	if (tofQl < 0.8 || tofQl > 1.6) return false;
	Hist::Head("hCOM_Cutflow")->fill(8, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(8, weight);
	
	if (fTof->numOfInTimeCluster > 4) return false;
	Hist::Head("hCOM_Cutflow")->fill(9, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(9, weight);
  
	// preselection (Mini-TRD Requirement)
	if (fTrd->numOfVertexWithTrTrack > 10) return false;
	if (fTrd->numOfCluster < 8 || fTrd->numOfCluster > 100) return false;
	Hist::Head("hCOM_Cutflow")->fill(10, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(10, weight);

	// preselection (Mini-Shower Requirement)
	ShowerInfo * shower = (fEcal->showers.size() >= 1) ? (&fEcal->showers.at(0)) : nullptr;

	// preselection (Min-Track Pattern Requirement) 
	if (trHitL2) Hist::Head("hCOMi_TrL2Q")->fill(refIRig, trL2Q, weight);
	if (trHitL2 && refSign>0) Hist::Head("hCOMp_TrL2Q")->fill(refNRig, trL2Q, weight);
	if (trHitL2 && refSign<0) Hist::Head("hCOMn_TrL2Q")->fill(refNRig, trL2Q, weight);
	if (trHitL2 && tvBin!=0) Hist::Head(StrFmt("hCOMi_TrL2Q%03d", tvBin))->fill(refIRig, trL2Q, weight);
	if (trHitL2 && tvBin!=0 && refSign>0) Hist::Head(StrFmt("hCOMp_TrL2Q%03d", tvBin))->fill(refNRig, trL2Q, weight);
	if (trHitL2 && tvBin!=0 && refSign<0) Hist::Head(StrFmt("hCOMn_TrL2Q%03d", tvBin))->fill(refNRig, trL2Q, weight);

	if (trHitL1) Hist::Head("hCOMi_TrL1Q")->fill(refIRig, trL1Q, weight);
	if (trHitL1 && refSign>0) Hist::Head("hCOMp_TrL1Q")->fill(refNRig, trL1Q, weight);
	if (trHitL1 && refSign<0) Hist::Head("hCOMn_TrL1Q")->fill(refNRig, trL1Q, weight);
	if (trHitL1 && tvBin!=0) Hist::Head(StrFmt("hCOMi_TrL1Q%03d", tvBin))->fill(refIRig, trL1Q, weight);
	if (trHitL1 && tvBin!=0 && refSign>0) Hist::Head(StrFmt("hCOMp_TrL1Q%03d", tvBin))->fill(refNRig, trL1Q, weight);
	if (trHitL1 && tvBin!=0 && refSign<0) Hist::Head(StrFmt("hCOMn_TrL1Q%03d", tvBin))->fill(refNRig, trL1Q, weight);

	if (trHitL9) Hist::Head("hCOMi_TrL9Q")->fill(refIRig, trL9Q, weight);
	if (trHitL9 && refSign>0) Hist::Head("hCOMp_TrL9Q")->fill(refNRig, trL9Q, weight);
	if (trHitL9 && refSign<0) Hist::Head("hCOMn_TrL9Q")->fill(refNRig, trL9Q, weight);
	if (trHitL9 && tvBin!=0) Hist::Head(StrFmt("hCOMi_TrL9Q%03d", tvBin))->fill(refIRig, trL9Q, weight);
	if (trHitL9 && tvBin!=0 && refSign>0) Hist::Head(StrFmt("hCOMp_TrL9Q%03d", tvBin))->fill(refNRig, trL9Q, weight);
	if (trHitL9 && tvBin!=0 && refSign<0) Hist::Head(StrFmt("hCOMn_TrL9Q%03d", tvBin))->fill(refNRig, trL9Q, weight);
	
	Bool_t isTrIn = (track->status[0][2] && (track->bitPatt&  9)==  9 && track->status[0][0] && track->status[0][1]); // Inner Y + L2 XY
	Bool_t isTrL1 = (track->status[0][3] && (track->bitPatt& 33)== 33 && (trL1Q > 0.7 && trL1Q < 1.6)); // Inner Y + L1 XY
	Bool_t isTrL9 = (track->status[0][4] && (track->bitPatt&137)==137 && (trL9Q > 0.7 && trL9Q < 1.6)); // Inner Y + L2 XY + L9 XY
	Bool_t isTrFs = (track->status[0][5] && (track->bitPatt&161)==161 && (trL1Q > 0.7 && trL1Q < 1.6) && (trL9Q > 0.7 && trL9Q < 1.6)); // Inner Y + L1 XY + L9 XY
	Bool_t isTrPt = (isTrIn || isTrL1 || isTrL9 || isTrFs);
	
	if (!isTrPt) return false;
	Hist::Head("hCOM_Cutflow")->fill(11, weight);
	if (tvBin!=0) Hist::Head(StrFmt("hCOM_Cutflow%03d", tvBin))->fill(11, weight);

#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] L Selection [][]\n";
#endif
	/**** Low Energy Region ****/
	while (true) {
		// TRK
		Int_t trPt = 2;
		std::string trNm = "In";
		if (!isTrIn) break;
		Float_t trig       = track->rigidity[0][trPt];
		Int_t   sign       = MgntNum::Compare(trig);
		Float_t nrig       = std::fabs(trig);
		Float_t irig       = 1. / trig;
		Float_t tbtaPr     = 1. / std::sqrt(1.+(massPr/nrig)*(massPr/nrig));
		Float_t tbtaPi     = 1. / std::sqrt(1.+(massPi/nrig)*(massPi/nrig));
		Float_t tbtaEl     = 1. / std::sqrt(1.+(massEl/nrig)*(massEl/nrig));
		
		Hist::Head("hLi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);

		Int_t   patAgl = -1;
		Float_t cosAgl = std::acos(-track->stateL1[0][trPt][5]) * TMath::RadToDeg();
		if      (cosAgl < 25.0) patAgl = 0; 
		else if (cosAgl < 30.0) patAgl = 1; 
		else if (cosAgl < 35.0) patAgl = 2; 
		else if (cosAgl < 40.0) patAgl = 3; 
		if (patAgl == -1) break;

		// Cutoff
		if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
			Double_t cfRig = DST::CfStableFT * fRti->cutoffIGRF[patAgl];
			//Double_t cfRig = DST::CfStableFT * fRti->cutoffIGRF[3];
			if (nrig < cfRig) break;
		}
		
		Hist::Head("hLi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		// Charge Confusion
		Float_t lchixIn = std::log(track->chisq[0][2][0]);
		Float_t lchiyIn = std::log(track->chisq[0][2][1]);
		
		Hist::Head("hLi_LchixIn")->fill(irig, lchixIn, weight);
		if (sign>0) Hist::Head("hLp_LchixIn")->fill(nrig, lchixIn, weight);
		if (sign<0) Hist::Head("hLn_LchixIn")->fill(nrig, lchixIn, weight);
		Hist::Head("hLi_LchiyIn")->fill(irig, lchiyIn, weight);
		if (sign>0) Hist::Head("hLp_LchiyIn")->fill(nrig, lchiyIn, weight);
		if (sign<0) Hist::Head("hLn_LchiyIn")->fill(nrig, lchiyIn, weight);
		
		if (tvBin!=0) Hist::Head(StrFmt("hLi_LchixIn%03d", tvBin))->fill(irig, lchixIn, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_LchixIn%03d", tvBin))->fill(nrig, lchixIn, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_LchixIn%03d", tvBin))->fill(nrig, lchixIn, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_LchiyIn%03d", tvBin))->fill(irig, lchiyIn, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_LchiyIn%03d", tvBin))->fill(nrig, lchiyIn, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_LchiyIn%03d", tvBin))->fill(nrig, lchiyIn, weight);

		const Float_t lasymSGM = 1.8;
		const Float_t lasymLMT = 1.0e-8;
		Float_t trSgmIU = GetRigSgm(0, nrig, massPr);
		Float_t trSgmIL = GetRigSgm(1, nrig, massPr);
		Float_t trSgmIn = GetRigSgm(2, nrig, massPr);

		Float_t crrParUL[6] = { 7.63668e-01, 2.13220e-01, -1.40775e+00, 9.60481e-01, 3.50012e+00, 1.49281e+00 };
		Float_t crrSgmUL    = crrParUL[0]*(std::erf(crrParUL[1]*(std::log(nrig)+crrParUL[2]))+1.0) + TMath::Landau(crrParUL[3], crrParUL[4], crrParUL[5]);
		Float_t trSgmUL     = std::sqrt(trSgmIU * trSgmIU + trSgmIL * trSgmIL) * crrSgmUL;
		Float_t asymUL      = (1.0/track->rigidity[0][0] - 1.0/track->rigidity[0][1]) / trSgmUL;
		Float_t lasymUL     = std::log(asymUL * asymUL + lasymLMT) / lasymSGM;
		
		Hist::Head("hLi_LasymUL")->fill(irig, lasymUL, weight);
		if (sign>0) Hist::Head("hLp_LasymUL")->fill(nrig, lasymUL, weight);
		if (sign<0) Hist::Head("hLn_LasymUL")->fill(nrig, lasymUL, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_LasymUL%03d", tvBin))->fill(irig, lasymUL, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_LasymUL%03d", tvBin))->fill(nrig, lasymUL, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_LasymUL%03d", tvBin))->fill(nrig, lasymUL, weight);
	
		if (lchixIn > 3.5) break;
		if (lchiyIn > 2.3) break;
		if (lasymUL > 1.8) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);

		// Charge Confusion
		Float_t ccest = 0.;
		if (trPt == 2) {
			Float_t ccelm[2] = { lchiyIn, lasymUL };
			Float_t ccwgt[2] = { GetSqrMDR(2), 0.5*(GetSqrMDR(0)+GetSqrMDR(1)) };
			Float_t ccsum = 0.0, ccsgm = 0.0;
			for (Int_t it = 0; it < 2; ++it) { ccsum += ccelm[it] * ccwgt[it]; ccsgm += ccwgt[it]; }
			ccest = ccsum / ccsgm;
		}
		else ccest = lchiyIn;
		
		Hist::Head("hLi_CCestIn")->fill(irig, ccest, weight);
		if (sign>0) Hist::Head("hLp_CCestIn")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hLn_CCestIn")->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_CCestIn%03d", tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_CCestIn%03d", tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_CCestIn%03d", tvBin))->fill(nrig, ccest, weight);
	
		Float_t CCcutIn = 1.75 * std::erf((std::log(40.0) - std::log(nrig)) / 1.7); // 50
		if (trPt == 2 && ccest > CCcutIn) break;

		Hist::Head("hLi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);

		// TRD		
		Int_t trdPt = 0;
		Bool_t  hastrd  = (fTrd->statusKCls[trdPt] && fTrd->LLR_nhit[trdPt] >= 8);
		if (!hastrd) break;
		
		Float_t trdl    = (fTrd->LLR[trdPt][0]);
		Bool_t  istrdHe = (fTrd->LLR[trdPt][2] > 0.3);
		Bool_t  istrdPr = (trdl > (0.8 - (0.1/nrig)));
		if (istrdHe) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 4, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 4, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 4, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 4, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		
		Hist::Head("hLi_Trdl")->fill(irig, trdl, weight);
		if (sign>0) Hist::Head("hLp_Trdl")->fill(nrig, trdl, weight);
		if (sign<0) Hist::Head("hLn_Trdl")->fill(nrig, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		
		if (!istrdPr) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 5, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 5, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 5, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 5, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		
		// RICH
		Int_t   richph    = fRich->numOfHit;
		Float_t richelph  = (fRich->kindOfRad == -1) ? 0.0 : fRich->numOfExpPE[0];
		Float_t richprph  = (fRich->kindOfRad == -1) ? 0.0 : fRich->numOfExpPE[3];
		
		Bool_t  hasrich   = (fRich->status && fRich->kindOfRad != -1);
		Float_t richbPr   = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbtaPr);
		Float_t richbPi   = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbtaPi);
		Float_t richbEl   = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbtaEl);
		Float_t tofrichb  = (!hasrich) ? 0.0 : (fTof->betaH - fRich->beta);

		if (fRich->kindOfRad == -1) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 6, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 6, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 6, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 6, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 6, weight);

		std::string richnm = (fRich->kindOfRad == 0) ? "Agl" : "Naf";
		Float_t richbth   = (fRich->kindOfRad == 0) ? 0.003 : 0.007;
		Float_t richbsh   = std::sqrt(1.0+(0.4/nrig)*(0.4/nrig))-1.0;
		Bool_t  isrichNo  = (!fRich->status && fRich->kindOfRad == 0 && MgntNum::Compare(richelph, 2.0f) > 0 && MgntNum::Compare(richprph, 0.5f) < 0 && richph <= 2);
		Bool_t  isrichPr  = ( fRich->status && std::fabs(richbPr) < (richbth + richbsh));
		Bool_t  isrichPi  = ( fRich->status && std::fabs(richbPi) < richbth);
		Bool_t  isrichEl  = ( fRich->status && std::fabs(richbEl) < richbth);
		Bool_t  cutrich   = (fRich->status) ? (std::fabs(richbPr) > (richbth + richbsh) || ((fRich->kindOfRad == 0) ? richph >= 25 : richph >= 25)) : 
		                                      ((fRich->kindOfRad == 0) ? !(MgntNum::Compare(richelph, 2.0f) > 0 && richph <= 10) : !(MgntNum::Compare(richelph, 0.8f) > 0 && richph <= 6));
		
		std::string patnm = (fRich->status) ? "HV" : "NO";

		Hist::Head(StrFmt("hLi_%sph%s", richnm.c_str(), patnm.c_str()))->fill(irig, richph, weight);
		if (sign>0) Hist::Head(StrFmt("hLp_%sph%s", richnm.c_str(), patnm.c_str()))->fill(nrig, richph, weight);
		if (sign<0) Hist::Head(StrFmt("hLn_%sph%s", richnm.c_str(), patnm.c_str()))->fill(nrig, richph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_%sph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(irig, richph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_%sph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(nrig, richph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_%sph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(nrig, richph, weight);
		
		Hist::Head(StrFmt("hLi_%sprph%s", richnm.c_str(), patnm.c_str()))->fill(irig, richprph, weight);
		if (sign>0) Hist::Head(StrFmt("hLp_%sprph%s", richnm.c_str(), patnm.c_str()))->fill(nrig, richprph, weight);
		if (sign<0) Hist::Head(StrFmt("hLn_%sprph%s", richnm.c_str(), patnm.c_str()))->fill(nrig, richprph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_%sprph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(irig, richprph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_%sprph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(nrig, richprph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_%sprph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(nrig, richprph, weight);
		
		Hist::Head(StrFmt("hLi_%selph%s", richnm.c_str(), patnm.c_str()))->fill(irig, richelph, weight);
		if (sign>0) Hist::Head(StrFmt("hLp_%selph%s", richnm.c_str(), patnm.c_str()))->fill(nrig, richelph, weight);
		if (sign<0) Hist::Head(StrFmt("hLn_%selph%s", richnm.c_str(), patnm.c_str()))->fill(nrig, richelph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_%selph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(irig, richelph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_%selph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(nrig, richelph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_%selph%s%03d", richnm.c_str(), patnm.c_str(), tvBin))->fill(nrig, richelph, weight);
		
		if (hasrich) Hist::Head(StrFmt("hLi_%sb", richnm.c_str()))->fill(irig, richbPr, weight);
		if (hasrich && sign>0) Hist::Head(StrFmt("hLp_%sb", richnm.c_str()))->fill(nrig, richbPr, weight);
		if (hasrich && sign<0) Hist::Head(StrFmt("hLn_%sb", richnm.c_str()))->fill(nrig, richbPr, weight);
		if (tvBin!=0 && hasrich) Hist::Head(StrFmt("hLi_%sb%03d", richnm.c_str(), tvBin))->fill(irig, richbPr, weight);
		if (tvBin!=0 && hasrich && sign>0) Hist::Head(StrFmt("hLp_%sb%03d", richnm.c_str(), tvBin))->fill(nrig, richbPr, weight);
		if (tvBin!=0 && hasrich && sign<0) Hist::Head(StrFmt("hLn_%sb%03d", richnm.c_str(), tvBin))->fill(nrig, richbPr, weight);
		
		if (hasrich) Hist::Head(StrFmt("hLi_Tof%sb", richnm.c_str()))->fill(irig, tofrichb, weight);
		if (hasrich && sign>0) Hist::Head(StrFmt("hLp_Tof%sb", richnm.c_str()))->fill(nrig, tofrichb, weight);
		if (hasrich && sign<0) Hist::Head(StrFmt("hLn_Tof%sb", richnm.c_str()))->fill(nrig, tofrichb, weight);
		if (tvBin!=0 && hasrich) Hist::Head(StrFmt("hLi_Tof%sb%03d", richnm.c_str(), tvBin))->fill(irig, tofrichb, weight);
		if (tvBin!=0 && hasrich && sign>0) Hist::Head(StrFmt("hLp_Tof%sb%03d", richnm.c_str(), tvBin))->fill(nrig, tofrichb, weight);
		if (tvBin!=0 && hasrich && sign<0) Hist::Head(StrFmt("hLn_Tof%sb%03d", richnm.c_str(), tvBin))->fill(nrig, tofrichb, weight);
		
		// TOF
		Hist::Head("hLi_TofulChrg")->fill(irig, tofQul, weight);
		if (sign>0) Hist::Head("hLp_TofulChrg")->fill(nrig, tofQul, weight);
		if (sign<0) Hist::Head("hLn_TofulChrg")->fill(nrig, tofQul, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_TofulChrg%03d", tvBin))->fill(irig, tofQul, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_TofulChrg%03d", tvBin))->fill(nrig, tofQul, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_TofulChrg%03d", tvBin))->fill(nrig, tofQul, weight);

		Float_t cutTofQul = 0.15 * (TMath::Erf(2.0 * (nrig - 2.5)) + 1.0) + 0.15;
		if (std::fabs(tofQul) > cutTofQul) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 7, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 7, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 7, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 7, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 7, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 7, weight);

		Float_t tofb = (fTof->betaH - tbtaPr);

		Bool_t isSignal     = (sign > 0 && (isrichNo || isrichPr));
		Bool_t isBackground = (sign < 0 && (isrichPi || isrichEl));
		if (isSignal)     Hist::Head("hLs_Tofb")->fill(nrig, tofb, weight);
		if (isBackground) Hist::Head("hLb_Tofb")->fill(nrig, tofb, weight);
		if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hLs_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hLb_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		
		if (cutrich) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 8, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 8, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 8, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 8, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 8, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 8, weight);
		
		Hist::Head("hLi_Tofb")->fill(irig, tofb, weight);
		if (sign>0) Hist::Head("hLp_Tofb")->fill(nrig, tofb, weight);
		if (sign<0) Hist::Head("hLn_Tofb")->fill(nrig, tofb, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Tofb%03d", tvBin))->fill(irig, tofb, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Tofb%03d", tvBin))->fill(nrig, tofb, weight);

		Hist::Head("hLi_Evt")->fill(irig, weight);
		if (sign>0) Hist::Head("hLp_Evt")->fill(nrig, weight);
		if (sign<0) Hist::Head("hLn_Evt")->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Evt%03d", tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Evt%03d", tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Evt%03d", tvBin))->fill(nrig, weight);
		
		if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
			Hist::Head("hLi_MCTrRso")->fill(MCIRig, irig, weight);
			if (MCTRig>0.0) Hist::Head("hLp_MCTrRso")->fill(MCNRig, (irig - MCIRig), weight);
			if (MCTRig<0.0) Hist::Head("hLn_MCTrRso")->fill(MCNRig, (irig - MCIRig), weight);
		
			Hist::Head("hLi_MCEvt")->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head("hLp_MCEvt")->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head("hLn_MCEvt")->fill(MCNRig, weight);
		}
	
		break;
	}


#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] I Selection [][]\n";
#endif
	/**** Intermedia Energy Region ****/
	while (true) {
		// TRK
		Int_t trPt = -1;
		std::string trNm = "";
		if      (track->status[0][5] && (track->bitPatt&161)==161) { trPt = 5; trNm = "Fs"; }
		//else if (track->status[0][4] && (track->bitPatt&129)==129) { trPt = 4; trNm = "L9"; }
		else if (track->status[0][4] && (track->bitPatt&137)==137) { trPt = 4; trNm = "L9"; }
		else if (track->status[0][3] && (track->bitPatt& 33)== 33) { trPt = 3; trNm = "L1"; }
		else if (track->status[0][2] && (track->bitPatt& 10)== 10 && track->status[0][0] && track->status[0][1]) { trPt = 2; trNm = "In"; }
		//if (track->status[0][2] && (track->bitPatt&  10)==10 && track->status[0][0] && track->status[0][1]) { trPt = 2; trNm = "In"; }
		if (trPt == -1) break;
		Float_t trig       = track->rigidity[0][trPt];
		Int_t   sign       = MgntNum::Compare(trig);
		Float_t nrig       = std::fabs(trig);
		Float_t irig       = 1. / trig;
		Float_t tbtaPr     = 1. / std::sqrt(1.+(massPr/nrig)*(massPr/nrig));
		Float_t tbtaPi     = 1. / std::sqrt(1.+(massPi/nrig)*(massPi/nrig));
		Float_t tbtaEl     = 1. / std::sqrt(1.+(massEl/nrig)*(massEl/nrig));
		
		Hist::Head("hIi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		
		Int_t   patAgl = -1;
		Float_t cosAgl = std::acos(-track->stateL1[0][trPt][5]) * TMath::RadToDeg();
		if      (cosAgl < 25.0) patAgl = 0; 
		else if (cosAgl < 30.0) patAgl = 1; 
		else if (cosAgl < 35.0) patAgl = 2; 
		else if (cosAgl < 40.0) patAgl = 3; 
		if (patAgl == -1) break;
		
		// Cutoff
		if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
			Double_t cfRig = DST::CfStableFT * fRti->cutoffIGRF[patAgl];
			//Double_t cfRig = DST::CfStableFT * fRti->cutoffIGRF[3];
			if (nrig < cfRig) break;
		}
		
		Hist::Head("hIi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		// TOF
		Float_t cutTofQul = 0.15 * (TMath::Erf(2.0 * (nrig - 2.5)) + 1.0) + 0.15;
		if (std::fabs(tofQul) > cutTofQul) break;
			
		// TRD
		Int_t trdPt = 1;
		Bool_t  hastrd  = (fTrd->statusKCls[trdPt] && fTrd->LLR_nhit[trdPt] >= 8);
		if (!hastrd) break;
		
		Float_t trdl    = (fTrd->LLR[trdPt][0]);
		Bool_t  istrdHe = (fTrd->LLR[trdPt][2] > 0.3);
		if (istrdHe) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		
		// ECAL
		Bool_t  hasecal  = (shower != nullptr);
		Float_t ecalbdt  = (!hasecal) ? -2.0 : shower->PisaBDT;
		Bool_t  isecalPr = (!hasecal) ? false : (ecalbdt < -0.8);
		Bool_t  isecalEl = (!hasecal) ? false : (ecalbdt >  0.6);
		
		if (hasecal) Hist::Head("hIi_Ecalbdt")->fill(irig, ecalbdt, weight);
		if (hasecal && sign>0) Hist::Head("hIp_Ecalbdt")->fill(nrig, ecalbdt, weight);
		if (hasecal && sign<0) Hist::Head("hIn_Ecalbdt")->fill(nrig, ecalbdt, weight);
		if (tvBin!=0 && hasecal) Hist::Head(StrFmt("hIi_Ecalbdt%03d", tvBin))->fill(irig, ecalbdt, weight);
		if (tvBin!=0 && hasecal && sign>0) Hist::Head(StrFmt("hIp_Ecalbdt%03d", tvBin))->fill(nrig, ecalbdt, weight);
		if (tvBin!=0 && hasecal && sign<0) Hist::Head(StrFmt("hIn_Ecalbdt%03d", tvBin))->fill(nrig, ecalbdt, weight);
		
		// Charge Confusion
		const Float_t lasymSGM = 1.8;
		const Float_t lasymLMT = 1.0e-8;
		Float_t trSgmIU = GetRigSgm(0, nrig, massPr);
		Float_t trSgmIL = GetRigSgm(1, nrig, massPr);
		Float_t trSgmIn = GetRigSgm(2, nrig, massPr);
		Float_t trSgmL1 = GetRigSgm(3, nrig, massPr);
		Float_t trSgmL9 = GetRigSgm(4, nrig, massPr);
		Float_t trSgmFs = GetRigSgm(5, nrig, massPr);
		
		Float_t trSgmPt = 1.0 / nrig;
		if      (trPt == 2) trSgmPt = trSgmIn;
		else if (trPt == 3) trSgmPt = trSgmL1;
		else if (trPt == 4) trSgmPt = trSgmL9;
		else if (trPt == 5) trSgmPt = trSgmFs;

		Bool_t  hasasymUL   = (track->status[0][0] && track->status[0][1]);
		Float_t crrParUL[6] = { 7.63668e-01, 2.13220e-01, -1.40775e+00, 9.60481e-01, 3.50012e+00, 1.49281e+00 };
		Float_t crrSgmUL    = crrParUL[0]*(std::erf(crrParUL[1]*(std::log(nrig)+crrParUL[2]))+1.0) + TMath::Landau(crrParUL[3], crrParUL[4], crrParUL[5]);
		Float_t trSgmUL     = (!hasasymUL) ? 0.0 : std::sqrt(trSgmIU * trSgmIU + trSgmIL * trSgmIL) * crrSgmUL;
		Float_t asymUL      = (!hasasymUL) ? 0.0 : (1.0/track->rigidity[0][0] - 1.0/track->rigidity[0][1]) / trSgmUL;
		Float_t lasymUL     = (!hasasymUL) ? 0.0 : std::log(asymUL * asymUL + lasymLMT) / lasymSGM;

		Bool_t  hasasym1I   = (track->status[0][2] && track->status[0][3]);
		Float_t crrPar1I[5] = { 4.93251e-01, 1.90417e+00, 9.18502e+01, 4.86389e+00, 3.29747e-01 };
		Float_t crrSgm1I    = crrPar1I[0]*std::erf(crrPar1I[1]*(std::log(nrig+crrPar1I[2])-crrPar1I[3]))+crrPar1I[4];
		Float_t trSgm1I     = (!hasasym1I) ? 0.0 : std::sqrt(trSgmIn * trSgmIn + trSgmL1 * trSgmL1) * crrSgm1I;
		Float_t asym1I      = (!hasasym1I) ? 0.0 : (1.0/track->rigidity[0][2] - 1.0/track->rigidity[0][3]) / trSgm1I; 
		Float_t lasym1I     = (!hasasym1I) ? 0.0 : std::log(asym1I * asym1I + lasymLMT) / lasymSGM;

		Bool_t  hasasym9I   = (track->status[0][2] && track->status[0][4]);
		Float_t crrPar9I[5] = { 4.17294e-01, 1.15603e+00, 1.62194e+01, 3.81780e+00, 3.92404e-01 };
		Float_t crrSgm9I    = crrPar9I[0]*std::erf(crrPar9I[1]*(std::log(nrig+crrPar9I[2])-crrPar9I[3]))+crrPar9I[4];
		Float_t trSgm9I     = (!hasasym9I) ? 0.0 : std::sqrt(trSgmIn * trSgmIn + trSgmL9 * trSgmL9) * crrSgm9I;
		Float_t asym9I      = (!hasasym9I) ? 0.0 : (1.0/track->rigidity[0][2] - 1.0/track->rigidity[0][4]) / trSgm9I;
		Float_t lasym9I     = (!hasasym9I) ? 0.0 : std::log(asym9I * asym9I + lasymLMT) / lasymSGM;

		Bool_t  hasasym91   = (track->status[0][3] && track->status[0][4]);
		Float_t crrPar91[5] = { 5.25352e-01, 9.88059e-01, 1.21200e+01, 3.72027e+00, 4.99424e-01 };
		Float_t crrSgm91    = crrPar91[0]*std::erf(crrPar91[1]*(std::log(nrig+crrPar91[2])-crrPar91[3]))+crrPar91[4];
		Float_t trSgm91     = (!hasasym91) ? 0.0 : std::sqrt(trSgmL1 * trSgmL1 + trSgmL9 * trSgmL9) * crrSgm91;
		Float_t asym91      = (!hasasym91) ? 0.0 : (1.0/track->rigidity[0][3] - 1.0/track->rigidity[0][4]) / trSgm91;
		Float_t lasym91     = (!hasasym91) ? 0.0 : std::log(asym91 * asym91 + lasymLMT) / lasymSGM;
	
		if (hasasymUL) Hist::Head("hIi_LasymUL")->fill(irig, lasymUL, weight);
		if (hasasymUL && sign>0) Hist::Head("hIp_LasymUL")->fill(nrig, lasymUL, weight);
		if (hasasymUL && sign<0) Hist::Head("hIn_LasymUL")->fill(nrig, lasymUL, weight);
		if (hasasym1I) Hist::Head("hIi_Lasym1I")->fill(irig, lasym1I, weight);
		if (hasasym1I && sign>0) Hist::Head("hIp_Lasym1I")->fill(nrig, lasym1I, weight);
		if (hasasym1I && sign<0) Hist::Head("hIn_Lasym1I")->fill(nrig, lasym1I, weight);
		if (hasasym9I) Hist::Head("hIi_Lasym9I")->fill(irig, lasym9I, weight);
		if (hasasym9I && sign>0) Hist::Head("hIp_Lasym9I")->fill(nrig, lasym9I, weight);
		if (hasasym9I && sign<0) Hist::Head("hIn_Lasym9I")->fill(nrig, lasym9I, weight);
		if (hasasym91) Hist::Head("hIi_Lasym91")->fill(irig, lasym91, weight);
		if (hasasym91 && sign>0) Hist::Head("hIp_Lasym91")->fill(nrig, lasym91, weight);
		if (hasasym91 && sign<0) Hist::Head("hIn_Lasym91")->fill(nrig, lasym91, weight);
		
		if (tvBin!=0 && hasasymUL) Hist::Head(StrFmt("hIi_LasymUL%03d", tvBin))->fill(irig, lasymUL, weight);
		if (tvBin!=0 && hasasymUL && sign>0) Hist::Head(StrFmt("hIp_LasymUL%03d", tvBin))->fill(nrig, lasymUL, weight);
		if (tvBin!=0 && hasasymUL && sign<0) Hist::Head(StrFmt("hIn_LasymUL%03d", tvBin))->fill(nrig, lasymUL, weight);
		if (tvBin!=0 && hasasym1I) Hist::Head(StrFmt("hIi_Lasym1I%03d", tvBin))->fill(irig, lasym1I, weight);
		if (tvBin!=0 && hasasym1I && sign>0) Hist::Head(StrFmt("hIp_Lasym1I%03d", tvBin))->fill(nrig, lasym1I, weight);
		if (tvBin!=0 && hasasym1I && sign<0) Hist::Head(StrFmt("hIn_Lasym1I%03d", tvBin))->fill(nrig, lasym1I, weight);
		if (tvBin!=0 && hasasym9I) Hist::Head(StrFmt("hIi_Lasym9I%03d", tvBin))->fill(irig, lasym9I, weight);
		if (tvBin!=0 && hasasym9I && sign>0) Hist::Head(StrFmt("hIp_Lasym9I%03d", tvBin))->fill(nrig, lasym9I, weight);
		if (tvBin!=0 && hasasym9I && sign<0) Hist::Head(StrFmt("hIn_Lasym9I%03d", tvBin))->fill(nrig, lasym9I, weight);
		if (tvBin!=0 && hasasym91) Hist::Head(StrFmt("hIi_Lasym91%03d", tvBin))->fill(irig, lasym91, weight);
		if (tvBin!=0 && hasasym91 && sign>0) Hist::Head(StrFmt("hIp_Lasym91%03d", tvBin))->fill(nrig, lasym91, weight);
		if (tvBin!=0 && hasasym91 && sign<0) Hist::Head(StrFmt("hIn_Lasym91%03d", tvBin))->fill(nrig, lasym91, weight);
			
		// Charge Confusion
		Float_t lchixIn = (!track->status[0][2]) ? 0.0 : std::log(track->chisq[0][2][0]);
		Float_t lchiyIn = (!track->status[0][2]) ? 0.0 : std::log(track->chisq[0][2][1]);
		Float_t lchixL1 = (!track->status[0][3]) ? 0.0 : std::log(track->chisq[0][3][0]);
		Float_t lchiyL1 = (!track->status[0][3]) ? 0.0 : std::log(track->chisq[0][3][1]);
		Float_t lchixL9 = (!track->status[0][4]) ? 0.0 : std::log(track->chisq[0][4][0]);
		Float_t lchiyL9 = (!track->status[0][4]) ? 0.0 : std::log(track->chisq[0][4][1]);
		Float_t lchixFs = (!track->status[0][5]) ? 0.0 : std::log(track->chisq[0][5][0]);
		Float_t lchiyFs = (!track->status[0][5]) ? 0.0 : std::log(track->chisq[0][5][1]);
		
		if (track->status[0][2]) Hist::Head("hIi_LchixIn")->fill(irig, lchixIn, weight);
		if (track->status[0][2] && sign>0) Hist::Head("hIp_LchixIn")->fill(nrig, lchixIn, weight);
		if (track->status[0][2] && sign<0) Hist::Head("hIn_LchixIn")->fill(nrig, lchixIn, weight);
		if (track->status[0][2]) Hist::Head("hIi_LchiyIn")->fill(irig, lchiyIn, weight);
		if (track->status[0][2] && sign>0) Hist::Head("hIp_LchiyIn")->fill(nrig, lchiyIn, weight);
		if (track->status[0][2] && sign<0) Hist::Head("hIn_LchiyIn")->fill(nrig, lchiyIn, weight);
		if (track->status[0][3]) Hist::Head("hIi_LchixL1")->fill(irig, lchixL1, weight);
		if (track->status[0][3] && sign>0) Hist::Head("hIp_LchixL1")->fill(nrig, lchixL1, weight);
		if (track->status[0][3] && sign<0) Hist::Head("hIn_LchixL1")->fill(nrig, lchixL1, weight);
		if (track->status[0][3]) Hist::Head("hIi_LchiyL1")->fill(irig, lchiyL1, weight);
		if (track->status[0][3] && sign>0) Hist::Head("hIp_LchiyL1")->fill(nrig, lchiyL1, weight);
		if (track->status[0][3] && sign<0) Hist::Head("hIn_LchiyL1")->fill(nrig, lchiyL1, weight);
		if (track->status[0][4]) Hist::Head("hIi_LchixL9")->fill(irig, lchixL9, weight);
		if (track->status[0][4] && sign>0) Hist::Head("hIp_LchixL9")->fill(nrig, lchixL9, weight);
		if (track->status[0][4] && sign<0) Hist::Head("hIn_LchixL9")->fill(nrig, lchixL9, weight);
		if (track->status[0][4]) Hist::Head("hIi_LchiyL9")->fill(irig, lchiyL9, weight);
		if (track->status[0][4] && sign>0) Hist::Head("hIp_LchiyL9")->fill(nrig, lchiyL9, weight);
		if (track->status[0][4] && sign<0) Hist::Head("hIn_LchiyL9")->fill(nrig, lchiyL9, weight);
		if (track->status[0][5]) Hist::Head("hIi_LchixFs")->fill(irig, lchixFs, weight);
		if (track->status[0][5] && sign>0) Hist::Head("hIp_LchixFs")->fill(nrig, lchixFs, weight);
		if (track->status[0][5] && sign<0) Hist::Head("hIn_LchixFs")->fill(nrig, lchixFs, weight);
		if (track->status[0][5]) Hist::Head("hIi_LchiyFs")->fill(irig, lchiyFs, weight);
		if (track->status[0][5] && sign>0) Hist::Head("hIp_LchiyFs")->fill(nrig, lchiyFs, weight);
		if (track->status[0][5] && sign<0) Hist::Head("hIn_LchiyFs")->fill(nrig, lchiyFs, weight);
		
		if (tvBin!=0 && track->status[0][2]) Hist::Head(StrFmt("hIi_LchixIn%03d", tvBin))->fill(irig, lchixIn, weight);
		if (tvBin!=0 && track->status[0][2] && sign>0) Hist::Head(StrFmt("hIp_LchixIn%03d", tvBin))->fill(nrig, lchixIn, weight);
		if (tvBin!=0 && track->status[0][2] && sign<0) Hist::Head(StrFmt("hIn_LchixIn%03d", tvBin))->fill(nrig, lchixIn, weight);
		if (tvBin!=0 && track->status[0][2]) Hist::Head(StrFmt("hIi_LchiyIn%03d", tvBin))->fill(irig, lchiyIn, weight);
		if (tvBin!=0 && track->status[0][2] && sign>0) Hist::Head(StrFmt("hIp_LchiyIn%03d", tvBin))->fill(nrig, lchiyIn, weight);
		if (tvBin!=0 && track->status[0][2] && sign<0) Hist::Head(StrFmt("hIn_LchiyIn%03d", tvBin))->fill(nrig, lchiyIn, weight);
		if (tvBin!=0 && track->status[0][3]) Hist::Head(StrFmt("hIi_LchixL1%03d", tvBin))->fill(irig, lchixL1, weight);
		if (tvBin!=0 && track->status[0][3] && sign>0) Hist::Head(StrFmt("hIp_LchixL1%03d", tvBin))->fill(nrig, lchixL1, weight);
		if (tvBin!=0 && track->status[0][3] && sign<0) Hist::Head(StrFmt("hIn_LchixL1%03d", tvBin))->fill(nrig, lchixL1, weight);
		if (tvBin!=0 && track->status[0][3]) Hist::Head(StrFmt("hIi_LchiyL1%03d", tvBin))->fill(irig, lchiyL1, weight);
		if (tvBin!=0 && track->status[0][3] && sign>0) Hist::Head(StrFmt("hIp_LchiyL1%03d", tvBin))->fill(nrig, lchiyL1, weight);
		if (tvBin!=0 && track->status[0][3] && sign<0) Hist::Head(StrFmt("hIn_LchiyL1%03d", tvBin))->fill(nrig, lchiyL1, weight);
		if (tvBin!=0 && track->status[0][4]) Hist::Head(StrFmt("hIi_LchixL9%03d", tvBin))->fill(irig, lchixL9, weight);
		if (tvBin!=0 && track->status[0][4] && sign>0) Hist::Head(StrFmt("hIp_LchixL9%03d", tvBin))->fill(nrig, lchixL9, weight);
		if (tvBin!=0 && track->status[0][4] && sign<0) Hist::Head(StrFmt("hIn_LchixL9%03d", tvBin))->fill(nrig, lchixL9, weight);
		if (tvBin!=0 && track->status[0][4]) Hist::Head(StrFmt("hIi_LchiyL9%03d", tvBin))->fill(irig, lchiyL9, weight);
		if (tvBin!=0 && track->status[0][4] && sign>0) Hist::Head(StrFmt("hIp_LchiyL9%03d", tvBin))->fill(nrig, lchiyL9, weight);
		if (tvBin!=0 && track->status[0][4] && sign<0) Hist::Head(StrFmt("hIn_LchiyL9%03d", tvBin))->fill(nrig, lchiyL9, weight);
		if (tvBin!=0 && track->status[0][5]) Hist::Head(StrFmt("hIi_LchixFs%03d", tvBin))->fill(irig, lchixFs, weight);
		if (tvBin!=0 && track->status[0][5] && sign>0) Hist::Head(StrFmt("hIp_LchixFs%03d", tvBin))->fill(nrig, lchixFs, weight);
		if (tvBin!=0 && track->status[0][5] && sign<0) Hist::Head(StrFmt("hIn_LchixFs%03d", tvBin))->fill(nrig, lchixFs, weight);
		if (tvBin!=0 && track->status[0][5]) Hist::Head(StrFmt("hIi_LchiyFs%03d", tvBin))->fill(irig, lchiyFs, weight);
		if (tvBin!=0 && track->status[0][5] && sign>0) Hist::Head(StrFmt("hIp_LchiyFs%03d", tvBin))->fill(nrig, lchiyFs, weight);
		if (tvBin!=0 && track->status[0][5] && sign<0) Hist::Head(StrFmt("hIn_LchiyFs%03d", tvBin))->fill(nrig, lchiyFs, weight);
			
		// Charge Confusion
		if (trPt == 5 && (lchixL1 > 3.5 || lchixL9 > 3.5 || lchixFs > 3.5)) break;
		if (trPt == 5 && (lchiyL1 > 3.5 || lchiyL9 > 2.5 || lchiyFs > 2.5)) break; // lchiyFs 2.5   new (2.2  1.9)
		if (trPt == 5 && (lasym1I > 2.0 || lasym9I > 2.0 || lasym91 > 2.0)) break; // lasym91 2.0   new (1.8  1.5)
		if (trPt == 4 && (lchixL9 > 3.5)) break;
		if (trPt == 4 && (lchiyL9 > 1.8)) break; // 2.5
		if (trPt == 4 && (lasym9I > 1.5)) break; // 2.0
		if (trPt == 3 && (lchixL1 > 3.5)) break;
		if (trPt == 3 && (lchiyL1 > 2.4)) break; // 3.5
		if (trPt == 3 && (lasym1I > 1.6)) break; // 2.0
		if (trPt == 2 && (lchixIn > 3.5)) break;
		if (trPt == 2 && (lchiyIn > 2.3)) break; // 3.4
		if (trPt == 2 && (lasymUL > 1.8)) break; // 1.8
		
		Hist::Head("hIi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		
		// Charge Confusion
		Float_t ccest = 0.;
		if (trPt == 5) {
			Float_t ccelm[6] = { lchiyL1, lchiyL9, lchiyFs, lasym1I, lasym9I, lasym91 };
			Float_t ccwgt[6] = { GetSqrMDR(3), GetSqrMDR(4), GetSqrMDR(5), 0.5*(GetSqrMDR(2)+GetSqrMDR(3)), 0.5*(GetSqrMDR(2)+GetSqrMDR(4)), 0.5*(GetSqrMDR(4)+GetSqrMDR(5)) };
			Float_t ccsum = 0.0, ccsgm = 0.0;
			for (Int_t it = 0; it < 6; ++it) { ccsum += ccelm[it] * ccwgt[it]; ccsgm += ccwgt[it]; }
			ccest = ccsum / ccsgm;
		}
		if (trPt == 4) {
			Float_t ccelm[2] = { lchiyL9, lasym9I };
			Float_t ccwgt[2] = { GetSqrMDR(4), 0.5*(GetSqrMDR(2)+GetSqrMDR(4)) };
			Float_t ccsum = 0.0, ccsgm = 0.0;
			for (Int_t it = 0; it < 2; ++it) { ccsum += ccelm[it] * ccwgt[it]; ccsgm += ccwgt[it]; }
			ccest = ccsum / ccsgm;
		}
		if (trPt == 3) {
			Float_t ccelm[2] = { lchiyL1, lasym1I };
			Float_t ccwgt[2] = { GetSqrMDR(3), 0.5*(GetSqrMDR(2)+GetSqrMDR(3)) };
			Float_t ccsum = 0.0, ccsgm = 0.0;
			for (Int_t it = 0; it < 2; ++it) { ccsum += ccelm[it] * ccwgt[it]; ccsgm += ccwgt[it]; }
			ccest = ccsum / ccsgm;
		}
		if (trPt == 2) {
			Float_t ccelm[2] = { lchiyIn, lasymUL };
			Float_t ccwgt[2] = { GetSqrMDR(2), 0.5*(GetSqrMDR(0)+GetSqrMDR(1)) };
			Float_t ccsum = 0.0, ccsgm = 0.0;
			for (Int_t it = 0; it < 2; ++it) { ccsum += ccelm[it] * ccwgt[it]; ccsgm += ccwgt[it]; }
			ccest = ccsum / ccsgm;
		}
		
		Hist::Head(StrFmt("hIi_CCest%s", trNm.c_str()))->fill(irig, ccest, weight);
		if (sign>0) Hist::Head(StrFmt("hIp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hIn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_CCest%s%03d", trNm.c_str(), tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_CCest%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_CCest%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		
		Float_t CCcutFs = 1.75 * std::erf((std::log( 100.0) - std::log(nrig)) / 1.7); // 100.   80.
		Float_t CCcutL9 = 1.75 * std::erf((std::log(  75.0) - std::log(nrig)) / 1.7); // 75.    65.
		Float_t CCcutL1 = 1.75 * std::erf((std::log(  75.0) - std::log(nrig)) / 1.7); // 75.    65.
		Float_t CCcutIn = 1.75 * std::erf((std::log(  40.0) - std::log(nrig)) / 1.7); // 50.    40.
		if (trPt == 5 && ccest > CCcutFs) break;
		if (trPt == 4 && ccest > CCcutL9) break;
		if (trPt == 3 && ccest > CCcutL1) break;
		if (trPt == 2 && ccest > CCcutIn) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 4, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 4, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 4, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 4, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		
		// RICH
		Int_t   richph    = fRich->numOfHit;
		Float_t richelph  = (fRich->kindOfRad == -1) ? 0.0 : fRich->numOfExpPE[0];
		Float_t richprph  = (fRich->kindOfRad == -1) ? 0.0 : fRich->numOfExpPE[3];
		
		Bool_t  hasrich   = (fRich->status && fRich->kindOfRad == 0);
		Float_t richbPr   = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbtaPr);
		Float_t richbPi   = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbtaPi);
		Float_t richbEl   = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbtaEl);
		Float_t tofrichb  = (!hasrich) ? 0.0 : (fTof->betaH - fRich->beta);

		if (fRich->kindOfRad == -1) break;
		
		Float_t richbth   = (fRich->kindOfRad == 0) ? 0.003 : 0.005;
		Float_t richbsh   = std::sqrt(1.0+(0.4/nrig)*(0.4/nrig))-1.0;
		Bool_t  isrichNo  = (!fRich->status && fRich->kindOfRad == 0 && MgntNum::Compare(richelph, 2.0f) > 0 && MgntNum::Compare(richprph, 0.5f) < 0 && richph <= 2);
		Bool_t  isrichPr  = ( fRich->status && std::fabs(richbPr) < (richbth + richbsh));
		Bool_t  isrichPi  = ( fRich->status && std::fabs(richbPi) < richbth);
		Bool_t  isrichEl  = ( fRich->status && std::fabs(richbEl) < richbth);

		Hist::Head("hIi_Cutflow")->fill(irig, 5, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 5, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 5, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 5, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		
		if (hasrich) Hist::Head("hIi_Aglb")->fill(irig, richbPr, weight);
		if (hasrich && sign>0) Hist::Head("hIp_Aglb")->fill(nrig, richbPr, weight);
		if (hasrich && sign<0) Hist::Head("hIn_Aglb")->fill(nrig, richbPr, weight);
		if (hasrich && tvBin!=0) Hist::Head(StrFmt("hIi_Aglb%03d", tvBin))->fill(irig, richbPr, weight);
		if (hasrich && tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Aglb%03d", tvBin))->fill(nrig, richbPr, weight);
		if (hasrich && tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Aglb%03d", tvBin))->fill(nrig, richbPr, weight);
		
		if (hasrich && !isrichPr) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 6, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 6, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 6, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 6, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 6, weight);

		const Double_t lmtPiEl = 9.26;
		if (nrig < lmtPiEl && !hasrich) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 7, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 7, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 7, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 7, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 7, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 7, weight);

		// Sp Rig	
		Float_t trigSp  = track->rigidity[0][2];
		Int_t   signSp  = MgntNum::Compare(trigSp);
		Float_t nrigSp  = std::fabs(trigSp);
		Float_t irigSp  = 1. / trigSp;
		
		// ECAL
		const Double_t lmtPiEl2 = 9.26;
		Bool_t isSignal     = (sign > 0 && isecalPr);
		Bool_t isBackground = (sign < 0 && isecalEl);
		if (isSignal)     Hist::Head("hIs_Trdl")->fill(nrig, trdl, weight);
		if (isBackground) Hist::Head("hIb_Trdl")->fill(nrig, trdl, weight);
		if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		
		if (isSignal)     Hist::Head(StrFmt("hIs_Trdl%s", trNm.c_str()))->fill(nrig, trdl, weight);
		if (isBackground) Hist::Head(StrFmt("hIb_Trdl%s", trNm.c_str()))->fill(nrig, trdl, weight);
		if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_Trdl%s%03d", trNm.c_str(), tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_Trdl%s%03d", trNm.c_str(), tvBin))->fill(nrig, trdl, weight);
		
		if (isSignal)     Hist::Head("hIs_TrdlSp")->fill(nrigSp, trdl, weight);
		if (isBackground) Hist::Head("hIb_TrdlSp")->fill(nrigSp, trdl, weight);
		if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_TrdlSp%03d", tvBin))->fill(nrigSp, trdl, weight);
		if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_TrdlSp%03d", tvBin))->fill(nrigSp, trdl, weight);

		if (hasecal && !isecalPr) break;

		Hist::Head("hIi_Cutflow")->fill(irig, 8, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 8, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 8, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 8, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 8, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 8, weight);

		Hist::Head("hIi_Trdl")->fill(irig, trdl, weight);
		if (sign>0) Hist::Head("hIp_Trdl")->fill(nrig, trdl, weight);
		if (sign<0) Hist::Head("hIn_Trdl")->fill(nrig, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		
		Hist::Head(StrFmt("hIi_Trdl%s", trNm.c_str()))->fill(irig, trdl, weight);
		if (sign>0) Hist::Head(StrFmt("hIp_Trdl%s", trNm.c_str()))->fill(nrig, trdl, weight);
		if (sign<0) Hist::Head(StrFmt("hIn_Trdl%s", trNm.c_str()))->fill(nrig, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Trdl%s%03d", trNm.c_str(), tvBin))->fill(irig, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Trdl%s%03d", trNm.c_str(), tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Trdl%s%03d", trNm.c_str(), tvBin))->fill(nrig, trdl, weight);
		
		Hist::Head("hIi_TrdlSp")->fill(irigSp, trdl, weight);
		if (sign>0) Hist::Head("hIp_TrdlSp")->fill(nrigSp, trdl, weight);
		if (sign<0) Hist::Head("hIn_TrdlSp")->fill(nrigSp, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_TrdlSp%03d", tvBin))->fill(irigSp, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_TrdlSp%03d", tvBin))->fill(nrigSp, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_TrdlSp%03d", tvBin))->fill(nrigSp, trdl, weight);
		
		Hist::Head("hIi_Evt")->fill(irig, weight);
		if (sign>0) Hist::Head("hIp_Evt")->fill(nrig, weight);
		if (sign<0) Hist::Head("hIn_Evt")->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Evt%03d", tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Evt%03d", tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Evt%03d", tvBin))->fill(nrig, weight);
		
		Hist::Head(StrFmt("hIi_Evt%s", trNm.c_str()))->fill(irig, weight);
		if (sign>0) Hist::Head(StrFmt("hIp_Evt%s", trNm.c_str()))->fill(nrig, weight);
		if (sign<0) Hist::Head(StrFmt("hIn_Evt%s", trNm.c_str()))->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Evt%s%03d", trNm.c_str(), tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Evt%s%03d", trNm.c_str(), tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Evt%s%03d", trNm.c_str(), tvBin))->fill(nrig, weight);
		
		Hist::Head("hIi_EvtSp")->fill(irigSp, weight);
		if (sign>0) Hist::Head("hIp_EvtSp")->fill(nrigSp, weight);
		if (sign<0) Hist::Head("hIn_EvtSp")->fill(nrigSp, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_EvtSp%03d", tvBin))->fill(irigSp, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_EvtSp%03d", tvBin))->fill(nrigSp, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_EvtSp%03d", tvBin))->fill(nrigSp, weight);
		
		if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
			Hist::Head("hIi_MCTrRso")->fill(MCIRig, irig, weight);
			if (MCTRig>0.0) Hist::Head("hIp_MCTrRso")->fill(MCNRig, (irig - MCIRig), weight);
			if (MCTRig<0.0) Hist::Head("hIn_MCTrRso")->fill(MCNRig, (irig - MCIRig), weight);
			
			Hist::Head(StrFmt("hIi_MCTrRso%s", trNm.c_str()))->fill(MCIRig, irig, weight);
			if (MCTRig>0.0) Hist::Head(StrFmt("hIp_MCTrRso%s", trNm.c_str()))->fill(MCNRig, (irig - MCIRig), weight);
			if (MCTRig<0.0) Hist::Head(StrFmt("hIn_MCTrRso%s", trNm.c_str()))->fill(MCNRig, (irig - MCIRig), weight);
		
			Hist::Head("hIi_MCTrRsoSp")->fill(MCIRig, irigSp, weight);
			if (MCTRig>0.0) Hist::Head("hIp_MCTrRsoSp")->fill(MCNRig, (irigSp - MCIRig), weight);
			if (MCTRig<0.0) Hist::Head("hIn_MCTrRsoSp")->fill(MCNRig, (irigSp - MCIRig), weight);
			
			Hist::Head("hIi_MCEvt")->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head("hIp_MCEvt")->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head("hIn_MCEvt")->fill(MCNRig, weight);
			
			Hist::Head(StrFmt("hIi_MCEvt%s", trNm.c_str()))->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head(StrFmt("hIp_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head(StrFmt("hIn_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
			
			Hist::Head("hIi_MCEvtSp")->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head("hIp_MCEvtSp")->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head("hIn_MCEvtSp")->fill(MCNRig, weight);
		}

		break;
	}

#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] H Selection [][]\n";
#endif
	/**** High Energy Region ****/
	while (true) {
		// TRK
		Int_t trPt = -1;
		std::string trNm = "";
		if      (track->status[0][5] && (track->bitPatt&161)==161) { trPt = 5; trNm = "Fs"; }
		//else if (track->status[0][4] && (track->bitPatt&129)==129) { trPt = 4; trNm = "L9"; }
		else if (track->status[0][4] && (track->bitPatt&137)==137) { trPt = 4; trNm = "L9"; }
		else if (track->status[0][3] && (track->bitPatt& 33)== 33) { trPt = 3; trNm = "L1"; }
		if (trPt == -1) break;
		Float_t trig = track->rigidity[0][trPt];
		Int_t   sign = MgntNum::Compare(trig);
		Float_t nrig = std::fabs(trig);
		Float_t irig = 1. / trig;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		
		Int_t   patAgl = -1;
		Float_t cosAgl = std::acos(-track->stateL1[0][trPt][5]) * TMath::RadToDeg();
		if      (cosAgl < 25.0) patAgl = 0; 
		else if (cosAgl < 30.0) patAgl = 1; 
		else if (cosAgl < 35.0) patAgl = 2; 
		else if (cosAgl < 40.0) patAgl = 3; 
		if (patAgl == -1) break;
		
		// Cutoff
		if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
			Double_t cfRig = DST::CfStableFT * fRti->cutoffIGRF[patAgl];
			//Double_t cfRig = DST::CfStableFT * fRti->cutoffIGRF[3];
			if (nrig < cfRig) break;
		}
		
		Hist::Head("hHi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		// TRD
		Int_t trdPt = 1;
		Bool_t  hastrd  = (fTrd->statusKCls[trdPt] && fTrd->LLR_nhit[trdPt] >= 8);
		if (!hastrd) break;
		
		Float_t trdl    = (fTrd->LLR[trdPt][0]);
		Bool_t  istrdHe = (fTrd->LLR[trdPt][2] > 0.3);
		Bool_t  istrdPr = (trdl > 0.75);
		if (istrdHe) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		
		Hist::Head("hHi_Trdl")->fill(irig, trdl, weight);
		if (sign>0) Hist::Head("hHp_Trdl")->fill(nrig, trdl, weight);
		if (sign<0) Hist::Head("hHn_Trdl")->fill(nrig, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		
		if (!istrdPr) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);

		// ECAL
		Bool_t  hasecal  = (shower != nullptr);
		Float_t ecalbdt  = (!hasecal) ? -2.0 : shower->PisaBDT;
		Bool_t  isecalPr = (!hasecal) ? false : (ecalbdt < -0.8);
		Bool_t  isecalEl = (!hasecal) ? false : (ecalbdt >  0.6);

		if (hasecal) Hist::Head("hHi_Ecalbdt")->fill(irig, ecalbdt, weight);
		if (hasecal && sign>0) Hist::Head("hHp_Ecalbdt")->fill(nrig, ecalbdt, weight);
		if (hasecal && sign<0) Hist::Head("hHn_Ecalbdt")->fill(nrig, ecalbdt, weight);
		if (tvBin!=0 && hasecal) Hist::Head(StrFmt("hHi_Ecalbdt%03d", tvBin))->fill(irig, ecalbdt, weight);
		if (tvBin!=0 && hasecal && sign>0) Hist::Head(StrFmt("hHp_Ecalbdt%03d", tvBin))->fill(nrig, ecalbdt, weight);
		if (tvBin!=0 && hasecal && sign<0) Hist::Head(StrFmt("hHn_Ecalbdt%03d", tvBin))->fill(nrig, ecalbdt, weight);
		
		if (hasecal && !isecalPr) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 4, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 4, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 4, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 4, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		
		// Charge Confusion
		const Float_t lasymLMT = 1.0e-8;
		const Float_t lasymSGM = 1.8;
		Float_t trSgmIn = GetRigSgm(2, nrig, massPr);
		Float_t trSgmL1 = GetRigSgm(3, nrig, massPr);
		Float_t trSgmL9 = GetRigSgm(4, nrig, massPr);
		Float_t trSgmFs = GetRigSgm(5, nrig, massPr);
		
		Float_t trSgmPt = 1.0 / nrig;
		if      (trPt == 3) trSgmPt = trSgmL1;
		else if (trPt == 4) trSgmPt = trSgmL9;
		else if (trPt == 5) trSgmPt = trSgmFs;
		
		Bool_t  hasasym1I   = (track->status[0][2] && track->status[0][3]);
		Float_t crrPar1I[5] = { 4.93251e-01, 1.90417e+00, 9.18502e+01, 4.86389e+00, 3.29747e-01 };
		Float_t crrSgm1I    = crrPar1I[0]*std::erf(crrPar1I[1]*(std::log(nrig+crrPar1I[2])-crrPar1I[3]))+crrPar1I[4];
		Float_t trSgm1I     = (!hasasym1I) ? 0.0 : std::sqrt(trSgmIn * trSgmIn + trSgmL1 * trSgmL1) * crrSgm1I;
		Float_t asym1I      = (!hasasym1I) ? 0.0 : (1.0/track->rigidity[0][2] - 1.0/track->rigidity[0][3]) / trSgm1I; 
		Float_t lasym1I     = (!hasasym1I) ? 0.0 : std::log(asym1I * asym1I + lasymLMT) / lasymSGM;

		Bool_t  hasasym9I   = (track->status[0][2] && track->status[0][4]);
		Float_t crrPar9I[5] = { 4.17294e-01, 1.15603e+00, 1.62194e+01, 3.81780e+00, 3.92404e-01 };
		Float_t crrSgm9I    = crrPar9I[0]*std::erf(crrPar9I[1]*(std::log(nrig+crrPar9I[2])-crrPar9I[3]))+crrPar9I[4];
		Float_t trSgm9I     = (!hasasym9I) ? 0.0 : std::sqrt(trSgmIn * trSgmIn + trSgmL9 * trSgmL9) * crrSgm9I;
		Float_t asym9I      = (!hasasym9I) ? 0.0 : (1.0/track->rigidity[0][2] - 1.0/track->rigidity[0][4]) / trSgm9I;
		Float_t lasym9I     = (!hasasym9I) ? 0.0 : std::log(asym9I * asym9I + lasymLMT) / lasymSGM;

		Bool_t  hasasym91   = (track->status[0][3] && track->status[0][4]);
		Float_t crrPar91[5] = { 5.25352e-01, 9.88059e-01, 1.21200e+01, 3.72027e+00, 4.99424e-01 };
		Float_t crrSgm91    = crrPar91[0]*std::erf(crrPar91[1]*(std::log(nrig+crrPar91[2])-crrPar91[3]))+crrPar91[4];
		Float_t trSgm91     = (!hasasym91) ? 0.0 : std::sqrt(trSgmL1 * trSgmL1 + trSgmL9 * trSgmL9) * crrSgm91;
		Float_t asym91      = (!hasasym91) ? 0.0 : (1.0/track->rigidity[0][3] - 1.0/track->rigidity[0][4]) / trSgm91;
		Float_t lasym91     = (!hasasym91) ? 0.0 : std::log(asym91 * asym91 + lasymLMT) / lasymSGM;

		if (hasasym1I) Hist::Head("hHi_Lasym1I")->fill(irig, lasym1I, weight);
		if (hasasym1I && sign>0) Hist::Head("hHp_Lasym1I")->fill(nrig, lasym1I, weight);
		if (hasasym1I && sign<0) Hist::Head("hHn_Lasym1I")->fill(nrig, lasym1I, weight);
		if (hasasym9I) Hist::Head("hHi_Lasym9I")->fill(irig, lasym9I, weight);
		if (hasasym9I && sign>0) Hist::Head("hHp_Lasym9I")->fill(nrig, lasym9I, weight);
		if (hasasym9I && sign<0) Hist::Head("hHn_Lasym9I")->fill(nrig, lasym9I, weight);
		if (hasasym91) Hist::Head("hHi_Lasym91")->fill(irig, lasym91, weight);
		if (hasasym91 && sign>0) Hist::Head("hHp_Lasym91")->fill(nrig, lasym91, weight);
		if (hasasym91 && sign<0) Hist::Head("hHn_Lasym91")->fill(nrig, lasym91, weight);
		
		if (tvBin!=0 && hasasym1I) Hist::Head(StrFmt("hHi_Lasym1I%03d", tvBin))->fill(irig, lasym1I, weight);
		if (tvBin!=0 && hasasym1I && sign>0) Hist::Head(StrFmt("hHp_Lasym1I%03d", tvBin))->fill(nrig, lasym1I, weight);
		if (tvBin!=0 && hasasym1I && sign<0) Hist::Head(StrFmt("hHn_Lasym1I%03d", tvBin))->fill(nrig, lasym1I, weight);
		if (tvBin!=0 && hasasym9I) Hist::Head(StrFmt("hHi_Lasym9I%03d", tvBin))->fill(irig, lasym9I, weight);
		if (tvBin!=0 && hasasym9I && sign>0) Hist::Head(StrFmt("hHp_Lasym9I%03d", tvBin))->fill(nrig, lasym9I, weight);
		if (tvBin!=0 && hasasym9I && sign<0) Hist::Head(StrFmt("hHn_Lasym9I%03d", tvBin))->fill(nrig, lasym9I, weight);
		if (tvBin!=0 && hasasym91) Hist::Head(StrFmt("hHi_Lasym91%03d", tvBin))->fill(irig, lasym91, weight);
		if (tvBin!=0 && hasasym91 && sign>0) Hist::Head(StrFmt("hHp_Lasym91%03d", tvBin))->fill(nrig, lasym91, weight);
		if (tvBin!=0 && hasasym91 && sign<0) Hist::Head(StrFmt("hHn_Lasym91%03d", tvBin))->fill(nrig, lasym91, weight);
			
		Float_t lchixL1 = (!track->status[0][3]) ? 0.0 : std::log(track->chisq[0][3][0]);
		Float_t lchiyL1 = (!track->status[0][3]) ? 0.0 : std::log(track->chisq[0][3][1]);
		Float_t lchixL9 = (!track->status[0][4]) ? 0.0 : std::log(track->chisq[0][4][0]);
		Float_t lchiyL9 = (!track->status[0][4]) ? 0.0 : std::log(track->chisq[0][4][1]);
		Float_t lchixFs = (!track->status[0][5]) ? 0.0 : std::log(track->chisq[0][5][0]);
		Float_t lchiyFs = (!track->status[0][5]) ? 0.0 : std::log(track->chisq[0][5][1]);
		
		if (track->status[0][3]) Hist::Head("hHi_LchixL1")->fill(irig, lchixL1, weight);
		if (track->status[0][3] && sign>0) Hist::Head("hHp_LchixL1")->fill(nrig, lchixL1, weight);
		if (track->status[0][3] && sign<0) Hist::Head("hHn_LchixL1")->fill(nrig, lchixL1, weight);
		if (track->status[0][3]) Hist::Head("hHi_LchiyL1")->fill(irig, lchiyL1, weight);
		if (track->status[0][3] && sign>0) Hist::Head("hHp_LchiyL1")->fill(nrig, lchiyL1, weight);
		if (track->status[0][3] && sign<0) Hist::Head("hHn_LchiyL1")->fill(nrig, lchiyL1, weight);
		if (track->status[0][4]) Hist::Head("hHi_LchixL9")->fill(irig, lchixL9, weight);
		if (track->status[0][4] && sign>0) Hist::Head("hHp_LchixL9")->fill(nrig, lchixL9, weight);
		if (track->status[0][4] && sign<0) Hist::Head("hHn_LchixL9")->fill(nrig, lchixL9, weight);
		if (track->status[0][4]) Hist::Head("hHi_LchiyL9")->fill(irig, lchiyL9, weight);
		if (track->status[0][4] && sign>0) Hist::Head("hHp_LchiyL9")->fill(nrig, lchiyL9, weight);
		if (track->status[0][4] && sign<0) Hist::Head("hHn_LchiyL9")->fill(nrig, lchiyL9, weight);
		if (track->status[0][5]) Hist::Head("hHi_LchixFs")->fill(irig, lchixFs, weight);
		if (track->status[0][5] && sign>0) Hist::Head("hHp_LchixFs")->fill(nrig, lchixFs, weight);
		if (track->status[0][5] && sign<0) Hist::Head("hHn_LchixFs")->fill(nrig, lchixFs, weight);
		if (track->status[0][5]) Hist::Head("hHi_LchiyFs")->fill(irig, lchiyFs, weight);
		if (track->status[0][5] && sign>0) Hist::Head("hHp_LchiyFs")->fill(nrig, lchiyFs, weight);
		if (track->status[0][5] && sign<0) Hist::Head("hHn_LchiyFs")->fill(nrig, lchiyFs, weight);
		
		if (tvBin!=0 && track->status[0][3]) Hist::Head(StrFmt("hHi_LchixL1%03d", tvBin))->fill(irig, lchixL1, weight);
		if (tvBin!=0 && track->status[0][3] && sign>0) Hist::Head(StrFmt("hHp_LchixL1%03d", tvBin))->fill(nrig, lchixL1, weight);
		if (tvBin!=0 && track->status[0][3] && sign<0) Hist::Head(StrFmt("hHn_LchixL1%03d", tvBin))->fill(nrig, lchixL1, weight);
		if (tvBin!=0 && track->status[0][3]) Hist::Head(StrFmt("hHi_LchiyL1%03d", tvBin))->fill(irig, lchiyL1, weight);
		if (tvBin!=0 && track->status[0][3] && sign>0) Hist::Head(StrFmt("hHp_LchiyL1%03d", tvBin))->fill(nrig, lchiyL1, weight);
		if (tvBin!=0 && track->status[0][3] && sign<0) Hist::Head(StrFmt("hHn_LchiyL1%03d", tvBin))->fill(nrig, lchiyL1, weight);
		if (tvBin!=0 && track->status[0][4]) Hist::Head(StrFmt("hHi_LchixL9%03d", tvBin))->fill(irig, lchixL9, weight);
		if (tvBin!=0 && track->status[0][4] && sign>0) Hist::Head(StrFmt("hHp_LchixL9%03d", tvBin))->fill(nrig, lchixL9, weight);
		if (tvBin!=0 && track->status[0][4] && sign<0) Hist::Head(StrFmt("hHn_LchixL9%03d", tvBin))->fill(nrig, lchixL9, weight);
		if (tvBin!=0 && track->status[0][4]) Hist::Head(StrFmt("hHi_LchiyL9%03d", tvBin))->fill(irig, lchiyL9, weight);
		if (tvBin!=0 && track->status[0][4] && sign>0) Hist::Head(StrFmt("hHp_LchiyL9%03d", tvBin))->fill(nrig, lchiyL9, weight);
		if (tvBin!=0 && track->status[0][4] && sign<0) Hist::Head(StrFmt("hHn_LchiyL9%03d", tvBin))->fill(nrig, lchiyL9, weight);
		if (tvBin!=0 && track->status[0][5]) Hist::Head(StrFmt("hHi_LchixFs%03d", tvBin))->fill(irig, lchixFs, weight);
		if (tvBin!=0 && track->status[0][5] && sign>0) Hist::Head(StrFmt("hHp_LchixFs%03d", tvBin))->fill(nrig, lchixFs, weight);
		if (tvBin!=0 && track->status[0][5] && sign<0) Hist::Head(StrFmt("hHn_LchixFs%03d", tvBin))->fill(nrig, lchixFs, weight);
		if (tvBin!=0 && track->status[0][5]) Hist::Head(StrFmt("hHi_LchiyFs%03d", tvBin))->fill(irig, lchiyFs, weight);
		if (tvBin!=0 && track->status[0][5] && sign>0) Hist::Head(StrFmt("hHp_LchiyFs%03d", tvBin))->fill(nrig, lchiyFs, weight);
		if (tvBin!=0 && track->status[0][5] && sign<0) Hist::Head(StrFmt("hHn_LchiyFs%03d", tvBin))->fill(nrig, lchiyFs, weight);
		
		// Res
		const Float_t lresLMT = 1.0e-6;
		const Float_t lresSGM = 1.6;
		
		Float_t trL1SclxL1  = 1.0e4 / (2.800752e-01 * std::sqrt(4.52677e+08/nrig/nrig + 1.0));
		Float_t trL1SclyL1  = 1.0e4 / (5.664780e-00 * (std::sqrt(1.12235e+06/nrig/nrig + 1.0) * std::exp(-7.89344e-03*std::pow(nrig,1.07759e+00)) + 1.0));
		Float_t trL1ResxL1  = (trPt == 3) ? trL1SclxL1 * (track->stateL1[0][trPt][0] - trHitL1->coo[0]) : 0.0; 
		Float_t trL1ResyL1  = (trPt == 3) ? trL1SclyL1 * (track->stateL1[0][trPt][1] - trHitL1->coo[1]) : 0.0;
		Float_t trL1LresxL1 = (trPt == 3) ? std::log(trL1ResxL1 * trL1ResxL1 + lresLMT) / lresSGM : 0.0; 
		Float_t trL1LresyL1 = (trPt == 3) ? std::log(trL1ResyL1 * trL1ResyL1 + lresLMT) / lresSGM : 0.0;
		
		Float_t trL9SclxL9  = 1.0e4 / (3.834576e+00 * std::sqrt(1.39170e+06/nrig/nrig + 1.0));
		Float_t trL9SclyL9  = 1.0e4 / (5.010260e-00 * (std::sqrt(1.31972e+06/nrig/nrig + 1.0) * std::exp(-5.23536e-02*std::pow(nrig,7.89620e-01)) + 1.0));
		Float_t trL9ResxL9  = (trPt == 4) ? trL9SclxL9 * (track->stateL9[0][trPt][0] - trHitL9->coo[0]) : 0.0;
		Float_t trL9ResyL9  = (trPt == 4) ? trL9SclyL9 * (track->stateL9[0][trPt][1] - trHitL9->coo[1]) : 0.0;
		Float_t trL9LresxL9 = (trPt == 4) ? std::log(trL9ResxL9 * trL9ResxL9 + lresLMT) / lresSGM : 0.0;
		Float_t trL9LresyL9 = (trPt == 4) ? std::log(trL9ResyL9 * trL9ResyL9 + lresLMT) / lresSGM : 0.0;
		
		Float_t trFsSclxL1  = 1.0e4 / (2.264057e+01 * std::sqrt(6.45000e+04/nrig/nrig + 1.0));
		Float_t trFsSclyL1  = 1.0e4 / (7.280530e-00 * (std::sqrt(4.86353e+05/nrig/nrig + 1.0) * std::exp(-1.40328e-03*std::pow(nrig,1.32389e+00)) + 1.0));
		Float_t trFsResxL1  = (trPt == 5) ? trFsSclxL1 * (track->stateL1[0][trPt][0] - trHitL1->coo[0]) : 0.0;
		Float_t trFsResyL1  = (trPt == 5) ? trFsSclyL1 * (track->stateL1[0][trPt][1] - trHitL1->coo[1]) : 0.0;
		Float_t trFsLresxL1 = (trPt == 5) ? std::log(trFsResxL1 * trFsResxL1 + lresLMT) / lresSGM : 0.0;
		Float_t trFsLresyL1 = (trPt == 5) ? std::log(trFsResyL1 * trFsResyL1 + lresLMT) / lresSGM : 0.0;
		
		Float_t trFsSclxL9  = 1.0e4 / (1.949518e+01 * std::sqrt(5.14464e+04/nrig/nrig + 1.0));
		Float_t trFsSclyL9  = 1.0e4 / (5.966690e-00 * (std::sqrt(6.82955e+05/nrig/nrig + 1.0) * std::exp(-4.23847e-02*std::pow(nrig,7.56334e-01)) + 1.0));
		Float_t trFsResxL9  = (trPt == 5) ? trFsSclxL9 * (track->stateL9[0][trPt][0] - trHitL9->coo[0]) : 0.0;
		Float_t trFsResyL9  = (trPt == 5) ? trFsSclyL9 * (track->stateL9[0][trPt][1] - trHitL9->coo[1]) : 0.0;
		Float_t trFsLresxL9 = (trPt == 5) ? std::log(trFsResxL9 * trFsResxL9 + lresLMT) / lresSGM : 0.0;
		Float_t trFsLresyL9 = (trPt == 5) ? std::log(trFsResyL9 * trFsResyL9 + lresLMT) / lresSGM : 0.0;
		
		if (trPt==3) Hist::Head("hHi_L1LresxL1")->fill(irig, trL1LresxL1, weight);
		if (trPt==3 && sign>0) Hist::Head("hHp_L1LresxL1")->fill(nrig, trL1LresxL1, weight);
		if (trPt==3 && sign<0) Hist::Head("hHn_L1LresxL1")->fill(nrig, trL1LresxL1, weight);
		if (trPt==3 && tvBin!=0) Hist::Head(StrFmt("hHi_L1LresxL1%03d", tvBin))->fill(irig, trL1LresxL1, weight);
		if (trPt==3 && tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_L1LresxL1%03d", tvBin))->fill(nrig, trL1LresxL1, weight);
		if (trPt==3 && tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_L1LresxL1%03d", tvBin))->fill(nrig, trL1LresxL1, weight);
		
		if (trPt==3) Hist::Head("hHi_L1LresyL1")->fill(irig, trL1LresyL1, weight);
		if (trPt==3 && sign>0) Hist::Head("hHp_L1LresyL1")->fill(nrig, trL1LresyL1, weight);
		if (trPt==3 && sign<0) Hist::Head("hHn_L1LresyL1")->fill(nrig, trL1LresyL1, weight);
		if (trPt==3 && tvBin!=0) Hist::Head(StrFmt("hHi_L1LresyL1%03d", tvBin))->fill(irig, trL1LresyL1, weight);
		if (trPt==3 && tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_L1LresyL1%03d", tvBin))->fill(nrig, trL1LresyL1, weight);
		if (trPt==3 && tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_L1LresyL1%03d", tvBin))->fill(nrig, trL1LresyL1, weight);
		
		if (trPt==4) Hist::Head("hHi_L9LresxL9")->fill(irig, trL9LresxL9, weight);
		if (trPt==4 && sign>0) Hist::Head("hHp_L9LresxL9")->fill(nrig, trL9LresxL9, weight);
		if (trPt==4 && sign<0) Hist::Head("hHn_L9LresxL9")->fill(nrig, trL9LresxL9, weight);
		if (trPt==4 && tvBin!=0) Hist::Head(StrFmt("hHi_L9LresxL9%03d", tvBin))->fill(irig, trL9LresxL9, weight);
		if (trPt==4 && tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_L9LresxL9%03d", tvBin))->fill(nrig, trL9LresxL9, weight);
		if (trPt==4 && tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_L9LresxL9%03d", tvBin))->fill(nrig, trL9LresxL9, weight);
		
		if (trPt==4) Hist::Head("hHi_L9LresyL9")->fill(irig, trL9LresyL9, weight);
		if (trPt==4 && sign>0) Hist::Head("hHp_L9LresyL9")->fill(nrig, trL9LresyL9, weight);
		if (trPt==4 && sign<0) Hist::Head("hHn_L9LresyL9")->fill(nrig, trL9LresyL9, weight);
		if (trPt==4 && tvBin!=0) Hist::Head(StrFmt("hHi_L9LresyL9%03d", tvBin))->fill(irig, trL9LresyL9, weight);
		if (trPt==4 && tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_L9LresyL9%03d", tvBin))->fill(nrig, trL9LresyL9, weight);
		if (trPt==4 && tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_L9LresyL9%03d", tvBin))->fill(nrig, trL9LresyL9, weight);
		
		if (trPt==5) Hist::Head("hHi_FsLresxL1")->fill(irig, trFsLresxL1, weight);
		if (trPt==5 && sign>0) Hist::Head("hHp_FsLresxL1")->fill(nrig, trFsLresxL1, weight);
		if (trPt==5 && sign<0) Hist::Head("hHn_FsLresxL1")->fill(nrig, trFsLresxL1, weight);
		if (trPt==5 && tvBin!=0) Hist::Head(StrFmt("hHi_FsLresxL1%03d", tvBin))->fill(irig, trFsLresxL1, weight);
		if (trPt==5 && tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_FsLresxL1%03d", tvBin))->fill(nrig, trFsLresxL1, weight);
		if (trPt==5 && tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_FsLresxL1%03d", tvBin))->fill(nrig, trFsLresxL1, weight);
		
		if (trPt==5) Hist::Head("hHi_FsLresyL1")->fill(irig, trFsLresyL1, weight);
		if (trPt==5 && sign>0) Hist::Head("hHp_FsLresyL1")->fill(nrig, trFsLresyL1, weight);
		if (trPt==5 && sign<0) Hist::Head("hHn_FsLresyL1")->fill(nrig, trFsLresyL1, weight);
		if (trPt==5 && tvBin!=0) Hist::Head(StrFmt("hHi_FsLresyL1%03d", tvBin))->fill(irig, trFsLresyL1, weight);
		if (trPt==5 && tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_FsLresyL1%03d", tvBin))->fill(nrig, trFsLresyL1, weight);
		if (trPt==5 && tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_FsLresyL1%03d", tvBin))->fill(nrig, trFsLresyL1, weight);
		
		if (trPt==5) Hist::Head("hHi_FsLresxL9")->fill(irig, trFsLresxL9, weight);
		if (trPt==5 && sign>0) Hist::Head("hHp_FsLresxL9")->fill(nrig, trFsLresxL9, weight);
		if (trPt==5 && sign<0) Hist::Head("hHn_FsLresxL9")->fill(nrig, trFsLresxL9, weight);
		if (trPt==5 && tvBin!=0) Hist::Head(StrFmt("hHi_FsLresxL9%03d", tvBin))->fill(irig, trFsLresxL9, weight);
		if (trPt==5 && tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_FsLresxL9%03d", tvBin))->fill(nrig, trFsLresxL9, weight);
		if (trPt==5 && tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_FsLresxL9%03d", tvBin))->fill(nrig, trFsLresxL9, weight);
		
		if (trPt==5) Hist::Head("hHi_FsLresyL9")->fill(irig, trFsLresyL9, weight);
		if (trPt==5 && sign>0) Hist::Head("hHp_FsLresyL9")->fill(nrig, trFsLresyL9, weight);
		if (trPt==5 && sign<0) Hist::Head("hHn_FsLresyL9")->fill(nrig, trFsLresyL9, weight);
		if (trPt==5 && tvBin!=0) Hist::Head(StrFmt("hHi_FsLresyL9%03d", tvBin))->fill(irig, trFsLresyL9, weight);
		if (trPt==5 && tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_FsLresyL9%03d", tvBin))->fill(nrig, trFsLresyL9, weight);
		if (trPt==5 && tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_FsLresyL9%03d", tvBin))->fill(nrig, trFsLresyL9, weight);
		
		// New
		if (trPt == 5 && (lchixL1 > 3.5 || lchixL9 > 3.5 || lchixFs > 3.5 || trFsLresxL1 > 2.0 || trFsLresxL9 > 2.0)) break;
		if (trPt == 4 && (lchixL9 > 3.5 || trL9LresxL9 > 2.0)) break;
		if (trPt == 3 && (lchixL1 > 3.5 || trL1LresxL1 > 2.0)) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 5, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 5, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 5, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 5, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 5, weight);

		// Charge Confusion
		Float_t ccest = 0.;
		if (trPt == 5) {
			//const Int_t ncc = 8;
			//Float_t ccelm[ncc] = { lchiyL1, lchiyL9, lchiyFs, lasym1I, lasym9I, lasym91, trFsLresyL1, trFsLresyL9 };
			//Float_t ccwgt[ncc] = { GetSqrMDR(3), GetSqrMDR(4), GetSqrMDR(5), 0.5*(GetSqrMDR(2)+GetSqrMDR(3)), 0.5*(GetSqrMDR(2)+GetSqrMDR(4)), 0.5*(GetSqrMDR(4)+GetSqrMDR(5)), 0.5*GetSqrMDR(5), 0.5*GetSqrMDR(5) };
			const Int_t ncc = 6;
			Float_t ccelm[ncc] = { lchiyL1, lchiyL9, lchiyFs, lasym1I, lasym9I, lasym91 };
			Float_t ccwgt[ncc] = { GetSqrMDR(3), GetSqrMDR(4), GetSqrMDR(5), 0.5*(GetSqrMDR(2)+GetSqrMDR(3)), 0.5*(GetSqrMDR(2)+GetSqrMDR(4)), 0.5*(GetSqrMDR(4)+GetSqrMDR(5)) };
			Float_t ccsum = 0.0, ccsgm = 0.0;
			for (Int_t it = 0; it < ncc; ++it) { ccsum += ccelm[it] * ccwgt[it]; ccsgm += ccwgt[it]; }
			ccest = ccsum / ccsgm;
		}
		if (trPt == 4) {
			//const Int_t ncc = 3;
			//Float_t ccelm[ncc] = { lchiyL9, lasym9I, trL9LresyL9 };
			//Float_t ccwgt[ncc] = { GetSqrMDR(4), 0.5*(GetSqrMDR(2)+GetSqrMDR(4)), 0.5*GetSqrMDR(4) };
			const Int_t ncc = 2;
			Float_t ccelm[ncc] = { lchiyL9, lasym9I };
			Float_t ccwgt[ncc] = { GetSqrMDR(4), 0.5*(GetSqrMDR(2)+GetSqrMDR(4)) };
			Float_t ccsum = 0.0, ccsgm = 0.0;
			for (Int_t it = 0; it < ncc; ++it) { ccsum += ccelm[it] * ccwgt[it]; ccsgm += ccwgt[it]; }
			ccest = ccsum / ccsgm;
		}
		if (trPt == 3) {
			//const Int_t ncc = 3;
			//Float_t ccelm[ncc] = { lchiyL1, lasym1I, trL1LresyL1 };
			//Float_t ccwgt[ncc] = { GetSqrMDR(3), 0.5*(GetSqrMDR(2)+GetSqrMDR(3)), 0.5*GetSqrMDR(3) };
			const Int_t ncc = 2;
			Float_t ccelm[ncc] = { lchiyL1, lasym1I };
			Float_t ccwgt[ncc] = { GetSqrMDR(3), 0.5*(GetSqrMDR(2)+GetSqrMDR(3)) };
			Float_t ccsum = 0.0, ccsgm = 0.0;
			for (Int_t it = 0; it < ncc; ++it) { ccsum += ccelm[it] * ccwgt[it]; ccsgm += ccwgt[it]; }
			ccest = ccsum / ccsgm;
		}
		
		// New
		if (trPt == 5 && (lchiyL1 > 3.0 || lchiyL9 > 3.0 || lchiyFs > 2.5 || trFsLresyL1 > 2.0 || trFsLresyL9 > 2.0)) break;
		if (trPt == 4 && (lchiyL9 > 2.5 || trL9LresyL9 > 2.0)) break;
		if (trPt == 3 && (lchiyL1 > 2.5 || trL1LresyL1 > 2.0)) break;
		
		Hist::Head("hHi_CCestQ")->fill(irig, ccest, weight);
		if (sign>0) Hist::Head("hHp_CCestQ")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hHn_CCestQ")->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_CCestQ%03d", tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_CCestQ%03d", tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_CCestQ%03d", tvBin))->fill(nrig, ccest, weight);
	
		Hist::Head(StrFmt("hHi_CCestQ%s", trNm.c_str()))->fill(irig, ccest, weight);
		if (sign>0) Hist::Head(StrFmt("hHp_CCestQ%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_CCestQ%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_CCestQ%s%03d", trNm.c_str(), tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_CCestQ%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_CCestQ%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		
		// New
		if (trPt == 5 && (lasym1I > 2.1 || lasym9I > 1.6 || lasym91 > 1.8)) break;
		if (trPt == 4 && (lasym9I > 1.5)) break;
		if (trPt == 3 && (lasym1I > 1.9)) break;
		
		// Org
		//if (trPt == 5 && (lchixL1 > 3.5 || lchixL9 > 3.5 || lchixFs > 3.5 || trFsLresxL1 > 2.0 || trFsLresxL9 > 2.0)) break;
		//if (trPt == 5 && (lchiyL1 > 3.5 || lchiyL9 > 2.5 || lchiyFs > 2.5 || trFsLresyL1 > 2.0 || trFsLresyL9 > 2.0)) break; // lchiyFs 2.5   new (2.2  1.9)
		//if (trPt == 5 && (lasym1I > 2.0 || lasym9I > 2.0 || lasym91 > 2.0)) break; // lasym91 2.0   new (1.8  1.5)
		//if (trPt == 4 && (lchixL9 > 3.5 || trL9LresxL9 > 2.0)) break;
		//if (trPt == 4 && (lchiyL9 > 1.8 || trL9LresyL9 > 2.0)) break; // 2.5
		//if (trPt == 4 && (lasym9I > 1.5)) break; // 2.0
		//if (trPt == 3 && (lchixL1 > 3.5 || trL1LresxL1 > 2.0)) break;
		//if (trPt == 3 && (lchiyL1 > 2.4 || trL1LresyL1 > 2.0)) break; // 3.5
		//if (trPt == 3 && (lasym1I > 1.6)) break; // 2.0
		
		Hist::Head("hHi_Cutflow")->fill(irig, 6, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 6, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 6, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 6, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
		
		Hist::Head("hHi_CCest")->fill(irig, ccest, weight);
		if (sign>0) Hist::Head("hHp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hHn_CCest")->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_CCest%03d", tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_CCest%03d", tvBin))->fill(nrig, ccest, weight);
	
		Hist::Head(StrFmt("hHi_CCest%s", trNm.c_str()))->fill(irig, ccest, weight);
		if (sign>0) Hist::Head(StrFmt("hHp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_CCest%s%03d", trNm.c_str(), tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_CCest%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_CCest%s%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		
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
			Hist::Head("hHi_MCTrRso")->fill(MCIRig, irig, weight);
			if (MCTRig>0.0) Hist::Head("hHp_MCTrRso")->fill(MCNRig, (irig - MCIRig), weight);
			if (MCTRig<0.0) Hist::Head("hHn_MCTrRso")->fill(MCNRig, (irig - MCIRig), weight);
			
			Hist::Head(StrFmt("hHi_MCTrRso%s", trNm.c_str()))->fill(MCIRig, irig, weight);
			if (MCTRig>0.0) Hist::Head(StrFmt("hHp_MCTrRso%s", trNm.c_str()))->fill(MCNRig, (irig - MCIRig), weight);
			if (MCTRig<0.0) Hist::Head(StrFmt("hHn_MCTrRso%s", trNm.c_str()))->fill(MCNRig, (irig - MCIRig), weight);
		
			Hist::Head("hHi_MCEvt")->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head("hHp_MCEvt")->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head("hHn_MCEvt")->fill(MCNRig, weight);
			
			Hist::Head(StrFmt("hHi_MCEvt%s", trNm.c_str()))->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head(StrFmt("hHp_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head(StrFmt("hHn_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
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
