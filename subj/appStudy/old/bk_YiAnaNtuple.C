#ifndef __YiAnaNtuple_C__
#define __YiAnaNtuple_C__

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.h"
#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.tcc"

#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2016Nov23/src/ClassDef.h"
#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2016Nov23/src/ClassDef.C"

using namespace MgntROOT;

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
			
			TString sbn = "/afs/cern.ch/user/h/hchou/public/DATABASE/physics/binning/antipp3_rig.root";
			TFile * fBin = TFile::Open(sbn);
			TH1D  * hBn0 = (TH1D*)fBin->Get("hbin0");
			TH1D  * hBn1 = (TH1D*)fBin->Get("hbin1");
			AXnr = Axis(hBn0, Axis::kX);
      AXir = Axis(hBn1, Axis::kX);
			static Axis AXrso("Rso", 400, -2.5, 2.5);
			fBin->Close();
			file->cd();
	
			// 6 month
			AXtme = Axis("Time", 
				{ 1312416000, 1326412800, 1340409600, 1354406400, 1368403200, 
				  1382400000, 1396396800, 1410393600 ,1424390400, 1438387200, 
					1452384000 , 1464299612 } );
			
			// 3 month
			//AXtme = Axis("Time", 
			//	{ 1326412800, 1333411200, 1340409600, 1347408000, 1354406400, 
			//	  1361404800, 1368403200, 1375401600, 1382400000, 1389398400, 
			//		1396396800, 1403395200, 1410393600, 1417392000, 1424390400, 
			//		1431388800, 1438387200, 1445385600, 1452384000, 1459382400 } );

			//TString sbnPH = "/afs/cern.ch/user/h/hchou/public/DATABASE/physics/binning/phe_bin2.root";
			//TFile * fPHBin = TFile::Open(sbnPH);
			//TH1D  * hPHBn0 = (TH1D*)fPHBin->Get("hist2");
			//Axis AXPHnr("Rigidity [GV]", hPHBn0, Axis::kX);
      //Axis AXPHir = MgntROOT::Axis::Invert("Inverse Rigidity [1/GV]", AXPHnr);
			//Axis AXnr = AXPHnr;
			//Axis AXir = AXPHir;
			//fPHBin->Close();
			//file->cd();

			// Cutoff and ExpTime
			//Hist::New("hGeoR25", "Geomagnetic Cutoff (25 degree)", AXir, AXnr);
			//Hist::New("hGeoR30", "Geomagnetic Cutoff (30 degree)", AXir, AXnr);
			//Hist::New("hGeoR35", "Geomagnetic Cutoff (35 degree)", AXir, AXnr);
			//Hist::New("hGeoR40", "Geomagnetic Cutoff (40 degree)", AXir, AXnr);
			//Hist::New("hExpT25", "Exposure Time (25 degree)", AXnr);
			//Hist::New("hExpT30", "Exposure Time (30 degree)", AXnr);
			//Hist::New("hExpT35", "Exposure Time (35 degree)", AXnr);
			//Hist::New("hExpT40", "Exposure Time (40 degree)", AXnr);


			//----  Low Energy  ----//
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Axis AXLcutoff("Cutoff", 400, 0., 35.);
				Hist::New("hLi_Cutoff", "", AXir, AXLcutoff);
				Hist::New("hLp_Cutoff", "", AXnr, AXLcutoff);
				Hist::New("hLn_Cutoff", "", AXnr, AXLcutoff);
				for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
					Hist::New(StrFmt("hLi_Cutoff%03d", it), "", AXir, AXLcutoff);
					Hist::New(StrFmt("hLp_Cutoff%03d", it), "", AXnr, AXLcutoff);
					Hist::New(StrFmt("hLn_Cutoff%03d", it), "", AXnr, AXLcutoff);
				}
			}
			
			Axis AXLcutflow("Cutflow", 8, 0., 8.);
			Hist::New("hLi_Cutflow", "", AXir, AXLcutflow);
			Hist::New("hLp_Cutflow", "", AXnr, AXLcutflow);
			Hist::New("hLn_Cutflow", "", AXnr, AXLcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Cutflow%03d", it), "", AXir, AXLcutflow);
				Hist::New(StrFmt("hLp_Cutflow%03d", it), "", AXnr, AXLcutflow);
				Hist::New(StrFmt("hLn_Cutflow%03d", it), "", AXnr, AXLcutflow);
			}
			
			Axis AXLtrdcls("Trdcls", 150, 0., 150.);
			Hist::New("hLi_Trdcls", "", AXir, AXLtrdcls);
			Hist::New("hLp_Trdcls", "", AXnr, AXLtrdcls);
			Hist::New("hLn_Trdcls", "", AXnr, AXLtrdcls);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Trdcls%03d", it), "", AXir, AXLtrdcls);
				Hist::New(StrFmt("hLp_Trdcls%03d", it), "", AXnr, AXLtrdcls);
				Hist::New(StrFmt("hLn_Trdcls%03d", it), "", AXnr, AXLtrdcls);
			}

			Axis AXLtrdvtx("Trdvtx", 15, 0., 15.);
			Hist::New("hLi_Trdvtx", "", AXir, AXLtrdvtx);
			Hist::New("hLp_Trdvtx", "", AXnr, AXLtrdvtx);
			Hist::New("hLn_Trdvtx", "", AXnr, AXLtrdvtx);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Trdvtx%03d", it), "", AXir, AXLtrdvtx);
				Hist::New(StrFmt("hLp_Trdvtx%03d", it), "", AXnr, AXLtrdvtx);
				Hist::New(StrFmt("hLn_Trdvtx%03d", it), "", AXnr, AXLtrdvtx);
			}
			
			Axis AXLtrdl("Trdl", 400, 0.2, 1.6);
			Hist::New("hLi_Trdl", "", AXir, AXLtrdl);
			Hist::New("hLp_Trdl", "", AXnr, AXLtrdl);
			Hist::New("hLn_Trdl", "", AXnr, AXLtrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Trdl%03d", it), "", AXir, AXLtrdl);
				Hist::New(StrFmt("hLp_Trdl%03d", it), "", AXnr, AXLtrdl);
				Hist::New(StrFmt("hLn_Trdl%03d", it), "", AXnr, AXLtrdl);
			}

			Axis AXLasym("Asym", 400, -8., 8.);
			Hist::New("hLi_AsymUL", "", AXir, AXLasym);
			Hist::New("hLp_AsymUL", "", AXnr, AXLasym);
			Hist::New("hLn_AsymUL", "", AXnr, AXLasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_AsymUL%03d", it), "", AXir, AXLasym);
				Hist::New(StrFmt("hLp_AsymUL%03d", it), "", AXnr, AXLasym);
				Hist::New(StrFmt("hLn_AsymUL%03d", it), "", AXnr, AXLasym);
			}
			
			Axis AXLlasym("Lasym", 400, -8., 8.);
			Hist::New("hLi_LasymUL", "", AXir, AXLlasym);
			Hist::New("hLp_LasymUL", "", AXnr, AXLlasym);
			Hist::New("hLn_LasymUL", "", AXnr, AXLlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_LasymUL%03d", it), "", AXir, AXLlasym);
				Hist::New(StrFmt("hLp_LasymUL%03d", it), "", AXnr, AXLlasym);
				Hist::New(StrFmt("hLn_LasymUL%03d", it), "", AXnr, AXLlasym);
			}
			
			Axis AXLlchi("Lchi", 400, -8., 8.);
			Hist::New("hLi_LchixIU", "", AXir, AXLlchi);
			Hist::New("hLp_LchixIU", "", AXnr, AXLlchi);
			Hist::New("hLn_LchixIU", "", AXnr, AXLlchi);
			Hist::New("hLi_LchiyIU", "", AXir, AXLlchi);
			Hist::New("hLp_LchiyIU", "", AXnr, AXLlchi);
			Hist::New("hLn_LchiyIU", "", AXnr, AXLlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_LchixIU%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_LchixIU%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_LchixIU%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLi_LchiyIU%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_LchiyIU%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_LchiyIU%03d", it), "", AXnr, AXLlchi);
			}
			
			Hist::New("hLi_LchixIL", "", AXir, AXLlchi);
			Hist::New("hLp_LchixIL", "", AXnr, AXLlchi);
			Hist::New("hLn_LchixIL", "", AXnr, AXLlchi);
			Hist::New("hLi_LchiyIL", "", AXir, AXLlchi);
			Hist::New("hLp_LchiyIL", "", AXnr, AXLlchi);
			Hist::New("hLn_LchiyIL", "", AXnr, AXLlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_LchixIL%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_LchixIL%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_LchixIL%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLi_LchiyIL%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_LchiyIL%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_LchiyIL%03d", it), "", AXnr, AXLlchi);
			}
			
			Hist::New("hLi_LchixIn", "", AXir, AXLlchi);
			Hist::New("hLp_LchixIn", "", AXnr, AXLlchi);
			Hist::New("hLn_LchixIn", "", AXnr, AXLlchi);
			Hist::New("hLi_LchiyIn", "", AXir, AXLlchi);
			Hist::New("hLp_LchiyIn", "", AXnr, AXLlchi);
			Hist::New("hLn_LchiyIn", "", AXnr, AXLlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_LchixIn%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_LchixIn%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_LchixIn%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLi_LchiyIn%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_LchiyIn%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_LchiyIn%03d", it), "", AXnr, AXLlchi);
			}

			Hist::New("hLi_Lchix", "", AXir, AXLlchi);
			Hist::New("hLp_Lchix", "", AXnr, AXLlchi);
			Hist::New("hLn_Lchix", "", AXnr, AXLlchi);
			Hist::New("hLi_Lchiy", "", AXir, AXLlchi);
			Hist::New("hLp_Lchiy", "", AXnr, AXLlchi);
			Hist::New("hLn_Lchiy", "", AXnr, AXLlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Lchix%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_Lchix%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_Lchix%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLi_Lchiy%03d", it), "", AXir, AXLlchi);
				Hist::New(StrFmt("hLp_Lchiy%03d", it), "", AXnr, AXLlchi);
				Hist::New(StrFmt("hLn_Lchiy%03d", it), "", AXnr, AXLlchi);
			}
			
			Axis AXLccest("CCest", 400, -8., 8.);
			Hist::New("hLi_CCest", "", AXir, AXLccest);
			Hist::New("hLp_CCest", "", AXnr, AXLccest);
			Hist::New("hLn_CCest", "", AXnr, AXLccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_CCest%03d", it), "", AXir, AXLccest);
				Hist::New(StrFmt("hLp_CCest%03d", it), "", AXnr, AXLccest);
				Hist::New(StrFmt("hLn_CCest%03d", it), "", AXnr, AXLccest);
			}
			
			Axis AXLaglprph("Aglprph", 400, 0., 10.);
			Hist::New("hLi_Aglprph", "", AXir, AXLaglprph);
			Hist::New("hLp_Aglprph", "", AXnr, AXLaglprph);
			Hist::New("hLn_Aglprph", "", AXnr, AXLaglprph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Aglprph%03d", it), "", AXir, AXLaglprph);
				Hist::New(StrFmt("hLp_Aglprph%03d", it), "", AXnr, AXLaglprph);
				Hist::New(StrFmt("hLn_Aglprph%03d", it), "", AXnr, AXLaglprph);
			}
			
			Axis AXLaglelph("Aglelph", 400, 0., 10.);
			Hist::New("hLi_Aglelph", "", AXir, AXLaglelph);
			Hist::New("hLp_Aglelph", "", AXnr, AXLaglelph);
			Hist::New("hLn_Aglelph", "", AXnr, AXLaglelph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Aglelph%03d", it), "", AXir, AXLaglelph);
				Hist::New(StrFmt("hLp_Aglelph%03d", it), "", AXnr, AXLaglelph);
				Hist::New(StrFmt("hLn_Aglelph%03d", it), "", AXnr, AXLaglelph);
			}
			
			Axis AXLaglph("Aglph", 50, 0., 50.);
			Hist::New("hLi_Aglph", "", AXir, AXLaglph);
			Hist::New("hLp_Aglph", "", AXnr, AXLaglph);
			Hist::New("hLn_Aglph", "", AXnr, AXLaglph);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Aglph%03d", it), "", AXir, AXLaglph);
				Hist::New(StrFmt("hLp_Aglph%03d", it), "", AXnr, AXLaglph);
				Hist::New(StrFmt("hLn_Aglph%03d", it), "", AXnr, AXLaglph);
			}
			
			Axis AXLaglb("Aglb", 400, -0.04, 0.02);
			Hist::New("hLi_Aglb", "", AXir, AXLaglb);
			Hist::New("hLp_Aglb", "", AXnr, AXLaglb);
			Hist::New("hLn_Aglb", "", AXnr, AXLaglb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Aglb%03d", it), "", AXir, AXLaglb);
				Hist::New(StrFmt("hLp_Aglb%03d", it), "", AXnr, AXLaglb);
				Hist::New(StrFmt("hLn_Aglb%03d", it), "", AXnr, AXLaglb);
			}
			
			Axis AXLtofaglb("TofAglb", 400, -0.16, 0.16);
			Hist::New("hLi_TofAglb", "", AXir, AXLtofaglb);
			Hist::New("hLp_TofAglb", "", AXnr, AXLtofaglb);
			Hist::New("hLn_TofAglb", "", AXnr, AXLtofaglb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_TofAglb%03d", it), "", AXir, AXLtofaglb);
				Hist::New(StrFmt("hLp_TofAglb%03d", it), "", AXnr, AXLtofaglb);
				Hist::New(StrFmt("hLn_TofAglb%03d", it), "", AXnr, AXLtofaglb);
			}
			
			Axis AXLtofb("Tofb", 50, -0.16, 0.4);
			Hist::New("hLi_Tofb", "", AXir, AXLtofb);
			Hist::New("hLp_Tofb", "", AXnr, AXLtofb);
			Hist::New("hLn_Tofb", "", AXnr, AXLtofb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Tofb%03d", it), "", AXir, AXLtofb);
				Hist::New(StrFmt("hLp_Tofb%03d", it), "", AXnr, AXLtofb);
				Hist::New(StrFmt("hLn_Tofb%03d", it), "", AXnr, AXLtofb);
			}
		
			Hist::New("hLs_Tofb", "", AXnr, AXLtofb);
			Hist::New("hLb_Tofb", "", AXnr, AXLtofb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLs_Tofb%03d", it), "", AXnr, AXLtofb);
				Hist::New(StrFmt("hLb_Tofb%03d", it), "", AXnr, AXLtofb);
			}

			Hist::New("hLi_Evt", "", AXir);
			Hist::New("hLp_Evt", "", AXnr);
			Hist::New("hLn_Evt", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hLi_Evt%03d", it), "", AXir);
				Hist::New(StrFmt("hLp_Evt%03d", it), "", AXnr);
				Hist::New(StrFmt("hLn_Evt%03d", it), "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hLi_MCTrRso", "", AXir, AXir);
				Hist::New("hLp_MCTrRso", "", AXnr, AXrso);
				Hist::New("hLn_MCTrRso", "", AXnr, AXrso);
				
				Hist::New("hLi_MCEvt", "", AXir);
				Hist::New("hLp_MCEvt", "", AXnr);
				Hist::New("hLn_MCEvt", "", AXnr);
			}
			
			//----  Intermedia Energy  ----//
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Axis AXIcutoff("Cutoff", 400, 0., 35.);
				Hist::New("hIi_Cutoff", "", AXir, AXIcutoff);
				Hist::New("hIp_Cutoff", "", AXnr, AXIcutoff);
				Hist::New("hIn_Cutoff", "", AXnr, AXIcutoff);
				for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
					Hist::New(StrFmt("hIi_Cutoff%03d", it), "", AXir, AXIcutoff);
					Hist::New(StrFmt("hIp_Cutoff%03d", it), "", AXnr, AXIcutoff);
					Hist::New(StrFmt("hIn_Cutoff%03d", it), "", AXnr, AXIcutoff);
				}
			}
			
			Axis AXIcutflow("Cutflow", 6, 0., 6.);
			Hist::New("hIi_Cutflow", "", AXir, AXIcutflow);
			Hist::New("hIp_Cutflow", "", AXnr, AXIcutflow);
			Hist::New("hIn_Cutflow", "", AXnr, AXIcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Cutflow%03d", it), "", AXir, AXIcutflow);
				Hist::New(StrFmt("hIp_Cutflow%03d", it), "", AXnr, AXIcutflow);
				Hist::New(StrFmt("hIn_Cutflow%03d", it), "", AXnr, AXIcutflow);
			}
			
			Axis AXIecalbdt("Ecalbdt", 400, -1.0, 1.0);
			Hist::New("hIi_Ecalbdt", "", AXir, AXIecalbdt);
			Hist::New("hIp_Ecalbdt", "", AXnr, AXIecalbdt);
			Hist::New("hIn_Ecalbdt", "", AXnr, AXIecalbdt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Ecalbdt%03d", it), "", AXir, AXIecalbdt);
				Hist::New(StrFmt("hIp_Ecalbdt%03d", it), "", AXnr, AXIecalbdt);
				Hist::New(StrFmt("hIn_Ecalbdt%03d", it), "", AXnr, AXIecalbdt);
			}

			Axis AXIaglb("Aglb", 400, -0.04, 0.02);
			Hist::New("hIi_Aglb", "", AXir, AXIaglb);
			Hist::New("hIp_Aglb", "", AXnr, AXIaglb);
			Hist::New("hIn_Aglb", "", AXnr, AXIaglb);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Aglb%03d", it), "", AXir, AXIaglb);
				Hist::New(StrFmt("hIp_Aglb%03d", it), "", AXnr, AXIaglb);
				Hist::New(StrFmt("hIn_Aglb%03d", it), "", AXnr, AXIaglb);
			}
			
			Axis AXIasym("Asym", 400, -8., 8.);
			Hist::New("hIi_AsymUL", "", AXir, AXIasym);
			Hist::New("hIp_AsymUL", "", AXnr, AXIasym);
			Hist::New("hIn_AsymUL", "", AXnr, AXIasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_AsymUL%03d", it), "", AXir, AXIasym);
				Hist::New(StrFmt("hIp_AsymUL%03d", it), "", AXnr, AXIasym);
				Hist::New(StrFmt("hIn_AsymUL%03d", it), "", AXnr, AXIasym);
			}
			
			Hist::New("hIi_Asym1I", "", AXir, AXIasym);
			Hist::New("hIp_Asym1I", "", AXnr, AXIasym);
			Hist::New("hIn_Asym1I", "", AXnr, AXIasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Asym1I%03d", it), "", AXir, AXIasym);
				Hist::New(StrFmt("hIp_Asym1I%03d", it), "", AXnr, AXIasym);
				Hist::New(StrFmt("hIn_Asym1I%03d", it), "", AXnr, AXIasym);
			}
			
			Hist::New("hIi_Asym9I", "", AXir, AXIasym);
			Hist::New("hIp_Asym9I", "", AXnr, AXIasym);
			Hist::New("hIn_Asym9I", "", AXnr, AXIasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Asym9I%03d", it), "", AXir, AXIasym);
				Hist::New(StrFmt("hIp_Asym9I%03d", it), "", AXnr, AXIasym);
				Hist::New(StrFmt("hIn_Asym9I%03d", it), "", AXnr, AXIasym);
			}
			
			Hist::New("hIi_Asym91", "", AXir, AXIasym);
			Hist::New("hIp_Asym91", "", AXnr, AXIasym);
			Hist::New("hIn_Asym91", "", AXnr, AXIasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Asym91%03d", it), "", AXir, AXIasym);
				Hist::New(StrFmt("hIp_Asym91%03d", it), "", AXnr, AXIasym);
				Hist::New(StrFmt("hIn_Asym91%03d", it), "", AXnr, AXIasym);
			}
			
			Axis AXIlasym("Lasym", 400, -8., 8.);
			Hist::New("hIi_LasymUL", "", AXir, AXIlasym);
			Hist::New("hIp_LasymUL", "", AXnr, AXIlasym);
			Hist::New("hIn_LasymUL", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_LasymUL%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_LasymUL%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_LasymUL%03d", it), "", AXnr, AXIlasym);
			}
			
			Hist::New("hIi_Lasym1I", "", AXir, AXIlasym);
			Hist::New("hIp_Lasym1I", "", AXnr, AXIlasym);
			Hist::New("hIn_Lasym1I", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Lasym1I%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_Lasym1I%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_Lasym1I%03d", it), "", AXnr, AXIlasym);
			}
			
			Hist::New("hIi_Lasym9I", "", AXir, AXIlasym);
			Hist::New("hIp_Lasym9I", "", AXnr, AXIlasym);
			Hist::New("hIn_Lasym9I", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Lasym9I%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_Lasym9I%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_Lasym9I%03d", it), "", AXnr, AXIlasym);
			}
			
			Hist::New("hIi_Lasym91", "", AXir, AXIlasym);
			Hist::New("hIp_Lasym91", "", AXnr, AXIlasym);
			Hist::New("hIn_Lasym91", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Lasym91%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_Lasym91%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_Lasym91%03d", it), "", AXnr, AXIlasym);
			}
			
			Hist::New("hIi_Lasym", "", AXir, AXIlasym);
			Hist::New("hIp_Lasym", "", AXnr, AXIlasym);
			Hist::New("hIn_Lasym", "", AXnr, AXIlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Lasym%03d", it), "", AXir, AXIlasym);
				Hist::New(StrFmt("hIp_Lasym%03d", it), "", AXnr, AXIlasym);
				Hist::New(StrFmt("hIn_Lasym%03d", it), "", AXnr, AXIlasym);
			}
			
			Axis AXIlchi("Lchi", 400, -8., 8.);
			Hist::New("hIi_LchixIU", "", AXir, AXIlchi);
			Hist::New("hIp_LchixIU", "", AXnr, AXIlchi);
			Hist::New("hIn_LchixIU", "", AXnr, AXIlchi);
			Hist::New("hIi_LchiyIU", "", AXir, AXIlchi);
			Hist::New("hIp_LchiyIU", "", AXnr, AXIlchi);
			Hist::New("hIn_LchiyIU", "", AXnr, AXIlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_LchixIU%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchixIU%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchixIU%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIi_LchiyIU%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchiyIU%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchiyIU%03d", it), "", AXnr, AXIlchi);
			}
			
			Hist::New("hIi_LchixIL", "", AXir, AXIlchi);
			Hist::New("hIp_LchixIL", "", AXnr, AXIlchi);
			Hist::New("hIn_LchixIL", "", AXnr, AXIlchi);
			Hist::New("hIi_LchiyIL", "", AXir, AXIlchi);
			Hist::New("hIp_LchiyIL", "", AXnr, AXIlchi);
			Hist::New("hIn_LchiyIL", "", AXnr, AXIlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_LchixIL%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchixIL%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchixIL%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIi_LchiyIL%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchiyIL%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchiyIL%03d", it), "", AXnr, AXIlchi);
			}
			
			Hist::New("hIi_LchixIn", "", AXir, AXIlchi);
			Hist::New("hIp_LchixIn", "", AXnr, AXIlchi);
			Hist::New("hIn_LchixIn", "", AXnr, AXIlchi);
			Hist::New("hIi_LchiyIn", "", AXir, AXIlchi);
			Hist::New("hIp_LchiyIn", "", AXnr, AXIlchi);
			Hist::New("hIn_LchiyIn", "", AXnr, AXIlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
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
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
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
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
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
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_LchixFs%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchixFs%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchixFs%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIi_LchiyFs%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_LchiyFs%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_LchiyFs%03d", it), "", AXnr, AXIlchi);
			}
			
			Hist::New("hIi_Lchix", "", AXir, AXIlchi);
			Hist::New("hIp_Lchix", "", AXnr, AXIlchi);
			Hist::New("hIn_Lchix", "", AXnr, AXIlchi);
			Hist::New("hIi_Lchiy", "", AXir, AXIlchi);
			Hist::New("hIp_Lchiy", "", AXnr, AXIlchi);
			Hist::New("hIn_Lchiy", "", AXnr, AXIlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Lchix%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_Lchix%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_Lchix%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIi_Lchiy%03d", it), "", AXir, AXIlchi);
				Hist::New(StrFmt("hIp_Lchiy%03d", it), "", AXnr, AXIlchi);
				Hist::New(StrFmt("hIn_Lchiy%03d", it), "", AXnr, AXIlchi);
			}
			
			Axis AXIccest("CCest", 400, -8., 8.);
			Hist::New("hIi_CCest", "", AXir, AXIccest);
			Hist::New("hIp_CCest", "", AXnr, AXIccest);
			Hist::New("hIn_CCest", "", AXnr, AXIccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_CCest%03d", it), "", AXir, AXIccest);
				Hist::New(StrFmt("hIp_CCest%03d", it), "", AXnr, AXIccest);
				Hist::New(StrFmt("hIn_CCest%03d", it), "", AXnr, AXIccest);
			}

			Axis AXItrdl("Trdl", 50, 0.2, 1.6);
			Hist::New("hIi_Trdl", "", AXir, AXItrdl);
			Hist::New("hIp_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIn_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Trdl%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIs_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIb_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIs_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIb_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_In_Trdl", "", AXir, AXItrdl);
			Hist::New("hIp_In_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIn_In_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_In_Trdl%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_In_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_In_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_L1_Trdl", "", AXir, AXItrdl);
			Hist::New("hIp_L1_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIn_L1_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_L1_Trdl%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_L1_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_L1_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_L9_Trdl", "", AXir, AXItrdl);
			Hist::New("hIp_L9_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIn_L9_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_L9_Trdl%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_L9_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_L9_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_Fs_Trdl", "", AXir, AXItrdl);
			Hist::New("hIp_Fs_Trdl", "", AXnr, AXItrdl);
			Hist::New("hIn_Fs_Trdl", "", AXnr, AXItrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Fs_Trdl%03d", it), "", AXir, AXItrdl);
				Hist::New(StrFmt("hIp_Fs_Trdl%03d", it), "", AXnr, AXItrdl);
				Hist::New(StrFmt("hIn_Fs_Trdl%03d", it), "", AXnr, AXItrdl);
			}
			
			Hist::New("hIi_Evt", "", AXir);
			Hist::New("hIp_Evt", "", AXnr);
			Hist::New("hIn_Evt", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hIi_Evt%03d", it), "", AXir);
				Hist::New(StrFmt("hIp_Evt%03d", it), "", AXnr);
				Hist::New(StrFmt("hIn_Evt%03d", it), "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hIi_MCTrRso", "", AXir, AXir);
				Hist::New("hIp_MCTrRso", "", AXnr, AXrso);
				Hist::New("hIn_MCTrRso", "", AXnr, AXrso);
				
				Hist::New("hIi_MCEvt", "", AXir);
				Hist::New("hIp_MCEvt", "", AXnr);
				Hist::New("hIn_MCEvt", "", AXnr);
			}
			
			//----  High Energy  ----//
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Axis AXHcutoff("Cutoff", 400, 0., 35.);
				Hist::New("hHi_Cutoff", "", AXir, AXHcutoff);
				Hist::New("hHp_Cutoff", "", AXnr, AXHcutoff);
				Hist::New("hHn_Cutoff", "", AXnr, AXHcutoff);
				for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
					Hist::New(StrFmt("hHi_Cutoff%03d", it), "", AXir, AXHcutoff);
					Hist::New(StrFmt("hHp_Cutoff%03d", it), "", AXnr, AXHcutoff);
					Hist::New(StrFmt("hHn_Cutoff%03d", it), "", AXnr, AXHcutoff);
				}
			}
			
			Axis AXHcutflow("Cutflow", 5, 0., 5.);
			Hist::New("hHi_Cutflow", "", AXir, AXHcutflow);
			Hist::New("hHp_Cutflow", "", AXnr, AXHcutflow);
			Hist::New("hHn_Cutflow", "", AXnr, AXHcutflow);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Cutflow%03d", it), "", AXir, AXHcutflow);
				Hist::New(StrFmt("hHp_Cutflow%03d", it), "", AXnr, AXHcutflow);
				Hist::New(StrFmt("hHn_Cutflow%03d", it), "", AXnr, AXHcutflow);
			}
			
			Axis AXHtrdl("Trdl", 400, 0.2, 1.6);
			Hist::New("hHi_Trdl", "", AXir, AXHtrdl);
			Hist::New("hHp_Trdl", "", AXnr, AXHtrdl);
			Hist::New("hHn_Trdl", "", AXnr, AXHtrdl);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Trdl%03d", it), "", AXir, AXHtrdl);
				Hist::New(StrFmt("hHp_Trdl%03d", it), "", AXnr, AXHtrdl);
				Hist::New(StrFmt("hHn_Trdl%03d", it), "", AXnr, AXHtrdl);
			}
			
			Axis AXHecalbdt("Ecalbdt", 400, -1.0, 1.0);
			Hist::New("hHi_Ecalbdt", "", AXir, AXHecalbdt);
			Hist::New("hHp_Ecalbdt", "", AXnr, AXHecalbdt);
			Hist::New("hHn_Ecalbdt", "", AXnr, AXHecalbdt);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Ecalbdt%03d", it), "", AXir, AXHecalbdt);
				Hist::New(StrFmt("hHp_Ecalbdt%03d", it), "", AXnr, AXHecalbdt);
				Hist::New(StrFmt("hHn_Ecalbdt%03d", it), "", AXnr, AXHecalbdt);
			}

			Axis AXHasym("Asym", 400, -6., 6.);
			Hist::New("hHi_Asym1I", "", AXir, AXHasym);
			Hist::New("hHp_Asym1I", "", AXnr, AXHasym);
			Hist::New("hHn_Asym1I", "", AXnr, AXHasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Asym1I%03d", it), "", AXir, AXHasym);
				Hist::New(StrFmt("hHp_Asym1I%03d", it), "", AXnr, AXHasym);
				Hist::New(StrFmt("hHn_Asym1I%03d", it), "", AXnr, AXHasym);
			}
			
			Hist::New("hHi_Asym9I", "", AXir, AXHasym);
			Hist::New("hHp_Asym9I", "", AXnr, AXHasym);
			Hist::New("hHn_Asym9I", "", AXnr, AXHasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Asym9I%03d", it), "", AXir, AXHasym);
				Hist::New(StrFmt("hHp_Asym9I%03d", it), "", AXnr, AXHasym);
				Hist::New(StrFmt("hHn_Asym9I%03d", it), "", AXnr, AXHasym);
			}
			
			Hist::New("hHi_Asym91", "", AXir, AXHasym);
			Hist::New("hHp_Asym91", "", AXnr, AXHasym);
			Hist::New("hHn_Asym91", "", AXnr, AXHasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Asym91%03d", it), "", AXir, AXHasym);
				Hist::New(StrFmt("hHp_Asym91%03d", it), "", AXnr, AXHasym);
				Hist::New(StrFmt("hHn_Asym91%03d", it), "", AXnr, AXHasym);
			}
			
			Axis AXHlasym("Lasym", 400, -6., 6.);
			Hist::New("hHi_Lasym1I", "", AXir, AXHlasym);
			Hist::New("hHp_Lasym1I", "", AXnr, AXHlasym);
			Hist::New("hHn_Lasym1I", "", AXnr, AXHlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Lasym1I%03d", it), "", AXir, AXHlasym);
				Hist::New(StrFmt("hHp_Lasym1I%03d", it), "", AXnr, AXHlasym);
				Hist::New(StrFmt("hHn_Lasym1I%03d", it), "", AXnr, AXHlasym);
			}
			
			Hist::New("hHi_Lasym9I", "", AXir, AXHlasym);
			Hist::New("hHp_Lasym9I", "", AXnr, AXHlasym);
			Hist::New("hHn_Lasym9I", "", AXnr, AXHlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Lasym9I%03d", it), "", AXir, AXHlasym);
				Hist::New(StrFmt("hHp_Lasym9I%03d", it), "", AXnr, AXHlasym);
				Hist::New(StrFmt("hHn_Lasym9I%03d", it), "", AXnr, AXHlasym);
			}
			
			Hist::New("hHi_Lasym91", "", AXir, AXHlasym);
			Hist::New("hHp_Lasym91", "", AXnr, AXHlasym);
			Hist::New("hHn_Lasym91", "", AXnr, AXHlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Lasym91%03d", it), "", AXir, AXHlasym);
				Hist::New(StrFmt("hHp_Lasym91%03d", it), "", AXnr, AXHlasym);
				Hist::New(StrFmt("hHn_Lasym91%03d", it), "", AXnr, AXHlasym);
			}
			
			Hist::New("hHi_Lasym", "", AXir, AXHlasym);
			Hist::New("hHp_Lasym", "", AXnr, AXHlasym);
			Hist::New("hHn_Lasym", "", AXnr, AXHlasym);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Lasym%03d", it), "", AXir, AXHlasym);
				Hist::New(StrFmt("hHp_Lasym%03d", it), "", AXnr, AXHlasym);
				Hist::New(StrFmt("hHn_Lasym%03d", it), "", AXnr, AXHlasym);
			}
			
			Axis AXHlchi("Lchi", 400, -6., 6.);
			Hist::New("hHi_LchixL1", "", AXir, AXHlchi);
			Hist::New("hHp_LchixL1", "", AXnr, AXHlchi);
			Hist::New("hHn_LchixL1", "", AXnr, AXHlchi);
			Hist::New("hHi_LchiyL1", "", AXir, AXHlchi);
			Hist::New("hHp_LchiyL1", "", AXnr, AXHlchi);
			Hist::New("hHn_LchiyL1", "", AXnr, AXHlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
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
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
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
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_LchixFs%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_LchixFs%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_LchixFs%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHi_LchiyFs%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_LchiyFs%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_LchiyFs%03d", it), "", AXnr, AXHlchi);
			}
			
			Hist::New("hHi_Lchix", "", AXir, AXHlchi);
			Hist::New("hHp_Lchix", "", AXnr, AXHlchi);
			Hist::New("hHn_Lchix", "", AXnr, AXHlchi);
			Hist::New("hHi_Lchiy", "", AXir, AXHlchi);
			Hist::New("hHp_Lchiy", "", AXnr, AXHlchi);
			Hist::New("hHn_Lchiy", "", AXnr, AXHlchi);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Lchix%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_Lchix%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_Lchix%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHi_Lchiy%03d", it), "", AXir, AXHlchi);
				Hist::New(StrFmt("hHp_Lchiy%03d", it), "", AXnr, AXHlchi);
				Hist::New(StrFmt("hHn_Lchiy%03d", it), "", AXnr, AXHlchi);
			}
			
			Axis AXHccest("CCest", 50, -4., 4.);
			Hist::New("hHi_CCest", "", AXir, AXHccest);
			Hist::New("hHp_CCest", "", AXnr, AXHccest);
			Hist::New("hHn_CCest", "", AXnr, AXHccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_CCest%03d", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_CCest%03d", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_CCest%03d", it), "", AXnr, AXHccest);
			}
			
			Hist::New("hHi_L1_CCest", "", AXir, AXHccest);
			Hist::New("hHp_L1_CCest", "", AXnr, AXHccest);
			Hist::New("hHn_L1_CCest", "", AXnr, AXHccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_L1_CCest%03d", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_L1_CCest%03d", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_L1_CCest%03d", it), "", AXnr, AXHccest);
			}
			
			Hist::New("hHi_L9_CCest", "", AXir, AXHccest);
			Hist::New("hHp_L9_CCest", "", AXnr, AXHccest);
			Hist::New("hHn_L9_CCest", "", AXnr, AXHccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_L9_CCest%03d", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_L9_CCest%03d", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_L9_CCest%03d", it), "", AXnr, AXHccest);
			}
			
			Hist::New("hHi_Fs_CCest", "", AXir, AXHccest);
			Hist::New("hHp_Fs_CCest", "", AXnr, AXHccest);
			Hist::New("hHn_Fs_CCest", "", AXnr, AXHccest);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Fs_CCest%03d", it), "", AXir, AXHccest);
				Hist::New(StrFmt("hHp_Fs_CCest%03d", it), "", AXnr, AXHccest);
				Hist::New(StrFmt("hHn_Fs_CCest%03d", it), "", AXnr, AXHccest);
			}
			
			Hist::New("hHi_Evt", "", AXir);
			Hist::New("hHp_Evt", "", AXnr);
			Hist::New("hHn_Evt", "", AXnr);
			for (Int_t it = 1; (it <= AXtme.nbin()) && YiNtuple::CheckEventMode(YiNtuple::ISS); ++it) {
				Hist::New(StrFmt("hHi_Evt%03d", it), "", AXir);
				Hist::New(StrFmt("hHp_Evt%03d", it), "", AXnr);
				Hist::New(StrFmt("hHn_Evt%03d", it), "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hHi_MCTrRso", "", AXir, AXir);
				Hist::New("hHp_MCTrRso", "", AXnr, AXrso);
				Hist::New("hHn_MCTrRso", "", AXnr, AXrso);
				
				Hist::New("hHi_MCEvt", "", AXir);
				Hist::New("hHp_MCEvt", "", AXnr);
				Hist::New("hHn_MCEvt", "", AXnr);
			}
			
			std::cout << "\n<<  init DST End  >>\n";
			file->cd();
			return;
		}

		void resetDST() {
		}

		void fillDST() {
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
		static Axis AXtme;

		static std::map<std::string, Axis> gAX_Map;
		static std::map<std::string, std::vector<TH1 *> > gCut_H1DMap;
		static std::map<std::string, std::vector<TH2 *> > gCut_H2DMap;
		static std::map<std::string, std::vector<TH3 *> > gCut_H3DMap;
};

Bool_t DST::gSaveDST = true;
Bool_t DST::gSaveDSTTree = false;
Bool_t DST::gSaveDSTClone = false;
UInt_t DST::gUTime[2] = {0, 0};

Axis DST::AXnr;
Axis DST::AXir;
Axis DST::AXtme;
std::map<std::string, Axis> DST::gAX_Map;
std::map<std::string, std::vector<TH1 *> > DST::gCut_H1DMap;
std::map<std::string, std::vector<TH2 *> > DST::gCut_H2DMap;
std::map<std::string, std::vector<TH3 *> > DST::gCut_H3DMap;


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
	if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
		//----  RTI Cut  ----//
		// TODO: To be fixed. NOTE:Is Not Work

		if (gUTime[0] != fRti->uTime) {
			gUTime[0] = (gUTime[1] == 0) ? gUTime[1] : fRti->uTime;
			gUTime[1] = fRti->uTime;

			/*
			// Expourse Time
			const Float_t stableFact = 1.2;
			const Int_t degree[4] = { 25, 30, 35, 40 };
			for (Int_t iAgl = 0; iAgl < 4; ++iAgl) {
				Hist * hist = Hist::Head(StrFmt("hExpT%d", degree[iAgl]));
				Float_t cutoff = stableFact * fRti->cutoffStormer[iAgl];
				Int_t initBin = (hist->axisX())->find(cutoff) + 1;
				for (Int_t ibin = initBin; ibin <= (hist->axisX())->nbin(); ++ibin) {
					Float_t cenVal = (hist->axisX())->center(ibin);
					hist->fill(cenVal, fRti->liveTime);
				}
			}
			*/
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
	Double_t runID   = fList->runID;
	Double_t eventID = fList->eventID;
	Double_t entryID = fList->entryID;
	Double_t weight  = fList->weight;

	// Constant
	const Float_t mass = 0.938272297; 
	const Float_t massPi = 0.139570180; 
	const Float_t massEl = 0.000510999; 

	// Monte Carlo
	Float_t MCTRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
	                 0.0 : (fG4mc->primPart.mom / fG4mc->primPart.chrg);
	Float_t MCNRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ? 
	                 0.0 : std::fabs(MCTRig);
	Float_t MCIRig = (!YiNtuple::CheckEventMode(YiNtuple::MC)) ?
	                 0.0 : (1. / MCTRig);

	// Time
	Bool_t tvStudy = YiNtuple::CheckEventMode(YiNtuple::ISS) && (fRti->uTime > AXtme.min() && fRti->uTime < AXtme.max());
	Int_t tvBin = (!tvStudy) ? 0 : AXtme.find(fRti->uTime);

	// preselection (based on RTI)
	const Float_t CfStableFT = 1.2;
	if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
		if (!fRti->flagRun) return false;
		//if (!fRti->isGoodSecond) return false;
		if (fRti->isInSAA) return false;
		if (fRti->isInShadow) return false;
	}

	// preselection (based on Trigger)
	if ((fTrg->bit&8) != 8) return false;

	// preselection (Mini-Track Requirement)
	TrackInfo * track = (fTrk->tracks.size() == 1) ? (&fTrk->tracks.at(0)) : nullptr;
  if (track == nullptr) return false;
	//if ((track->bitPatt&2) != 2) return false;

	// preselection (Mini-Shower Requirement)
	ShowerInfo * shower = (fEcal->showers.size() >= 1) ? (&fEcal->showers.at(0)) : nullptr;

  // preselection (based on BetaH)
	if (!fTof->statusBetaH) return false;
	if (fTof->betaHPatt != 15) return false;
  if (fTof->normChisqT > 10) return false;
  if (fTof->normChisqC > 10) return false;
	if (fTof->betaH < 0.5 || fTof->betaH > 1.3) return false;
	//if (fTof->numOfInTimeCluster > 4) return false;
	
  // preselection (based on Charge)
	if (track->Qinner > 1.4 || track->Qinner < 0.7) return false;
	Float_t tofQu = 0.5 * (fTof->Q[0] + fTof->Q[1]);
	Float_t tofQl = 0.5 * (fTof->Q[2] + fTof->Q[3]);
	if (tofQu < 0.7 || tofQu > 1.6) return false;
	if (tofQl < 0.7 || tofQl > 1.6) return false;
	
	/**** Low Energy Region ****/
	while (true) {
		// TRK
		Int_t trPt = 2;
		std::string trNm = "In";
		if (!track->status[0][trPt]) break;
		if (!(track->status[0][0] && track->status[0][1])) break;
		Float_t trig       = track->rigidity[0][trPt];
		Int_t   sign       = MgntNum::Compare(trig);
		Float_t nrig       = std::fabs(trig);
		Float_t irig       = 1. / trig;
		Float_t tbta       = 1. / std::sqrt(1.+(mass/nrig)*(mass/nrig));
		Float_t tbtaPi     = 1. / std::sqrt(1.+(massPi/nrig)*(massPi/nrig));
		Float_t tbtaEl     = 1. / std::sqrt(1.+(massEl/nrig)*(massEl/nrig));
		
		// Cutoff
		if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
			Int_t icf = 0;
			Hist::Head("hLi_Cutoff")->fill(irig, fRti->cutoffIGRF[icf], weight);
			if (sign>0) Hist::Head("hLp_Cutoff")->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (sign<0) Hist::Head("hLn_Cutoff")->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutoff%03d", tvBin))->fill(irig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutoff%03d", tvBin))->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutoff%03d", tvBin))->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (nrig < CfStableFT * fRti->cutoffIGRF[icf]) break;
		}
		
		Hist::Head("hLi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		
		// TRD
		Int_t ntrdcls = fTrd->numOfCluster;
		Int_t ntrdvtx = fTrd->numOfVertexWithTrTrack;
		
		Hist::Head("hLi_Trdcls")->fill(irig, ntrdcls, weight);
		if (sign>0) Hist::Head("hLp_Trdcls")->fill(nrig, ntrdcls, weight);
		if (sign<0) Hist::Head("hLn_Trdcls")->fill(nrig, ntrdcls, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Trdcls%03d", tvBin))->fill(irig, ntrdcls, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Trdcls%03d", tvBin))->fill(nrig, ntrdcls, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Trdcls%03d", tvBin))->fill(nrig, ntrdcls, weight);
		
		Hist::Head("hLi_Trdvtx")->fill(irig, ntrdvtx, weight);
		if (sign>0) Hist::Head("hLp_Trdvtx")->fill(nrig, ntrdvtx, weight);
		if (sign<0) Hist::Head("hLn_Trdvtx")->fill(nrig, ntrdvtx, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Trdvtx%03d", tvBin))->fill(irig, ntrdvtx, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Trdvtx%03d", tvBin))->fill(nrig, ntrdvtx, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Trdvtx%03d", tvBin))->fill(nrig, ntrdvtx, weight);
	
		if (ntrdvtx >= 5) break;
		if (ntrdcls >= 60) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		Int_t trdPt = 0;
		Bool_t  hastrd  = (fTrd->statusKCls[trdPt] && fTrd->LLR_nhit[trdPt] >= 8);
		if (!hastrd) break;
		
		Float_t trdl    = (fTrd->LLR[trdPt][0]);
		Bool_t  istrdHe = (fTrd->LLR[trdPt][2] > 0.3);
		Bool_t  istrdPr = (trdl > 0.8);
		if (istrdHe) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		
		Hist::Head("hLi_Trdl")->fill(irig, trdl, weight);
		if (sign>0) Hist::Head("hLp_Trdl")->fill(nrig, trdl, weight);
		if (sign<0) Hist::Head("hLn_Trdl")->fill(nrig, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		
		if (!istrdPr) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);

		// Charge Confusion ---- TrSigma
		Float_t trParIU[3] = { 1.44615e-02, 1.58178e-03, 5.11515e-05 };
		Float_t trParIL[3] = { 9.60391e-03, 7.54520e-04, 3.81112e-05 };
		Float_t trParIn[3] = { 8.39123e-03, 3.36288e-04, 1.20139e-05 };
		Float_t trSgmIU = ((!track->status[0][0]) ? 1.0 : 
		                  std::sqrt(trParIU[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIU[1] * nrig + trParIU[2] * nrig * nrig)) / nrig;
		Float_t trSgmIL = ((!track->status[0][1]) ? 1.0 : 
		                  std::sqrt(trParIL[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIL[1] * nrig + trParIL[2] * nrig * nrig)) / nrig;
		Float_t trSgmIn = ((!track->status[0][2]) ? 1.0 : 
		                  std::sqrt(trParIn[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIn[1] * nrig + trParIn[2] * nrig * nrig)) / nrig;
		
		Float_t trSgmPt = trSgmIn;
		Float_t MDRFT = std::erfc(5. * nrig * std::sqrt(trParIn[2]));
	
		Float_t trParIUSgm[4] = { trParIU[0] * mass * mass, trParIU[0], trParIU[1], trParIU[2] };
		Float_t trParILSgm[4] = { trParIL[0] * mass * mass, trParIL[0], trParIL[1], trParIL[2] };
		Float_t trParInSgm[4] = { trParIn[0] * mass * mass, trParIn[0], trParIn[1], trParIn[2] };
		Float_t trSgmIUSgm = (!track->status[0][0]) ? 0.0 : 
		                     (0.5 / trSgmIn) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParIUSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParIUSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParIUSgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmILSgm = (!track->status[0][1]) ? 0.0 : 
		                     (0.5 / trSgmIn) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParILSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParILSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParILSgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmInSgm = (!track->status[0][2]) ? 0.0 : 
		                     (0.5 / trSgmIn) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParInSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParInSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParInSgm[2] * trSgmPt) )
		                     ); 

		if (track->status[0][0]) trSgmIU = std::sqrt(trSgmIU * trSgmIU + trSgmIUSgm * trSgmIUSgm);
		if (track->status[0][1]) trSgmIL = std::sqrt(trSgmIL * trSgmIL + trSgmILSgm * trSgmILSgm);
		if (track->status[0][2]) trSgmIn = std::sqrt(trSgmIn * trSgmIn + trSgmInSgm * trSgmInSgm);

		// Charge Confusion ---- Asym
		const Float_t lasymLMT = 1.0e-8;
		const Float_t lasymSGM = 1.8;
	
		Bool_t  hasasymUL   = (track->status[0][0] && track->status[0][1]);
		Float_t crrParUL[6] = { 7.63668e-01, 2.13220e-01, -1.40775e+00, 9.60481e-01, 3.50012e+00, 1.49281e+00 };
		Float_t crrSgmUL    = crrParUL[0]*(std::erf(crrParUL[1]*(std::log(nrig)+crrParUL[2]))+1.0) + TMath::Landau(crrParUL[3], crrParUL[4], crrParUL[5]);
		Float_t trSgmUL     = (!hasasymUL) ? 0.0 : std::sqrt(trSgmIU * trSgmIU + trSgmIL * trSgmIL) * crrSgmUL;
		Float_t asymUL      = (!hasasymUL) ? 0.0 : (1.0/track->rigidity[0][0] - 1.0/track->rigidity[0][1]) / trSgmUL;
		Float_t lasymUL     = (!hasasymUL) ? 0.0 : std::log(asymUL * asymUL + lasymLMT) / lasymSGM;
	
		Float_t lasym = lasymUL;

		if (hasasymUL) Hist::Head("hLi_AsymUL")->fill(irig, asymUL, weight);
		if (hasasymUL && sign>0) Hist::Head("hLp_AsymUL")->fill(nrig, asymUL, weight);
		if (hasasymUL && sign<0) Hist::Head("hLn_AsymUL")->fill(nrig, asymUL, weight);
		if (tvBin!=0 && hasasymUL) Hist::Head(StrFmt("hLi_AsymUL%03d", tvBin))->fill(irig, asymUL, weight);
		if (tvBin!=0 && hasasymUL && sign>0) Hist::Head(StrFmt("hLp_AsymUL%03d", tvBin))->fill(nrig, asymUL, weight);
		if (tvBin!=0 && hasasymUL && sign<0) Hist::Head(StrFmt("hLn_AsymUL%03d", tvBin))->fill(nrig, asymUL, weight);
		
		if (hasasymUL) Hist::Head("hLi_LasymUL")->fill(irig, lasymUL, weight);
		if (hasasymUL && sign>0) Hist::Head("hLp_LasymUL")->fill(nrig, lasymUL, weight);
		if (hasasymUL && sign<0) Hist::Head("hLn_LasymUL")->fill(nrig, lasymUL, weight);
		if (tvBin!=0 && hasasymUL) Hist::Head(StrFmt("hLi_LasymUL%03d", tvBin))->fill(irig, lasymUL, weight);
		if (tvBin!=0 && hasasymUL && sign>0) Hist::Head(StrFmt("hLp_LasymUL%03d", tvBin))->fill(nrig, lasymUL, weight);
		if (tvBin!=0 && hasasymUL && sign<0) Hist::Head(StrFmt("hLn_LasymUL%03d", tvBin))->fill(nrig, lasymUL, weight);
		
		// Charge Confusion ---- ChisqCC
		Float_t lchixIU = (!track->status[0][0]) ? 0.0 : std::log(track->chisq[0][0][0]);
		Float_t lchiyIU = (!track->status[0][0]) ? 0.0 : std::log(track->chisq[0][0][1]);
		Float_t lchixIL = (!track->status[0][1]) ? 0.0 : std::log(track->chisq[0][1][0]);
		Float_t lchiyIL = (!track->status[0][1]) ? 0.0 : std::log(track->chisq[0][1][1]);
		Float_t lchixIn = (!track->status[0][2]) ? 0.0 : std::log(track->chisq[0][2][0]);
		Float_t lchiyIn = (!track->status[0][2]) ? 0.0 : std::log(track->chisq[0][2][1]);
		
		std::vector<Float_t> chixvec;
		if (track->status[0][0]) { chixvec.push_back(lchixIU); }
		if (track->status[0][1]) { chixvec.push_back(lchixIL); }
		if (track->status[0][2]) { chixvec.push_back(lchixIn); }
		Float_t lchix = (chixvec.size()==0) ? 0.0 : *std::max_element(std::begin(chixvec), std::end(chixvec));
		
		std::vector<Float_t> chiyvec;
		if (track->status[0][0]) { chiyvec.push_back(lchiyIU); }
		if (track->status[0][1]) { chiyvec.push_back(lchiyIL); }
		if (track->status[0][2]) { chiyvec.push_back(lchiyIn); }
		Float_t lchiy = (chiyvec.size()==0) ? 0.0 : *std::max_element(std::begin(chiyvec), std::end(chiyvec));

		if (track->status[0][0]) Hist::Head("hLi_LchixIU")->fill(irig, lchixIU, weight);
		if (track->status[0][0] && sign>0) Hist::Head("hLp_LchixIU")->fill(nrig, lchixIU, weight);
		if (track->status[0][0] && sign<0) Hist::Head("hLn_LchixIU")->fill(nrig, lchixIU, weight);
		if (track->status[0][0]) Hist::Head("hLi_LchiyIU")->fill(irig, lchiyIU, weight);
		if (track->status[0][0] && sign>0) Hist::Head("hLp_LchiyIU")->fill(nrig, lchiyIU, weight);
		if (track->status[0][0] && sign<0) Hist::Head("hLn_LchiyIU")->fill(nrig, lchiyIU, weight);
		if (track->status[0][1]) Hist::Head("hLi_LchixIL")->fill(irig, lchixIL, weight);
		if (track->status[0][1] && sign>0) Hist::Head("hLp_LchixIL")->fill(nrig, lchixIL, weight);
		if (track->status[0][1] && sign<0) Hist::Head("hLn_LchixIL")->fill(nrig, lchixIL, weight);
		if (track->status[0][1]) Hist::Head("hLi_LchiyIL")->fill(irig, lchiyIL, weight);
		if (track->status[0][1] && sign>0) Hist::Head("hLp_LchiyIL")->fill(nrig, lchiyIL, weight);
		if (track->status[0][1] && sign<0) Hist::Head("hLn_LchiyIL")->fill(nrig, lchiyIL, weight);
		if (track->status[0][2]) Hist::Head("hLi_LchixIn")->fill(irig, lchixIn, weight);
		if (track->status[0][2] && sign>0) Hist::Head("hLp_LchixIn")->fill(nrig, lchixIn, weight);
		if (track->status[0][2] && sign<0) Hist::Head("hLn_LchixIn")->fill(nrig, lchixIn, weight);
		if (track->status[0][2]) Hist::Head("hLi_LchiyIn")->fill(irig, lchiyIn, weight);
		if (track->status[0][2] && sign>0) Hist::Head("hLp_LchiyIn")->fill(nrig, lchiyIn, weight);
		if (track->status[0][2] && sign<0) Hist::Head("hLn_LchiyIn")->fill(nrig, lchiyIn, weight);
		
		if (tvBin!=0 && track->status[0][0]) Hist::Head(StrFmt("hLi_LchixIU%03d", tvBin))->fill(irig, lchixIU, weight);
		if (tvBin!=0 && track->status[0][0] && sign>0) Hist::Head(StrFmt("hLp_LchixIU%03d", tvBin))->fill(nrig, lchixIU, weight);
		if (tvBin!=0 && track->status[0][0] && sign<0) Hist::Head(StrFmt("hLn_LchixIU%03d", tvBin))->fill(nrig, lchixIU, weight);
		if (tvBin!=0 && track->status[0][0]) Hist::Head(StrFmt("hLi_LchiyIU%03d", tvBin))->fill(irig, lchiyIU, weight);
		if (tvBin!=0 && track->status[0][0] && sign>0) Hist::Head(StrFmt("hLp_LchiyIU%03d", tvBin))->fill(nrig, lchiyIU, weight);
		if (tvBin!=0 && track->status[0][0] && sign<0) Hist::Head(StrFmt("hLn_LchiyIU%03d", tvBin))->fill(nrig, lchiyIU, weight);
		if (tvBin!=0 && track->status[0][1]) Hist::Head(StrFmt("hLi_LchixIL%03d", tvBin))->fill(irig, lchixIL, weight);
		if (tvBin!=0 && track->status[0][1] && sign>0) Hist::Head(StrFmt("hLp_LchixIL%03d", tvBin))->fill(nrig, lchixIL, weight);
		if (tvBin!=0 && track->status[0][1] && sign<0) Hist::Head(StrFmt("hLn_LchixIL%03d", tvBin))->fill(nrig, lchixIL, weight);
		if (tvBin!=0 && track->status[0][1]) Hist::Head(StrFmt("hLi_LchiyIL%03d", tvBin))->fill(irig, lchiyIL, weight);
		if (tvBin!=0 && track->status[0][1] && sign>0) Hist::Head(StrFmt("hLp_LchiyIL%03d", tvBin))->fill(nrig, lchiyIL, weight);
		if (tvBin!=0 && track->status[0][1] && sign<0) Hist::Head(StrFmt("hLn_LchiyIL%03d", tvBin))->fill(nrig, lchiyIL, weight);
		if (tvBin!=0 && track->status[0][2]) Hist::Head(StrFmt("hLi_LchixIn%03d", tvBin))->fill(irig, lchixIn, weight);
		if (tvBin!=0 && track->status[0][2] && sign>0) Hist::Head(StrFmt("hLp_LchixIn%03d", tvBin))->fill(nrig, lchixIn, weight);
		if (tvBin!=0 && track->status[0][2] && sign<0) Hist::Head(StrFmt("hLn_LchixIn%03d", tvBin))->fill(nrig, lchixIn, weight);
		if (tvBin!=0 && track->status[0][2]) Hist::Head(StrFmt("hLi_LchiyIn%03d", tvBin))->fill(irig, lchiyIn, weight);
		if (tvBin!=0 && track->status[0][2] && sign>0) Hist::Head(StrFmt("hLp_LchiyIn%03d", tvBin))->fill(nrig, lchiyIn, weight);
		if (tvBin!=0 && track->status[0][2] && sign<0) Hist::Head(StrFmt("hLn_LchiyIn%03d", tvBin))->fill(nrig, lchiyIn, weight);
		
		Hist::Head("hLi_Lchix")->fill(irig, lchix, weight);
		if (sign>0) Hist::Head("hLp_Lchix")->fill(nrig, lchix, weight);
		if (sign<0) Hist::Head("hLn_Lchix")->fill(nrig, lchix, weight);
		Hist::Head("hLi_Lchiy")->fill(irig, lchiy, weight);
		if (sign>0) Hist::Head("hLp_Lchiy")->fill(nrig, lchiy, weight);
		if (sign<0) Hist::Head("hLn_Lchiy")->fill(nrig, lchiy, weight);

		if (tvBin!=0) Hist::Head(StrFmt("hLi_Lchix%03d", tvBin))->fill(irig, lchix, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Lchix%03d", tvBin))->fill(nrig, lchix, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Lchix%03d", tvBin))->fill(nrig, lchix, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Lchiy%03d", tvBin))->fill(irig, lchiy, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Lchiy%03d", tvBin))->fill(nrig, lchiy, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Lchiy%03d", tvBin))->fill(nrig, lchiy, weight);

		// Charge Confusion
		Float_t ccest = std::max(lasym, lchiy);
		Hist::Head("hLi_CCest")->fill(irig, ccest, weight);
		if (sign>0) Hist::Head("hLp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hLn_CCest")->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_CCest%03d", tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_CCest%03d", tvBin))->fill(nrig, ccest, weight);
	
		//if (lchixIU > 2.0) break;
		//if (lchiyIU > 2.0 * MDRFT) break;
		//if (lchixIL > 2.0) break;
		//if (lchiyIL > 2.0 * MDRFT) break;
		if (lchixIn > 2.0) break;
		if (lchiyIn > 2.0 * MDRFT) break;
		if (lasymUL > 1.5 * MDRFT) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 4, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 4, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 4, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 4, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 4, weight);

		// RICH
		Bool_t  richgood  = (fRich->kindOfRad == 0 && fRich->isGoodTile && fRich->isInFiducialVolume);
		Int_t   richph    = fRich->numOfHit;
		Float_t richelph  = (!richgood) ? 0.0 : fRich->numOfExpPE[0];
		Float_t richprph  = (!richgood) ? 0.0 : fRich->numOfExpPE[3];

		Bool_t  hasrich   = (fRich->status);
		Float_t richb     = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbta);
		Float_t richbPi   = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbtaPi);
		Float_t richbEl   = (!hasrich) ? 0.0 : (1./fRich->beta - 1./tbtaEl);
		Float_t tofrichb  = (!hasrich) ? 0.0 : (fTof->betaH - fRich->beta);

		Bool_t  isrichPr   = ((richph <= 2) || (hasrich && richph >= 4 && std::fabs(richb) < 0.003)); 
		Bool_t  isrichPiEl = (hasrich && richph >= 4 && (richbEl+0.003) > 0. && (richbPi-0.003) < 0.);
		
		if (!richgood) break;

		Hist::Head("hLi_Cutflow")->fill(irig, 5, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 5, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 5, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 5, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		
		Hist::Head("hLi_Aglprph")->fill(irig, richprph, weight);
		if (sign>0) Hist::Head("hLp_Aglprph")->fill(nrig, richprph, weight);
		if (sign<0) Hist::Head("hLn_Aglprph")->fill(nrig, richprph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Aglprph%03d", tvBin))->fill(irig, richprph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Aglprph%03d", tvBin))->fill(nrig, richprph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Aglprph%03d", tvBin))->fill(nrig, richprph, weight);
		
		Hist::Head("hLi_Aglelph")->fill(irig, richelph, weight);
		if (sign>0) Hist::Head("hLp_Aglelph")->fill(nrig, richelph, weight);
		if (sign<0) Hist::Head("hLn_Aglelph")->fill(nrig, richelph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Aglelph%03d", tvBin))->fill(irig, richelph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Aglelph%03d", tvBin))->fill(nrig, richelph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Aglelph%03d", tvBin))->fill(nrig, richelph, weight);
		
		Hist::Head("hLi_Aglph")->fill(irig, richph, weight);
		if (sign>0) Hist::Head("hLp_Aglph")->fill(nrig, richph, weight);
		if (sign<0) Hist::Head("hLn_Aglph")->fill(nrig, richph, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Aglph%03d", tvBin))->fill(irig, richph, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Aglph%03d", tvBin))->fill(nrig, richph, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Aglph%03d", tvBin))->fill(nrig, richph, weight);
		
		if (hasrich) Hist::Head("hLi_Aglb")->fill(irig, richb, weight);
		if (hasrich && sign>0) Hist::Head("hLp_Aglb")->fill(nrig, richb, weight);
		if (hasrich && sign<0) Hist::Head("hLn_Aglb")->fill(nrig, richb, weight);
		if (tvBin!=0 && hasrich) Hist::Head(StrFmt("hLi_Aglb%03d", tvBin))->fill(irig, richb, weight);
		if (tvBin!=0 && hasrich && sign>0) Hist::Head(StrFmt("hLp_Aglb%03d", tvBin))->fill(nrig, richb, weight);
		if (tvBin!=0 && hasrich && sign<0) Hist::Head(StrFmt("hLn_Aglb%03d", tvBin))->fill(nrig, richb, weight);
		
		if (hasrich) Hist::Head("hLi_TofAglb")->fill(irig, tofrichb, weight);
		if (hasrich && sign>0) Hist::Head("hLp_TofAglb")->fill(nrig, tofrichb, weight);
		if (hasrich && sign<0) Hist::Head("hLn_TofAglb")->fill(nrig, tofrichb, weight);
		if (tvBin!=0 && hasrich) Hist::Head(StrFmt("hLi_TofAglb%03d", tvBin))->fill(irig, tofrichb, weight);
		if (tvBin!=0 && hasrich && sign>0) Hist::Head(StrFmt("hLp_TofAglb%03d", tvBin))->fill(nrig, tofrichb, weight);
		if (tvBin!=0 && hasrich && sign<0) Hist::Head(StrFmt("hLn_TofAglb%03d", tvBin))->fill(nrig, tofrichb, weight);

		if (richelph < 2.0) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 6, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 6, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 6, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 6, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 6, weight);
			
		// TOF
		Float_t tofb = (fTof->betaH - tbta);

		Bool_t isSignal     = (sign > 0 && isrichPr);
		Bool_t isBackground = (sign < 0 && isrichPiEl);
		if (isSignal)     Hist::Head("hLs_Tofb")->fill(nrig, tofb, weight);
		if (isBackground) Hist::Head("hLb_Tofb")->fill(nrig, tofb, weight);
		if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hLs_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hLb_Tofb%03d", tvBin))->fill(nrig, tofb, weight);
		
		if (!isrichPr) break;
		
		Hist::Head("hLi_Cutflow")->fill(irig, 7, weight);
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 7, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 7, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hLi_Cutflow%03d", tvBin))->fill(irig, 7, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hLp_Cutflow%03d", tvBin))->fill(nrig, 7, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hLn_Cutflow%03d", tvBin))->fill(nrig, 7, weight);
		
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
			if (MCTRig>0.0) Hist::Head("hLp_MCTrRso")->fill(MCNRig, ((MCTRig*irig)-1.0), weight);
			if (MCTRig<0.0) Hist::Head("hLn_MCTrRso")->fill(MCNRig, ((MCTRig*irig)-1.0), weight);
		
			Hist::Head("hLi_MCEvt")->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head("hLp_MCEvt")->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head("hLn_MCEvt")->fill(MCNRig, weight);
		}
	
		break;
	}


	/**** Intermedia Energy Region ****/
	while (true) {
		// TRK
		Int_t trPt = -1;
		std::string trNm = "";
		if      (track->status[0][5] && ((track->bitPatt& 162)==162)) { trPt = 5; trNm = "Fs"; }
		else if (track->status[0][4] && ((track->bitPatt& 130)==130)) { trPt = 4; trNm = "L9"; }
		else if (track->status[0][3] && ((track->bitPatt&  34)== 34)) { trPt = 3; trNm = "L1"; }
		else if (track->status[0][2] && track->status[0][0] && track->status[0][1]) { trPt = 2; trNm = "In"; }
		if (trPt == -1) break;
		if (!track->status[0][trPt]) break;
		Float_t trig       = track->rigidity[0][trPt];
		Int_t   sign       = MgntNum::Compare(trig);
		Float_t nrig       = std::fabs(trig);
		Float_t irig       = 1. / trig;
		Float_t tbta       = 1. / std::sqrt(1.+(mass/nrig)*(mass/nrig));
		Float_t tbtaPi     = 1. / std::sqrt(1.+(massPi/nrig)*(massPi/nrig));
		Float_t tbtaEl     = 1. / std::sqrt(1.+(massEl/nrig)*(massEl/nrig));
		
		// TODO
		//if (trPt == 2) break;
		//if (trPt != 2) break;
		
		// Cutoff
		if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
			Int_t icf = 0;
			Hist::Head("hIi_Cutoff")->fill(irig, fRti->cutoffIGRF[icf], weight);
			if (sign>0) Hist::Head("hIp_Cutoff")->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (sign<0) Hist::Head("hIn_Cutoff")->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutoff%03d", tvBin))->fill(irig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutoff%03d", tvBin))->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutoff%03d", tvBin))->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (nrig < CfStableFT * fRti->cutoffIGRF[icf]) break;
		}
		
		Hist::Head("hIi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		
		// TRD
		Int_t trdPt = 1;
		Bool_t  hastrd  = (fTrd->statusKCls[trdPt] && fTrd->LLR_nhit[trdPt] >= 8);
		if (!hastrd) break;
		
		Float_t trdl    = (fTrd->LLR[trdPt][0]);
		Bool_t  istrdHe = (fTrd->LLR[trdPt][2] > 0.3);
		if (istrdHe) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		// ECAL
		Bool_t  hasecal  = (shower != nullptr);
		Float_t ecalbdt  = (!hasecal) ? -2.0 : shower->PisaBDT;
		Bool_t  isecalPr = (!hasecal) ? false : (ecalbdt < -0.8);
		Bool_t  isecalEl = (!hasecal) ? false : (ecalbdt >  0.5);
		
		if (hasecal) Hist::Head("hIi_Ecalbdt")->fill(irig, ecalbdt, weight);
		if (hasecal && sign>0) Hist::Head("hIp_Ecalbdt")->fill(nrig, ecalbdt, weight);
		if (hasecal && sign<0) Hist::Head("hIn_Ecalbdt")->fill(nrig, ecalbdt, weight);
		if (tvBin!=0 && hasecal) Hist::Head(StrFmt("hIi_Ecalbdt%03d", tvBin))->fill(irig, ecalbdt, weight);
		if (tvBin!=0 && hasecal && sign>0) Hist::Head(StrFmt("hIp_Ecalbdt%03d", tvBin))->fill(nrig, ecalbdt, weight);
		if (tvBin!=0 && hasecal && sign<0) Hist::Head(StrFmt("hIn_Ecalbdt%03d", tvBin))->fill(nrig, ecalbdt, weight);
		
		// RICH
		Bool_t  richgood  = (fRich->status && 
		                     fRich->kindOfRad == 0 && 
												 fRich->isGoodTile && fRich->isInFiducialVolume);
		Float_t richb     = (!richgood) ? 0.0 : (1./fRich->beta - 1./tbta);
		Float_t richbPi   = (!richgood) ? 0.0 : (1./fRich->beta - 1./tbtaPi);
		Float_t richbEl   = (!richgood) ? 0.0 : (1./fRich->beta - 1./tbtaEl);
		Bool_t  isrichPr  = (!richgood) ? 0.0 : (std::fabs(richb) < 0.003);
		if (!richgood) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		
		Hist::Head("hIi_Aglb")->fill(irig, richb, weight);
		if (sign>0) Hist::Head("hIp_Aglb")->fill(nrig, richb, weight);
		if (sign<0) Hist::Head("hIn_Aglb")->fill(nrig, richb, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Aglb%03d", tvBin))->fill(irig, richb, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Aglb%03d", tvBin))->fill(nrig, richb, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Aglb%03d", tvBin))->fill(nrig, richb, weight);
		
		// Charge Confusion ---- TrSigma
		Float_t trParIU[3] = { 1.44615e-02, 1.58178e-03, 5.11515e-05 };
		Float_t trParIL[3] = { 9.60391e-03, 7.54520e-04, 3.81112e-05 };
		Float_t trParIn[3] = { 8.39123e-03, 3.36288e-04, 1.20139e-05 };
		Float_t trParL1[3] = { 7.64252e-03, 5.34561e-04, 1.36103e-06 };
		Float_t trParL9[3] = { 8.23647e-03, 2.97335e-04, 8.57550e-07 };
		Float_t trParFs[3] = { 9.29537e-03, 1.14609e-04, 2.50000e-07 };
		Float_t trSgmIU = ((!track->status[0][0]) ? 1.0 : 
		                  std::sqrt(trParIU[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIU[1] * nrig + trParIU[2] * nrig * nrig)) / nrig;
		Float_t trSgmIL = ((!track->status[0][1]) ? 1.0 : 
		                  std::sqrt(trParIL[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIL[1] * nrig + trParIL[2] * nrig * nrig)) / nrig;
		Float_t trSgmIn = ((!track->status[0][2]) ? 1.0 : 
		                  std::sqrt(trParIn[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIn[1] * nrig + trParIn[2] * nrig * nrig)) / nrig;
		Float_t trSgmL1 = ((!track->status[0][3]) ? 1.0 : 
		                  std::sqrt(trParL1[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParL1[1] * nrig + trParL1[2] * nrig * nrig)) / nrig;
		Float_t trSgmL9 = ((!track->status[0][4]) ? 1.0 : 
		                  std::sqrt(trParL9[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParL9[1] * nrig + trParL9[2] * nrig * nrig)) / nrig;
		Float_t trSgmFs = ((!track->status[0][5]) ? 1.0 : 
		                  std::sqrt(trParFs[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParFs[1] * nrig + trParFs[2] * nrig * nrig)) / nrig;
		
		Float_t trSgmPt = 1.0 / nrig;
		if      (trPt == 2) trSgmPt = trSgmIn;
		else if (trPt == 3) trSgmPt = trSgmL1;
		else if (trPt == 4) trSgmPt = trSgmL9;
		else if (trPt == 5) trSgmPt = trSgmFs;

		Float_t trMDRPt = nrig;
		if      (trPt == 2) trMDRPt = std::sqrt(trParIn[2]);
		else if (trPt == 3) trMDRPt = std::sqrt(trParL1[2]);
		else if (trPt == 4) trMDRPt = std::sqrt(trParL9[2]);
		else if (trPt == 5) trMDRPt = std::sqrt(trParFs[2]);
		Float_t MDRFT = std::erfc(5. * nrig * trMDRPt);

		Float_t trParIUSgm[4] = { trParIU[0] * mass * mass, trParIU[0], trParIU[1], trParIU[2] };
		Float_t trParILSgm[4] = { trParIL[0] * mass * mass, trParIL[0], trParIL[1], trParIL[2] };
		Float_t trParInSgm[4] = { trParIn[0] * mass * mass, trParIn[0], trParIn[1], trParIn[2] };
		Float_t trParL1Sgm[4] = { trParL1[0] * mass * mass, trParL1[0], trParL1[1], trParL1[2] };
		Float_t trParL9Sgm[4] = { trParL9[0] * mass * mass, trParL9[0], trParL9[1], trParL9[2] };
		Float_t trParFsSgm[4] = { trParFs[0] * mass * mass, trParFs[0], trParFs[1], trParFs[2] };
		Float_t trSgmIUSgm = (!track->status[0][0]) ? 0.0 : 
		                     (0.5 / trSgmIn) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParIUSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParIUSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParIUSgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmILSgm = (!track->status[0][1]) ? 0.0 : 
		                     (0.5 / trSgmIn) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParILSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParILSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParILSgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmInSgm = (!track->status[0][2]) ? 0.0 : 
		                     (0.5 / trSgmIn) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParInSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParInSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParInSgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmL1Sgm = (!track->status[0][3]) ? 0.0 : 
		                     (0.5 / trSgmL1) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParL1Sgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParL1Sgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParL1Sgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmL9Sgm = (!track->status[0][4]) ? 0.0 : 
		                     (0.5 / trSgmL9) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParL9Sgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParL9Sgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParL9Sgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmFsSgm = (!track->status[0][5]) ? 0.0 : 
		                     (0.5 / trSgmFs) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParFsSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParFsSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParFsSgm[2] * trSgmPt) )
		                     ); 
	
		if (track->status[0][0]) trSgmIU = std::sqrt(trSgmIU * trSgmIU + trSgmIUSgm * trSgmIUSgm);
		if (track->status[0][1]) trSgmIL = std::sqrt(trSgmIL * trSgmIL + trSgmILSgm * trSgmILSgm);
		if (track->status[0][2]) trSgmIn = std::sqrt(trSgmIn * trSgmIn + trSgmInSgm * trSgmInSgm);
		if (track->status[0][3]) trSgmL1 = std::sqrt(trSgmL1 * trSgmL1 + trSgmL1Sgm * trSgmL1Sgm);
		if (track->status[0][4]) trSgmL9 = std::sqrt(trSgmL9 * trSgmL9 + trSgmL9Sgm * trSgmL9Sgm);
		if (track->status[0][5]) trSgmFs = std::sqrt(trSgmFs * trSgmFs + trSgmFsSgm * trSgmFsSgm);

		// Charge Confusion ---- Asym
		const Float_t lasymLMT = 1.0e-8;
		const Float_t lasymSGM = 1.8;
	
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
	
		std::vector<Float_t> asymvec;
		if (hasasymUL) { asymvec.push_back(lasymUL);  }
		if (hasasym1I) { asymvec.push_back(lasym1I);  }
		if (hasasym9I) { asymvec.push_back(lasym9I);  }
		if (hasasym91) { asymvec.push_back(lasym91);  }
		Float_t lasym = (asymvec.size()==0) ? 0.0 : *std::max_element(std::begin(asymvec), std::end(asymvec));

		if (hasasymUL) Hist::Head("hIi_AsymUL")->fill(irig, asymUL, weight);
		if (hasasymUL && sign>0) Hist::Head("hIp_AsymUL")->fill(nrig, asymUL, weight);
		if (hasasymUL && sign<0) Hist::Head("hIn_AsymUL")->fill(nrig, asymUL, weight);
		if (hasasym1I) Hist::Head("hIi_Asym1I")->fill(irig, asym1I, weight);
		if (hasasym1I && sign>0) Hist::Head("hIp_Asym1I")->fill(nrig, asym1I, weight);
		if (hasasym1I && sign<0) Hist::Head("hIn_Asym1I")->fill(nrig, asym1I, weight);
		if (hasasym9I) Hist::Head("hIi_Asym9I")->fill(irig, asym9I, weight);
		if (hasasym9I && sign>0) Hist::Head("hIp_Asym9I")->fill(nrig, asym9I, weight);
		if (hasasym9I && sign<0) Hist::Head("hIn_Asym9I")->fill(nrig, asym9I, weight);
		if (hasasym91) Hist::Head("hIi_Asym91")->fill(irig, asym91, weight);
		if (hasasym91 && sign>0) Hist::Head("hIp_Asym91")->fill(nrig, asym91, weight);
		if (hasasym91 && sign<0) Hist::Head("hIn_Asym91")->fill(nrig, asym91, weight);
		
		if (tvBin!=0 && hasasymUL) Hist::Head(StrFmt("hIi_AsymUL%03d", tvBin))->fill(irig, asymUL, weight);
		if (tvBin!=0 && hasasymUL && sign>0) Hist::Head(StrFmt("hIp_AsymUL%03d", tvBin))->fill(nrig, asymUL, weight);
		if (tvBin!=0 && hasasymUL && sign<0) Hist::Head(StrFmt("hIn_AsymUL%03d", tvBin))->fill(nrig, asymUL, weight);
		if (tvBin!=0 && hasasym1I) Hist::Head(StrFmt("hIi_Asym1I%03d", tvBin))->fill(irig, asym1I, weight);
		if (tvBin!=0 && hasasym1I && sign>0) Hist::Head(StrFmt("hIp_Asym1I%03d", tvBin))->fill(nrig, asym1I, weight);
		if (tvBin!=0 && hasasym1I && sign<0) Hist::Head(StrFmt("hIn_Asym1I%03d", tvBin))->fill(nrig, asym1I, weight);
		if (tvBin!=0 && hasasym9I) Hist::Head(StrFmt("hIi_Asym9I%03d", tvBin))->fill(irig, asym9I, weight);
		if (tvBin!=0 && hasasym9I && sign>0) Hist::Head(StrFmt("hIp_Asym9I%03d", tvBin))->fill(nrig, asym9I, weight);
		if (tvBin!=0 && hasasym9I && sign<0) Hist::Head(StrFmt("hIn_Asym9I%03d", tvBin))->fill(nrig, asym9I, weight);
		if (tvBin!=0 && hasasym91) Hist::Head(StrFmt("hIi_Asym91%03d", tvBin))->fill(irig, asym91, weight);
		if (tvBin!=0 && hasasym91 && sign>0) Hist::Head(StrFmt("hIp_Asym91%03d", tvBin))->fill(nrig, asym91, weight);
		if (tvBin!=0 && hasasym91 && sign<0) Hist::Head(StrFmt("hIn_Asym91%03d", tvBin))->fill(nrig, asym91, weight);
		
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
			
		Hist::Head("hIi_Lasym")->fill(irig, lasym, weight);
		if (sign>0) Hist::Head("hIp_Lasym")->fill(nrig, lasym, weight);
		if (sign<0) Hist::Head("hIn_Lasym")->fill(nrig, lasym, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Lasym", tvBin))->fill(irig, lasym, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Lasym", tvBin))->fill(nrig, lasym, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Lasym", tvBin))->fill(nrig, lasym, weight);

		// Charge Confusion ---- ChisqCC
		Float_t lchixIU = (!track->status[0][0]) ? 0.0 : std::log(track->chisq[0][0][0]);
		Float_t lchiyIU = (!track->status[0][0]) ? 0.0 : std::log(track->chisq[0][0][1]);
		Float_t lchixIL = (!track->status[0][1]) ? 0.0 : std::log(track->chisq[0][1][0]);
		Float_t lchiyIL = (!track->status[0][1]) ? 0.0 : std::log(track->chisq[0][1][1]);
		Float_t lchixIn = (!track->status[0][2]) ? 0.0 : std::log(track->chisq[0][2][0]);
		Float_t lchiyIn = (!track->status[0][2]) ? 0.0 : std::log(track->chisq[0][2][1]);
		Float_t lchixL1 = (!track->status[0][3]) ? 0.0 : std::log(track->chisq[0][3][0]);
		Float_t lchiyL1 = (!track->status[0][3]) ? 0.0 : std::log(track->chisq[0][3][1]);
		Float_t lchixL9 = (!track->status[0][4]) ? 0.0 : std::log(track->chisq[0][4][0]);
		Float_t lchiyL9 = (!track->status[0][4]) ? 0.0 : std::log(track->chisq[0][4][1]);
		Float_t lchixFs = (!track->status[0][5]) ? 0.0 : std::log(track->chisq[0][5][0]);
		Float_t lchiyFs = (!track->status[0][5]) ? 0.0 : std::log(track->chisq[0][5][1]);
		
		std::vector<Float_t> chixvec;
		if (track->status[0][0]) { chixvec.push_back(lchixIU); }
		if (track->status[0][1]) { chixvec.push_back(lchixIL); }
		if (track->status[0][2]) { chixvec.push_back(lchixIn); }
		if (track->status[0][3]) { chixvec.push_back(lchixL1); }
		if (track->status[0][4]) { chixvec.push_back(lchixL9); }
		if (track->status[0][5]) { chixvec.push_back(lchixFs); }
		Float_t lchix = (chixvec.size()==0) ? 0.0 : *std::max_element(std::begin(chixvec), std::end(chixvec));
		
		std::vector<Float_t> chiyvec;
		if (track->status[0][0]) { chiyvec.push_back(lchiyIU); }
		if (track->status[0][1]) { chiyvec.push_back(lchiyIL); }
		if (track->status[0][2]) { chiyvec.push_back(lchiyIn); }
		if (track->status[0][3]) { chiyvec.push_back(lchiyL1); }
		if (track->status[0][4]) { chiyvec.push_back(lchiyL9); }
		if (track->status[0][5]) { chiyvec.push_back(lchiyFs); }
		Float_t lchiy = (chiyvec.size()==0) ? 0.0 : *std::max_element(std::begin(chiyvec), std::end(chiyvec));

		if (track->status[0][0]) Hist::Head("hIi_LchixIU")->fill(irig, lchixIU, weight);
		if (track->status[0][0] && sign>0) Hist::Head("hIp_LchixIU")->fill(nrig, lchixIU, weight);
		if (track->status[0][0] && sign<0) Hist::Head("hIn_LchixIU")->fill(nrig, lchixIU, weight);
		if (track->status[0][0]) Hist::Head("hIi_LchiyIU")->fill(irig, lchiyIU, weight);
		if (track->status[0][0] && sign>0) Hist::Head("hIp_LchiyIU")->fill(nrig, lchiyIU, weight);
		if (track->status[0][0] && sign<0) Hist::Head("hIn_LchiyIU")->fill(nrig, lchiyIU, weight);
		if (track->status[0][1]) Hist::Head("hIi_LchixIL")->fill(irig, lchixIL, weight);
		if (track->status[0][1] && sign>0) Hist::Head("hIp_LchixIL")->fill(nrig, lchixIL, weight);
		if (track->status[0][1] && sign<0) Hist::Head("hIn_LchixIL")->fill(nrig, lchixIL, weight);
		if (track->status[0][1]) Hist::Head("hIi_LchiyIL")->fill(irig, lchiyIL, weight);
		if (track->status[0][1] && sign>0) Hist::Head("hIp_LchiyIL")->fill(nrig, lchiyIL, weight);
		if (track->status[0][1] && sign<0) Hist::Head("hIn_LchiyIL")->fill(nrig, lchiyIL, weight);
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
		
		if (tvBin!=0 && track->status[0][0]) Hist::Head(StrFmt("hIi_LchixIU%03d", tvBin))->fill(irig, lchixIU, weight);
		if (tvBin!=0 && track->status[0][0] && sign>0) Hist::Head(StrFmt("hIp_LchixIU%03d", tvBin))->fill(nrig, lchixIU, weight);
		if (tvBin!=0 && track->status[0][0] && sign<0) Hist::Head(StrFmt("hIn_LchixIU%03d", tvBin))->fill(nrig, lchixIU, weight);
		if (tvBin!=0 && track->status[0][0]) Hist::Head(StrFmt("hIi_LchiyIU%03d", tvBin))->fill(irig, lchiyIU, weight);
		if (tvBin!=0 && track->status[0][0] && sign>0) Hist::Head(StrFmt("hIp_LchiyIU%03d", tvBin))->fill(nrig, lchiyIU, weight);
		if (tvBin!=0 && track->status[0][0] && sign<0) Hist::Head(StrFmt("hIn_LchiyIU%03d", tvBin))->fill(nrig, lchiyIU, weight);
		if (tvBin!=0 && track->status[0][1]) Hist::Head(StrFmt("hIi_LchixIL%03d", tvBin))->fill(irig, lchixIL, weight);
		if (tvBin!=0 && track->status[0][1] && sign>0) Hist::Head(StrFmt("hIp_LchixIL%03d", tvBin))->fill(nrig, lchixIL, weight);
		if (tvBin!=0 && track->status[0][1] && sign<0) Hist::Head(StrFmt("hIn_LchixIL%03d", tvBin))->fill(nrig, lchixIL, weight);
		if (tvBin!=0 && track->status[0][1]) Hist::Head(StrFmt("hIi_LchiyIL%03d", tvBin))->fill(irig, lchiyIL, weight);
		if (tvBin!=0 && track->status[0][1] && sign>0) Hist::Head(StrFmt("hIp_LchiyIL%03d", tvBin))->fill(nrig, lchiyIL, weight);
		if (tvBin!=0 && track->status[0][1] && sign<0) Hist::Head(StrFmt("hIn_LchiyIL%03d", tvBin))->fill(nrig, lchiyIL, weight);
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
			
		Hist::Head("hIi_Lchix")->fill(irig, lchix, weight);
		if (sign>0) Hist::Head("hIp_Lchix")->fill(nrig, lchix, weight);
		if (sign<0) Hist::Head("hIn_Lchix")->fill(nrig, lchix, weight);
		Hist::Head("hIi_Lchiy")->fill(irig, lchiy, weight);
		if (sign>0) Hist::Head("hIp_Lchiy")->fill(nrig, lchiy, weight);
		if (sign<0) Hist::Head("hIn_Lchiy")->fill(nrig, lchiy, weight);
		
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Lchix%03d", tvBin))->fill(irig, lchix, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Lchix%03d", tvBin))->fill(nrig, lchix, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Lchix%03d", tvBin))->fill(nrig, lchix, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Lchiy%03d", tvBin))->fill(irig, lchiy, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Lchiy%03d", tvBin))->fill(nrig, lchiy, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Lchiy%03d", tvBin))->fill(nrig, lchiy, weight);
		
		// Charge Confusion
		Float_t ccest = std::max(lasym, lchiy);
		Hist::Head("hIi_CCest")->fill(irig, ccest, weight);
		if (sign>0) Hist::Head("hIp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hIn_CCest")->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_CCest%03d", tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_CCest%03d", tvBin))->fill(nrig, ccest, weight);

		if (trPt == 5 && (lchixL1 > 2.0 || lchixL9 > 2.0 || lchixFs > 2.0)) break;
		if (trPt == 5 && (lchiyL1 > 2.0 * MDRFT || lchiyL9 > 2.0 * MDRFT || lchiyFs > 2.0 * MDRFT)) break;
		if (trPt == 5 && (lasym1I > 1.5 * MDRFT || lasym9I > 1.5 * MDRFT || lasym91 > 1.5 * MDRFT)) break;
		if (trPt == 4 && (lchixIn > 2.0 || lchixL9 > 2.0)) break;
		if (trPt == 4 && (lchiyIn > 2.0 * MDRFT || lchiyL9 > 2.0 * MDRFT)) break;
		if (trPt == 4 && (lasymUL > 1.5 * MDRFT || lasym9I > 1.5 * MDRFT)) break;
		if (trPt == 3 && (lchixIn > 2.0 || lchixL1 > 2.0)) break;
		if (trPt == 3 && (lchiyIn > 2.0 * MDRFT || lchiyL1 > 2.0 * MDRFT)) break;
		if (trPt == 3 && (lasymUL > 1.5 * MDRFT || lasym1I > 1.5 * MDRFT)) break;
		if (trPt == 2 && (lchixIU > 2.0 || lchixIL > 2.0 || lchixIn > 2.0)) break;
		if (trPt == 2 && (lchiyIU > 2.0 * MDRFT || lchiyIL > 2.0 * MDRFT || lchiyIn > 2.0 * MDRFT)) break;
		if (trPt == 2 && (lasymUL > 1.5 * MDRFT)) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);

		// ECAL
		Bool_t isSignal     = (sign > 0 && isecalPr);
		Bool_t isBackground = (sign < 0 && isecalEl);
		if (isSignal)     Hist::Head("hIs_Trdl")->fill(nrig, trdl, weight);
		if (isBackground) Hist::Head("hIb_Trdl")->fill(nrig, trdl, weight);
		if (tvBin!=0 && isSignal)     Hist::Head(StrFmt("hIs_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && isBackground) Hist::Head(StrFmt("hIb_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		
		if (hasecal && !isecalPr) break;
		
		Hist::Head("hIi_Cutflow")->fill(irig, 4, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 4, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 4, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 4, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		
		if (!isrichPr) break;

		Hist::Head("hIi_Cutflow")->fill(irig, 5, weight);
		if (sign>0) Hist::Head("hIp_Cutflow")->fill(nrig, 5, weight);
		if (sign<0) Hist::Head("hIn_Cutflow")->fill(nrig, 5, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Cutflow%03d", tvBin))->fill(irig, 5, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Cutflow%03d", tvBin))->fill(nrig, 5, weight);
		
		Hist::Head("hIi_Trdl")->fill(irig, trdl, weight);
		if (sign>0) Hist::Head("hIp_Trdl")->fill(nrig, trdl, weight);
		if (sign<0) Hist::Head("hIn_Trdl")->fill(nrig, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		
		Hist::Head(StrFmt("hIi_%s_Trdl", trNm.c_str()))->fill(irig, trdl, weight);
		if (sign>0) Hist::Head(StrFmt("hIp_%s_Trdl", trNm.c_str()))->fill(nrig, trdl, weight);
		if (sign<0) Hist::Head(StrFmt("hIn_%s_Trdl", trNm.c_str()))->fill(nrig, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_%s_Trdl%03d", trNm.c_str(), tvBin))->fill(irig, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_%s_Trdl%03d", trNm.c_str(), tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_%s_Trdl%03d", trNm.c_str(), tvBin))->fill(nrig, trdl, weight);
		
		if (tvBin!=0) Hist::Head("hIi_Evt")->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head("hIp_Evt")->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head("hIn_Evt")->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hIi_Evt%03d", tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hIp_Evt%03d", tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hIn_Evt%03d", tvBin))->fill(nrig, weight);
		
		if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
			Hist::Head("hIi_MCTrRso")->fill(MCIRig, irig, weight);
			if (MCTRig>0.0) Hist::Head("hIp_MCTrRso")->fill(MCNRig, ((MCTRig*irig)-1.0), weight);
			if (MCTRig<0.0) Hist::Head("hIn_MCTrRso")->fill(MCNRig, ((MCTRig*irig)-1.0), weight);
		
			Hist::Head("hIi_MCEvt")->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head("hIp_MCEvt")->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head("hIn_MCEvt")->fill(MCNRig, weight);
		}

		break;
	}

	/**** High Energy Region ****/
	while (true) {
		// TRK
		Int_t trPt = -1;
		std::string trNm = "";
		if      (track->status[0][5] && ((track->bitPatt& 162)==162)) { trPt = 5; trNm = "Fs"; }
		else if (track->status[0][4] && ((track->bitPatt& 130)==130)) { trPt = 4; trNm = "L9"; }
		else if (track->status[0][3] && ((track->bitPatt&  34)== 34)) { trPt = 3; trNm = "L1"; }
		if (trPt == -1) break;
		if (!track->status[0][trPt]) break;
		Float_t trig = track->rigidity[0][trPt];
		Int_t   sign = MgntNum::Compare(trig);
		Float_t nrig = std::fabs(trig);
		Float_t irig = 1. / trig;

		// TODO
		// (93.00, 108.00000000, 125.00000000, 147.00000000, 175.00000000, 211.00000000, 259.00000000, 450.00000000)
		//const Float_t cutRigL1 = 93.00;
		//const Float_t cutRigL9 = 93.00;
		//if (trPt == 3 && nrig > cutRigL1) break;
		//if (trPt == 4 && nrig > cutRigL9) break;
		//if (trPt != 5) break;	
		
		// Cutoff
		if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
			Int_t icf = 0;
			Hist::Head("hHi_Cutoff")->fill(irig, fRti->cutoffIGRF[icf], weight);
			if (sign>0) Hist::Head("hHp_Cutoff")->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (sign<0) Hist::Head("hHn_Cutoff")->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutoff%03d", tvBin))->fill(irig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutoff%03d", tvBin))->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutoff%03d", tvBin))->fill(nrig, fRti->cutoffIGRF[icf], weight);
			if (nrig < CfStableFT * fRti->cutoffIGRF[icf]) break;
		}
		
		Hist::Head("hHi_Cutflow")->fill(irig, 0, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 0, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 0, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 0, weight);
		
		// TRD
		Int_t trdPt = 1;
		Bool_t  hastrd  = (fTrd->statusKCls[trdPt] && fTrd->LLR_nhit[trdPt] >= 8);
		if (!hastrd) break;
		
		Float_t trdl    = (fTrd->LLR[trdPt][0]);
		Bool_t  istrdHe = (fTrd->LLR[trdPt][2] > 0.3);
		Bool_t  istrdPr = (trdl > 0.75);
		if (istrdHe) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 1, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 1, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 1, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 1, weight);
		
		Hist::Head("hHi_Trdl")->fill(irig, trdl, weight);
		if (sign>0) Hist::Head("hHp_Trdl")->fill(nrig, trdl, weight);
		if (sign<0) Hist::Head("hHn_Trdl")->fill(nrig, trdl, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Trdl%03d", tvBin))->fill(irig, trdl, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Trdl%03d", tvBin))->fill(nrig, trdl, weight);
		
		if (!istrdPr) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 2, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 2, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 2, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 2, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 2, weight);

		// ECAL
		Bool_t  hasecal  = (shower != nullptr);
		Float_t ecalbdt  = (!hasecal) ? -2.0 : shower->PisaBDT;
		Bool_t  isecalPr = (!hasecal) ? false : (ecalbdt < -0.8);
		Bool_t  isecalEl = (!hasecal) ? false : (ecalbdt >  0.5);

		if (hasecal) Hist::Head("hHi_Ecalbdt")->fill(irig, ecalbdt, weight);
		if (hasecal && sign>0) Hist::Head("hHp_Ecalbdt")->fill(nrig, ecalbdt, weight);
		if (hasecal && sign<0) Hist::Head("hHn_Ecalbdt")->fill(nrig, ecalbdt, weight);
		if (tvBin!=0 && hasecal) Hist::Head(StrFmt("hHi_Ecalbdt%03d", tvBin))->fill(irig, ecalbdt, weight);
		if (tvBin!=0 && hasecal && sign>0) Hist::Head(StrFmt("hHp_Ecalbdt%03d", tvBin))->fill(nrig, ecalbdt, weight);
		if (tvBin!=0 && hasecal && sign<0) Hist::Head(StrFmt("hHn_Ecalbdt%03d", tvBin))->fill(nrig, ecalbdt, weight);
		
		if (hasecal && !isecalPr) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 3, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 3, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 3, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 3, weight);
		
		// Charge Confusion ---- TrSigma
		Float_t trParIn[3] = { 8.39123e-03, 3.36288e-04, 1.20139e-05 };
		Float_t trParL1[3] = { 7.64252e-03, 5.34561e-04, 1.36103e-06 };
		Float_t trParL9[3] = { 8.23647e-03, 2.97335e-04, 8.57550e-07 };
		Float_t trParFs[3] = { 9.29537e-03, 1.14609e-04, 2.50000e-07 };
		Float_t trSgmIn = ((!track->status[0][2]) ? 1.0 : 
		                  std::sqrt(trParIn[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParIn[1] * nrig + trParIn[2] * nrig * nrig)) / nrig;
		Float_t trSgmL1 = ((!track->status[0][3]) ? 1.0 : 
		                  std::sqrt(trParL1[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParL1[1] * nrig + trParL1[2] * nrig * nrig)) / nrig;
		Float_t trSgmL9 = ((!track->status[0][4]) ? 1.0 : 
		                  std::sqrt(trParL9[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParL9[1] * nrig + trParL9[2] * nrig * nrig)) / nrig;
		Float_t trSgmFs = ((!track->status[0][5]) ? 1.0 : 
		                  std::sqrt(trParFs[0] * (1.0 + mass*mass/nrig/nrig) + 
											trParFs[1] * nrig + trParFs[2] * nrig * nrig)) / nrig;
		
		Float_t trSgmPt = 1.0 / nrig;
		if      (trPt == 2) trSgmPt = trSgmIn;
		else if (trPt == 3) trSgmPt = trSgmL1;
		else if (trPt == 4) trSgmPt = trSgmL9;
		else if (trPt == 5) trSgmPt = trSgmFs;

		Float_t trMDRPt = nrig;
		if      (trPt == 2) trMDRPt = std::sqrt(trParIn[2]);
		else if (trPt == 3) trMDRPt = std::sqrt(trParL1[2]);
		else if (trPt == 4) trMDRPt = std::sqrt(trParL9[2]);
		else if (trPt == 5) trMDRPt = std::sqrt(trParFs[2]);
		Float_t MDRFT = std::erfc(5. * nrig * trMDRPt);
		
		Float_t trParInSgm[4] = { trParIn[0] * mass * mass, trParIn[0], trParIn[1], trParIn[2] };
		Float_t trParL1Sgm[4] = { trParL1[0] * mass * mass, trParL1[0], trParL1[1], trParL1[2] };
		Float_t trParL9Sgm[4] = { trParL9[0] * mass * mass, trParL9[0], trParL9[1], trParL9[2] };
		Float_t trParFsSgm[4] = { trParFs[0] * mass * mass, trParFs[0], trParFs[1], trParFs[2] };
		Float_t trSgmInSgm = (!track->status[0][2]) ? 0.0 : 
		                     (0.5 / trSgmIn) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParInSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParInSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParInSgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmL1Sgm = (!track->status[0][3]) ? 0.0 : 
		                     (0.5 / trSgmL1) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParL1Sgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParL1Sgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParL1Sgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmL9Sgm = (!track->status[0][4]) ? 0.0 : 
		                     (0.5 / trSgmL9) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParL9Sgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParL9Sgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParL9Sgm[2] * trSgmPt) )
		                     ); 
		Float_t trSgmFsSgm = (!track->status[0][5]) ? 0.0 : 
		                     (0.5 / trSgmFs) * std::sqrt(
		                     ( MgntNum::Square(4.0 * trParFsSgm[0] * trSgmPt / nrig / nrig / nrig) +
												   MgntNum::Square(2.0 * trParFsSgm[1] * trSgmPt / nrig) +
													 MgntNum::Square(1.0 * trParFsSgm[2] * trSgmPt) )
		                     ); 
	
		if (track->status[0][2]) trSgmIn = std::sqrt(trSgmIn * trSgmIn + trSgmInSgm * trSgmInSgm);
		if (track->status[0][3]) trSgmL1 = std::sqrt(trSgmL1 * trSgmL1 + trSgmL1Sgm * trSgmL1Sgm);
		if (track->status[0][4]) trSgmL9 = std::sqrt(trSgmL9 * trSgmL9 + trSgmL9Sgm * trSgmL9Sgm);
		if (track->status[0][5]) trSgmFs = std::sqrt(trSgmFs * trSgmFs + trSgmFsSgm * trSgmFsSgm);

		// Charge Confusion ---- Asym
		const Float_t lasymLMT = 1.0e-8;
		const Float_t lasymSGM = 1.8;

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

		Float_t asymwgt = 0.0;
		Float_t lasym   = 0.0;
		if (hasasym1I) { lasym += (lasym1I/trSgm1I/trSgm1I); asymwgt += (1.0/trSgm1I/trSgm1I); }
		if (hasasym9I) { lasym += (lasym9I/trSgm9I/trSgm9I); asymwgt += (1.0/trSgm9I/trSgm9I); }
		if (hasasym91) { lasym += (lasym91/trSgm91/trSgm91); asymwgt += (1.0/trSgm91/trSgm91); }
		lasym = (lasym / asymwgt);

		if (hasasym1I) Hist::Head("hHi_Asym1I")->fill(irig, asym1I, weight);
		if (hasasym1I && sign>0) Hist::Head("hHp_Asym1I")->fill(nrig, asym1I, weight);
		if (hasasym1I && sign<0) Hist::Head("hHn_Asym1I")->fill(nrig, asym1I, weight);
		if (hasasym9I) Hist::Head("hHi_Asym9I")->fill(irig, asym9I, weight);
		if (hasasym9I && sign>0) Hist::Head("hHp_Asym9I")->fill(nrig, asym9I, weight);
		if (hasasym9I && sign<0) Hist::Head("hHn_Asym9I")->fill(nrig, asym9I, weight);
		if (hasasym91) Hist::Head("hHi_Asym91")->fill(irig, asym91, weight);
		if (hasasym91 && sign>0) Hist::Head("hHp_Asym91")->fill(nrig, asym91, weight);
		if (hasasym91 && sign<0) Hist::Head("hHn_Asym91")->fill(nrig, asym91, weight);
		
		if (tvBin!=0 && hasasym1I) Hist::Head(StrFmt("hHi_Asym1I%03d", tvBin))->fill(irig, asym1I, weight);
		if (tvBin!=0 && hasasym1I && sign>0) Hist::Head(StrFmt("hHp_Asym1I%03d", tvBin))->fill(nrig, asym1I, weight);
		if (tvBin!=0 && hasasym1I && sign<0) Hist::Head(StrFmt("hHn_Asym1I%03d", tvBin))->fill(nrig, asym1I, weight);
		if (tvBin!=0 && hasasym9I) Hist::Head(StrFmt("hHi_Asym9I%03d", tvBin))->fill(irig, asym9I, weight);
		if (tvBin!=0 && hasasym9I && sign>0) Hist::Head(StrFmt("hHp_Asym9I%03d", tvBin))->fill(nrig, asym9I, weight);
		if (tvBin!=0 && hasasym9I && sign<0) Hist::Head(StrFmt("hHn_Asym9I%03d", tvBin))->fill(nrig, asym9I, weight);
		if (tvBin!=0 && hasasym91) Hist::Head(StrFmt("hHi_Asym91%03d", tvBin))->fill(irig, asym91, weight);
		if (tvBin!=0 && hasasym91 && sign>0) Hist::Head(StrFmt("hHp_Asym91%03d", tvBin))->fill(nrig, asym91, weight);
		if (tvBin!=0 && hasasym91 && sign<0) Hist::Head(StrFmt("hHn_Asym91%03d", tvBin))->fill(nrig, asym91, weight);
		
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
			
		Hist::Head("hHi_Lasym")->fill(irig, lasym, weight);
		if (sign>0) Hist::Head("hHp_Lasym")->fill(nrig, lasym, weight);
		if (sign<0) Hist::Head("hHn_Lasym")->fill(nrig, lasym, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Lasym%03d", tvBin))->fill(irig, lasym, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Lasym%03d", tvBin))->fill(nrig, lasym, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Lasym%03d", tvBin))->fill(nrig, lasym, weight);

		// Charge Confusion ---- ChisqCC
		Float_t lchixL1 = (!track->status[0][3]) ? 0.0 : std::log(track->chisq[0][3][0]);
		Float_t lchiyL1 = (!track->status[0][3]) ? 0.0 : std::log(track->chisq[0][3][1]);
		Float_t lchixL9 = (!track->status[0][4]) ? 0.0 : std::log(track->chisq[0][4][0]);
		Float_t lchiyL9 = (!track->status[0][4]) ? 0.0 : std::log(track->chisq[0][4][1]);
		Float_t lchixFs = (!track->status[0][5]) ? 0.0 : std::log(track->chisq[0][5][0]);
		Float_t lchiyFs = (!track->status[0][5]) ? 0.0 : std::log(track->chisq[0][5][1]);
		
		Float_t chixwgt = 0.0;
		Float_t lchix   = 0.0;
		if (track->status[0][3]) { lchix += (lchixL1/trSgmL1/trSgmL1); chixwgt += (1.0/trSgmL1/trSgmL1); }
		if (track->status[0][4]) { lchix += (lchixL9/trSgmL9/trSgmL9); chixwgt += (1.0/trSgmL9/trSgmL9); }
		if (track->status[0][5]) { lchix += (lchixFs/trSgmFs/trSgmFs); chixwgt += (1.0/trSgmFs/trSgmFs); }
		lchix = (lchix / chixwgt);

		Float_t chiywgt = 0.0;
		Float_t lchiy   = 0.0;
		if (track->status[0][3]) { lchiy += (lchiyL1/trSgmL1/trSgmL1); chiywgt += (1.0/trSgmL1/trSgmL1); }
		if (track->status[0][4]) { lchiy += (lchiyL9/trSgmL9/trSgmL9); chiywgt += (1.0/trSgmL9/trSgmL9); }
		if (track->status[0][5]) { lchiy += (lchiyFs/trSgmFs/trSgmFs); chiywgt += (1.0/trSgmFs/trSgmFs); }
		lchiy = (lchiy / chiywgt);

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
			
		Hist::Head("hHi_Lchix")->fill(irig, lchix, weight);
		if (sign>0) Hist::Head("hHp_Lchix")->fill(nrig, lchix, weight);
		if (sign<0) Hist::Head("hHn_Lchix")->fill(nrig, lchix, weight);
		Hist::Head("hHi_Lchiy")->fill(irig, lchiy, weight);
		if (sign>0) Hist::Head("hHp_Lchiy")->fill(nrig, lchiy, weight);
		if (sign<0) Hist::Head("hHn_Lchiy")->fill(nrig, lchiy, weight);
		
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Lchix%03d", tvBin))->fill(irig, lchix, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Lchix%03d", tvBin))->fill(nrig, lchix, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Lchix%03d", tvBin))->fill(nrig, lchix, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Lchiy%03d", tvBin))->fill(irig, lchiy, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Lchiy%03d", tvBin))->fill(nrig, lchiy, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Lchiy%03d", tvBin))->fill(nrig, lchiy, weight);
		
		if (trPt == 5 && (lchixL1 > 2.0 || lchixL9 > 2.0 || lchixFs > 2.0)) break;
		if (trPt == 5 && (lchiyL1 > 2.0 || lchiyL9 > 2.0 || lchiyFs > 2.0)) break;
		if (trPt == 5 && (lasym1I > 1.5 || lasym9I > 1.5 || lasym91 > 1.5)) break;
		if (trPt == 4 && (lchixL9 > 2.0)) break;
		if (trPt == 4 && (lchiyL9 > 2.0)) break;
		if (trPt == 4 && (lasym9I > 1.5)) break;
		if (trPt == 3 && (lchixL1 > 2.0)) break;
		if (trPt == 3 && (lchiyL1 > 2.0)) break;
		if (trPt == 3 && (lasym1I > 1.5)) break;
		
		Hist::Head("hHi_Cutflow")->fill(irig, 4, weight);
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 4, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 4, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Cutflow%03d", tvBin))->fill(irig, 4, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Cutflow%03d", tvBin))->fill(nrig, 4, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Cutflow%03d", tvBin))->fill(nrig, 4, weight);

		// Charge Confusion
		Float_t ccwgt = (asymwgt + chiywgt);
		Float_t ccest = (lasym*asymwgt + lchiy*chiywgt) / ccwgt;
		Hist::Head("hHi_CCest")->fill(irig, ccest, weight);
		if (sign>0) Hist::Head("hHp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hHn_CCest")->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_CCest%03d", tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_CCest%03d", tvBin))->fill(nrig, ccest, weight);
		
		Hist::Head(StrFmt("hHi_%s_CCest", trNm.c_str()))->fill(irig, ccest, weight);
		if (sign>0) Hist::Head(StrFmt("hHp_%s_CCest", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_%s_CCest", trNm.c_str()))->fill(nrig, ccest, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_%s_CCest%03d", trNm.c_str(), tvBin))->fill(irig, ccest, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_%s_CCest%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_%s_CCest%03d", trNm.c_str(), tvBin))->fill(nrig, ccest, weight);
		
		Hist::Head("hHi_Evt")->fill(irig, weight);
		if (sign>0) Hist::Head("hHp_Evt")->fill(nrig, weight);
		if (sign<0) Hist::Head("hHn_Evt")->fill(nrig, weight);
		if (tvBin!=0) Hist::Head(StrFmt("hHi_Evt%03d", tvBin))->fill(irig, weight);
		if (tvBin!=0 && sign>0) Hist::Head(StrFmt("hHp_Evt%03d", tvBin))->fill(nrig, weight);
		if (tvBin!=0 && sign<0) Hist::Head(StrFmt("hHn_Evt%03d", tvBin))->fill(nrig, weight);
		
		if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
			Hist::Head("hHi_MCTrRso")->fill(MCIRig, irig, weight);
			if (MCTRig>0.0) Hist::Head("hHp_MCTrRso")->fill(MCNRig, ((MCTRig*irig)-1.0), weight);
			if (MCTRig<0.0) Hist::Head("hHn_MCTrRso")->fill(MCNRig, ((MCTRig*irig)-1.0), weight);
		
			Hist::Head("hHi_MCEvt")->fill(MCIRig, weight);
			if (MCTRig>0.0) Hist::Head("hHp_MCEvt")->fill(MCNRig, weight);
			if (MCTRig<0.0) Hist::Head("hHn_MCEvt")->fill(MCNRig, weight);
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
