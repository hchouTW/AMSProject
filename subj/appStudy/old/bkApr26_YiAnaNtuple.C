#ifndef __YiAnaNtuple_C__
#define __YiAnaNtuple_C__

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.h"
#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.tcc"

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/src/minDST.h"

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
		
			//TTree * DSTRun = (runChain) ? runChain->CopyTree("") : 0;
			
			if (runChain != nullptr) {
				UInt_t trgEVSum = 0;
				UInt_t selEVSum = 0;
				TTree * DSTRun = new TTree("DSTRun", "DST Run Info");
				DSTRun->Branch("trgEV", &trgEVSum);
				DSTRun->Branch("selEV", &selEVSum);
				
				UInt_t trgEV = 0;
				UInt_t selEV = 0;
				runChain->SetBranchAddress("trgEV", &trgEV);
				runChain->SetBranchAddress("selEV", &selEV);
				for (Long64_t it = 0; it < runChain->GetEntries(); ++it) {
					runChain->GetEntry(it);
					trgEVSum += trgEV;
					selEVSum += selEV;
				}
				DSTRun->Fill();
			} 

			if (dataChain != 0 && gSaveDSTClone) fDST = dataChain->CloneTree(0);
			if (fDST == 0      && gSaveDSTTree)  {
				fDST = new TTree("data", "MinDST data");
				fDST->Branch("MCRig",  &fMCRig);
				fDST->Branch("TRRig",  &fTRRig);
				fDST->Branch("OptL" ,  &fOptL);
				fDST->Branch("OptIV" , &fOptI);
				fDST->Branch("OptIW" , &fOptM);
				fDST->Branch("OptHL1", &fOptHL1);
				fDST->Branch("OptHL9", &fOptHL9);
				fDST->Branch("OptHFs", &fOptHFs);
			}
			file->cd();
			
			Axis AXrso("Rigidity Resolution", 400, -1.5, 1.5);

			// rigidity binning
			//AXnr = Axis("|Rigidity| [GV]",
			//	{   1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
			//	    3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
			//		  8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
			//		 19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
			//		 41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
			//		 93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } );
			//AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
			//DST::WithPub = true;

			AXnr = Axis("|Rigidity| [GV]",
				{   1.00,   2.00,   3.00,   4.12,   5.00,   
				    6.00,   7.10,   8.30,   9.62,  11.04,  
					 12.59,  14.25,  16.05,  17.98,  20.04,  
					 22.25,  24.62,  27.25,  30.21  } );
			AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
			
			// 3 month
			AXtme = Axis("Time", 
				{ 1305417600, 1312416000, 1319414400, 1326412800, 1333411200,
				  1340409600, 1347408000, 1354406400, 1361404800, 1368403200,
					1375401600, 1382400000, 1389398400, 1396396800, 1403395200,
					1410393600, 1417392000, 1424390400, 1431388800, 1438387200,
					1445385600, 1452384000, 1459382400, 1466380800, 1473379200, 
					1480377600 } );

			//----  Low Energy  ----//
			Axis AXLaglPh("Ph", 50, 0, 50.);
			Hist::New("hLp_AglPh", "", AXnr, AXLaglPh);
			Hist::New("hLn_AglPh", "", AXnr, AXLaglPh);

			Axis AXLnafPh("Ph", 25, 0, 25.);
			Hist::New("hLp_NafPh", "", AXnr, AXLnafPh);
			Hist::New("hLn_NafPh", "", AXnr, AXLnafPh);

			Axis AXLaglM("Mass Estimator", 100, -6., 6.);
			Hist::New("hLp_AglM", "", AXnr, AXLaglM);
			Hist::New("hLn_AglM", "", AXnr, AXLaglM);
			
			Axis AXLnafM("Mass Estimator", 100, -8., 8.);
			Hist::New("hLp_NafM", "", AXnr, AXLnafM);
			Hist::New("hLn_NafM", "", AXnr, AXLnafM);
			
			Axis AXLlchi("lchi", 200, -6., 8.);
			Hist::New("hLp_Lchix", "", AXnr, AXLlchi);
			Hist::New("hLn_Lchix", "", AXnr, AXLlchi);
			Hist::New("hLp_Lchiy", "", AXnr, AXLlchi);
			Hist::New("hLn_Lchiy", "", AXnr, AXLlchi);

			Axis AXLlasym("Lasym", 200, -4., 5.);
			Hist::New("hLp_Lasym", "", AXnr, AXLlasym);
			Hist::New("hLn_Lasym", "", AXnr, AXLlasym);
			
			Axis AXLccest("CCest", 200, -5, 5);
			Hist::New("hLp_CCest", "", AXnr, AXLccest);
			Hist::New("hLn_CCest", "", AXnr, AXLccest);

			Axis AXLtrdEst("TRD Estimator", 200, 0.2, 1.6);
			Hist::New("hLp_TrdEst", "", AXnr, AXLtrdEst);
			Hist::New("hLn_TrdEst", "", AXnr, AXLtrdEst);
			
			Axis AXLnvtx("NVtx", 50, 0., 50.);
			Hist::New("hLp_TrdNVtx", "", AXnr, AXLnvtx);
			Hist::New("hLn_TrdNVtx", "", AXnr, AXLnvtx);
			
			Axis AXLagl("AglPh", 100, 0, 12.);
			Hist::New("hLp_AglPhEl", "", AXnr, AXLagl);
			Hist::New("hLn_AglPhEl", "", AXnr, AXLagl);
			Hist::New("hLp_AglPhPr", "", AXnr, AXLagl);
			Hist::New("hLn_AglPhPr", "", AXnr, AXLagl);
			
			Axis AXLnaf("NafPh", 100, 0, 12.);
			Hist::New("hLp_NafPhEl", "", AXnr, AXLnaf);
			Hist::New("hLn_NafPhEl", "", AXnr, AXLnaf);
			Hist::New("hLp_NafPhPr", "", AXnr, AXLnaf);
			Hist::New("hLn_NafPhPr", "", AXnr, AXLnaf);
		
			Axis AXLcutflow("Cutflow", 4, 0., 4.);
			Axis AXLtofM("Mass Estimator", 100, -4.5, 4.5);
	
			const Int_t N_TRDL = 1;
			const Int_t N_VETOL = 8;
			std::vector< std::vector<Int_t> > SetL;
			for (Int_t iTRD = 1; iTRD <= N_TRDL; ++iTRD) {
			for (Int_t iVETO = 1; iVETO <= N_VETOL; ++iVETO) {
				std::vector<Int_t> iSet( { iTRD, iVETO } );
				SetL.push_back(iSet);
			}}

			for (Int_t iSet = 0; iSet <= SetL.size(); ++iSet) {
				Int_t iTRD  = (iSet==SetL.size()) ? 0 : SetL.at(iSet).at(0);
				Int_t iVETO = (iSet==SetL.size()) ? 0 : SetL.at(iSet).at(1);
				std::string nmSet = (iSet==SetL.size()) ? "" : StrFmt("%d%d", iTRD, iVETO);

				hLpCutflow.push_back(Hist::New(StrFmt("hL%sp_Cutflow", nmSet.c_str()), "", AXnr, AXLcutflow));
				hLnCutflow.push_back(Hist::New(StrFmt("hL%sn_Cutflow", nmSet.c_str()), "", AXnr, AXLcutflow));
			
				hLsTofM.push_back(Hist::New(StrFmt("hL%ss_TofM", nmSet.c_str()), "", AXnr, AXLtofM));
				hLbTofM.push_back(Hist::New(StrFmt("hL%sb_TofM", nmSet.c_str()), "", AXnr, AXLtofM));
				
				hLpTofM.push_back(Hist::New(StrFmt("hL%sp_TofM", nmSet.c_str()), "", AXnr, AXLtofM));
				hLnTofM.push_back(Hist::New(StrFmt("hL%sn_TofM", nmSet.c_str()), "", AXnr, AXLtofM));

				hLpEvt.push_back(Hist::New(StrFmt("hL%sp_Evt", nmSet.c_str()), "", AXnr));
				hLnEvt.push_back(Hist::New(StrFmt("hL%sn_Evt", nmSet.c_str()), "", AXnr));
			
				if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
					hLpMCEvt.push_back(Hist::New(StrFmt("hL%sp_MCEvt", nmSet.c_str()), "", AXnr));
					hLnMCEvt.push_back(Hist::New(StrFmt("hL%sn_MCEvt", nmSet.c_str()), "", AXnr));
				}
			
				if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
					hTLpCutflow.push_back(Hist::New(StrFmt("hTL%sp_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXLcutflow));
					hTLnCutflow.push_back(Hist::New(StrFmt("hTL%sn_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXLcutflow));
					
					hTLsTofM.push_back(Hist::New(StrFmt("hTL%ss_TofM", nmSet.c_str()), "", AXtme, AXnr, AXLtofM));
					hTLbTofM.push_back(Hist::New(StrFmt("hTL%sb_TofM", nmSet.c_str()), "", AXtme, AXnr, AXLtofM));
          
					hTLpTofM.push_back(Hist::New(StrFmt("hTL%sp_TofM", nmSet.c_str()), "", AXtme, AXnr, AXLtofM));
					hTLnTofM.push_back(Hist::New(StrFmt("hTL%sn_TofM", nmSet.c_str()), "", AXtme, AXnr, AXLtofM));
					
					hTLpEvt.push_back(Hist::New(StrFmt("hTL%sp_Evt", nmSet.c_str()), "", AXtme, AXnr));
					hTLnEvt.push_back(Hist::New(StrFmt("hTL%sn_Evt", nmSet.c_str()), "", AXtme, AXnr));
				}
			}


			//----  Intermedia Energy  ----//
			Axis AXIlchi("lchi", 200, -6., 8.);
			Hist::New("hIp_Lchix", "", AXnr, AXIlchi);
			Hist::New("hIn_Lchix", "", AXnr, AXIlchi);
			Hist::New("hIp_Lchiy", "", AXnr, AXIlchi);
			Hist::New("hIn_Lchiy", "", AXnr, AXIlchi);

			Axis AXIlasym("Lasym", 200, -4., 5.);
			Hist::New("hIp_Lasym", "", AXnr, AXIlasym);
			Hist::New("hIn_Lasym", "", AXnr, AXIlasym);
			
			Axis AXIccest("CCest", 200, -5, 5);
			Hist::New("hIp_CCest", "", AXnr, AXIccest);
			Hist::New("hIn_CCest", "", AXnr, AXIccest);
				
			Axis AXIcutflow("Cutflow", 3, 0., 3.);
			Axis AXItrdEst("TRD Estimator", 100, 0.2, 1.6);
			
			const Int_t N_MASSI = 1;
			std::vector< std::vector<Int_t> > SetI;
			for (Int_t iMASS = 1; iMASS <= N_MASSI; ++iMASS) {
				std::vector<Int_t> iSet( { iMASS } );
				SetI.push_back(iSet);
			}

			for (Int_t iSet = 0; iSet <= SetI.size(); ++iSet) {
				Int_t iMASS = (iSet==SetI.size()) ? 0 : SetI.at(iSet).at(0);
				std::string nmSet = (iSet==SetI.size()) ? "" : StrFmt("%d", iMASS);

				hIpCutflow.push_back(Hist::New(StrFmt("hI%sp_Cutflow", nmSet.c_str()), "", AXnr, AXIcutflow));
				hInCutflow.push_back(Hist::New(StrFmt("hI%sn_Cutflow", nmSet.c_str()), "", AXnr, AXIcutflow));

				hIsTrdL.push_back(Hist::New(StrFmt("hI%ss_TrdEst", nmSet.c_str()), "", AXnr, AXItrdEst));
				hIbTrdL.push_back(Hist::New(StrFmt("hI%sb_TrdEst", nmSet.c_str()), "", AXnr, AXItrdEst));
				
				hIpTrdL.push_back(Hist::New(StrFmt("hI%sp_TrdEst", nmSet.c_str()), "", AXnr, AXItrdEst));
				hInTrdL.push_back(Hist::New(StrFmt("hI%sn_TrdEst", nmSet.c_str()), "", AXnr, AXItrdEst));
				
				hIpEvt.push_back(Hist::New(StrFmt("hI%sp_Evt", nmSet.c_str()), "", AXnr));
				hInEvt.push_back(Hist::New(StrFmt("hI%sn_Evt", nmSet.c_str()), "", AXnr));
				
				if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
					hIpMCEvt.push_back(Hist::New(StrFmt("hI%sp_MCEvt", nmSet.c_str()), "", AXnr));
					hInMCEvt.push_back(Hist::New(StrFmt("hI%sn_MCEvt", nmSet.c_str()), "", AXnr));
				}
				
				if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
					hTIpCutflow.push_back(Hist::New(StrFmt("hTI%sp_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXIcutflow));
					hTInCutflow.push_back(Hist::New(StrFmt("hTI%sn_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXIcutflow));

					hTIsTrdL.push_back(Hist::New(StrFmt("hTI%ss_TrdEst", nmSet.c_str()), "", AXtme, AXnr, AXItrdEst));
					hTIbTrdL.push_back(Hist::New(StrFmt("hTI%sb_TrdEst", nmSet.c_str()), "", AXtme, AXnr, AXItrdEst));
					
					hTIpTrdL.push_back(Hist::New(StrFmt("hTI%sp_TrdEst", nmSet.c_str()), "", AXtme, AXnr, AXItrdEst));
					hTInTrdL.push_back(Hist::New(StrFmt("hTI%sn_TrdEst", nmSet.c_str()), "", AXtme, AXnr, AXItrdEst));
					
					hTIpEvt.push_back(Hist::New(StrFmt("hTI%sp_Evt", nmSet.c_str()), "", AXtme, AXnr));
					hTInEvt.push_back(Hist::New(StrFmt("hTI%sn_Evt", nmSet.c_str()), "", AXtme, AXnr));
				}
			}	
			
			Axis AXMlchi("lchi", 200, -6., 8.);
			Hist::New("hMp_Lchix", "", AXnr, AXMlchi);
			Hist::New("hMn_Lchix", "", AXnr, AXMlchi);
			Hist::New("hMp_Lchiy", "", AXnr, AXMlchi);
			Hist::New("hMn_Lchiy", "", AXnr, AXMlchi);

			Axis AXMlasym("Lasym", 200, -4., 5.);
			Hist::New("hMp_Lasym", "", AXnr, AXMlasym);
			Hist::New("hMn_Lasym", "", AXnr, AXMlasym);
			
			Axis AXMccest("CCest", 200, -5, 5);
			Hist::New("hMp_CCest", "", AXnr, AXMccest);
			Hist::New("hMn_CCest", "", AXnr, AXMccest);
				
			Axis AXMcutflow("Cutflow", 3, 0., 3.);
			Axis AXMtrdEst("TRD Estimator", 100, 0.2, 1.6);
			
			const Int_t N_MASSM = 1;
			std::vector< std::vector<Int_t> > SetM;
			for (Int_t iMASS = 1; iMASS <= N_MASSM; ++iMASS) {
				std::vector<Int_t> iSet( { iMASS } );
				SetM.push_back(iSet);
			}

			for (Int_t iSet = 0; iSet <= SetM.size(); ++iSet) {
				Int_t iMASS = (iSet==SetM.size()) ? 0 : SetM.at(iSet).at(0);
				std::string nmSet = (iSet==SetM.size()) ? "" : StrFmt("%d", iMASS);

				hMpCutflow.push_back(Hist::New(StrFmt("hM%sp_Cutflow", nmSet.c_str()), "", AXnr, AXMcutflow));
				hMnCutflow.push_back(Hist::New(StrFmt("hM%sn_Cutflow", nmSet.c_str()), "", AXnr, AXMcutflow));

				hMsTrdL.push_back(Hist::New(StrFmt("hM%ss_TrdEst", nmSet.c_str()), "", AXnr, AXMtrdEst));
				hMbTrdL.push_back(Hist::New(StrFmt("hM%sb_TrdEst", nmSet.c_str()), "", AXnr, AXMtrdEst));
				
				hMpTrdL.push_back(Hist::New(StrFmt("hM%sp_TrdEst", nmSet.c_str()), "", AXnr, AXMtrdEst));
				hMnTrdL.push_back(Hist::New(StrFmt("hM%sn_TrdEst", nmSet.c_str()), "", AXnr, AXMtrdEst));
				
				hMpEvt.push_back(Hist::New(StrFmt("hM%sp_Evt", nmSet.c_str()), "", AXnr));
				hMnEvt.push_back(Hist::New(StrFmt("hM%sn_Evt", nmSet.c_str()), "", AXnr));
				
				if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
					hMpMCEvt.push_back(Hist::New(StrFmt("hM%sp_MCEvt", nmSet.c_str()), "", AXnr));
					hMnMCEvt.push_back(Hist::New(StrFmt("hM%sn_MCEvt", nmSet.c_str()), "", AXnr));
				}
				
				if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
					hTMpCutflow.push_back(Hist::New(StrFmt("hTM%sp_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXMcutflow));
					hTMnCutflow.push_back(Hist::New(StrFmt("hTM%sn_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXMcutflow));

					hTMsTrdL.push_back(Hist::New(StrFmt("hTM%ss_TrdEst", nmSet.c_str()), "", AXtme, AXnr, AXMtrdEst));
					hTMbTrdL.push_back(Hist::New(StrFmt("hTM%sb_TrdEst", nmSet.c_str()), "", AXtme, AXnr, AXMtrdEst));
					
					hTMpTrdL.push_back(Hist::New(StrFmt("hTM%sp_TrdEst", nmSet.c_str()), "", AXtme, AXnr, AXMtrdEst));
					hTMnTrdL.push_back(Hist::New(StrFmt("hTM%sn_TrdEst", nmSet.c_str()), "", AXtme, AXnr, AXMtrdEst));
					
					hTMpEvt.push_back(Hist::New(StrFmt("hTM%sp_Evt", nmSet.c_str()), "", AXtme, AXnr));
					hTMnEvt.push_back(Hist::New(StrFmt("hTM%sn_Evt", nmSet.c_str()), "", AXtme, AXnr));
				}
			}	

/*
			//----  High Energy  ----//
			Axis AXHL1cutflow("Cutflow", 3, 0., 3.);
			Axis AXHL1cc("CC Estimator", 100, 0.2, 1.6);
			
			const Int_t N_CCHL1 = 2;
			std::vector< std::vector<Int_t> > SetHL1;
			for (Int_t iCC = 1; iCC <= N_CCHL1; ++iCC) {
				std::vector<Int_t> iSet( { iCC } );
				SetHL1.push_back(iSet);
			}

			for (Int_t iSet = 0; iSet <= SetHL1.size(); ++iSet) {
				Int_t iTRD  = (iSet==SetHL1.size()) ? 0 : SetHL1.at(iSet).at(0);
				std::string nmSet = (iSet==SetHL1.size()) ? "" : StrFmt("%d", iTRD);

				hHL1pCutflow.push_back(Hist::New(StrFmt("hHL1%sp_Cutflow", nmSet.c_str()), "", AXnr, AXHL1cutflow));
				hHL1nCutflow.push_back(Hist::New(StrFmt("hHL1%sn_Cutflow", nmSet.c_str()), "", AXnr, AXHL1cutflow));

				hHL1sTrdL.push_back(Hist::New(StrFmt("hHL1%ss_CC", nmSet.c_str()), "", AXnr, AXHL1cc));
				hHL1bTrdL.push_back(Hist::New(StrFmt("hHL1%sb_CC", nmSet.c_str()), "", AXnr, AXHL1cc));
				
				hHL1pTrdL.push_back(Hist::New(StrFmt("hHL1%sp_CC", nmSet.c_str()), "", AXnr, AXHL1cc));
				hHL1nTrdL.push_back(Hist::New(StrFmt("hHL1%sn_CC", nmSet.c_str()), "", AXnr, AXHL1cc));
				
				hHL1pEvt.push_back(Hist::New(StrFmt("hHL1%sp_Evt", nmSet.c_str()), "", AXnr));
				hHL1nEvt.push_back(Hist::New(StrFmt("hHL1%sn_Evt", nmSet.c_str()), "", AXnr));
				
				if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
					hHL1pMCEvt.push_back(Hist::New(StrFmt("hHL1%sp_MCEvt", nmSet.c_str()), "", AXnr));
					hHL1nMCEvt.push_back(Hist::New(StrFmt("hHL1%sn_MCEvt", nmSet.c_str()), "", AXnr));
				}
				
				if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
					hTHL1pCutflow.push_back(Hist::New(StrFmt("hTHL1%sp_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXHL1cutflow));
					hTHL1nCutflow.push_back(Hist::New(StrFmt("hTHL1%sn_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXHL1cutflow));

					hTHL1sTrdL.push_back(Hist::New(StrFmt("hTHL1%ss_CC", nmSet.c_str()), "", AXtme, AXnr, AXHL1cc));
					hTHL1bTrdL.push_back(Hist::New(StrFmt("hTHL1%sb_CC", nmSet.c_str()), "", AXtme, AXnr, AXHL1cc));
					
					hTHL1pTrdL.push_back(Hist::New(StrFmt("hTHL1%sp_CC", nmSet.c_str()), "", AXtme, AXnr, AXHL1cc));
					hTHL1nTrdL.push_back(Hist::New(StrFmt("hTHL1%sn_CC", nmSet.c_str()), "", AXtme, AXnr, AXHL1cc));
					
					hTHL1pEvt.push_back(Hist::New(StrFmt("hTHL1%sp_Evt", nmSet.c_str()), "", AXtme, AXnr));
					hTHL1nEvt.push_back(Hist::New(StrFmt("hTHL1%sn_Evt", nmSet.c_str()), "", AXtme, AXnr));
				}
			}	













*/




























			//----  High Energy  ----//
			/*
			Axis AXHcutflow("Cutflow", 4, 0., 4.);
			Hist::New("hHp_CutflowL1", "", AXnr, AXHcutflow);
			Hist::New("hHn_CutflowL1", "", AXnr, AXHcutflow);
			Hist::New("hHp_CutflowL9", "", AXnr, AXHcutflow);
			Hist::New("hHn_CutflowL9", "", AXnr, AXHcutflow);
			Hist::New("hHp_CutflowFs", "", AXnr, AXHcutflow);
			Hist::New("hHn_CutflowFs", "", AXnr, AXHcutflow);
			
			Axis AXHccest("CCest", 50, -4, 4);
			Hist::New("hHp_CCestL1", "", AXnr, AXHccest);
			Hist::New("hHn_CCestL1", "", AXnr, AXHccest);
			Hist::New("hHp_CCestL9", "", AXnr, AXHccest);
			Hist::New("hHn_CCestL9", "", AXnr, AXHccest);
			Hist::New("hHp_CCestFs", "", AXnr, AXHccest);
			Hist::New("hHn_CCestFs", "", AXnr, AXHccest);
			
			Hist::New("hHp_EvtL1", "", AXnr);
			Hist::New("hHn_EvtL1", "", AXnr);
			Hist::New("hHp_EvtL9", "", AXnr);
			Hist::New("hHn_EvtL9", "", AXnr);
			Hist::New("hHp_EvtFs", "", AXnr);
			Hist::New("hHn_EvtFs", "", AXnr);
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hHp_MCEvtL1", "", AXnr);
				Hist::New("hHn_MCEvtL1", "", AXnr);
				
				Hist::New("hHp_MCEvtL9", "", AXnr);
				Hist::New("hHn_MCEvtL9", "", AXnr);
				
				Hist::New("hHp_MCEvtFs", "", AXnr);
				Hist::New("hHn_MCEvtFs", "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Hist::New("hTHp_CutflowL1", "", AXtme, AXnr, AXHcutflow);
				Hist::New("hTHn_CutflowL1", "", AXtme, AXnr, AXHcutflow);
				Hist::New("hTHp_CutflowL9", "", AXtme, AXnr, AXHcutflow);
				Hist::New("hTHn_CutflowL9", "", AXtme, AXnr, AXHcutflow);
				Hist::New("hTHp_CutflowFs", "", AXtme, AXnr, AXHcutflow);
				Hist::New("hTHn_CutflowFs", "", AXtme, AXnr, AXHcutflow);
				
				Axis AXTHccest("CCest", 50, -4, 4);
				Hist::New("hTHp_CCestL1", "", AXtme, AXnr, AXTHccest);
				Hist::New("hTHn_CCestL1", "", AXtme, AXnr, AXTHccest);
				Hist::New("hTHp_CCestL9", "", AXtme, AXnr, AXTHccest);
				Hist::New("hTHn_CCestL9", "", AXtme, AXnr, AXTHccest);
				Hist::New("hTHp_CCestFs", "", AXtme, AXnr, AXTHccest);
				Hist::New("hTHn_CCestFs", "", AXtme, AXnr, AXTHccest);
				
				Hist::New("hTHp_EvtL1", "", AXtme, AXnr);
				Hist::New("hTHn_EvtL1", "", AXtme, AXnr);
				Hist::New("hTHp_EvtL9", "", AXtme, AXnr);
				Hist::New("hTHn_EvtL9", "", AXtme, AXnr);
				Hist::New("hTHp_EvtFs", "", AXtme, AXnr);
				Hist::New("hTHn_EvtFs", "", AXtme, AXnr);
			}
			*/

			std::cout << "\n<<  init DST End  >>\n";
			file->cd();
			return;
		}

		void resetDST() {
#if Debug == true
	std::cerr << "DST::resetDST()\n";
#endif
			fMCRig = 0;
			fTRRig = 0;
			fOptL   = false;
			fOptI   = false;
			fOptM   = false;
			fOptHL1 = false;
			fOptHL9 = false;
			fOptHFs = false;
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
		TTree *  fDST;
	
	public :
		Float_t fMCRig;
		Float_t fTRRig;
		Bool_t  fOptL;
		Bool_t  fOptI;
		Bool_t  fOptM;
		Bool_t  fOptHL1;
		Bool_t  fOptHL9;
		Bool_t  fOptHFs;

	public :
		static Bool_t gSaveDST;
		static Bool_t gSaveDSTTree;
		static Bool_t gSaveDSTClone;
		static UInt_t gUTime[2]; // (pre, cur)

	public :
		static Axis AXnr;
		static Axis AXir;
		static Axis AXtme;

	public :
		static Bool_t WithPub;

	public :
		static std::vector<Hist *> hLpCutflow;
		static std::vector<Hist *> hLnCutflow;
		static std::vector<Hist *> hLsTofM;
		static std::vector<Hist *> hLbTofM;
		static std::vector<Hist *> hLpTofM;
		static std::vector<Hist *> hLnTofM;
		static std::vector<Hist *> hLpEvt;
		static std::vector<Hist *> hLnEvt;
		static std::vector<Hist *> hLpMCEvt;
		static std::vector<Hist *> hLnMCEvt;
		
		static std::vector<Hist *> hTLpCutflow;
		static std::vector<Hist *> hTLnCutflow;
		static std::vector<Hist *> hTLsTofM;
		static std::vector<Hist *> hTLbTofM;
		static std::vector<Hist *> hTLpTofM;
		static std::vector<Hist *> hTLnTofM;
		static std::vector<Hist *> hTLpEvt;
		static std::vector<Hist *> hTLnEvt;
		
		static std::vector<Hist *> hIpCutflow;
		static std::vector<Hist *> hInCutflow;
		static std::vector<Hist *> hIsTrdL;
		static std::vector<Hist *> hIbTrdL;
		static std::vector<Hist *> hIpTrdL;
		static std::vector<Hist *> hInTrdL;
		static std::vector<Hist *> hIpEvt;
		static std::vector<Hist *> hInEvt;
		static std::vector<Hist *> hIpMCEvt;
		static std::vector<Hist *> hInMCEvt;
		
		static std::vector<Hist *> hTIpCutflow;
		static std::vector<Hist *> hTInCutflow;
		static std::vector<Hist *> hTIsTrdL;
		static std::vector<Hist *> hTIbTrdL;
		static std::vector<Hist *> hTIpTrdL;
		static std::vector<Hist *> hTInTrdL;
		static std::vector<Hist *> hTIpEvt;
		static std::vector<Hist *> hTInEvt;
		
		static std::vector<Hist *> hMpCutflow;
		static std::vector<Hist *> hMnCutflow;
		static std::vector<Hist *> hMsTrdL;
		static std::vector<Hist *> hMbTrdL;
		static std::vector<Hist *> hMpTrdL;
		static std::vector<Hist *> hMnTrdL;
		static std::vector<Hist *> hMpEvt;
		static std::vector<Hist *> hMnEvt;
		static std::vector<Hist *> hMpMCEvt;
		static std::vector<Hist *> hMnMCEvt;
		
		static std::vector<Hist *> hTMpCutflow;
		static std::vector<Hist *> hTMnCutflow;
		static std::vector<Hist *> hTMsTrdL;
		static std::vector<Hist *> hTMbTrdL;
		static std::vector<Hist *> hTMpTrdL;
		static std::vector<Hist *> hTMnTrdL;
		static std::vector<Hist *> hTMpEvt;
		static std::vector<Hist *> hTMnEvt;
		
		static std::vector<Hist *> hHL1pCutflow;
		static std::vector<Hist *> hHL1nCutflow;
		static std::vector<Hist *> hHL1sTrdL;
		static std::vector<Hist *> hHL1bTrdL;
		static std::vector<Hist *> hHL1pTrdL;
		static std::vector<Hist *> hHL1nTrdL;
		static std::vector<Hist *> hHL1pEvt;
		static std::vector<Hist *> hHL1nEvt;
		static std::vector<Hist *> hHL1pMCEvt;
		static std::vector<Hist *> hHL1nMCEvt;
		
		static std::vector<Hist *> hTHL1pCutflow;
		static std::vector<Hist *> hTHL1nCutflow;
		static std::vector<Hist *> hTHL1sTrdL;
		static std::vector<Hist *> hTHL1bTrdL;
		static std::vector<Hist *> hTHL1pTrdL;
		static std::vector<Hist *> hTHL1nTrdL;
		static std::vector<Hist *> hTHL1pEvt;
		static std::vector<Hist *> hTHL1nEvt;
};

Bool_t DST::gSaveDST = true;
Bool_t DST::gSaveDSTTree = false;
Bool_t DST::gSaveDSTClone = false;
UInt_t DST::gUTime[2] = {0, 0};

Axis   DST::AXnr;
Axis   DST::AXir;
Axis   DST::AXtme;
Bool_t DST::WithPub = false;
		
std::vector<Hist *> DST::hLpCutflow;
std::vector<Hist *> DST::hLnCutflow;
std::vector<Hist *> DST::hLsTofM;
std::vector<Hist *> DST::hLbTofM;
std::vector<Hist *> DST::hLpTofM;
std::vector<Hist *> DST::hLnTofM;
std::vector<Hist *> DST::hLpEvt;
std::vector<Hist *> DST::hLnEvt;
std::vector<Hist *> DST::hLpMCEvt;
std::vector<Hist *> DST::hLnMCEvt;

std::vector<Hist *> DST::hTLpCutflow;
std::vector<Hist *> DST::hTLnCutflow;
std::vector<Hist *> DST::hTLsTofM;
std::vector<Hist *> DST::hTLbTofM;
std::vector<Hist *> DST::hTLpTofM;
std::vector<Hist *> DST::hTLnTofM;
std::vector<Hist *> DST::hTLpEvt;
std::vector<Hist *> DST::hTLnEvt;

std::vector<Hist *> DST::hIpCutflow;
std::vector<Hist *> DST::hInCutflow;
std::vector<Hist *> DST::hIsTrdL;
std::vector<Hist *> DST::hIbTrdL;
std::vector<Hist *> DST::hIpTrdL;
std::vector<Hist *> DST::hInTrdL;
std::vector<Hist *> DST::hIpEvt;
std::vector<Hist *> DST::hInEvt;
std::vector<Hist *> DST::hIpMCEvt;
std::vector<Hist *> DST::hInMCEvt;

std::vector<Hist *> DST::hTIpCutflow;
std::vector<Hist *> DST::hTInCutflow;
std::vector<Hist *> DST::hTIsTrdL;
std::vector<Hist *> DST::hTIbTrdL;
std::vector<Hist *> DST::hTIpTrdL;
std::vector<Hist *> DST::hTInTrdL;
std::vector<Hist *> DST::hTIpEvt;
std::vector<Hist *> DST::hTInEvt;

std::vector<Hist *> DST::hMpCutflow;
std::vector<Hist *> DST::hMnCutflow;
std::vector<Hist *> DST::hMsTrdL;
std::vector<Hist *> DST::hMbTrdL;
std::vector<Hist *> DST::hMpTrdL;
std::vector<Hist *> DST::hMnTrdL;
std::vector<Hist *> DST::hMpEvt;
std::vector<Hist *> DST::hMnEvt;
std::vector<Hist *> DST::hMpMCEvt;
std::vector<Hist *> DST::hMnMCEvt;

std::vector<Hist *> DST::hTMpCutflow;
std::vector<Hist *> DST::hTMnCutflow;
std::vector<Hist *> DST::hTMsTrdL;
std::vector<Hist *> DST::hTMbTrdL;
std::vector<Hist *> DST::hTMpTrdL;
std::vector<Hist *> DST::hTMnTrdL;
std::vector<Hist *> DST::hTMpEvt;
std::vector<Hist *> DST::hTMnEvt;

std::vector<Hist *> DST::hHL1pCutflow;
std::vector<Hist *> DST::hHL1nCutflow;
std::vector<Hist *> DST::hHL1sTrdL;
std::vector<Hist *> DST::hHL1bTrdL;
std::vector<Hist *> DST::hHL1pTrdL;
std::vector<Hist *> DST::hHL1nTrdL;
std::vector<Hist *> DST::hHL1pEvt;
std::vector<Hist *> DST::hHL1nEvt;
std::vector<Hist *> DST::hHL1pMCEvt;
std::vector<Hist *> DST::hHL1nMCEvt;

std::vector<Hist *> DST::hTHL1pCutflow;
std::vector<Hist *> DST::hTHL1nCutflow;
std::vector<Hist *> DST::hTHL1sTrdL;
std::vector<Hist *> DST::hTHL1bTrdL;
std::vector<Hist *> DST::hTHL1pTrdL;
std::vector<Hist *> DST::hTHL1nTrdL;
std::vector<Hist *> DST::hTHL1pEvt;
std::vector<Hist *> DST::hTHL1nEvt;



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
		MinDST * fMDst;
};

YiAna::YiAna() {
	DST::fDST = 0;
	YiNtuple::fRunChain = 0;
	YiNtuple::fDataChain = 0;
	YiNtuple::init();
	YiNtuple::fTimer.start();
	fMDst = 0;
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
	fMDst = new MinDST;
	
	SetBranchAddressMinDST(fDataChain, *fMDst);
	/*******************/
	/**               **/
	/*******************/

	std::cout << "\n<<  Set Branch Address End  >>\n";
}

void YiAna::freeBranchAddress() {
	if (fMDst) { delete fMDst; fMDst = 0; }
}

void YiAna::analyzeEvent() {
	std::cout << "\n<<  Analyze Event Info  >>\n";

	/*******************/
	/**  USER DEFINE  **/
	/*******************/
	//if (fileList.at(0).find("MC_Pr_PL1_0_510") != std::string::npos) DST::gSaveDSTTree = true;
	//if (fileList.at(0).find("MC_Pr_PL1_1800") != std::string::npos) DST::gSaveDSTTree = true;
	//if (fileList.at(0).find("MC_Ap_PL1_0_510") != std::string::npos) DST::gSaveDSTTree = true;
	//if (fileList.at(0).find("MC_Ap_PL1_1800") != std::string::npos) DST::gSaveDSTTree = true;
	//if (fileList.at(0).find("Pr0_510") != std::string::npos) DST::gSaveDSTTree = true;
	//if (fileList.at(0).find("Pr1800") != std::string::npos) DST::gSaveDSTTree = true;
	//if (fileList.at(0).find("Ap0_510") != std::string::npos) DST::gSaveDSTTree = true;
	//if (fileList.at(0).find("Ap1800") != std::string::npos) DST::gSaveDSTTree = true;
	
	initDST(fFile, fRunChain, fDataChain);
	fFile->cd();
	/*******************/
	/**               **/
	/*******************/

	Long64_t nentries = fDataChain->GetEntries();
	Long64_t npassed = 0;
	Long64_t nprocessed = 0;
	Long64_t printRate = nentries / 100;
	if (printRate < 50000) printRate = 50000;
	if (printRate > 500000) printRate = 500000;

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
		if (!analyzeRunInfo()) continue;
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
	
	// Event
	Float_t weight = fMDst->weight;
	
	// MC
	Bool_t  IsMonteCarlo = (YiNtuple::CheckEventMode(YiNtuple::MC) && (fMDst->mcSign != 0));
	Short_t MCSign       = fMDst->mcSign;
	Float_t MCNRig       = fMDst->mcNRig;
	
	// RTI
	Bool_t  IsISS = (YiNtuple::CheckEventMode(YiNtuple::ISS) && (fMDst->mcSign == 0));
	UInt_t  uTime = fMDst->uTime;

	// <<Pub>> 4-years < May 19, 2011 ~ May 26, 2015 >
	UInt_t period[2] = { 1305763200, 1432684800 };
	if (DST::WithPub && IsISS && (uTime < period[0] || uTime >= period[1])) return false;

	// TRD
	Float_t trdEst  = fMDst->trdEst;

	// ECAL
	Bool_t isEcalPr    = (!fMDst->hasEcal) ? false : (fMDst->ecalEst < -0.8); 
	Bool_t isEcalEl    = (!fMDst->hasEcal) ? false : (fMDst->ecalEst >  0.6); 
	Bool_t isNotEcalPr = (!fMDst->hasEcal) ? false : (fMDst->ecalEst > -0.8); 

	// IsProtonLike (TRD ECAL RICH)
	Bool_t isPrLike = ( (fMDst->trdEst > 0.75) && 
	                    (fMDst->hasEcal ? (fMDst->ecalEst < -0.8) : true) && 
										  (fMDst->hasRich ? (fMDst->richMPr > -2.0) : true) );

	// Charge Estimator
	Float_t lchiyTh[4] = { 2.90, 2.20, 2.30, 1.90 }; // ~98%
	Float_t lasymTh[4] = { 1.00, 1.60, 1.40, 1.40 }; // ~98%
	Float_t lchiyIn = (fMDst->lchiy[0] / lchiyTh[0]);
	Float_t lchiyL1 = (fMDst->lchiy[1] / lchiyTh[1]);
	Float_t lchiyL9 = (fMDst->lchiy[2] / lchiyTh[2]);
	Float_t lchiyFs = (fMDst->lchiy[3] / lchiyTh[3]);
	Float_t lasymUL = ((2.0 * fMDst->lasym[0] - fMDst->lchiy[0]) / lasymTh[0]);
	Float_t lasym1I = ((2.0 * fMDst->lasym[1] - fMDst->lchiy[1]) / lasymTh[1]);
	Float_t lasym9I = ((2.0 * fMDst->lasym[2] - fMDst->lchiy[2]) / lasymTh[2]);
	Float_t lasym91 = ((2.0 * fMDst->lasym[3] - fMDst->lchiy[3]) / lasymTh[3]);

	Float_t ccestSf[4] = { 0.20, 0.27, 0.21, 0.27 };
	Float_t ccestSg[4] = { 0.40, 0.46, 0.47, 0.50 };
	Float_t ccestIn = ((((lchiyIn>0&&lasymUL>0) ? std::sqrt(lchiyIn*lchiyIn+lasymUL*lasymUL) : std::max(lchiyIn, lasymUL)) - ccestSf[0]) / ccestSg[0]);
	Float_t ccestL1 = ((((lchiyL1>0&&lasym1I>0) ? std::sqrt(lchiyL1*lchiyL1+lasym1I*lasym1I) : std::max(lchiyL1, lasym1I)) - ccestSf[1]) / ccestSg[1]);
	Float_t ccestL9 = ((((lchiyL9>0&&lasym9I>0) ? std::sqrt(lchiyL9*lchiyL9+lasym9I*lasym9I) : std::max(lchiyL9, lasym9I)) - ccestSf[2]) / ccestSg[2]);
	Float_t ccestFs = ((((lchiyFs>0&&lasym91>0) ? std::sqrt(lchiyFs*lchiyFs+lasym91*lasym91) : std::max(lchiyFs, lasym91)) - ccestSf[3]) / ccestSg[3]);
	
	// Charge-Confusion Estimator
	//Float_t ccestTh[4] = { (std::erfc(std::log(50.0)-std::log(1+std::fabs(fMDst->trTRigIn)))+2.5), 2.5, 3.2, 2.7 }; // ~98%
	Float_t ccestTh[4] = { 4.2, 2.5, 3.2, 2.7 }; // ~98%
	Float_t ccestPt[4] = { ccestIn, ccestL1, ccestL9, ccestFs };
	Float_t tuneCCPnt[4] = { 20.0, 25.0, 30.0, 40.0 };

	/**** Low Energy Region ****/
	while (true) {
		Short_t     trPt  = 2;
		std::string trNm  = "In";
		Short_t     sign  = fMDst->sign[trPt-2];
		Float_t     nrig  = fMDst->nrig[trPt-2];
		Float_t     lchix = fMDst->lchix[trPt-2];
		Float_t     lchiy = fMDst->lchiy[trPt-2];
		Float_t     lasym = fMDst->lasym[trPt-2];
		Float_t     ccest = ccestPt[trPt-2];
		if (sign == 0) break;
		
		// RTI
		if ((nrig / fMDst->cfRig) < 1.2) break;
			
		// ECAL
		if (isNotEcalPr) break;
	
		// RICH
		if (fMDst->richRad == -1) break;
			
		// RICH NO-RING
		Bool_t isAglNORing = 
			(fMDst->richRad == 0 && !fMDst->hasRich && 
			 fMDst->richPhEl > 2.0 && fMDst->richPhPr < 2.0);
		if (isAglNORing) {
			if (sign>0) Hist::Head("hLp_AglPh")->fill(nrig, fMDst->richNHit, weight);
			if (sign<0) Hist::Head("hLn_AglPh")->fill(nrig, fMDst->richNHit, weight);
		}


		Bool_t isNafNORing = 
			(fMDst->richRad == 1 && !fMDst->hasRich && 
			 fMDst->richPhEl > 0.5 && fMDst->richPhPr < 0.5);
		if (isNafNORing) {
			if (sign>0) Hist::Head("hLp_NafPh")->fill(nrig, fMDst->richNHit, weight);
			if (sign<0) Hist::Head("hLn_NafPh")->fill(nrig, fMDst->richNHit, weight);
		}
	
		//if (isAglNORing && fMDst->richNHit >= 12) break;  // 5 -> 15
		//if (isNafNORing && fMDst->richNHit >=  5) break;  // 4 -> 6
		
		// RICH RING
		Bool_t isAglRing = (fMDst->richRad == 0 && fMDst->hasRich);
		if (isAglRing) {
			if (sign>0) Hist::Head("hLp_AglM")->fill(nrig, fMDst->richMPr, weight);
			if (sign<0) Hist::Head("hLn_AglM")->fill(nrig, fMDst->richMPi, weight);
		}
		Bool_t isNafRing = (fMDst->richRad == 1 && fMDst->hasRich);
		if (isNafRing) {
			if (sign>0) Hist::Head("hLp_NafM")->fill(nrig, fMDst->richMPr, weight);
			if (sign<0) Hist::Head("hLn_NafM")->fill(nrig, fMDst->richMPi, weight);
		}

		// Charge Confusion
		if (sign>0) Hist::Head("hLp_Lchix")->fill(nrig, lchix, weight);
		if (sign<0) Hist::Head("hLn_Lchix")->fill(nrig, lchix, weight);
		
		if (lchix > 2.7) break;
	
		if (sign>0) Hist::Head("hLp_Lchiy")->fill(nrig, lchiy, weight);
		if (sign<0) Hist::Head("hLn_Lchiy")->fill(nrig, lchiy, weight);
		
		if (lchiy > 2.4) break;
		
		if (sign>0) Hist::Head("hLp_Lasym")->fill(nrig, lasym, weight);
		if (sign<0) Hist::Head("hLn_Lasym")->fill(nrig, lasym, weight);
		
		if (lasym > 1.3) break;
		
		if (sign>0) Hist::Head("hLp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hLn_CCest")->fill(nrig, ccest, weight);
	
		// TRD
		if (sign>0) Hist::Head("hLp_TrdEst")->fill(nrig, trdEst, weight);
		if (sign<0) Hist::Head("hLn_TrdEst")->fill(nrig, trdEst, weight);
		
		
		if (sign>0) Hist::Head("hLp_TrdNVtx")->fill(nrig, fMDst->trdNVtx, weight);
		if (sign<0) Hist::Head("hLn_TrdNVtx")->fill(nrig, fMDst->trdNVtx, weight);
		
		if (fMDst->richRad == 0 && !fMDst->hasRich) {
			if (sign>0) Hist::Head("hLp_AglPhEl")->fill(nrig, fMDst->richPhEl, weight);
			if (sign<0) Hist::Head("hLn_AglPhEl")->fill(nrig, fMDst->richPhEl, weight);
			if (fMDst->richPhEl > 2.0) {
				if (sign>0) Hist::Head("hLp_AglPhPr")->fill(nrig, fMDst->richPhPr, weight);
				if (sign<0) Hist::Head("hLn_AglPhPr")->fill(nrig, fMDst->richPhPr, weight);
			}
		}
		
		if (fMDst->richRad == 1 && !fMDst->hasRich) {
			if (sign>0) Hist::Head("hLp_NafPhEl")->fill(nrig, fMDst->richPhEl, weight);
			if (sign<0) Hist::Head("hLn_NafPhEl")->fill(nrig, fMDst->richPhEl, weight);
			if (fMDst->richPhEl > 0.5) {
				if (sign>0) Hist::Head("hLp_NafPhPr")->fill(nrig, fMDst->richPhPr, weight);
				if (sign<0) Hist::Head("hLn_NafPhPr")->fill(nrig, fMDst->richPhPr, weight);
			}
		}
	
	
		Float_t tofQCut = ( 1.15 + 0.125 * std::erfc( 4.0 * (std::log(3.5)-std::log(1.0+nrig)) ) );
		if (fMDst->tofQ  < 0.8 || fMDst->tofQ > tofQCut) break;
		
		const Float_t tofQulLmt = 0.10; 
		const Float_t tofQuTh   = 0.40;
		const Float_t tofQlTh   = 0.30;
		Float_t tofQuCut =  (tofQulLmt + 0.5 * (tofQuTh-tofQulLmt) * std::erfc( 4.0 * (std::log(4.0)-std::log(1.0+nrig)) ) );
		Float_t tofQlCut = -(tofQulLmt + 0.5 * (tofQlTh-tofQulLmt) * std::erfc( 4.0 * (std::log(3.5)-std::log(1.0+nrig)) ) );
		Float_t tofQulCut = (std::fabs(tofQlCut) + fabs(tofQuCut)) * 0.5;
		if (std::fabs(fMDst->tofQul) > tofQulCut) break;
		//if (tofQul < tofQlCut || tofQul > tofQuCut) return false;
	
		// Mass Estimator (Scale)
		Int_t   nrigBin = AXnr.find(nrig);
		Float_t nrigScl = (nrigBin <= AXnr.nbin()) ? AXnr.bins(nrigBin) : AXnr.max();
		Float_t massScl = (4.0*std::exp(-0.2*nrigScl*nrigScl) + 1.0);
		Float_t tofM    = fMDst->tofM / massScl;
		
		const Int_t   N_TRDL = 1;
		const Float_t V_TRD[N_TRDL+1] = 
			{ 0.850,
		    0.850 };
		
		const Int_t   N_VETOL = 8;
		const Float_t V_VETOL[N_VETOL+1] = 
			{ 10,
			  3, 4, 5, 6, 7, 8, 9, 10 };
	
		std::vector< std::vector<Int_t> > Set;
		for (Int_t iTRD = 1; iTRD <= N_TRDL; ++iTRD) {
		for (Int_t iVETO = 1; iVETO <= N_VETOL; ++iVETO) {
			std::vector<Int_t> iSet( { iTRD, iVETO } );
			Set.push_back(iSet);
		}}

		for (Int_t iSet = 0; iSet <= Set.size(); ++iSet) {
			Int_t iTRD  = (iSet==Set.size()) ? 0 : Set.at(iSet).at(0);
			Int_t iVETO = (iSet==Set.size()) ? 0 : Set.at(iSet).at(1);
			std::string nmSet = (iSet==Set.size()) ? "" : StrFmt("%d%d", iTRD, iVETO);

			if (sign>0) hLpCutflow.at(iSet)->fill(nrig, 0, weight);
			if (sign<0) hLnCutflow.at(iSet)->fill(nrig, 0, weight);
			if (IsISS && sign>0) hTLpCutflow.at(iSet)->fill(uTime, nrig, 0, weight);
			if (IsISS && sign<0) hTLnCutflow.at(iSet)->fill(uTime, nrig, 0, weight);
		
			// TRD
			Bool_t isTrdPr = (fMDst->trdEst > V_TRD[iTRD]);
			
			// TRD
			if (!isTrdPr) continue;
			
			// RICH
			//if (fMDst->richRad != 0) continue;
	
			// RICH VETO
			if (isAglNORing && fMDst->richNHit >= V_VETOL[iVETO]) continue;  // 5 -> 15
			if (isNafNORing && fMDst->richNHit >= 2) continue;  // 4 -> 6

			// RICH RING
			Bool_t isNafDE = (isNafRing && fMDst->richMPr > 2.5);
			Bool_t isAglDE = (isAglRing && fMDst->richMPr > 2.5);
			Bool_t isNafPE = (isNafRing && fMDst->richMPr < -1.5 && fMDst->richMPi < 2.0);
			Bool_t isAglPE = (isAglRing && fMDst->richMPr < -1.5 && fMDst->richMPi < 2.0);
		
			// RICH NAF RING
			if (isNafDE || isAglDE) continue;
			if (isNafPE) continue;

			if (sign>0) hLpCutflow.at(iSet)->fill(nrig, 1, weight);
			if (sign<0) hLnCutflow.at(iSet)->fill(nrig, 1, weight);
			if (IsISS && sign>0) hTLpCutflow.at(iSet)->fill(uTime, nrig, 1, weight);
			if (IsISS && sign<0) hTLnCutflow.at(iSet)->fill(uTime, nrig, 1, weight);

			// Template Fit
			Bool_t isSignal     = (sign > 0) && (isTrdPr && !isAglPE);
			Bool_t isBackground = (sign < 0) && (isAglPE);
	
			if (isSignal)     hLsTofM.at(iSet)->fill(nrig, tofM, weight);
			if (isBackground) hLbTofM.at(iSet)->fill(nrig, tofM, weight);
			if (IsISS && isSignal)     hTLsTofM.at(iSet)->fill(uTime, nrig, tofM, weight);
			if (IsISS && isBackground) hTLbTofM.at(iSet)->fill(uTime, nrig, tofM, weight);
			
			// RICH AGL RING
			if (isAglPE) continue;
			
			if (sign>0) hLpCutflow.at(iSet)->fill(nrig, 2, weight);
			if (sign<0) hLnCutflow.at(iSet)->fill(nrig, 2, weight);
			if (IsISS && sign>0) hTLpCutflow.at(iSet)->fill(uTime, nrig, 2, weight);
			if (IsISS && sign<0) hTLnCutflow.at(iSet)->fill(uTime, nrig, 2, weight);

			if (sign>0) hLpCutflow.at(iSet)->fill(nrig, 3, weight);
			if (sign<0) hLnCutflow.at(iSet)->fill(nrig, 3, weight);
			if (IsISS && sign>0) hTLpCutflow.at(iSet)->fill(uTime, nrig, 3, weight);
			if (IsISS && sign<0) hTLnCutflow.at(iSet)->fill(uTime, nrig, 3, weight);

			if (sign>0) hLpTofM.at(iSet)->fill(nrig, tofM, weight);
			if (sign<0) hLnTofM.at(iSet)->fill(nrig, tofM, weight);
			if (IsISS && sign>0) hTLpTofM.at(iSet)->fill(uTime, nrig, tofM, weight);
			if (IsISS && sign<0) hTLnTofM.at(iSet)->fill(uTime, nrig, tofM, weight);

			if (sign>0) hLpEvt.at(iSet)->fill(nrig, weight);
			if (sign<0) hLnEvt.at(iSet)->fill(nrig, weight);
			if (IsISS && sign>0) hTLpEvt.at(iSet)->fill(uTime, nrig, weight);
			if (IsISS && sign<0) hTLnEvt.at(iSet)->fill(uTime, nrig, weight);
			
			if (IsMonteCarlo) {
				if (MCSign>0) hLpMCEvt.at(iSet)->fill(MCNRig, weight);
				if (MCSign<0) hLnMCEvt.at(iSet)->fill(MCNRig, weight);
			}
		
			if (iSet==Set.size()) fOptL = true;
		}
		break;	
	}	


	/**** Intermedia Energy Region ****/
	while (true) {
		Short_t     trPt  = 2;
		std::string trNm  = "In";
		Short_t     sign  = fMDst->sign[trPt-2];
		Float_t     nrig  = fMDst->nrig[trPt-2];
		Float_t     lchix = fMDst->lchix[trPt-2];
		Float_t     lchiy = fMDst->lchiy[trPt-2];
		Float_t     lasym = fMDst->lasym[trPt-2];
		Float_t     ccest = ccestPt[trPt-2];
		if (sign == 0) break;
		
		// RTI
		if ((nrig / fMDst->cfRig) < 1.2) break;
		
		// RICH
		Bool_t isAglRing  = (fMDst->richRad == 0 && fMDst->hasRich);
		if (!isAglRing) break;
		
		// Charge Confusion
		if (sign>0) Hist::Head("hIp_Lchix")->fill(nrig, lchix, weight);
		if (sign<0) Hist::Head("hIn_Lchix")->fill(nrig, lchix, weight);
		
		if (lchix > 2.7) break;
	
		if (sign>0) Hist::Head("hIp_Lchiy")->fill(nrig, lchiy, weight);
		if (sign<0) Hist::Head("hIn_Lchiy")->fill(nrig, lchiy, weight);
		
		if (lchiy > 2.4) break;
		
		if (sign>0) Hist::Head("hIp_Lasym")->fill(nrig, lasym, weight);
		if (sign<0) Hist::Head("hIn_Lasym")->fill(nrig, lasym, weight);
		
		if (lasym > 1.3) break;
		
		if (sign>0) Hist::Head("hIp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hIn_CCest")->fill(nrig, ccest, weight);
	
		const Int_t   N_MASSI = 1;
		const Float_t V_MASS[N_MASSI+1] = 
			{ 1.5,
			  1.5 };
	
		std::vector< std::vector<Int_t> > Set;
		for (Int_t iMASS = 1; iMASS <= N_MASSI; ++iMASS) {
			std::vector<Int_t> iSet( { iMASS } );
			Set.push_back(iSet);
		}

		for (Int_t iSet = 0; iSet <= Set.size(); ++iSet) {
			Int_t iMASS = (iSet==Set.size()) ? 0 : Set.at(iSet).at(0);
			std::string nmSet = (iSet==Set.size()) ? "" : StrFmt("%d", iMASS);

			if (sign>0) hIpCutflow.at(iSet)->fill(nrig, 0, weight);
			if (sign<0) hInCutflow.at(iSet)->fill(nrig, 0, weight);
			if (IsISS && sign>0) hTIpCutflow.at(iSet)->fill(uTime, nrig, 0, weight);
			if (IsISS && sign<0) hTInCutflow.at(iSet)->fill(uTime, nrig, 0, weight);
		
			// RICH
			Bool_t isAglDE = (isAglRing && fMDst->richMPr > 2.5);
			Bool_t isAglPE = (isAglRing && fMDst->richMPr < -V_MASS[iMASS] && fMDst->richMPi < 2.5);

			if (isAglDE) continue;

			// Template Fit
			Bool_t isSignal     = (sign > 0 && isEcalPr && !isAglPE);
			Bool_t isBackground = (sign < 0 && isEcalEl);
		
			if (isSignal)     hIsTrdL.at(iSet)->fill(nrig, trdEst, weight);
			if (isBackground) hIbTrdL.at(iSet)->fill(nrig, trdEst, weight);
			if (IsISS && isSignal)     hTIsTrdL.at(iSet)->fill(uTime, nrig, trdEst, weight);
			if (IsISS && isBackground) hTIbTrdL.at(iSet)->fill(uTime, nrig, trdEst, weight);

			// ECAL
			if (fMDst->hasEcal && !isEcalPr) continue;
			
			if (sign>0) hIpCutflow.at(iSet)->fill(nrig, 1, weight);
			if (sign<0) hInCutflow.at(iSet)->fill(nrig, 1, weight);
			if (IsISS && sign>0) hTIpCutflow.at(iSet)->fill(uTime, nrig, 1, weight);
			if (IsISS && sign<0) hTInCutflow.at(iSet)->fill(uTime, nrig, 1, weight);

			// RICH AGL RING
			if (isAglPE) continue; 
			
			if (sign>0) hIpCutflow.at(iSet)->fill(nrig, 2, weight);
			if (sign<0) hInCutflow.at(iSet)->fill(nrig, 2, weight);
			if (IsISS && sign>0) hTIpCutflow.at(iSet)->fill(uTime, nrig, 2, weight);
			if (IsISS && sign<0) hTInCutflow.at(iSet)->fill(uTime, nrig, 2, weight);
				
			if (sign>0) hIpTrdL.at(iSet)->fill(nrig, trdEst, weight);
			if (sign<0) hInTrdL.at(iSet)->fill(nrig, trdEst, weight);
			if (IsISS && sign>0) hTIpTrdL.at(iSet)->fill(uTime, nrig, trdEst, weight);
			if (IsISS && sign<0) hTInTrdL.at(iSet)->fill(uTime, nrig, trdEst, weight);
			
			if (sign>0) hIpEvt.at(iSet)->fill(nrig, weight);
			if (sign<0) hInEvt.at(iSet)->fill(nrig, weight);
			if (IsISS && sign>0) hTIpEvt.at(iSet)->fill(uTime, nrig, weight);
			if (IsISS && sign<0) hTInEvt.at(iSet)->fill(uTime, nrig, weight);

			if (IsMonteCarlo) {
				if (MCSign>0) hIpMCEvt.at(iSet)->fill(MCNRig, weight);
				if (MCSign<0) hInMCEvt.at(iSet)->fill(MCNRig, weight);
			}
			
			if (iSet==Set.size()) fOptI = true;
		}
		break;
	}


	while (true) {
		Short_t     trPt  = 2;
		std::string trNm  = "In";
		Short_t     sign  = fMDst->sign[trPt-2];
		Float_t     nrig  = fMDst->nrig[trPt-2];
		Float_t     lchix = fMDst->lchix[trPt-2];
		Float_t     lchiy = fMDst->lchiy[trPt-2];
		Float_t     lasym = fMDst->lasym[trPt-2];
		Float_t     ccest = ccestPt[trPt-2];
		if (sign == 0) break;
		
		// RTI
		if ((nrig / fMDst->cfRig) < 1.2) break;
		
		// Charge Confusion
		if (sign>0) Hist::Head("hMp_Lchix")->fill(nrig, lchix, weight);
		if (sign<0) Hist::Head("hMn_Lchix")->fill(nrig, lchix, weight);
		
		if (lchix > 2.7) break;
	
		if (sign>0) Hist::Head("hMp_Lchiy")->fill(nrig, lchiy, weight);
		if (sign<0) Hist::Head("hMn_Lchiy")->fill(nrig, lchiy, weight);
		
		if (lchiy > 2.4) break;
		
		if (sign>0) Hist::Head("hMp_Lasym")->fill(nrig, lasym, weight);
		if (sign<0) Hist::Head("hMn_Lasym")->fill(nrig, lasym, weight);
		
		if (lasym > 1.3) break;
		
		if (sign>0) Hist::Head("hMp_CCest")->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head("hMn_CCest")->fill(nrig, ccest, weight);
		
		// RICH
		Bool_t isAglRing  = (fMDst->richRad == 0 && fMDst->hasRich);

		
		const Int_t   N_MASSM = 1;
		const Float_t V_MASS[N_MASSM+1] = 
			{ 1.5, 
			  1.5 };
	
		std::vector< std::vector<Int_t> > Set;
		for (Int_t iMASS = 1; iMASS <= N_MASSM; ++iMASS) {
			std::vector<Int_t> iSet( { iMASS } );
			Set.push_back(iSet);
		}

		for (Int_t iSet = 0; iSet <= Set.size(); ++iSet) {
			Int_t iMASS = (iSet==Set.size()) ? 0 : Set.at(iSet).at(0);
			std::string nmSet = (iSet==Set.size()) ? "" : StrFmt("%d", iMASS);

			if (sign>0) hMpCutflow.at(iSet)->fill(nrig, 0, weight);
			if (sign<0) hMnCutflow.at(iSet)->fill(nrig, 0, weight);
			if (IsISS && sign>0) hTMpCutflow.at(iSet)->fill(uTime, nrig, 0, weight);
			if (IsISS && sign<0) hTMnCutflow.at(iSet)->fill(uTime, nrig, 0, weight);
		
			// RICH
			Bool_t isAglDE = (isAglRing && fMDst->richMPr > 2.5);
			Bool_t isAglPE = (isAglRing && fMDst->richMPr < -V_MASS[iMASS] && fMDst->richMPi < 3.0);

			if (isAglDE) continue;

			// Template Fit
			Bool_t isSignal     = (sign > 0 && isEcalPr && !isAglPE);
			Bool_t isBackground = (sign < 0 && isEcalEl);
		
			if (isSignal)     hMsTrdL.at(iSet)->fill(nrig, trdEst, weight);
			if (isBackground) hMbTrdL.at(iSet)->fill(nrig, trdEst, weight);
			if (IsISS && isSignal)     hTMsTrdL.at(iSet)->fill(uTime, nrig, trdEst, weight);
			if (IsISS && isBackground) hTMbTrdL.at(iSet)->fill(uTime, nrig, trdEst, weight);

			// ECAL
			if (fMDst->hasEcal && !isEcalPr) continue;
			
			if (sign>0) hMpCutflow.at(iSet)->fill(nrig, 1, weight);
			if (sign<0) hMnCutflow.at(iSet)->fill(nrig, 1, weight);
			if (IsISS && sign>0) hTMpCutflow.at(iSet)->fill(uTime, nrig, 1, weight);
			if (IsISS && sign<0) hTMnCutflow.at(iSet)->fill(uTime, nrig, 1, weight);

			// RICH AGL RING
			if (isAglPE) continue;
			
			if (sign>0) hMpCutflow.at(iSet)->fill(nrig, 2, weight);
			if (sign<0) hMnCutflow.at(iSet)->fill(nrig, 2, weight);
			if (IsISS && sign>0) hTMpCutflow.at(iSet)->fill(uTime, nrig, 2, weight);
			if (IsISS && sign<0) hTMnCutflow.at(iSet)->fill(uTime, nrig, 2, weight);
				
			if (sign>0) hMpTrdL.at(iSet)->fill(nrig, trdEst, weight);
			if (sign<0) hMnTrdL.at(iSet)->fill(nrig, trdEst, weight);
			if (IsISS && sign>0) hTMpTrdL.at(iSet)->fill(uTime, nrig, trdEst, weight);
			if (IsISS && sign<0) hTMnTrdL.at(iSet)->fill(uTime, nrig, trdEst, weight);
			
			if (sign>0) hMpEvt.at(iSet)->fill(nrig, weight);
			if (sign<0) hMnEvt.at(iSet)->fill(nrig, weight);
			if (IsISS && sign>0) hTMpEvt.at(iSet)->fill(uTime, nrig, weight);
			if (IsISS && sign<0) hTMnEvt.at(iSet)->fill(uTime, nrig, weight);

			if (IsMonteCarlo) {
				if (MCSign>0) hMpMCEvt.at(iSet)->fill(MCNRig, weight);
				if (MCSign<0) hMnMCEvt.at(iSet)->fill(MCNRig, weight);
			}
			
			if (iSet==Set.size()) fOptM = true;
		}
		break;
	}

	
	/**** High Energy Region ****/
	while (false) {
		Short_t     trPt  = fMDst->trPt;
		std::string trNm  = "In";
		if      (trPt == 2) trNm = "In";
		else if (trPt == 3) trNm = "L1";
		else if (trPt == 4) trNm = "L9";
		else if (trPt == 5) trNm = "Fs";
		Short_t     sign  = fMDst->trSign[trPt-2];
		Float_t     nrig  = fMDst->trNRig[trPt-2];
		Float_t     lchix = fMDst->trLchix[trPt-2];
		Float_t     ccest = ccestPt[trPt-2];
		Float_t     ccCut = 1.0 * (std::erf(std::log(tuneCCPnt[trPt-2]) - std::log(nrig)) + 1.0);
		Float_t     cfScl = (nrig / fMDst->cfRig);
		
		if (sign == 0) break;

		if (cfScl < 1.2) break;
		
		// TRK
		if (lchix > 3.) break;
		if (trPt == 2) break; 

		if (sign>0) Hist::Head(StrFmt("hHp_Cutflow%s", trNm.c_str()))->fill(nrig, 0, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_Cutflow%s", trNm.c_str()))->fill(nrig, 0, weight);
		if (IsISS && sign>0) Hist::Head(StrFmt("hTHp_Cutflow%s", trNm.c_str()))->fill(uTime, nrig, 0, weight);
		if (IsISS && sign<0) Hist::Head(StrFmt("hTHn_Cutflow%s", trNm.c_str()))->fill(uTime, nrig, 0, weight);
	
		// ECAL
		if (fMDst->hasEcal && !isEcalPr) break;

		if (sign>0) Hist::Head(StrFmt("hHp_Cutflow%s", trNm.c_str()))->fill(nrig, 1, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_Cutflow%s", trNm.c_str()))->fill(nrig, 1, weight);
		if (IsISS && sign>0) Hist::Head(StrFmt("hTHp_Cutflow%s", trNm.c_str()))->fill(uTime, nrig, 1, weight);
		if (IsISS && sign<0) Hist::Head(StrFmt("hTHn_Cutflow%s", trNm.c_str()))->fill(uTime, nrig, 1, weight);
		
		// TRD
		Bool_t isTrdPr = (fMDst->trdEst > 0.75);
		if (!isTrdPr) break;
		if (sign>0) Hist::Head(StrFmt("hHp_Cutflow%s", trNm.c_str()))->fill(nrig, 2, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_Cutflow%s", trNm.c_str()))->fill(nrig, 2, weight);
		if (IsISS && sign>0) Hist::Head(StrFmt("hTHp_Cutflow%s", trNm.c_str()))->fill(uTime, nrig, 2, weight);
		if (IsISS && sign<0) Hist::Head(StrFmt("hTHn_Cutflow%s", trNm.c_str()))->fill(uTime, nrig, 2, weight);



		// TRK
		Bool_t ccInCut2 = (ccestIn > 0.90 * ccestTh[0]); 
		Bool_t ccL1Cut2 = (ccestL1 > 0.90 * ccestTh[1]); 
		Bool_t ccL9Cut2 = (ccestL9 > 0.90 * ccestTh[2]); 
		Bool_t ccFsCut2 = (ccestFs > 0.90 * ccestTh[3]); 
		
		if (trPt == 3 && (ccInCut2)) break; 
		if (trPt == 4 && (ccInCut2)) break; 
		if (trPt == 5 && (ccL1Cut2 || ccL9Cut2)) break; 
		
		if (sign>0) Hist::Head(StrFmt("hHp_Cutflow%s", trNm.c_str()))->fill(nrig, 3, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_Cutflow%s", trNm.c_str()))->fill(nrig, 3, weight);
		if (IsISS && sign>0) Hist::Head(StrFmt("hTHp_Cutflow%s", trNm.c_str()))->fill(uTime, nrig, 3, weight);
		if (IsISS && sign<0) Hist::Head(StrFmt("hTHn_Cutflow%s", trNm.c_str()))->fill(uTime, nrig, 3, weight);
	
		// CCest
		if (sign>0) Hist::Head(StrFmt("hHp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		
		if (IsISS && sign>0) Hist::Head(StrFmt("hTHp_CCest%s", trNm.c_str()))->fill(uTime, nrig, ccest, weight);
		if (IsISS && sign<0) Hist::Head(StrFmt("hTHn_CCest%s", trNm.c_str()))->fill(uTime, nrig, ccest, weight);
		
		if (sign>0) Hist::Head(StrFmt("hHp_Evt%s", trNm.c_str()))->fill(nrig, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_Evt%s", trNm.c_str()))->fill(nrig, weight);
		
		if (IsISS && sign>0) Hist::Head(StrFmt("hTHp_Evt%s", trNm.c_str()))->fill(uTime, nrig, weight);
		if (IsISS && sign<0) Hist::Head(StrFmt("hTHn_Evt%s", trNm.c_str()))->fill(uTime, nrig, weight);
		
		if (IsMonteCarlo) {
			if (MCSign>0) Hist::Head(StrFmt("hHp_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
			if (MCSign<0) Hist::Head(StrFmt("hHn_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
		}

		if (trPt == 3) fOptHL1 = true;
		if (trPt == 4) fOptHL9 = true;
		if (trPt == 5) fOptHFs = true;
		break;
	}


	Bool_t isSave = (fOptL || fOptI || fOptM || fOptHL1 || fOptHL9 || fOptHFs);
	if (!isSave) return false;

	fMCRig = MCSign * MCNRig;
	fTRRig = fMDst->trSign[0] * fMDst->trNRig[0];

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
