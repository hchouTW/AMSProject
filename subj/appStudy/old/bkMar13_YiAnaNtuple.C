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
				fDST->Branch("MCRig", &fMCRig);
				fDST->Branch("TRRig", &fTRRig);
				fDST->Branch("OptL" , &fOptL);
				fDST->Branch("OptIV" , &fOptIV);
				fDST->Branch("OptIW" , &fOptIW);
				fDST->Branch("OptH" , &fOptH);
			}
			file->cd();
			
			Axis AXrso("Rigidity Resolution", 400, -1.5, 1.5);

			// rigidity binning
			AXnr = Axis("|Rigidity| [GV]",
				{   1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
				    3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
					  8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
					 19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
					 41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
					 93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } );
			AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
			
			//AXnr = Axis("|Rigidity| [GV]",
			//	{   1.00,   2.00,   3.00,   4.12, 
			//	    5.00,   6.00,   7.10,   8.30,   9.62, 
			//	   11.04,  12.59,  14.25,  16.05,  17.98, 
			//	   20.04,  22.25,  24.62,  27.25,  30.21  } );
			//AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
			
			// 3 month
			AXtme = Axis("Time", 
				{ 1305417600, 1312416000, 1319414400, 1326412800, 1333411200,
				  1340409600, 1347408000, 1354406400, 1361404800, 1368403200,
					1375401600, 1382400000, 1389398400, 1396396800, 1403395200,
					1410393600, 1417392000, 1424390400, 1431388800, 1438387200,
					1445385600, 1452384000, 1459382400, 1466380800, 1473379200, 
					1480377600 } );
	
			//----  All Energy  ----//
			//Axis AXCOMtrQ("Chrg", 100, 0.5, 3.0);
			//Hist::New("hCOMp_TrMxQ", "", AXnr, AXCOMtrQ);
			//Hist::New("hCOMn_TrMxQ", "", AXnr, AXCOMtrQ);
			//Hist::New("hCOMp_TrL2Q", "", AXnr, AXCOMtrQ);
			//Hist::New("hCOMn_TrL2Q", "", AXnr, AXCOMtrQ);
			//Hist::New("hCOMp_TrL1Q", "", AXnr, AXCOMtrQ);
			//Hist::New("hCOMn_TrL1Q", "", AXnr, AXCOMtrQ);
			//Hist::New("hCOMp_TrL9Q", "", AXnr, AXCOMtrQ);
			//Hist::New("hCOMn_TrL9Q", "", AXnr, AXCOMtrQ);
			
			Axis AXCOMlchi("lchi", 100, -4., 4.);
			Hist::New("hCOMp_LchiyIn", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyIn", "", AXnr, AXCOMlchi);
			Hist::New("hCOMp_LchiyL1", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyL1", "", AXnr, AXCOMlchi);
			Hist::New("hCOMp_LchiyL9", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyL9", "", AXnr, AXCOMlchi);
			Hist::New("hCOMp_LchiyFs", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyFs", "", AXnr, AXCOMlchi);

			Axis AXCOMlasym("Lasym", 100, -4., 3.);
			Hist::New("hCOMp_LasymUL", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_LasymUL", "", AXnr, AXCOMlasym);
			Hist::New("hCOMp_Lasym1I", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym1I", "", AXnr, AXCOMlasym);
			Hist::New("hCOMp_Lasym9I", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym9I", "", AXnr, AXCOMlasym);
			Hist::New("hCOMp_Lasym91", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym91", "", AXnr, AXCOMlasym);
			
			Hist::New("hCOMp_CCIn", "", AXCOMlchi, AXCOMlasym);
			Hist::New("hCOMn_CCIn", "", AXCOMlchi, AXCOMlasym);
			Hist::New("hCOMp_CCL1", "", AXCOMlchi, AXCOMlasym);
			Hist::New("hCOMn_CCL1", "", AXCOMlchi, AXCOMlasym);
			Hist::New("hCOMp_CCL9", "", AXCOMlchi, AXCOMlasym);
			Hist::New("hCOMn_CCL9", "", AXCOMlchi, AXCOMlasym);
			Hist::New("hCOMp_CCFs", "", AXCOMlchi, AXCOMlasym);
			Hist::New("hCOMn_CCFs", "", AXCOMlchi, AXCOMlasym);

			Axis AXCOMccest("CCest", 400, -6.0, 8.0);
			Hist::New("hCOMp_CCestIn", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestIn", "", AXnr, AXCOMccest);
			Hist::New("hCOMp_CCestL1", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestL1", "", AXnr, AXCOMccest);
			Hist::New("hCOMp_CCestL9", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestL9", "", AXnr, AXCOMccest);
			Hist::New("hCOMp_CCestFs", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestFs", "", AXnr, AXCOMccest);

			//----  Low Energy  ----//
			Axis AXLcutflow("Cutflow", 8, 0., 8.);
			Hist::New("hLp_Cutflow", "", AXnr, AXLcutflow);
			Hist::New("hLn_Cutflow", "", AXnr, AXLcutflow);
	
			Axis AXLaglM("Mass Estimator", 400, -0.08, 0.08);
			Hist::New("hLp_AglMPr", "", AXnr, AXLaglM);
			Hist::New("hLn_AglMPr", "", AXnr, AXLaglM);
			
			Hist::New("hLp_AglMPi", "", AXnr, AXLaglM);
			Hist::New("hLn_AglMPi", "", AXnr, AXLaglM);
			
			Hist::New("hLp_AglMEl", "", AXnr, AXLaglM);
			Hist::New("hLn_AglMEl", "", AXnr, AXLaglM);
			
			Axis AXLnafM("Mass Estimator", 400, -0.12, 0.12);
			Hist::New("hLp_NafMPr", "", AXnr, AXLnafM);
			Hist::New("hLn_NafMPr", "", AXnr, AXLnafM);
			
			Hist::New("hLp_NafMPi", "", AXnr, AXLnafM);
			Hist::New("hLn_NafMPi", "", AXnr, AXLnafM);
			
			Hist::New("hLp_NafMEl", "", AXnr, AXLnafM);
			Hist::New("hLn_NafMEl", "", AXnr, AXLnafM);
			
			Axis AXLaglPrPh("PrPh", 200, 0., 15.);
			Hist::New("hLp_AglPrPh", "", AXnr, AXLaglPrPh);
			Hist::New("hLn_AglPrPh", "", AXnr, AXLaglPrPh);
			
			Axis AXLnafPrPh("PrPh", 100, 0., 5.);
			Hist::New("hLp_NafPrPh", "", AXnr, AXLnafPrPh);
			Hist::New("hLn_NafPrPh", "", AXnr, AXLnafPrPh);
			
			Axis AXLaglElPh("ElPh", 200, 0., 15.);
			Hist::New("hLp_AglElPh", "", AXnr, AXLaglElPh);
			Hist::New("hLn_AglElPh", "", AXnr, AXLaglElPh);
			
			Axis AXLnafElPh("ElPh", 100, 0., 5.);
			Hist::New("hLp_NafElPh", "", AXnr, AXLnafElPh);
			Hist::New("hLn_NafElPh", "", AXnr, AXLnafElPh);
			
			Axis AXLaglPh("Ph", 50, 0., 50.);
			Hist::New("hLp_AglPh", "", AXnr, AXLaglPh);
			Hist::New("hLn_AglPh", "", AXnr, AXLaglPh);
			
			Axis AXLnafPh("Ph", 25, 0., 25.);
			Hist::New("hLp_NafPh", "", AXnr, AXLnafPh);
			Hist::New("hLn_NafPh", "", AXnr, AXLnafPh);
			
			Axis AXLtofM("Mass Estimator", 200, -1.1, 0.35); // best 200
			Hist::New("hLs_TofM", "", AXnr, AXLtofM);
			Hist::New("hLb_TofM", "", AXnr, AXLtofM);
			
			Hist::New("hLp_TofM", "", AXnr, AXLtofM);
			Hist::New("hLn_TofM", "", AXnr, AXLtofM);
			
			Axis AXLtofMScl("Mass Estimator (Scale)", 200, -4.0, 4.0); // best 200
			Hist::New("hLs_TofMScl", "", AXnr, AXLtofMScl);
			Hist::New("hLb_TofMScl", "", AXnr, AXLtofMScl);
			
			Hist::New("hLp_TofMScl", "", AXnr, AXLtofMScl);
			Hist::New("hLn_TofMScl", "", AXnr, AXLtofMScl);

			Hist::New("hLp_Evt", "", AXnr);
			Hist::New("hLn_Evt", "", AXnr);
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hLp_MCEvt", "", AXnr);
				Hist::New("hLn_MCEvt", "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Hist::New("hTLp_Cutflow", "", AXtme, AXnr, AXLcutflow);
				Hist::New("hTLn_Cutflow", "", AXtme, AXnr, AXLcutflow);
				
				Hist::New("hTLs_TofM", "", AXtme, AXnr, AXLtofM);
				Hist::New("hTLb_TofM", "", AXtme, AXnr, AXLtofM);
				
				Hist::New("hTLp_TofM", "", AXtme, AXnr, AXLtofM);
				Hist::New("hTLn_TofM", "", AXtme, AXnr, AXLtofM);
				
				Hist::New("hTLs_TofMScl", "", AXtme, AXnr, AXLtofMScl);
				Hist::New("hTLb_TofMScl", "", AXtme, AXnr, AXLtofMScl);
				
				Hist::New("hTLp_TofMScl", "", AXtme, AXnr, AXLtofMScl);
				Hist::New("hTLn_TofMScl", "", AXtme, AXnr, AXLtofMScl);
				
				Hist::New("hTLp_Evt", "", AXtme, AXnr);
				Hist::New("hTLn_Evt", "", AXtme, AXnr);
			}
			
			//----  Intermedia Energy  ----//
			Axis AXIcutflow("Cutflow", 3, 0., 3.);
			Hist::New("hIp_CutflowV", "", AXnr, AXIcutflow);
			Hist::New("hIn_CutflowV", "", AXnr, AXIcutflow);
			
			Hist::New("hIp_CutflowW", "", AXnr, AXIcutflow);
			Hist::New("hIn_CutflowW", "", AXnr, AXIcutflow);
			
			Axis AXItrdEst("TrdEst", 100, 0.2, 1.6);
			Hist::New("hIs_TrdEstV", "", AXnr, AXItrdEst);
			Hist::New("hIb_TrdEstV", "", AXnr, AXItrdEst);
			
			Hist::New("hIp_TrdEstV", "", AXnr, AXItrdEst);
			Hist::New("hIn_TrdEstV", "", AXnr, AXItrdEst);
			
			Hist::New("hIs_TrdEstW", "", AXnr, AXItrdEst);
			Hist::New("hIb_TrdEstW", "", AXnr, AXItrdEst);
			
			Hist::New("hIp_TrdEstW", "", AXnr, AXItrdEst);
			Hist::New("hIn_TrdEstW", "", AXnr, AXItrdEst);
			
			Hist::New("hIp_EvtV", "", AXnr);
			Hist::New("hIn_EvtV", "", AXnr);
			
			Hist::New("hIp_EvtW", "", AXnr);
			Hist::New("hIn_EvtW", "", AXnr);
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hIp_MCEvtV", "", AXnr);
				Hist::New("hIn_MCEvtV", "", AXnr);
				
				Hist::New("hIp_MCEvtW", "", AXnr);
				Hist::New("hIn_MCEvtW", "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Hist::New("hTIp_CutflowV", "", AXtme, AXnr, AXIcutflow);
				Hist::New("hTIn_CutflowV", "", AXtme, AXnr, AXIcutflow);
				
				Hist::New("hTIp_CutflowW", "", AXtme, AXnr, AXIcutflow);
				Hist::New("hTIn_CutflowW", "", AXtme, AXnr, AXIcutflow);
				
				Hist::New("hTIs_TrdEstV", "", AXtme, AXnr, AXItrdEst);
				Hist::New("hTIb_TrdEstV", "", AXtme, AXnr, AXItrdEst);
				
				Hist::New("hTIp_TrdEstV", "", AXtme, AXnr, AXItrdEst);
				Hist::New("hTIn_TrdEstV", "", AXtme, AXnr, AXItrdEst);
				
				Hist::New("hTIs_TrdEstW", "", AXtme, AXnr, AXItrdEst);
				Hist::New("hTIb_TrdEstW", "", AXtme, AXnr, AXItrdEst);
				
				Hist::New("hTIp_TrdEstW", "", AXtme, AXnr, AXItrdEst);
				Hist::New("hTIn_TrdEstW", "", AXtme, AXnr, AXItrdEst);
				
				Hist::New("hTIp_EvtV", "", AXtme, AXnr);
				Hist::New("hTIn_EvtV", "", AXtme, AXnr);
			
				Hist::New("hTIp_EvtW", "", AXtme, AXnr);
				Hist::New("hTIn_EvtW", "", AXtme, AXnr);
			}

			//----  High Energy  ----//
			Axis AXHcutflow("Cutflow", 4, 0., 4.);
			Hist::New("hHp_Cutflow", "", AXnr, AXHcutflow);
			Hist::New("hHn_Cutflow", "", AXnr, AXHcutflow);
			
			Axis AXHccest("CCest", 50, -4, 4);
			Hist::New("hHp_CCest", "", AXnr, AXHccest);
			Hist::New("hHn_CCest", "", AXnr, AXHccest);
			Hist::New("hHp_CCestL1", "", AXnr, AXHccest);
			Hist::New("hHn_CCestL1", "", AXnr, AXHccest);
			Hist::New("hHp_CCestL9", "", AXnr, AXHccest);
			Hist::New("hHn_CCestL9", "", AXnr, AXHccest);
			Hist::New("hHp_CCestFs", "", AXnr, AXHccest);
			Hist::New("hHn_CCestFs", "", AXnr, AXHccest);
			
			Hist::New("hHp_Evt", "", AXnr);
			Hist::New("hHn_Evt", "", AXnr);
			Hist::New("hHp_EvtL1", "", AXnr);
			Hist::New("hHn_EvtL1", "", AXnr);
			Hist::New("hHp_EvtL9", "", AXnr);
			Hist::New("hHn_EvtL9", "", AXnr);
			Hist::New("hHp_EvtFs", "", AXnr);
			Hist::New("hHn_EvtFs", "", AXnr);
			
			if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
				Hist::New("hHp_MCEvt", "", AXnr);
				Hist::New("hHn_MCEvt", "", AXnr);
				
				Hist::New("hHp_MCEvtL1", "", AXnr);
				Hist::New("hHn_MCEvtL1", "", AXnr);
				
				Hist::New("hHp_MCEvtL9", "", AXnr);
				Hist::New("hHn_MCEvtL9", "", AXnr);
				
				Hist::New("hHp_MCEvtFs", "", AXnr);
				Hist::New("hHn_MCEvtFs", "", AXnr);
			}
			
			if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
				Hist::New("hTHp_Cutflow", "", AXtme, AXnr, AXHcutflow);
				Hist::New("hTHn_Cutflow", "", AXtme, AXnr, AXHcutflow);
				
				Hist::New("hTHp_CCest", "", AXtme, AXnr, AXHccest);
				Hist::New("hTHn_CCest", "", AXtme, AXnr, AXHccest);
				Hist::New("hTHp_CCestL1", "", AXtme, AXnr, AXHccest);
				Hist::New("hTHn_CCestL1", "", AXtme, AXnr, AXHccest);
				Hist::New("hTHp_CCestL9", "", AXtme, AXnr, AXHccest);
				Hist::New("hTHn_CCestL9", "", AXtme, AXnr, AXHccest);
				Hist::New("hTHp_CCestFs", "", AXtme, AXnr, AXHccest);
				Hist::New("hTHn_CCestFs", "", AXtme, AXnr, AXHccest);
				
				Hist::New("hTHp_Evt", "", AXtme, AXnr);
				Hist::New("hTHn_Evt", "", AXtme, AXnr);
				Hist::New("hTHp_EvtL1", "", AXtme, AXnr);
				Hist::New("hTHn_EvtL1", "", AXtme, AXnr);
				Hist::New("hTHp_EvtL9", "", AXtme, AXnr);
				Hist::New("hTHn_EvtL9", "", AXtme, AXnr);
				Hist::New("hTHp_EvtFs", "", AXtme, AXnr);
				Hist::New("hTHn_EvtFs", "", AXtme, AXnr);
			}

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
			fOptL = false;
			fOptIV = false;
			fOptIW = false;
			fOptH = false;
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
		Bool_t  fOptIV;
		Bool_t  fOptIW;
		Bool_t  fOptH;

	public :
		static Bool_t gSaveDST;
		static Bool_t gSaveDSTTree;
		static Bool_t gSaveDSTClone;
		static UInt_t gUTime[2]; // (pre, cur)

	public :
		static Axis AXnr;
		static Axis AXir;
		static Axis AXtme;
		static const Float_t CfStableFT = 1.2;

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
Axis DST::AXtme;

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
	if (fileList.at(0).find("MC_Pr_PL1_0_510") != std::string::npos) DST::gSaveDSTTree = true;
	if (fileList.at(0).find("MC_Pr_PL1_1800") != std::string::npos) DST::gSaveDSTTree = true;
	if (fileList.at(0).find("MC_Ap_PL1_0_510") != std::string::npos) DST::gSaveDSTTree = true;
	if (fileList.at(0).find("MC_Ap_PL1_1800") != std::string::npos) DST::gSaveDSTTree = true;
	
	initDST(fFile, fRunChain, fDataChain);
	fFile->cd();
	/*******************/
	/**               **/
	/*******************/

	Long64_t nentries = fDataChain->GetEntries();
	Long64_t npassed = 0;
	Long64_t nprocessed = 0;
	Long64_t printRate = nentries / 50;
	if (printRate < 500000) printRate = 500000;
	if (printRate > 5000000) printRate = 5000000;

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
	const Float_t MPr = 0.938272297; 
	const Float_t MPi = 0.139570180; 
	const Float_t MEl = 0.000510999; 

	// Event
	Float_t weight = fMDst->weight;
	
	// RTI
	Bool_t  IsISS = (YiNtuple::CheckEventMode(YiNtuple::ISS) && (fMDst->mcSign == 0));
	UInt_t  uTime = fMDst->uTime;
	if (IsISS && fMDst->cfScl < 1.2) return false;
		
	// MC
	Bool_t  IsMonteCarlo = (YiNtuple::CheckEventMode(YiNtuple::MC) && (fMDst->mcSign != 0));
	Short_t MCSign       = fMDst->mcSign;
	Float_t MCNRig       = fMDst->mcNRig;
	
	// TRK
	std::string trNm = "";
	Short_t     trPt = fMDst->trPt;
	if      (trPt == 5) trNm = "Fs"; 
	else if (trPt == 4) trNm = "L9"; 
	else if (trPt == 3) trNm = "L1"; 
	else if (trPt == 2) trNm = "In"; 

	Short_t sign = fMDst->trSign;
	Float_t nrig = fMDst->trNRig; 
	if (fMDst->trLchix > 3.0) return false; 
	
	// TRD
	Float_t trdEst  = fMDst->trdEst;
	Bool_t  isTrdPr = (fMDst->trdEst > (0.50 + 0.25 / std::sqrt(1.+(MPr/nrig)*(MPr/nrig))));

	// ECAL
	Bool_t isEcalPr = (!fMDst->hasEcal) ? false : (fMDst->ecalEst < -0.8); 
	Bool_t isEcalEl = (!fMDst->hasEcal) ? false : (fMDst->ecalEst >  0.6); 

	// RICH
	Float_t richMSgm  = (fMDst->richRad == 0) ? 0.0025 : 0.0075;
	Float_t richMPrTh = richMSgm * (1.0 + 100.0*(std::sqrt(1.0+(MPr/nrig)*(MPr/nrig)) - 1.0));
	Float_t richMPiTh = richMSgm * (1.0 + 100.0*(std::sqrt(1.0+(MPi/nrig)*(MPi/nrig)) - 1.0));
	Float_t richMElTh = richMSgm * (1.0 + 100.0*(std::sqrt(1.0+(MEl/nrig)*(MEl/nrig)) - 1.0));
	
	Bool_t  isRichPr    = (fMDst->hasRich && std::fabs(fMDst->richMPr) < 2.0 * richMPrTh);
	Bool_t  isRichPE    = (fMDst->hasRich && fMDst->richMPi < 1.5 * richMPiTh);
	Bool_t  isNotPrRing = (fMDst->hasRich && fMDst->richMPr < -2.0 * richMPrTh);

	// IsProtonLike (TRD ECAL RICH)
	Bool_t isPrLike = (isTrdPr && (fMDst->hasEcal ? isEcalPr : true) && (fMDst->hasRich ? isRichPr : true));

	// TRK Charge
	Bool_t hasMxQ = (MgntNum::Compare(fMDst->trMxQ) > 0);
	Bool_t hasL2Q = (MgntNum::Compare(fMDst->trL2Q) > 0);
	Bool_t hasL1Q = (MgntNum::Compare(fMDst->trL1Q) > 0);
	Bool_t hasL9Q = (MgntNum::Compare(fMDst->trL9Q) > 0);
	Float_t trMxQ = (hasMxQ) ? fMDst->trMxQ : -1.0;
	Float_t trL2Q = (hasL2Q) ? fMDst->trL2Q : -1.0;
	Float_t trL1Q = (hasL1Q) ? fMDst->trL1Q : -1.0;
	Float_t trL9Q = (hasL9Q) ? fMDst->trL9Q : -1.0;

	Bool_t trMxQCut = (hasMxQ && (trMxQ < 1.0 || trMxQ > 2.3)); // (1.0 2.3) 
	Bool_t trL2QCut = (hasL2Q && (trL2Q < 0.8 || trL2Q > 2.3)); // (0.8 2.3)
	Bool_t trL1QCut = (hasL1Q && (trL1Q < 0.8 || trL1Q > 2.3)); // (0.8 2.3)
	Bool_t trL9QCut = (hasL9Q && (trL9Q < 0.8 || trL9Q > 2.3)); // (0.8 2.3)
	
	//if (isPrLike) {
	//	if (hasMxQ && sign>0) Hist::Head("hCOMp_TrMxQ")->fill(nrig, trMxQ, weight);
	//	if (hasMxQ && sign<0) Hist::Head("hCOMn_TrMxQ")->fill(nrig, trMxQ, weight);
	//	if (hasL2Q && sign>0) Hist::Head("hCOMp_TrL2Q")->fill(nrig, trL2Q, weight);
	//	if (hasL2Q && sign<0) Hist::Head("hCOMn_TrL2Q")->fill(nrig, trL2Q, weight);
	//	if (hasL1Q && sign>0) Hist::Head("hCOMp_TrL1Q")->fill(nrig, trL1Q, weight);
	//	if (hasL1Q && sign<0) Hist::Head("hCOMn_TrL1Q")->fill(nrig, trL1Q, weight);
	//	if (hasL9Q && sign>0) Hist::Head("hCOMp_TrL9Q")->fill(nrig, trL9Q, weight);
	//	if (hasL9Q && sign<0) Hist::Head("hCOMn_TrL9Q")->fill(nrig, trL9Q, weight);
	//}

	if (trMxQCut) return false;
	if (trL2QCut) return false;
	if (trL1QCut) return false;
	if (trL9QCut) return false;
	
	// Pre-CC Cut
	Float_t lchiyTh[4] = { 2.90, 2.20, 2.30, 1.90 }; // ~98%
	Float_t lasymTh[4] = { 1.00, 1.60, 1.40, 1.40 }; // ~98%
	Float_t lchiyIn = (fMDst->trLchiy[0] / lchiyTh[0]);
	Float_t lchiyL1 = (fMDst->trLchiy[1] / lchiyTh[1]);
	Float_t lchiyL9 = (fMDst->trLchiy[2] / lchiyTh[2]);
	Float_t lchiyFs = (fMDst->trLchiy[3] / lchiyTh[3]);
	Float_t lasymUL = ((2.0 * fMDst->trLasym[0] - fMDst->trLchiy[0]) / lasymTh[0]);
	Float_t lasym1I = ((2.0 * fMDst->trLasym[1] - fMDst->trLchiy[1]) / lasymTh[1]);
	Float_t lasym9I = ((2.0 * fMDst->trLasym[2] - fMDst->trLchiy[2]) / lasymTh[2]);
	Float_t lasym91 = ((2.0 * fMDst->trLasym[3] - fMDst->trLchiy[3]) / lasymTh[3]);
	
	if (isPrLike) {
		if (trPt==2 && sign>0) Hist::Head("hCOMp_LchiyIn")->fill(nrig, lchiyIn, weight);
		if (trPt==2 && sign<0) Hist::Head("hCOMn_LchiyIn")->fill(nrig, lchiyIn, weight);
		if (trPt==2 && sign>0) Hist::Head("hCOMp_LasymUL")->fill(nrig, lasymUL, weight);
		if (trPt==2 && sign<0) Hist::Head("hCOMn_LasymUL")->fill(nrig, lasymUL, weight);
	
		if (trPt==2 && sign>0 && nrig > 20) Hist::Head("hCOMp_CCIn")->fill(lchiyIn, lasymUL, weight);
		if (trPt==2 && sign<0 && nrig > 20) Hist::Head("hCOMn_CCIn")->fill(lchiyIn, lasymUL, weight);
		
		if (trPt==3 && sign>0) Hist::Head("hCOMp_LchiyL1")->fill(nrig, lchiyL1, weight);
		if (trPt==3 && sign<0) Hist::Head("hCOMn_LchiyL1")->fill(nrig, lchiyL1, weight);
		if (trPt==3 && sign>0) Hist::Head("hCOMp_Lasym1I")->fill(nrig, lasym1I, weight);
		if (trPt==3 && sign<0) Hist::Head("hCOMn_Lasym1I")->fill(nrig, lasym1I, weight);
		
		if (trPt==3 && sign>0 && nrig > 40) Hist::Head("hCOMp_CCL1")->fill(lchiyL1, lasym1I, weight);
		if (trPt==3 && sign<0 && nrig > 40) Hist::Head("hCOMn_CCL1")->fill(lchiyL1, lasym1I, weight);
		
		if (trPt==4 && sign>0) Hist::Head("hCOMp_LchiyL9")->fill(nrig, lchiyL9, weight);
		if (trPt==4 && sign<0) Hist::Head("hCOMn_LchiyL9")->fill(nrig, lchiyL9, weight);
		if (trPt==4 && sign>0) Hist::Head("hCOMp_Lasym9I")->fill(nrig, lasym9I, weight);
		if (trPt==4 && sign<0) Hist::Head("hCOMn_Lasym9I")->fill(nrig, lasym9I, weight);
		
		if (trPt==4 && sign>0 && nrig > 50) Hist::Head("hCOMp_CCL9")->fill(lchiyL9, lasym9I, weight);
		if (trPt==4 && sign<0 && nrig > 50) Hist::Head("hCOMn_CCL9")->fill(lchiyL9, lasym9I, weight);
		
		if (trPt==5 && sign>0) Hist::Head("hCOMp_LchiyFs")->fill(nrig, lchiyFs, weight);
		if (trPt==5 && sign<0) Hist::Head("hCOMn_LchiyFs")->fill(nrig, lchiyFs, weight);
		if (trPt==5 && sign>0) Hist::Head("hCOMp_Lasym91")->fill(nrig, lasym91, weight);
		if (trPt==5 && sign<0) Hist::Head("hCOMn_Lasym91")->fill(nrig, lasym91, weight);
		
		if (trPt==5 && sign>0 && nrig > 80) Hist::Head("hCOMp_CCFs")->fill(lchiyFs, lasym91, weight);
		if (trPt==5 && sign<0 && nrig > 80) Hist::Head("hCOMn_CCFs")->fill(lchiyFs, lasym91, weight);
	}

	Float_t ccestSf[4] = { 0.20, 0.27, 0.21, 0.27 };
	Float_t ccestSg[4] = { 0.40, 0.46, 0.47, 0.50 };
	Float_t ccestIn = ((((lchiyIn>0&&lasymUL>0) ? std::sqrt(lchiyIn*lchiyIn+lasymUL*lasymUL) : std::max(lchiyIn, lasymUL)) - ccestSf[0]) / ccestSg[0]);
	Float_t ccestL1 = ((((lchiyL1>0&&lasym1I>0) ? std::sqrt(lchiyL1*lchiyL1+lasym1I*lasym1I) : std::max(lchiyL1, lasym1I)) - ccestSf[1]) / ccestSg[1]);
	Float_t ccestL9 = ((((lchiyL9>0&&lasym9I>0) ? std::sqrt(lchiyL9*lchiyL9+lasym9I*lasym9I) : std::max(lchiyL9, lasym9I)) - ccestSf[2]) / ccestSg[2]);
	Float_t ccestFs = ((((lchiyFs>0&&lasym91>0) ? std::sqrt(lchiyFs*lchiyFs+lasym91*lasym91) : std::max(lchiyFs, lasym91)) - ccestSf[3]) / ccestSg[3]);

	Float_t ccestTh[4] = { (0.9*std::erfc(std::log(50.0)-std::log(1+nrig))+2.4), 2.4, 2.4, 2.3 }; // ~98%
	Bool_t ccInCut = (ccestIn > ccestTh[0]); 
	Bool_t ccL1Cut = (ccestL1 > ccestTh[1]); 
	Bool_t ccL9Cut = (ccestL9 > ccestTh[2]); 
	Bool_t ccFsCut = (ccestFs > ccestTh[3]); 
	
	//if (trPt == 3 && (ccInCut)) return false; 
	//if (trPt == 4 && (ccInCut)) return false; 
	//if (trPt == 5 && (ccL1Cut || ccL9Cut)) return false; 
	
	// Charge-Confusion Estimator
	Float_t ccest = 0.0;
	if      (trPt == 2) { ccest = ccestIn; }
	else if (trPt == 3) { ccest = ccestL1; }
	else if (trPt == 4) { ccest = ccestL9; }
	else if (trPt == 5) { ccest = ccestFs; }

	if (isPrLike) {
		if (sign>0) Hist::Head(StrFmt("hCOMp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hCOMn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
	}
	
	// CC Cut
	Float_t tuneCCPnt = 0.0;
	if      (trPt == 2) tuneCCPnt = 20.0;
	else if (trPt == 3) tuneCCPnt = 25.0;
	else if (trPt == 4) tuneCCPnt = 30.0;
	else if (trPt == 5) tuneCCPnt = 40.0;
	Float_t ccestCut = 1.0 * (std::erf(std::log(tuneCCPnt) - std::log(nrig)) + 1.0);

	//Float_t lchiyTh[4] = { 4.70, 3.90, 2.70, 2.50 }; // ~99%
	//Float_t lasymTh[4] = { 1.70, 1.90, 1.90, 1.90 }; // ~99%
	//Bool_t ccInCut = (lchiyIn > lchiyTh[0] || lasymUL > lasymTh[0]) && 
	//                 (lchiyIn > 0 && lasymUL > 0 && 
	//								  (lchiyIn*lchiyIn/lchiyTh[0]/lchiyTh[0] + lasymUL*lasymUL/lasymTh[0]/lasymTh[0]) > 1);
	//Bool_t ccL1Cut = (lchiyL1 > lchiyTh[1] || lasym1I > lasymTh[1]);
	//                 (lchiyL1 > 0 && lasym1I > 0 && 
	//								  (lchiyL1*lchiyL1/lchiyTh[1]/lchiyTh[1] + lasym1I*lasym1I/lasymTh[1]/lasymTh[1]) > 1);
	//Bool_t ccL9Cut = (lchiyL9 > lchiyTh[2] || lasym9I > lasymTh[2]);
	//                 (lchiyL9 > 0 && lasym9I > 0 && 
	//								  (lchiyL9*lchiyL9/lchiyTh[2]/lchiyTh[2] + lasym9I*lasym9I/lasymTh[2]/lasymTh[2]) > 1);
	//
	//if (trPt == 3 && (ccInCut)) return false; 
	//if (trPt == 4 && (ccInCut)) return false; 
	//if (trPt == 5 && (ccL1Cut || ccL9Cut)) return false; 
	

	// Charge-Confusion Estimator
	//Float_t lchiyCut98 = 0.0; // ~98%
	//if      (trPt == 2) lchiyCut98 = 2.90;
	//else if (trPt == 3) lchiyCut98 = 2.20;
	//else if (trPt == 4) lchiyCut98 = 2.30;
	//else if (trPt == 5) lchiyCut98 = 1.90;
	//
	//Float_t lasymCut98 = 0.0; // ~98%
	//if      (trPt == 2) lasymCut98 = 1.00;
	//else if (trPt == 3) lasymCut98 = 1.60;
	//else if (trPt == 4) lasymCut98 = 1.40;
	//else if (trPt == 5) lasymCut98 = 1.40;

	//Float_t lchiyPre = fMDst->trLchiy[trPt-2] / lchiyCut98;
	//Float_t lasymPre = (2.0 * fMDst->trLasym[trPt-2] - fMDst->trLchiy[trPt-2]) / lasymCut98;
	//Float_t ccestPre = (lchiyPre>0&&lasymPre>0) ? std::sqrt(lchiyPre*lchiyPre+lasymPre*lasymPre) : std::max(lchiyPre, lasymPre);

	//Float_t ccestSft = 0.0;
	//Float_t ccestSgm = 1.0;
	//if      (trPt == 2) { ccestSft = 0.20; ccestSgm = 0.39; }
	//else if (trPt == 3) { ccestSft = 0.28; ccestSgm = 0.48; }
	//else if (trPt == 4) { ccestSft = 0.20; ccestSgm = 0.47; }
	//else if (trPt == 5) { ccestSft = 0.26; ccestSgm = 0.50; }
	//Float_t ccest = (ccestPre - ccestSft) / ccestSgm;

	//if (isPrLike) {
	//	if (sign>0) Hist::Head(StrFmt("hCOMp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
	//	if (sign<0) Hist::Head(StrFmt("hCOMn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
	//}



	/**** Low Energy Region ****/
	while (true) {
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 0, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Cutflow")->fill(uTime, nrig, 0, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Cutflow")->fill(uTime, nrig, 0, weight);

		// TRK
		if (ccest > ccestCut) break;
		
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 1, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Cutflow")->fill(uTime, nrig, 1, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Cutflow")->fill(uTime, nrig, 1, weight);

		// TRD
		if (!isTrdPr) break;

		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 2, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Cutflow")->fill(uTime, nrig, 2, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Cutflow")->fill(uTime, nrig, 2, weight);

		// ECAL
		if (fMDst->hasEcal && !isEcalPr) break;
	
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 3, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Cutflow")->fill(uTime, nrig, 3, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Cutflow")->fill(uTime, nrig, 3, weight);

		// RICH
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 4, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 4, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Cutflow")->fill(uTime, nrig, 4, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Cutflow")->fill(uTime, nrig, 4, weight);
	
		// TOF
		Float_t tofM = fMDst->tofM;
		Float_t CutTofQul = 0.3 * (TMath::Erf(6.0 * (std::sqrt(nrig) - 1.4)) + 1.0) + 0.10;
		if (std::fabs(fMDst->tofQul) > CutTofQul) return false;

		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 5, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 5, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Cutflow")->fill(uTime, nrig, 5, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Cutflow")->fill(uTime, nrig, 5, weight);

		//if (!fMDst->hasRich) {
		//	if (fMDst->richRad==0 && sign>0) Hist::Head("hLp_AglPrPh")->fill(nrig, fMDst->richPrPh, weight);
		//	if (fMDst->richRad==0 && sign<0) Hist::Head("hLn_AglPrPh")->fill(nrig, fMDst->richPrPh, weight);
		//	if (fMDst->richRad==1 && sign>0) Hist::Head("hLp_NafPrPh")->fill(nrig, fMDst->richPrPh, weight);
		//	if (fMDst->richRad==1 && sign<0) Hist::Head("hLn_NafPrPh")->fill(nrig, fMDst->richPrPh, weight);
		//	
		//	if (fMDst->richRad==0 && sign>0) Hist::Head("hLp_AglElPh")->fill(nrig, fMDst->richElPh, weight);
		//	if (fMDst->richRad==0 && sign<0) Hist::Head("hLn_AglElPh")->fill(nrig, fMDst->richElPh, weight);
		//	if (fMDst->richRad==1 && sign>0) Hist::Head("hLp_NafElPh")->fill(nrig, fMDst->richElPh, weight);
		//	if (fMDst->richRad==1 && sign<0) Hist::Head("hLn_NafElPh")->fill(nrig, fMDst->richElPh, weight);
		//}
		//Bool_t isNoRing = (!fMDst->hasRich) && 
		//                  ((fMDst->richRad==0) ? 
		//									 (fMDst->richElPh > 1.0 && fMDst->richPrPh < 0.5) :
		//									 (fMDst->richElPh > 0.5 && fMDst->richPrPh < 0.5));
		//Bool_t isNoRingCut = (isNoRing && fMDst->richPh > ((fMDst->richRad==0) ? std::min(10+Int_t(5.0*(nrig-1.0)), 15) : 5));

		//if (isNoRing) {
		//	if (fMDst->richRad==0 && sign>0) Hist::Head("hLp_AglPh")->fill(nrig, fMDst->richPh, weight);
		//	if (fMDst->richRad==0 && sign<0) Hist::Head("hLn_AglPh")->fill(nrig, fMDst->richPh, weight);
		//	if (fMDst->richRad==1 && sign>0) Hist::Head("hLp_NafPh")->fill(nrig, fMDst->richPh, weight);
		//	if (fMDst->richRad==1 && sign<0) Hist::Head("hLn_NafPh")->fill(nrig, fMDst->richPh, weight);
		//}

		//if (isNoRingCut) break;
		
		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 6, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 6, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Cutflow")->fill(uTime, nrig, 6, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Cutflow")->fill(uTime, nrig, 6, weight);

		Bool_t isAglRing = (fMDst->richRad == 0 && fMDst->hasRich);
		if (isAglRing) {
			if (sign>0) Hist::Head("hLp_AglMPr")->fill(nrig, fMDst->richMPr, weight);
			if (sign<0) Hist::Head("hLn_AglMPr")->fill(nrig, fMDst->richMPr, weight);
			if (sign>0) Hist::Head("hLp_AglMPi")->fill(nrig, fMDst->richMPi, weight);
			if (sign<0) Hist::Head("hLn_AglMPi")->fill(nrig, fMDst->richMPi, weight);
			if (sign>0) Hist::Head("hLp_AglMEl")->fill(nrig, fMDst->richMEl, weight);
			if (sign<0) Hist::Head("hLn_AglMEl")->fill(nrig, fMDst->richMEl, weight);
		}
		Bool_t isNafRing = (fMDst->richRad == 1 && fMDst->hasRich);
		if (isNafRing) {
			if (sign>0) Hist::Head("hLp_NafMPr")->fill(nrig, fMDst->richMPr, weight);
			if (sign<0) Hist::Head("hLn_NafMPr")->fill(nrig, fMDst->richMPr, weight);
			if (sign>0) Hist::Head("hLp_NafMPi")->fill(nrig, fMDst->richMPi, weight);
			if (sign<0) Hist::Head("hLn_NafMPi")->fill(nrig, fMDst->richMPi, weight);
			if (sign>0) Hist::Head("hLp_NafMEl")->fill(nrig, fMDst->richMEl, weight);
			if (sign<0) Hist::Head("hLn_NafMEl")->fill(nrig, fMDst->richMEl, weight);
		}

		
		Int_t   nrigBin = AXnr.find(nrig);
		Float_t nrigScl = (nrigBin <= AXnr.nbin()) ? AXnr.bins(nrigBin) : AXnr.max();
		Float_t massRso = 0.08;
		Float_t massScl = (1.0/(nrigScl*nrigScl+10.0*massRso) + massRso);
		Float_t tofMScl = tofM / massScl;

		Bool_t isSignal     = (sign > 0) && (!isNotPrRing);
		Bool_t isBackground = (sign < 0) && (fMDst->richRad == 0 && isRichPE && isNotPrRing);

		if (isSignal)     Hist::Head("hLs_TofM")->fill(nrig, tofM, weight);
		if (isBackground) Hist::Head("hLb_TofM")->fill(nrig, tofM, weight);
		if (IsISS && isSignal)     Hist::Head("hTLs_TofM")->fill(uTime, nrig, tofM, weight);
		if (IsISS && isBackground) Hist::Head("hTLb_TofM")->fill(uTime, nrig, tofM, weight);

		if (isSignal)     Hist::Head("hLs_TofMScl")->fill(nrig, tofMScl, weight);
		if (isBackground) Hist::Head("hLb_TofMScl")->fill(nrig, tofMScl, weight);
		if (IsISS && isSignal)     Hist::Head("hTLs_TofMScl")->fill(uTime, nrig, tofMScl, weight);
		if (IsISS && isBackground) Hist::Head("hTLb_TofMScl")->fill(uTime, nrig, tofMScl, weight);
		
		if (isNotPrRing) break;

		if (sign>0) Hist::Head("hLp_Cutflow")->fill(nrig, 7, weight);
		if (sign<0) Hist::Head("hLn_Cutflow")->fill(nrig, 7, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Cutflow")->fill(uTime, nrig, 7, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Cutflow")->fill(uTime, nrig, 7, weight);

		if (sign>0) Hist::Head("hLp_TofM")->fill(nrig, tofM, weight);
		if (sign<0) Hist::Head("hLn_TofM")->fill(nrig, tofM, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_TofM")->fill(uTime, nrig, tofM, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_TofM")->fill(uTime, nrig, tofM, weight);
		
		if (sign>0) Hist::Head("hLp_TofMScl")->fill(nrig, tofMScl, weight);
		if (sign<0) Hist::Head("hLn_TofMScl")->fill(nrig, tofMScl, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_TofMScl")->fill(uTime, nrig, tofMScl, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_TofMScl")->fill(uTime, nrig, tofMScl, weight);
			
		if (sign>0) Hist::Head("hLp_Evt")->fill(nrig, weight);
		if (sign<0) Hist::Head("hLn_Evt")->fill(nrig, weight);
		if (IsISS && sign>0) Hist::Head("hTLp_Evt")->fill(uTime, nrig, weight);
		if (IsISS && sign<0) Hist::Head("hTLn_Evt")->fill(uTime, nrig, weight);
		
		if (IsMonteCarlo) {
			if (MCSign>0) Hist::Head("hLp_MCEvt")->fill(MCNRig, weight);
			if (MCSign<0) Hist::Head("hLn_MCEvt")->fill(MCNRig, weight);
		}

		fOptL = true;
		break;	
	}	
	
	/**** Intermedia Energy Region ****/
	while (true) {
		// TRK
		if (ccest > ccestCut) break;
		
		// Version V
		while (true) {
			if (!(fMDst->hasRich && fMDst->richRad == 0)) break;
			
			if (sign>0) Hist::Head("hIp_CutflowV")->fill(nrig, 0, weight);
			if (sign<0) Hist::Head("hIn_CutflowV")->fill(nrig, 0, weight);
			if (IsISS && sign>0) Hist::Head("hTIp_CutflowV")->fill(uTime, nrig, 0, weight);
			if (IsISS && sign<0) Hist::Head("hTIn_CutflowV")->fill(uTime, nrig, 0, weight);
			
			if (sign>0) Hist::Head("hIp_CutflowV")->fill(nrig, 1, weight);
			if (sign<0) Hist::Head("hIn_CutflowV")->fill(nrig, 1, weight);
			if (IsISS && sign>0) Hist::Head("hTIp_CutflowV")->fill(uTime, nrig, 1, weight);
			if (IsISS && sign<0) Hist::Head("hTIn_CutflowV")->fill(uTime, nrig, 1, weight);

			Bool_t isSignal     = (sign > 0 && isEcalPr);
			Bool_t isBackground = (sign < 0 && isEcalEl);
			if (isSignal)     Hist::Head("hIs_TrdEstV")->fill(nrig, trdEst, weight);
			if (isBackground) Hist::Head("hIb_TrdEstV")->fill(nrig, trdEst, weight);

			if (!isRichPr) break;
			if (fMDst->hasEcal && !isEcalPr) break;
			
			if (sign>0) Hist::Head("hIp_CutflowV")->fill(nrig, 2, weight);
			if (sign<0) Hist::Head("hIn_CutflowV")->fill(nrig, 2, weight);
			if (IsISS && sign>0) Hist::Head("hTIp_CutflowV")->fill(uTime, nrig, 2, weight);
			if (IsISS && sign<0) Hist::Head("hTIn_CutflowV")->fill(uTime, nrig, 2, weight);
			
			if (sign>0) Hist::Head("hIp_TrdEstV")->fill(nrig, trdEst, weight);
			if (sign<0) Hist::Head("hIn_TrdEstV")->fill(nrig, trdEst, weight);
		
			if (sign>0) Hist::Head("hIp_EvtV")->fill(nrig, weight);
			if (sign<0) Hist::Head("hIn_EvtV")->fill(nrig, weight);
			if (IsISS && sign>0) Hist::Head("hTIp_EvtV")->fill(uTime, nrig, weight);
			if (IsISS && sign<0) Hist::Head("hTIn_EvtV")->fill(uTime, nrig, weight);

			if (IsMonteCarlo) {
				if (MCSign>0) Hist::Head("hIp_MCEvtV")->fill(MCNRig, weight);
				if (MCSign<0) Hist::Head("hIn_MCEvtV")->fill(MCNRig, weight);
			}
			
			fOptIV = true;
			break;
		}	
		
		// Version W
		while (true) {
			if (sign>0) Hist::Head("hIp_CutflowW")->fill(nrig, 0, weight);
			if (sign<0) Hist::Head("hIn_CutflowW")->fill(nrig, 0, weight);
			if (IsISS && sign>0) Hist::Head("hTIp_CutflowW")->fill(uTime, nrig, 0, weight);
			if (IsISS && sign<0) Hist::Head("hTIn_CutflowW")->fill(uTime, nrig, 0, weight);
			
			if (sign>0) Hist::Head("hIp_CutflowW")->fill(nrig, 1, weight);
			if (sign<0) Hist::Head("hIn_CutflowW")->fill(nrig, 1, weight);
			if (IsISS && sign>0) Hist::Head("hTIp_CutflowW")->fill(uTime, nrig, 1, weight);
			if (IsISS && sign<0) Hist::Head("hTIn_CutflowW")->fill(uTime, nrig, 1, weight);

			Bool_t isSignal     = (sign > 0 && isEcalPr);
			Bool_t isBackground = (sign < 0 && isEcalEl);
			if (isSignal)     Hist::Head("hIs_TrdEstW")->fill(nrig, trdEst, weight);
			if (isBackground) Hist::Head("hIb_TrdEstW")->fill(nrig, trdEst, weight);

			if (fMDst->hasRich && !isRichPr) break;
			if (fMDst->hasEcal && !isEcalPr) break;
			
			if (sign>0) Hist::Head("hIp_CutflowW")->fill(nrig, 2, weight);
			if (sign<0) Hist::Head("hIn_CutflowW")->fill(nrig, 2, weight);
			if (IsISS && sign>0) Hist::Head("hTIp_CutflowW")->fill(uTime, nrig, 2, weight);
			if (IsISS && sign<0) Hist::Head("hTIn_CutflowW")->fill(uTime, nrig, 2, weight);
			
			if (sign>0) Hist::Head("hIp_TrdEstW")->fill(nrig, trdEst, weight);
			if (sign<0) Hist::Head("hIn_TrdEstW")->fill(nrig, trdEst, weight);
		
			if (sign>0) Hist::Head("hIp_EvtW")->fill(nrig, weight);
			if (sign<0) Hist::Head("hIn_EvtW")->fill(nrig, weight);
			if (IsISS && sign>0) Hist::Head("hTIp_EvtW")->fill(uTime, nrig, weight);
			if (IsISS && sign<0) Hist::Head("hTIn_EvtW")->fill(uTime, nrig, weight);

			if (IsMonteCarlo) {
				if (MCSign>0) Hist::Head("hIp_MCEvtW")->fill(MCNRig, weight);
				if (MCSign<0) Hist::Head("hIn_MCEvtW")->fill(MCNRig, weight);
			}
			
			fOptIW = true;
			break;
		}	
		
		break;
	}
	
	/**** High Energy Region ****/
	while (true) {
		// TRK
		if (trPt == 2) break; 
		Bool_t isTrMgPt = true;
		if (trPt == 4 && nrig >= 175.00) isTrMgPt = false;
		if (trPt == 3 && nrig >= 147.00) isTrMgPt = false;

		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 0, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 0, weight);
		if (IsISS && sign>0) Hist::Head("hTHp_Cutflow")->fill(uTime, nrig, 0, weight);
		if (IsISS && sign<0) Hist::Head("hTHn_Cutflow")->fill(uTime, nrig, 0, weight);
		
		if (fMDst->hasEcal && !isEcalPr) break;

		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 1, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 1, weight);
		if (IsISS && sign>0) Hist::Head("hTHp_Cutflow")->fill(uTime, nrig, 1, weight);
		if (IsISS && sign<0) Hist::Head("hTHn_Cutflow")->fill(uTime, nrig, 1, weight);
		
		// TRD
		if (!isTrdPr) break;
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 2, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 2, weight);
		if (IsISS && sign>0) Hist::Head("hTHp_Cutflow")->fill(uTime, nrig, 2, weight);
		if (IsISS && sign<0) Hist::Head("hTHn_Cutflow")->fill(uTime, nrig, 2, weight);
		
		// TRK
		if (sign>0) Hist::Head("hHp_Cutflow")->fill(nrig, 3, weight);
		if (sign<0) Hist::Head("hHn_Cutflow")->fill(nrig, 3, weight);
		if (IsISS && sign>0) Hist::Head("hTHp_Cutflow")->fill(uTime, nrig, 3, weight);
		if (IsISS && sign<0) Hist::Head("hTHn_Cutflow")->fill(uTime, nrig, 3, weight);
	
		// CCest
		if (sign>0 && isTrMgPt) Hist::Head("hHp_CCest")->fill(nrig, ccest, weight);
		if (sign<0 && isTrMgPt) Hist::Head("hHn_CCest")->fill(nrig, ccest, weight);
		if (sign>0) Hist::Head(StrFmt("hHp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
		
		if (IsISS && sign>0 && isTrMgPt) Hist::Head("hTHp_CCest")->fill(uTime, nrig, ccest, weight);
		if (IsISS && sign<0 && isTrMgPt) Hist::Head("hTHn_CCest")->fill(uTime, nrig, ccest, weight);
		if (IsISS && sign>0) Hist::Head(StrFmt("hTHp_CCest%s", trNm.c_str()))->fill(uTime, nrig, ccest, weight);
		if (IsISS && sign<0) Hist::Head(StrFmt("hTHn_CCest%s", trNm.c_str()))->fill(uTime, nrig, ccest, weight);
		
		if (sign>0 && isTrMgPt) Hist::Head("hHp_Evt")->fill(nrig, weight);
		if (sign<0 && isTrMgPt) Hist::Head("hHn_Evt")->fill(nrig, weight);
		if (sign>0) Hist::Head(StrFmt("hHp_Evt%s", trNm.c_str()))->fill(nrig, weight);
		if (sign<0) Hist::Head(StrFmt("hHn_Evt%s", trNm.c_str()))->fill(nrig, weight);
		
		if (IsISS && sign>0 && isTrMgPt) Hist::Head("hTHp_Evt")->fill(uTime, nrig, weight);
		if (IsISS && sign<0 && isTrMgPt) Hist::Head("hTHn_Evt")->fill(uTime, nrig, weight);
		if (IsISS && sign>0) Hist::Head(StrFmt("hTHp_Evt%s", trNm.c_str()))->fill(uTime, nrig, weight);
		if (IsISS && sign<0) Hist::Head(StrFmt("hTHn_Evt%s", trNm.c_str()))->fill(uTime, nrig, weight);
		
		if (IsMonteCarlo) {
			if (MCSign>0 && isTrMgPt) Hist::Head("hHp_MCEvt")->fill(MCNRig, weight);
			if (MCSign<0 && isTrMgPt) Hist::Head("hHn_MCEvt")->fill(MCNRig, weight);
			
			if (MCSign>0) Hist::Head(StrFmt("hHp_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
			if (MCSign<0) Hist::Head(StrFmt("hHn_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
		}

		if (isTrMgPt) fOptH = true;
		break;
	}

	fMCRig = MCSign * MCNRig;
	fTRRig = sign * nrig;

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
