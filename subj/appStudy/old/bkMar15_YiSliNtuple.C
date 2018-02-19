#ifndef __YiAnaNtuple_C__
#define __YiAnaNtuple_C__

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.h"
#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.tcc"

#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jan22/src/ClassDef.h"
#include "/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jan22/src/ClassDef.C"

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

			if (runChain != nullptr) {
				UInt_t runID = 0;
				UInt_t trgEV = 0;
				UInt_t selEV = 0;
				TTree * DSTRun = new TTree("runTag", "DST Run Info");
				DSTRun->Branch("runID", &runID);
				DSTRun->Branch("trgEV", &trgEV);
				DSTRun->Branch("selEV", &selEV);
				
				RunTagInfo * runTag = new RunTagInfo;
				runChain->SetBranchAddress("runTag", &runTag);
				for (Long64_t it = 0; it < runChain->GetEntries(); ++it) {
					runChain->GetEntry(it);
					runID = runTag->runID;
					trgEV = runTag->numOfTrgEvent;
					selEV = runTag->numOfSelEvent;
					DSTRun->Fill();
				}
				delete runTag;
			}

			if (dataChain != 0 && gSaveDSTClone) fDST = dataChain->CloneTree(0);
			if (fDST == 0      && gSaveDSTTree)  {
				fDST = new TTree("data", "MinDST data");
				BranchMinDST(fDST, fMDst);
			}
			file->cd();
			
			Axis AXrso("Rigidity Resolution", 400, -1.5, 1.5);

			// rigidity binning
			AXnr = Axis("Rigidity [GV]",
				{   1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97, 
				    3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
					  8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
					 19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
					 41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
					 93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } ); // extern bin
			AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);
		
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
				const Long64_t TmeReg[2] = { 1305417600, 1480377600 };
				Long64_t ntme = Long64_t(((TmeReg[1] - TmeReg[0]) / 86400.) / 7.) + 1;
				Axis AXRTItme = Axis("Time", ntme, TmeReg[0], TmeReg[1]);
				
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
			}
			
			Axis AXCOMcut("Cutflow", 8, 0., 8.);
			Hist::New("hCOM_Cutflow", "", AXCOMcut);
			
			Axis AXCOMcutflow("Cutflow", 18, 0., 18.);
			Hist::New("hCOMi_Cutflow", "", AXir, AXCOMcutflow);
			Hist::New("hCOMp_Cutflow", "", AXnr, AXCOMcutflow);
			Hist::New("hCOMn_Cutflow", "", AXnr, AXCOMcutflow);
			
			Axis AXCOMtrInQ("Chrg", 400, 0.5, 2.6);
			Hist::New("hCOMi_TrInQ", "", AXir, AXCOMtrInQ);
			Hist::New("hCOMp_TrInQ", "", AXnr, AXCOMtrInQ);
			Hist::New("hCOMn_TrInQ", "", AXnr, AXCOMtrInQ);
			
			Axis AXCOMtrExQ("Chrg", 400, 0.5, 3.2);
			Hist::New("hCOMi_TrL2Q", "", AXir, AXCOMtrExQ);
			Hist::New("hCOMp_TrL2Q", "", AXnr, AXCOMtrExQ);
			Hist::New("hCOMn_TrL2Q", "", AXnr, AXCOMtrExQ);

			Hist::New("hCOMi_TrL1Q", "", AXir, AXCOMtrExQ);
			Hist::New("hCOMp_TrL1Q", "", AXnr, AXCOMtrExQ);
			Hist::New("hCOMn_TrL1Q", "", AXnr, AXCOMtrExQ);
			
			Hist::New("hCOMi_TrL9Q", "", AXir, AXCOMtrExQ);
			Hist::New("hCOMp_TrL9Q", "", AXnr, AXCOMtrExQ);
			Hist::New("hCOMn_TrL9Q", "", AXnr, AXCOMtrExQ);
			
			Hist::New("hCOMi_TrMxQ", "", AXir, AXCOMtrExQ);
			Hist::New("hCOMp_TrMxQ", "", AXnr, AXCOMtrExQ);
			Hist::New("hCOMn_TrMxQ", "", AXnr, AXCOMtrExQ);
	
			// 
			// dx = 0.5 * (sqrt(x*x + 4 * dA / TWO_PI) - x)
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

			Axis AXCOMtofQ("Chrg", 400, 0.0, 2.6);
			Hist::New("hCOMi_TofQ", "", AXir, AXCOMtofQ);
			Hist::New("hCOMp_TofQ", "", AXnr, AXCOMtofQ);
			Hist::New("hCOMn_TofQ", "", AXnr, AXCOMtofQ);
			
			Hist::New("hCOMi_TofQu", "", AXir, AXCOMtofQ);
			Hist::New("hCOMp_TofQu", "", AXnr, AXCOMtofQ);
			Hist::New("hCOMn_TofQu", "", AXnr, AXCOMtofQ);
			
			Hist::New("hCOMi_TofQl", "", AXir, AXCOMtofQ);
			Hist::New("hCOMp_TofQl", "", AXnr, AXCOMtofQ);
			Hist::New("hCOMn_TofQl", "", AXnr, AXCOMtofQ);
			
			Axis AXCOMdltTofQ("Delta Chrg", 400, -1.4, 1.4);
			Hist::New("hCOMi_TofQul", "", AXir, AXCOMdltTofQ);
			Hist::New("hCOMp_TofQul", "", AXnr, AXCOMdltTofQ);
			Hist::New("hCOMn_TofQul", "", AXnr, AXCOMdltTofQ);
			
			Axis AXCOMtofM("Mass Estimator", 400, -1.3, 1.3);
			Hist::New("hCOMi_TofM", "", AXir, AXCOMtofM);
			Hist::New("hCOMp_TofM", "", AXnr, AXCOMtofM);
			Hist::New("hCOMn_TofM", "", AXnr, AXCOMtofM);
			
			Axis AXCOMnacc("NAcc", 6, 0., 6.);
			Hist::New("hCOMi_NAcc", "", AXir, AXCOMnacc);
			Hist::New("hCOMp_NAcc", "", AXnr, AXCOMnacc);
			Hist::New("hCOMn_NAcc", "", AXnr, AXCOMnacc);
			
			Axis AXCOMtrdl("Trdl", 200, 0.2, 1.6);
			Hist::New("hCOMi_Trdl", "", AXir, AXCOMtrdl);
			Hist::New("hCOMp_Trdl", "", AXnr, AXCOMtrdl);
			Hist::New("hCOMn_Trdl", "", AXnr, AXCOMtrdl);
			
			Axis AXCOMtrdQ("TrdQ", 200, 0.0, 2.2);
			Hist::New("hCOMi_TrdQ", "", AXir, AXCOMtrdQ);
			Hist::New("hCOMp_TrdQ", "", AXnr, AXCOMtrdQ);
			Hist::New("hCOMn_TrdQ", "", AXnr, AXCOMtrdQ);
			
			Axis AXCOMbdt("EcalBDT", 400, -1.0, 1.0);
			Hist::New("hCOMi_EcalBDT", "", AXir, AXCOMbdt);
			Hist::New("hCOMp_EcalBDT", "", AXnr, AXCOMbdt);
			Hist::New("hCOMn_EcalBDT", "", AXnr, AXCOMbdt);
			
			Axis AXCOMaglbPr("AglbPr", 400, -0.04, 0.02);
			Hist::New("hCOMi_AglbPr", "", AXir, AXCOMaglbPr);
			Hist::New("hCOMp_AglbPr", "", AXnr, AXCOMaglbPr);
			Hist::New("hCOMn_AglbPr", "", AXnr, AXCOMaglbPr);
			
			Axis AXCOMnafbPr("NafbPr", 400, -0.06, 0.03);
			Hist::New("hCOMi_NafbPr", "", AXir, AXCOMnafbPr);
			Hist::New("hCOMp_NafbPr", "", AXnr, AXCOMnafbPr);
			Hist::New("hCOMn_NafbPr", "", AXnr, AXCOMnafbPr);
			
			Axis AXCOMaglbEl("AglbEl", 400, -0.02, 0.04);
			Hist::New("hCOMi_AglbEl", "", AXir, AXCOMaglbEl);
			Hist::New("hCOMp_AglbEl", "", AXnr, AXCOMaglbEl);
			Hist::New("hCOMn_AglbEl", "", AXnr, AXCOMaglbEl);
			
			Axis AXCOMnafbEl("NafbEl", 400, -0.03, 0.06);
			Hist::New("hCOMi_NafbEl", "", AXir, AXCOMnafbEl);
			Hist::New("hCOMp_NafbEl", "", AXnr, AXCOMnafbEl);
			Hist::New("hCOMn_NafbEl", "", AXnr, AXCOMnafbEl);
		
			Axis AXCOMaglMPr("AglMPr", 200, -0.03, 0.05);
			Hist::New("hCOMi_AglMPr", "", AXir, AXCOMaglMPr);
			Hist::New("hCOMp_AglMPr", "", AXnr, AXCOMaglMPr);
			Hist::New("hCOMn_AglMPr", "", AXnr, AXCOMaglMPr);
			
			Axis AXCOMnafMPr("NafMPr", 200, -0.06, 0.10);
			Hist::New("hCOMi_NafMPr", "", AXir, AXCOMnafMPr);
			Hist::New("hCOMp_NafMPr", "", AXnr, AXCOMnafMPr);
			Hist::New("hCOMn_NafMPr", "", AXnr, AXCOMnafMPr);
			
			Axis AXCOMaglMEl("AglMEl", 200, -0.03, 0.05);
			Hist::New("hCOMi_AglMEl", "", AXir, AXCOMaglMEl);
			Hist::New("hCOMp_AglMEl", "", AXnr, AXCOMaglMEl);
			Hist::New("hCOMn_AglMEl", "", AXnr, AXCOMaglMEl);
			
			Axis AXCOMnafMEl("NafMEl", 200, -0.06, 0.10);
			Hist::New("hCOMi_NafMEl", "", AXir, AXCOMnafMEl);
			Hist::New("hCOMp_NafMEl", "", AXnr, AXCOMnafMEl);
			Hist::New("hCOMn_NafMEl", "", AXnr, AXCOMnafMEl);
			
			Axis AXCOMlchi("lchi", 400, -7., 7.);
			
			Hist::New("hCOMi_LchixIn", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchixIn", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchixIn", "", AXnr, AXCOMlchi);
			
			Hist::New("hCOMi_LchixL1", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchixL1", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchixL1", "", AXnr, AXCOMlchi);
			
			Hist::New("hCOMi_LchixL9", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchixL9", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchixL9", "", AXnr, AXCOMlchi);
			
			Hist::New("hCOMi_LchixFs", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchixFs", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchixFs", "", AXnr, AXCOMlchi);
		
			Hist::New("hCOMi_LchiyIn", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchiyIn", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyIn", "", AXnr, AXCOMlchi);

			Hist::New("hCOMi_LchiyL1", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchiyL1", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyL1", "", AXnr, AXCOMlchi);

			Hist::New("hCOMi_LchiyL9", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchiyL9", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyL9", "", AXnr, AXCOMlchi);

			Hist::New("hCOMi_LchiyFs", "", AXir, AXCOMlchi);
			Hist::New("hCOMp_LchiyFs", "", AXnr, AXCOMlchi);
			Hist::New("hCOMn_LchiyFs", "", AXnr, AXCOMlchi);

			Axis AXCOMlasym("Lasym", 400, -6., 4.);
			
			Hist::New("hCOMi_LasymUL", "", AXir, AXCOMlasym);
			Hist::New("hCOMp_LasymUL", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_LasymUL", "", AXnr, AXCOMlasym);
			
			Hist::New("hCOMi_Lasym1I", "", AXir, AXCOMlasym);
			Hist::New("hCOMp_Lasym1I", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym1I", "", AXnr, AXCOMlasym);
			
			Hist::New("hCOMi_Lasym9I", "", AXir, AXCOMlasym);
			Hist::New("hCOMp_Lasym9I", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym9I", "", AXnr, AXCOMlasym);
			
			Hist::New("hCOMi_Lasym91", "", AXir, AXCOMlasym);
			Hist::New("hCOMp_Lasym91", "", AXnr, AXCOMlasym);
			Hist::New("hCOMn_Lasym91", "", AXnr, AXCOMlasym);

			Axis AXCOMccest("CCest", 200, -5.0, 3.0);
			
			Hist::New("hCOMi_CCestIn", "", AXir, AXCOMccest);
			Hist::New("hCOMp_CCestIn", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestIn", "", AXnr, AXCOMccest);
			
			Hist::New("hCOMi_CCestL1", "", AXir, AXCOMccest);
			Hist::New("hCOMp_CCestL1", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestL1", "", AXnr, AXCOMccest);
			
			Hist::New("hCOMi_CCestL9", "", AXir, AXCOMccest);
			Hist::New("hCOMp_CCestL9", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestL9", "", AXnr, AXCOMccest);
			
			Hist::New("hCOMi_CCestFs", "", AXir, AXCOMccest);
			Hist::New("hCOMp_CCestFs", "", AXnr, AXCOMccest);
			Hist::New("hCOMn_CCestFs", "", AXnr, AXCOMccest);
		
			Hist::New("hCOMi_EvtIn", "", AXir);
			Hist::New("hCOMp_EvtIn", "", AXnr);
			Hist::New("hCOMn_EvtIn", "", AXnr);
			
			Hist::New("hCOMi_EvtL1", "", AXir);
			Hist::New("hCOMp_EvtL1", "", AXnr);
			Hist::New("hCOMn_EvtL1", "", AXnr);
			
			Hist::New("hCOMi_EvtL9", "", AXir);
			Hist::New("hCOMp_EvtL9", "", AXnr);
			Hist::New("hCOMn_EvtL9", "", AXnr);
			
			Hist::New("hCOMi_EvtFs", "", AXir);
			Hist::New("hCOMp_EvtFs", "", AXnr);
			Hist::New("hCOMn_EvtFs", "", AXnr);

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

			std::cout << "\n<<  init DST End  >>\n";
			file->cd();
			return;
		}

		void resetDST() {
#if Debug == true
	std::cerr << "DST::resetDST()\n";
#endif
			InitMinDST(fMDst);
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
		MinDST   fMDst;
	
	public :
		static Bool_t gSaveDST;
		static Bool_t gSaveDSTTree;
		static Bool_t gSaveDSTClone;
		static UInt_t gUTime[2]; // (pre, cur)

	public :
		static Axis AXnr;
		static Axis AXir;
		static const Float_t CfStableFT = 1.0; // 1.2

	public :
		static Float_t GetStdMDR(Int_t ipatt);
		static Float_t GetStdSqrMDR(Int_t ipatt);
		static Float_t GetMDR(Int_t ipatt);
		static Float_t GetSqrMDR(Int_t ipatt);
		static Float_t GetRigSgm(Int_t ipatt, Float_t rig, Float_t mass);
};

Bool_t DST::gSaveDST = true;
Bool_t DST::gSaveDSTTree = true;
Bool_t DST::gSaveDSTClone = false;
UInt_t DST::gUTime[2] = {0, 0};

Axis DST::AXnr;
Axis DST::AXir;

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

#if Debug == true
	std::cerr << "YiAna::analyzeCore()   [][] COM Selection [][]\n";
#endif
	
	// Sample
	Hist::Head("hCOM_Cutflow")->fill(0, weight);

	// preselection (Mini-RTI Requirement)
	if (!YiAna::analyzeRunInfo()) return false;
	
	Hist::Head("hCOM_Cutflow")->fill(1, weight);
	
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

	if (isTrIn) Hist::Head("hCOM_Cutflow")->fill(4, weight);
	if (isTrL1) Hist::Head("hCOM_Cutflow")->fill(5, weight);
	if (isTrL9) Hist::Head("hCOM_Cutflow")->fill(6, weight);
	if (isTrFs) Hist::Head("hCOM_Cutflow")->fill(7, weight);

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
	Bool_t  isOverCf = true;
	Float_t cfRig = 0.0;
	Float_t cfScl = 0.0;
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
			if (MgntNum::Compare(minRig, DST::CfStableFT * fRti->cutoffIGRF[patAgl]) >= 0) isOverCf = true;
			cfRig = fRti->cutoffIGRF[patAgl];
			cfScl = minRig / fRti->cutoffIGRF[patAgl];
		}
		else isOverCf = false;
	}
	if (!isOverCf) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 0, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 0, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 0, weight);
	
	Hist::Head("hCOMi_TrInQ")->fill(irig, track->Qinner, weight);
	if (sign>0) Hist::Head("hCOMp_TrInQ")->fill(nrig, track->Qinner, weight);
	if (sign<0) Hist::Head("hCOMn_TrInQ")->fill(nrig, track->Qinner, weight);
	
	// preselection (Min-Track Charge Requirement) 
	if (track->Qinner < 0.8 || track->Qinner > 1.4) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 1, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 1, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 1, weight);
	
	// preselection (Mini-FriducialVolume Requirement)
	Hist::Head("hCOMi_TrRadius")->fill(irig, trRad, weight);
	if (sign>0) Hist::Head("hCOMp_TrRadius")->fill(nrig, trRad, weight);
	if (sign<0) Hist::Head("hCOMn_TrRadius")->fill(nrig, trRad, weight);
	
	Hist::Head("hCOMi_TrAngle")->fill(irig, trAgl, weight);
	if (sign>0) Hist::Head("hCOMp_TrAngle")->fill(nrig, trAgl, weight);
	if (sign<0) Hist::Head("hCOMn_TrAngle")->fill(nrig, trAgl, weight);

	const Float_t trFiducialVolume = 47.0;
	if (trRad > trFiducialVolume) return false;

	Hist::Head("hCOMi_Cutflow")->fill(irig, 2, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 2, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 2, weight);
	
	// preselection (Mini-Trigger Requirement)
	if ((fTrg->bit&8) != 8) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 3, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 3, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 3, weight);
	
	// preselection (Mini-TOF Requirement)
	if (!fTof->statusBetaH) return false;
	if (fTof->betaHPatt != 15) return false;
	if (fTof->betaH < 0.3 || fTof->betaH > 1.3) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 4, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 4, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 4, weight);
	
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
	
	if (trHitMx) Hist::Head("hCOMi_TrMxQ")->fill(irig, trMxQ, weight);
	if (trHitMx && sign>0) Hist::Head("hCOMp_TrMxQ")->fill(nrig, trMxQ, weight);
	if (trHitMx && sign<0) Hist::Head("hCOMn_TrMxQ")->fill(nrig, trMxQ, weight);
	
	if (trHitL2) Hist::Head("hCOMi_TrL2Q")->fill(irig, trL2Q, weight);
	if (trHitL2 && sign>0) Hist::Head("hCOMp_TrL2Q")->fill(nrig, trL2Q, weight);
	if (trHitL2 && sign<0) Hist::Head("hCOMn_TrL2Q")->fill(nrig, trL2Q, weight);

	//if (trHitMx && (trMxQ < 0.90 || trMxQ > 2.3)) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 5, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 5, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 5, weight);
	
	//if (trHitL2 && (trL2Q < 0.75 || trL2Q > 2.1)) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 6, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 6, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 6, weight);
	
	if (trHitL1) Hist::Head("hCOMi_TrL1Q")->fill(irig, trL1Q, weight);
	if (trHitL1 && sign>0) Hist::Head("hCOMp_TrL1Q")->fill(nrig, trL1Q, weight);
	if (trHitL1 && sign<0) Hist::Head("hCOMn_TrL1Q")->fill(nrig, trL1Q, weight);

	if (trHitL9) Hist::Head("hCOMi_TrL9Q")->fill(irig, trL9Q, weight);
	if (trHitL9 && sign>0) Hist::Head("hCOMp_TrL9Q")->fill(nrig, trL9Q, weight);
	if (trHitL9 && sign<0) Hist::Head("hCOMn_TrL9Q")->fill(nrig, trL9Q, weight);
	
	//if (trHitL1 && (trL1Q < 0.8 || trL1Q > 2.1)) return false;
	//if (trHitL9 && (trL9Q < 0.8 || trL9Q > 2.1)) return false;

	Hist::Head("hCOMi_Cutflow")->fill(irig, 7, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 7, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 7, weight);
	
	// preselection (Mini-TOF Requirement)
	Float_t tofb   = (fTof->betaH - tbtaPr);
	Float_t tofM   = ((1.0/fTof->betaH/fTof->betaH - 1.0) - (massPr/nrig)*(massPr/nrig));
	Float_t tofQ   = fTof->Qall;
	Float_t tofQu  = 0.50 * (fTof->Q[0] + fTof->Q[1]);
	Float_t tofQl  = 0.50 * (fTof->Q[2] + fTof->Q[3]);
	Float_t tofQul = tofQu - tofQl;

	if (fTof->normChisqT > 10) return false;
	if (fTof->normChisqC > 10) return false;

	Hist::Head("hCOMi_Cutflow")->fill(irig, 8, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 8, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 8, weight);
	
	if (fTof->numOfInTimeCluster > 4) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 9, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 9, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 9, weight);
	
	Hist::Head("hCOMi_TofQ")->fill(irig, tofQ, weight);
	if (sign>0) Hist::Head("hCOMp_TofQ")->fill(nrig, tofQ, weight);
	if (sign<0) Hist::Head("hCOMn_TofQ")->fill(nrig, tofQ, weight);

	if (tofQ  < 0.80 || tofQ  > 1.60) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 10, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 10, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 10, weight);
	
	Hist::Head("hCOMi_TofQu")->fill(irig, tofQu, weight);
	if (sign>0) Hist::Head("hCOMp_TofQu")->fill(nrig, tofQu, weight);
	if (sign<0) Hist::Head("hCOMn_TofQu")->fill(nrig, tofQu, weight);
	
	Hist::Head("hCOMi_TofQl")->fill(irig, tofQl, weight);
	if (sign>0) Hist::Head("hCOMp_TofQl")->fill(nrig, tofQl, weight);
	if (sign<0) Hist::Head("hCOMn_TofQl")->fill(nrig, tofQl, weight);
	
	Hist::Head("hCOMi_TofM")->fill(irig, tofM, weight);
	if (sign>0) Hist::Head("hCOMp_TofM")->fill(nrig, tofM, weight);
	if (sign<0) Hist::Head("hCOMn_TofM")->fill(nrig, tofM, weight);
	
	if (tofQu < 0.75 || tofQu > 1.8) return false;
	if (tofQl < 0.75 || tofQl > 1.8) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 11, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 11, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 11, weight);
	
	Float_t CutTofQul = 0.30 * (TMath::Erf(1.5 * (std::sqrt(nrig) - 1.8)) + 1.0) + 0.10;
	
	Hist::Head("hCOMi_TofQul")->fill(irig, tofQul, weight);
	if (sign>0) Hist::Head("hCOMp_TofQul")->fill(nrig, tofQul, weight);
	if (sign<0) Hist::Head("hCOMn_TofQul")->fill(nrig, tofQul, weight);

	//if (std::fabs(tofQul) > CutTofQul) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 12, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 12, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 12, weight);
	
	// preselection (Mini-Acc Requirement)
	Hist::Head("hCOMi_NAcc")->fill(irig, fAcc->numOfCluster, weight);
	if (sign>0) Hist::Head("hCOMp_NAcc")->fill(nrig, fAcc->numOfCluster, weight);
	if (sign<0) Hist::Head("hCOMn_NAcc")->fill(nrig, fAcc->numOfCluster, weight);
	
	Int_t nAccCut = 4;
	if (fAcc->numOfCluster > nAccCut) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 13, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 13, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 13, weight);
	
	// preselection (Mini-TRD Requirement)
	Bool_t isCleanTrdTrack = 
		(fTrd->numOfTrack <= 2 && fTrd->numOfHTrack <= 2) && 
		(fTrd->numOfTrack == 1 || fTrd->numOfHTrack == 1);
	if (!isCleanTrdTrack) return false;
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 14, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 14, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 14, weight);
	
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
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 15, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 15, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 15, weight);
	
	Hist::Head("hCOMi_TrdQ")->fill(irig, trdQ, weight);
	if (sign>0) Hist::Head("hCOMp_TrdQ")->fill(nrig, trdQ, weight);
	if (sign<0) Hist::Head("hCOMn_TrdQ")->fill(nrig, trdQ, weight);
	
	if (istrdHe || istrdSQ) return false;

	Hist::Head("hCOMi_Trdl")->fill(irig, trdl, weight);
	if (sign>0) Hist::Head("hCOMp_Trdl")->fill(nrig, trdl, weight);
	if (sign<0) Hist::Head("hCOMn_Trdl")->fill(nrig, trdl, weight);
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 16, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 16, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 16, weight);

	// preselection (Mini-Shower Requirement)
	ShowerInfo * shower = (fEcal->showers.size() >= 1) ? (&fEcal->showers.at(0)) : nullptr;
	Bool_t  hasEcal  = (shower != nullptr);
	Float_t ecalbdt  = (!hasEcal) ? -2.0  : shower->PisaBDT;
	Bool_t  isecalPr = (!hasEcal) ? false : (ecalbdt < -0.8);
	Bool_t  isecalEl = (!hasEcal) ? false : (ecalbdt >  0.6);
	
	if (hasEcal) Hist::Head("hCOMi_EcalBDT")->fill(irig, ecalbdt, weight);
	if (hasEcal && sign>0) Hist::Head("hCOMp_EcalBDT")->fill(nrig, ecalbdt, weight);
	if (hasEcal && sign<0) Hist::Head("hCOMn_EcalBDT")->fill(nrig, ecalbdt, weight);

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
	Bool_t  isrichNo  = (fRich->kindOfRad == -1) ? false : (!fRich->status && fRich->kindOfRad == 0 && MgntNum::Compare(richelph, 2.0f) > 0 && MgntNum::Compare(richprph, 0.5f) < 0 && richph == 0);
	Bool_t  isrichPr  = (fRich->kindOfRad == -1) ? false : ( fRich->status && std::fabs(richbPr) < (richbTh + richbSh));
	Bool_t  isrichPi  = (fRich->kindOfRad == -1) ? false : ( fRich->status && std::fabs(richbPi) < (richbTh));
	Bool_t  isrichEl  = (fRich->kindOfRad == -1) ? false : ( fRich->status && std::fabs(richbEl) < (richbTh));

	if (hasRich && fRich->kindOfRad == 0) Hist::Head("hCOMi_AglbPr")->fill(irig, richbPr, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head("hCOMp_AglbPr")->fill(nrig, richbPr, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head("hCOMn_AglbPr")->fill(nrig, richbPr, weight);
	
	if (hasRich && fRich->kindOfRad == 1) Hist::Head("hCOMi_NafbPr")->fill(irig, richbPr, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head("hCOMp_NafbPr")->fill(nrig, richbPr, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head("hCOMn_NafbPr")->fill(nrig, richbPr, weight);
	
	if (hasRich && fRich->kindOfRad == 0) Hist::Head("hCOMi_AglbEl")->fill(irig, richbEl, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head("hCOMp_AglbEl")->fill(nrig, richbEl, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head("hCOMn_AglbEl")->fill(nrig, richbEl, weight);

	if (hasRich && fRich->kindOfRad == 1) Hist::Head("hCOMi_NafbEl")->fill(irig, richbEl, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head("hCOMp_NafbEl")->fill(nrig, richbEl, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head("hCOMn_NafbEl")->fill(nrig, richbEl, weight);
	
	if (hasRich && fRich->kindOfRad == 0) Hist::Head("hCOMi_AglMPr")->fill(irig, richMPr, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head("hCOMp_AglMPr")->fill(nrig, richMPr, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head("hCOMn_AglMPr")->fill(nrig, richMPr, weight);
	
	if (hasRich && fRich->kindOfRad == 1) Hist::Head("hCOMi_NafMPr")->fill(irig, richMPr, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head("hCOMp_NafMPr")->fill(nrig, richMPr, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head("hCOMn_NafMPr")->fill(nrig, richMPr, weight);
	
	if (hasRich && fRich->kindOfRad == 0) Hist::Head("hCOMi_AglMEl")->fill(irig, richMEl, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign>0) Hist::Head("hCOMp_AglMEl")->fill(nrig, richMEl, weight);
	if (hasRich && fRich->kindOfRad == 0 && sign<0) Hist::Head("hCOMn_AglMEl")->fill(nrig, richMEl, weight);

	if (hasRich && fRich->kindOfRad == 1) Hist::Head("hCOMi_NafMEl")->fill(irig, richMEl, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign>0) Hist::Head("hCOMp_NafMEl")->fill(nrig, richMEl, weight);
	if (hasRich && fRich->kindOfRad == 1 && sign<0) Hist::Head("hCOMn_NafMEl")->fill(nrig, richMEl, weight);
	
	// Charge Confusion
	Float_t tuneCCPnt = 0.0;
	if      (trPt == 2) tuneCCPnt = 40.0;
	else if (trPt == 3) tuneCCPnt = 40.0;
	else if (trPt == 4) tuneCCPnt = 60.0;
	else if (trPt == 5) tuneCCPnt = 75.0;
	
	// Charge Confusion ---- Chisq
	Float_t lchix    = std::log(track->chisq[0][trPt][0]);
	Float_t lchiy    = std::log(track->chisq[0][trPt][1]);
	Float_t lchiWgt  = GetSqrMDR(trPt);
	Float_t lchixCut   = 3.0; // ~98%
	Float_t lchiyCutTh[4] = { 3.6, 3.4, 2.7, 2.6 }; // ~99%
	Float_t lchiyCut98 = 0.0; // ~98%
	Float_t lchiyCut95 = 0.0; // ~95%
	if      (trPt == 2) { lchiyCut98 = 2.35; lchiyCut95 = 1.95; } 
	else if (trPt == 3) { lchiyCut98 = 2.20; lchiyCut95 = 1.80; }
	else if (trPt == 4) { lchiyCut98 = 2.00; lchiyCut95 = 1.65; }
	else if (trPt == 5) { lchiyCut98 = 2.20; lchiyCut95 = 1.75; }
	
	Float_t lchiyIn = (track->status[0][2]) ? std::log(track->chisq[0][2][1]) : 0.0;
	Float_t lchiyL1 = (track->status[0][3]) ? std::log(track->chisq[0][3][1]) : 0.0;
	Float_t lchiyL9 = (track->status[0][4]) ? std::log(track->chisq[0][4][1]) : 0.0;
	Float_t lchiyFs = (track->status[0][5]) ? std::log(track->chisq[0][5][1]) : 0.0;

	// Charge Confusion ---- AsymRig
	std::string lasymNm    = "";
	Float_t     lasym      = 0.0;
	Float_t     lasymWgt   = 0.0;
	Float_t     lasymCutTh[4] = { 1.75, 1.90, 1.60, 1.90 }; // ~99%
	Float_t     lasymCut98 = 0.0; // ~98%
	Float_t     lasymCut95 = 0.0; // ~95%
	if      (trPt == 2) { lasymCut98 = 1.25; lasymCut95 = 0.95; }
	else if (trPt == 3) { lasymCut98 = 1.55; lasymCut95 = 1.20; }
	else if (trPt == 4) { lasymCut98 = 1.35; lasymCut95 = 1.05; }
	else if (trPt == 5) { lasymCut98 = 1.70; lasymCut95 = 1.30; }

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
	Float_t lasymUL     = (track->status[0][0] && track->status[0][1]) ? std::log(asymUL * asymUL + lasymLMT) / lasymSGM : 0.0;
	if (trPt == 2) {
		lasymNm  = "UL";
		lasym    = lasymUL;
		lasymWgt = GetStdSqrMDR(2);
	}

	Float_t crrPar1I[5] = { 4.93251e-01, 1.90417e+00, 9.18502e+01, 4.86389e+00, 3.29747e-01 };
	Float_t crrSgm1I    = crrPar1I[0]*std::erf(crrPar1I[1]*(std::log(nrig+crrPar1I[2])-crrPar1I[3]))+crrPar1I[4];
	Float_t trSgm1I     = std::sqrt(trSgmIn * trSgmIn + trSgmL1 * trSgmL1) * crrSgm1I;
	Float_t asym1I      = (1.0/track->rigidity[0][2] - 1.0/track->rigidity[0][3]) / trSgm1I; 
	Float_t lasym1I     = (track->status[0][2] && track->status[0][3]) ? std::log(asym1I * asym1I + lasymLMT) / lasymSGM : 0.0;
	if (trPt == 3) {
		lasymNm  = "1I";
		lasym    = lasym1I;
		lasymWgt = GetStdSqrMDR(3);
	}
	
	Float_t crrPar9I[5] = { 4.17294e-01, 1.15603e+00, 1.62194e+01, 3.81780e+00, 3.92404e-01 };
	Float_t crrSgm9I    = crrPar9I[0]*std::erf(crrPar9I[1]*(std::log(nrig+crrPar9I[2])-crrPar9I[3]))+crrPar9I[4];
	Float_t trSgm9I     = std::sqrt(trSgmIn * trSgmIn + trSgmL9 * trSgmL9) * crrSgm9I;
	Float_t asym9I      = (1.0/track->rigidity[0][2] - 1.0/track->rigidity[0][4]) / trSgm9I;
	Float_t lasym9I     = (track->status[0][2] && track->status[0][4]) ? std::log(asym9I * asym9I + lasymLMT) / lasymSGM : 0.0;
	if (trPt == 4) {
		lasymNm  = "9I";
		lasym    = lasym9I;
		lasymWgt = GetStdSqrMDR(4);
	}

	Float_t crrPar91[5] = { 5.25352e-01, 9.88059e-01, 1.21200e+01, 3.72027e+00, 4.99424e-01 };
	Float_t crrSgm91    = crrPar91[0]*std::erf(crrPar91[1]*(std::log(nrig+crrPar91[2])-crrPar91[3]))+crrPar91[4];
	Float_t trSgm91     = std::sqrt(trSgmL1 * trSgmL1 + trSgmL9 * trSgmL9) * crrSgm91;
	Float_t asym91      = (1.0/track->rigidity[0][3] - 1.0/track->rigidity[0][4]) / trSgm91;
	Float_t lasym91     = (track->status[0][3] && track->status[0][4]) ? std::log(asym91 * asym91 + lasymLMT) / lasymSGM : 0.0;
	if (trPt == 5) {
		lasymNm  = "91";
		lasym    = lasym91;
		lasymWgt = GetStdSqrMDR(5);
	}
	
	// Charge Confusion ---- CC Estimator
	Float_t ccest    = ((lchiy * lchiWgt + lasym * lasymWgt) / (lchiWgt + lasymWgt));

	Bool_t ccIn = (lchiyIn > lchiyCutTh[0] || lasymUL > lasymCutTh[0] ||
	               (lchiyIn > 0 && lasymUL > 0 && (lchiyIn*lchiyIn/lchiyCutTh[0]/lchiyCutTh[0] + lasymUL*lasymUL/lasymCutTh[0]/lasymCutTh[0]) > 1) );
	Bool_t ccL1 = (lchiyL1 > lchiyCutTh[1] || lasym1I > lasymCutTh[1] ||
	               (lchiyL1 > 0 && lasym1I > 0 && (lchiyL1*lchiyL1/lchiyCutTh[1]/lchiyCutTh[1] + lasym1I*lasym1I/lasymCutTh[1]/lasymCutTh[1]) > 1) );
	Bool_t ccL9 = (lchiyL9 > lchiyCutTh[2] || lasym9I > lasymCutTh[2] ||
	               (lchiyL9 > 0 && lasym9I > 0 && (lchiyL9*lchiyL9/lchiyCutTh[2]/lchiyCutTh[2] + lasym9I*lasym9I/lasymCutTh[2]/lasymCutTh[2]) > 1) );
	Bool_t ccFs = (lchiyFs > lchiyCutTh[3] || lasym91 > lasymCutTh[3] ||
	               (lchiyFs > 0 && lasym91 > 0 && (lchiyFs*lchiyFs/lchiyCutTh[3]/lchiyCutTh[3] + lasym91*lasym91/lasymCutTh[3]/lasymCutTh[3]) > 1) );

	//If (trPt == 3 && ccIn) return false;
	//If (trPt == 4 && ccIn) return false;
	//If (trPt == 5 && (ccL1 || ccL9)) return false;
	
	Bool_t isPrLike = (istrdPr && (hasEcal ? isecalPr : true) && (hasRich ? isrichPr : true));
	
	Hist::Head("hCOMi_Cutflow")->fill(irig, 17, weight);
	if (sign>0) Hist::Head("hCOMp_Cutflow")->fill(nrig, 17, weight);
	if (sign<0) Hist::Head("hCOMn_Cutflow")->fill(nrig, 17, weight);
	
	Hist::Head(StrFmt("hCOMi_Lchix%s", trNm.c_str()))->fill(irig, lchix, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Lchix%s", trNm.c_str()))->fill(nrig, lchix, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Lchix%s", trNm.c_str()))->fill(nrig, lchix, weight);
	
	Hist::Head(StrFmt("hCOMi_Lchiy%s", trNm.c_str()))->fill(irig, lchiy, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Lchiy%s", trNm.c_str()))->fill(nrig, lchiy, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Lchiy%s", trNm.c_str()))->fill(nrig, lchiy, weight);
	
	Hist::Head(StrFmt("hCOMi_Lasym%s", lasymNm.c_str()))->fill(irig, lasym, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Lasym%s", lasymNm.c_str()))->fill(nrig, lasym, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Lasym%s", lasymNm.c_str()))->fill(nrig, lasym, weight);
	
	Hist::Head(StrFmt("hCOMi_CCest%s", trNm.c_str()))->fill(irig, ccest, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
	
	Hist::Head(StrFmt("hCOMi_Evt%s", trNm.c_str()))->fill(irig, weight);
	if (sign>0) Hist::Head(StrFmt("hCOMp_Evt%s", trNm.c_str()))->fill(nrig, weight);
	if (sign<0) Hist::Head(StrFmt("hCOMn_Evt%s", trNm.c_str()))->fill(nrig, weight);

	if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
		Hist::Head(StrFmt("hCOMi_MCTrRso%s", trNm.c_str()))->fill(MCIRig, irig, weight);
		if (MCSign >0) Hist::Head(StrFmt("hCOMp_MCTrRso%s", trNm.c_str()))->fill(MCNRig, MCNRigCen * (irig - MCIRig), weight);
		if (MCSign <0) Hist::Head(StrFmt("hCOMn_MCTrRso%s", trNm.c_str()))->fill(MCNRig, MCNRigCen * (irig - MCIRig), weight);
	
		Hist::Head(StrFmt("hCOMi_MCEvt%s", trNm.c_str()))->fill(MCIRig, weight);
		if (MCSign >0) Hist::Head(StrFmt("hCOMp_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
		if (MCSign <0) Hist::Head(StrFmt("hCOMn_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
	}

	// MinDST
	fMDst.run    = runID;
	fMDst.event  = eventID;
	fMDst.weight = weight;

	if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
	fMDst.uTime    = fRti->uTime;
	fMDst.UTCyr    = fRti->UTCyr;
	fMDst.UTCyd    = fRti->UTCyd;
	fMDst.UTChr    = fRti->UTChr;
	fMDst.cfRig    = cfRig;
	fMDst.cfScl    = cfScl;
	fMDst.lvtme    = fRti->liveTime;
	}

	if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
	fMDst.mcSign = MCSign;
	fMDst.mcNRig = MCNRig;
	}

	fMDst.trPt    = trPt;
	fMDst.trSign  = sign;
	fMDst.trNRig  = nrig;
	
	fMDst.trLchix = lchix;
	fMDst.trLchiy[0]  = lchiyIn;
	fMDst.trLchiy[1]  = lchiyL1;
	fMDst.trLchiy[2]  = lchiyL9;
	fMDst.trLchiy[3]  = lchiyFs;
	fMDst.trLasym[0]  = lasymUL;
	fMDst.trLasym[1]  = lasym1I;
	fMDst.trLasym[2]  = lasym9I;
	fMDst.trLasym[3]  = lasym91;
	fMDst.trCCest = ccest;

	//fMDst.trLchix = lchix;
	//fMDst.trLchiy = lchiy;
	//fMDst.trLasym = lasym;
	//fMDst.trCCest = ccest;
	fMDst.trMxQ   = trMxQ;
	fMDst.trL2Q   = trL2Q;
	fMDst.trL1Q   = trL1Q;
	fMDst.trL9Q   = trL9Q;

	fMDst.tofQul = tofQul;
	//fMDst.tofb   = tofb;
	fMDst.tofM   = tofM;
	
	fMDst.trdEst = trdl;
	
	if (hasEcal) {
	fMDst.hasEcal = hasEcal;
	fMDst.ecalEst = ecalbdt;
	}

	fMDst.richRad  = fRich->kindOfRad;
	fMDst.richPh   = richph;
	fMDst.richPrPh = richprph;
	fMDst.richElPh = richelph;
	
	if (hasRich) {
	fMDst.hasRich  = hasRich;
	//fMDst.richbPr  = richbPr;
	//fMDst.richbPi  = richbPi;
	//fMDst.richbEl  = richbEl;
	fMDst.richMPr  = richMPr;
	fMDst.richMPi  = richMPi;
	fMDst.richMEl  = richMEl;
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
