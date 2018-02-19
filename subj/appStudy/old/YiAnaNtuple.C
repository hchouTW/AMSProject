#ifndef __YiAnaNtuple_C__
#define __YiAnaNtuple_C__

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.h"
#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/template/YiAnaNtuple.tcc"

#include "/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/src/minDST.h"

using namespace MGROOT;
#include <TGraphAsymmErrors.h>
#include <TF1.h>


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
                fDST->Branch("OptI" ,  &fOptI);
                fDST->Branch("OptM" ,  &fOptM);
                fDST->Branch("OptHL1", &fOptHL1);
                fDST->Branch("OptHL9", &fOptHL9);
                fDST->Branch("OptHFs", &fOptHFs);
            }
            file->cd();

            Axis AXrso("Rigidity Resolution", 400, -1.5, 1.5);

            DST::UnfoldRatio->SetParameters(-6.29727e+02, 3.49266e-02, 3.26913e+01, -5.17749e+01);

            // rigidity binning
            AXnr = Axis("|Rigidity| [GV]", BinList(
                        {   1.00,   1.16,   1.33,   1.51,   1.71,   1.92,   2.15,   2.40,   2.67,   2.97,
                        3.29,   3.64,   4.02,   4.43,   4.88,   5.37,   5.90,   6.47,   7.09,   7.76,
                        8.48,   9.26,  10.10,  11.00,  12.00,  13.00,  14.10,  15.30,  16.60,  18.00, 
                        19.50,  21.10,  22.80,  24.70,  26.70,  28.80,  31.10,  33.50,  36.10,  38.90, 
                        41.90,  45.10,  48.50,  52.20,  56.10,  60.30,  64.80,  69.70,  74.90,  80.50, 
                        93.00, 108.00, 125.00, 147.00, 175.00, 211.00, 259.00, 450.00 } ) );
            AXir = Axis("1/Rigidity [1/GV]", AXnr, 1, true);
            //DST::WithPub = true;

            //AXnr = Axis("|Rigidity| [GV]", BinList(
            //	{   1.00,   2.00,   3.00,   4.12,   5.00,   
            //	    6.00,   7.10,   8.30,   9.62,  11.04,  
            //		12.59,  14.25,  16.05,  17.98,  20.04,  
            //		22.25,  24.62,  27.25,  30.21  } ) );
            //AXir = Axis::Invert("1/Rigidity [1/GV]", AXnr);

            // 3 month
            AXtme = Axis("Time", BinList( 
                        { 1305417600, 1312416000, 1319414400, 1326412800, 1333411200,
                        1340409600, 1347408000, 1354406400, 1361404800, 1368403200,
                        1375401600, 1382400000, 1389398400, 1396396800, 1403395200,
                        1410393600, 1417392000, 1424390400, 1431388800, 1438387200,
                        1445385600, 1452384000, 1459382400, 1466380800, 1473379200, 
                        1480377600 } ) );

            //----  Low Energy  ----//
            Axis AXLlchi("lchi", 200, -6., 8.);
            Hist::New("hLp_LchixIn", "", HistAxis(AXnr, AXLlchi));
            Hist::New("hLn_LchixIn", "", HistAxis(AXnr, AXLlchi));
            Hist::New("hLp_LchiyIn", "", HistAxis(AXnr, AXLlchi));
            Hist::New("hLn_LchiyIn", "", HistAxis(AXnr, AXLlchi));

            Hist::New("hLp_LchixL1", "", HistAxis(AXnr, AXLlchi));
            Hist::New("hLn_LchixL1", "", HistAxis(AXnr, AXLlchi));
            Hist::New("hLp_LchiyL1", "", HistAxis(AXnr, AXLlchi));
            Hist::New("hLn_LchiyL1", "", HistAxis(AXnr, AXLlchi));

            Axis AXLlasym("Lasym", 200, -2.5, 5.0);
            Hist::New("hLp_Lasym1I", "", HistAxis(AXnr, AXLlasym));
            Hist::New("hLn_Lasym1I", "", HistAxis(AXnr, AXLlasym));

            Axis AXLccest("CCest", 200, -4, 6);
            Hist::New("hLp_CCestL1", "", HistAxis(AXnr, AXLccest));
            Hist::New("hLn_CCestL1", "", HistAxis(AXnr, AXLccest));

            Axis AXLtrest("TrEst", 200, -4, 4);
            Hist::New("hLp_TrEstL1", "", HistAxis(AXnr, AXLtrest));
            Hist::New("hLn_TrEstL1", "", HistAxis(AXnr, AXLtrest));

            Axis AXLaglM("Mass Estimator", 100, -6., 6.);
            Hist::New("hLp_AglM", "", HistAxis(AXnr, AXLaglM));
            Hist::New("hLn_AglM", "", HistAxis(AXnr, AXLaglM));

            Axis AXLnafM("Mass Estimator", 100, -8., 8.);
            Hist::New("hLp_NafM", "", HistAxis(AXnr, AXLnafM));
            Hist::New("hLn_NafM", "", HistAxis(AXnr, AXLnafM));

            Axis AXLtrdEst("TRD Estimator", 200, 0.2, 1.6);
            Hist::New("hLp_TrdEst", "", HistAxis(AXnr, AXLtrdEst));
            Hist::New("hLn_TrdEst", "", HistAxis(AXnr, AXLtrdEst));

            Axis AXLtofM("Mass Estimator", 100, -4.5, 4.5);
            Hist::New("hLs_TofM", "", HistAxis(AXnr, AXLtofM));
            Hist::New("hLb_TofM", "", HistAxis(AXnr, AXLtofM));

            Hist::New("hLp_TofM", "", HistAxis(AXnr, AXLtofM));
            Hist::New("hLn_TofM", "", HistAxis(AXnr, AXLtofM));

            for (Int_t icf = 12; icf <= 14; ++icf) {
                Hist::New(STR_FMT("hLs_TofM_CF%02d", icf), "", HistAxis(AXnr, AXLtofM));
                Hist::New(STR_FMT("hLb_TofM_CF%02d", icf), "", HistAxis(AXnr, AXLtofM));
                Hist::New(STR_FMT("hLp_TofM_CF%02d", icf), "", HistAxis(AXnr, AXLtofM));
                Hist::New(STR_FMT("hLn_TofM_CF%02d", icf), "", HistAxis(AXnr, AXLtofM));
            }

            for (Int_t irs = 0; irs <= 1; ++irs) {
                Hist::New(STR_FMT("hLs_TofM_RS%d", irs), "", HistAxis(AXnr, AXLtofM));
                Hist::New(STR_FMT("hLb_TofM_RS%d", irs), "", HistAxis(AXnr, AXLtofM));
                Hist::New(STR_FMT("hLp_TofM_RS%d", irs), "", HistAxis(AXnr, AXLtofM));
                Hist::New(STR_FMT("hLn_TofM_RS%d", irs), "", HistAxis(AXnr, AXLtofM));
            }

            Axis AXLcutflow("Cutflow", 4, 0., 4.);
            Hist::New("hLp_Cutflow", "", HistAxis(AXnr, AXLcutflow));
            Hist::New("hLn_Cutflow", "", HistAxis(AXnr, AXLcutflow));

            Hist::New("hLp_Evt", "", HistAxis(AXnr));
            Hist::New("hLn_Evt", "", HistAxis(AXnr));

            if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
                Hist::New("hTLs_TofM", "", HistAxis(AXtme, AXnr, AXLtofM));
                Hist::New("hTLb_TofM", "", HistAxis(AXtme, AXnr, AXLtofM));

                Hist::New("hTLp_TofM", "", HistAxis(AXtme, AXnr, AXLtofM));
                Hist::New("hTLn_TofM", "", HistAxis(AXtme, AXnr, AXLtofM));

                Hist::New("hTLp_Cutflow", "", HistAxis(AXtme, AXnr, AXLcutflow));
                Hist::New("hTLn_Cutflow", "", HistAxis(AXtme, AXnr, AXLcutflow));

                Hist::New("hTLp_Evt", "", HistAxis(AXtme, AXnr));
                Hist::New("hTLn_Evt", "", HistAxis(AXtme, AXnr));
            }

            if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
                Hist::New("hLp_MCEvt", "", HistAxis(AXnr));
                Hist::New("hLn_MCEvt", "", HistAxis(AXnr));

                Hist::New("hLp_Evt_RWGT", "", HistAxis(AXnr));
                Hist::New("hLn_Evt_RWGT", "", HistAxis(AXnr));
                Hist::New("hLp_MCEvt_RWGT", "", HistAxis(AXnr));
                Hist::New("hLn_MCEvt_RWGT", "", HistAxis(AXnr));
            }


            //----  Intermedia Energy  ----//
            Axis AXIlchi("lchi", 200, -6., 8.);
            Hist::New("hIp_LchixIn", "", HistAxis(AXnr, AXIlchi));
            Hist::New("hIn_LchixIn", "", HistAxis(AXnr, AXIlchi));
            Hist::New("hIp_LchiyIn", "", HistAxis(AXnr, AXIlchi));
            Hist::New("hIn_LchiyIn", "", HistAxis(AXnr, AXIlchi));

            Hist::New("hIp_LchixL1", "", HistAxis(AXnr, AXIlchi));
            Hist::New("hIn_LchixL1", "", HistAxis(AXnr, AXIlchi));
            Hist::New("hIp_LchiyL1", "", HistAxis(AXnr, AXIlchi));
            Hist::New("hIn_LchiyL1", "", HistAxis(AXnr, AXIlchi));

            Axis AXIlasym("Lasym", 200, -2.5, 5.0);
            Hist::New("hIp_Lasym1I", "", HistAxis(AXnr, AXIlasym));
            Hist::New("hIn_Lasym1I", "", HistAxis(AXnr, AXIlasym));

            Axis AXIccest("CCest", 200, -4, 6);
            Hist::New("hIp_CCestL1", "", HistAxis(AXnr, AXIccest));
            Hist::New("hIn_CCestL1", "", HistAxis(AXnr, AXIccest));

            Axis AXItrest("TrEst", 200, -4, 4);
            Hist::New("hIp_TrEstL1", "", HistAxis(AXnr, AXItrest));
            Hist::New("hIn_TrEstL1", "", HistAxis(AXnr, AXItrest));

            Axis AXItrdEst("TRD Estimator", 100, 0.2, 1.6);
            Hist::New("hIs_TrdEst", "", HistAxis(AXnr, AXItrdEst));
            Hist::New("hIb_TrdEst", "", HistAxis(AXnr, AXItrdEst));

            Hist::New("hIp_TrdEst", "", HistAxis(AXnr, AXItrdEst));
            Hist::New("hIn_TrdEst", "", HistAxis(AXnr, AXItrdEst));

            for (Int_t icf = 12; icf <= 14; ++icf) {
                Hist::New(STR_FMT("hIs_TrdEst_CF%02d", icf), "", HistAxis(AXnr, AXItrdEst));
                Hist::New(STR_FMT("hIb_TrdEst_CF%02d", icf), "", HistAxis(AXnr, AXItrdEst));
                Hist::New(STR_FMT("hIp_TrdEst_CF%02d", icf), "", HistAxis(AXnr, AXItrdEst));
                Hist::New(STR_FMT("hIn_TrdEst_CF%02d", icf), "", HistAxis(AXnr, AXItrdEst));
            }

            for (Int_t irs = 0; irs <= 1; ++irs) {
                Hist::New(STR_FMT("hIs_TrdEst_RS%d", irs), "", HistAxis(AXnr, AXItrdEst));
                Hist::New(STR_FMT("hIb_TrdEst_RS%d", irs), "", HistAxis(AXnr, AXItrdEst));
                Hist::New(STR_FMT("hIp_TrdEst_RS%d", irs), "", HistAxis(AXnr, AXItrdEst));
                Hist::New(STR_FMT("hIn_TrdEst_RS%d", irs), "", HistAxis(AXnr, AXItrdEst));
            }

            Axis AXIcutflow("Cutflow", 3, 0., 3.);
            Hist::New("hIp_Cutflow", "", HistAxis(AXnr, AXIcutflow));
            Hist::New("hIn_Cutflow", "", HistAxis(AXnr, AXIcutflow));

            Hist::New("hIp_Evt", "", HistAxis(AXnr));
            Hist::New("hIn_Evt", "", HistAxis(AXnr));

            if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
                Hist::New("hTIs_TrdEst", "", HistAxis(AXtme, AXnr, AXItrdEst));
                Hist::New("hTIb_TrdEst", "", HistAxis(AXtme, AXnr, AXItrdEst));

                Hist::New("hTIp_TrdEst", "", HistAxis(AXtme, AXnr, AXItrdEst));
                Hist::New("hTIn_TrdEst", "", HistAxis(AXtme, AXnr, AXItrdEst));

                Hist::New("hTIp_Cutflow", "", HistAxis(AXtme, AXnr, AXIcutflow));
                Hist::New("hTIn_Cutflow", "", HistAxis(AXtme, AXnr, AXIcutflow));

                Hist::New("hTIp_Evt", "", HistAxis(AXtme, AXnr));
                Hist::New("hTIn_Evt", "", HistAxis(AXtme, AXnr));
            }	

            if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
                Hist::New("hIp_MCEvt", "", HistAxis(AXnr));
                Hist::New("hIn_MCEvt", "", HistAxis(AXnr));

                Hist::New("hIp_Evt_RWGT", "", HistAxis(AXnr));
                Hist::New("hIn_Evt_RWGT", "", HistAxis(AXnr));
                Hist::New("hIp_MCEvt_RWGT", "", HistAxis(AXnr));
                Hist::New("hIn_MCEvt_RWGT", "", HistAxis(AXnr));
            }


            //----  Intermedia Energy  ----//
            Axis AXMlchi("lchi", 200, -6., 8.);
            Hist::New("hMp_LchixIn", "", HistAxis(AXnr, AXMlchi));
            Hist::New("hMn_LchixIn", "", HistAxis(AXnr, AXMlchi));
            Hist::New("hMp_LchiyIn", "", HistAxis(AXnr, AXMlchi));
            Hist::New("hMn_LchiyIn", "", HistAxis(AXnr, AXMlchi));

            Hist::New("hMp_LchixL1", "", HistAxis(AXnr, AXMlchi));
            Hist::New("hMn_LchixL1", "", HistAxis(AXnr, AXMlchi));
            Hist::New("hMp_LchiyL1", "", HistAxis(AXnr, AXMlchi));
            Hist::New("hMn_LchiyL1", "", HistAxis(AXnr, AXMlchi));

            Axis AXMlasym("Lasym", 200, -2.5, 5.0);
            Hist::New("hMp_Lasym1I", "", HistAxis(AXnr, AXMlasym));
            Hist::New("hMn_Lasym1I", "", HistAxis(AXnr, AXMlasym));

            Axis AXMccest("CCest", 200, -4, 6);
            Hist::New("hMp_CCestL1", "", HistAxis(AXnr, AXMccest));
            Hist::New("hMn_CCestL1", "", HistAxis(AXnr, AXMccest));

            Axis AXMtrest("TrEst", 200, -4, 4);
            Hist::New("hMp_TrEstL1", "", HistAxis(AXnr, AXMtrest));
            Hist::New("hMn_TrEstL1", "", HistAxis(AXnr, AXMtrest));

            Axis AXMtrdEst("TRD Estimator", 100, 0.2, 1.6);
            Hist::New("hMs_TrdEst", "", HistAxis(AXnr, AXMtrdEst));
            Hist::New("hMb_TrdEst", "", HistAxis(AXnr, AXMtrdEst));

            Hist::New("hMp_TrdEst", "", HistAxis(AXnr, AXMtrdEst));
            Hist::New("hMn_TrdEst", "", HistAxis(AXnr, AXMtrdEst));

            for (Int_t icf = 12; icf <= 14; ++icf) {
                Hist::New(STR_FMT("hMs_TrdEst_CF%02d", icf), "", HistAxis(AXnr, AXMtrdEst));
                Hist::New(STR_FMT("hMb_TrdEst_CF%02d", icf), "", HistAxis(AXnr, AXMtrdEst));
                Hist::New(STR_FMT("hMp_TrdEst_CF%02d", icf), "", HistAxis(AXnr, AXMtrdEst));
                Hist::New(STR_FMT("hMn_TrdEst_CF%02d", icf), "", HistAxis(AXnr, AXMtrdEst));
            }

            for (Int_t irs = 0; irs <= 1; ++irs) {
                Hist::New(STR_FMT("hMs_TrdEst_RS%d", irs), "", HistAxis(AXnr, AXMtrdEst));
                Hist::New(STR_FMT("hMb_TrdEst_RS%d", irs), "", HistAxis(AXnr, AXMtrdEst));
                Hist::New(STR_FMT("hMp_TrdEst_RS%d", irs), "", HistAxis(AXnr, AXMtrdEst));
                Hist::New(STR_FMT("hMn_TrdEst_RS%d", irs), "", HistAxis(AXnr, AXMtrdEst));
            }

            Axis AXMcutflow("Cutflow", 3, 0., 3.);
            Hist::New("hMp_Cutflow", "", HistAxis(AXnr, AXMcutflow));
            Hist::New("hMn_Cutflow", "", HistAxis(AXnr, AXMcutflow));

            Hist::New("hMp_Evt", "", HistAxis(AXnr));
            Hist::New("hMn_Evt", "", HistAxis(AXnr));

            if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
                Hist::New("hTMs_TrdEst", "", HistAxis(AXtme, AXnr, AXMtrdEst));
                Hist::New("hTMb_TrdEst", "", HistAxis(AXtme, AXnr, AXMtrdEst));

                Hist::New("hTMp_TrdEst", "", HistAxis(AXtme, AXnr, AXMtrdEst));
                Hist::New("hTMn_TrdEst", "", HistAxis(AXtme, AXnr, AXMtrdEst));

                Hist::New("hTMp_Cutflow", "", HistAxis(AXtme, AXnr, AXMcutflow));
                Hist::New("hTMn_Cutflow", "", HistAxis(AXtme, AXnr, AXMcutflow));

                Hist::New("hTMp_Evt", "", HistAxis(AXtme, AXnr));
                Hist::New("hTMn_Evt", "", HistAxis(AXtme, AXnr));
            }	

            if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
                Hist::New("hMp_MCEvt", "", HistAxis(AXnr));
                Hist::New("hMn_MCEvt", "", HistAxis(AXnr));

                Hist::New("hMp_Evt_RWGT", "", HistAxis(AXnr));
                Hist::New("hMn_Evt_RWGT", "", HistAxis(AXnr));
                Hist::New("hMp_MCEvt_RWGT", "", HistAxis(AXnr));
                Hist::New("hMn_MCEvt_RWGT", "", HistAxis(AXnr));
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
            std::string nmSet = (iSet==SetHL1.size()) ? "" : STR_FMT("%d", iTRD);

            hHL1pCutflow.push_back(Hist::New(STR_FMT("hHL1%sp_Cutflow", nmSet.c_str()), "", AXnr, AXHL1cutflow));
            hHL1nCutflow.push_back(Hist::New(STR_FMT("hHL1%sn_Cutflow", nmSet.c_str()), "", AXnr, AXHL1cutflow));

            hHL1sTrdL.push_back(Hist::New(STR_FMT("hHL1%ss_CC", nmSet.c_str()), "", AXnr, AXHL1cc));
            hHL1bTrdL.push_back(Hist::New(STR_FMT("hHL1%sb_CC", nmSet.c_str()), "", AXnr, AXHL1cc));

            hHL1pTrdL.push_back(Hist::New(STR_FMT("hHL1%sp_CC", nmSet.c_str()), "", AXnr, AXHL1cc));
            hHL1nTrdL.push_back(Hist::New(STR_FMT("hHL1%sn_CC", nmSet.c_str()), "", AXnr, AXHL1cc));

            hHL1pEvt.push_back(Hist::New(STR_FMT("hHL1%sp_Evt", nmSet.c_str()), "", AXnr));
            hHL1nEvt.push_back(Hist::New(STR_FMT("hHL1%sn_Evt", nmSet.c_str()), "", AXnr));

            if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
            hHL1pMCEvt.push_back(Hist::New(STR_FMT("hHL1%sp_MCEvt", nmSet.c_str()), "", AXnr));
            hHL1nMCEvt.push_back(Hist::New(STR_FMT("hHL1%sn_MCEvt", nmSet.c_str()), "", AXnr));
            }

            if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
            hTHL1pCutflow.push_back(Hist::New(STR_FMT("hTHL1%sp_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXHL1cutflow));
            hTHL1nCutflow.push_back(Hist::New(STR_FMT("hTHL1%sn_Cutflow", nmSet.c_str()), "", AXtme, AXnr, AXHL1cutflow));

            hTHL1sTrdL.push_back(Hist::New(STR_FMT("hTHL1%ss_CC", nmSet.c_str()), "", AXtme, AXnr, AXHL1cc));
            hTHL1bTrdL.push_back(Hist::New(STR_FMT("hTHL1%sb_CC", nmSet.c_str()), "", AXtme, AXnr, AXHL1cc));

            hTHL1pTrdL.push_back(Hist::New(STR_FMT("hTHL1%sp_CC", nmSet.c_str()), "", AXtme, AXnr, AXHL1cc));
            hTHL1nTrdL.push_back(Hist::New(STR_FMT("hTHL1%sn_CC", nmSet.c_str()), "", AXtme, AXnr, AXHL1cc));

            hTHL1pEvt.push_back(Hist::New(STR_FMT("hTHL1%sp_Evt", nmSet.c_str()), "", AXtme, AXnr));
            hTHL1nEvt.push_back(Hist::New(STR_FMT("hTHL1%sn_Evt", nmSet.c_str()), "", AXtme, AXnr));
            }
            }	

*/

            //----  High Energy  ----//
            Axis AXHtrdEst("TRD Estimator", 100, 0.2, 1.6);
            Hist::New("hHp_TrdEstL1", "", HistAxis(AXnr, AXHtrdEst));
            Hist::New("hHn_TrdEstL1", "", HistAxis(AXnr, AXHtrdEst));
            Hist::New("hHp_TrdEstL9", "", HistAxis(AXnr, AXHtrdEst));
            Hist::New("hHn_TrdEstL9", "", HistAxis(AXnr, AXHtrdEst));
            Hist::New("hHp_TrdEstFs", "", HistAxis(AXnr, AXHtrdEst));
            Hist::New("hHn_TrdEstFs", "", HistAxis(AXnr, AXHtrdEst));

            Axis AXHlchi("lchi", 200, -6., 8.);
            Hist::New("hHp_LchixIn", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHn_LchixIn", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHp_LchixL1", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHn_LchixL1", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHp_LchixL9", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHn_LchixL9", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHp_LchixFs", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHn_LchixFs", "", HistAxis(AXnr, AXHlchi));

            Hist::New("hHp_LchiyIn", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHn_LchiyIn", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHp_LchiyL1", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHn_LchiyL1", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHp_LchiyL9", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHn_LchiyL9", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHp_LchiyFs", "", HistAxis(AXnr, AXHlchi));
            Hist::New("hHn_LchiyFs", "", HistAxis(AXnr, AXHlchi));

            Axis AXHlasym("Lasym", 200, -2.5, 5.0);
            Hist::New("hHp_Lasym1I", "", HistAxis(AXnr, AXHlasym));
            Hist::New("hHn_Lasym1I", "", HistAxis(AXnr, AXHlasym));
            Hist::New("hHp_Lasym9I", "", HistAxis(AXnr, AXHlasym));
            Hist::New("hHn_Lasym9I", "", HistAxis(AXnr, AXHlasym));
            Hist::New("hHp_Lasym91", "", HistAxis(AXnr, AXHlasym));
            Hist::New("hHn_Lasym91", "", HistAxis(AXnr, AXHlasym));

            Hist::New("hHp_LCAL1", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHn_LCAL1", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHp_LCAL9", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHn_LCAL9", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHp_LCAFs", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHn_LCAFs", "", HistAxis(AXHlchi, AXHlasym));

            Axis AXHccest("CCest", 50, -4, 6);
            Hist::New("hHp_CCestL1", "", HistAxis(AXnr, AXHccest));
            Hist::New("hHn_CCestL1", "", HistAxis(AXnr, AXHccest));
            Hist::New("hHp_CCestL9", "", HistAxis(AXnr, AXHccest));
            Hist::New("hHn_CCestL9", "", HistAxis(AXnr, AXHccest));
            Hist::New("hHp_CCestFs", "", HistAxis(AXnr, AXHccest));
            Hist::New("hHn_CCestFs", "", HistAxis(AXnr, AXHccest));

            Axis AXHtrest("TrEst", 40, -4, 4);
            Hist::New("hHp_TrEstL1", "", HistAxis(AXnr, AXHtrest));
            Hist::New("hHn_TrEstL1", "", HistAxis(AXnr, AXHtrest));
            Hist::New("hHp_TrEstL9", "", HistAxis(AXnr, AXHtrest));
            Hist::New("hHn_TrEstL9", "", HistAxis(AXnr, AXHtrest));
            Hist::New("hHp_TrEstFs", "", HistAxis(AXnr, AXHtrest));
            Hist::New("hHn_TrEstFs", "", HistAxis(AXnr, AXHtrest));

            for (Int_t icf = 12; icf <= 14; ++icf) {
                Hist::New(STR_FMT("hHp_TrEstL1_CF%02d", icf), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHn_TrEstL1_CF%02d", icf), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHp_TrEstL9_CF%02d", icf), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHn_TrEstL9_CF%02d", icf), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHp_TrEstFs_CF%02d", icf), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHn_TrEstFs_CF%02d", icf), "", HistAxis(AXnr, AXHtrest));
            }

            for (Int_t irs = 0; irs <= 1; ++irs) {
                Hist::New(STR_FMT("hHp_TrEstL1_RS%d", irs), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHn_TrEstL1_RS%d", irs), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHp_TrEstL9_RS%d", irs), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHn_TrEstL9_RS%d", irs), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHp_TrEstFs_RS%d", irs), "", HistAxis(AXnr, AXHtrest));
                Hist::New(STR_FMT("hHn_TrEstFs_RS%d", irs), "", HistAxis(AXnr, AXHtrest));
            }

            Axis AXHccest2("CCest", 200, -6, 6);
            Hist::New("hHp_LCCFs", "", HistAxis(AXHccest2, AXHccest2));
            Hist::New("hHn_LCCFs", "", HistAxis(AXHccest2, AXHccest2));

            Hist::New("hHp_LCAFsL1", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHn_LCAFsL1", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHp_LCAFsL9", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHn_LCAFsL9", "", HistAxis(AXHlchi, AXHlasym));
            Hist::New("hHp_LCAFsFs", "", HistAxis(AXHlasym, AXHlasym));
            Hist::New("hHn_LCAFsFs", "", HistAxis(AXHlasym, AXHlasym));

            Axis AXHcutflow("Cutflow", 4, 0., 4.);
            Hist::New("hHp_CutflowL1", "", HistAxis(AXnr, AXHcutflow));
            Hist::New("hHn_CutflowL1", "", HistAxis(AXnr, AXHcutflow));
            Hist::New("hHp_CutflowL9", "", HistAxis(AXnr, AXHcutflow));
            Hist::New("hHn_CutflowL9", "", HistAxis(AXnr, AXHcutflow));
            Hist::New("hHp_CutflowFs", "", HistAxis(AXnr, AXHcutflow));
            Hist::New("hHn_CutflowFs", "", HistAxis(AXnr, AXHcutflow));

            Hist::New("hHp_EvtL1", "", HistAxis(AXnr));
            Hist::New("hHn_EvtL1", "", HistAxis(AXnr));
            Hist::New("hHp_EvtL9", "", HistAxis(AXnr));
            Hist::New("hHn_EvtL9", "", HistAxis(AXnr));
            Hist::New("hHp_EvtFs", "", HistAxis(AXnr));
            Hist::New("hHn_EvtFs", "", HistAxis(AXnr));


            if (YiNtuple::CheckEventMode(YiNtuple::ISS)) {
                Hist::New("hTHp_CCestL1", "", HistAxis(AXtme, AXnr, AXHccest));
                Hist::New("hTHn_CCestL1", "", HistAxis(AXtme, AXnr, AXHccest));
                Hist::New("hTHp_CCestL9", "", HistAxis(AXtme, AXnr, AXHccest));
                Hist::New("hTHn_CCestL9", "", HistAxis(AXtme, AXnr, AXHccest));
                Hist::New("hTHp_CCestFs", "", HistAxis(AXtme, AXnr, AXHccest));
                Hist::New("hTHn_CCestFs", "", HistAxis(AXtme, AXnr, AXHccest));

                Hist::New("hTHp_TrEstL1", "", HistAxis(AXtme, AXnr, AXHtrest));
                Hist::New("hTHn_TrEstL1", "", HistAxis(AXtme, AXnr, AXHtrest));
                Hist::New("hTHp_TrEstL9", "", HistAxis(AXtme, AXnr, AXHtrest));
                Hist::New("hTHn_TrEstL9", "", HistAxis(AXtme, AXnr, AXHtrest));
                Hist::New("hTHp_TrEstFs", "", HistAxis(AXtme, AXnr, AXHtrest));
                Hist::New("hTHn_TrEstFs", "", HistAxis(AXtme, AXnr, AXHtrest));

                Hist::New("hTHp_CutflowL1", "", HistAxis(AXtme, AXnr, AXHcutflow));
                Hist::New("hTHn_CutflowL1", "", HistAxis(AXtme, AXnr, AXHcutflow));
                Hist::New("hTHp_CutflowL9", "", HistAxis(AXtme, AXnr, AXHcutflow));
                Hist::New("hTHn_CutflowL9", "", HistAxis(AXtme, AXnr, AXHcutflow));
                Hist::New("hTHp_CutflowFs", "", HistAxis(AXtme, AXnr, AXHcutflow));
                Hist::New("hTHn_CutflowFs", "", HistAxis(AXtme, AXnr, AXHcutflow));

                Hist::New("hTHp_EvtL1", "", HistAxis(AXtme, AXnr));
                Hist::New("hTHn_EvtL1", "", HistAxis(AXtme, AXnr));
                Hist::New("hTHp_EvtL9", "", HistAxis(AXtme, AXnr));
                Hist::New("hTHn_EvtL9", "", HistAxis(AXtme, AXnr));
                Hist::New("hTHp_EvtFs", "", HistAxis(AXtme, AXnr));
                Hist::New("hTHn_EvtFs", "", HistAxis(AXtme, AXnr));
            }

            if (YiNtuple::CheckEventMode(YiNtuple::MC)) {
                Hist::New("hHp_MCEvtL1", "", HistAxis(AXnr));
                Hist::New("hHn_MCEvtL1", "", HistAxis(AXnr));
                Hist::New("hHp_MCEvtL9", "", HistAxis(AXnr));
                Hist::New("hHn_MCEvtL9", "", HistAxis(AXnr));
                Hist::New("hHp_MCEvtFs", "", HistAxis(AXnr));
                Hist::New("hHn_MCEvtFs", "", HistAxis(AXnr));
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
        static TF1 * LasymPDFunc;
        static constexpr double LasymCen = 6.826895e-01;

    public :
        static constexpr double RigSclFact = (1. / 26000.);

    public :
        static TGraphAsymmErrors * PrFlux;
        static TGraphAsymmErrors * AppRatio;
        static TF1 * UnfoldRatio;
};

Bool_t DST::gSaveDST = true;
Bool_t DST::gSaveDSTTree = false;
Bool_t DST::gSaveDSTClone = false;
UInt_t DST::gUTime[2] = {0, 0};

Axis   DST::AXnr;
Axis   DST::AXir;
Axis   DST::AXtme;
Bool_t DST::WithPub = false;

TF1 * DST::LasymPDFunc = new TF1("LasymPDFunc", "TMath::Exp( -0.5 * (-x + TMath::Exp(x) ) ) / 2.50663e+00", -20, 20);

TGraphAsymmErrors * DST::PrFlux   = (TGraphAsymmErrors*)(TFile::Open("/afs/cern.ch/work/h/hchou/public/DATABASE/physics/database_antipp/database_pr.root")->Get("gr_exp2"));
//TGraphAsymmErrors * DST::AppRatio = (TGraphAsymmErrors*)(TFile::Open("/afs/cern.ch/work/h/hchou/public/DATABASE/physics/database_antipp/database_ap.root")->Get("gr_exp1"));
TGraphAsymmErrors * DST::AppRatio = new TGraphAsymmErrors((TH1D*)(TFile::Open("/afs/cern.ch/user/h/hchou/private/YiService/analysis/project/antipp5/temp/YiAna.root")->Get("hC_RatStat")));
//TF1               * DST::UnfoldRatio = new TF1("UnfoldRatio", "1.0+[0]*TMath::Power([1]*(x+[2]), [3])", 1, 1000);
TF1               * DST::UnfoldRatio = new TF1("UnfoldRatio", "1.0", 1, 1000);

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
    YiNtuple::fStopwatch.start();
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
    //if (fileList.at(0).find("Pr1800") != std::string::npos) DST::gSaveDSTTree = true;
    //if (fileList.at(0).find("MC_Ap_PL1_0_510") != std::string::npos) DST::gSaveDSTTree = true;
    //if (fileList.at(0).find("Ap1800") != std::string::npos) DST::gSaveDSTTree = true;
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
    Long64_t printRate = nentries / 50;
    if (printRate < 200000) printRate = 200000;
    if (printRate > 2000000) printRate = 2000000;

    for (Long64_t ientry = 0; ientry < nentries; ientry++) {
        if (nprocessed%printRate == 0) {
            const UInt_t MemSize = 1024;
            ProcInfo_t procinfo;
            gSystem->GetProcInfo(&procinfo);
            Long64_t memRes = procinfo.fMemResident / MemSize;
            Long64_t memVrl = procinfo.fMemVirtual  / MemSize;
            fStopwatch.stop();

            std::cout << Form("Info :: %lf %\n", 100. * float(nprocessed)/float(nentries));
            std::cout << Form("        Processed       : %ld / %ld\n", nprocessed, nentries);
            std::cout << Form("        Passed          : %ld / %ld\n", npassed, nprocessed);
            std::cout << Form("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
            std::cout << Form("        Real Time       : %9.2f (second)\n", fStopwatch.time());
            std::cout << Form("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fStopwatch.time());
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
        fStopwatch.stop();

        std::cout << Form("Info :: %lf %\n", 100. * float(nprocessed)/float(nentries));
        std::cout << Form("        Processed       : %ld / %ld\n", nprocessed, nentries);
        std::cout << Form("        Passed          : %ld / %ld\n", npassed, nprocessed);
        std::cout << Form("        Passed Ratio    : %lf %\n", ((nprocessed == 0) ? 0. : (100. * float(npassed)/float(nprocessed))));
        std::cout << Form("        Real Time       : %9.2f (second)\n", fStopwatch.time());
        std::cout << Form("        Processed Rate  : %8.2f (Hz)\n", nprocessed / fStopwatch.time());
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

    const Float_t IGRFFact = 1.2;

    // <<Pub>> 4-years < May 19, 2011 ~ May 26, 2015 >
    UInt_t period[2] = { 1305763200, 1432684800 };
    if (DST::WithPub && IsISS && (uTime < period[0] || uTime >= period[1])) return false;

    // ECAL
    Bool_t isEcalPr    = (!fMDst->hasEcal) ? false : (fMDst->ecalEst < -0.8); 
    Bool_t isEcalEl    = (!fMDst->hasEcal) ? false : (fMDst->ecalEst >  0.6); 
    Bool_t isNotEcalPr = (!fMDst->hasEcal) ? false : (fMDst->ecalEst > -0.8); 

    // Charge Estimator
    //Float_t lchiyTh[3] = { 2.20, 2.30, 1.90 }; // ~98%
    Float_t lchixL1 = (fMDst->lchix[1]);
    Float_t lchixL9 = (fMDst->lchix[2]);
    Float_t lchixFs = (fMDst->lchix[3]);
    Float_t lchixPt[3] = { lchixL1, lchixL9, lchixFs };

    Float_t lchiyL1 = (fMDst->lchiy[1]);
    Float_t lchiyL9 = (fMDst->lchiy[2]);
    Float_t lchiyFs = (fMDst->lchiy[3]);
    Float_t lchiyPt[3] = { lchiyL1, lchiyL9, lchiyFs };

    Float_t lasymTh[3] = { 0.64, 0.68, 0.62 }; // ~98%
    Float_t lasym1I = (fMDst->lasym1I / lasymTh[0]);
    Float_t lasym9I = (fMDst->lasym9I / lasymTh[1]);
    Float_t lasym91 = (fMDst->lasym91 / lasymTh[2]);

    Float_t ccestLasym1I = LasymCen + ((lasym1I >= 0) ? LasymPDFunc->Integral(0.0, (lasym1I>100?100:lasym1I)) : -LasymPDFunc->Integral((lasym1I<-100?-100:lasym1I), 0.0));
    Float_t ccestLasym9I = LasymCen + ((lasym9I >= 0) ? LasymPDFunc->Integral(0.0, (lasym9I>100?100:lasym9I)) : -LasymPDFunc->Integral((lasym9I<-100?-100:lasym9I), 0.0));
    Float_t ccestLasym91 = LasymCen + ((lasym91 >= 0) ? LasymPDFunc->Integral(0.0, (lasym91>100?100:lasym91)) : -LasymPDFunc->Integral((lasym91<-100?-100:lasym91), 0.0));
    lasym1I = TMath::ErfInverse(2.0*ccestLasym1I-1.0);
    lasym9I = TMath::ErfInverse(2.0*ccestLasym9I-1.0);
    lasym91 = TMath::ErfInverse(2.0*ccestLasym91-1.0);
    Float_t lasymPt[3] = { lasym1I, lasym9I, lasym91 };

    //Float_t ccestLchiyL1 = 0.5 * (TMath::Erf(lchiyL1) + 1.0);
    //Float_t ccestL1 = TMath::ErfInverse( (ccestLasym1I + ccestLchiyL1) - 1.0 );
    Float_t ccestL1 = (lchiyL1 + lasym1I) / std::sqrt(2.);

    //Float_t ccestLchiyL9 = 0.5 * (TMath::Erf(lchiyL9) + 1.0);
    //Float_t ccestL9 = TMath::ErfInverse( (ccestLasym9I + ccestLchiyL9) - 1.0 );
    Float_t ccestL9 = (lchiyL9 + lasym9I) / std::sqrt(2.);

    //Float_t ccestLchiyFs = 0.5 * (TMath::Erf(lchiyFs) + 1.0);
    //Float_t ccestFs = TMath::ErfInverse( (ccestLasym91 + ccestLchiyFs) - 1.0 );
    Float_t ccestFs = (lchiyFs + lasym91) / std::sqrt(2.);

    Float_t ccestPt[3] = { ccestL1, ccestL9, ccestFs };


    //Float_t ccestSf[3] = { 0.27, 0.21, 0.27 };
    //Float_t ccestSg[3] = { 0.46, 0.47, 0.50 };
    //Float_t ccestL1 = ((((lchiyL1>0&&lasym1I>0) ? std::sqrt(lchiyL1*lchiyL1+lasym1I*lasym1I) : std::max(lchiyL1, lasym1I)) - ccestSf[0]) / ccestSg[0]) / 2.0;
    //Float_t ccestL9 = ((((lchiyL9>0&&lasym9I>0) ? std::sqrt(lchiyL9*lchiyL9+lasym9I*lasym9I) : std::max(lchiyL9, lasym9I)) - ccestSf[1]) / ccestSg[1]) / 2.0;
    //Float_t ccestFs = ((((lchiyFs>0&&lasym91>0) ? std::sqrt(lchiyFs*lchiyFs+lasym91*lasym91) : std::max(lchiyFs, lasym91)) - ccestSf[2]) / ccestSg[2]) / 2.0;
    //Float_t ccestPt[3]   = { ccestL1, ccestL9, ccestFs };
    //Float_t tuneCCPnt[3] = { 25.0, 30.0, 40.0 };

    // Charge-Confusion Estimator
    //Float_t ccestTh[4] = { (std::erfc(std::log(50.0)-std::log(1+std::fabs(fMDst->trTRigIn)))+2.5), 2.5, 3.2, 2.7 }; // ~98%
    //Float_t ccestTh[4] = { 4.2, 2.5, 3.2, 2.7 }; // ~98%
    //Float_t ccestPt[4] = { 0, ccestL1, ccestL9, ccestFs };
    //Float_t tuneCCPnt[4] = { 20.0, 25.0, 30.0, 40.0 };


    /**** Low Energy Region ****/
    while (true) {
        Short_t     trPt  = ((fMDst->sign[1]!=0) ? 1 : 0);
        std::string trNm  = (trPt == 1 ? "L1" : "In");
        Short_t     sign  = fMDst->sign[trPt];
        Float_t     nrig  = fMDst->nrig[trPt];
        Float_t     lchix = fMDst->lchix[trPt];
        Float_t     lchiy = fMDst->lchiy[trPt];
        Float_t     lasym = (trPt == 1 ? lasym1I : 0);
        Float_t     ccest = (trPt == 1 ? ccestPt[0] : 0);
        if (sign == 0) break;

        // TESTCODE ////////////////////////////////
        //if (fMDst->tofCN[0] >= 2 || fMDst->tofCN[1] >= 2) break;
        if (fMDst->richCh[0]+fMDst->richCh[1] >= 15) break;
        //if (fMDst->tofCN[0]==1 && fMDst->tofCN[1] == 1) break;
        //if (fMDst->tofCN[0] > 0 || fMDst->tofCN[1] > 0) break;
        if (fMDst->tofCN[0]+fMDst->tofCN[1] > 1) break;
        //if (fMDst->richRhEl[2] > 8) break;
        //if (fMDst->richRhPr[1] > 12) break;
        if (fMDst->trNHL1 >= 2) break;
        if (fMDst->trNHL9 >= 2) break;
        //if (fMDst->tofCN[0] != 0 || fMDst->tofCN[1] != 0)
            ///std::cerr << Form("Q %8.2f %8.2f\n", fMDst->tofCQ[0][0], fMDst->tofCQ[1][0]);
        //if (fMDst->tofCQ[0][0] > 0.9) break;
        //if (fMDst->tofCQ[1][0] > 0.9) break;
        //if (fMDst->tofCN[0] != 0 || fMDst->tofCN[1] != 0)
        //    std::cerr << Form("T %8.2f %8.2f\n", fMDst->tofCT[0][0], fMDst->tofCT[1][0]);
        ///////////////////////////////////////////


        // Rigidity Scale	
        Float_t irigRS[2] = { (1./(sign*nrig)-RigSclFact),  (1./(sign*nrig)+RigSclFact) };
        Short_t signRS[2] = { ((irigRS[0]>0)?1:-1), ((irigRS[1]>0)?1:-1) };
        Float_t nrigRS[2] = { 1./std::fabs(irigRS[0]), 1./std::fabs(irigRS[1]) };

        // RTI
        Float_t cfRig = (nrig / fMDst->cfRig);
        if (IsISS && cfRig < IGRFFact) break;

        // ECAL
        if (isNotEcalPr) break;

        // TRD
        Bool_t isTrdPr = (fMDst->trdEst > 0.9);

        // Charge Confusion
        if (trPt==1 && (fMDst->lchix[0] > 3.0)) break;
        if (trPt==1 && (fMDst->lchiy[0] > 3.0)) break;

        if (trPt==0 && sign>0 && isTrdPr) Hist::Head("hLp_LchixIn")->fill(nrig, lchix, weight);
        if (trPt==0 && sign<0 && isTrdPr) Hist::Head("hLn_LchixIn")->fill(nrig, lchix, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hLp_LchixL1")->fill(nrig, lchix, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hLn_LchixL1")->fill(nrig, lchix, weight);

        if (lchix > 2.5) break;

        if (trPt==0 && sign>0 && isTrdPr) Hist::Head("hLp_LchiyIn")->fill(nrig, lchiy, weight);
        if (trPt==0 && sign<0 && isTrdPr) Hist::Head("hLn_LchiyIn")->fill(nrig, lchiy, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hLp_LchiyL1")->fill(nrig, lchiy, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hLn_LchiyL1")->fill(nrig, lchiy, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hLp_CCestL1")->fill(nrig, ccest, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hLn_CCestL1")->fill(nrig, ccest, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hLp_Lasym1I")->fill(nrig, lasym, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hLn_Lasym1I")->fill(nrig, lasym, weight);

        Bool_t CCestCut = (lasym > (0.5*lchiy+1.1) || ccest > 2.5);
        if (trPt == 1 && CCestCut) break;

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hLp_TrEstL1")->fill(nrig, lchiy, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hLn_TrEstL1")->fill(nrig, lchiy, weight);

        Float_t lchiyInCut = 2.0 * std::erf(std::log(30.0) - std::log(nrig));
        Float_t lchiyL1Cut = 2.0 * std::erf(std::log(60.0) - std::log(nrig));
        if (trPt==0 && (lchiy > lchiyInCut)) break;
        if (trPt==1 && (lchiy > lchiyL1Cut)) break;

        // RICH
        if (fMDst->richRad == -1) break;

        // RICH RING
        Bool_t isAglRing = (fMDst->richRad == 0 && fMDst->hasRich);
        Bool_t isAglPE = (isAglRing && fMDst->richMPr < -1.5 && fMDst->richMPi < 2.0);
        Bool_t isAglDE = (isAglRing && fMDst->richMPr > 2.5);

        if (isAglRing) {
            if (sign>0) Hist::Head("hLp_AglM")->fill(nrig, fMDst->richMPr, weight);
            if (sign<0) Hist::Head("hLn_AglM")->fill(nrig, fMDst->richMPi, weight);
        }

        Bool_t isNafRing = (fMDst->richRad == 1 && fMDst->hasRich);
        Bool_t isNafDE = (isNafRing && fMDst->richMPr > 2.5);
        Bool_t isNafPE = (isNafRing && fMDst->richMPr < -1.5 && fMDst->richMPi < 2.0);

        if (isNafRing) {
            if (sign>0) Hist::Head("hLp_NafM")->fill(nrig, fMDst->richMPr, weight);
            if (sign<0) Hist::Head("hLn_NafM")->fill(nrig, fMDst->richMPi, weight);
        }

        if (isNafDE) break;
        if (isNafPE) break;
        if (isAglDE) break;

        // RICH NO-RING
        Bool_t isAglNORing = (fMDst->richRad == 0 && !fMDst->hasRich && fMDst->richPhEl > 2.0);
        Bool_t isNafNORing = (fMDst->richRad == 1 && !fMDst->hasRich);

        // RICH
        if (fMDst->richRad != 0) break;

        // TRD
        if (sign>0) Hist::Head("hLp_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (sign<0) Hist::Head("hLn_TrdEst")->fill(nrig, fMDst->trdEst, weight);

        if (!isTrdPr) break;

        // Mass Estimator (Scale)
        Int_t   nrigBin = AXnr.find(nrig);
        Float_t nrigScl = (nrigBin <= AXnr.nbin()) ? AXnr(nrigBin) : AXnr.max();
        Float_t massScl = (4.0*std::exp(-0.2*nrigScl*nrigScl) + 1.0);
        Float_t tofM    = fMDst->tofM / massScl;

        // Template Fit
        Bool_t isSignal     = (sign > 0);
        Bool_t isBackground = (sign < 0) && (isAglPE);

        if (isSignal)     Hist::Head("hLs_TofM")->fill(nrig, tofM, weight);
        if (isBackground) Hist::Head("hLb_TofM")->fill(nrig, tofM, weight);
        if (IsISS && isSignal)     Hist::Head("hTLs_TofM")->fill(uTime, nrig, tofM, weight);
        if (IsISS && isBackground) Hist::Head("hTLb_TofM")->fill(uTime, nrig, tofM, weight);

        for (Int_t icf = 12; icf <= 14; ++icf) {
            if (IsISS && cfRig < (0.1*icf)) continue;
            if (isSignal)     Hist::Head(STR_FMT("hLs_TofM_CF%02d", icf))->fill(nrig, tofM, weight);
            if (isBackground) Hist::Head(STR_FMT("hLb_TofM_CF%02d", icf))->fill(nrig, tofM, weight);
        }

        for (Int_t irs = 0; irs <= 1; ++irs) {
            if (sign != signRS[irs]) continue;
            if (isSignal)     Hist::Head(STR_FMT("hLs_TofM_RS%d", irs))->fill(nrigRS[irs], tofM, weight);
            if (isBackground) Hist::Head(STR_FMT("hLb_TofM_RS%d", irs))->fill(nrigRS[irs], tofM, weight);
        }

        // RICH AGL RING
        if (isAglPE) break;

        if (sign>0) Hist::Head("hLp_TofM")->fill(nrig, tofM, weight);
        if (sign<0) Hist::Head("hLn_TofM")->fill(nrig, tofM, weight);
        if (IsISS && sign>0) Hist::Head("hTLp_TofM")->fill(uTime, nrig, tofM, weight);
        if (IsISS && sign<0) Hist::Head("hTLn_TofM")->fill(uTime, nrig, tofM, weight);

        for (Int_t icf = 12; icf <= 14; ++icf) {
            if (IsISS && cfRig < (0.1*icf)) continue;
            if (sign>0) Hist::Head(STR_FMT("hLp_TofM_CF%02d", icf))->fill(nrig, tofM, weight);
            if (sign<0) Hist::Head(STR_FMT("hLn_TofM_CF%02d", icf))->fill(nrig, tofM, weight);
        }

        for (Int_t irs = 0; irs <= 1; ++irs) {
            if (sign != signRS[irs]) continue;
            if (sign>0) Hist::Head(STR_FMT("hLp_TofM_RS%d", irs))->fill(nrigRS[irs], tofM, weight);
            if (sign<0) Hist::Head(STR_FMT("hLn_TofM_RS%d", irs))->fill(nrigRS[irs], tofM, weight);
        }

        if (sign>0) Hist::Head("hLp_Evt")->fill(nrig, weight);
        if (sign<0) Hist::Head("hLn_Evt")->fill(nrig, weight);
        if (IsISS && sign>0) Hist::Head("hTLp_Evt")->fill(uTime, nrig, weight);
        if (IsISS && sign<0) Hist::Head("hTLn_Evt")->fill(uTime, nrig, weight);

        if (IsMonteCarlo) {
            if (MCSign>0) Hist::Head("hLp_MCEvt")->fill(MCNRig, weight);
            if (MCSign<0) Hist::Head("hLn_MCEvt")->fill(MCNRig, weight);

            Float_t prwgt = std::pow(MCNRig, -1.7) * (DST::PrFlux->Eval(MCNRig) * 1e-4);
            Float_t apwgt = std::pow(MCNRig, -1.7) * (DST::PrFlux->Eval(MCNRig) * DST::AppRatio->Eval(MCNRig)) * DST::UnfoldRatio->Eval(MCNRig);
            //std::cout << DST::UnfoldRatio->Eval(MCNRig) << std::endl;
            if (MCSign>0 && sign>0) Hist::Head("hLp_Evt_RWGT")->fill(nrig, weight * prwgt);
            if (MCSign<0 && sign<0) Hist::Head("hLn_Evt_RWGT")->fill(nrig, weight * apwgt);
            if (MCSign>0) Hist::Head("hLp_MCEvt_RWGT")->fill(MCNRig, weight * prwgt);
            if (MCSign<0) Hist::Head("hLn_MCEvt_RWGT")->fill(MCNRig, weight * apwgt);
        }

        fOptL = true;
        break;	
    }	


    /**** Intermedia Energy Region ****/
    while (true) {
        Short_t     trPt  = ((fMDst->sign[1]!=0) ? 1 : 0);
        std::string trNm  = ((fMDst->sign[1]!=0) ? "L1" : "In");
        Short_t     sign  = fMDst->sign[trPt];
        Float_t     nrig  = fMDst->nrig[trPt];
        Float_t     lchix = fMDst->lchix[trPt];
        Float_t     lchiy = fMDst->lchiy[trPt];
        Float_t     lasym = (trPt == 1 ? lasym1I : 0);
        Float_t     ccest = (trPt == 1 ? ccestPt[0] : 0);
        if (sign == 0) break;

        // TESTCODE ////////////////////////////////
        //if (fMDst->tofCN[0] >= 2 || fMDst->tofCN[1] >= 2) break;
        if (fMDst->richCh[0]+fMDst->richCh[1] >= 15) break;
        //if (fMDst->tofCN[0]==1 && fMDst->tofCN[1] == 1) break;
        //if (fMDst->tofCN[0] > 0 || fMDst->tofCN[1] > 0) break;
        if (fMDst->tofCN[0]+fMDst->tofCN[1] > 1) break;
        //if (fMDst->richRhEl[2] > 8) break;
        //if (fMDst->richRhPr[1] > 12) break;
        if (fMDst->trNHL1 >= 2) break;
        if (fMDst->trNHL9 >= 2) break;
        //if (fMDst->tofCQ[0][0] > 0.9) break;
        //if (fMDst->tofCQ[1][0] > 0.9) break;
        ///////////////////////////////////////////


        // Rigidity Scale	
        Float_t irigRS[2] = { (1./(sign*nrig)-RigSclFact),  (1./(sign*nrig)+RigSclFact) };
        Short_t signRS[2] = { ((irigRS[0]>0)?1:-1), ((irigRS[1]>0)?1:-1) };
        Float_t nrigRS[2] = { 1./std::fabs(irigRS[0]), 1./std::fabs(irigRS[1]) };

        // RTI
        Float_t cfRig = (nrig / fMDst->cfRig);
        if (IsISS && cfRig < IGRFFact) break;

        // TRD
        Bool_t isTrdPr = (fMDst->trdEst > 0.9);

        // RICH
        Bool_t isAglRing  = (fMDst->richRad == 0 && fMDst->hasRich);
        if (!isAglRing) break;

        // Charge Confusion
        if (trPt==1 && (fMDst->lchix[0] > 3.0)) break;
        if (trPt==1 && (fMDst->lchiy[0] > 3.0)) break;

        if (trPt==0 && sign>0 && isTrdPr) Hist::Head("hIp_LchixIn")->fill(nrig, lchix, weight);
        if (trPt==0 && sign<0 && isTrdPr) Hist::Head("hIn_LchixIn")->fill(nrig, lchix, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hIp_LchixL1")->fill(nrig, lchix, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hIn_LchixL1")->fill(nrig, lchix, weight);

        if (lchix > 2.5) break;

        if (trPt==0 && sign>0 && isTrdPr) Hist::Head("hIp_LchiyIn")->fill(nrig, lchiy, weight);
        if (trPt==0 && sign<0 && isTrdPr) Hist::Head("hIn_LchiyIn")->fill(nrig, lchiy, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hIp_LchiyL1")->fill(nrig, lchiy, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hIn_LchiyL1")->fill(nrig, lchiy, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hIp_CCestL1")->fill(nrig, ccest, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hIn_CCestL1")->fill(nrig, ccest, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hIp_Lasym1I")->fill(nrig, lasym, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hIn_Lasym1I")->fill(nrig, lasym, weight);

        Bool_t CCestCut = (lasym > (0.5*lchiy+1.1) || ccest > 2.5);
        if (trPt == 1 && CCestCut) break;

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hIp_TrEstL1")->fill(nrig, lchiy, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hIn_TrEstL1")->fill(nrig, lchiy, weight);

        Float_t lchiyInCut = 2.0 * std::erf(std::log(30.0) - std::log(nrig));
        Float_t lchiyL1Cut = 2.0 * std::erf(std::log(60.0) - std::log(nrig));
        if (trPt==0 && (lchiy > lchiyInCut)) break;
        if (trPt==1 && (lchiy > lchiyL1Cut)) break;

        // RICH		
        Bool_t isAglDE = (isAglRing && fMDst->richMPr > 2.5);
        Bool_t isAglPE = (isAglRing && fMDst->richMPr < -1.5 && fMDst->richMPi < 2.0);

        if (isAglDE) break;

        // Template Fit
        Bool_t isSignal     = (sign > 0 && isEcalPr && !isAglPE);
        Bool_t isBackground = (sign < 0 && isEcalEl);

        if (isSignal)     Hist::Head("hIs_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (isBackground) Hist::Head("hIb_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (IsISS && isSignal)     Hist::Head("hTIs_TrdEst")->fill(uTime, nrig, fMDst->trdEst, weight);
        if (IsISS && isBackground) Hist::Head("hTIb_TrdEst")->fill(uTime, nrig, fMDst->trdEst, weight);

        for (Int_t icf = 12; icf <= 14; ++icf) {	
            if (IsISS && cfRig < (0.1*icf)) continue;
            if (isSignal)     Hist::Head(STR_FMT("hIs_TrdEst_CF%02d", icf))->fill(nrig, fMDst->trdEst, weight);
            if (isBackground) Hist::Head(STR_FMT("hIb_TrdEst_CF%02d", icf))->fill(nrig, fMDst->trdEst, weight);
        }

        for (Int_t irs = 0; irs <= 1; ++irs) {
            if (sign != signRS[irs]) continue;
            if (isSignal)     Hist::Head(STR_FMT("hIs_TrdEst_RS%d", irs))->fill(nrigRS[irs], fMDst->trdEst, weight);
            if (isBackground) Hist::Head(STR_FMT("hIb_TrdEst_RS%d", irs))->fill(nrigRS[irs], fMDst->trdEst, weight);
        }

        // ECAL
        if (fMDst->hasEcal && !isEcalPr) break;

        // RICH AGL RING
        if (isAglPE) break;

        if (sign>0) Hist::Head("hIp_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (sign<0) Hist::Head("hIn_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (IsISS && sign>0) Hist::Head("hTIp_TrdEst")->fill(uTime, nrig, fMDst->trdEst, weight);
        if (IsISS && sign<0) Hist::Head("hTIn_TrdEst")->fill(uTime, nrig, fMDst->trdEst, weight);

        for (Int_t icf = 12; icf <= 14; ++icf) {	
            if (IsISS && cfRig < (0.1*icf)) continue;
            if (sign>0) Hist::Head(STR_FMT("hIp_TrdEst_CF%02d", icf))->fill(nrig, fMDst->trdEst, weight);
            if (sign<0) Hist::Head(STR_FMT("hIn_TrdEst_CF%02d", icf))->fill(nrig, fMDst->trdEst, weight);
        }

        for (Int_t irs = 0; irs <= 1; ++irs) {
            if (sign != signRS[irs]) continue;
            if (sign>0) Hist::Head(STR_FMT("hIp_TrdEst_RS%d", irs))->fill(nrigRS[irs], fMDst->trdEst, weight);
            if (sign<0) Hist::Head(STR_FMT("hIn_TrdEst_RS%d", irs))->fill(nrigRS[irs], fMDst->trdEst, weight);
        }

        if (sign>0) Hist::Head("hIp_Evt")->fill(nrig, weight);
        if (sign<0) Hist::Head("hIn_Evt")->fill(nrig, weight);
        if (IsISS && sign>0) Hist::Head("hTIp_Evt")->fill(uTime, nrig, weight);
        if (IsISS && sign<0) Hist::Head("hTIn_Evt")->fill(uTime, nrig, weight);

        if (IsMonteCarlo) {
            if (MCSign>0) Hist::Head("hIp_MCEvt")->fill(MCNRig, weight);
            if (MCSign<0) Hist::Head("hIn_MCEvt")->fill(MCNRig, weight);

            Float_t prwgt = std::pow(MCNRig, -1.7) * (DST::PrFlux->Eval(MCNRig) * 1e-4);
            Float_t apwgt = std::pow(MCNRig, -1.7) * (DST::PrFlux->Eval(MCNRig) * DST::AppRatio->Eval(MCNRig)) * DST::UnfoldRatio->Eval(MCNRig);
            if (MCSign>0 && sign>0) Hist::Head("hIp_Evt_RWGT")->fill(nrig, weight * prwgt);
            if (MCSign<0 && sign<0) Hist::Head("hIn_Evt_RWGT")->fill(nrig, weight * apwgt);
            if (MCSign>0) Hist::Head("hIp_MCEvt_RWGT")->fill(MCNRig, weight * prwgt);
            if (MCSign<0) Hist::Head("hIn_MCEvt_RWGT")->fill(MCNRig, weight * apwgt);
        }

        fOptI = true;
        break;
    }


    //----  Intermedia Energy  ----//
    while (true) {
        Short_t     trPt  = ((fMDst->sign[1]!=0) ? 1 : 0);
        std::string trNm  = ((fMDst->sign[1]!=0) ? "L1" : "In");
        Short_t     sign  = fMDst->sign[trPt];
        Float_t     nrig  = fMDst->nrig[trPt];
        Float_t     lchix = fMDst->lchix[trPt];
        Float_t     lchiy = fMDst->lchiy[trPt];
        Float_t     lasym = (trPt == 1 ? lasym1I : 0);
        Float_t     ccest = (trPt == 1 ? ccestPt[0] : 0);
        if (sign == 0) break;

        // TESTCODE ////////////////////////////////
        //if (fMDst->tofCN[0] >= 2 || fMDst->tofCN[1] >= 2) break;
        if (fMDst->richCh[0]+fMDst->richCh[1] >= 15) break;
        //if (fMDst->tofCN[0]==1 && fMDst->tofCN[1] == 1) break;
        //if (fMDst->tofCN[0] > 0 || fMDst->tofCN[1] > 0) break;
        if (fMDst->tofCN[0]+fMDst->tofCN[1] > 1) break;
        //if (fMDst->richRhEl[2] > 8) break;
        //if (fMDst->richRhPr[1] > 12) break;
        if (fMDst->trNHL1 >= 2) break;
        if (fMDst->trNHL9 >= 2) break;
        //if (fMDst->tofCQ[0][0] > 0.9) break;
        //if (fMDst->tofCQ[1][0] > 0.9) break;
        ///////////////////////////////////////////


        // Rigidity Scale	
        Float_t irigRS[2] = { (1./(sign*nrig)-RigSclFact),  (1./(sign*nrig)+RigSclFact) };
        Short_t signRS[2] = { ((irigRS[0]>0)?1:-1), ((irigRS[1]>0)?1:-1) };
        Float_t nrigRS[2] = { 1./std::fabs(irigRS[0]), 1./std::fabs(irigRS[1]) };

        // RTI
        Float_t cfRig = (nrig / fMDst->cfRig);
        if (IsISS && cfRig < IGRFFact) break;

        // TRD
        Bool_t isTrdPr = (fMDst->trdEst > 0.9);

        // Charge Confusion
        if (trPt==1 && (fMDst->lchix[0] > 3.0)) break;
        if (trPt==1 && (fMDst->lchiy[0] > 3.0)) break;

        if (trPt==0 && sign>0 && isTrdPr) Hist::Head("hMp_LchixIn")->fill(nrig, lchix, weight);
        if (trPt==0 && sign<0 && isTrdPr) Hist::Head("hMn_LchixIn")->fill(nrig, lchix, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hMp_LchixL1")->fill(nrig, lchix, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hMn_LchixL1")->fill(nrig, lchix, weight);

        if (lchix > 2.5) break;

        if (trPt==0 && sign>0 && isTrdPr) Hist::Head("hMp_LchiyIn")->fill(nrig, lchiy, weight);
        if (trPt==0 && sign<0 && isTrdPr) Hist::Head("hMn_LchiyIn")->fill(nrig, lchiy, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hMp_LchiyL1")->fill(nrig, lchiy, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hMn_LchiyL1")->fill(nrig, lchiy, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hMp_Lasym1I")->fill(nrig, lasym, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hMn_Lasym1I")->fill(nrig, lasym, weight);

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hMp_CCestL1")->fill(nrig, ccest, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hMn_CCestL1")->fill(nrig, ccest, weight);

        Bool_t CCestCut = (lasym > (0.5*lchiy+1.1) || ccest > 2.5);
        if (trPt == 1 && CCestCut) break;

        if (trPt==1 && sign>0 && isTrdPr) Hist::Head("hMp_TrEstL1")->fill(nrig, lchiy, weight);
        if (trPt==1 && sign<0 && isTrdPr) Hist::Head("hMn_TrEstL1")->fill(nrig, lchiy, weight);

        Float_t lchiyInCut = 2.0 * std::erf(std::log(30.0) - std::log(nrig));
        Float_t lchiyL1Cut = 2.0 * std::erf(std::log(60.0) - std::log(nrig));
        if (trPt==0 && (lchiy > lchiyInCut)) break;
        if (trPt==1 && (lchiy > lchiyL1Cut)) break;

        // RICH
        Bool_t isAglRing  = (fMDst->richRad == 0 && fMDst->hasRich);
        Bool_t isAglDE = (isAglRing && fMDst->richMPr > 2.5);
        Bool_t isAglPE = (isAglRing && fMDst->richMPr < -1.5 && fMDst->richMPi < 2.0);

        if (isAglDE) break;

        // Template Fit
        Bool_t isSignal     = (sign > 0 && isEcalPr && !isAglPE);
        Bool_t isBackground = (sign < 0 && isEcalEl);

        if (isSignal)     Hist::Head("hMs_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (isBackground) Hist::Head("hMb_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (IsISS && isSignal)     Hist::Head("hTMs_TrdEst")->fill(uTime, nrig, fMDst->trdEst, weight);
        if (IsISS && isBackground) Hist::Head("hTMb_TrdEst")->fill(uTime, nrig, fMDst->trdEst, weight);

        for (Int_t icf = 12; icf <= 14; ++icf) {	
            if (IsISS && cfRig < (0.1*icf)) continue;
            if (isSignal)     Hist::Head(STR_FMT("hMs_TrdEst_CF%02d", icf))->fill(nrig, fMDst->trdEst, weight);
            if (isBackground) Hist::Head(STR_FMT("hMb_TrdEst_CF%02d", icf))->fill(nrig, fMDst->trdEst, weight);
        }

        for (Int_t irs = 0; irs <= 1; ++irs) {
            if (sign != signRS[irs]) continue;
            if (isSignal)     Hist::Head(STR_FMT("hMs_TrdEst_RS%d", irs))->fill(nrigRS[irs], fMDst->trdEst, weight);
            if (isBackground) Hist::Head(STR_FMT("hMb_TrdEst_RS%d", irs))->fill(nrigRS[irs], fMDst->trdEst, weight);
        }

        // ECAL
        if (fMDst->hasEcal && !isEcalPr) break;

        // RICH AGL RING
        if (isAglPE) break;

        if (sign>0) Hist::Head("hMp_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (sign<0) Hist::Head("hMn_TrdEst")->fill(nrig, fMDst->trdEst, weight);
        if (IsISS && sign>0) Hist::Head("hTMp_TrdEst")->fill(uTime, nrig, fMDst->trdEst, weight);
        if (IsISS && sign<0) Hist::Head("hTMn_TrdEst")->fill(uTime, nrig, fMDst->trdEst, weight);

        for (Int_t icf = 12; icf <= 14; ++icf) {	
            if (IsISS && cfRig < (0.1*icf)) continue;
            if (sign>0) Hist::Head(STR_FMT("hMp_TrdEst_CF%02d", icf))->fill(nrig, fMDst->trdEst, weight);
            if (sign<0) Hist::Head(STR_FMT("hMn_TrdEst_CF%02d", icf))->fill(nrig, fMDst->trdEst, weight);
        }

        for (Int_t irs = 0; irs <= 1; ++irs) {
            if (sign != signRS[irs]) continue;
            if (sign>0) Hist::Head(STR_FMT("hMp_TrdEst_RS%d", irs))->fill(nrigRS[irs], fMDst->trdEst, weight);
            if (sign<0) Hist::Head(STR_FMT("hMn_TrdEst_RS%d", irs))->fill(nrigRS[irs], fMDst->trdEst, weight);
        }

        if (sign>0) Hist::Head("hMp_Evt")->fill(nrig, weight);
        if (sign<0) Hist::Head("hMn_Evt")->fill(nrig, weight);
        if (IsISS && sign>0) Hist::Head("hTMp_Evt")->fill(uTime, nrig, weight);
        if (IsISS && sign<0) Hist::Head("hTMn_Evt")->fill(uTime, nrig, weight);

        if (IsMonteCarlo) {
            if (MCSign>0) Hist::Head("hMp_MCEvt")->fill(MCNRig, weight);
            if (MCSign<0) Hist::Head("hMn_MCEvt")->fill(MCNRig, weight);

            Float_t prwgt = std::pow(MCNRig, -1.7) * (DST::PrFlux->Eval(MCNRig) * 1e-4);
            Float_t apwgt = std::pow(MCNRig, -1.7) * (DST::PrFlux->Eval(MCNRig) * DST::AppRatio->Eval(MCNRig)) * DST::UnfoldRatio->Eval(MCNRig);
            if (MCSign>0 && sign>0) Hist::Head("hMp_Evt_RWGT")->fill(nrig, weight * prwgt);
            if (MCSign<0 && sign<0) Hist::Head("hMn_Evt_RWGT")->fill(nrig, weight * apwgt);
            if (MCSign>0) Hist::Head("hMp_MCEvt_RWGT")->fill(MCNRig, weight * prwgt);
            if (MCSign<0) Hist::Head("hMn_MCEvt_RWGT")->fill(MCNRig, weight * apwgt);
        }

        fOptM = true;
        break;
    }


    /**** High Energy Region ****/
    while (true) {
        if (fMDst->trPt == 0) break; 
        Short_t   trPt   = fMDst->trPt;
        std::string trNm = "In";
        std::string trNs = "UL";
        if      (trPt == 1) { trNm = "L1"; trNs = "1I"; }
        else if (trPt == 2) { trNm = "L9"; trNs = "9I"; }
        else if (trPt == 3) { trNm = "Fs"; trNs = "91"; }
        Short_t     sign  = fMDst->sign[trPt];
        Float_t     nrig  = fMDst->nrig[trPt];
        Float_t     lchix = fMDst->lchix[trPt];
        Float_t     lchiy = lchiyPt[trPt-1];
        Float_t     lasym = lasymPt[trPt-1];
        Float_t     ccest = ccestPt[trPt-1];
        if (sign == 0) break;

        // TESTCODE ////////////////////////////////
        //if (fMDst->tofCN[0] >= 2 || fMDst->tofCN[1] >= 2) break;
        if (fMDst->richCh[0]+fMDst->richCh[1] >= 15) break;
        //if (fMDst->tofCN[0]==1 && fMDst->tofCN[1] == 1) break;
        //if (fMDst->tofCN[0] > 0 || fMDst->tofCN[1] > 0) break;
        if (fMDst->tofCN[0]+fMDst->tofCN[1] > 1) break;
        //if (fMDst->richRhEl[2] > 8) break;
        //if (fMDst->richRhPr[1] > 12) break;
        if (fMDst->trNHL1 >= 2) break;
        if (fMDst->trNHL9 >= 2) break;
        //if (fMDst->tofCQ[0][0] > 0.9) break;
        //if (fMDst->tofCQ[1][0] > 0.9) break;
        ///////////////////////////////////////////


        // Rigidity Scale	
        Float_t irigRS[2] = { (1./(sign*nrig)-RigSclFact),  (1./(sign*nrig)+RigSclFact) };
        Short_t signRS[2] = { ((irigRS[0]>0)?1:-1), ((irigRS[1]>0)?1:-1) };
        Float_t nrigRS[2] = { 1./std::fabs(irigRS[0]), 1./std::fabs(irigRS[1]) };

        // RTI
        Float_t cfRig = (nrig / fMDst->cfRig);
        if (IsISS && cfRig < IGRFFact) break;

        // ECAL
        if (fMDst->hasEcal && !isEcalPr) break;

        // TRD
        Float_t trdCut = (0.8 + 0.05 * std::erf(std::log(50.) - std::log(nrig)) );
        Bool_t isTrdPr = (fMDst->trdEst > trdCut);

        // TRK
        if (isTrdPr && sign>0) Hist::Head("hHp_LchixIn")->fill(nrig, fMDst->lchix[0], weight);
        if (isTrdPr && sign<0) Hist::Head("hHn_LchixIn")->fill(nrig, fMDst->lchix[0], weight);

        if (fMDst->lchix[0] > 3.0) break;

        if (isTrdPr && sign>0) Hist::Head("hHp_LchiyIn")->fill(nrig, fMDst->lchiy[0], weight);
        if (isTrdPr && sign<0) Hist::Head("hHn_LchiyIn")->fill(nrig, fMDst->lchiy[0], weight);

        if (fMDst->lchiy[0] > 3.0) break;

        if (isTrdPr && sign>0) Hist::Head(STR_FMT("hHp_Lchix%s", trNm.c_str()))->fill(nrig, lchix, weight);
        if (isTrdPr && sign<0) Hist::Head(STR_FMT("hHn_Lchix%s", trNm.c_str()))->fill(nrig, lchix, weight);

        if (lchix > 2.5) break;

        // TRD	
        if (sign>0) Hist::Head(STR_FMT("hHp_TrdEst%s", trNm.c_str()))->fill(nrig, fMDst->trdEst, weight);
        if (sign<0) Hist::Head(STR_FMT("hHn_TrdEst%s", trNm.c_str()))->fill(nrig, fMDst->trdEst, weight);
        if (!isTrdPr) break;

        if (trPt == 3) {	
            const Float_t ccestFsCut = 2.5;
            Float_t ccestFsX = (ccestPt[0] - ccestPt[1]) / sqrt(2.) / 0.6;
            Float_t ccestFsY = (ccestPt[0] + ccestPt[1]) / sqrt(2.) / 1.3;
            if (sign>0 && nrig > 80 && nrig < 450) Hist::Head("hHp_LCCFs")->fill(ccestFsX, ccestFsY, weight);
            if (sign<0 && nrig > 80 && nrig < 450) Hist::Head("hHn_LCCFs")->fill(ccestFsX, ccestFsY, weight);

            //if (std::sqrt(ccestFsX*ccestFsX+ccestFsY*ccestFsY) > ccestFsCut) break;
            if (((ccestFsY > 0) ? 
                        std::sqrt(ccestFsX*ccestFsX+ccestFsY*ccestFsY) > ccestFsCut : 
                        std::fabs(ccestFsX) > ccestFsCut)) break;

            if (sign>0 && nrig > 80 && nrig < 450) Hist::Head(STR_FMT("hHp_LCAFsL1", trNm.c_str()))->fill(lchiyPt[0], lasymPt[0], weight);
            if (sign<0 && nrig > 80 && nrig < 450) Hist::Head(STR_FMT("hHn_LCAFsL1", trNm.c_str()))->fill(lchiyPt[0], lasymPt[0], weight);
            if (sign>0 && nrig > 80 && nrig < 450) Hist::Head(STR_FMT("hHp_LCAFsL9", trNm.c_str()))->fill(lchiyPt[1], lasymPt[1], weight);
            if (sign<0 && nrig > 80 && nrig < 450) Hist::Head(STR_FMT("hHn_LCAFsL9", trNm.c_str()))->fill(lchiyPt[1], lasymPt[1], weight);

            Bool_t LasymCutL1 = (lasymPt[0] > (0.5*lchiyPt[0]+1.1));
            Bool_t LasymCutL9 = (lasymPt[1] > (0.5*lchiyPt[1]+1.1));
            if (LasymCutL1) break;
            if (LasymCutL9) break;

            const Float_t lasymFsCut = 2.0;
            Float_t lasymFsX = (lasymPt[0] - lasymPt[1]) / sqrt(2.) / 0.6;
            Float_t lasymFsY = (lasymPt[0] + lasymPt[1]) / sqrt(2.) / 0.8; // only Y has meaning

            if (sign>0 && nrig > 80 && nrig < 450) Hist::Head(STR_FMT("hHp_LCAFsFs", trNm.c_str()))->fill(lasymFsX, lasymFsY, weight);
            if (sign<0 && nrig > 80 && nrig < 450) Hist::Head(STR_FMT("hHn_LCAFsFs", trNm.c_str()))->fill(lasymFsX, lasymFsY, weight);

            //if (((lasymFsY > 0) ? 
            //	 std::sqrt(lasymFsX*lasymFsX+lasymFsY*lasymFsY) > lasymFsCut : 
            //	 std::fabs(lasymFsX) > lasymFsCut)) break;
        }

        if (sign>0) Hist::Head(STR_FMT("hHp_Lchiy%s", trNm.c_str()))->fill(nrig, lchiy, weight);
        if (sign<0) Hist::Head(STR_FMT("hHn_Lchiy%s", trNm.c_str()))->fill(nrig, lchiy, weight);

        if (sign>0) Hist::Head(STR_FMT("hHp_Lasym%s", trNs.c_str()))->fill(nrig, lasym, weight);
        if (sign<0) Hist::Head(STR_FMT("hHn_Lasym%s", trNs.c_str()))->fill(nrig, lasym, weight);

        if (sign>0 && nrig > 80 && nrig < 450) Hist::Head(STR_FMT("hHp_LCA%s", trNm.c_str()))->fill(lchiy, lasym, weight);
        if (sign<0 && nrig > 80 && nrig < 450) Hist::Head(STR_FMT("hHn_LCA%s", trNm.c_str()))->fill(lchiy, lasym, weight);

        // Charge Confusion
        Bool_t LasymCut = (lasym > (0.5*lchiy+1.1));
        if (LasymCut) break;

        if (sign>0) Hist::Head(STR_FMT("hHp_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
        if (sign<0) Hist::Head(STR_FMT("hHn_CCest%s", trNm.c_str()))->fill(nrig, ccest, weight);
        if (IsISS && sign>0) Hist::Head(STR_FMT("hTHp_CCest%s", trNm.c_str()))->fill(uTime, nrig, ccest, weight);
        if (IsISS && sign<0) Hist::Head(STR_FMT("hTHn_CCest%s", trNm.c_str()))->fill(uTime, nrig, ccest, weight);

        // Best One, Now
        Bool_t CCestCut = (ccest > 2.0);
        //Bool_t CCestCut = (ccest > 2.5);
        if (CCestCut) break;

        if (sign>0) Hist::Head(STR_FMT("hHp_TrEst%s", trNm.c_str()))->fill(nrig, lchiy, weight);
        if (sign<0) Hist::Head(STR_FMT("hHn_TrEst%s", trNm.c_str()))->fill(nrig, lchiy, weight);

        for (Int_t icf = 12; icf <= 14; ++icf) {
            if (IsISS && cfRig < (0.1*icf)) continue;
            if (sign>0) Hist::Head(STR_FMT("hHp_TrEst%s_CF%02d", trNm.c_str(), icf))->fill(nrig, lchiy, weight);
            if (sign<0) Hist::Head(STR_FMT("hHn_TrEst%s_CF%02d", trNm.c_str(), icf))->fill(nrig, lchiy, weight);
        }

        for (Int_t irs = 0; irs <= 1; ++irs) {
            if (sign != signRS[irs]) continue;
            if (sign>0) Hist::Head(STR_FMT("hHp_TrEst%s_RS%d", trNm.c_str(), irs))->fill(nrigRS[irs], lchiy, weight);
            if (sign<0) Hist::Head(STR_FMT("hHn_TrEst%s_RS%d", trNm.c_str(), irs))->fill(nrigRS[irs], lchiy, weight);
        }

        if (IsISS && sign>0) Hist::Head(STR_FMT("hTHp_TrEst%s", trNm.c_str()))->fill(uTime, nrig, lchiy, weight);
        if (IsISS && sign<0) Hist::Head(STR_FMT("hTHn_TrEst%s", trNm.c_str()))->fill(uTime, nrig, lchiy, weight);

        if (sign>0) Hist::Head(STR_FMT("hHp_Evt%s", trNm.c_str()))->fill(nrig, weight);
        if (sign<0) Hist::Head(STR_FMT("hHn_Evt%s", trNm.c_str()))->fill(nrig, weight);

        if (IsISS && sign>0) Hist::Head(STR_FMT("hTHp_Evt%s", trNm.c_str()))->fill(uTime, nrig, weight);
        if (IsISS && sign<0) Hist::Head(STR_FMT("hTHn_Evt%s", trNm.c_str()))->fill(uTime, nrig, weight);

        if (IsMonteCarlo) {
            if (MCSign>0) Hist::Head(STR_FMT("hHp_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
            if (MCSign<0) Hist::Head(STR_FMT("hHn_MCEvt%s", trNm.c_str()))->fill(MCNRig, weight);
        }

        if (trPt == 1) fOptHL1 = true;
        if (trPt == 2) fOptHL9 = true;
        if (trPt == 3) fOptHFs = true;
        break;
    }

    Bool_t isSave = (fOptL || fOptI || fOptM || fOptHL1 || fOptHL9 || fOptHFs);
    if (!isSave) return false;

    fMCRig = MCSign * MCNRig;
    fTRRig = fMDst->sign[0] * fMDst->nrig[0];

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
    Style::LoadDefaultEnvironment();

    std::cout << std::endl << std::endl;
    std::cout << "Usage : YiAnaNtuple event_mode file_list group_th group_size (path)\n";
    std::cout << "    Parameters : \n";
    std::cout << "    event_mode [ISS BT MC]\n";
    std::cout << "    file_list\n";
    std::cout << "    group_th\n";
    std::cout << "    group_size\n";
    std::cout << "    (path)\n";
    std::cout << std::endl << std::endl;

    if (argc != 5 && argc != 6) MGSys::ShowErrorAndExit(LOC_ADDR(), MGSys::Message("Number of argument is not conform! Exiting ..."));

    std::string event_mode = argv[1];
    std::string file_list = argv[2];
    Long64_t group_th = atol(argv[3]);
    Long64_t group_size = atol(argv[4]);

    std::use_facet<std::ctype<char> >(std::locale()).toupper(&event_mode[0], &event_mode[0] + event_mode.size());
    if (event_mode == "ISS") YiNtuple::SetEventMode(YiNtuple::ISS);
    else if (event_mode == "BT") YiNtuple::SetEventMode(YiNtuple::BT);
    else if (event_mode == "MC") YiNtuple::SetEventMode(YiNtuple::MC);
    else MGSys::ShowErrorAndExit(LOC_ADDR(), MGSys::Message("Can't find event mode (ISS, BT, MC)! Exiting ..."));

    std::string outputFile = STR_FMT("YiAnalytics_%s.%07ld.root", event_mode.c_str(), group_th);
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
