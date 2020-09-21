#ifndef __Variables_H__
#define __Variables_H__

#include "TTree.h"
#include "TChain.h"

#ifdef __USE_TMVA_READER__
#include <TMVA/Reader.h>
#endif // __USE_TMVA_READER__

class Vars {
    public :
        Vars() {
            init_base();
            tree = nullptr;
            chain = nullptr;
#ifdef __USE_TMVA_READER__
            tmva_reader = nullptr;
#endif // __USE_TMVA_READER__
        }
        ~Vars() {}

        TTree* get_tree() { return ((tree != nullptr) ? tree : nullptr); }
        TChain* get_chain() { return ((chain != nullptr) ? chain : nullptr); }
        
#ifdef __USE_TMVA_READER__
        TMVA::Reader* get_tmva_reader() { return ((tmva_reader != nullptr) ? tmva_reader : nullptr); }
        Double_t get_tmva_value() { return ((tmva_reader != nullptr) ? tmva_reader->EvaluateMVA(tmva_name.c_str()) : 0); }
#endif // __USE_TMVA_READER__

        void init_base() {
            run = 0;
            evt = 0;
            ut  = 0;
            wgt = 1;

            mc = 0;
            mc_rig = 0;
            mc_w10 = 1;
            mc_w27 = 1;

            cfsec = 0;
            cfevt = 0;
            lv    = 1;
            
            trg   = 0;
        }

        virtual void branch() = 0;
        void branch_base() {
            if (!tree) return;
            tree->Branch("run"   , &run);
            tree->Branch("evt"   , &evt);
            tree->Branch("ut"    , &ut);
            tree->Branch("wgt"   , &wgt);
            tree->Branch("mc"    , &mc);
            tree->Branch("mc_rig", &mc_rig);
            tree->Branch("mc_w10", &mc_w10);
            tree->Branch("mc_w27", &mc_w27);
            tree->Branch("cfsec" , &cfsec);
            tree->Branch("cfevt" , &cfevt);
            tree->Branch("lv"    , &lv);
            tree->Branch("trg"   , &trg);
        }
        
        void set_tree(const char* name, const char* title = "") {
            if (tree) return;
            tree = new TTree(name, title);
            branch_base();
            branch();
        }
        
        virtual void set_branch_address() = 0;
        void set_branch_address_base() {
            if (!chain) return;
            chain->SetBranchAddress("run"   , &run);
            chain->SetBranchAddress("evt"   , &evt);
            chain->SetBranchAddress("ut"    , &ut);
            chain->SetBranchAddress("wgt"   , &wgt);
            chain->SetBranchAddress("mc"    , &mc);
            chain->SetBranchAddress("mc_rig", &mc_rig);
            chain->SetBranchAddress("mc_w10", &mc_w10);
            chain->SetBranchAddress("mc_w27", &mc_w27);
            chain->SetBranchAddress("cfsec" , &cfsec);
            chain->SetBranchAddress("cfevt" , &cfevt);
            chain->SetBranchAddress("lv"    , &lv);
            chain->SetBranchAddress("trg"   , &trg);
        }
        
        void set_tree(TChain* ref) {
            if (chain || !ref) return;
            chain = ref;
            set_branch_address_base();
            set_branch_address();
        }
        
#ifdef __USE_TMVA_READER__
        virtual void set_tmva_vars() = 0;
        void set_tmva_vars_base() {
            if (!tmva_reader) return;
        }

        void set_tmva_reader(const char* name, const char* path) {
            if (tmva_reader) return;
            tmva_name = name;
            tmva_path = path;
            tmva_reader = new TMVA::Reader();
            set_tmva_vars_base();
            set_tmva_vars();
            tmva_reader->BookMVA(name, path);
        }
#endif // __USE_TMVA_READER__
        
        void fill() {
            if (!tree) return;
            tree->Fill();
        }

    protected :
        TTree*  tree;
        TChain* chain;

#ifdef __USE_TMVA_READER__
        std::string tmva_name;
        std::string tmva_path;
        TMVA::Reader* tmva_reader;
#endif // __USE_TMVA_READER__
        
    public :
        UInt_t  run;
        UInt_t  evt;
        UInt_t  ut;
        Float_t wgt;
       
        Short_t mc;
        Float_t mc_rig;
        Float_t mc_w10;
        Float_t mc_w27;

        Float_t cfsec;
        Float_t cfevt;
        Float_t lv;
        
        Short_t trg;
};


class VarsBT : public Vars {
    public :
        VarsBT() : Vars() { init(); }
        ~VarsBT() {}
        
        void init() {
            init_base();

            L2  = 0;
            nhx = 0;
            nhy = 0;
           
            bt_rig = 0;
            
            tdQ = 0;
            tdN = 0;
            std::fill_n(tdLQ, 20, -1);
            std::fill_n(tdLL, 20, 0);

            tfN = 0;
            tfQ = 0;
            std::fill_n(tfL , 4, 0);
            std::fill_n(tfLQ, 4, -1);
            
            tkQIn = 0;
            std::fill_n(tkL  , 7, 0);
            std::fill_n(tkLQ , 7, -1);
            std::fill_n(tkLQx, 7, -1);
            std::fill_n(tkLQy, 7, -1);
            
            std::fill_n(tkHitL, 9, 0);
            std::fill_n(tkHitX, 9, 0);
            std::fill_n(tkHitY, 9, 0);
            std::fill_n(tkHitZ, 9, 0);
           
            std::fill_n(ck        , 4, 0);
            std::fill_n(ck_rig    , 4, 0);
            std::fill_n(ck_nchi[0], 4*2, 0);
            std::fill_n(ck_cpu    , 4, 0);
            
            std::fill_n(kf        , 4, 0);
            std::fill_n(kf_rig    , 4, 0);
            std::fill_n(kf_nchi[0], 4*2, 0);
            std::fill_n(kf_cpu    , 4, 0);
            
            std::fill_n(hc        , 4, 0);
            std::fill_n(hc_rig    , 4, 0);
            std::fill_n(hc_nchi[0], 4*2, 0);
            std::fill_n(hc_cpu    , 4, 0);
            
            std::fill_n(new_hc        , 4, 0);
            std::fill_n(new_hc_rig    , 4, 0);
            std::fill_n(new_hc_nchi[0], 4*2, 0);

            new_tdQ      = 0;
            new_tdQ_nchi = 0;
            
            new_tfQ      = 0;
            new_tfQ_nchi = 0;
            
            new_tkQ      = 0;
            new_tkQ_nchi = 0;
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("L2" , &L2);
            tree->Branch("nhx", &nhx);
            tree->Branch("nhy", &nhy);
            
            tree->Branch("bt_rig", &bt_rig);
            
            tree->Branch("tdQ" , &tdQ);
            tree->Branch("tdN" , &tdN);
            tree->Branch("tdLQ",  tdLQ, "tdLQ[20]/F");
            tree->Branch("tdLL",  tdLL, "tdLL[20]/F");
            
            tree->Branch("tfN" , &tfN);
            tree->Branch("tfQ" , &tfQ);
            tree->Branch("tfL" ,  tfL,  "tfL[4]/S");
            tree->Branch("tfLQ",  tfLQ, "tfLQ[4]/F");
            
            tree->Branch("tkQIn", &tkQIn);
            tree->Branch("tkL"  ,  tkL  , "tkL[7]/S");
            tree->Branch("tkLQ" ,  tkLQ , "tkLQ[7]/F");
            tree->Branch("tkLQx",  tkLQx, "tkLQx[7]/F");
            tree->Branch("tkLQy",  tkLQy, "tkLQy[7]/F");
            
            tree->Branch("tkHitL", tkHitL, "tkHitL[9]/S");
            tree->Branch("tkHitX", tkHitX, "tkHitX[9]/D");
            tree->Branch("tkHitY", tkHitY, "tkHitY[9]/D");
            tree->Branch("tkHitZ", tkHitZ, "tkHitZ[9]/D");
            
            tree->Branch("ck"     , ck     , "ck[4]/S");
            tree->Branch("ck_rig" , ck_rig , "ck_rig[4]/F");
            tree->Branch("ck_nchi", ck_nchi, "ck_nchi[4][2]/F");
            tree->Branch("ck_cpu" , ck_cpu , "ck_cpu[4]/F");
            
            tree->Branch("kf"     , kf     , "kf[4]/S");
            tree->Branch("kf_rig" , kf_rig , "kf_rig[4]/F");
            tree->Branch("kf_nchi", kf_nchi, "kf_nchi[4][2]/F");
            tree->Branch("kf_cpu" , kf_cpu , "kf_cpu[4]/F");
            
            tree->Branch("hc"     , hc     , "hc[4]/S");
            tree->Branch("hc_rig" , hc_rig , "hc_rig[4]/F");
            tree->Branch("hc_nchi", hc_nchi, "hc_nchi[4][2]/F");
            tree->Branch("hc_cpu" , hc_cpu , "hc_cpu[4]/F");
            
            tree->Branch("new_hc"     , new_hc     , "new_hc[4]/S");
            tree->Branch("new_hc_rig" , new_hc_rig , "new_hc_rig[4]/F");
            tree->Branch("new_hc_nchi", new_hc_nchi, "new_hc_nchi[4][2]/F");
            
            tree->Branch("new_tdQ"     , &new_tdQ);
            tree->Branch("new_tdQ_nchi", &new_tdQ_nchi);
            
            tree->Branch("new_tfQ"     , &new_tfQ);
            tree->Branch("new_tfQ_nchi", &new_tfQ_nchi);
            
            tree->Branch("new_tkQ"     , &new_tkQ);
            tree->Branch("new_tkQ_nchi", &new_tkQ_nchi);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("L2" , &L2);
            chain->SetBranchAddress("nhx", &nhx);
            chain->SetBranchAddress("nhy", &nhy);
            
            chain->SetBranchAddress("bt_rig", &bt_rig);
            
            chain->SetBranchAddress("tdQ" , &tdQ);
            chain->SetBranchAddress("tdN" , &tdN);
            chain->SetBranchAddress("tdLQ",  tdLQ);
            chain->SetBranchAddress("tdLL",  tdLL);
            
            chain->SetBranchAddress("tfN" , &tfN);
            chain->SetBranchAddress("tfQ" , &tfQ);
            chain->SetBranchAddress("tfL" ,  tfL);
            chain->SetBranchAddress("tfLQ",  tfLQ);
            
            chain->SetBranchAddress("tkQIn", &tkQIn);
            chain->SetBranchAddress("tkL"  ,  tkL);
            chain->SetBranchAddress("tkLQ" ,  tkLQ);
            chain->SetBranchAddress("tkLQx",  tkLQx);
            chain->SetBranchAddress("tkLQy",  tkLQy);
            
            chain->SetBranchAddress("tkHitL", tkHitL);
            chain->SetBranchAddress("tkHitX", tkHitX);
            chain->SetBranchAddress("tkHitY", tkHitY);
            chain->SetBranchAddress("tkHitZ", tkHitZ);
            
            chain->SetBranchAddress("ck"     , ck     );
            chain->SetBranchAddress("ck_rig" , ck_rig );
            chain->SetBranchAddress("ck_nchi", ck_nchi);
            chain->SetBranchAddress("ck_cpu" , ck_cpu );
            
            chain->SetBranchAddress("kf"     , kf     );
            chain->SetBranchAddress("kf_rig" , kf_rig );
            chain->SetBranchAddress("kf_nchi", kf_nchi);
            chain->SetBranchAddress("kf_cpu" , kf_cpu );
            
            chain->SetBranchAddress("hc"     , hc     );
            chain->SetBranchAddress("hc_rig" , hc_rig );
            chain->SetBranchAddress("hc_nchi", hc_nchi);
            chain->SetBranchAddress("hc_cpu" , hc_cpu );
            
            chain->SetBranchAddress("new_hc"     , new_hc     );
            chain->SetBranchAddress("new_hc_rig" , new_hc_rig );
            chain->SetBranchAddress("new_hc_nchi", new_hc_nchi);
            
            chain->SetBranchAddress("new_tdQ"     , &new_tdQ);
            chain->SetBranchAddress("new_tdQ_nchi", &new_tdQ_nchi);
            
            chain->SetBranchAddress("new_tfQ"     , &new_tfQ);
            chain->SetBranchAddress("new_tfQ_nchi", &new_tfQ_nchi);
            
            chain->SetBranchAddress("new_tkQ"     , &new_tkQ);
            chain->SetBranchAddress("new_tkQ_nchi", &new_tkQ_nchi);
        }

#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t L2;
        Short_t nhx;
        Short_t nhy;
        
        Float_t bt_rig;

        Float_t tdQ;
        Short_t tdN;
        Float_t tdLQ[20];
        Float_t tdLL[20];
        
        Short_t tfN;
        Float_t tfQ;
        Short_t tfL[4];
        Float_t tfLQ[4];

        Float_t tkQIn;
        Short_t tkL[7];
        Float_t tkLQ[7];
        Float_t tkLQx[7];
        Float_t tkLQy[7];

        Short_t  tkHitL[9];
        Double_t tkHitX[9];
        Double_t tkHitY[9];
        Double_t tkHitZ[9];

        Short_t ck[4];
        Float_t ck_rig[4];
        Float_t ck_nchi[4][2];
        Float_t ck_cpu[4];
        
        Short_t kf[4];
        Float_t kf_rig[4];
        Float_t kf_nchi[4][2];
        Float_t kf_cpu[4];
        
        Short_t hc[4];
        Float_t hc_rig[4];
        Float_t hc_nchi[4][2];
        Float_t hc_cpu[4];
        
        Short_t new_hc[4];
        Float_t new_hc_rig[4];
        Float_t new_hc_nchi[4][2];
        
        Float_t new_tdQ;
        Float_t new_tdQ_nchi;
        
        Float_t new_tfQ;
        Float_t new_tfQ_nchi;
        
        Float_t new_tkQ;
        Float_t new_tkQ_nchi;
};


#endif // __Variables_H__
