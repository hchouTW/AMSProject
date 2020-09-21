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
            mc_z = 0;
            mc_rig = 0;
            mc_igb = 0;
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
            tree->Branch("mc_z"  , &mc_z);
            tree->Branch("mc_rig", &mc_rig);
            tree->Branch("mc_igb", &mc_igb);
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
            chain->SetBranchAddress("mc_z"  , &mc_z);
            chain->SetBranchAddress("mc_rig", &mc_rig);
            chain->SetBranchAddress("mc_igb", &mc_igb);
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
        Short_t mc_z;
        Float_t mc_rig;
        Float_t mc_igb;
        Float_t mc_w10;
        Float_t mc_w27;

        Float_t cfsec;
        Float_t cfevt;
        Float_t lv;
        
        Short_t trg;
};


class VarsDB : public Vars {
    public :
        VarsDB() : Vars() { init(); }
        ~VarsDB() {}
        
        void init() {
            init_base();

            L2  = 0;
            nhx = 0;
            nhy = 0;
            
            std::fill_n(ck        , 4, 0);
            std::fill_n(ck_rig    , 4, 0);
            std::fill_n(ck_nchi[0], 4*2, 0);
           
            tfN = 0;
            tfQ = 0;
            std::fill_n(tfL , 4, 0);
            std::fill_n(tfLQ, 4, -1);
            
            tfQ_nchi_c = 0;
            tfQ_nchi_t = 0;
            
            tkQhl = 0;
            tkQyj = 0;
            
            tkQ = 0;
            std::fill_n(tkL  , 9, 0);
            std::fill_n(tkLQ , 9, -1);
            std::fill_n(tkLQx, 9, -1);
            std::fill_n(tkLQy, 9, -1);
            
            std::fill_n(tkHitL, 9, 0);
            std::fill_n(tkHitX, 9, 0);
            std::fill_n(tkHitY, 9, 0);
            std::fill_n(tkHitZ, 9, 0);
            
            new_tfQ      = 0;
            new_tfQ_nchi = 0;
            new_tfQ_nhit = 0;
            new_tfQ_ncls = 0;
            
            new_tkQ      = 0;
            new_tkQ_nchi = 0;
            new_tkQ_nhit = 0;
            new_tkQ_ncls = 0;
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("L2" , &L2);
            tree->Branch("nhx", &nhx);
            tree->Branch("nhy", &nhy);
            
            tree->Branch("ck"     , ck     , "ck[4]/S");
            tree->Branch("ck_rig" , ck_rig , "ck_rig[4]/F");
            tree->Branch("ck_nchi", ck_nchi, "ck_nchi[4][2]/F");
            
            tree->Branch("tfN" , &tfN);
            tree->Branch("tfQ" , &tfQ);
            tree->Branch("tfL" ,  tfL,  "tfL[4]/S");
            tree->Branch("tfLQ",  tfLQ, "tfLQ[4]/F");
            
            tree->Branch("tfQ_nchi_c" , &tfQ_nchi_c);
            tree->Branch("tfQ_nchi_t" , &tfQ_nchi_t);
            
            tree->Branch("tkQhl"  , &tkQhl);
            tree->Branch("tkQyj"  , &tkQyj);
            
            tree->Branch("tkQ"  , &tkQ);
            tree->Branch("tkL"  ,  tkL  , "tkL[9]/S");
            tree->Branch("tkLQ" ,  tkLQ , "tkLQ[9]/F");
            tree->Branch("tkLQx",  tkLQx, "tkLQx[9]/F");
            tree->Branch("tkLQy",  tkLQy, "tkLQy[9]/F");
            
            tree->Branch("tkHitL", tkHitL, "tkHitL[9]/S");
            tree->Branch("tkHitX", tkHitX, "tkHitX[9]/D");
            tree->Branch("tkHitY", tkHitY, "tkHitY[9]/D");
            tree->Branch("tkHitZ", tkHitZ, "tkHitZ[9]/D");
            
            tree->Branch("new_tfQ"     , &new_tfQ);
            tree->Branch("new_tfQ_nchi", &new_tfQ_nchi);
            tree->Branch("new_tfQ_nhit", &new_tfQ_nhit);
            tree->Branch("new_tfQ_ncls", &new_tfQ_ncls);
            
            tree->Branch("new_tkQ"     , &new_tkQ);
            tree->Branch("new_tkQ_nchi", &new_tkQ_nchi);
            tree->Branch("new_tkQ_nhit", &new_tkQ_nhit);
            tree->Branch("new_tkQ_ncls", &new_tkQ_ncls);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("L2" , &L2);
            chain->SetBranchAddress("nhx", &nhx);
            chain->SetBranchAddress("nhy", &nhy);
            
            chain->SetBranchAddress("ck"     , ck     );
            chain->SetBranchAddress("ck_rig" , ck_rig );
            chain->SetBranchAddress("ck_nchi", ck_nchi);
            
            chain->SetBranchAddress("tfN" , &tfN);
            chain->SetBranchAddress("tfQ" , &tfQ);
            chain->SetBranchAddress("tfL" ,  tfL);
            chain->SetBranchAddress("tfLQ",  tfLQ);
            
            chain->SetBranchAddress("tfQ_nchi_c" , &tfQ_nchi_c);
            chain->SetBranchAddress("tfQ_nchi_t" , &tfQ_nchi_t);
            
            chain->SetBranchAddress("tkQhl", &tkQhl);
            chain->SetBranchAddress("tkQyj", &tkQyj);
            
            chain->SetBranchAddress("tkQ"  , &tkQ);
            chain->SetBranchAddress("tkL"  ,  tkL);
            chain->SetBranchAddress("tkLQ" ,  tkLQ);
            chain->SetBranchAddress("tkLQx",  tkLQx);
            chain->SetBranchAddress("tkLQy",  tkLQy);
            
            chain->SetBranchAddress("tkHitL", tkHitL);
            chain->SetBranchAddress("tkHitX", tkHitX);
            chain->SetBranchAddress("tkHitY", tkHitY);
            chain->SetBranchAddress("tkHitZ", tkHitZ);
            
            chain->SetBranchAddress("new_tfQ"     , &new_tfQ);
            chain->SetBranchAddress("new_tfQ_nchi", &new_tfQ_nchi);
            chain->SetBranchAddress("new_tfQ_nhit", &new_tfQ_nhit);
            chain->SetBranchAddress("new_tfQ_ncls", &new_tfQ_ncls);
            
            chain->SetBranchAddress("new_tkQ"     , &new_tkQ);
            chain->SetBranchAddress("new_tkQ_nchi", &new_tkQ_nchi);
            chain->SetBranchAddress("new_tkQ_nhit", &new_tkQ_nhit);
            chain->SetBranchAddress("new_tkQ_ncls", &new_tkQ_ncls);
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
        
        Short_t ck[4];
        Float_t ck_rig[4];
        Float_t ck_nchi[4][2];
        
        Short_t tfN;
        Float_t tfQ;
        Short_t tfL[4];
        Float_t tfLQ[4];
        
        Float_t tfQ_nchi_c;
        Float_t tfQ_nchi_t;
        
        Float_t tkQhl;
        Float_t tkQyj;

        Float_t tkQ;
        Short_t tkL[9];
        Float_t tkLQ[9];
        Float_t tkLQx[9];
        Float_t tkLQy[9];
        
        Short_t  tkHitL[9];
        Double_t tkHitX[9];
        Double_t tkHitY[9];
        Double_t tkHitZ[9];

        Float_t new_tfQ;
        Float_t new_tfQ_nchi;
        Short_t new_tfQ_nhit;
        Short_t new_tfQ_ncls;
        
        Float_t new_tkQ;
        Float_t new_tkQ_nchi;
        Short_t new_tkQ_nhit;
        Short_t new_tkQ_ncls;
};


#endif // __Variables_H__
