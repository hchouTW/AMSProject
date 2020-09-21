#ifndef __Variables_H__
#define __Variables_H__

#include "TTree.h"
#include "TChain.h"

class Vars {
    public :
        Vars() { init_base(); tree = nullptr; chain = nullptr; }
        ~Vars() {}
        
        void init_base() {
            run = 0;
            evt = 0;
            ut  = 0;
            wgt = 1;

            mc = false;
            mc_rig = 0;
            mc_w10 = 1;
            mc_w27 = 1;

            cfsec = 0;
            cfevt = 0;
            lv    = 1;
        }

        void set_tree_base(const char* name, const char* title = "") {
            if (tree) return;
            tree = new TTree(name, title);
        }
        
        void set_tree_base(TChain* ref) {
            if (chain) return;
            chain = ref;
        }

        void branch_base() {
            if (!tree) return;

            tree->Branch("run", &run);
            tree->Branch("evt", &evt);
            tree->Branch("ut",  &ut);
            tree->Branch("wgt", &wgt);
            
            tree->Branch("mc", &mc);
            tree->Branch("mc_rig", &mc_rig);
            tree->Branch("mc_w10", &mc_w10);
            tree->Branch("mc_w27", &mc_w27);
            
            tree->Branch("cfsec", &cfsec);
            tree->Branch("cfevt", &cfevt);
            tree->Branch("lv",    &lv);
        }

        void set_branch_address_base() {
            if (!chain) return;
            chain->SetBranchAddress("run", &run);
            chain->SetBranchAddress("evt", &evt);
            chain->SetBranchAddress("ut",  &ut);
            chain->SetBranchAddress("wgt", &wgt);

            chain->SetBranchAddress("mc", &mc);
            chain->SetBranchAddress("mc_rig", &mc_rig);
            chain->SetBranchAddress("mc_w10", &mc_w10);
            chain->SetBranchAddress("mc_w27", &mc_w27);
            
            chain->SetBranchAddress("cfsec", &cfsec);
            chain->SetBranchAddress("cfevt", &cfevt);
            chain->SetBranchAddress("lv",    &lv);
        }
        
        virtual void branch() = 0;
        void set_tree(const char* name, const char* title = "") {
            set_tree_base(name, title);
            branch();
        }
        
        virtual void set_branch_address() = 0;
        void set_tree(TChain* ref) {
            set_tree_base(ref);
            set_branch_address();
        }
        
        void fill() {
            if (!tree) return;
            tree->Fill();
        }

    protected :
        TTree*  tree;
        TChain* chain;

    public :
        UInt_t  run;
        UInt_t  evt;
        UInt_t  ut;
        Float_t wgt;
       
        Bool_t  mc;
        Float_t mc_rig;
        Float_t mc_w10;
        Float_t mc_w27;

        Float_t cfsec;
        Float_t cfevt;
        Float_t lv;
};


class VarsL : public Vars {
    public :
        VarsL() : Vars() { init(); }
        ~VarsL() {}
        
        void init() {
            init_base();

            sign    = 0;
            rig     = 0;
            bta     = 0;
            sqrm    = 0;
            llr     = 0;
            rich    = 0;
            rich_pr = 0;
            rich_el = 0;
            tdvtx_x = 0;
            tdvtx_y = 0;
        }
        
        void branch() {
            if (!tree) return;
            branch_base();

            tree->Branch("sign",    &sign);
            tree->Branch("rig",     &rig);
            tree->Branch("bta",     &bta);
            tree->Branch("sqrm",    &sqrm);
            tree->Branch("llr",     &llr);
            tree->Branch("rich",    &rich);
            tree->Branch("rich_pr", &rich_pr);
            tree->Branch("rich_el", &rich_el);
            tree->Branch("tdvtx_x", &tdvtx_x);
            tree->Branch("tdvtx_y", &tdvtx_y);
        }
        
        void set_branch_address() {
            if (!chain) return;
            set_branch_address_base();

            chain->SetBranchAddress("sign",    &sign);
            chain->SetBranchAddress("rig",     &rig);
            chain->SetBranchAddress("bta",     &bta);
            chain->SetBranchAddress("sqrm",    &sqrm);
            chain->SetBranchAddress("llr",     &llr);
            chain->SetBranchAddress("rich",    &rich);
            chain->SetBranchAddress("rich_pr", &rich_pr);
            chain->SetBranchAddress("rich_el", &rich_el);
            chain->SetBranchAddress("tdvtx_x", &tdvtx_x);
            chain->SetBranchAddress("tdvtx_y", &tdvtx_y);
        }

    public :
        Short_t sign;
        Float_t rig;
        Float_t bta;
        Float_t sqrm;
        Float_t llr;
        Short_t rich;
        Short_t rich_pr;
        Short_t rich_el;
        Short_t tdvtx_x;
        Short_t tdvtx_y;
};


class VarsM : public Vars {
    public :
        VarsM() : Vars() { init(); }
        ~VarsM() {}
        
        void init() {
            init_base();

            sign = 0;
            rig  = 0;
            bta  = 0;
            sqrm = 0;
            llr  = 0;
        }

        void branch() {
            if (!tree) return;
            branch_base();

            tree->Branch("sign", &sign);
            tree->Branch("rig",  &rig);
            tree->Branch("bta",  &bta);
            tree->Branch("sqrm", &sqrm);
            tree->Branch("llr",  &llr);
        }
        
        void set_branch_address() {
            if (!tree) return;
            set_branch_address_base();

            tree->SetBranchAddress("sign", &sign);
            tree->SetBranchAddress("rig",  &rig);
            tree->SetBranchAddress("bta",  &bta);
            tree->SetBranchAddress("sqrm", &sqrm);
            tree->SetBranchAddress("llr",  &llr);
        }

    public :
        Short_t sign;
        Float_t rig;
        Float_t bta;
        Float_t sqrm;
        Float_t llr;
};


class VarsI : public Vars {
    public :
        VarsI() : Vars() { init(); }
        ~VarsI() {}
        
        void init() {
            init_base();

            sign   = 0;
            rig    = 0;
            llr    = 0;
            ecal   = 0;
            mvaBDT = 0;
        }
        
        void branch() {
            if (!tree) return;
            branch_base();

            tree->Branch("sign",   &sign);
            tree->Branch("rig",    &rig);
            tree->Branch("llr",    &llr);
            tree->Branch("ecal",   &ecal);
            tree->Branch("mvaBDT", &mvaBDT);
        }
        
        void set_branch_address() {
            if (!chain) return;
            set_branch_address_base();

            chain->SetBranchAddress("sign",   &sign);
            chain->SetBranchAddress("rig",    &rig);
            chain->SetBranchAddress("llr",    &llr);
            chain->SetBranchAddress("ecal",   &ecal);
            chain->SetBranchAddress("mvaBDT", &mvaBDT);
        }

    public :
        Short_t sign;
        Float_t rig;
        Float_t llr;
        Short_t ecal;
        Float_t mvaBDT;
};


class VarsL1 : public Vars {
    public :
        VarsL1() : Vars() { init(); }
        ~VarsL1() {}

        void init() {
            init_base();

            sign  = 0;
            rig   = 0;
            llr   = 0;
            lxin  = 0;
            lyin  = 0;
            lxl1  = 0;
            lyl1  = 0;
            lrl1  = 0;
            nrmlx = 0;
            nrmly = 0;
            nrmms = 0;
        }
        
        void branch() {
            if (!tree) return;
            branch_base();

            tree->Branch("sign",  &sign);
            tree->Branch("rig",   &rig);
            tree->Branch("llr",   &llr);
            tree->Branch("lxin",  &lxin);
            tree->Branch("lyin",  &lyin);
            tree->Branch("lxl1",  &lxl1);
            tree->Branch("lyl1",  &lyl1);
            tree->Branch("lrl1",  &lrl1);
            tree->Branch("nrmlx", &nrmlx);
            tree->Branch("nrmly", &nrmly);
            tree->Branch("nrmms", &nrmms);
        }
        
        void set_branch_address() {
            if (!chain) return;
            set_branch_address_base();

            chain->SetBranchAddress("sign",  &sign);
            chain->SetBranchAddress("rig",   &rig);
            chain->SetBranchAddress("llr",   &llr);
            chain->SetBranchAddress("lxin",  &lxin);
            chain->SetBranchAddress("lyin",  &lyin);
            chain->SetBranchAddress("lxl1",  &lxl1);
            chain->SetBranchAddress("lyl1",  &lyl1);
            chain->SetBranchAddress("lrl1",  &lrl1);
            chain->SetBranchAddress("nrmlx", &nrmlx);
            chain->SetBranchAddress("nrmly", &nrmly);
            chain->SetBranchAddress("nrmms", &nrmms);
        }

    public :
        Short_t sign;
        Float_t rig;
        Float_t llr;
        Float_t lxin;
        Float_t lyin;
        Float_t lxl1;
        Float_t lyl1;
        Float_t lrl1;
        Float_t nrmlx;
        Float_t nrmly;
        Float_t nrmms;
};


class VarsL9 : public Vars {
    public :
        VarsL9() : Vars() { init(); }
        ~VarsL9() {}

        void init() {
            init_base();

            sign  = 0;
            rig   = 0;
            llr   = 0;
            lxin  = 0;
            lyin  = 0;
            lxl9  = 0;
            lyl9  = 0;
            lrl9  = 0;
            nrmlx = 0;
            nrmly = 0;
            nrmms = 0;
        }
        
        void branch() {
            if (!tree) return;
            branch_base();

            tree->Branch("sign",  &sign);
            tree->Branch("rig",   &rig);
            tree->Branch("llr",   &llr);
            tree->Branch("lxin",  &lxin);
            tree->Branch("lyin",  &lyin);
            tree->Branch("lxl9",  &lxl9);
            tree->Branch("lyl9",  &lyl9);
            tree->Branch("lrl9",  &lrl9);
            tree->Branch("nrmlx", &nrmlx);
            tree->Branch("nrmly", &nrmly);
            tree->Branch("nrmms", &nrmms);
        }
        
        void set_branch_address() {
            if (!chain) return;
            set_branch_address_base();

            chain->SetBranchAddress("sign",  &sign);
            chain->SetBranchAddress("rig",   &rig);
            chain->SetBranchAddress("llr",   &llr);
            chain->SetBranchAddress("lxin",  &lxin);
            chain->SetBranchAddress("lyin",  &lyin);
            chain->SetBranchAddress("lxl9",  &lxl9);
            chain->SetBranchAddress("lyl9",  &lyl9);
            chain->SetBranchAddress("lrl9",  &lrl9);
            chain->SetBranchAddress("nrmlx", &nrmlx);
            chain->SetBranchAddress("nrmly", &nrmly);
            chain->SetBranchAddress("nrmms", &nrmms);
        }

    public :
        Short_t sign;
        Float_t rig;
        Float_t llr;
        Float_t lxin;
        Float_t lyin;
        Float_t lxl9;
        Float_t lyl9;
        Float_t lrl9;
        Float_t nrmlx;
        Float_t nrmly;
        Float_t nrmms;
};


class VarsFS : public Vars {
    public :
        VarsFS() : Vars() { init(); }
        ~VarsFS() {}

        void init() {
            init_base();

            sign  = 0;
            rig   = 0;
            llr   = 0;
            lxin  = 0;
            lyin  = 0;
            lxl1  = 0;
            lyl1  = 0;
            lxl9  = 0;
            lyl9  = 0;
            lxfs  = 0;
            lyfs  = 0;
            lrl1  = 0;
            lrl9  = 0;
            lrfs  = 0;
            nrmlx = 0;
            nrmly = 0;
            nrmms = 0;
        }
        
        void branch() {
            if (!tree) return;
            branch_base();

            tree->Branch("sign",  &sign);
            tree->Branch("rig",   &rig);
            tree->Branch("llr",   &llr);
            tree->Branch("lxin",  &lxin);
            tree->Branch("lyin",  &lyin);
            tree->Branch("lxl1",  &lxl1);
            tree->Branch("lyl1",  &lyl1);
            tree->Branch("lxl9",  &lxl9);
            tree->Branch("lyl9",  &lyl9);
            tree->Branch("lxfs",  &lxfs);
            tree->Branch("lyfs",  &lyfs);
            tree->Branch("lrl1",  &lrl1);
            tree->Branch("lrl9",  &lrl9);
            tree->Branch("lrfs",  &lrfs);
            tree->Branch("nrmlx", &nrmlx);
            tree->Branch("nrmly", &nrmly);
            tree->Branch("nrmms", &nrmms);
        }
        
        void set_branch_address() {
            if (!chain) return;
            set_branch_address_base();

            chain->SetBranchAddress("sign",  &sign);
            chain->SetBranchAddress("rig",   &rig);
            chain->SetBranchAddress("llr",   &llr);
            chain->SetBranchAddress("lxin",  &lxin);
            chain->SetBranchAddress("lyin",  &lyin);
            chain->SetBranchAddress("lxl1",  &lxl1);
            chain->SetBranchAddress("lyl1",  &lyl1);
            chain->SetBranchAddress("lxl9",  &lxl9);
            chain->SetBranchAddress("lyl9",  &lyl9);
            chain->SetBranchAddress("lxfs",  &lxfs);
            chain->SetBranchAddress("lyfs",  &lyfs);
            chain->SetBranchAddress("lrl1",  &lrl1);
            chain->SetBranchAddress("lrl9",  &lrl9);
            chain->SetBranchAddress("lrfs",  &lrfs);
            chain->SetBranchAddress("nrmlx", &nrmlx);
            chain->SetBranchAddress("nrmly", &nrmly);
            chain->SetBranchAddress("nrmms", &nrmms);
        }

    public :
        Short_t sign;
        Float_t rig;
        Float_t llr;
        Float_t lxin;
        Float_t lyin;
        Float_t lxl1;
        Float_t lyl1;
        Float_t lxl9;
        Float_t lyl9;
        Float_t lxfs;
        Float_t lyfs;
        Float_t lrl1;
        Float_t lrl9;
        Float_t lrfs;
        Float_t nrmlx;
        Float_t nrmly;
        Float_t nrmms;
};



















class Variables {
    public :
        Variables() { init(); tree = nullptr; }
        ~Variables() {}

        void init() {
            run = 0;
            ev  = 0;
            std::fill_n(pt, 4, 0);
            w10 = 1;
            w27 = 1;
            rig = 0;
            std::fill_n(chix, 4, 0);
            std::fill_n(chiy, 4, 0);
            std::fill_n(rvar, 3, 0);
            std::fill_n(norm, 3, 0);
        }

        void fill() {
            if (tree == nullptr) return;
            tree->Fill();
        }

        void branch(TTree* reftree) {
            if (reftree == nullptr) return;
            tree = reftree;
            tree->Branch("run", &run);
            tree->Branch("ev" , &ev);
            tree->Branch("pt" ,  pt, "pt[4]/S");
            tree->Branch("w10", &w10);
            tree->Branch("w27", &w27);
            tree->Branch("rig", &rig);
            tree->Branch("chix", chix, "chix[4]/F");
            tree->Branch("chiy", chiy, "chiy[4]/F");
            tree->Branch("rvar", rvar, "rvar[3]/F");
            tree->Branch("norm", norm, "norm[3]/F");
        }
       
    public :
        TTree*  tree;

        UInt_t  run;
        UInt_t  ev;
        Short_t pt[4];
        Float_t w10;
        Float_t w27;
        Float_t rig;
        Float_t chix[4];
        Float_t chiy[4];
        Float_t rvar[3];
        Float_t norm[3];
};

#endif // __Variables_H__
