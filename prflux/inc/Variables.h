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
};


class VarsFLX : public Vars {
    public :
        VarsFLX() : Vars() { init(); }
        ~VarsFLX() {}
        
        void init() {
            init_base();

            trg = 0;

            tdllr = 0;
            tkllr = 0;
            
            tkL1 = 0;
            tkL9 = 0;
            tkL2 = 0;
           
            sign  = 0;
            rig   = 0;
            rigc  = 0;
            lxfs  = 0;
            lyfs  = 0;
            rigIn = 0;
            rigL1 = 0;
            rigL9 = 0;
            
            signC  = 0;
            rigC   = 0;
            rigcC  = 0;
            lxfsC  = 0;
            lyfsC  = 0;
            rigInC = 0;
            rigL1C = 0;
            rigL9C = 0;

            tf = 0;
            tf_sqrm = 0;
            tf_bta = 0;
            tf_rig = 0;
            tf_lb = 0;
            tf_lx = 0;
            tf_ly = 0;
            
            rh = 0;
            rh_sqrm = 0;
            rh_bta = 0;
            rh_rig = 0;
            rh_lb = 0;
            rh_lx = 0;
            rh_ly = 0;

            ecal = 0;
            engE = 0;
            engD = 0;
            engH = 0;
            apxL = -1;
            mvaBDT = -2;
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("trg"  , &trg);
            tree->Branch("tdllr", &tdllr);
            tree->Branch("tkllr", &tkllr);
            
            tree->Branch("tkL1", &tkL1);
            tree->Branch("tkL9", &tkL9);
            tree->Branch("tkL2", &tkL2);

            tree->Branch("sign" , &sign);
            tree->Branch("rig"  , &rig);
            tree->Branch("rigc" , &rigc);
            tree->Branch("lxfs" , &lxfs);
            tree->Branch("lyfs" , &lyfs);
            tree->Branch("rigIn", &rigIn);
            tree->Branch("rigL1", &rigL1);
            tree->Branch("rigL9", &rigL9);

            tree->Branch("signC" , &signC);
            tree->Branch("rigC"  , &rigC);
            tree->Branch("rigcC" , &rigcC);
            tree->Branch("lxfsC" , &lxfsC);
            tree->Branch("lyfsC" , &lyfsC);
            tree->Branch("rigInC", &rigInC);
            tree->Branch("rigL1C", &rigL1C);
            tree->Branch("rigL9C", &rigL9C);

            tree->Branch("tf"     , &tf);
            tree->Branch("tf_sqrm", &tf_sqrm);
            tree->Branch("tf_bta" , &tf_bta);
            tree->Branch("tf_rig" , &tf_rig);
            tree->Branch("tf_lb"  , &tf_lb);
            tree->Branch("tf_lx"  , &tf_lx);
            tree->Branch("tf_ly"  , &tf_ly);

            tree->Branch("rh"     , &rh);
            tree->Branch("rh_sqrm", &rh_sqrm);
            tree->Branch("rh_bta" , &rh_bta);
            tree->Branch("rh_rig" , &rh_rig);
            tree->Branch("rh_lb"  , &rh_lb);
            tree->Branch("rh_lx"  , &rh_lx);
            tree->Branch("rh_ly"  , &rh_ly);

            tree->Branch("ecal"  , &ecal);
            tree->Branch("engE"  , &engE);
            tree->Branch("engD"  , &engD);
            tree->Branch("engH"  , &engH);
            tree->Branch("apxL"  , &apxL);
            tree->Branch("mvaBDT", &mvaBDT);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("trg"  , &trg);
            chain->SetBranchAddress("tdllr", &tdllr);
            chain->SetBranchAddress("tkllr", &tkllr);
            
            chain->SetBranchAddress("tkL1", &tkL1);
            chain->SetBranchAddress("tkL9", &tkL9);
            chain->SetBranchAddress("tkL2", &tkL2);

            chain->SetBranchAddress("sign" , &sign);
            chain->SetBranchAddress("rig"  , &rig);
            chain->SetBranchAddress("rigc" , &rigc);
            chain->SetBranchAddress("lxfs" , &lxfs);
            chain->SetBranchAddress("lyfs" , &lyfs);
            chain->SetBranchAddress("rigIn", &rigIn);
            chain->SetBranchAddress("rigL1", &rigL1);
            chain->SetBranchAddress("rigL9", &rigL9);

            chain->SetBranchAddress("signC" , &signC);
            chain->SetBranchAddress("rigC"  , &rigC);
            chain->SetBranchAddress("rigcC" , &rigcC);
            chain->SetBranchAddress("lxfsC" , &lxfsC);
            chain->SetBranchAddress("lyfsC" , &lyfsC);
            chain->SetBranchAddress("rigInC", &rigInC);
            chain->SetBranchAddress("rigL1C", &rigL1C);
            chain->SetBranchAddress("rigL9C", &rigL9C);

            chain->SetBranchAddress("tf"     , &tf);
            chain->SetBranchAddress("tf_sqrm", &tf_sqrm);
            chain->SetBranchAddress("tf_bta" , &tf_bta);
            chain->SetBranchAddress("tf_rig" , &tf_rig);
            chain->SetBranchAddress("tf_lb"  , &tf_lb);
            chain->SetBranchAddress("tf_lx"  , &tf_lx);
            chain->SetBranchAddress("tf_ly"  , &tf_ly);

            chain->SetBranchAddress("rh"     , &rh);
            chain->SetBranchAddress("rh_sqrm", &rh_sqrm);
            chain->SetBranchAddress("rh_bta" , &rh_bta);
            chain->SetBranchAddress("rh_rig" , &rh_rig);
            chain->SetBranchAddress("rh_lb"  , &rh_lb);
            chain->SetBranchAddress("rh_lx"  , &rh_lx);
            chain->SetBranchAddress("rh_ly"  , &rh_ly);

            chain->SetBranchAddress("ecal"  , &ecal);
            chain->SetBranchAddress("engE"  , &engE);
            chain->SetBranchAddress("engD"  , &engD);
            chain->SetBranchAddress("engH"  , &engH);
            chain->SetBranchAddress("apxL"  , &apxL);
            chain->SetBranchAddress("mvaBDT", &mvaBDT);
        }

#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t trg;

        Float_t tdllr;
        Float_t tkllr;
        
        Short_t tkL1;
        Short_t tkL9;
        Short_t tkL2;
    
        Short_t sign;
        Float_t rig;
        Float_t rigc;
        Float_t lxfs;
        Float_t lyfs;
        Float_t rigIn;
        Float_t rigL1;
        Float_t rigL9;
        
        Short_t signC;
        Float_t rigC;
        Float_t rigcC;
        Float_t lxfsC;
        Float_t lyfsC;
        Float_t rigInC;
        Float_t rigL1C;
        Float_t rigL9C;

        Short_t tf;
        Float_t tf_sqrm;
        Float_t tf_bta;
        Float_t tf_rig;
        Float_t tf_lb;
        Float_t tf_lx;
        Float_t tf_ly;
        
        Short_t rh;
        Float_t rh_sqrm;
        Float_t rh_bta;
        Float_t rh_rig;
        Float_t rh_lb;
        Float_t rh_lx;
        Float_t rh_ly;
        
        Short_t ecal;
        Float_t engE;
        Float_t engD;
        Float_t engH;
        Short_t apxL;
        Float_t mvaBDT;
};


class VarsLTF : public Vars {
    public :
        VarsLTF() : Vars() { init(); }
        ~VarsLTF() {}
        
        void init() {
            init_base();

            sign    = 0;
            rig     = 0;
            bta     = 0;
            sqrm    = 0;
            llr     = 0;
            
            lxin    = 0;
            lyin    = 0;
            
            rich    = 0;
            rich_pr = 0;
            rich_el = 0;
            nhinn   = 0;
            nhout   = 0;

            std::fill_n(nvtxx, 2, 0);
            std::fill_n(nvtxy, 2, 0);
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("sign",    &sign);
            tree->Branch("rig",     &rig);
            tree->Branch("bta",     &bta);
            tree->Branch("sqrm",    &sqrm);
            tree->Branch("llr",     &llr);
            
            tree->Branch("lxin",    &lxin);
            tree->Branch("lyin",    &lyin);
            
            tree->Branch("rich",    &rich);
            tree->Branch("rich_pr", &rich_pr);
            tree->Branch("rich_el", &rich_el);
            tree->Branch("nvtxx",   nvtxx,  "nvtxx[2]/S");
            tree->Branch("nvtxy",   nvtxy,  "nvtxy[2]/S");
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("sign",    &sign);
            chain->SetBranchAddress("rig",     &rig);
            chain->SetBranchAddress("bta",     &bta);
            chain->SetBranchAddress("sqrm",    &sqrm);
            chain->SetBranchAddress("llr",     &llr);
            
            chain->SetBranchAddress("lxin",    &lxin);
            chain->SetBranchAddress("lyin",    &lyin);
            
            chain->SetBranchAddress("rich",    &rich);
            chain->SetBranchAddress("rich_pr", &rich_pr);
            chain->SetBranchAddress("rich_el", &rich_el);
            chain->SetBranchAddress("nhinn",   &nhinn);
            chain->SetBranchAddress("nhout",   &nhout);
            
            chain->SetBranchAddress("nvtxx",   nvtxx);
            chain->SetBranchAddress("nvtxy",   nvtxy);
        }

#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t sign;
        Float_t rig;
        Float_t bta;
        Float_t sqrm;
        Float_t llr;
        
        Float_t lxin;
        Float_t lyin;
        
        Short_t rich;
        Short_t rich_pr;
        Short_t rich_el;
        Short_t nhinn;
        Short_t nhout;
        
        Short_t nvtxx[2];
        Short_t nvtxy[2];
};


class VarsLRH : public Vars {
    public :
        VarsLRH() : Vars() { init(); }
        ~VarsLRH() {}
        
        void init() {
            init_base();
            
            sign  = 0;
            rig   = 0;
            bta   = 0;
            sqrm  = 0;
            llr   = 0;
            
            lxin  = 0;
            lyin  = 0;
        }

        void branch() {
            if (!tree) return;
            tree->Branch("sign",  &sign);
            tree->Branch("rig",   &rig);
            tree->Branch("bta",   &bta);
            tree->Branch("sqrm",  &sqrm);
            tree->Branch("llr",   &llr);
            
            tree->Branch("lxin",  &lxin);
            tree->Branch("lyin",  &lyin);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("sign",  &sign);
            chain->SetBranchAddress("rig",   &rig);
            chain->SetBranchAddress("bta",   &bta);
            chain->SetBranchAddress("sqrm",  &sqrm);
            chain->SetBranchAddress("llr",   &llr);
            
            chain->SetBranchAddress("lxin",  &lxin);
            chain->SetBranchAddress("lyin",  &lyin);
        }
        
#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t sign;
        Float_t rig;
        Float_t bta;
        Float_t sqrm;
        Float_t llr;
        
        Float_t lxin;
        Float_t lyin;
};


class VarsIIN : public Vars {
    public :
        VarsIIN() : Vars() { init(); }
        ~VarsIIN() {}
        
        void init() {
            init_base();

            sign   = 0;
            rig    = 0;
            llr    = 0;
            ecal   = 0;
            engE   = 0;
            engD   = 0;
            mvaBDT = -2;
            
            lxin  = 0;
            lyin  = 0;
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("sign",   &sign);
            tree->Branch("rig",    &rig);
            tree->Branch("llr",    &llr);
            tree->Branch("ecal",   &ecal);
            tree->Branch("engE",   &engE);
            tree->Branch("engD",   &engD);
            tree->Branch("mvaBDT", &mvaBDT);
            
            tree->Branch("lxin",  &lxin);
            tree->Branch("lyin",  &lyin);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("sign",   &sign);
            chain->SetBranchAddress("rig",    &rig);
            chain->SetBranchAddress("llr",    &llr);
            chain->SetBranchAddress("ecal",   &ecal);
            chain->SetBranchAddress("engE",   &engE);
            chain->SetBranchAddress("engD",   &engD);
            chain->SetBranchAddress("mvaBDT", &mvaBDT);
            
            chain->SetBranchAddress("lxin",  &lxin);
            chain->SetBranchAddress("lyin",  &lyin);
        }
        
#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t sign;
        Float_t rig;
        Float_t llr;
        Short_t ecal;
        Float_t engE;
        Float_t engD;
        Float_t mvaBDT;
        
        Float_t lxin;
        Float_t lyin;
};


class VarsIEX : public Vars {
    public :
        VarsIEX() : Vars() { init(); }
        ~VarsIEX() {}
        
        void init() {
            init_base();

            sign   = 0;
            rig    = 0;
            llr    = 0;
            ecal   = 0;
            engE   = 0;
            engD   = 0;
            mvaBDT = -2;
            
            patt   = 0;
            lxin   = 0;
            lyin   = 0;
            lxex   = 0;
            lyex   = 0;
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("sign",   &sign);
            tree->Branch("rig",    &rig);
            tree->Branch("llr",    &llr);
            tree->Branch("ecal",   &ecal);
            tree->Branch("engE",   &engE);
            tree->Branch("engD",   &engD);
            tree->Branch("mvaBDT", &mvaBDT);
            
            tree->Branch("patt",   &patt);
            tree->Branch("lxin",   &lxin);
            tree->Branch("lyin",   &lyin);
            tree->Branch("lxex",   &lxex);
            tree->Branch("lyex",   &lyex);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("sign",   &sign);
            chain->SetBranchAddress("rig",    &rig);
            chain->SetBranchAddress("llr",    &llr);
            chain->SetBranchAddress("ecal",   &ecal);
            chain->SetBranchAddress("engE",   &engE);
            chain->SetBranchAddress("engD",   &engD);
            chain->SetBranchAddress("mvaBDT", &mvaBDT);
            
            chain->SetBranchAddress("patt",   &patt);
            chain->SetBranchAddress("lxin",   &lxin);
            chain->SetBranchAddress("lyin",   &lyin);
            chain->SetBranchAddress("lxex",   &lxex);
            chain->SetBranchAddress("lyex",   &lyex);
        }
        
#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t sign;
        Float_t rig;
        Float_t llr;
        Short_t ecal;
        Float_t engE;
        Float_t engD;
        Float_t mvaBDT;
        
        Short_t patt;
        Float_t lxin;
        Float_t lyin;
        Float_t lxex;
        Float_t lyex;
};


class VarsHEX : public Vars {
    public :
        VarsHEX() : Vars() { init(); }
        ~VarsHEX() {}

        void init() {
            init_base();

            sign   = 0;
            rig    = 0;
            llr    = 0;
            ecal   = 0;
            engE   = 0;
            engD   = 0;
            mvaBDT = -2;
            
            patt   = 0;
            lxin   = 0;
            lyin   = 0;
            lxex   = 0;
            lyex   = 0;
            
            nlxin  = 0;
            nlyin  = 0;
            nlxex  = 0;
            nlyex  = 0;
            
            lxinC  = 0;
            lyinC  = 0;
            lxexC  = 0;
            lyexC  = 0;
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("sign",   &sign);
            tree->Branch("rig",    &rig);
            tree->Branch("llr",    &llr);
            tree->Branch("ecal",   &ecal);
            tree->Branch("engE",   &engE);
            tree->Branch("engD",   &engD);
            tree->Branch("mvaBDT", &mvaBDT);
            
            tree->Branch("patt",   &patt);
            tree->Branch("lxin",   &lxin);
            tree->Branch("lyin",   &lyin);
            tree->Branch("lxex",   &lxex);
            tree->Branch("lyex",   &lyex);
            
            tree->Branch("nlxin",  &nlxin);
            tree->Branch("nlyin",  &nlyin);
            tree->Branch("nlxex",  &nlxex);
            tree->Branch("nlyex",  &nlyex);
            
            tree->Branch("lxinC",  &lxinC);
            tree->Branch("lyinC",  &lyinC);
            tree->Branch("lxexC",  &lxexC);
            tree->Branch("lyexC",  &lyexC);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("sign",   &sign);
            chain->SetBranchAddress("rig",    &rig);
            chain->SetBranchAddress("llr",    &llr);
            chain->SetBranchAddress("ecal",   &ecal);
            chain->SetBranchAddress("engE",   &engE);
            chain->SetBranchAddress("engD",   &engD);
            chain->SetBranchAddress("mvaBDT", &mvaBDT);
            
            chain->SetBranchAddress("patt",   &patt);
            chain->SetBranchAddress("lxin",   &lxin);
            chain->SetBranchAddress("lyin",   &lyin);
            chain->SetBranchAddress("lxex",   &lxex);
            chain->SetBranchAddress("lyex",   &lyex);
         
            chain->SetBranchAddress("nlxin",  &nlxin);
            chain->SetBranchAddress("nlyin",  &nlyin);
            chain->SetBranchAddress("nlxex",  &nlxex);
            chain->SetBranchAddress("nlyex",  &nlyex);
            
            chain->SetBranchAddress("lxinC",  &lxinC);
            chain->SetBranchAddress("lyinC",  &lyinC);
            chain->SetBranchAddress("lxexC",  &lxexC);
            chain->SetBranchAddress("lyexC",  &lyexC);
        }
        
#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
            tmva_reader->AddVariable("lxin",  &lxin);
            tmva_reader->AddVariable("lyin",  &lyin);
            tmva_reader->AddVariable("lxex",  &lxex);
            tmva_reader->AddVariable("lyex",  &lyex);
            
            tmva_reader->AddVariable("nlxin", &nlxin);
            tmva_reader->AddVariable("nlyin", &nlyin);
            tmva_reader->AddVariable("nlxex", &nlxex);
            tmva_reader->AddVariable("nlyex", &nlyex);
            
            tmva_reader->AddVariable("lxinC",  &lxinC);
            tmva_reader->AddVariable("lyinC",  &lyinC);
            tmva_reader->AddVariable("lxexC",  &lxexC);
            tmva_reader->AddVariable("lyexC",  &lyexC);
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t sign;
        Float_t rig;
        Float_t llr;
        Short_t ecal;
        Float_t engE;
        Float_t engD;
        Float_t mvaBDT;
        
        Short_t patt;
        Float_t lxin;
        Float_t lyin;
        Float_t lxex;
        Float_t lyex;
        
        Float_t nlxin;
        Float_t nlyin;
        Float_t nlxex;
        Float_t nlyex;
        
        Float_t lxinC;
        Float_t lyinC;
        Float_t lxexC;
        Float_t lyexC;
};


class VarsHFS : public Vars {
    public :
        VarsHFS() : Vars() { init(); }
        ~VarsHFS() {}

        void init() {
            init_base();

            sign  = 0;
            rig   = 0;
            llr   = 0;
            ecal   = 0;
            engE   = 0;
            engD   = 0;
            mvaBDT = -2;
            
            lxin  = 0;
            lyin  = 0;
            lxl1  = 0;
            lyl1  = 0;
            lxl9  = 0;
            lyl9  = 0;
            lxfs  = 0;
            lyfs  = 0;
            
            nlxin = 0;
            nlyin = 0;
            nlxl1 = 0;
            nlyl1 = 0;
            nlxl9 = 0;
            nlyl9 = 0;
            nlxfs = 0;
            nlyfs = 0;
            
            lxinC  = 0;
            lyinC  = 0;
            lxl1C  = 0;
            lyl1C  = 0;
            lxl9C  = 0;
            lyl9C  = 0;
            lxfsC  = 0;
            lyfsC  = 0;

            lrin = 0;
            lrfs = 0;
            lrinC = 0;
            lrfsC = 0;
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("sign",  &sign);
            tree->Branch("rig",   &rig);
            tree->Branch("llr",   &llr);
            tree->Branch("ecal",   &ecal);
            tree->Branch("engE",   &engE);
            tree->Branch("engD",   &engD);
            tree->Branch("mvaBDT", &mvaBDT);
            
            tree->Branch("lxin",  &lxin);
            tree->Branch("lyin",  &lyin);
            tree->Branch("lxl1",  &lxl1);
            tree->Branch("lyl1",  &lyl1);
            tree->Branch("lxl9",  &lxl9);
            tree->Branch("lyl9",  &lyl9);
            tree->Branch("lxfs",  &lxfs);
            tree->Branch("lyfs",  &lyfs);
            
            tree->Branch("nlxin", &nlxin);
            tree->Branch("nlyin", &nlyin);
            tree->Branch("nlxl1", &nlxl1);
            tree->Branch("nlyl1", &nlyl1);
            tree->Branch("nlxl9", &nlxl9);
            tree->Branch("nlyl9", &nlyl9);
            tree->Branch("nlxfs", &nlxfs);
            tree->Branch("nlyfs", &nlyfs);
            
            tree->Branch("lxinC",  &lxinC);
            tree->Branch("lyinC",  &lyinC);
            tree->Branch("lxl1C",  &lxl1C);
            tree->Branch("lyl1C",  &lyl1C);
            tree->Branch("lxl9C",  &lxl9C);
            tree->Branch("lyl9C",  &lyl9C);
            tree->Branch("lxfsC",  &lxfsC);
            tree->Branch("lyfsC",  &lyfsC);
            
            tree->Branch("lrin",  &lrin);
            tree->Branch("lrfs",  &lrfs);
            tree->Branch("lrinC",  &lrinC);
            tree->Branch("lrfsC",  &lrfsC);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("sign",  &sign);
            chain->SetBranchAddress("rig",   &rig);
            chain->SetBranchAddress("llr",   &llr);
            chain->SetBranchAddress("ecal",   &ecal);
            chain->SetBranchAddress("engE",   &engE);
            chain->SetBranchAddress("engD",   &engD);
            chain->SetBranchAddress("mvaBDT", &mvaBDT);
            
            chain->SetBranchAddress("lxin",  &lxin);
            chain->SetBranchAddress("lyin",  &lyin);
            chain->SetBranchAddress("lxl1",  &lxl1);
            chain->SetBranchAddress("lyl1",  &lyl1);
            chain->SetBranchAddress("lxl9",  &lxl9);
            chain->SetBranchAddress("lyl9",  &lyl9);
            chain->SetBranchAddress("lxfs",  &lxfs);
            chain->SetBranchAddress("lyfs",  &lyfs);
            
            chain->SetBranchAddress("nlxin", &nlxin);
            chain->SetBranchAddress("nlyin", &nlyin);
            chain->SetBranchAddress("nlxl1", &nlxl1);
            chain->SetBranchAddress("nlyl1", &nlyl1);
            chain->SetBranchAddress("nlxl9", &nlxl9);
            chain->SetBranchAddress("nlyl9", &nlyl9);
            chain->SetBranchAddress("nlxfs", &nlxfs);
            chain->SetBranchAddress("nlyfs", &nlyfs);
            
            chain->SetBranchAddress("lxinC",  &lxinC);
            chain->SetBranchAddress("lyinC",  &lyinC);
            chain->SetBranchAddress("lxl1C",  &lxl1C);
            chain->SetBranchAddress("lyl1C",  &lyl1C);
            chain->SetBranchAddress("lxl9C",  &lxl9C);
            chain->SetBranchAddress("lyl9C",  &lyl9C);
            chain->SetBranchAddress("lxfsC",  &lxfsC);
            chain->SetBranchAddress("lyfsC",  &lyfsC);
            
            chain->SetBranchAddress("lrin",  &lrin);
            chain->SetBranchAddress("lrfs",  &lrfs);
            chain->SetBranchAddress("lrinC",  &lrinC);
            chain->SetBranchAddress("lrfsC",  &lrfsC);
        }
        
#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
            tmva_reader->AddVariable("lxin",  &lxin);
            tmva_reader->AddVariable("lyin",  &lyin);
            tmva_reader->AddVariable("lxl1",  &lxl1);
            tmva_reader->AddVariable("lyl1",  &lyl1);
            tmva_reader->AddVariable("lxl9",  &lxl9);
            tmva_reader->AddVariable("lyl9",  &lyl9);
            tmva_reader->AddVariable("lxfs",  &lxfs);
            tmva_reader->AddVariable("lyfs",  &lyfs);
            
            tmva_reader->AddVariable("nlxin", &nlxin);
            tmva_reader->AddVariable("nlyin", &nlyin);
            tmva_reader->AddVariable("nlxl1", &nlxl1);
            tmva_reader->AddVariable("nlyl1", &nlyl1);
            tmva_reader->AddVariable("nlxl9", &nlxl9);
            tmva_reader->AddVariable("nlyl9", &nlyl9);
            tmva_reader->AddVariable("nlxfs", &nlxfs);
            tmva_reader->AddVariable("nlyfs", &nlyfs);
            
            tmva_reader->AddVariable("lxinC",  &lxinC);
            tmva_reader->AddVariable("lyinC",  &lyinC);
            tmva_reader->AddVariable("lxl1C",  &lxl1C);
            tmva_reader->AddVariable("lyl1C",  &lyl1C);
            tmva_reader->AddVariable("lxl9C",  &lxl9C);
            tmva_reader->AddVariable("lyl9C",  &lyl9C);
            tmva_reader->AddVariable("lxfsC",  &lxfsC);
            tmva_reader->AddVariable("lyfsC",  &lyfsC);
            
            tmva_reader->AddVariable("lrin",  &lrin);
            tmva_reader->AddVariable("lrfs",  &lrfs);
            tmva_reader->AddVariable("lrinC",  &lrinC);
            tmva_reader->AddVariable("lrfsC",  &lrfsC);
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t sign;
        Float_t rig;
        Float_t llr;
        Short_t ecal;
        Float_t engE;
        Float_t engD;
        Float_t mvaBDT;
        
        Float_t lxin;
        Float_t lyin;
        Float_t lxl1;
        Float_t lyl1;
        Float_t lxl9;
        Float_t lyl9;
        Float_t lxfs;
        Float_t lyfs;
        
        Float_t nlxin;
        Float_t nlyin;
        Float_t nlxl1;
        Float_t nlyl1;
        Float_t nlxl9;
        Float_t nlyl9;
        Float_t nlxfs;
        Float_t nlyfs;

        Float_t lxinC;
        Float_t lyinC;
        Float_t lxl1C;
        Float_t lyl1C;
        Float_t lxl9C;
        Float_t lyl9C;
        Float_t lxfsC;
        Float_t lyfsC;
        
        Float_t lrin;
        Float_t lrfs;
        Float_t lrinC;
        Float_t lrfsC;
};

#endif // __Variables_H__
