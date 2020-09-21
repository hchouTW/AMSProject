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
            br  = 0;
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
            tree->Branch("br"    , &br);
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
            chain->SetBranchAddress("br"    , &br);
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
        UInt_t  br;
        Float_t wgt;
       
        Short_t mc;
        Float_t mc_rig;
        Float_t mc_w10;
        Float_t mc_w27;

        Float_t cfsec;
        Float_t cfevt;
        Float_t lv;
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
            
            ecal   = 0;
            mvaBDT = -2;
            
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
            
            tree->Branch("ecal"  , &ecal);
            tree->Branch("mvaBDT", &mvaBDT);
            
            tree->Branch("lxin",    &lxin);
            tree->Branch("lyin",    &lyin);
            
            tree->Branch("rich",    &rich);
            tree->Branch("rich_pr", &rich_pr);
            tree->Branch("rich_el", &rich_el);
            tree->Branch("nhinn",   &nhinn);
            tree->Branch("nhout",   &nhout);
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
            
            chain->SetBranchAddress("ecal"  , &ecal);
            chain->SetBranchAddress("mvaBDT", &mvaBDT);
            
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
        
        Short_t ecal;
        Float_t mvaBDT;
        
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
           
            ecal   = 0;
            mvaBDT = -2;

            lxin  = 0;
            lyin  = 0;
            
            nhinn   = 0;
            nhout   = 0;
        }

        void branch() {
            if (!tree) return;
            tree->Branch("sign",  &sign);
            tree->Branch("rig",   &rig);
            tree->Branch("bta",   &bta);
            tree->Branch("sqrm",  &sqrm);
            tree->Branch("llr",   &llr);
            
            tree->Branch("ecal"  , &ecal);
            tree->Branch("mvaBDT", &mvaBDT);
            
            tree->Branch("lxin",  &lxin);
            tree->Branch("lyin",  &lyin);
            
            tree->Branch("nhinn", &nhinn);
            tree->Branch("nhout", &nhout);
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("sign",  &sign);
            chain->SetBranchAddress("rig",   &rig);
            chain->SetBranchAddress("bta",   &bta);
            chain->SetBranchAddress("sqrm",  &sqrm);
            chain->SetBranchAddress("llr",   &llr);
            
            chain->SetBranchAddress("ecal"  , &ecal);
            chain->SetBranchAddress("mvaBDT", &mvaBDT);
            
            chain->SetBranchAddress("lxin",  &lxin);
            chain->SetBranchAddress("lyin",  &lyin);
            
            chain->SetBranchAddress("nhinn", &nhinn);
            chain->SetBranchAddress("nhout", &nhout);
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
        
        Short_t ecal;
        Float_t mvaBDT;
        
        Float_t lxin;
        Float_t lyin;
        
        Short_t nhinn;
        Short_t nhout;
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
            
            patt  = 0;
            optL1 = 0;
            optL9 = 0;
            optL2 = 0;
            nhL1 = 0;
            nhL9 = 0;
            nhL2 = 0;
            
            lxin = 0;
            lyin = 0;
            lxex = 0;
            lyex = 0;
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
            tree->Branch("optL1",  &optL1);
            tree->Branch("optL9",  &optL9);
            tree->Branch("optL2",  &optL2);
            tree->Branch("nhL1" ,  &nhL1);
            tree->Branch("nhL9" ,  &nhL9);
            tree->Branch("nhL2" ,  &nhL2);
            
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
            
            chain->SetBranchAddress("patt" , &patt);
            chain->SetBranchAddress("optL1", &optL1);
            chain->SetBranchAddress("optL9", &optL9);
            chain->SetBranchAddress("optL2", &optL2);
            chain->SetBranchAddress("nhL1" , &nhL1);
            chain->SetBranchAddress("nhL9" , &nhL9);
            chain->SetBranchAddress("nhL2" , &nhL2);
            
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
        Short_t optL1;
        Short_t optL9;
        Short_t optL2;
        Short_t nhL1;
        Short_t nhL9;
        Short_t nhL2;
        
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
            
            patt  = 0;
            optL1 = 0;
            optL9 = 0;
            optL2 = 0;
            nhL1 = 0;
            nhL9 = 0;
            nhL2 = 0;
            
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
            tree->Branch("optL1",  &optL1);
            tree->Branch("optL9",  &optL9);
            tree->Branch("optL2",  &optL2);
            tree->Branch("nhL1" ,  &nhL1);
            tree->Branch("nhL9" ,  &nhL9);
            tree->Branch("nhL2" ,  &nhL2);
            
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
            
            chain->SetBranchAddress("patt" , &patt);
            chain->SetBranchAddress("optL1", &optL1);
            chain->SetBranchAddress("optL9", &optL9);
            chain->SetBranchAddress("optL2", &optL2);
            chain->SetBranchAddress("nhL1" , &nhL1);
            chain->SetBranchAddress("nhL9" , &nhL9);
            chain->SetBranchAddress("nhL2" , &nhL2);
            
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
        Short_t optL1;
        Short_t optL9;
        Short_t optL2;
        Short_t nhL1;
        Short_t nhL9;
        Short_t nhL2;
        
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
            
            optL1 = 0;
            optL9 = 0;
            optL2 = 0;
            nhL1 = 0;
            nhL9 = 0;
            nhL2 = 0;
            
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
        
            tree->Branch("optL1", &optL1);
            tree->Branch("optL9", &optL9);
            tree->Branch("optL2", &optL2);
            tree->Branch("nhL1" , &nhL1);
            tree->Branch("nhL9" , &nhL9);
            tree->Branch("nhL2" , &nhL2);
            
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
            
            chain->SetBranchAddress("optL1", &optL1);
            chain->SetBranchAddress("optL9", &optL9);
            chain->SetBranchAddress("optL2", &optL2);
            chain->SetBranchAddress("nhL1" , &nhL1);
            chain->SetBranchAddress("nhL9" , &nhL9);
            chain->SetBranchAddress("nhL2" , &nhL2);
            
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
            
            tmva_reader->AddVariable("lrin" , &lrin);
            tmva_reader->AddVariable("lrfs" , &lrfs);
            tmva_reader->AddVariable("lrinC", &lrinC);
            tmva_reader->AddVariable("lrfsC", &lrfsC);
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
    
        Short_t optL1;
        Short_t optL9;
        Short_t optL2;
        Short_t nhL1;
        Short_t nhL9;
        Short_t nhL2;
        
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
