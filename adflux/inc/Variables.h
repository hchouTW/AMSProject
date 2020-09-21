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


class VarsTF : public Vars {
    public :
        VarsTF() : Vars() { init(); }
        ~VarsTF() {}
        
        void init() {
            init_base();

            chrg    = 0;
            sign    = 0;
            rig     = 0;
            bta     = 0;
            sqrm    = 0;
            llr     = 0;
           
            L2  = 0;
            nhx = 0;
            nhy = 0;
            
            ftL2  = 0;
            ftSL1 = 0;
            ftSL2 = 0;
            ftSL3 = 0;

            rich    = 0;
            veto    = 0;
            ncls    = 0;
            nhinn   = 0;
            nhout   = 0;
            nhinn2   = 0;
            nhout2   = 0;

            std::fill_n(nvtxx, 2, 0);
            std::fill_n(nvtxy, 2, 0);

            ext = 0;
            extlx = 0;
            extly = 0;

            geom_lx = 0; 
            geom_ly = 0; 
            
            std::fill_n(tf_bta, 2, 0); 
            
            std::fill_n(vel_lb , 2, 0); 
            std::fill_n(mutr_lx, 2, 0); 
            std::fill_n(mutr_ly, 2, 0); 
            std::fill_n(mutr_lb, 2, 0); 
            std::fill_n(a_phys_lx, 2, 0); 
            std::fill_n(a_phys_ly, 2, 0); 
            std::fill_n(a_phys_lb, 2, 0); 
            std::fill_n(b_phys_lx, 2, 0); 
            std::fill_n(b_phys_ly, 2, 0); 
            std::fill_n(b_phys_lb, 2, 0); 
            
            std::fill_n(a_top_rig, 2, 0);
            std::fill_n(a_top_bta, 2, 0);
            std::fill_n(b_top_rig, 2, 0);
            std::fill_n(b_top_bta, 2, 0);
            
            std::fill_n(a_cen_rig, 2, 0);
            std::fill_n(a_cen_bta, 2, 0);
            std::fill_n(b_cen_rig, 2, 0);
            std::fill_n(b_cen_bta, 2, 0);
        }
        
        void branch() {
            if (!tree) return;
            tree->Branch("chrg",    &chrg);
            tree->Branch("sign",    &sign);
            tree->Branch("rig",     &rig);
            tree->Branch("bta",     &bta);
            tree->Branch("sqrm",    &sqrm);
            tree->Branch("llr",     &llr);
            
            tree->Branch("L2" , &L2);
            tree->Branch("nhx", &nhx);
            tree->Branch("nhy", &nhy);
            
            tree->Branch("ftL2" , &ftL2 );
            tree->Branch("ftSL1", &ftSL1);
            tree->Branch("ftSL2", &ftSL2);
            tree->Branch("ftSL3", &ftSL3);
            
            tree->Branch("rich",    &rich);
            tree->Branch("veto",    &veto);
            tree->Branch("ncls",    &ncls);
            tree->Branch("nhinn",   &nhinn);
            tree->Branch("nhout",   &nhout);
            tree->Branch("nhinn2",   &nhinn2);
            tree->Branch("nhout2",   &nhout2);
            
            tree->Branch("nvtxx",   nvtxx,  "nvtxx[2]/S");
            tree->Branch("nvtxy",   nvtxy,  "nvtxy[2]/S");
            
            tree->Branch("ext", &ext);
            tree->Branch("extlx", &extlx);
            tree->Branch("extly", &extly);
            
            tree->Branch("geom_lx", &geom_lx);
            tree->Branch("geom_ly", &geom_ly);
            
            tree->Branch("tf_bta", tf_bta, "tf_bta[2]/F");

            tree->Branch("vel_lb" , vel_lb , "vel_lb[2]/F" );
            tree->Branch("mutr_lx", mutr_lx, "mutr_lx[2]/F");
            tree->Branch("mutr_ly", mutr_ly, "mutr_ly[2]/F");
            tree->Branch("mutr_lb", mutr_lb, "mutr_lb[2]/F");
            tree->Branch("a_phys_lx", a_phys_lx, "a_phys_lx[2]/F");
            tree->Branch("a_phys_ly", a_phys_ly, "a_phys_ly[2]/F");
            tree->Branch("a_phys_lb", a_phys_lb, "a_phys_lb[2]/F");
            tree->Branch("b_phys_lx", b_phys_lx, "b_phys_lx[2]/F");
            tree->Branch("b_phys_ly", b_phys_ly, "b_phys_ly[2]/F");
            tree->Branch("b_phys_lb", b_phys_lb, "b_phys_lb[2]/F");
            
            tree->Branch("a_top_rig", a_top_rig, "a_top_rig[2]/F");
            tree->Branch("a_top_bta", a_top_bta, "a_top_bta[2]/F");
            tree->Branch("b_top_rig", b_top_rig, "b_top_rig[2]/F");
            tree->Branch("b_top_bta", b_top_bta, "b_top_bta[2]/F");
            
            tree->Branch("a_cen_rig", a_cen_rig, "a_cen_rig[2]/F");
            tree->Branch("a_cen_bta", a_cen_bta, "a_cen_bta[2]/F");
            tree->Branch("b_cen_rig", b_cen_rig, "b_cen_rig[2]/F");
            tree->Branch("b_cen_bta", b_cen_bta, "b_cen_bta[2]/F");
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("chrg",    &chrg);
            chain->SetBranchAddress("sign",    &sign);
            chain->SetBranchAddress("rig",     &rig);
            chain->SetBranchAddress("bta",     &bta);
            chain->SetBranchAddress("sqrm",    &sqrm);
            chain->SetBranchAddress("llr",     &llr);
            
            chain->SetBranchAddress("L2" , &L2);
            chain->SetBranchAddress("nhx", &nhx);
            chain->SetBranchAddress("nhy", &nhy);
            
            chain->SetBranchAddress("ftL2" , &ftL2 );
            chain->SetBranchAddress("ftSL1", &ftSL1);
            chain->SetBranchAddress("ftSL2", &ftSL2);
            chain->SetBranchAddress("ftSL3", &ftSL3);
            
            chain->SetBranchAddress("rich",    &rich);
            chain->SetBranchAddress("veto",    &veto);
            chain->SetBranchAddress("ncls",    &ncls);
            chain->SetBranchAddress("nhinn",   &nhinn);
            chain->SetBranchAddress("nhout",   &nhout);
            chain->SetBranchAddress("nhinn2",   &nhinn2);
            chain->SetBranchAddress("nhout2",   &nhout2);
            
            chain->SetBranchAddress("nvtxx",   nvtxx);
            chain->SetBranchAddress("nvtxy",   nvtxy);
            
            chain->SetBranchAddress("ext", &ext);
            chain->SetBranchAddress("extlx", &extlx);
            chain->SetBranchAddress("extly", &extly);
            
            chain->SetBranchAddress("geom_lx", &geom_lx);
            chain->SetBranchAddress("geom_ly", &geom_ly);
            
            chain->SetBranchAddress("tf_bta", tf_bta);
            
            chain->SetBranchAddress("vel_lb" , vel_lb );
            chain->SetBranchAddress("mutr_lx", mutr_lx);
            chain->SetBranchAddress("mutr_ly", mutr_ly);
            chain->SetBranchAddress("mutr_lb", mutr_lb);
            chain->SetBranchAddress("a_phys_lx", a_phys_lx);
            chain->SetBranchAddress("a_phys_ly", a_phys_ly);
            chain->SetBranchAddress("a_phys_lb", a_phys_lb);
            chain->SetBranchAddress("b_phys_lx", b_phys_lx);
            chain->SetBranchAddress("b_phys_ly", b_phys_ly);
            chain->SetBranchAddress("b_phys_lb", b_phys_lb);
            
            chain->SetBranchAddress("a_top_rig", a_top_rig);
            chain->SetBranchAddress("a_top_bta", a_top_bta);
            chain->SetBranchAddress("b_top_rig", b_top_rig);
            chain->SetBranchAddress("b_top_bta", b_top_bta);
            
            chain->SetBranchAddress("a_cen_rig", a_cen_rig);
            chain->SetBranchAddress("a_cen_bta", a_cen_bta);
            chain->SetBranchAddress("b_cen_rig", b_cen_rig);
            chain->SetBranchAddress("b_cen_bta", b_cen_bta);
        }

#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t chrg;
        Short_t sign;
        Float_t rig;
        Float_t bta;
        Float_t sqrm;
        Float_t llr;
        
        Short_t L2;
        Short_t nhx;
        Short_t nhy;
        
        Float_t ftL2;
        Float_t ftSL1;
        Float_t ftSL2;
        Float_t ftSL3;
        
        Short_t rich;
        Short_t veto;
        Short_t ncls;
        Short_t nhinn;
        Short_t nhout;
        Short_t nhinn2;
        Short_t nhout2;
        
        Short_t nvtxx[2];
        Short_t nvtxy[2];
        
        Short_t ext;
        Float_t extlx;
        Float_t extly;
        
        Float_t geom_lx;
        Float_t geom_ly;
        
        Float_t tf_bta[2];
        
        Float_t vel_lb [2];
        Float_t mutr_lx[2];
        Float_t mutr_ly[2];
        Float_t mutr_lb[2];
        Float_t a_phys_lx[2];
        Float_t a_phys_ly[2];
        Float_t a_phys_lb[2];
        Float_t b_phys_lx[2];
        Float_t b_phys_ly[2];
        Float_t b_phys_lb[2];
        
        Float_t a_top_rig[2];
        Float_t a_top_bta[2];
        Float_t b_top_rig[2];
        Float_t b_top_bta[2];
        
        Float_t a_cen_rig[2];
        Float_t a_cen_bta[2];
        Float_t b_cen_rig[2];
        Float_t b_cen_bta[2];
};


class VarsRH : public Vars {
    public :
        VarsRH() : Vars() { init(); }
        ~VarsRH() {}
        
        void init() {
            init_base();
           
            chrg  = 0;
            sign  = 0;
            rig   = 0;
            bta   = 0;
            sqrm  = 0;
            llr   = 0;
            
            L2  = 0;
            nhx = 0;
            nhy = 0;
          
            ftL2  = 0;
            ftSL1 = 0;
            ftSL2 = 0;
            ftSL3 = 0;

            nhinn   = 0;
            nhout   = 0;
            nhinn2  = 0;
            nhout2  = 0;
            
            nclsZ1  = 0;
            selfqq  = 0;
            
            std::fill_n(rhbta, 2, 0); 
            
            std::fill_n(nvtxx, 2, 0);
            std::fill_n(nvtxy, 2, 0);
            
            ext = 0;
            extlx = 0;
            extly = 0;
            
            geom_lx = 0; 
            geom_ly = 0; 
            
            std::fill_n(tf_bta, 2, 0); 
            std::fill_n(rh_bta, 2, 0); 
            
            std::fill_n(vel_lb , 4, 0); 
            std::fill_n(mutr_lx, 4, 0); 
            std::fill_n(mutr_ly, 4, 0); 
            std::fill_n(mutr_lb, 4, 0); 
            std::fill_n(a_phys_lx, 4, 0); 
            std::fill_n(a_phys_ly, 4, 0); 
            std::fill_n(a_phys_lb, 4, 0); 
            std::fill_n(b_phys_lx, 4, 0); 
            std::fill_n(b_phys_ly, 4, 0); 
            std::fill_n(b_phys_lb, 4, 0); 
            
            std::fill_n(a_top_rig, 4, 0);
            std::fill_n(a_top_bta, 4, 0);
            std::fill_n(b_top_rig, 4, 0);
            std::fill_n(b_top_bta, 4, 0);
            
            std::fill_n(a_cen_rig, 4, 0);
            std::fill_n(a_cen_bta, 4, 0);
            std::fill_n(b_cen_rig, 4, 0);
            std::fill_n(b_cen_bta, 4, 0);
        }

        void branch() {
            if (!tree) return;
            tree->Branch("chrg",  &chrg);
            tree->Branch("sign",  &sign);
            tree->Branch("rig",   &rig);
            tree->Branch("bta",   &bta);
            tree->Branch("sqrm",  &sqrm);
            tree->Branch("llr",   &llr);
            
            tree->Branch("L2" , &L2);
            tree->Branch("nhx", &nhx);
            tree->Branch("nhy", &nhy);
            
            tree->Branch("ftL2" , &ftL2 );
            tree->Branch("ftSL1", &ftSL1);
            tree->Branch("ftSL2", &ftSL2);
            tree->Branch("ftSL3", &ftSL3);
            
            tree->Branch("nhinn", &nhinn);
            tree->Branch("nhout", &nhout);
            tree->Branch("nhinn2", &nhinn2);
            tree->Branch("nhout2", &nhout2);
            
            tree->Branch("nclsZ1", &nclsZ1);
            tree->Branch("selfqq", &selfqq);
            
            tree->Branch("rhbta", rhbta, "rhbta[2]/F");
            
            tree->Branch("nvtxx",   nvtxx,  "nvtxx[2]/S");
            tree->Branch("nvtxy",   nvtxy,  "nvtxy[2]/S");
            
            tree->Branch("ext", &ext);
            tree->Branch("extlx", &extlx);
            tree->Branch("extly", &extly);
            
            tree->Branch("geom_lx", &geom_lx);
            tree->Branch("geom_ly", &geom_ly);
            
            tree->Branch("tf_bta", tf_bta, "tf_bta[2]/F");
            tree->Branch("rh_bta", rh_bta, "rh_bta[2]/F");
            
            tree->Branch("vel_lb" , vel_lb , "vel_lb[4]/F" );
            tree->Branch("mutr_lx", mutr_lx, "mutr_lx[4]/F");
            tree->Branch("mutr_ly", mutr_ly, "mutr_ly[4]/F");
            tree->Branch("mutr_lb", mutr_lb, "mutr_lb[4]/F");
            tree->Branch("a_phys_lx", a_phys_lx, "a_phys_lx[4]/F");
            tree->Branch("a_phys_ly", a_phys_ly, "a_phys_ly[4]/F");
            tree->Branch("a_phys_lb", a_phys_lb, "a_phys_lb[4]/F");
            tree->Branch("b_phys_lx", b_phys_lx, "b_phys_lx[4]/F");
            tree->Branch("b_phys_ly", b_phys_ly, "b_phys_ly[4]/F");
            tree->Branch("b_phys_lb", b_phys_lb, "b_phys_lb[4]/F");
            
            tree->Branch("a_top_rig", a_top_rig, "a_top_rig[4]/F");
            tree->Branch("a_top_bta", a_top_bta, "a_top_bta[4]/F");
            tree->Branch("b_top_rig", b_top_rig, "b_top_rig[4]/F");
            tree->Branch("b_top_bta", b_top_bta, "b_top_bta[4]/F");
            
            tree->Branch("a_cen_rig", a_cen_rig, "a_cen_rig[4]/F");
            tree->Branch("a_cen_bta", a_cen_bta, "a_cen_bta[4]/F");
            tree->Branch("b_cen_rig", b_cen_rig, "b_cen_rig[4]/F");
            tree->Branch("b_cen_bta", b_cen_bta, "b_cen_bta[4]/F");
        }
        
        void set_branch_address() {
            if (!chain) return;
            chain->SetBranchAddress("chrg",  &chrg);
            chain->SetBranchAddress("sign",  &sign);
            chain->SetBranchAddress("rig",   &rig);
            chain->SetBranchAddress("bta",   &bta);
            chain->SetBranchAddress("sqrm",  &sqrm);
            chain->SetBranchAddress("llr",   &llr);
            
            chain->SetBranchAddress("L2" , &L2);
            chain->SetBranchAddress("nhx", &nhx);
            chain->SetBranchAddress("nhy", &nhy);
            
            chain->SetBranchAddress("ftL2" , &ftL2 );
            chain->SetBranchAddress("ftSL1", &ftSL1);
            chain->SetBranchAddress("ftSL2", &ftSL2);
            chain->SetBranchAddress("ftSL3", &ftSL3);
            
            chain->SetBranchAddress("nhinn", &nhinn);
            chain->SetBranchAddress("nhout", &nhout);
            chain->SetBranchAddress("nhinn2", &nhinn2);
            chain->SetBranchAddress("nhout2", &nhout2);
            
            chain->SetBranchAddress("nclsZ1", &nclsZ1);
            chain->SetBranchAddress("selfqq", &selfqq);
            
            chain->SetBranchAddress("rhbta", rhbta);
            
            chain->SetBranchAddress("nvtxx",   nvtxx);
            chain->SetBranchAddress("nvtxy",   nvtxy);
            
            chain->SetBranchAddress("ext", &ext);
            chain->SetBranchAddress("extlx", &extlx);
            chain->SetBranchAddress("extly", &extly);
            
            chain->SetBranchAddress("geom_lx", &geom_lx);
            chain->SetBranchAddress("geom_ly", &geom_ly);
            
            chain->SetBranchAddress("tf_bta", tf_bta);
            chain->SetBranchAddress("rh_bta", rh_bta);
            
            chain->SetBranchAddress("vel_lb" , vel_lb );
            chain->SetBranchAddress("mutr_lx", mutr_lx);
            chain->SetBranchAddress("mutr_ly", mutr_ly);
            chain->SetBranchAddress("mutr_lb", mutr_lb);
            chain->SetBranchAddress("a_phys_lx", a_phys_lx);
            chain->SetBranchAddress("a_phys_ly", a_phys_ly);
            chain->SetBranchAddress("a_phys_lb", a_phys_lb);
            chain->SetBranchAddress("b_phys_lx", b_phys_lx);
            chain->SetBranchAddress("b_phys_ly", b_phys_ly);
            chain->SetBranchAddress("b_phys_lb", b_phys_lb);
            
            chain->SetBranchAddress("a_top_rig", a_top_rig);
            chain->SetBranchAddress("a_top_bta", a_top_bta);
            chain->SetBranchAddress("b_top_rig", b_top_rig);
            chain->SetBranchAddress("b_top_bta", b_top_bta);
            
            chain->SetBranchAddress("a_cen_rig", a_cen_rig);
            chain->SetBranchAddress("a_cen_bta", a_cen_bta);
            chain->SetBranchAddress("b_cen_rig", b_cen_rig);
            chain->SetBranchAddress("b_cen_bta", b_cen_bta);
        }
        
#ifdef __USE_TMVA_READER__
        void set_tmva_vars() {
            if (!tmva_reader) return;
        }
#endif // __USE_TMVA_READER__

    public :
        Short_t chrg;
        Short_t sign;
        Float_t rig;
        Float_t bta;
        Float_t sqrm;
        Float_t llr;
        
        Short_t L2;
        Short_t nhx;
        Short_t nhy;
        
        Float_t ftL2;
        Float_t ftSL1;
        Float_t ftSL2;
        Float_t ftSL3;
        
        Short_t nhinn;
        Short_t nhout;
        Short_t nhinn2;
        Short_t nhout2;
        
        Short_t nclsZ1;
        Float_t selfqq;
       
        Float_t rhbta[2];

        Short_t nvtxx[2];
        Short_t nvtxy[2];
       
        Short_t ext;
        Float_t extlx;
        Float_t extly;

        Float_t geom_lx;
        Float_t geom_ly;
        
        Float_t tf_bta[2];
        Float_t rh_bta[2];
       
        Float_t vel_lb[4];
        Float_t mutr_lx[4];
        Float_t mutr_ly[4];
        Float_t mutr_lb[4];
        Float_t a_phys_lx[4];
        Float_t a_phys_ly[4];
        Float_t a_phys_lb[4];
        Float_t b_phys_lx[4];
        Float_t b_phys_ly[4];
        Float_t b_phys_lb[4];

        Float_t a_top_rig[4];
        Float_t a_top_bta[4];
        Float_t b_top_rig[4];
        Float_t b_top_bta[4];
        
        Float_t a_cen_rig[4];
        Float_t a_cen_bta[4];
        Float_t b_cen_rig[4];
        Float_t b_cen_bta[4];
};


#endif // __Variables_H__
