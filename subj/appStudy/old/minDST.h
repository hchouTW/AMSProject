#ifndef __MinDST_H__
#define __MinDST_H__

struct MinDST {
    public :
        // Event
        UInt_t  run, event;
        Float_t weight;
        // RTI
        UInt_t  uTime;
        Float_t cfRig;
        // MC
        Short_t mcSign;
        Float_t mcNRig;
        Bool_t  mcInt;
        Float_t mcVtx[3];
        Float_t mcFrac;
        // TRK
        Short_t trPt;
        Short_t trNHL1;
        Short_t trNHL9;
        Short_t sign[4];
        Float_t nrig[4];
        Float_t lchix[4];
        Float_t lchiy[4];
        Float_t lasym1I;
        Float_t lasym9I;
        Float_t lasym91;
        // TOF
        Float_t tofM;
        Short_t tofCN[2];
        Short_t tofCL[2][2];
        Float_t tofCQ[2][2];
        Float_t tofCT[2][2];
        // TRD
        Float_t trdEst;
        Bool_t  trdVtx;
        Short_t trdSeg[2];
        // ECAL
        Bool_t  hasEcal;
        Float_t ecalEst;
        // RICH
        Short_t richRad;
        Float_t richPhEl, richPhPr;
        Short_t richCh[2];
        Short_t richRhEl[3];
        Short_t richRhPr[3];
        Bool_t  hasRich;
        Float_t richMPr, richMPi;
};

void InitMinDST(MinDST& mdst) {
    // Event
    mdst.run = 0; mdst.event = 0;
    mdst.weight = 1;
    // RTI
    mdst.uTime = 0;
    mdst.cfRig = 0;
    // MC
    mdst.mcSign = 0;
    mdst.mcNRig = 0;
    mdst.mcInt = false;
    std::fill_n(mdst.mcVtx, 3, 0);
    mdst.mcFrac = -1;
    // TRK
    mdst.trPt = -1;
    mdst.trNHL1 = 0;
    mdst.trNHL9 = 0;
    std::fill_n(mdst.sign , 4, 0);
    std::fill_n(mdst.nrig , 4, 0);
    std::fill_n(mdst.lchix, 4, 0);
    std::fill_n(mdst.lchiy, 4, 0);
    mdst.lasym1I = 0;
    mdst.lasym9I = 0;
    mdst.lasym91 = 0;
    // TOF
    mdst.tofM = 0;
    std::fill_n(mdst.tofCN, 2, 0);
    std::fill_n(mdst.tofCL[0], 2*2, -1);
    std::fill_n(mdst.tofCQ[0], 2*2, -1);
    std::fill_n(mdst.tofCT[0], 2*2, -10);
    // TRD
    mdst.trdEst = -1;
    mdst.trdVtx = false;
    std::fill_n(mdst.trdSeg, 2, 0);
    // ECAL
    mdst.hasEcal =  false;
    mdst.ecalEst = -2;
    // RICH
    mdst.richRad = -1;
    mdst.richPhEl = 0; mdst.richPhPr = 0;
    std::fill_n(mdst.richCh, 2, 0);
    std::fill_n(mdst.richRhEl, 3, 0);
    std::fill_n(mdst.richRhPr, 3, 0);
    mdst.hasRich = false;
    mdst.richMPr = 0; mdst.richMPi = 0;
}

void BranchMinDST(TTree * tree, MinDST& mdst) {
    InitMinDST(mdst);
    tree->Branch("run"     ,  &mdst.run     );
    tree->Branch("event"   ,  &mdst.event   );
    tree->Branch("weight"  ,  &mdst.weight  );
    tree->Branch("uTime"   ,  &mdst.uTime   );
    tree->Branch("cfRig"   ,  &mdst.cfRig   );
    tree->Branch("mcSign"  ,  &mdst.mcSign  );
    tree->Branch("mcNRig"  ,  &mdst.mcNRig  );
    tree->Branch("mcInt"   ,  &mdst.mcInt   );
    tree->Branch("mcVtx"   ,   mdst.mcVtx, "mcVtx[3]/F");
    tree->Branch("mcFrac"  ,  &mdst.mcFrac  );
    tree->Branch("trPt"    ,  &mdst.trPt    );
    tree->Branch("trNHL1"  ,  &mdst.trNHL1  );
    tree->Branch("trNHL9"  ,  &mdst.trNHL9  );
    tree->Branch("sign"    ,   mdst.sign , "sign[4]/S");
    tree->Branch("nrig"    ,   mdst.nrig , "nrig[4]/F");
    tree->Branch("lchix"   ,   mdst.lchix, "lchix[4]/F");
    tree->Branch("lchiy"   ,   mdst.lchiy, "lchiy[4]/F");
    tree->Branch("lasym1I" ,  &mdst.lasym1I );
    tree->Branch("lasym9I" ,  &mdst.lasym9I );
    tree->Branch("lasym91" ,  &mdst.lasym91 );
    tree->Branch("tofM"    ,  &mdst.tofM    );
    tree->Branch("tofCN"   ,   mdst.tofCN, "tofCN[2]/S");
    tree->Branch("tofCL"   ,   mdst.tofCL, "tofCL[2][2]/S");
    tree->Branch("tofCQ"   ,   mdst.tofCQ, "tofCQ[2][2]/F");
    tree->Branch("tofCT"   ,   mdst.tofCT, "tofCT[2][2]/F");
    tree->Branch("trdEst"  ,  &mdst.trdEst  );
    tree->Branch("trdVtx"  ,  &mdst.trdVtx  );
    tree->Branch("trdSeg"  ,   mdst.trdSeg, "trdSeg[2]/S");
    tree->Branch("hasEcal" ,  &mdst.hasEcal );
    tree->Branch("ecalEst" ,  &mdst.ecalEst );
    tree->Branch("richRad" ,  &mdst.richRad );
    tree->Branch("richPhEl",  &mdst.richPhEl);
    tree->Branch("richPhPr",  &mdst.richPhPr);
    tree->Branch("richCh"  ,   mdst.richCh, "richCh[2]/S");
    tree->Branch("richRhEl",   mdst.richRhEl, "richRhEl[3]/S");
    tree->Branch("richRhPr",   mdst.richRhPr, "richRhPr[3]/S");
    tree->Branch("hasRich" ,  &mdst.hasRich );
    tree->Branch("richMPr" ,  &mdst.richMPr );
    tree->Branch("richMPi" ,  &mdst.richMPi );
}

void SetBranchAddressMinDST(TTree * tree, MinDST& mdst) {
    InitMinDST(mdst);
    tree->SetBranchAddress("run"     ,  &mdst.run     );
    tree->SetBranchAddress("event"   ,  &mdst.event   );
    tree->SetBranchAddress("weight"  ,  &mdst.weight  );
    tree->SetBranchAddress("uTime"   ,  &mdst.uTime   );
    tree->SetBranchAddress("cfRig"   ,  &mdst.cfRig   );
    tree->SetBranchAddress("mcSign"  ,  &mdst.mcSign  );
    tree->SetBranchAddress("mcNRig"  ,  &mdst.mcNRig  );
    tree->SetBranchAddress("mcInt"   ,  &mdst.mcInt   );
    tree->SetBranchAddress("mcVtx"   ,   mdst.mcVtx   );
    tree->SetBranchAddress("mcFrac"  ,  &mdst.mcFrac  );
    tree->SetBranchAddress("trPt"    ,  &mdst.trPt    );
    tree->SetBranchAddress("trNHL1"  ,  &mdst.trNHL1  );
    tree->SetBranchAddress("trNHL9"  ,  &mdst.trNHL9  );
    tree->SetBranchAddress("sign"    ,   mdst.sign    );
    tree->SetBranchAddress("nrig"    ,   mdst.nrig    );
    tree->SetBranchAddress("lchix"   ,   mdst.lchix   );
    tree->SetBranchAddress("lchiy"   ,   mdst.lchiy   );
    tree->SetBranchAddress("lasym1I" ,  &mdst.lasym1I );
    tree->SetBranchAddress("lasym9I" ,  &mdst.lasym9I );
    tree->SetBranchAddress("lasym91" ,  &mdst.lasym91 );
    tree->SetBranchAddress("tofM"    ,  &mdst.tofM    );
    tree->SetBranchAddress("tofCN"   ,   mdst.tofCN   );
    tree->SetBranchAddress("tofCL"   ,   mdst.tofCL   );
    tree->SetBranchAddress("tofCQ"   ,   mdst.tofCQ   );
    tree->SetBranchAddress("tofCT"   ,   mdst.tofCT   );
    tree->SetBranchAddress("trdEst"  ,  &mdst.trdEst  );
    tree->SetBranchAddress("trdVtx"  ,  &mdst.trdVtx  );
    tree->SetBranchAddress("trdSeg"  ,   mdst.trdSeg  );
    tree->SetBranchAddress("hasEcal" ,  &mdst.hasEcal );
    tree->SetBranchAddress("ecalEst" ,  &mdst.ecalEst );
    tree->SetBranchAddress("richRad" ,  &mdst.richRad );
    tree->SetBranchAddress("richPhEl",  &mdst.richPhEl);
    tree->SetBranchAddress("richPhPr",  &mdst.richPhPr);
    tree->SetBranchAddress("richCh"  ,   mdst.richCh  );
    tree->SetBranchAddress("richRhEl",   mdst.richRhEl);
    tree->SetBranchAddress("richRhPr",   mdst.richRhPr);
    tree->SetBranchAddress("hasRich" ,  &mdst.hasRich );
    tree->SetBranchAddress("richMPr" ,  &mdst.richMPr );
    tree->SetBranchAddress("richMPi" ,  &mdst.richMPi );
}

#endif // __MinDST_H__
