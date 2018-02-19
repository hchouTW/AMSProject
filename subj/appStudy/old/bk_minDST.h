#ifndef __MinDST_H__
#define __MinDST_H__

struct MinDST {
	public :
		// Event
		UInt_t  run, event;
		Float_t weight;
		// RTI
		UInt_t  uTime;
		Short_t UTCyr, UTCyd, UTChr, UTCmn;
		Float_t cfRig[4];
		Float_t lvtme;
		// MC
		Short_t mcSign;
		Float_t mcNRig;
		// TRK
		Short_t trPt;
		Float_t trAgl;
		Float_t trCooX;
		Float_t trCooY;
		Float_t trTRigIn;
		Short_t trSign[4];
		Float_t trNRig[4];
		Float_t trLchix[4];
		Float_t trLchiy[4];
		Float_t trLasym[4];
		Float_t trMxQ, trL2Q, trL1Q, trL9Q;
		// TOF
		Float_t tofQul;
		Float_t tofM;
		// TRD
		Float_t trdEst;
		Float_t trdQ;
		// ECAL
		Bool_t  hasEcal;
		Float_t ecalEst;
		// RICH
		Short_t richRad;
		Short_t richPh;
		Float_t richPrPh, richElPh;
		Bool_t  hasRich;
		Float_t richMPr, richMPi;
};

void InitMinDST(MinDST& mdst) {
	// Event
	mdst.run = 0; mdst.event = 0;
	mdst.weight = 1;
	// RTI
	mdst.uTime = 0;
	mdst.UTCyr = 0; mdst.UTCyd = 0; mdst.UTChr = 0; mdst.UTCmn = 0;
	std::fill_n(mdst.cfRig , 4, 0);
	mdst.lvtme = 0;
	// MC
	mdst.mcSign = 0;
	mdst.mcNRig = 0;
	// TRK
	mdst.trPt   = -1;
	mdst.trAgl  = 0;
	mdst.trCooX = 0;
	mdst.trCooY = 0;
	mdst.trTRigIn = 0;
	std::fill_n(mdst.trSign , 4, 0);
	std::fill_n(mdst.trNRig , 4, 0);
	std::fill_n(mdst.trLchix, 4, 0);
	std::fill_n(mdst.trLchiy, 4, 0);
	std::fill_n(mdst.trLasym, 4, 0);
	mdst.trMxQ = -1; mdst.trL2Q = -1; mdst.trL1Q = -1; mdst.trL9Q = -1;
	// TOF
	mdst.tofQul = 0;
	mdst.tofM = 0;
	// TRD
	mdst.trdEst = -1;
	mdst.trdQ   = -1;
	// ECAL
	mdst.hasEcal =  false;
	mdst.ecalEst = -2;
	// RICH
	mdst.richRad = -1;
	mdst.richPh  = 0;
	mdst.richPrPh = -1; mdst.richElPh = -1;
	mdst.hasRich = false;
	mdst.richMPr = 0; mdst.richMPi = 0;
}

void BranchMinDST(TTree * tree, MinDST& mdst) {
	InitMinDST(mdst);
	tree->Branch("run"     ,  &mdst.run     );
	tree->Branch("event"   ,  &mdst.event   );
	tree->Branch("weight"  ,  &mdst.weight  );
	tree->Branch("uTime"   ,  &mdst.uTime   );
	tree->Branch("UTCyr"   ,  &mdst.UTCyr   );
	tree->Branch("UTCyd"   ,  &mdst.UTCyd   );
	tree->Branch("UTChr"   ,  &mdst.UTChr   );
	tree->Branch("UTCmn"   ,  &mdst.UTCmn   );
	tree->Branch("cfRig"   ,   mdst.cfRig   , "cfRig[4]/F");
	tree->Branch("lvtme"   ,  &mdst.lvtme   );
	tree->Branch("mcSign"  ,  &mdst.mcSign  );
	tree->Branch("mcNRig"  ,  &mdst.mcNRig  );
	tree->Branch("trPt"    ,  &mdst.trPt    );
	tree->Branch("trAgl"   ,  &mdst.trAgl   );
	tree->Branch("trCooX"  ,  &mdst.trCooX  );
	tree->Branch("trCooY"  ,  &mdst.trCooY  );
	tree->Branch("trTRigIn",  &mdst.trTRigIn);
	tree->Branch("trSign"  ,   mdst.trSign  , "trSign[4]/S");
	tree->Branch("trNRig"  ,   mdst.trNRig  , "trNRig[4]/F");
	tree->Branch("trLchix" ,   mdst.trLchix , "trLchix[4]/F");
	tree->Branch("trLchiy" ,   mdst.trLchiy , "trLchiy[4]/F");
	tree->Branch("trLasym" ,   mdst.trLasym , "trLasym[4]/F");
	tree->Branch("trMxQ"   ,  &mdst.trMxQ   );
	tree->Branch("trL2Q"   ,  &mdst.trL2Q   );
	tree->Branch("trL1Q"   ,  &mdst.trL1Q   );
	tree->Branch("trL9Q"   ,  &mdst.trL9Q   );
	tree->Branch("tofQul"  ,  &mdst.tofQul  );
	tree->Branch("tofM"    ,  &mdst.tofM    );
	tree->Branch("trdEst"  ,  &mdst.trdEst  );
	tree->Branch("trdQ"    ,  &mdst.trdQ    );
	tree->Branch("hasEcal" ,  &mdst.hasEcal );
	tree->Branch("ecalEst" ,  &mdst.ecalEst );
	tree->Branch("richRad" ,  &mdst.richRad );
	tree->Branch("richPh"  ,  &mdst.richPh  );
	tree->Branch("richPrPh",  &mdst.richPrPh);
	tree->Branch("richElPh",  &mdst.richElPh);
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
	tree->SetBranchAddress("UTCyr"   ,  &mdst.UTCyr   );
	tree->SetBranchAddress("UTCyd"   ,  &mdst.UTCyd   );
	tree->SetBranchAddress("UTChr"   ,  &mdst.UTChr   );
	tree->SetBranchAddress("UTCmn"   ,  &mdst.UTCmn   );
	tree->SetBranchAddress("cfRig"   ,   mdst.cfRig   );
	tree->SetBranchAddress("lvtme"   ,  &mdst.lvtme   );
	tree->SetBranchAddress("mcSign"  ,  &mdst.mcSign  );
	tree->SetBranchAddress("mcNRig"  ,  &mdst.mcNRig  );
	tree->SetBranchAddress("trPt"    ,  &mdst.trPt    );
	tree->SetBranchAddress("trAgl"   ,  &mdst.trAgl   );
	tree->SetBranchAddress("trCooX"  ,  &mdst.trCooX  );
	tree->SetBranchAddress("trCooY"  ,  &mdst.trCooY  );
	tree->SetBranchAddress("trTRigIn",  &mdst.trTRigIn);
	tree->SetBranchAddress("trSign"  ,   mdst.trSign  );
	tree->SetBranchAddress("trNRig"  ,   mdst.trNRig  );
	tree->SetBranchAddress("trLchix" ,   mdst.trLchix );
	tree->SetBranchAddress("trLchiy" ,   mdst.trLchiy );
	tree->SetBranchAddress("trLasym" ,   mdst.trLasym );
	tree->SetBranchAddress("trMxQ"   ,  &mdst.trMxQ   );
	tree->SetBranchAddress("trL2Q"   ,  &mdst.trL2Q   );
	tree->SetBranchAddress("trL1Q"   ,  &mdst.trL1Q   );
	tree->SetBranchAddress("trL9Q"   ,  &mdst.trL9Q   );
	tree->SetBranchAddress("tofQul"  ,  &mdst.tofQul  );
	tree->SetBranchAddress("tofM"    ,  &mdst.tofM    );
	tree->SetBranchAddress("trdEst"  ,  &mdst.trdEst  );
	tree->SetBranchAddress("trdQ"    ,  &mdst.trdQ    );
	tree->SetBranchAddress("hasEcal" ,  &mdst.hasEcal );
	tree->SetBranchAddress("ecalEst" ,  &mdst.ecalEst );
	tree->SetBranchAddress("richRad" ,  &mdst.richRad );
	tree->SetBranchAddress("richPh"  ,  &mdst.richPh  );
	tree->SetBranchAddress("richPrPh",  &mdst.richPrPh);
	tree->SetBranchAddress("richElPh",  &mdst.richElPh);
	tree->SetBranchAddress("hasRich" ,  &mdst.hasRich );
	tree->SetBranchAddress("richMPr" ,  &mdst.richMPr );
	tree->SetBranchAddress("richMPi" ,  &mdst.richMPi );
}

#endif // __MinDST_H__
