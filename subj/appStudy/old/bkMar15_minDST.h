#ifndef __MinDST_H__
#define __MinDST_H__

struct MinDST {
	public :
		// Event
		UInt_t  run, event;
		Float_t weight;
		// RTI
		UInt_t  uTime;
		Short_t UTCyr, UTCyd, UTChr;
		Float_t cfRig;
		Float_t cfScl;
		Float_t lvtme;
		// MC
		Short_t mcSign;
		Float_t mcNRig;
		// TRK
		Short_t trPt;
		Short_t trSign;
		Float_t trNRig;
		
		Float_t trLchix;
		Float_t trLchiy[4];
		Float_t trLasym[4];
		Float_t trCCest;
		
		//Float_t trLchix, trLchiy, trLasym, trCCest;
		Float_t trMxQ, trL2Q, trL1Q, trL9Q;
		// TOF
		Float_t tofQul;
		//Float_t tofb;
		Float_t tofM;
		// TRD
		Float_t trdEst;
		// ECAL
		Bool_t  hasEcal;
		Float_t ecalEst;
		// RICH
		Short_t richRad;
		Short_t richPh;
		Float_t richPrPh, richElPh;
		Bool_t  hasRich;
		//Float_t richbPr, richbPi, richbEl;
		Float_t richMPr, richMPi, richMEl;
};

void InitMinDST(MinDST& mdst) {
	// Event
	mdst.run = 0; mdst.event = 0;
	mdst.weight = 1;
	// RTI
	mdst.uTime = 0;
	mdst.UTCyr = 0; mdst.UTCyd = 0; mdst.UTChr = 0;
	mdst.cfRig = 0;
	mdst.cfScl = 0;
	mdst.lvtme = 0;
	// MC
	mdst.mcSign = 0;
	mdst.mcNRig = 0;
	// TRK
	mdst.trPt = -1;
	mdst.trSign = 0;
	mdst.trNRig = 0;
	mdst.trLchix = 0;
	std::fill_n(mdst.trLchiy , 4, 0);
	std::fill_n(mdst.trLasym , 4, 0);
	mdst.trCCest = 0;
	//mdst.trLchix = 0; mdst.trLchiy = 0, mdst.trLasym = 0, mdst.trCCest = 0;
	mdst.trMxQ = -1; mdst.trL2Q = -1; mdst.trL1Q = -1; mdst.trL9Q = -1;
	// TOF
	mdst.tofQul = 0;
	//mdst.tofb = 0;
	mdst.tofM = 0;
	// TRD
	mdst.trdEst = -1;
	// ECAL
	mdst.hasEcal =  false;
	mdst.ecalEst = -2;
	// RICH
	mdst.richRad = -1;
	mdst.richPh = 0;
	mdst.richPrPh = 0; mdst.richElPh = 0;
	mdst.hasRich = false;
	//mdst.richbPr = 0; mdst.richbPi = 0; mdst.richbEl = 0;
	mdst.richMPr = 0; mdst.richMPi = 0; mdst.richMEl = 0;
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
	tree->Branch("cfRig"   ,  &mdst.cfRig   );
	tree->Branch("cfScl"   ,  &mdst.cfScl   );
	tree->Branch("lvtme"   ,  &mdst.lvtme   );
	tree->Branch("mcSign"  ,  &mdst.mcSign  );
	tree->Branch("mcNRig"  ,  &mdst.mcNRig  );
	tree->Branch("trPt"    ,  &mdst.trPt    );
	tree->Branch("trSign"  ,  &mdst.trSign  );
	tree->Branch("trNRig"  ,  &mdst.trNRig  );
	
	tree->Branch("trLchix" ,  &mdst.trLchix );
	tree->Branch("trLchiy" ,   mdst.trLchiy , "trLchiy[4]/F");
	tree->Branch("trLasym" ,   mdst.trLasym , "trLasym[4]/F");
	tree->Branch("trCCest" ,  &mdst.trCCest );
	
	//tree->Branch("trLchix" ,  &mdst.trLchix );
	//tree->Branch("trLchiy" ,  &mdst.trLchiy );
	//tree->Branch("trLasym" ,  &mdst.trLasym );
	//tree->Branch("trCCest" ,  &mdst.trCCest );
	tree->Branch("trMxQ"   ,  &mdst.trMxQ   );
	tree->Branch("trL2Q"   ,  &mdst.trL2Q   );
	tree->Branch("trL1Q"   ,  &mdst.trL1Q   );
	tree->Branch("trL9Q"   ,  &mdst.trL9Q   );
	tree->Branch("tofQul"  ,  &mdst.tofQul  );
	//tree->Branch("tofb"    ,  &mdst.tofb    );
	tree->Branch("tofM"    ,  &mdst.tofM    );
	tree->Branch("trdEst"  ,  &mdst.trdEst  );
	tree->Branch("hasEcal" ,  &mdst.hasEcal );
	tree->Branch("ecalEst" ,  &mdst.ecalEst );
	tree->Branch("richRad" ,  &mdst.richRad );
	tree->Branch("richPh"  ,  &mdst.richPh  );
	tree->Branch("richPrPh",  &mdst.richPrPh);
	tree->Branch("richElPh",  &mdst.richElPh);
	tree->Branch("hasRich" ,  &mdst.hasRich );
	//tree->Branch("richbPr" ,  &mdst.richbPr );
	//tree->Branch("richbPi" ,  &mdst.richbPi );
	//tree->Branch("richbEl" ,  &mdst.richbEl );
	tree->Branch("richMPr" ,  &mdst.richMPr );
	tree->Branch("richMPi" ,  &mdst.richMPi );
	tree->Branch("richMEl" ,  &mdst.richMEl );
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
	tree->SetBranchAddress("cfRig"   ,  &mdst.cfRig   );
	tree->SetBranchAddress("cfScl"   ,  &mdst.cfScl   );
	tree->SetBranchAddress("lvtme"   ,  &mdst.lvtme   );
	tree->SetBranchAddress("mcSign"  ,  &mdst.mcSign  );
	tree->SetBranchAddress("mcNRig"  ,  &mdst.mcNRig  );
	tree->SetBranchAddress("trPt"    ,  &mdst.trPt    );
	tree->SetBranchAddress("trSign"  ,  &mdst.trSign  );
	tree->SetBranchAddress("trNRig"  ,  &mdst.trNRig  );
	
	tree->SetBranchAddress("trLchix" ,  &mdst.trLchix );
	tree->SetBranchAddress("trLchiy" ,   mdst.trLchiy );
	tree->SetBranchAddress("trLasym" ,   mdst.trLasym );
	tree->SetBranchAddress("trCCest" ,  &mdst.trCCest );
	
	//tree->SetBranchAddress("trLchix" ,  &mdst.trLchix );
	//tree->SetBranchAddress("trLchiy" ,  &mdst.trLchiy );
	//tree->SetBranchAddress("trLasym" ,  &mdst.trLasym );
	//tree->SetBranchAddress("trCCest" ,  &mdst.trCCest );
	tree->SetBranchAddress("trMxQ"   ,  &mdst.trMxQ   );
	tree->SetBranchAddress("trL2Q"   ,  &mdst.trL2Q   );
	tree->SetBranchAddress("trL1Q"   ,  &mdst.trL1Q   );
	tree->SetBranchAddress("trL9Q"   ,  &mdst.trL9Q   );
	tree->SetBranchAddress("tofQul"  ,  &mdst.tofQul  );
	//tree->SetBranchAddress("tofb"    ,  &mdst.tofb    );
	tree->SetBranchAddress("tofM"    ,  &mdst.tofM    );
	tree->SetBranchAddress("trdEst"  ,  &mdst.trdEst  );
	tree->SetBranchAddress("hasEcal" ,  &mdst.hasEcal );
	tree->SetBranchAddress("ecalEst" ,  &mdst.ecalEst );
	tree->SetBranchAddress("richRad" ,  &mdst.richRad );
	tree->SetBranchAddress("richPh"  ,  &mdst.richPh  );
	tree->SetBranchAddress("richPrPh",  &mdst.richPrPh);
	tree->SetBranchAddress("richElPh",  &mdst.richElPh);
	tree->SetBranchAddress("hasRich" ,  &mdst.hasRich );
	//tree->SetBranchAddress("richbPr" ,  &mdst.richbPr );
	//tree->SetBranchAddress("richbPi" ,  &mdst.richbPi );
	//tree->SetBranchAddress("richbEl" ,  &mdst.richbEl );
	tree->SetBranchAddress("richMPr" ,  &mdst.richMPr );
	tree->SetBranchAddress("richMPi" ,  &mdst.richMPi );
	tree->SetBranchAddress("richMEl" ,  &mdst.richMEl );
}

#endif // __MinDST_H__
