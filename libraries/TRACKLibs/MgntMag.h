#ifndef __MgntMag_H__
#define __MgntMag_H__
#include <iostream>
#include <cmath>
#include "TString.h"
#include "TMath.h"

//---- Magntic ---//
class Magnetic {
	public :
		Magnetic() : fMagVal(0.) {}
		~Magnetic() {}

		void init() {
			fMagVal = 0.;
			std::fill(fCoo.begin(), fCoo.end(), 0.);
			std::fill(fMagDir.begin(), fMagDir.end(), 0.);
			std::fill(fDevMag.begin(), fDevMag.end(), 0.);
		}

		void print() {
			std::cout << Form("**** Magnetic ****\n");
			std::cout << Form("Coord               : (%14.8f %14.8f %14.8f) [cm]\n", fCoo(0), fCoo(1), fCoo(2));
			std::cout << Form("Magnetic Field      : %14.8f [kGauss]\n", fMagVal);
			std::cout << Form("                      (%14.8f %14.8f %14.8f)\n", fMagDir(0), fMagDir(1), fMagDir(2));
			std::cout << Form("Magnetic Dev Field  : (dF/dP)\n");
			std::cout << Form("                       %14.8f %14.8f %14.8f\n", fDevMag(0, 0), fDevMag(1, 0), fDevMag(2, 0));
			std::cout << Form("                       %14.8f %14.8f %14.8f\n", fDevMag(0, 1), fDevMag(1, 1), fDevMag(2, 1));
			std::cout << Form("                       %14.8f %14.8f %14.8f\n", fDevMag(0, 2), fDevMag(1, 2), fDevMag(2, 2));
			std::cout << Form("\n");
		}

	public :
		SVecD<3> fCoo;    // coord [cm]
		Double_t        fMagVal; // magnetic field [10^3Gauss = 0.1Tesla]
		SVecD<3> fMagDir; // magnetic field unit vector [1]
		SMtxD<3> fDevMag; // derivative magnetic field [x, y, z][derivative x, y, z]
};


//---- MgntMag ----//
#include "MagField.h"
class MgntMag {
	public :
		MgntMag() {}
		~MgntMag() {}
		
		static Bool_t Load();

    static Magnetic & GetMag() { return MagVec; }
    static Magnetic & GetMag(Float_t coo[3], Bool_t isWithDev = false);
		static Magnetic & GetMag(SVecD<3> & coo, Bool_t isWithDev = false);
	
		static SVecD<3> & GetMagFast() { return MagVecFast; }
		static SVecD<3> & GetMagFast(Float_t coo[3]);
		static SVecD<3> & GetMagFast(SVecD<3> & coo);

		static SVecD<3> & Coo() { return MagVec.fCoo; }
		static SVecD<3> & MagDir() { return MagVec.fMagDir; }
		static SMtxD<3> & DevMag() { return MagVec.fDevMag; }
		static Double_t Coo(Short_t i) { return MagVec.fCoo(i); }
		static Double_t MagDir(Short_t i) { return MagVec.fMagDir(i); }
		static Double_t DevMag(Short_t i, Short_t j) { return MagVec.fDevMag(i, j); }

		static Double_t Mag()  { return MagVec.fMagVal; }
    static Double_t MagX() { return (MagVec.fMagVal * MagVec.fMagDir(0)); }
    static Double_t MagY() { return (MagVec.fMagVal * MagVec.fMagDir(1)); }
    static Double_t MagZ() { return (MagVec.fMagVal * MagVec.fMagDir(2)); }
    static Double_t MagDirX() { return MagVec.fMagDir(0); }
    static Double_t MagDirY() { return MagVec.fMagDir(1); }
    static Double_t MagDirZ() { return MagVec.fMagDir(2); }
		static Double_t DevMagXX() { return MagVec.fDevMag(0, 0); }
		static Double_t DevMagXY() { return MagVec.fDevMag(0, 1); }
		static Double_t DevMagXZ() { return MagVec.fDevMag(0, 2); }
		static Double_t DevMagYX() { return MagVec.fDevMag(1, 0); }
		static Double_t DevMagYY() { return MagVec.fDevMag(1, 1); }
		static Double_t DevMagYZ() { return MagVec.fDevMag(1, 2); }
		static Double_t DevMagZX() { return MagVec.fDevMag(2, 0); }
		static Double_t DevMagZY() { return MagVec.fDevMag(2, 1); }
		static Double_t DevMagZZ() { return MagVec.fDevMag(2, 2); }

		static Double_t GetMagFunc(Double_t z);
		static Double_t GetMagFuncINT1(Double_t z);
		static Double_t GetMagFuncINT2(Double_t z);

	protected :
		static constexpr Double_t MagSqrt2      = 1.41421356237309515e+00;
		static constexpr Double_t MagInvSqrt2   = 7.07106781186547462e-01;
		static constexpr Double_t MagInvSqrt2Pi = 3.98942280401432703e-01;
		static constexpr Double_t MagInvSqrtPi  = 5.64189583547756279e-01;
		static constexpr Double_t MagParam[4]   = { 123.21342990, 37.01735081, 20.80309443, 79.25552614 };

    static Bool_t          LoadStatus;
		static MagField *      MagFld;
		static Magnetic        MagVec;
		static SVecD<3> MagVecFast;
};

Bool_t          MgntMag::LoadStatus = false;
MagField *      MgntMag::MagFld = 0;
Magnetic        MgntMag::MagVec;
SVecD<3> MgntMag::MagVecFast;


Bool_t MgntMag::Load() {
  if (LoadStatus) return LoadStatus;
	MagFld = MagField::GetPtr();
	if (!MagFld->GetMap()) {
		bool isWork = true;
		static int magerr = 0;
		if (!MagFld->GetMap() && !magerr) {
			std::string filePath = STR_FMT("%s/v5.00/MagneticFieldMapPermanent_NEW.bin", MGSys::GetEnv("AMSDataDir").c_str());
			//std::string filePath = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/MagneticFieldMapPermanent_NEW.bin";
			if ((MagFld->Read(filePath.c_str())) < 0) {
				std::cerr << "Magnetic Field map not found : " << filePath.c_str() << std::endl;
				magerr = -1;
				isWork = false;
			}
			else {
				std::cout << "MgntMag::Load() Open file : "
					<< filePath.c_str() << std::endl;
			}
			MagFld->SetMagstat(1);
			MagFld->SetScale(1);
		}
		if (!isWork) MagFld = 0;
	}
  if (MagFld != 0) LoadStatus = true;

  return LoadStatus;
}

Magnetic & MgntMag::GetMag(Float_t coo[3], Bool_t isWithDev) {
  if (!Load()) return MagVec;
	MagVec.init();
	Float_t magVec[3] = {0, 0, 0};
	MagFld->GuFld(coo, magVec);
	for (Int_t idir = 0; idir < 3; ++idir) {
		MagVec.fCoo(idir) = coo[idir];	
		MagVec.fMagDir(idir) = magVec[idir];
	}
	MagVec.fMagVal = LA::Mag(MagVec.fMagDir);
	Bool_t vacuum = MGNumc::EqualToZero(MagVec.fMagVal);
	if (!vacuum) MagVec.fMagDir.Unit();

	if (isWithDev && !vacuum) {
		Float_t devMag[2][3] = {0};
		MagFld->TkFld(coo, devMag);
		for (Int_t idir = 0; idir < 3; ++idir) {
			MagVec.fDevMag(idir, 0) = (devMag[0][idir] / MagVec.fMagVal);
			MagVec.fDevMag(idir, 1) = (devMag[1][idir] / MagVec.fMagVal);
		}
	}

	return MagVec;
}

Magnetic & MgntMag::GetMag(SVecD<3> & coo, Bool_t isWithDev) {
	Float_t _coo[3] = { Float_t(coo(0)), Float_t(coo(1)), Float_t(coo(2)) };
	return MgntMag::GetMag(_coo, isWithDev);
}

SVecD<3> & MgntMag::GetMagFast(Float_t coo[3]) {
  if (!Load()) return MagVecFast;
	Float_t magVecFast[3] = {0, 0, 0};
	MagFld->GuFld(coo, magVecFast);
	MagVecFast(0) = magVecFast[0];
	MagVecFast(1) = magVecFast[1];
	MagVecFast(2) = magVecFast[2];
	return MagVecFast;
}

SVecD<3> & MgntMag::GetMagFast(SVecD<3> & coo) {
	Float_t _coo[3] = { Float_t(coo(0)), Float_t(coo(1)), Float_t(coo(2)) };
	return MgntMag::GetMagFast(_coo);
}

Double_t MgntMag::GetMagFunc(Double_t z) {
	Double_t NEG = -1.;
	Double_t normz1 = MagInvSqrt2 * z / MagParam[1];
	Double_t normz2 = MagInvSqrt2 * z / MagParam[3];
	Double_t mag = 
		MagParam[0] * (MagInvSqrt2Pi / MagParam[1]) * TMath::Exp(NEG * normz1 * normz1) +
		MagParam[2] * (MagInvSqrt2Pi / MagParam[3]) * TMath::Exp(NEG * normz2 * normz2);
	return mag;
}

Double_t MgntMag::GetMagFuncINT1(Double_t z) { // reference Z = 0
	Double_t NEG = -1.;
	Double_t HALF = 0.5;
	Double_t normz1 = MagInvSqrt2 * z / MagParam[1];
	Double_t normz2 = MagInvSqrt2 * z / MagParam[3];
	Double_t magINT1 =
		MagParam[0] * HALF * TMath::Erf(normz1) +
		MagParam[2] * HALF * TMath::Erf(normz2);
	return magINT1;
}

Double_t MgntMag::GetMagFuncINT2(Double_t z) { // reference Z = 0
	Double_t NEG = -1.;
	Double_t ONE = 1.;
	Double_t HALF = 0.5;
	Double_t factz1 = MagSqrt2 * MagParam[1];
	Double_t factz2 = MagSqrt2 * MagParam[3];
	Double_t normz1 = z / factz1;
	Double_t normz2 = z / factz2;
	Double_t magINT2 =
		MagParam[0] * (HALF * factz1) * ((normz1 * TMath::Erf(normz1)) + ((TMath::Exp(NEG * normz1 * normz1) - ONE) / MagInvSqrtPi)) + 
		MagParam[2] * (HALF * factz2) * ((normz2 * TMath::Erf(normz2)) + ((TMath::Exp(NEG * normz2 * normz2) - ONE) / MagInvSqrtPi)); 
	return magINT2;
}

#endif // __MgntMag_H__
