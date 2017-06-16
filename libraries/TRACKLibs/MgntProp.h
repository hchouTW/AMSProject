#ifndef __MgntProp_H__
#define __MgntProp_H__

#include "MgntMat.h"
#include "MgntMag.h"
#include "MgntPhySt.h"

// Official TrProp
#include "TrFit.h"

/********************************
 **    MgntProp based on DS    **
 ********************************/

//---- MatPropParam ----//
class MatPropParam {
	public :
		MatPropParam(Bool_t optMscat = false, Bool_t optEngls = false) { setOption(optMscat, optEngls); }
		~MatPropParam() {}

		void init() { setOption(); }

		void setOption(Bool_t optMscat = false, Bool_t optEngls = false) {
			fOptMscat = optMscat;
			fOptEngls = optEngls;
			fVacuum = !(fOptMscat || fOptEngls);

			fMscatDV = ZERO;
			fMscatDW = ZERO;
			fEnglsIN = ZERO;
			fEnglsBR = ZERO;
		}

		void setMscat(Double_t mscatDV = 0, Double_t mscatDW = 0) {
			if (!fOptMscat) return;
			fMscatDV = (!MGNumc::Valid(mscatDV)) ? ZERO : mscatDV;
			fMscatDW = (!MGNumc::Valid(mscatDW)) ? ZERO : mscatDW;
		}
		
		void setEngls(Double_t englsIN = 0, Double_t englsBR = 0) {
			if (!fOptEngls) return;
			fEnglsIN  = (!MGNumc::Valid(englsIN)) ? ZERO : englsIN;
			fEnglsBR  = (!MGNumc::Valid(englsBR) || MGNumc::Compare(englsBR) <= 0) ? ZERO : englsBR;
		}

		void random(MatPhyCalParam & matCal, Bool_t optMscat = false, Bool_t optEngls = false) {
			if (matCal.fVacuum) { init(); return; }
			else setOption(optMscat, optEngls);

			if (fOptMscat) {
				fMscatDV = MGRndm::NormalGaussian();
				fMscatDW = MGRndm::NormalGaussian();
			}

			if (fOptEngls) {
				Double_t englsINLMT = (NEG * matCal.fEnglsIMPV / matCal.fEnglsISGM);
				if (!MGNumc::Valid(englsINLMT)) fEnglsIN = 0.;
				else {
					constexpr Short_t niter = 3;
					Short_t iter = 0;
					fEnglsIN = MGROOT::Rndm::Landau(0.0, 1.0);
					while (fEnglsIN < englsINLMT && iter < niter) { fEnglsIN = MGROOT::Rndm::Landau(0.0, 1.0); iter++; }
					if (iter == niter) fEnglsIN = englsINLMT;
					//--- Testing Gaus ----//
					//fEnglsIN = MGRndm::NormalGaussian();
					//while (fEnglsIN < englsINLMT && iter < niter) { fEnglsIN = MGRndm::NormalGaussian(); iter++; }
					//if (iter == niter) fEnglsIN = englsINLMT;
					//--------//
				}
				Double_t bremslen = matCal.fNumRadLen * Bremsstrahlung_1oln2;
				fEnglsBR = (MGNumc::Compare(bremslen)<=0) ? ZERO : MGRndm::Gamma(bremslen, 1./bremslen)();
			}
		}
		
		void print() {
			std::cout << Form("**** MatPropParam ****\n");
			std::cout << Form("Vacuum : %d\n",     fVacuum);
			if (fVacuum) { std::cout << Form("\n"); return; }
			if (fOptMscat) {
				std::cout << Form("MscatDV : %14.8f\n", fMscatDV);
				std::cout << Form("MscatDW : %14.8f\n", fMscatDW);
			}
			if (fOptEngls) {
				std::cout << Form("EnglsIN : %14.8f\n", fEnglsIN);
				std::cout << Form("EnglsBR : %14.8f\n", fEnglsBR);
			}
			std::cout << Form("\n");
		}

	public :
		Bool_t Vacuum() { return fVacuum; }
		Bool_t OptMscat()  { return fOptMscat; }
		Bool_t OptEngls()  { return fOptEngls; }
	
		Double_t MscatDV() { return fMscatDV; }
		Double_t MscatDW() { return fMscatDW; }
		Double_t EnglsIN() { return fEnglsIN; }
		Double_t EnglsBR() { return fEnglsBR; }

		SVecD<4> Param() { 
			SVecD<4> param;
			if (fVacuum) return param;
			if (fOptMscat) { param(0) = fMscatDV; param(1) = fMscatDW; }
			if (fOptEngls) { param(2) = fEnglsIN; param(3) = fEnglsBR; }
			return param;
		}

		SVecD<4> Chisq(Double_t numRadLen = ZERO) {
			SVecD<4> chisq;
			if (fVacuum) return chisq;
			if (fOptMscat) {
				chisq(0) = (fMscatDV * fMscatDV);
				chisq(1) = (fMscatDW * fMscatDW);
			}
			if (fOptEngls) {
				Double_t bremslen = numRadLen * Bremsstrahlung_1oln2;
				Double_t englsBR = (MGNumc::Compare(fEnglsBR, LIMIT) > 0) ? fEnglsBR : LIMIT;
				chisq(2) = (fEnglsIN + std::exp(-fEnglsIN) - ONE);
				chisq(3) = ((englsBR) * bremslen + (ONE - bremslen) * (std::log(englsBR) - std::log(LIMIT)));
			}
			return chisq;
		}

		SVecD<4> Grad(Double_t numRadLen = ZERO) {
			if (numRadLen < ZERO || !MGNumc::Valid(numRadLen)) numRadLen = ZERO;
			SVecD<4> grad;
			if (fVacuum) return grad;
			if (fOptMscat) { 
				grad(0) = fMscatDV; 
				grad(1) = fMscatDW; 
			}
			if (fOptEngls) { 
				Double_t bremslen = numRadLen * Bremsstrahlung_1oln2;
				Double_t englsBR = (MGNumc::Compare(fEnglsBR, LIMIT) > 0) ? fEnglsBR : LIMIT;
				grad(2) = HALF * (ONE - std::exp(-fEnglsIN)); 
				//--- Testing Gaus ----//
				//grad(2) = fEnglsIN; 
				//--------//
				grad(3) = (bremslen + (ONE - bremslen) / englsBR);
			}
			return grad;
		}
		
		SMtxSymD<4> Cov(Double_t numRadLen = ZERO) {
			if (numRadLen < ZERO || !MGNumc::Valid(numRadLen)) numRadLen = ZERO;
			SMtxSymD<4> cov;
			if (fVacuum) return cov;
			if (fOptMscat) { 
				cov(0, 0) = ONE; 
				cov(1, 1) = ONE; 
			}
			if (fOptEngls) { 
				Double_t bremslen = numRadLen * Bremsstrahlung_1oln2;
				Double_t englsBR = (MGNumc::Compare(fEnglsBR, LIMIT) > 0) ? fEnglsBR : LIMIT;
				cov(2, 2) = HALF * std::exp(-fEnglsIN); 
				//--- Testing Gaus ----//
				//cov(2, 2) = ONE; 
				//--------//
				cov(3, 3) = std::fabs((bremslen - ONE) / englsBR / englsBR);
			}
			return cov;
		}

	protected :
		static constexpr Double_t ZERO = 0.;
		static constexpr Double_t ONE  = 1.;
		static constexpr Double_t HALF = 0.5;
		static constexpr Double_t NEG  = -1;
		static constexpr Double_t LIMIT = 1e-7;
		// Bremsstrahlung
		static constexpr Double_t Bremsstrahlung_1oln2 = 1.44269504088896339e+00;

	protected :
		Bool_t fVacuum;
		Bool_t fOptMscat;
		Bool_t fOptEngls;
	
	public :	
		Double_t fMscatDV;
		Double_t fMscatDW;
		Double_t fEnglsIN;
		Double_t fEnglsBR;
};


//---- PhyJb ----//
class PhyJb {
	public :
		enum MatrixKind {
			Zero = 0, Identity = 1
		};
	
		static constexpr Short_t GDim = 5;
		static constexpr Short_t LDim = 4;

	public :
		PhyJb(MatrixKind kind = PhyJb::Zero, Bool_t vacuum = true) { init(kind, vacuum); }
		~PhyJb() {}

		inline void init(MatrixKind kind = PhyJb::Zero, Bool_t vacuum = true) { 
			fVacuum = vacuum; 
			fJacbL = SMtxD<5, 4>(); 
			switch (kind) {
				case Zero :
					fJacbG = SMtxD<5>(); break;
				case Identity :
					fJacbG = SMtxId(); break;
				default :
					break;
			}
			fDeltaZ = 0.;
			fNumRadLen = 0.;
		}

		inline Bool_t Vacuum() { return fVacuum; }
		inline SMtxD<5>    & G() { return fJacbG; }
		inline SMtxD<5, 4> & L() { return fJacbL; }
    inline Double_t & G(Short_t i, Short_t j) { return fJacbG(i, j); }
		inline Double_t & L(Short_t i, Short_t j) { return fJacbL(i, j); }
		inline Double_t & DeltaZ() { return fDeltaZ; }
		inline Double_t & NumRadLen() { return fNumRadLen; }

		inline SMtxD<2, 5> SubJbGXY() { return this->fJacbG.Sub<SMtxD<2, 5>>(0, 0); }
		inline SMtxD<2, 4> SubJbLXY() { return this->fJacbL.Sub<SMtxD<2, 4>>(0, 0); }
		inline SMtxD<2>    SubJbLXYWithMscat() { return this->fJacbL.Sub<SMtxD<2>>(0, 0); }
		inline SMtxD<2>    SubJbLXYWithEngls() { return this->fJacbL.Sub<SMtxD<2>>(0, 2); }

		void multiply(PhyJb & jb, Bool_t onlyG = false) {
			this->fDeltaZ += jb.DeltaZ();
			if (!onlyG) {
				if (!jb.Vacuum()) this->L() += (this->G() * jb.L());
				this->fVacuum = (this->Vacuum() && jb.Vacuum());
				this->fNumRadLen += jb.NumRadLen();
			}
			this->G() *= jb.G();
		}

		void multiplied(PhyJb & jb, Bool_t onlyG = false) {
			this->fDeltaZ += jb.DeltaZ();
			if (!onlyG) {
				if (!this->Vacuum()) this->L() = (jb.G() * this->L());
				if (!jb.Vacuum()) this->L() += jb.L();
				this->fVacuum = (this->Vacuum() && jb.Vacuum());
				this->fNumRadLen += jb.NumRadLen();
			}
			this->G() = (jb.G() * this->G());
		}

		PhySt transfer(PhySt & phySt, MatPropParam * matPar = nullptr) { // loss one term --- ion-loss peak
			SVecD<5> stOLD;
			stOLD(0) = phySt.X();    stOLD(1) = phySt.Y();
			stOLD(2) = phySt.DirX(); stOLD(3) = phySt.DirY();
			stOLD(4) = phySt.InvEta();
			SVecD<5> && stNEW = this->G() * stOLD;
			if (matPar != nullptr && !matPar->Vacuum() && !this->Vacuum()) stNEW += (this->L() * matPar->Param());
			PhySt state = phySt;
			state.setSpatialWithCos(stNEW(0), stNEW(1), phySt.Z() + fDeltaZ, stNEW(2), stNEW(3), phySt.DirZ(), false);
			state.setInvEta(stNEW(4));
			return state;
		}
		
		void print() {
			std::cout << Form("**** PhyJb ****\n");
			std::cout << Form("DeltaZ    %14.8f\n", fDeltaZ);
			std::cout << Form("NumRadLen %14.8f\n", fNumRadLen);
			TString var[9] = { "x", "y", "ux", "uy", "1/eta", "mscatDV", "mscatDW", "englsIN", "englsBR" };
			if (fVacuum)
				std::cout << Form("%14s %14s %14s %14s %14s %14s\n", "", 
				var[0].Data(), var[1].Data(), var[2].Data(), var[3].Data(), var[4].Data()); 
			else
				std::cout << Form("%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n", "", 
					var[0].Data(), var[1].Data(), var[2].Data(), var[3].Data(), var[4].Data(), 
					var[5].Data(), var[6].Data(), var[7].Data(), var[8].Data());

			for (Int_t row = 0; row < PhyJb::GDim; ++row) {
				std::cout << Form("%14s ", var[row].Data());
				for (Int_t col = 0; col < PhyJb::GDim; ++col)
					std::cout << Form("%14.8f ", fJacbG(row, col));
				if (!fVacuum)
					for (Int_t col = 0; col < PhyJb::LDim; ++col)
						std::cout << Form("%14.8f ", fJacbL(row, col));
				std::cout << Form("\n");
			}
			std::cout << Form("\n");
		}

	private :
		Bool_t             fVacuum;
		SMtxD<5>    fJacbG;
		SMtxD<5, 4> fJacbL;
		Double_t           fDeltaZ;
		Double_t           fNumRadLen;
};


//---- MgntProp ----//
class MgntProp {
	public :
    enum MethodOption {
			kEuler = 1, // first order
			kEulerHeun = 2, // second order
			kRungeKuttaNystrom = 4 // fourth order
		};
    
		static void   SetMethod(MgntProp::MethodOption method = MgntProp::kRungeKuttaNystrom) { Method = method; }
    static Bool_t CheckMethod(MgntProp::MethodOption method) { return (method == Method); }

	private :
    enum EnvironmentOption {
      kVacuum = 0,
      kMaterial = 1
    };
		
		static void   SetEnvironment(MgntProp::EnvironmentOption env = MgntProp::kVacuum) { Environment = env; }
    static Bool_t CheckEnvironment(MgntProp::EnvironmentOption env) { return (env == Environment); }

	public :
		MgntProp() {}
		~MgntProp() {}

		static Bool_t Prop(Double_t step, PhySt & phySt, PhyJb * phyJb = nullptr, MatPropParam * matPar = nullptr);
		static Bool_t PropToZ(Double_t zcoo, PhySt & phySt, PhyJb * phyJb = nullptr, MatPropParam * matPar = nullptr);

		// Toy Monte Carlo
		static Bool_t Prop_MonteCarlo(Double_t step, PhySt & phySt, Bool_t optMscat = false, Bool_t optEngls = false);
		static Bool_t PropToZ_MonteCarlo(Double_t zcoo, PhySt & phySt, Bool_t optMscat = false, Bool_t optEngls = false);
	
		// Official
		static Bool_t PropToZ_Official(Double_t zcoo, PhySt & phySt);

		// Analytic
		static Bool_t PropToZ_Analytic(Double_t zcoo, PhySt & phySt, SMtxD<5> * phyJb = nullptr);

	private :
		// Step Length
		// step < 0, backward trace
		// step > 0, forward trace
		static Double_t GetPropStep(PhySt & phySt, Double_t sign);
		static Double_t GetStep(PhySt & phySt, Double_t resStep);
		static Double_t GetStepToZ(PhySt & phySt, Double_t resStepZ);

		static Bool_t Prop_EulerMethod(Double_t step, PhySt & phySt, PhyJb * phyJb = nullptr); 
		static Bool_t Prop_EulerHeunMethod(Double_t step, PhySt & phySt, PhyJb * phyJb = nullptr);
		static Bool_t Prop_RungeKuttaNystromMethod(Double_t step, PhySt & phySt, PhyJb * phyJb = nullptr);

	private :	
		static constexpr Short_t PX = 0;
		static constexpr Short_t PY = 1;
		static constexpr Short_t PZ = 2;
		static constexpr Short_t UX = 3;
		static constexpr Short_t UY = 4;
		static constexpr Short_t UZ = 5;
		static constexpr Short_t GB = 6;
		
		static constexpr Short_t JPX = 0;
		static constexpr Short_t JPY = 1;
		static constexpr Short_t JUX = 2;
		static constexpr Short_t JUY = 3;
		static constexpr Short_t JGB = 4;

		static constexpr Double_t NEG  = -1.;
		static constexpr Double_t ZERO = 0.;
		static constexpr Double_t ONE  = 1.;
		static constexpr Double_t TWO  = 2.;
		static constexpr Double_t THREE = 3.;
		static constexpr Double_t SIX   = 6.;
		static constexpr Double_t EIGHT = 8.;
		static constexpr Double_t HALF  = 0.5;
		static constexpr Double_t ONE_SIXTH = 1. / 6.;
		static constexpr Double_t ONE_EIGHTH = 1. / 8.;

		static constexpr Double_t TUNE_ETA  = 3.00; // (min-eta)
		static constexpr Double_t TUNE_STEP = 5e-3; // (angular threshold)
		static constexpr Double_t TUNE_MAGF = 5e-2; // (mini magnetic field) [kG]
		static constexpr Double_t TUNE_MAT  = 5e-2; // (radiation length)
		static constexpr Double_t MIN_STEP  = 8.;   // [cm]
		static constexpr Double_t INT_STEP  = 15.;  // [cm]
		static constexpr Double_t MAX_STEP  = 30.;  // [cm]
		static constexpr Double_t MIN_MOM   = 1e-4; // [GeV]
		static constexpr Double_t MIN_CONT  = 1e-4; // [cm]
		static constexpr Long64_t MAX_ITER  = 100;

	public :
		static constexpr Double_t PI = 3.14159265358979312;
		static constexpr Double_t PROP_FACT = 2.99792458e-04; // propragate factor [(GeV/c)(kG^-1)(cm^-1)]
		static constexpr Double_t BetaLimit = 0.3; // energy loss limit
		static constexpr Double_t IEtaLimit = 3.179797e+00;	
	
	private :
		static enum MethodOption      Method;
    static enum EnvironmentOption Environment;
		static MatPropParam           MatPar;
		static MatPhyCalParam         MatCal;
};

MgntProp::MethodOption      MgntProp::Method      = MgntProp::kRungeKuttaNystrom;
MgntProp::EnvironmentOption MgntProp::Environment = MgntProp::kVacuum;
MatPropParam                MgntProp::MatPar;
MatPhyCalParam              MgntProp::MatCal;


//---- OrthCoord ----//
// Orthogonal Coordinate based on seed vector (1, 0, 0)
// Taget is particle direction
// Axis T := V x W
class OrthCoord {
	public :
		OrthCoord(PhySt & phySt, Bool_t isWithDev = false);
		~OrthCoord() {}

		inline void init() { fTaget = SVecD<3>(); fAxisV = SVecD<3>(); fAxisW = SVecD<3>(); fAxisVDev = SMtxD<2>(); fAxisWDev = SMtxD<2>(); }
	
		inline SVecD<3> & T() { return fTaget; }
		inline SVecD<3> & V() { return fAxisV; }
		inline SVecD<3> & W() { return fAxisW; }
		inline SMtxD<2> & VDev() { return fAxisVDev; }
		inline SMtxD<2> & WDev() { return fAxisWDev; }
		
		inline Double_t & T(Short_t i) { return fTaget(i); }
		inline Double_t & V(Short_t i) { return fAxisV(i); }
		inline Double_t & W(Short_t i) { return fAxisW(i); }
		inline Double_t & VDev(Short_t i, Short_t j) { return fAxisVDev(i, j); }
		inline Double_t & WDev(Short_t i, Short_t j) { return fAxisWDev(i, j); }


	private :
		static const SVecD<3> ORTH_SEED; // orthogonal seed vector
		SVecD<3> fTaget;
		SVecD<3> fAxisV;
		SVecD<3> fAxisW;
		SMtxD<2> fAxisVDev;
		SMtxD<2> fAxisWDev;
		
		static constexpr Short_t  X = 0;
		static constexpr Short_t  Y = 1;
		static constexpr Short_t  Z = 2;
		static constexpr Double_t ZERO = 0.;
		static constexpr Double_t ONE  = 1.;
		static constexpr Double_t TWO  =  2.;
		static constexpr Double_t NEG  = -1.;
};

const SVecD<3> OrthCoord::ORTH_SEED(1, 0, 0);


//---- DevStatusFunc ----//
// Status (x, y, z, 1/eta) --- eta := ChrgSign * GammaBeta
// 1st Derivative (ds) of Status Function
// 2st Derivative (ds) of Status Function
// devFunc(x, y, z, ux, uy, uz, 1/eta)
class DevStatusFunc {
	public :
		DevStatusFunc(PhySt & phySt, MatPropParam * matPar = nullptr, MatPhyCalParam * matCal = nullptr);
		~DevStatusFunc() {}
	
		inline void init() { fVacuum = true; fDevFunc = SVecD<7>(); }
		inline SVecD<7> & operator() () { return fDevFunc; }
    inline Double_t & operator() (Short_t i) { return fDevFunc(i); }

		void print() {
			std::cout << Form("**** DevStatusFunc ****\n");
			TString var[7] = { "x", "y", "z", "ux", "uy", "uz", "1/eta" };
			std::cout << Form("%14s %14s %14s %14s %14s %14s %14s\n", 
				var[0].Data(), var[1].Data(), var[2].Data(), var[3].Data(), var[4].Data(), var[5].Data(), var[6].Data()); 
			for (Int_t col = 0; col < 7; ++col)
				std::cout << Form("%14.8f ", fDevFunc(col));
			std::cout << Form("\n\n");
		}

	private :
		Bool_t          fVacuum;
		SVecD<7>	fDevFunc;

		static constexpr Short_t  X = 0;
		static constexpr Short_t  Y = 1;
		static constexpr Short_t  Z = 2;
		static constexpr Double_t ZERO = 0.;
		static constexpr Double_t ONE  = 1.;
		static constexpr Double_t NEG  = -1.;
		static constexpr Double_t HALF = 0.5;
		static constexpr Double_t INV_SQRT_TWELVE = 2.88675134594812866e-01;
};


//---- DevParamFunc ----//
// Status (x, y, z, 1/eta) --- eta := ChrgSign * GammaBeta
// 1st Derivative of Status Function Relation with Parameters
// 2st Derivative of Status Function Relation with Parameters
// devMtxG(x, y, z, ux, uy, uz, 1/eta)
// devMtxL(mscatL, mscatD, mscatC, englsI)
class DevParamFunc {
	public :
		static constexpr Short_t GDim = 5;
		static constexpr Short_t LDim = 4;

	public :
		DevParamFunc(PhySt & phySt, MatPropParam * matPar = nullptr, MatPhyCalParam * matCal = nullptr);
		~DevParamFunc() {}
	
		inline void init() { fVacuum = true; fDevMtxG = SMtxD<5>(); fDevMtxL = SMtxD<5, 4>(); }
		inline SMtxD<5>    & G() { return fDevMtxG; }
		inline SMtxD<5, 4> & L() { return fDevMtxL; }
    inline Double_t & G(Short_t i, Short_t j) { return fDevMtxG(i, j); }
    inline Double_t & L(Short_t i, Short_t j) { return fDevMtxL(i, j); }

		void print() {
			std::cout << Form("**** DevParamFunc ****\n");
			TString var[9] = { "x", "y", "ux", "uy", "1/eta", "mscatDV", "mscatDW", "englsIN", "englsBR" };
			if (fVacuum)
				std::cout << Form("%14s %14s %14s %14s %14s %14s\n", "", 
				var[0].Data(), var[1].Data(), var[2].Data(), var[3].Data(), var[4].Data()); 
			else
				std::cout << Form("%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n", "", 
					var[0].Data(), var[1].Data(), var[2].Data(), var[3].Data(), var[4].Data(), 
					var[5].Data(), var[6].Data(), var[7].Data(), var[8].Data());

			for (Int_t row = 0; row < DevParamFunc::GDim; ++row) {
				std::cout << Form("%14s ", var[row].Data());
				for (Int_t col = 0; col < DevParamFunc::GDim; ++col) 
					std::cout << Form("%14.8f ", fDevMtxG(row, col));
				if (!fVacuum)
					for (Int_t col = 0; col < DevParamFunc::LDim; ++col)
						std::cout << Form("%14.8f ", fDevMtxL(row, col));
				std::cout << Form("\n");
			}
			std::cout << Form("\n");
		}

	private :
		Bool_t             fVacuum;
		SMtxD<5>    fDevMtxG;
		SMtxD<5, 4> fDevMtxL;
		
		static constexpr Short_t  X = 0;
		static constexpr Short_t  Y = 1;
		static constexpr Short_t  Z = 2;
		static constexpr Double_t ZERO  = 0.;
		static constexpr Double_t ONE   = 1.;
		static constexpr Double_t TWO   = 2.;
		static constexpr Double_t THREE = 3.;
		static constexpr Double_t FOUR  = 4.;
		static constexpr Double_t FIVE  = 5.;
		static constexpr Double_t NEG   = -1.;
		static constexpr Double_t HALF  = 0.5;
		static constexpr Double_t INV_SQRT_TWELVE = 2.88675134594812866e-01;
};


#endif // __MgntProp_H__
