#ifndef __ROOTLibs_MGMath_H__
#define __ROOTLibs_MGMath_H__


#include <Math/Functor.h>
#include <Math/BrentMinimizer1D.h>
#include <Math/GSLMinimizer1D.h>
#include <Math/GSLMinimizer.h>
#include <Math/GSLSimAnMinimizer.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>

namespace MGROOT {
	namespace Math {
		// One-Dimensional Minimization
		using ROOT::Math::Functor1D;
		using ROOT::Math::BrentMinimizer1D;
		using ROOT::Math::GSLMinimizer1D;

		// Multi-Dimensional Minimization
		using ROOT::Math::Functor;
		using ROOT::Math::Minimizer;
		using ROOT::Math::GSLMinimizer;
		using ROOT::Math::GSLSimAnMinimizer;
		using ROOT::Minuit2::Minuit2Minimizer;
	}
}


#include <TArray.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TArrayD.h>
#include <TVector.h>
#include <TMatrix.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
namespace MGROOT {
	// TArray 
	using TArrI = TArrayI;
	using TArrF = TArrayF;
	using TArrD = TArrayD;
	
	// TVector TMatrix Package
	using TVecO = TVectorT<Bool_t>;
	using TVecI = TVectorT<Int_t>;
	using TVecF = TVectorT<Float_t>;
	using TVecD = TVectorT<Double_t>;
	
	using TMtxO = TMatrixT<Bool_t>;
	using TMtxI = TMatrixT<Int_t>;
	using TMtxF = TMatrixT<Float_t>;
	using TMtxD = TMatrixT<Double_t>;
			
	using TMtxSymO = TMatrixTSym<Bool_t>;
	using TMtxSymI = TMatrixTSym<Int_t>;
	using TMtxSymF = TMatrixTSym<Float_t>;
	using TMtxSymD = TMatrixTSym<Double_t>;

	using TMtxSparseO = TMatrixTSparse<Bool_t>;
	using TMtxSparseI = TMatrixTSparse<Int_t>;
	using TMtxSparseF = TMatrixTSparse<Float_t>;
	using TMtxSparseD = TMatrixTSparse<Double_t>;
	
	// SVector SMatrix Package
	template <unsigned int D>
	using SVecO = ROOT::Math::SVector<Bool_t, D>;
	
    template <unsigned int D>
	using SVecI = ROOT::Math::SVector<Int_t, D>;
	
	template <unsigned int D>
	using SVecF = ROOT::Math::SVector<Float_t, D>;
	
	template <unsigned int D>
	using SVecD = ROOT::Math::SVector<Double_t, D>;
	
	template <unsigned int D1, unsigned int D2 = D1>
	using SMtxO = ROOT::Math::SMatrix<Bool_t, D1, D2, ROOT::Math::MatRepStd<Bool_t, D1, D2>>;
	
    template <unsigned int D1, unsigned int D2 = D1>
	using SMtxI = ROOT::Math::SMatrix<Int_t, D1, D2, ROOT::Math::MatRepStd<Int_t, D1, D2>>;
	
	template <unsigned int D1, unsigned int D2 = D1>
	using SMtxF = ROOT::Math::SMatrix<Float_t, D1, D2, ROOT::Math::MatRepStd<Float_t, D1, D2>>;
	
	template <unsigned int D1, unsigned int D2 = D1>
	using SMtxD = ROOT::Math::SMatrix<Double_t, D1, D2, ROOT::Math::MatRepStd<Double_t, D1, D2>>;
	
	template <unsigned int D>
	using SMtxSymO = ROOT::Math::SMatrix<Bool_t, D, D, ROOT::Math::MatRepSym<Bool_t, D>>;

	template <unsigned int D>
	using SMtxSymI = ROOT::Math::SMatrix<Int_t, D, D, ROOT::Math::MatRepSym<Int_t, D>>;
	
	template <unsigned int D>
	using SMtxSymF = ROOT::Math::SMatrix<Float_t, D, D, ROOT::Math::MatRepSym<Float_t, D>>;
	
	template <unsigned int D>
	using SMtxSymD = ROOT::Math::SMatrix<Double_t, D, D, ROOT::Math::MatRepSym<Double_t, D>>;
	
	using SMtxId = ROOT::Math::SMatrixIdentity;
    
    static const SMtxD<2> SMtxId2D = SMtxId();
    static const SMtxD<3> SMtxId3D = SMtxId();
    static const SMtxD<4> SMtxId4D = SMtxId();
    static const SMtxD<5> SMtxId5D = SMtxId();
    
    static const SMtxSymD<2> SMtxIdSym2D = SMtxId();
    static const SMtxSymD<3> SMtxIdSym3D = SMtxId();
    static const SMtxSymD<4> SMtxIdSym4D = SMtxId();
    static const SMtxSymD<5> SMtxIdSym5D = SMtxId();
	
	namespace LA { // Linear Algebra
		// SVector Template Function 
		using ROOT::Math::Cross;
		using ROOT::Math::Dot;
		using ROOT::Math::Mag2;
		using ROOT::Math::Mag;
		using ROOT::Math::Lmag2;
		using ROOT::Math::Lmag;
		using ROOT::Math::Unit;
		using ROOT::Math::TensorProd;
		using ROOT::Math::fabs;
		using ROOT::Math::sqr;
		using ROOT::Math::sqrt;
        using ROOT::Math::Unit;
        using ROOT::Math::TensorProd;
		
		// SMatrix Template Function 
		using ROOT::Math::Transpose;
		using ROOT::Math::Similarity;
		using ROOT::Math::SimilarityT;
		using ROOT::Math::fabs;
		using ROOT::Math::sqr;
		using ROOT::Math::sqrt;
	}
}


#include <TRandom3.h>
namespace MGROOT {
	namespace Rndm {
		TRandom3 rndmEngTR3(0);
		
		inline Double_t Uniform  (Double_t x1 = 0.0, Double_t x2 = 1.0)      { return rndmEngTR3.Uniform(x1, x2); }
		inline Double_t Gaus     (Double_t mean = 0.0, Double_t sigma = 1.0) { return rndmEngTR3.Gaus(mean, sigma); }
		inline Double_t Landau   (Double_t mpv = 0.0, Double_t sigma = 1.0)  { return rndmEngTR3.Landau(mpv, sigma); }
		inline Double_t Exp      (Double_t tau = 1.0)                        { return rndmEngTR3.Exp(tau); }
		inline Int_t    Binomial (Int_t ntot = 1, Double_t prob = 0.5)       { return rndmEngTR3.Binomial(ntot, prob); }
		inline Int_t    Poisson  (Double_t mean = 1.0)                       { return rndmEngTR3.Poisson(mean); }
	}
}


#endif // __ROOTLibs_MGMath_H__
