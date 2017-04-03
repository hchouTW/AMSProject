#ifndef __MgntPhySt_H__
#define __MgntPhySt_H__

#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "TString.h"

class PhyStSpace {
  public :
    PhyStSpace() { init(); }
    ~PhyStSpace() {}

    void init() { fPos = MtxLB::SVecD<3>(); fDir = MtxLB::SVecD<3>(); }

    void setSpatialPos(MtxLB::SVecD<3> & pos) { fPos = pos; }
    void setSpatialPos(Double_t x, Double_t y, Double_t z) { fPos(0) = x; fPos(1) = y; fPos(2) = z; }

    void setSpatialCos(MtxLB::SVecD<3> & dir) { fDir = dir; fDir.Unit(); }
    void setSpatialCos(Double_t ux, Double_t uy, Double_t uz, Bool_t isNorm = true) { 
      if (isNorm) {
				fDir(0) = ux; fDir(1) = uy; fDir(2) = uz;
			}
			else {
      	Double_t signz = NGNumc::Compare(uz);
				Double_t uz_Val = signz * std::sqrt(1. - ux * ux - uy * uy);
				if (!NGNumc::Valid(uz_Val)) uz_Val = uz;
				fDir(0) = ux; fDir(1) = uy; fDir(2) = uz_Val;
			}
      fDir.Unit();
		}
    
		void setSpatialTan(Double_t tx, Double_t ty, Double_t signz = -1) {
      Double_t norm = NGNumc::Compare(signz) / std::sqrt(tx * tx + ty * ty + 1.);
      fDir(0) = norm * tx; fDir(1) = norm * ty; fDir(2) = norm;
    }

    void setSpatial(MtxLB::SVecD<3> & pos, MtxLB::SVecD<3> & dir) { fPos = pos; fDir = dir; fDir.Unit(); }

    void setSpatialWithCos(Double_t px, Double_t py, Double_t pz, Double_t ux, Double_t uy, Double_t uz, Bool_t isNorm = true) {
			setSpatialPos(px, py, pz);
			setSpatialCos(ux, uy, uz, isNorm);
		}

    void setSpatialWithTan(Double_t px, Double_t py, Double_t pz, Double_t tx, Double_t ty, Double_t signz = -1) {
			setSpatialPos(px, py, pz);
			setSpatialTan(tx, ty, signz);
		}

    void setSpatial(Double_t state[6]) {
			fPos(0) = state[0]; fPos(1) = state[1]; fPos(2) = state[2];
			fDir(3) = state[3]; fDir(4) = state[4]; fDir(5) = state[5];
			fDir.Unit();
		}

    void setSpatial(Double_t pos[3], Double_t dir[3]) {
			fPos(0) = pos[0]; fPos(1) = pos[1]; fPos(2) = pos[2];
			fDir(0) = dir[0]; fDir(1) = dir[1]; fDir(2) = dir[2];
      fDir.Unit();
		}

	public :
		Bool_t vaild() {
			Bool_t vaild = (
				NGNumc::Valid(fPos(0)) && 
				NGNumc::Valid(fPos(1)) && 
				NGNumc::Valid(fPos(2)) && 
				NGNumc::Valid(fDir(0)) && 
				NGNumc::Valid(fDir(1)) && 
				NGNumc::Valid(fDir(2)) 
			);
			if (!vaild) init();
			return vaild;
		}
		
		void print() {
			std::cout << Form("**** PhyStSpace ****\n");
			std::cout << Form("Position  : (%14.8f %14.8f %14.8f) [cm]\n", fPos(0), fPos(1), fPos(2));
			std::cout << Form("Direction : (%14.8f %14.8f %14.8f)\n", fDir(0), fDir(1), fDir(2));
			std::cout << Form("\n");
		}

  public :
		MtxLB::SVecD<3> & Pos() { return fPos; }
		Double_t & Pos(Short_t idx) { return fPos(idx); }
		Double_t X() { return fPos(0); }
		Double_t Y() { return fPos(1); }
		Double_t Z() { return fPos(2); }

		MtxLB::SVecD<3> & Dir() { return fDir; }
		Double_t & Dir(Short_t idx) { return fDir(idx); }
		Double_t TanX() { return fDir(0)/fDir(2); }
		Double_t TanY() { return fDir(1)/fDir(2); }
		Double_t DirX() { return fDir(0); }
		Double_t DirY() { return fDir(1); }
		Double_t DirZ() { return fDir(2); }

  protected :
    MtxLB::SVecD<3> fPos; // position
		MtxLB::SVecD<3> fDir; // direction

	public :
		static void SortSpace(std::vector<PhyStSpace> & phyStSpace, Bool_t dirType = true) {
			if (phyStSpace.size() < 2) return;
			if (dirType) std::sort(phyStSpace.begin(), phyStSpace.end(), PhyStSpace::SortFromMaxToMin());
			else         std::sort(phyStSpace.begin(), phyStSpace.end(), PhyStSpace::SortFromMinToMax());
		}
	
	protected :
		struct SortFromMinToMax {
			Bool_t operator() (const PhyStSpace & physs1, const PhyStSpace & physs2) {
				if (physs1.fPos(2) < physs2.fPos(2)) return true;
				else return false;
			}
			Bool_t operator() (PhyStSpace * physs1, PhyStSpace * physs2) {
				if (physs1->fPos(2) < physs2->fPos(2)) return true;
				else return false;
			}
		};

		struct SortFromMaxToMin {
			Bool_t operator() (const PhyStSpace & physs1, const PhyStSpace & physs2) {
				if (physs1.fPos(2) > physs2.fPos(2)) return true;
				else return false;
			}
			Bool_t operator() (PhyStSpace * physs1, PhyStSpace * physs2) {
				if (physs1->fPos(2) > physs2->fPos(2)) return true;
				else return false;
			}
		};
};


class PhySt : public PhyStSpace {
  public :
    enum ParticleList {
      Photon, Electron, Positron, Muon, 
			PionPlus, PionMinus,KaonPlus, KaonMinus,
			Proton, Antiproton, Helium3, Antihelium3, Helium4, Antihelium4,
			Lithium6, Lithium7,
			Beryllium7, Beryllium9, Beryllium10,
			Boron10, Boron11,
			Carbon12, Carbon13, Carbon14,
			Nitrogen13, Nitrogen14, Nitrogen15,
			Oxygen16, Oxygen17, Oxygen18
    };

		static Bool_t CheckPart(enum ParticleList part, PhySt & phySt) { return (part == phySt.fPart); }

	public :
		PhySt(enum ParticleList part = PhySt::Proton) { init(part); }
		~PhySt() {}

		void init(enum ParticleList part = PhySt::Proton) { fPartName = ""; fStatus = MtxLB::SVecD<5>(); fPos = MtxLB::SVecD<3>(); fDir = MtxLB::SVecD<3>(); setChrgMass(part); }

    void setChrgMass(enum ParticleList part = PhySt::Proton) {
			// Atomic mass unit u = 0.931494095 GeV/c^2
      switch (part) {
        //case  : fPartName = ""; fPart = PhySt::; fStatus(0) = ; fStatus(1) = ; break;
        case Photon      : fPartName = "Photon"     ; fPart = PhySt::Photon     ; fStatus(0) =   0; fStatus(1) =  0.000000000; break;
        case Electron    : fPartName = "Electron"   ; fPart = PhySt::Electron   ; fStatus(0) =  -1; fStatus(1) =  0.000510999; break;
        case Positron    : fPartName = "Positron"   ; fPart = PhySt::Positron   ; fStatus(0) =   1; fStatus(1) =  0.000510999; break;
        case Muon        : fPartName = "Muon"       ; fPart = PhySt::Muon       ; fStatus(0) =  -1; fStatus(1) =  0.105658367; break;
        case PionPlus    : fPartName = "PionPlus"   ; fPart = PhySt::PionPlus   ; fStatus(0) =   1; fStatus(1) =  0.139570180; break;
        case PionMinus   : fPartName = "PionMinus"  ; fPart = PhySt::PionMinus  ; fStatus(0) =  -1; fStatus(1) =  0.139570180; break;
        case KaonPlus    : fPartName = "KaonPlus"   ; fPart = PhySt::KaonPlus   ; fStatus(0) =   1; fStatus(1) =  0.493677000; break;
        case KaonMinus   : fPartName = "KaonMinus"  ; fPart = PhySt::KaonMinus  ; fStatus(0) =  -1; fStatus(1) =  0.493677000; break;
        case Proton      : fPartName = "Proton"     ; fPart = PhySt::Proton     ; fStatus(0) =   1; fStatus(1) =  0.938272297; break;
        case Antiproton  : fPartName = "Antiproton" ; fPart = PhySt::Antiproton ; fStatus(0) =  -1; fStatus(1) =  0.938272297; break;
        case Helium3     : fPartName = "Helium3"    ; fPart = PhySt::Helium3    ; fStatus(0) =   2; fStatus(1) =  2.809413500; break;
        case Antihelium3 : fPartName = "Antihelium3"; fPart = PhySt::Antihelium3; fStatus(0) =  -2; fStatus(1) =  2.809413500; break;
        case Helium4     : fPartName = "Helium4"    ; fPart = PhySt::Helium4    ; fStatus(0) =   2; fStatus(1) =  3.727379240; break;
        case Antihelium4 : fPartName = "Antihelium4"; fPart = PhySt::Antihelium4; fStatus(0) =  -2; fStatus(1) =  3.727379240; break;
        case Lithium6    : fPartName = "Lithium6"   ; fPart = PhySt::Lithium6   ; fStatus(0) =   3; fStatus(1) =  5.603051363; break;
        case Lithium7    : fPartName = "Lithium7"   ; fPart = PhySt::Lithium7   ; fStatus(0) =   3; fStatus(1) =  6.535366807; break;
        case Beryllium7  : fPartName = "Beryllium7" ; fPart = PhySt::Beryllium7 ; fStatus(0) =   4; fStatus(1) =  6.536228700; break;
        case Beryllium9  : fPartName = "Beryllium9" ; fPart = PhySt::Beryllium9 ; fStatus(0) =   4; fStatus(1) =  8.394794503; break;
        case Beryllium10 : fPartName = "Beryllium10"; fPart = PhySt::Beryllium10; fStatus(0) =   4; fStatus(1) =  9.327547622; break;
        case Boron10     : fPartName = "Boron10"    ; fPart = PhySt::Boron10    ; fStatus(0) =   5; fStatus(1) =  9.326991682; break;
        case Boron11     : fPartName = "Boron11"    ; fPart = PhySt::Boron11    ; fStatus(0) =   5; fStatus(1) = 10.255102976; break;
        case Carbon12    : fPartName = "Carbon12"   ; fPart = PhySt::Carbon12   ; fStatus(0) =   6; fStatus(1) = 11.177929140; break;
        case Carbon13    : fPartName = "Carbon13"   ; fPart = PhySt::Carbon13   ; fStatus(0) =   6; fStatus(1) = 12.112548247; break;
        case Carbon14    : fPartName = "Carbon14"   ; fPart = PhySt::Carbon14   ; fStatus(0) =   6; fStatus(1) = 13.043937223; break;
        case Nitrogen13  : fPartName = "Nitrogen13" ; fPart = PhySt::Nitrogen13 ; fStatus(0) =   7; fStatus(1) = 12.114768715; break;
        case Nitrogen14  : fPartName = "Nitrogen14" ; fPart = PhySt::Nitrogen14 ; fStatus(0) =   7; fStatus(1) = 13.043780747; break;
        case Nitrogen15  : fPartName = "Nitrogen15" ; fPart = PhySt::Nitrogen15 ; fStatus(0) =   7; fStatus(1) = 13.972512863; break;
        case Oxygen16    : fPartName = "Oxygen16"   ; fPart = PhySt::Oxygen16   ; fStatus(0) =   8; fStatus(1) = 14.899168518; break;
        case Oxygen17    : fPartName = "Oxygen17"   ; fPart = PhySt::Oxygen17   ; fStatus(0) =   8; fStatus(1) = 15.834590801; break;
        case Oxygen18    : fPartName = "Oxygen18"   ; fPart = PhySt::Oxygen18   ; fStatus(0) =   8; fStatus(1) = 16.766112187; break;
        default :
          std::cerr << Form("(ERROR) PhySt::setChrgMass() : It is not in PhySt's ParticleList.\n");
      }
    }

		void setPart(enum ParticleList part = PhySt::Proton) { setChrgMass(part); }

		void setBeta(Double_t beta, Double_t chrgSign = 0.) {
			if (IsMassLess()) return;
			
			Double_t sign = 0;
			if (NGNumc::EqualToZero(chrgSign) && !IsChrgLess()) 
				sign = NGNumc::Compare(Chrg());
			else 
				sign = ((IsChrgLess()) ? 1. : NGNumc::Compare(chrgSign));
			
			beta = std::fabs(beta);
			fStatus(2) = ((beta>=(1.-1e-6)) ? (1.-1.e-6) : ((beta<=0) ?0:beta));
			fStatus(4) = sign * ((NGNumc::EqualToZero(fStatus(2))) ? 0 : std::sqrt((1. / fStatus(2) / fStatus(2)) - 1.)); 
			fStatus(3) = std::fabs(fStatus(1) * fStatus(4)); 
		}

		void setMom(Double_t mom, Double_t chrgSign = 0.) {
			Double_t sign = 0;
			if (NGNumc::EqualToZero(chrgSign) && !IsChrgLess()) 
				sign = NGNumc::Compare(Chrg());
			else 
				sign = ((IsChrgLess()) ? 1. : NGNumc::Compare(chrgSign));
			
			fStatus(3) = std::fabs(mom);
			fStatus(4) = sign * (IsMassLess() ? fStatus(3) : (fStatus(3) / fStatus(1)));
			fStatus(2) = 1. / std::sqrt(1. + (1. / fStatus(4) / fStatus(4)));
		}
    
		void setEng(Double_t eng, Double_t chrgSign = 0.) { 
			eng = std::fabs(eng);
			if (eng < fStatus(1)) return;
			
			Double_t sign = 0;
			if (NGNumc::EqualToZero(chrgSign) && !IsChrgLess()) 
				sign = NGNumc::Compare(Chrg());
			else 
				sign = ((IsChrgLess()) ? 1. : NGNumc::Compare(chrgSign));
			
			fStatus(3) = std::sqrt(eng * eng	- fStatus(1) * fStatus(1));
			fStatus(4) = sign * (IsMassLess() ? fStatus(3) : (fStatus(3) / fStatus(1)));
			fStatus(2) = 1. / std::sqrt(1. + (1. / fStatus(4) / fStatus(4)));
		}
		
		void setKEng(Double_t keng, Double_t chrgSign = 0.) {
			Double_t eng = std::fabs(keng) + fStatus(1);
			setEng(eng, chrgSign);
		}

		void setEta(Double_t eta) {
			fStatus(4) = eta;
			fStatus(3) = (IsMassLess() ? std::fabs(fStatus(4)) : std::fabs(fStatus(1) * fStatus(4)));
			fStatus(2) = (IsMassLess() ? 1. : (1. / std::sqrt(1. + (1. / fStatus(4) / fStatus(4)))));
		}

		void setInvEta(Double_t ieta) {
			Double_t eta = 1. / ieta;
			if (!NGNumc::Valid(eta)) eta = 0.;
			setEta(eta);
		}

		void setRig(Double_t rig) {
			if (IsMassLess() || IsChrgLess()) return;
			Double_t eta = rig * std::fabs(fStatus(0) / fStatus(1));
			setEta(eta);
		}
   
		void setInvRig(Double_t invRig) {
			Double_t rig = 1. / invRig;
			if (!NGNumc::Valid(rig)) rig = 0.;
			setRig(rig);
		}

		void setMomVec(Double_t px, Double_t py, Double_t pz, Double_t eng, Double_t chrgSign = 0.) {
			eng = std::fabs(eng);
			fDir(0) = px; fDir(1) = py; fDir(2) = pz;
			fDir.Unit();
			setEng(eng, chrgSign);
		}
    
		void setMomVec(MtxLB::SVecD<4> & momVec, Double_t chrgSign = 0.) {
			setMomVec(momVec(0), momVec(1), momVec(2), momVec(3), chrgSign);
		}

	public :
		Bool_t vaild() {
			Bool_t vaild = (
				NGNumc::Valid(fPos(0)) && 
				NGNumc::Valid(fPos(1)) && 
				NGNumc::Valid(fPos(2)) && 
				NGNumc::Valid(fDir(0)) && 
				NGNumc::Valid(fDir(1)) && 
				NGNumc::Valid(fDir(2)) &&
				NGNumc::Valid(fStatus(0)) &&
				NGNumc::Valid(fStatus(1)) &&
				NGNumc::Valid(fStatus(2)) &&
				NGNumc::Valid(fStatus(3)) &&
				NGNumc::Valid(fStatus(4)) 
			);
			if (!vaild) init(this->fPart);
			return vaild;
		}

		void limit() {
			if (std::fabs(fPos(0)) > 1e4) fPos(0) = ((fPos(0) > 0) ? 1e4 : -1e4);
			if (std::fabs(fPos(1)) > 1e4) fPos(1) = ((fPos(1) > 0) ? 1e4 : -1e4);
			if (std::fabs(fPos(2)) > 1e4) fPos(2) = ((fPos(2) > 0) ? 1e4 : -1e4);
			if (!IsChrgLess()) {
				Double_t rigv = this->Rig();
				Double_t absr = std::fabs(rigv);
				Double_t sign = NGNumc::Compare(rigv);
				if (absr < 1e-6) { Double_t lmtv = 1e-6; this->setRig(sign * lmtv); }
				if (absr > 1e+9) { Double_t lmtv = 1e+9; this->setRig(sign * lmtv); }
			}
		}

		void print() {
			std::cout << Form("**** PhySt ****\n");
			std::cout << Form("Particle  : %s\n",           fPartName.Data());
			std::cout << Form("Charge    : %14.8f\n",       fStatus(0));
			std::cout << Form("Mass      : %14.8f [GeV]\n", fStatus(1));
			std::cout << Form("Beta      : %14.8f\n",       fStatus(2));
			std::cout << Form("Momentum  : %14.8f [GeV]\n", fStatus(3));
			std::cout << Form("Eta       : %14.8f\n",       fStatus(4));
			std::cout << Form("GammaBeta : %14.8f\n",       GammaBeta());
			std::cout << Form("Energy    : %14.8f [GeV]\n", Eng());
			std::cout << Form("Rigidity  : %14.8f [GV]\n",  Rig());
			std::cout << Form("Position  : (%14.8f %14.8f %14.8f) [cm]\n", fPos(0), fPos(1), fPos(2));
			std::cout << Form("Direction : (%14.8f %14.8f %14.8f)\n", fDir(0), fDir(1), fDir(2));
			std::cout << Form("\n");
		}

  public :
		TString PartName() { return fPartName; }
		ParticleList Part() { return fPart; }
		MtxLB::SVecD<5> & Status() { return fStatus; }
		Double_t Status(Short_t idx) { return fStatus(idx); }
		Double_t Chrg() { return fStatus(0); }
		Double_t Mass() { return fStatus(1); }
		Double_t Beta() { return fStatus(2); }
		Double_t Mom()  { return fStatus(3); }
		Double_t Eta()  { return fStatus(4); }

		Double_t ChrgMass() { return (IsMassLess() ? 0 : (fStatus(0) / fStatus(1))); }
		Double_t GammaBeta() { return std::fabs(fStatus(4)); }
		Double_t Gamma() { return (GammaBeta() / Beta()); }
		Double_t Eng() { return std::sqrt(fStatus(3) * fStatus(3) + fStatus(1) * fStatus(1)); }
		Double_t KEng() { return (Eng() - Mass()); }

		Double_t InvEta() { Double_t ieta = (1. / fStatus(4)); return (NGNumc::Valid(ieta) ? ieta : 0.); }
		Double_t Rig() { return (IsChrgLess() ? 0 : (fStatus(4) * std::fabs(fStatus(1) / fStatus(0)))); }
		Double_t InvRig() { return (IsChrgLess() ? 0 : (1. / Rig())); }

    MtxLB::SVecD<4> MomVec() { return MtxLB::SVecD<4>(fStatus(3) * fDir(0), fStatus(3) * fDir(1), fStatus(3) * fDir(2), Eng()); }
      
	public :
		Bool_t IsMassLess() { return (fPart == PhySt::Photon); }
		Bool_t IsChrgLess() { return (fPart == PhySt::Photon); }

	protected :
		TString           fPartName;
    enum ParticleList fPart; // kind of particle
		MtxLB::SVecD<5>   fStatus; // (charge, mass, beta, momentum, eta := (gammaBeta * chrgSign) or momentum)
	
	public :
		static void SortSpace(std::vector<PhySt> & phySt, Bool_t dirType = true) {
			if (phySt.size() < 2) return;
			if (dirType) std::sort(phySt.begin(), phySt.end(), PhyStSpace::SortFromMaxToMin());
			else         std::sort(phySt.begin(), phySt.end(), PhyStSpace::SortFromMinToMax());
		}
};


#endif // __MgntPhySt_H__
