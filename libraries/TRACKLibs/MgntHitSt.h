#ifndef __MgntHitSt_H__
#define __MgntHitSt_H__

#include <iostream>
#include <vector>
#include <map>
#include <cmath>

class HitSt {
  public :
    HitSt() { init(); }
    ~HitSt() {}

    void init() { 
			fId = -1; fSide = 0; fChrg = 1.; 
			fCoo = SVecD<3>(); fErr = HitSt::DEF_ERR;
			fDetCov(0, 0) = fErr(0) * fErr(0); fDetCov(1, 1) = fErr(1) * fErr(1);
		}

		void setId(Int_t id) { fId = id; }
		void setSide(Bool_t sx, Bool_t sy, Bool_t sz) { fSide = (sx * 1 + sy * 2 + sz * 4); }
		void setIdWithSide(Int_t id, Bool_t sx, Bool_t sy, Bool_t sz) { fId = id; fSide = (sx * 1 + sy * 2 + sz * 4); }
		void setChrg(Double_t chrg = 1.) { fChrg = chrg; }

    void setSpatialCoo(SVecD<3> & coo) { fCoo = coo; }
    void setSpatialCoo(Double_t x, Double_t y, Double_t z) { fCoo(0) = x; fCoo(1) = y; fCoo(2) = z; }

    void setSpatialErr(SVecD<3> & err) { 
			fErr = err; 
			fDetCov(0, 0) = fErr(0) * fErr(0); fDetCov(1, 1) = fErr(1) * fErr(1);
		}
    void setSpatialErr(Double_t ex, Double_t ey, Double_t ez) { 
			fErr(0) = ex; fErr(1) = ey; fErr(2) = ez; 
			fDetCov(0, 0) = fErr(0) * fErr(0); fDetCov(1, 1) = fErr(1) * fErr(1);
		}
    
		void setSpatial(SVecD<3> & coo, SVecD<3> & err) { 
			fCoo = coo; fErr = err; 
			fDetCov(0, 0) = fErr(0) * fErr(0); fDetCov(1, 1) = fErr(1) * fErr(1);
		}

    void setSpatial(Double_t px, Double_t py, Double_t pz, Double_t ex, Double_t ey, Double_t ez) {
      fCoo(0) = px; fCoo(1) = py; fCoo(2) = pz;
      fErr(0) = ex; fErr(1) = ey; fErr(2) = ez;
			fDetCov(0, 0) = fErr(0) * fErr(0); fDetCov(1, 1) = fErr(1) * fErr(1);
		}

    void setSpatial(Double_t state[6]) {
			fCoo(0) = state[0]; fCoo(1) = state[1]; fCoo(2) = state[2];
			fErr(3) = state[3]; fErr(4) = state[4]; fErr(5) = state[5];
			fDetCov(0, 0) = fErr(0) * fErr(0); fDetCov(1, 1) = fErr(1) * fErr(1);
		}

    void setSpatial(Double_t coo[3], Double_t err[3]) {
			fCoo(0) = coo[0]; fCoo(1) = coo[1]; fCoo(2) = coo[2];
			fErr(0) = err[0]; fErr(1) = err[1]; fErr(2) = err[2];
			fDetCov(0, 0) = fErr(0) * fErr(0); fDetCov(1, 1) = fErr(1) * fErr(1);
		}
		
		void print() {
			std::cout << Form("**** HitSt ****\n");
			std::cout << Form("Id    : %d\n", fId);
			std::cout << Form("Side  : [%d] -> (%d %d %d)\n", fSide, SideX(), SideY(), SideZ());
			std::cout << Form("Coo   : (%14.8f %14.8f %14.8f) [cm]\n", fCoo(0), fCoo(1), fCoo(2));
			std::cout << Form("Err   : (%14.8f %14.8f %14.8f) [cm]\n", fErr(0), fErr(1), fErr(2));
			std::cout << Form("\n");
		}

  public :
		Int_t   Id() { return fId; }
		Short_t Side() { return fSide; }
		Bool_t  SideX() { return ((fSide&1) == 1); }
		Bool_t  SideY() { return ((fSide&2) == 2); }
		Bool_t  SideZ() { return ((fSide&4) == 4); }

		inline SVecD<3> & Coo() { return fCoo; }
		inline Double_t & Coo(Short_t idx) { return fCoo(idx); }
		Double_t CX() { return fCoo(0); }
		Double_t CY() { return fCoo(1); }
		Double_t CZ() { return fCoo(2); }

		inline SVecD<3> & Err() { return fErr; }
		Double_t & Err(Short_t idx) { return fErr(idx); }
		Double_t EX() { return fErr(0); }
		Double_t EY() { return fErr(1); }
		Double_t EZ() { return fErr(2); }

		inline SMtxSymD<2> & DetCov() { return fDetCov; }
		inline Double_t & DetCov(Int_t i, Int_t j) { return fDetCov(i, j); }

  protected :
		Int_t              fId;     // hit Id
		Short_t            fSide;   // side := (x * 1 + y * 2 + z * 4)
		Double_t           fChrg;   // chrg
    SVecD<3>    fCoo;    // coord
		SVecD<3>    fErr;    // error
		SMtxSymD<2> fDetCov; // error matrix
		
		static const SVecD<3> DEF_ERR; // default error

	public :
		static void SortHits(std::vector<HitSt> & hits, Bool_t dirType = true) {
			if (hits.size() < 2) return;
			if (dirType) std::sort(hits.begin(), hits.end(), HitSt::SortFromMaxToMin());
			else         std::sort(hits.begin(), hits.end(), HitSt::SortFromMinToMax());
		}
	
	protected :
		struct SortFromMinToMax {
			Bool_t operator() (const HitSt & hit1, const HitSt & hit2) {
				if (hit1.fCoo(2) < hit2.fCoo(2)) return true;
				else return false;
			}
			Bool_t operator() (HitSt * hit1, HitSt * hit2) {
				if (hit1->fCoo(2) < hit2->fCoo(2)) return true;
				else return false;
			}
		};

		struct SortFromMaxToMin {
			Bool_t operator() (const HitSt & hit1, const HitSt & hit2) {
				if (hit1.fCoo(2) > hit2.fCoo(2)) return true;
				else return false;
			}
			Bool_t operator() (HitSt * hit1, HitSt * hit2) {
				if (hit1->fCoo(2) > hit2->fCoo(2)) return true;
				else return false;
			}
		};
};

const SVecD<3> HitSt::DEF_ERR(24.0e-4, 10.0e-4, 300.0e-4);

#endif // __MgntHitSt_H__
