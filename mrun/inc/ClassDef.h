/**********************************************************
 * Author        : Hsin-Yi Chou
 * Email         : hchou@cern.ch
 * Last modified : 2015-07-22 16:29
 * Filename      : ClassDef.h
 * Description   :
 * *******************************************************/
#ifndef __ClassDef_H__
#define __ClassDef_H__

#include <iostream>
#include <vector>
#include <TObject.h>
#include <TString.h>


// LIST
class LIST : public TObject {
	public :
		LIST() { init(); }
		~LIST() {}

		void init() {
            file   = "";
			run    = 0;
			event  = 0;
            entry  = 0;
            utime  = 0;
		}

	public :
        TString file;
		UInt_t  run;
		UInt_t  event;
        UInt_t  entry;
        UInt_t  utime;

	ClassDef(LIST, 1)
};

// RTI
class RTI : public TObject {
	public :
		RTI() { init(); }
		~RTI() {}

		void init() {
			flag = true;
            good = true;
			zenith = 0;
            livetime = 0;

            std::fill_n(GTOD, 3, 0);
            std::fill_n(GM, 2, 0);
            std::fill_n(GAT, 2, 0);

            std::fill_n(Stoermer[0], 4*2, 0);
            std::fill_n(IGRF[0], 4*2, 0);
            
            min_Stoermer = 0;
            max_Stoermer = 0;
            min_IGRF = 0;
            max_IGRF = 0;
            
            std::fill_n(tk_align[0], 2*2, 0);
            tk_temp = 0;
            
            is_in_SAA = false;
		}

	public :
        Bool_t  flag;
        Bool_t  good;
		Float_t zenith;         // ams zenith (z-axis) angle [degrees]
		Float_t livetime;       // fraction of "non-busy" time
		
        Float_t GTOD[3];        // earth coordinate [m, rad] (radius, theta, phi)
		Float_t GM[2];          // geomagnetic [rad] (latitude, longitude)
		Float_t GAT[2];         // ams pointing galatic [rad] (latitude, longitude)
        
        Float_t Stoermer[4][2]; // Stoermer cutoff [GV]
        Float_t IGRF[4][2];     // IGRF cutoff [GV]

        Float_t min_Stoermer;   // Stoermer min cutoff [GV]
        Float_t max_Stoermer;   // Stoermer max cutoff [GV]
        Float_t min_IGRF;       // IGRF min cutoff [GV]
        Float_t max_IGRF;       // IGRF max cutoff [GV]

		Float_t tk_align[2][2]; // L1 x,y L9 x,y
		Float_t tk_temp;        // tracker temperature (Sensor A)
        
        Bool_t  is_in_SAA;      // true, if ams in south atlantic anomaly

	ClassDef(RTI, 1)
};


#endif // __ClassDef_H__
