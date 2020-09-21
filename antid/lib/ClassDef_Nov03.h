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
            weight = 1;
		}

	public :
        TString file;
		UInt_t  run;
		UInt_t  event;
        UInt_t  entry;
        UInt_t  utime;
        Float_t weight;

	ClassDef(LIST, 1)
};


// G4MC
class G4MC : public TObject {
	public :
		G4MC() { init(); }
		~G4MC() {}

		void init() {
            prm_chrg = 0;
            prm_mass = 0;
            prm_mom = 0;
            std::fill_n(prm_loc, 3, 0);
            std::fill_n(prm_dir, 3, 0);

            std::fill_n(tk, 9, false);
            std::fill_n(tk_mom, 9, 0);
            std::fill_n(tk_mes[0], 9*3, 0);
            std::fill_n(tk_loc[0], 9*3, 0);
            std::fill_n(tk_dir[0], 9*3, 0);
            std::fill_n(tk_edep, 9, 0);
            
            std::fill_n(tf, 4, false);
            std::fill_n(tf_beta, 4, 0);
            std::fill_n(tf_time, 4, 0);
            std::fill_n(tf_loc[0], 4*3, 0);
            
            std::fill_n(td, 20, false);
            std::fill_n(td_mom, 20, 0);
            std::fill_n(td_loc[0], 20*3, 0);
		}

	public :
        Short_t prm_chrg;
        Float_t prm_mass;
        Float_t prm_mom;
        Float_t prm_loc[3];
        Float_t prm_dir[3];

        Bool_t  tk[9];
        Float_t tk_mom[9];
        Float_t tk_mes[9][3];
        Float_t tk_loc[9][3];
        Float_t tk_dir[9][3];
        Float_t tk_edep[9];
        
        Bool_t  tf[4];
        Float_t tf_beta[4];
        Float_t tf_time[4];
        Float_t tf_loc[4][3];
        
        Bool_t  td[20];
        Float_t td_mom[20];
        Float_t td_loc[20][3];

	ClassDef(G4MC, 1)
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
            live_time = 0;

            std::fill_n(GTOD, 3, 0);
            std::fill_n(GM, 2, 0);
            std::fill_n(GAT, 2, 0);

            max_Stormer = 0;
            max_IGRF = 0;
            
            std::fill_n(tk_align[0], 2*2, 0);
            tk_temp = 0;
            
            is_in_SAA = false;
            is_in_shadow = false;
		}

	public :
        Bool_t  flag;
        Bool_t  good;
		Float_t zenith;         // ams zenith (z-axis) angle [degrees]
		Float_t live_time;      // fraction of "non-busy" time
		
        Float_t GTOD[3];        // earth coordinate [m, rad] (radius, theta, phi)
		Float_t GM[2];          // geomagnetic [rad] (latitude, longitude)
		Float_t GAT[2];         // ams pointing galatic [rad] (latitude, longitude)
        
        Float_t max_Stormer;    // Stormer max cutoff [GV]
        Float_t max_IGRF;       // IGRF max cutoff [GV]

		Float_t tk_align[2][2]; // L1 x,y L9 x,y
		Float_t tk_temp;        // tracker temperature (Sensor A)
        
        Bool_t  is_in_SAA;      // true, if ams in south atlantic anomaly

		Bool_t is_in_shadow;    // particle pass through the ISS Solar Array
		                        // return
								//   0, not in shadow
								//   1, in shadow

	ClassDef(RTI, 1)
};


// TRG
class TRG : public TObject {
	public :
		TRG() { init(); }
		~TRG() {}

		void init() {
			bit = 0;
			logic = 0;
			physics = 0;
		}

	public :
		Short_t bit;     // level 1 trigger : bit := extTrg * 1 + unBiasTrgTOF * 2 + unBiasTrgECAL * 4 + physTrg * 8;
		Short_t logic;   // level 1 trigger : logic
		Short_t physics; // level 1 trigger : physics

	ClassDef(TRG, 1)
};


// ACC
class ACC : public TObject {
	public :
		ACC() { init(); }
		~ACC() {}

		void init() {
			num_cls = 0;
		}

	public :
		Short_t num_cls;

	ClassDef(ACC, 1)
};


// TOF
class TOF : public TObject {
	public :
		TOF() { init(); }
		~TOF() {}

		void init() {
            num_cls = 0;
            num_beta = 0;
            
            status = 0;
            bit = 0;
            patt = 0;
            beta = 0;
   
            mc_beta = 0;

            nchi_t = 0;
            nchi_c = 0;
            
            Qall = 0;
            Zall = 0;

            std::fill_n(lay, 4, false);
            std::fill_n(loc[0], 4*3, 0);
            std::fill_n(T, 4, 0);
            std::fill_n(Q, 4, 0);

            std::fill_n(num_extcls, 4, 0);
            extcls_noise = 0;

            num_in_time_cls = 0;
		}

	public :
        Short_t num_cls;
        Short_t num_beta;

        Bool_t  status;
        Short_t bit;
        Short_t patt;
        Float_t beta;
        
        Float_t mc_beta;

        Float_t nchi_t;
        Float_t nchi_c;

        Float_t Qall;
        Short_t Zall;

        Bool_t  lay[4];
        Float_t loc[4][3];
        Float_t T[4];
        Float_t Q[4];
		
        // extern clusters
		Short_t num_extcls[4];
        Short_t extcls_noise; // noiseALL*1 + noiseUTOF*2 + noiseLTOF*4
		
        Short_t num_in_time_cls;

	ClassDef(TOF, 1)
};


// TRK
class TRK : public TObject {
	public :
		TRK() { init(); }
		~TRK() {}

		void init() {
            patt = 0;
			pattXY = 0;
			
            QInMin = 0;
            QIn = 0;
			QL2 = 0;
			QL1 = 0;
			QL9 = 0;
            
            std::fill_n(lay, 9, 0);
            std::fill_n(strip[0], 9*2, 0);
            std::fill_n(eta[0], 9*2, 0);
            std::fill_n(loc[0], 9*3, 0);
            std::fill_n(chrg_hl[0], 9*3, 0);
            std::fill_n(chrg_yj[0], 9*3, 0);

            std::fill_n(ext_num_hit, 9, 0);
            std::fill_n(ext_chrg_hl, 9, 0);
            std::fill_n(ext_chrg_yj, 9, 0);

            std::fill_n(ck_status, 4, false);
            std::fill_n(ck_ndof[0], 4*2, 0);
            std::fill_n(ck_nchi[0], 4*2, 0);
            std::fill_n(ck_rig, 4, 0);
            std::fill_n(ck_cen_state, 6, 0);
            std::fill_n(ck_top_state, 6, 0);
            
            std::fill_n(ck_cpu_time, 4, 0);
            
            std::fill_n(kf_status, 4, false);
            std::fill_n(kf_ndof[0], 4*2, 0);
            std::fill_n(kf_nchi[0], 4*2, 0);
            std::fill_n(kf_cen_rig, 4, 0);
            std::fill_n(kf_top_rig, 4, 0);
            std::fill_n(kf_cen_state, 6, 0);
            std::fill_n(kf_top_state, 6, 0);
            
            std::fill_n(kf_cpu_time, 4, 0);
		}

	public :
        Short_t num_track;

		// bitPatt   := hasInner   * 1 + hasL2   * 2 + hasL1   * 4 + hasL9   * 8
		// bitPattXY := hasInnerXY * 1 + hasL2XY * 2 + hasL1XY * 4 + hasL9XY * 8
		// Inner    (XY) := ((bitPattXY&  1)==  1)
		// InnerL1  (XY) := ((bitPattXY&  5)==  5)
		// InnerL9  (XY) := ((bitPattXY&  9)==  9)
		// FullSpan (XY) := ((bitPattXY& 13)== 13)
		Short_t patt;
		Short_t pattXY;

		// Track Charge
		Float_t QInMin;
		Float_t QIn;
		Float_t QL1;
		Float_t QL2;
		Float_t QL9;
        
        // Track Hits
        Short_t  lay[9]; // 0:no, 1:x, 2:y, 3:xy
        Short_t  strip[9][2];
        Float_t  eta[9][2];
        Double_t loc[9][3];
        Float_t  chrg_hl[9][3];
        Float_t  chrg_yj[9][3];
        
        // Track External Hits
        Short_t ext_num_hit[9];
        Float_t ext_chrg_hl[9];
        Float_t ext_chrg_yj[9];

        // Choutko [Inn InnL1 InnL9 FS]
        Bool_t  ck_status[4];
        Short_t ck_ndof[4][2];
        Float_t ck_nchi[4][2];
        Float_t ck_rig[4];
        Float_t ck_cen_state[6]; // [Inn]
        Float_t ck_top_state[6]; // [Inn]
        
        Float_t ck_cpu_time[4];

        Bool_t  kf_status[4];
        Short_t kf_ndof[4][2];
        Float_t kf_nchi[4][2];
        Float_t kf_cen_rig[4];
        Float_t kf_top_rig[4];
        Float_t kf_cen_state[6]; // [Inn]
        Float_t kf_top_state[6]; // [Inn]
        
        Float_t kf_cpu_time[4];

	ClassDef(TRK, 1)
};

// TRD
class TRD : public TObject {
	public :
		TRD() { init(); }
		~TRD() {}

		void init() {
            num_cls = 0;
            
            num_track  = 0;
            num_hit    = 0;

            num_Htrack = 0;
            num_Hhit   = 0;

            status = false;
            std::fill_n(state, 6, 0);

            num_extra_seg = 0;
            num_extra_hit = 0;
            std::fill_n(num_vtx, 2, 0);

            tdLLR_status = false;
            tdLLR_num_hit = 0;
            tdLLR_ep = 0;
            tdLLR_eh = 0;
            tdLLR_ph = 0;
            tdLL_el = 0;
            tdLL_pr = 0;
            tdLL_he = 0;
            
            num_tdHit = 0;
            tdHit_lay.clear();
            tdHit_len.clear();
            tdHit_amp.clear();
            tdHit_lz.clear();
            
            tkLLR_status = false;
            tkLLR_num_hit = 0;
            tkLLR_ep = 0;
            tkLLR_eh = 0;
            tkLLR_ph = 0;
            tkLL_el = 0;
            tkLL_pr = 0;
            tkLL_he = 0;
            
            num_tkHit = 0;
            tkHit_lay.clear();
            tkHit_len.clear();
            tkHit_amp.clear();
            tkHit_lz.clear();
		}

	public :
        Short_t num_cls;
        
        Short_t num_track;
        Short_t num_hit;

        Short_t num_Htrack;
        Short_t num_Hhit;
        
        // TrdTrack or TrdHTrack (first TrdH, second Trd)
		Bool_t  status;
		Float_t state[6]; // coo, dir
        
        Short_t num_extra_seg;
        Short_t num_extra_hit;
        Short_t num_vtx[2];

        Bool_t  tdLLR_status;
        Short_t tdLLR_num_hit;
        Float_t tdLLR_ep;
        Float_t tdLLR_eh;
        Float_t tdLLR_ph;
        Float_t tdLL_el;
        Float_t tdLL_pr;
        Float_t tdLL_he;
        
        Short_t num_tdHit;
        std::vector<Short_t> tdHit_lay;
        std::vector<Float_t> tdHit_len;
        std::vector<Float_t> tdHit_amp;
        std::vector<Float_t> tdHit_lz;

        Bool_t  tkLLR_status;
        Short_t tkLLR_num_hit;
        Float_t tkLLR_ep;
        Float_t tkLLR_eh;
        Float_t tkLLR_ph;
        Float_t tkLL_el;
        Float_t tkLL_pr;
        Float_t tkLL_he;
        
        Short_t num_tkHit;
        std::vector<Short_t> tkHit_lay;
        std::vector<Float_t> tkHit_len;
        std::vector<Float_t> tkHit_amp;
        std::vector<Float_t> tkHit_lz;


	ClassDef(TRD, 1)
};


// ECAL
class ECAL : public TObject {
	public :
		ECAL() { init(); }
		~ECAL() {}

		void init() {
            num_shower = 0;

            status = false;
            engD = 0;
            engE = 0;
            mvaBDT = 0;
            chrg = 0;

            hadron_apex = 0;
            hadron_eng = 0;

            num_hit = 0;
            hit_idx.clear();
            hit_edep.clear();
		}

	public :
        Short_t num_shower;

        Bool_t  status;
        Float_t engD;
        Float_t engE;
        Float_t mvaBDT;
        Float_t chrg;

        Short_t hadron_apex;
        Float_t hadron_eng;

        Short_t num_hit;
        std::vector<Short_t> hit_idx;  // plane * NCell + cell
        std::vector<Float_t> hit_edep;

	ClassDef(ECAL, 1)
};


// RICH
class RICH : public TObject {
	public :
		RICH() { init(); }
		~RICH() {}

		void init() {
            num_ring = 0;

            status = false;
            kind = 0;
            tile = 0;
            refz = 0;
            dist = 0;
            beta = 0;
            chrg = 0;

            isGood = false;
            nhit = 0;
            npmt = 0;
            prob = 0;
            cstcb = 0;
            cstcq = 0;

            num_clsPE = 0;
            num_expPE = 0;
            num_colPE = 0;
            eft_colPE = 0;
            num_clsZ1 = 0;

            self_status = false;
            self_kind = 0;
            self_tile = 0;
            self_index = 0;
            self_dist = 0;
            self_loc_lx = 0;
            self_loc_ly = 0;
            self_loc_tha = 0;
            self_loc_phi = 0;
            self_beta_crr = 1;
            self_is_good_geom = false;
            self_is_bad_tile = false;

            std::fill_n(self_rad_loc, 3, 0);
            std::fill_n(self_rad_dir, 3, 0);
            std::fill_n(self_pmt_loc, 3, 0);
            
            self_expnpe = 0;
            self_border = 0;
            self_trace  = 0;
            
            self_num_stone = 0;
            self_num_cloud = 0;
            self_num_tumor = 0;
            self_num_ghost = 0;
            
            self_nhit_total = 0;
            self_nhit_stone = 0;
            self_nhit_cloud = 0;
            self_nhit_tumor = 0;
            self_nhit_ghost = 0;
            self_nhit_other = 0;
            self_nhit_other_inn = 0;
            self_nhit_other_out = 0;
            
            self_npe_total = 0.0;
            self_npe_stone = 0.0;
            self_npe_cloud = 0.0;
            self_npe_tumor = 0.0;
            self_npe_ghost = 0.0;
            self_npe_other_inn = 0.0;
            self_npe_other_out = 0.0;
            
            self_stn_status = false;
            self_stn_nhit   = 0;
            self_stn_npmt   = 0;
            self_stn_lx     = 0;
            self_stn_ly     = 0;
            self_stn_npe    = 0;
            self_stn_dist   = 0;
            self_stn_nchi   = 0;
            self_stn_chic   = 0;
            
            self_cld_status   = false;
            self_cld_nhit     = 0;
            self_cld_npmt     = 0;
            self_cld_beta     = 0;
            self_cld_cbta     = 0;
            self_cld_npe      = 0;
            self_cld_nchi     = 0;
            self_cld_misjudge = 0;
            self_cld_expnpe   = 0;
            self_cld_border   = 0;
            self_cld_trace    = 0;
            self_cld_accuracy = 0;
            
            self_cldhit_chann.clear();
            self_cldhit_beta.clear();
            self_cldhit_phi.clear();
            self_cldhit_npe.clear();
            self_cldhit_lx.clear();
            self_cldhit_ly.clear();
            
            num_hit = 0;
            hit_chann.clear();
            hit_type.clear();
            hit_dbeta.clear();
            hit_rbetaA.clear();
            hit_rbetaB.clear();
            hit_npe.clear();
            hit_lx.clear();
            hit_ly.clear();
		}

	public :
        Short_t num_ring;

        // Official RichRingR
		Bool_t  status;
        Short_t kind;
        Short_t tile;
        Float_t refz;
        Float_t dist;
		Float_t beta;
		Float_t chrg;
        
		Bool_t  isGood;
        Short_t nhit;
        Short_t npmt;
        Float_t prob;
        Float_t cstcb; // consistency beta
        Float_t cstcq; // consistency charge
        Float_t num_clsPE;
        Float_t num_expPE;
        Float_t num_colPE;
        Float_t eft_colPE;
        Short_t num_clsZ1;

        // Self Ring
        Bool_t  self_status;
        Short_t self_kind;
        Short_t self_tile;
        Float_t self_index;
        Float_t self_dist;
        Float_t self_loc_lx;
        Float_t self_loc_ly;
        Float_t self_loc_tha;
        Float_t self_loc_phi;
        Float_t self_beta_crr;
        Bool_t  self_is_good_geom;
        Bool_t  self_is_bad_tile;

        Float_t self_rad_loc[3];
        Float_t self_rad_dir[3];
        Float_t self_pmt_loc[3];
        
        Float_t self_expnpe;
        Float_t self_border;
        Float_t self_trace;

        Short_t self_num_stone;
        Short_t self_num_cloud;
        Short_t self_num_tumor;
        Short_t self_num_ghost;

        Short_t self_nhit_total;
        Short_t self_nhit_stone;
        Short_t self_nhit_cloud;
        Short_t self_nhit_tumor;
        Short_t self_nhit_ghost;
        Short_t self_nhit_other;
        Short_t self_nhit_other_inn;
        Short_t self_nhit_other_out;

        Float_t self_npe_total;
        Float_t self_npe_stone;
        Float_t self_npe_cloud;
        Float_t self_npe_tumor;
        Float_t self_npe_ghost;
        Float_t self_npe_other;
        Float_t self_npe_other_inn;
        Float_t self_npe_other_out;

        Bool_t  self_stn_status;
        Short_t self_stn_nhit;
        Short_t self_stn_npmt;
        Float_t self_stn_lx;
        Float_t self_stn_ly;
        Float_t self_stn_npe;
        Float_t self_stn_dist;
        Float_t self_stn_nchi;
        Float_t self_stn_chic;

        Bool_t  self_cld_status;
        Short_t self_cld_nhit;
        Short_t self_cld_npmt;
        Float_t self_cld_beta;
        Float_t self_cld_cbta;
        Float_t self_cld_npe;
        Float_t self_cld_nchi;
        Float_t self_cld_misjudge;
        Float_t self_cld_expnpe;
        Float_t self_cld_border;
        Float_t self_cld_trace;
        Float_t self_cld_accuracy;

        std::vector<Short_t>  self_cldhit_chann;
        std::vector<Double_t> self_cldhit_beta;
        std::vector<Float_t>  self_cldhit_phi;
        std::vector<Float_t>  self_cldhit_npe;
        std::vector<Float_t>  self_cldhit_lx;
        std::vector<Float_t>  self_cldhit_ly;

        // Hit Information
        Short_t num_hit;
        std::vector<Short_t> hit_chann;
        std::vector<Short_t> hit_type;
        std::vector<Float_t> hit_dbeta;
        std::vector<Float_t> hit_rbetaA;
        std::vector<Float_t> hit_rbetaB;
        std::vector<Float_t> hit_npe;
        std::vector<Float_t> hit_lx;
        std::vector<Float_t> hit_ly;

    ClassDef(RICH, 1)
};


// HYC
class HYC : public TObject {
	public :
		HYC() { init(); }
		~HYC() {}

		void init() {
            chrg = 0;
            mass = 0;

            std::fill_n(geom_status, 4, false);
            std::fill_n(geom_ndof_x, 4, 0);
            std::fill_n(geom_ndof_y, 4, 0);
            std::fill_n(geom_nchi_x, 4, 0);
            std::fill_n(geom_nchi_y, 4, 0);
            
            std::fill_n(geom_cen_loc[0], 4*3, 0);
            std::fill_n(geom_cen_dir[0], 4*3, 0);
            std::fill_n(geom_cen_rig, 4, 0);
            
            std::fill_n(geom_top_loc[0], 4*3, 0);
            std::fill_n(geom_top_dir[0], 4*3, 0);
            std::fill_n(geom_top_rig, 4, 0);
            
            std::fill_n(geom_cpu_time, 4, 0);
            
            std::fill_n(vel_status, 3, false);
            std::fill_n(vel_ndof, 3, 0);
            std::fill_n(vel_nchi, 3, 0);
            
            std::fill_n(vel_cen_loc[0], 3*3, 0);
            std::fill_n(vel_cen_dir[0], 3*3, 0);
            std::fill_n(vel_cen_bta, 3, 0);
            
            std::fill_n(vel_top_loc[0], 3*3, 0);
            std::fill_n(vel_top_dir[0], 3*3, 0);
            std::fill_n(vel_top_bta, 3, 0);
            
            std::fill_n(vel_cpu_time, 3, 0);
            
            std::fill_n(phys_status[0], 2*3, false);
            std::fill_n(phys_ndof_x[0], 2*3, 0);
            std::fill_n(phys_ndof_y[0], 2*3, 0);
            std::fill_n(phys_ndof_b[0], 2*3, 0);
            std::fill_n(phys_nchi_x[0], 2*3, 0);
            std::fill_n(phys_nchi_y[0], 2*3, 0);
            std::fill_n(phys_nchi_b[0], 2*3, 0);
            
            std::fill_n(phys_cen_loc[0][0], 2*3*3, 0);
            std::fill_n(phys_cen_dir[0][0], 2*3*3, 0);
            std::fill_n(phys_cen_rig[0], 2*3, 0);
            std::fill_n(phys_cen_bta[0], 2*3, 0);
            
            std::fill_n(phys_top_loc[0][0], 2*3*3, 0);
            std::fill_n(phys_top_dir[0][0], 2*3*3, 0);
            std::fill_n(phys_top_rig[0], 2*3, 0);
            std::fill_n(phys_top_bta[0], 2*3, 0);
            
            std::fill_n(phys_cpu_time[0], 2*3, 0);
            
            std::fill_n(mutr_status, 3, false);
            std::fill_n(mutr_ndof_x, 3, 0);
            std::fill_n(mutr_ndof_y, 3, 0);
            std::fill_n(mutr_ndof_b, 3, 0);
            std::fill_n(mutr_nchi_x, 3, 0);
            std::fill_n(mutr_nchi_y, 3, 0);
            std::fill_n(mutr_nchi_b, 3, 0);
            
            std::fill_n(mutr_mass, 3, 0);
            std::fill_n(mutr_sqrm, 3, 0);
            
            std::fill_n(mutr_cen_loc[0], 3*3, 0);
            std::fill_n(mutr_cen_dir[0], 3*3, 0);
            std::fill_n(mutr_cen_rig, 3, 0);
            std::fill_n(mutr_cen_bta, 3, 0);
            
            std::fill_n(mutr_top_loc[0], 3*3, 0);
            std::fill_n(mutr_top_dir[0], 3*3, 0);
            std::fill_n(mutr_top_rig, 3, 0);
            std::fill_n(mutr_top_bta, 3, 0);
            
            std::fill_n(mutr_cpu_time, 3, 0);
        }

    public :
        Int_t   chrg;
        Float_t mass;

        // Geometry Fitting
        // inn, innL1, innL9, fullspan, innU, innL, innM
        // state (lx ly lz ux uy uz rig)
        Bool_t  geom_status[4];
        Short_t geom_ndof_x[4];
        Short_t geom_ndof_y[4];
        Float_t geom_nchi_x[4];
        Float_t geom_nchi_y[4];
        
        Float_t geom_cen_loc[4][3];
        Float_t geom_cen_dir[4][3];
        Float_t geom_cen_rig[4];

        Float_t geom_top_loc[4][3];
        Float_t geom_top_dir[4][3];
        Float_t geom_top_rig[4];

        Float_t geom_cpu_time[4];

        // Velocity Fitting
        // time, time-chrg, cherenkov
        Bool_t  vel_status[3];
        Short_t vel_ndof[3];
        Float_t vel_nchi[3];
        
        Float_t vel_cen_loc[3][3];
        Float_t vel_cen_dir[3][3];
        Float_t vel_cen_bta[3];
        
        Float_t vel_top_loc[3][3];
        Float_t vel_top_dir[3][3];
        Float_t vel_top_bta[3];
        
        Float_t vel_cpu_time[3];
        
        // Phys Fitting
        // particle type (prm, sec)
        // inn-time-chrg, inn-cherenkov
        // state (lx ly lz ux uy uz rig beta)
        Bool_t  phys_status[2][3];
        Short_t phys_ndof_x[2][3];
        Short_t phys_ndof_y[2][3];
        Short_t phys_ndof_b[2][3];
        Float_t phys_nchi_x[2][3];
        Float_t phys_nchi_y[2][3];
        Float_t phys_nchi_b[2][3];

        Float_t phys_cen_loc[2][3][3];
        Float_t phys_cen_dir[2][3][3];
        Float_t phys_cen_rig[2][3];
        Float_t phys_cen_bta[2][3];

        Float_t phys_top_loc[2][3][3];
        Float_t phys_top_dir[2][3][3];
        Float_t phys_top_rig[2][3];
        Float_t phys_top_bta[2][3];
        
        Float_t phys_cpu_time[2][3];

        // Mass Fitting
        // inn-time-chrg, inn-cherenkov
        // state (lx ly lz ux uy uz rig beta)
        Bool_t  mutr_status[3];
        Short_t mutr_ndof_x[3];
        Short_t mutr_ndof_y[3];
        Short_t mutr_ndof_b[3];
        Float_t mutr_nchi_x[3];
        Float_t mutr_nchi_y[3];
        Float_t mutr_nchi_b[3];
        
        Float_t mutr_sqrm[3];
        Float_t mutr_mass[3];
        
        Float_t mutr_cen_loc[3][3];
        Float_t mutr_cen_dir[3][3];
        Float_t mutr_cen_rig[3];
        Float_t mutr_cen_bta[3];
        
        Float_t mutr_top_loc[3][3];
        Float_t mutr_top_dir[3][3];
        Float_t mutr_top_rig[3];
        Float_t mutr_top_bta[3];
        
        Float_t mutr_cpu_time[3];

    ClassDef(HYC, 1)
};


#endif // __ClassDef_H__
