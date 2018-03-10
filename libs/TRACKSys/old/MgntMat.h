#ifndef __MgntMat_H__
#define __MgntMat_H__

#include <iostream>

#include <sys/mman.h>  
#include <unistd.h>  
#include <stdio.h>  
#include <fcntl.h>  
#include <sys/stat.h>  
#include <stdlib.h>  
#include <string.h>  
#include <errno.h> 
#include <sys/types.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TDirectory.h"

#include "MgntPhySt.h"


//---- Material ---//
class Material {
    public :
        Material() { init(); }
        ~Material() {}

        void init() {
            fVac = true;
            std::fill_n(fElm, Material::NumOfElm, false);
            std::fill_n(fDen, Material::NumOfElm, 0.);
        }

        void print() {
            std::cout << Form("**** Material ****\n");
            std::cout << Form("Vacuum             : %d\n", fVac);
            std::cout << Form("Element            : %7s %7s %7s %7s %7s %7s %7s %7s %7s\n", "H", "C", "Ne", "O", "F", "Na", "Al", "Si", "Pb");
            std::cout << Form("                   : %7d %7d %7d %7d %7d %7d %7d %7d %7d\n", fElm[0], fElm[1], fElm[2], fElm[3], fElm[4], fElm[5], fElm[6], fElm[7], fElm[8]);
            std::cout << Form("        [mol/cm^3] : %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n", fDen[0], fDen[1], fDen[2], fDen[3], fDen[4], fDen[5], fDen[6], fDen[7], fDen[8]);
            std::cout << Form("\n");
        }

    public :
        // 1 [H]     2 [C]   3 [N]     4 [O]   5 [F]     6 [Na]  7 [Al]    8 [Si]   9 [Pb]
        // Hydrogen, Carbon, Nitrogen, Oxygen, Fluorine, Sodium, Aluminum, Silicon, Lead
        static constexpr Short_t  NumOfElm = 9;
        static constexpr Double_t AtomChrg[Material::NumOfElm] = // [1]
        { 1, 6, 7, 8, 9, 11, 13, 14, 82 };
        static constexpr Double_t AtomMass[Material::NumOfElm] = // [g mol^-1]
        { 1.007947, 12.01078, 14.00672, 15.99943, 18.99840325, 22.989769282, 26.98153868, 28.08553, 207.21};

        Bool_t   fVac;
        Bool_t   fElm[Material::NumOfElm];
        Double_t fDen[Material::NumOfElm]; // [mol cm^-3]
};

//---- MgntMat ----//
class MgntMat {
    public:
        MgntMat() {}
        ~MgntMat() {}

        static Bool_t Load();

        static Material & GetMat() { return CurMat; }
        static Material & GetMat(SVecD<3> & coo, Bool_t isUsedTree = false) { return GetMat(coo(0), coo(1), coo(2), isUsedTree); } 
        static Material & GetMat(Double_t coo[3], Bool_t isUsedTree = false) { return GetMat(coo[0], coo[1], coo[2], isUsedTree); } 

        static Material & GetMat(SVecD<3> & sat, SVecD<3> & end, Bool_t isUsedTree = false, Bool_t isFastScan = true) { return GetMat(sat(0), sat(1), sat(2), end(0), end(1), end(2), isUsedTree, isFastScan); }
        static Material & GetMat(Double_t sat[3], Double_t end[3], Bool_t isUsedTree = false, Bool_t isFastScan = true) { return GetMat(sat[0], sat[1], sat[2], end[0], end[1], end[2], isUsedTree, isFastScan); }

    protected :
        static Bool_t  LoadStatus;
        static TFile * LoadFileROOT;
        static TTree * TreeElm;
        static Bool_t  Vac;
        static Float_t Coo[3];
        static Bool_t  Elm[Material::NumOfElm];
        static Float_t Den[Material::NumOfElm]; // [mol cm^-3]

        static constexpr Long64_t BinXSize = 560;
        static constexpr Double_t  BinXWidth = 0.5; // [cm]
        static constexpr Double_t  BinXRangeL = -140.0;
        static constexpr Double_t  BinXRangeU =  140.0;

        static constexpr Long64_t BinYSize = 780;
        static constexpr Double_t  BinYWidth = 0.5; // [cm]
        static constexpr Double_t  BinYRangeL = -195.0;
        static constexpr Double_t  BinYRangeU =  195.0;

        static constexpr Long64_t BinZSize = 780;
        static constexpr Double_t  BinZWidth = 0.5; // [cm]
        static constexpr Double_t  BinZRangeL = -195.0;
        static constexpr Double_t  BinZRangeU =  195.0;

        static constexpr Double_t  BinFastScanStep = 0.35; // [cm]
        static constexpr Double_t  BinScanStep = 0.15; // [cm]

        static Long64_t GetEntry(Double_t x, Double_t y, Double_t z, Bool_t isUsedTree = false);
        static Material & GetMat(Long64_t entry, Bool_t isUsedTree = false);
        static Material & GetMat(Double_t x, Double_t y, Double_t z, Bool_t isUsedTree = false); 
        static Material & GetMat(Double_t satx, Double_t saty, Double_t satz, Double_t endx, Double_t endy, Double_t endz, Bool_t isUsedTree = false, Bool_t isFastScan = true); 

        //---- MaterialBox ----//
        struct MaterialBox {
            Bool_t  fVac;
            Float_t fCoo[3];
            Bool_t  fElm[9];
            Float_t fDen[9];
        };

        static       void *        LoadFileBIN;
        static       MaterialBox * MatBox;
        static const size_t        MatBoxSize;
        static const Long64_t      MatDetBIdx[7];
        static const Long64_t      MatDetSize[7];
        static const Double_t      MatBinStep[3];
        static const Long64_t      MatBoundIndex[12]; 
        static const Double_t      MatBoundRange[13];
        static const Double_t      MatDetRange[7][3];
        static const Long64_t      MatDetBSize[7][3];


        static Material MergeMat;
        static Material NullMat;
        static Material CurMat;
        static Long64_t CurEntryID;
};

Bool_t  MgntMat::LoadStatus = false;
TFile * MgntMat::LoadFileROOT = 0;
TTree * MgntMat::TreeElm = 0;
Bool_t  MgntMat::Vac = false;
Float_t MgntMat::Coo[3] = {0};
Bool_t  MgntMat::Elm[Material::NumOfElm] = {false};
Float_t MgntMat::Den[Material::NumOfElm] = {0}; // [mol cm^-3]

void *                  MgntMat::LoadFileBIN = (void*)0;
MgntMat::MaterialBox *  MgntMat::MatBox = 0;
const size_t            MgntMat::MatBoxSize = 3520800 * sizeof(MgntMat::MaterialBox);
const Long64_t          MgntMat::MatDetBIdx[7] = {     0,  64800, 177696, 315680, 338720, 361760,  384800 };
const Long64_t          MgntMat::MatDetSize[7] = { 64800, 112896, 137984,  23040,  23040,  23040, 3136000 };
const Double_t          MgntMat::MatBinStep[3] = { 2.5, 2.5, 0.5 };
const Long64_t          MgntMat::MatBoundIndex[12] = {   6, -1,  5, -1, 4, -1,  3,   -1,   2,  -1,    1,    0 };
const Double_t          MgntMat::MatBoundRange[13] = { 175, 50, 30, 25, 2, -3, -25, -30, -54, -76, -122, -140, -165 };
const Double_t          MgntMat::MatDetRange[7][3] = { { -165, -140,  45 },   // ECAL
    { -140, -122,  70 },   // TRKL9 & RICH PMT
    {  -76,  -54,  70 },   // RICH RAD & LTOF
    {  -30,  -25,  60 },   // TRKL78
    {   -3,    2,  60 },   // TRKL56
    {   25,   30,  60 },   // TRKL34
    {   50,  175, 140 } }; // TRKL2 & UTOF & TRD & TRKL1 & RAD
const Long64_t          MgntMat::MatDetBSize[7][3] = { {  36,  36,  50 },    // ECAL
    {  56,  56,  36 },    // TRKL9 & RICH PMT
    {  56,  56,  44 },    // RICH RAD & LTOF
    {  48,  48,  10 },    // TRKL78
    {  48,  48,  10 },    // TRKL56
    {  48,  48,  10 },    // TRKL34
    { 112, 112, 250 } };  // TRKL2 & UTOF & TRD & TRKL1 & RAD;

Material MgntMat::MergeMat;
Material MgntMat::NullMat;
Material MgntMat::CurMat;
Long64_t MgntMat::CurEntryID = -1;

Bool_t MgntMat::Load() {
    if (LoadStatus) return LoadStatus;
    TString fileDirPath = "/afs/cern.ch/work/h/hchou/public/DATABASE/detector/";
    TString fileROOT = fileDirPath + "g4elmmap_full.root";
    TString fileBIN  = fileDirPath + "g4elmmap.bin";

    LoadFileROOT = TFile::Open(fileROOT.Data());
    Int_t fileDes = open(fileBIN, O_RDONLY);
    LoadFileBIN = (fileDes == -1)	? (void*)-1 : mmap(nullptr, MgntMat::MatBoxSize, PROT_READ, MAP_SHARED, fileDes, 0);

    if (LoadFileROOT && LoadFileBIN != (void*)-1) {
        std::cout << "MgntMat::Load() Open file : " << fileROOT.Data() << std::endl;
        TreeElm = (TTree *) LoadFileROOT->Get("elmmap");
        if (TreeElm) TreeElm->SetBranchAddress("vac", &Vac);
        if (TreeElm) TreeElm->SetBranchAddress("coo",  Coo);
        if (TreeElm) TreeElm->SetBranchAddress("elm",  Elm);
        if (TreeElm) TreeElm->SetBranchAddress("den",  Den);

        std::cout << "MgntMat::Load() Open file : " << fileBIN.Data() << std::endl;
        MatBox = (MgntMat::MaterialBox *) LoadFileBIN;

        LoadStatus = true;
    }
    else {
        if (!LoadFileROOT)            std::cerr << "\nMaterial map not found : " << fileROOT.Data() << std::endl;
        if (LoadFileBIN == (void*)-1) std::cerr << "\nMaterial map not found : " << fileBIN.Data() << std::endl;
    }

    return LoadStatus;
}

Long64_t MgntMat::GetEntry(Double_t x, Double_t y, Double_t z, Bool_t isUsedTree) {
    Long64_t entry = -1;

    if (isUsedTree) {
        if (x <= BinXRangeL || x >= BinXRangeU) return -1;
        if (y <= BinYRangeL || y >= BinYRangeU) return -1;
        if (z <= BinZRangeL || z >= BinZRangeU) return -1;
        Long64_t entryx = Long64_t((x - BinXRangeL) / BinXWidth); 
        Long64_t entryy = Long64_t((y - BinYRangeL) / BinYWidth); 
        Long64_t entryz = Long64_t((z - BinZRangeL) / BinZWidth);
        entry = (entryx) + (entryy * BinXSize) + (entryz * BinXSize * BinYSize);
    }
    else {
        if (z >= MgntMat::MatBoundRange[0] || z <= MgntMat::MatBoundRange[12]) return -1;
        Int_t range[2] = { 0, 12 };
        while ((range[1] - range[0]) != 1) {
            Int_t mid = (range[0] + range[1]) / 2;
            if (z > MgntMat::MatBoundRange[mid] || 
                    MGNumc::Equal(z, MgntMat::MatBoundRange[mid])) range[1] = mid;
            else                                                range[0] = mid;
        }
        Int_t    ibin = range[0];
        Long64_t idet = MgntMat::MatBoundIndex[ibin];
        if (idet == -1) return -1;

        Double_t maxxy = std::max(std::fabs(x), std::fabs(y));
        if (maxxy >= MgntMat::MatDetRange[idet][2]) return -1;
        Long64_t entryx = Long64_t((x + MgntMat::MatDetRange[idet][2]) / MgntMat::MatBinStep[0]); 
        Long64_t entryy = Long64_t((y + MgntMat::MatDetRange[idet][2]) / MgntMat::MatBinStep[1]); 
        Long64_t entryz = Long64_t((z - MgntMat::MatDetRange[idet][0]) / MgntMat::MatBinStep[2]); 
        entry = MgntMat::MatDetBIdx[idet] + 
            (entryx) + 
            (entryy * MgntMat::MatDetBSize[idet][0]) + 
            (entryz * MgntMat::MatDetBSize[idet][0] * MgntMat::MatDetBSize[idet][1]);
    }

    return entry;
}

Material & MgntMat::GetMat(Long64_t entry, Bool_t isUsedTree) {
    if (entry < 0) return NullMat;
    CurMat.init();
    CurEntryID = entry;
    if (isUsedTree) {
        TreeElm->GetEntry(CurEntryID);
        CurMat.fVac = Vac;
        for (Int_t ie = 0; ie < Material::NumOfElm; ++ie) {
            if (!Elm[ie]) continue;
            CurMat.fElm[ie] = true;
            CurMat.fDen[ie] = Den[ie]; // [mol cm^-3]
        }
    }
    else {
        CurMat.fVac = MatBox[CurEntryID].fVac;
        for (Int_t ie = 0; ie < Material::NumOfElm; ++ie) {
            if (!MatBox[CurEntryID].fElm[ie]) continue;
            CurMat.fElm[ie] = true;
            CurMat.fDen[ie] = MatBox[CurEntryID].fDen[ie]; // [mol cm^-3]
        }
    }
    return CurMat;
}

Material & MgntMat::GetMat(Double_t x, Double_t y, Double_t z, Bool_t isUsedTree) {
    if (!Load()) return NullMat;
    Long64_t entry = GetEntry(x, y, z, isUsedTree);
    Material & mat = GetMat(entry, isUsedTree);
    return mat;
}

Material & MgntMat::GetMat(Double_t satx, Double_t saty, Double_t satz, Double_t endx, Double_t endy, Double_t endz, Bool_t isUsedTree, Bool_t isFastScan) {
    if (!Load()) return NullMat;
    MergeMat.init();
    Double_t binScanStep = ((isFastScan) ? MgntMat::BinFastScanStep : MgntMat::BinScanStep);
    Double_t dir[3]  = { (endx - satx), (endy - saty), (endz - satz) };
    Double_t len     = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    Int_t   nstep   = Int_t(len / binScanStep);
    Double_t steplen = (len / Double_t(nstep));
    if (nstep == 0) {
        Double_t cen[3] = { Double_t(0.5) * (endx + satx), Double_t(0.5) * (endy + saty), Double_t(0.5) * (endz + satz) };
        MergeMat = GetMat(cen[0], cen[1], cen[2], isUsedTree); return MergeMat;
    }
    else {
        dir[0] *= (steplen / len);
        dir[1] *= (steplen / len);
        dir[2] *= (steplen / len);
    }

    Double_t wgtfrac = (1. / Double_t(nstep));
    for (Int_t istep = 0; istep < nstep; ++istep) {
        Double_t mag = (0.5 + istep);
        Double_t x = satx + mag * dir[0];
        Double_t y = saty + mag * dir[1];
        Double_t z = satz + mag * dir[2];
        Material & mat = GetMat(x, y, z, isUsedTree);
        if (mat.fVac) continue;
        MergeMat.fVac = false;
        for (Int_t ie = 0; ie < Material::NumOfElm; ++ie) {
            if (!mat.fElm[ie]) continue;
            MergeMat.fElm[ie] = true;
            MergeMat.fDen[ie] += wgtfrac * mat.fDen[ie]; // [mol cm^-3]
        }
    }
    return MergeMat;
}


//---- MatPhyCalParam ----//
class MatPhyCalParam {
    public :
        MatPhyCalParam() { init(); }
        ~MatPhyCalParam() {}

        void init() {
            fVacuum = true;
            fNumRadLen = 0;
            fInvRadLen = 0;
            fMscatL = 0;
            fMscatD = 0;
            fEnglsISGM = 0;
            fEnglsIMPV = 0;
            fEnglsBMEN = 0;
        }

        void print() {
            std::cout << Form("**** MatPhyCalParam ****\n");
            std::cout << Form("Vacuum    : %d\n",     fVacuum);
            std::cout << Form("NumRadLen : %14.8f\n", fNumRadLen);
            std::cout << Form("InvRadLen : %14.8f\n", fInvRadLen);
            std::cout << Form("MscatL    : %14.8f\n", fMscatL);
            std::cout << Form("MscatD    : %14.8f\n", fMscatD);
            std::cout << Form("EnglsISGM : %14.8f\n", fEnglsISGM);
            std::cout << Form("EnglsIMPV : %14.8f\n", fEnglsIMPV);
            std::cout << Form("EnglsBMEN : %14.8f\n", fEnglsBMEN);
            std::cout << Form("\n");
        }

    public :
        Bool_t Vacuum() { return fVacuum; }
        Double_t NumRadLen() { return fNumRadLen; }

    public :
        Bool_t   fVacuum;    // vacuum
        Double_t fNumRadLen; // number of radiation length [1]
        Double_t fInvRadLen; // inverse radiation length [cm^-1]
        Double_t fMscatL;    // multiple-scattering length [1]
        Double_t fMscatD;    // multiple-scattering direction [cm^-1]
        Double_t fEnglsISGM; // ionization-energy-loss SGM [cm^-1]
        Double_t fEnglsIMPV; // ionization-energy-loss MPV [cm^-1]
        Double_t fEnglsBMEN; // bremsstrahlung-energy-loss Mean [cm^-1]
};


//---- MgntMatPhyCal ----//
class MgntMatPhyCal {
    public :
        MgntMatPhyCal() {}
        ~MgntMatPhyCal() {}

        // Material Parameters
        static MatPhyCalParam GetMatPhyCalParam(PhySt & phySt, SVecD<3> & sat, SVecD<3> & end, Bool_t isUsedTree = false, Bool_t isFastScan = true);

        // Number of Radiation Length
        static Double_t GetNumRadLen(SVecD<3> & sat, SVecD<3> & end, Bool_t isUsedTree = false, Bool_t isFastScan = true);

    protected :
        // Inverse Radiation Length [cm^-1]
        static Double_t GetInverseRadiationLength(Material & mat);

        // Multiple-Scattering Length [1]
        static Double_t GetMultipleScatteringLength(PhySt & phySt, Material & mat, Double_t invRadLen, Double_t length);

        // Multiple-Scattering Direction [cm^-1]
        static Double_t GetMultipleScatteringDirection(PhySt & phySt, Material & mat, Double_t invRadLen, Double_t length);
        static Double_t GetMultipleScatteringDirection(Material & mat, Double_t mscatL, Double_t length);

        // Density effect correction [1]
        static void GetDensityEffectCorrection(PhySt & phySt, Material & mat, Double_t * delta);

        // Ionization Energy Loss [cm^-1] // (Normalized by Incident Particle Mass)
        static void GetIonizationEnergyLoss(PhySt & phySt, Material & mat, Double_t length, Double_t & englsSGM, Double_t & englsMPV);

        // Bremsstrahlung Energy Loss [cm^-1] // (Normalized by Incident Particle Mass)
        static Double_t GetBremsstrahlungEnergyLoss(PhySt & phySt, Material & mat, Double_t invRadLen);

    protected :
        // Inverse Radiation Length [g cm^-2]
        static const Double_t AtomRadLen[Material::NumOfElm];

        // Coulomb Multiple Scattering, the Highland-Lynch-Dahl equation
        // Sigma_plane_angle = (RydbergConstant / abs(beta * rigidity) *
        //                     sqrt( radiationLength ) *
        //                     (1. + 0.038 * log(radiationLength)) )
        static constexpr Double_t RygConst    = 0.0136; // [GeV]
        static constexpr Double_t CmsCorrFact = 0.0380; // [1]
        static constexpr Double_t LimitLenL   = 5e-4; // number of radiation length
        static constexpr Double_t LimitLenU   = 100.; // number of radiation length

        // Mean Excitation Energy I = 16eV * Z^0.9 [eV]
        static const Double_t MeanExEng[Material::NumOfElm]; // [MeV]
        static const Double_t NegLnMeanExEng[Material::NumOfElm]; // [MeV]

        // Density effect correction
        static const Double_t DenEffCorr_2ln10 = 4.60517e+00;
        static const Double_t DenEffCorr_C[Material::NumOfElm];
        static const Double_t DenEffCorr_X0[Material::NumOfElm];
        static const Double_t DenEffCorr_X1[Material::NumOfElm];
        static const Double_t DenEffCorr_A[Material::NumOfElm];
        static const Double_t DenEffCorr_M[Material::NumOfElm];
        static const Double_t DenEffCorr_Delta0[Material::NumOfElm];
        static const Double_t DenEffCorr_DeltaM[Material::NumOfElm];

        // Energy Loss from ionization, the Bethe-Bloch equation
        static constexpr Double_t BetheBlochK = 0.307075; // [MeV mol^-1 cm^2]
        static constexpr Double_t LandauEngLossCorr = 0.2;
        static constexpr Double_t MassElInMeV = 0.510999; // [MeV]
        static constexpr Double_t MassElInGeV = 0.000510999; // [GeV]
        static constexpr Double_t BetaLimit = 0.3;
        static constexpr Double_t EtaLimit  = 3.144855e-01;	

        // Bremsstrahlung
        static constexpr Double_t Bremsstrahlung_1oln2 = 1.44269504088896339e+00;

    public :
        // test
        static SMtxSymD<2> GetMscatCov(PhySt & phySt, Double_t zcoo, Bool_t varType = false); // varType := (0 tan) (1 cos)

        static Double_t GetSimplicityNumRadLen(SVecD<3> & sat, SVecD<3> & end);

    protected :
        static Int_t GetSimplicityEntry(Double_t zcoo, Bool_t ulEqOpt = false); // ulEqOpt := (0 up) (1 low)

    protected :
        static const Int_t        NumOfSec = 34;
        static const Bool_t       VacOfSec[NumOfSec];
        static const Double_t     RegOfSec[NumOfSec+1];
        static const Double_t     NrlOfSec[NumOfSec][4]; // ( SAT, END, NRL, DENNRL )

        static const Double_t   Inv3 = (1./3.);
        static const Double_t   InvSqrt3 = 5.77350269189625731e-01;
};

// Element  0"H", 1"C", 2"N", 3"O", 4"F", 5"Na", 6"Al", 7"Si", 8"Pb"
const Double_t MgntMatPhyCal::AtomRadLen[Material::NumOfElm] = 
{ 63.04, 42.70, 37.99, 34.24, 32.93, 27.74, 24.01, 21.82, 6.37 }; // [g cm^-2]

const Double_t MgntMatPhyCal::MeanExEng[Material::NumOfElm] = 
{  1.92e-05,  8.10e-05,  8.20e-05,  9.50e-05,  1.15e-04,  1.49e-04,  1.66e-04,  1.73e-04,  8.23e-04 }; // From NIST, ESTART
const Double_t MgntMatPhyCal::NegLnMeanExEng[Material::NumOfElm] = 
{ 1.086e+01, 9.421e+00, 9.409e+00, 9.262e+00, 9.071e+00, 8.812e+00, 8.704e+00, 8.662e+00, 7.103e+00 }; // From NIST, ESTART

const Double_t MgntMatPhyCal::DenEffCorr_C[Material::NumOfElm] = 
{ 3.2632,  2.9925, 10.5400, 10.7004, 10.9653, 5.0526, 4.2395, 4.4351, 6.2018 }; // Form Geant4
const Double_t MgntMatPhyCal::DenEffCorr_X0[Material::NumOfElm] = 
{ 0.4759, -0.0351,  1.7378,  1.7541,  1.8433, 0.2880, 0.1708, 0.2014, 0.3776 }; // Form Geant4
const Double_t MgntMatPhyCal::DenEffCorr_X1[Material::NumOfElm] = 
{ 1.9215,  2.4860,  4.1323,  4.3213,  4.4096, 3.1962, 3.0127, 2.8715, 3.8073 }; // Form Geant4
const Double_t MgntMatPhyCal::DenEffCorr_A[Material::NumOfElm] = 
{ 0.1348,  0.2024,  0.1535,  0.1178,  0.1108, 0.0777, 0.0802, 0.1492, 0.0936 }; // Form Geant4
const Double_t MgntMatPhyCal::DenEffCorr_M[Material::NumOfElm] = 
{ 5.6249,  3.0036,  3.2125,  3.2913,  3.2962, 3.6452, 3.6345, 3.2546, 3.1608 }; // Form Geant4
const Double_t MgntMatPhyCal::DenEffCorr_Delta0[Material::NumOfElm] = 
{   0.13,     0.1,    0.19,    0.11,    0.11,   0.08,   0.12,   0.14,   0.14 }; // Form Geant4
const Double_t MgntMatPhyCal::DenEffCorr_DeltaM[Material::NumOfElm] = 
{  0.021,   0.038,  0.086,  0.101,     0.121,  0.098,  0.061,  0.059,  0.019 }; // Form Geant4

const Bool_t MgntMatPhyCal::VacOfSec[MgntMatPhyCal::NumOfSec] = {
    // 01. RAD    02. VAC    03. TRKL1   04. VAC    05. TRDU    06. VAC    07. TRD    08. VAC    09. LTRD    10. VAC    
    0,         1,         0,          1,         0,          1,         0,         1,         0,          1,
    // 11. SSTR   12. UTOF   13. SSTR    14. VAC    15. SSTR    16. TRKL2  17. VAC    18. TRKL34 19. VAC     20. TRKL56 
    0,         0,         0,          1,         0,          0,         1,         0,         1,          0,
    // 21. VAC    22. TRKL78 23. VAC     24. SSTR   25. VAC     26. SSTR   27. LTOF   28. SSTR   29. VAC     30. RICH    
    1,         0,         1,          0,         1,          0,         0,         0,         1,          0,
    // 31. VAC    32. TRKL9  33. VAC     34. ECAL    
    1,         0,         1,          0
};

const Double_t MgntMatPhyCal::RegOfSec[MgntMatPhyCal::NumOfSec+1] = {
    //   00        01        02        03        04        05        06        07        08        09
    173.000,  165.000,  164.500,  158.500,  155.500,  146.000,  144.500,   84.500,   83.500,   79.000, // +00
    77.000,   66.500,   58.500,   58.000,   54.500,   54.000,   52.500,   29.500,   25.000,    2.000, // +10
    -2.500,  -25.000,  -29.500,  -54.000,  -54.500,  -58.000,  -58.500,  -66.500,  -72.000,  -73.000, // +20
    -76.000, -134.500, -136.500, -137.000, -164.000                                                    // +30
};

const Double_t MgntMatPhyCal::NrlOfSec[MgntMatPhyCal::NumOfSec][4] = { 
    //     SAT       END          NRL      DENNRL       Section
    {  173.000,  165.000,  0.02704704, 0.00338088 }, // 01. RAD    
    {  165.000,  164.500,  0.00000000, 0.00000000 }, // 02. VAC    
    {  164.500,  158.500,  0.01815394, 0.00302566 }, // 03. TRKL1   
    {  158.500,  155.500,  0.00000000, 0.00000000 }, // 04. VAC    
    {  155.500,  146.000,  0.02348655, 0.00247227 }, // 05. TRDU    
    {  146.000,  144.500,  0.00000000, 0.00000000 }, // 06. VAC    
    {  144.500,   84.500,  0.07611529, 0.00126859 }, // 07. TRD     
    {   84.500,   83.500,  0.00000000, 0.00000000 }, // 08. VAC    
    {   83.500,   79.000,  0.00887603, 0.00197245 }, // 09. LTRD    
    {   79.000,   77.000,  0.00000000, 0.00000000 }, // 10. VAC    
    {   77.000,   66.500,  0.03558093, 0.00338866 }, // 11. SSTR    
    {   66.500,   58.500,  0.05395483, 0.00674435 }, // 12. UTOF    
    {   58.500,   58.000,  0.00363748, 0.00727495 }, // 13. SSTR    
    {   58.000,   54.500,  0.00000000, 0.00000000 }, // 14. VAC    
    {   54.500,   54.000,  0.00290885, 0.00581771 }, // 15. SSTR    
    {   54.000,   52.500,  0.00532545, 0.00355030 }, // 16. TRKL2   
    {   52.500,   29.500,  0.00000000, 0.00000000 }, // 17. VAC    
    {   29.500,   25.000,  0.00979646, 0.00217699 }, // 18. TRKL34  
    {   25.000,    2.500,  0.00000000, 0.00000000 }, // 19. VAC    
    {    2.000,   -2.500,  0.01067055, 0.00237123 }, // 20. TRKL56  
    {   -2.500,  -25.000,  0.00000000, 0.00000000 }, // 21. VAC    
    {  -25.000,  -29.500,  0.00979555, 0.00217679 }, // 22. TRKL78  
    {  -29.500,  -54.000,  0.00000000, 0.00000000 }, // 23. VAC    
    {  -54.000,  -54.500,  0.00290885, 0.00581771 }, // 24. SSTR    
    {  -54.500,  -58.000,  0.00000000, 0.00000000 }, // 25. VAC    
    {  -58.000,  -58.500,  0.00290885, 0.00581771 }, // 26. SSTR    
    {  -58.500,  -66.500,  0.06242154, 0.00780269 }, // 27. LTOF    
    {  -66.500,  -72.000,  0.01784947, 0.00324536 }, // 28. SSTR    
    {  -72.000,  -73.000,  0.00000000, 0.00000000 }, // 29. VAC    
    {  -73.000,  -76.000,  0.01502663, 0.00500888 }, // 30. RICH    
    {  -76.000, -134.500,  0.00000000, 0.00000000 }, // 31. VAC    
    { -134.500, -136.500,  0.00847910, 0.00423955 }, // 32. TRKL9   
    { -136.500, -137.000,  0.00000000, 0.00000000 }, // 33. VAC     
    { -137.000, -164.000, 14.16256368, 0.52453940 }  // 34. ECAL    
};


// Material Parameters
MatPhyCalParam MgntMatPhyCal::GetMatPhyCalParam(PhySt & phySt, SVecD<3> & sat, SVecD<3> & end, Bool_t isUsedTree, Bool_t isFastScan) {
    MatPhyCalParam param;
    Material & mat = MgntMat::GetMat(sat, end, isUsedTree, isFastScan);
    if (mat.fVac) return param;
    param.fVacuum = false;
    Double_t length = LA::Mag(end - sat);
    param.fInvRadLen = GetInverseRadiationLength(mat);
    param.fNumRadLen = length * param.fInvRadLen;
    param.fMscatL = GetMultipleScatteringLength   (phySt, mat, param.fInvRadLen, length);
    param.fMscatD = GetMultipleScatteringDirection(mat, param.fMscatL, length);
    GetIonizationEnergyLoss(phySt, mat, length, param.fEnglsISGM, param.fEnglsIMPV);
    param.fEnglsBMEN = GetBremsstrahlungEnergyLoss(phySt, mat, param.fInvRadLen);

    // Tuning by H.Y.Chou
    param.fMscatL *= 0.00;
    param.fMscatD *= 1.11e+00;

    // Testing by H.Y.Chou
    //param.fMscatD *= 1.11e-3;
    //param.fEnglsIMPV *= 10.0;
    //param.fEnglsISGM *= 10.0;

    return param;
}

Double_t MgntMatPhyCal::GetNumRadLen(SVecD<3> & sat, SVecD<3> & end, Bool_t isUsedTree, Bool_t isFastScan) {
    Material & mat = MgntMat::GetMat(sat, end, isUsedTree, isFastScan);
    if (mat.fVac) return 0.;
    Double_t length = LA::Mag(end - sat);
    Double_t invRadLen = GetInverseRadiationLength(mat);
    Double_t numRadLen = length * invRadLen;
    if (!MGNumc::Valid(numRadLen) || numRadLen < 0) return 0;
    return numRadLen;
}

// Inverse Radiation Length [cm^-1]
Double_t MgntMatPhyCal::GetInverseRadiationLength(Material & mat) {
    if (mat.fVac) return 0;
    Double_t invRadLen = 0;
    for (Int_t ie = 0; ie < Material::NumOfElm; ++ie) {
        if (!mat.fElm[ie]) continue;
        invRadLen += (mat.fDen[ie] * Material::AtomMass[ie] / MgntMatPhyCal::AtomRadLen[ie]);
    }
    if (!MGNumc::Valid(invRadLen) || invRadLen < 0) return 0;
    return invRadLen;
}

// Multiple-Scattering Length [1]
Double_t MgntMatPhyCal::GetMultipleScatteringLength(PhySt & phySt, Material & mat, Double_t invRadLen, Double_t length) {
    if (mat.fVac) return 0;
    if (phySt.IsChrgLess() || phySt.IsMassLess()) return 0;
    if (MGNumc::Compare(invRadLen) <= 0) return 0;
    Double_t numRadLen = length * invRadLen;
    if (numRadLen < MgntMatPhyCal::LimitLenL || numRadLen > MgntMatPhyCal::LimitLenU) return 0;
    Double_t parfact = MgntMatPhyCal::RygConst * std::fabs(phySt.ChrgMass());
    Double_t radfact = std::sqrt(numRadLen) * (1. + MgntMatPhyCal::CmsCorrFact * std::log(numRadLen));
    Double_t mscatL = parfact * radfact;
    if (!MGNumc::Valid(mscatL) || mscatL < 0) return 0;
    return mscatL;
}

// Multiple-Scattering Direction [cm^-1]
Double_t MgntMatPhyCal::GetMultipleScatteringDirection(PhySt & phySt, Material & mat, Double_t invRadLen, Double_t length) {
    if (mat.fVac) return 0;
    if (phySt.IsChrgLess() || phySt.IsMassLess()) return 0;
    if (MGNumc::Compare(invRadLen) <= 0) return 0;
    Double_t numRadLen = length * invRadLen;
    if (numRadLen < MgntMatPhyCal::LimitLenL || numRadLen > MgntMatPhyCal::LimitLenU) return 0;
    Double_t parfact = MgntMatPhyCal::RygConst * std::fabs(phySt.ChrgMass());
    Double_t radfact = (std::sqrt(numRadLen) * (1. + MgntMatPhyCal::CmsCorrFact * std::log(numRadLen))) / length;
    Double_t mscatD = parfact * radfact;
    if (!MGNumc::Valid(mscatD) || mscatD < 0) return 0;
    return mscatD;
}

Double_t MgntMatPhyCal::GetMultipleScatteringDirection(Material & mat, Double_t mscatL, Double_t length) {
    if (mat.fVac) return 0;
    if (MGNumc::Compare(length) <= 0) return 0;
    if (MGNumc::Compare(mscatL) <= 0) return 0;
    Double_t mscatD = mscatL / length;
    if (!MGNumc::Valid(mscatD)) return 0;
    return mscatD;
}

// Density effect correction [1]
void MgntMatPhyCal::GetDensityEffectCorrection(PhySt & phySt, Material & mat, Double_t * delta) {
    if (delta == 0) return;
    std::fill_n(delta, Material::NumOfElm, 0.);
    if (mat.fVac || MGNumc::EqualToZero(phySt.GammaBeta())) return;

    Double_t gammabeta = (phySt.GammaBeta() < MgntMatPhyCal::EtaLimit) ? MgntMatPhyCal::EtaLimit : phySt.GammaBeta();
    Double_t logGB = std::log10(gammabeta);
    for (Int_t ie = 0; ie < Material::NumOfElm; ++ie) {
        if (!mat.fElm[ie]) continue;
        if (logGB < MgntMatPhyCal::DenEffCorr_X0[ie]) continue; // if nonconductors
        Double_t dlogGB = MgntMatPhyCal::DenEffCorr_X1[ie] - logGB;
        delta[ie] = MgntMatPhyCal::DenEffCorr_2ln10 * logGB - MgntMatPhyCal::DenEffCorr_C[ie];
        if (dlogGB < 0) continue;
        delta[ie] += MgntMatPhyCal::DenEffCorr_A[ie] * std::pow(dlogGB, MgntMatPhyCal::DenEffCorr_M[ie]);
    }

    for (Int_t ie = 0; ie < Material::NumOfElm; ++ie)
        if (!MGNumc::Valid(delta[ie]) || delta[ie] < 0) delta[ie] = 0;
}

// Ionization Energy Loss [cm^-1] // (Normalized by Incident Particle Mass)
void MgntMatPhyCal::GetIonizationEnergyLoss(PhySt & phySt, Material & mat, Double_t length, Double_t & englsSGM, Double_t & englsMPV) {
    englsSGM = 0;
    englsMPV = 0;
    if (mat.fVac) return;
    if (phySt.IsChrgLess() || phySt.IsMassLess()) return;
    Double_t gammabeta = (phySt.GammaBeta() < MgntMatPhyCal::EtaLimit) ? MgntMatPhyCal::EtaLimit : phySt.GammaBeta();
    Double_t beta = (phySt.Beta() < MgntMatPhyCal::BetaLimit) ? MgntMatPhyCal::BetaLimit : phySt.Beta();
    Double_t beta2 = beta  * beta;

    Double_t mass  = phySt.Mass() * 1.0e3; // [MeV]
    Double_t maxKEng = (2. * MgntMatPhyCal::MassElInMeV * gammabeta * gammabeta);
    Double_t sgmfact = (0.5 * MgntMatPhyCal::BetheBlochK) * (phySt.Chrg() * phySt.Chrg());

    Double_t dlt[Material::NumOfElm] = {0};
    GetDensityEffectCorrection(phySt, mat, dlt);

    Double_t elmwgt = 0;
    Double_t mexeng = 0;
    Double_t avgdlt = 0;
    for (Int_t ie = 0; ie < Material::NumOfElm; ++ie) {
        if (!mat.fElm[ie]) continue;
        Double_t wgt = Material::AtomChrg[ie] * mat.fDen[ie];
        mexeng += wgt * MgntMatPhyCal::NegLnMeanExEng[ie];
        avgdlt += wgt * dlt[ie];
        elmwgt += wgt;
    }
    mexeng /= elmwgt;
    avgdlt /= elmwgt;
    if (!MGNumc::Valid(elmwgt) || elmwgt < 0) elmwgt = 0;
    if (!MGNumc::Valid(mexeng) || mexeng < 0) mexeng = 0;
    if (!MGNumc::Valid(avgdlt) || avgdlt < 0) avgdlt = 0;

    Double_t sgm = sgmfact * elmwgt * length;
    //sgm *= (1.327457e-01*TMath::Erf((1.32475e+00*TMath::Log(gammabeta)-1.99679e+00))+5.052637e-01); // New Tuning by H.Y.Chou (Oct 28, 2016)
    //sgm *= (1.896367e-01*TMath::Erf((1.32475e+00*TMath::Log(gammabeta)-1.99679e+00))+4.294749e-01); // New Tuning by H.Y.Chou (Oct 28, 2016)

    Double_t mpv = sgm * (std::log(maxKEng * sgm / beta2) + (2. * mexeng) + MgntMatPhyCal::LandauEngLossCorr - beta2 - avgdlt);
    //mpv *= (5.64008e-01*TMath::Erfc((1.07651e+00*TMath::Log(gammabeta)-1.14087e+00))+1.07237e+00); // New Tuning by H.Y.Chou (Oct 28, 2016)
    //mpv *= (1.407764e+00*TMath::Erfc((1.07651e+00*TMath::Log(gammabeta)-1.14087e+00))+1.07237e+00); // New Tuning by H.Y.Chou (Oct 28, 2016)

    // Tuning by H.Y.Chou
    sgm *= 7.025551e-01;
    sgm /= (0.668148*(TMath::Erfc(0.971628*(TMath::Log(gammabeta)-0.939321))+1.51119));
    mpv *= 6.649170e-01;
    mpv /= (1.82142*(TMath::Erf(1.72868*(TMath::Power(gammabeta,0.246414)-0.30212))-0.454936));

    englsSGM = (sgm / length / mass);
    englsMPV = (mpv / length / mass);
    if (!MGNumc::Valid(englsSGM) || englsSGM < 0) englsSGM = 0;
    if (!MGNumc::Valid(englsMPV) || englsMPV < 0) englsMPV = 0;
}

// Bremsstrahlung Energy Loss [cm^-1]
Double_t MgntMatPhyCal::GetBremsstrahlungEnergyLoss(PhySt & phySt, Material & mat, Double_t invRadLen) {
    if (mat.fVac) return 0;	
    if (phySt.IsChrgLess() || phySt.IsMassLess()) return 0;
    if (MGNumc::Compare(invRadLen) <= 0) return 0;
    Double_t frac = (phySt.ChrgMass() * phySt.ChrgMass() * MgntMatPhyCal::MassElInGeV * MgntMatPhyCal::MassElInGeV);
    Double_t englsMEN = frac * MgntMatPhyCal::Bremsstrahlung_1oln2 * invRadLen; 
    if (!MGNumc::Valid(englsMEN) || englsMEN < 0) englsMEN = 0;
    return englsMEN;
}

// Multiple-Scattering Convariance
SMtxSymD<2> MgntMatPhyCal::GetMscatCov(PhySt & phySt, Double_t zcoo, Bool_t varType) { // varType := (0 tan) (1 cos)
    SMtxSymD<2> mscat;
    Double_t dlenz = (zcoo - phySt.Z());
    if (phySt.IsChrgLess() || MGNumc::EqualToZero(dlenz)) return mscat;
    Int_t sdir = ((MGNumc::Compare(dlenz) <= 0) ? 1 : -1); // (1 down-going) (-1 up-going)
    Int_t ssat = GetSimplicityEntry(phySt.Z(), (sdir==1?0:1));
    Int_t send = GetSimplicityEntry(zcoo,      (sdir==1?1:0));
    if (ssat == -1 || send == -1) return mscat;

    Double_t tx  = phySt.TanX();
    Double_t ty  = phySt.TanY();
    Double_t ux  = phySt.DirX();
    Double_t uy  = phySt.DirY();
    Double_t dsz = std::fabs(1. / phySt.DirZ());

    Bool_t   vacuum = true;
    Double_t parfact = MgntMatPhyCal::RygConst * std::fabs(phySt.ChrgMass() * phySt.InvEta());
    Double_t cxxfact = ((varType) ? (1. + tx * tx) : (1. - uy * uy)) * dsz * dsz;
    Double_t cyyfact = ((varType) ? (1. + ty * ty) : (1. - ux * ux)) * dsz * dsz;
    Double_t cxyfact = ((varType) ? (     tx * ty) : (     ux * uy)) * dsz * dsz;
    for (Int_t isec = ssat; ((sdir == 1) ? (isec <= send) : (isec >= send)); isec+=sdir) {
        if (isec == 0 || isec == MgntMatPhyCal::NumOfSec) continue;
        if (MgntMatPhyCal::VacOfSec[isec-1]) continue;
        Bool_t satInSec = (isec == ssat);
        Bool_t endInSec = (isec == send);
        Int_t idx = isec - 1;

        Double_t reg[2] = { 
            (satInSec ? phySt.Z() : MgntMatPhyCal::NrlOfSec[idx][(sdir==1?0:1)]), 
            (endInSec ? zcoo      : MgntMatPhyCal::NrlOfSec[idx][(sdir==1?1:0)]) 
        };

        Double_t nrl = NrlOfSec[idx][2];
        if (satInSec || endInSec) {
            Double_t stp = std::fabs(reg[1] - reg[0]);
            Double_t rat = NrlOfSec[idx][3];
            nrl = stp * rat;
        }
        nrl *= dsz;

        Double_t lenz2 = (reg[1] - reg[0]) * (reg[1] - reg[0]);
        Double_t exlz2 = (zcoo - reg[1]) * (zcoo - reg[1]);
        Double_t radfact = std::sqrt(nrl) * (1. + MgntMatPhyCal::CmsCorrFact * std::log(nrl));
        Double_t mscatFT = parfact * radfact;
        Double_t sigma2 = mscatFT * mscatFT * (lenz2 / 3. + exlz2);

        mscat(0, 0) += cxxfact * sigma2;
        mscat(1, 1) += cyyfact * sigma2;
        mscat(0, 1) += cxyfact * sigma2;

        vacuum = false;
    }
    if (vacuum) return mscat;

    if (!MGNumc::Valid(mscat(0, 0)) || 
            !MGNumc::Valid(mscat(1, 1)) ||
            !MGNumc::Valid(mscat(0, 1))) mscat = SMtxSymD<2>();

    return mscat;
}

Double_t MgntMatPhyCal::GetSimplicityNumRadLen(SVecD<3> & sat, SVecD<3> & end) {
    Double_t NRL = 0.;
    SVecD<3> len = end - sat;
    Double_t dsz = LA::Mag(len) / std::fabs(len(2));
    if (MGNumc::EqualToZero(len(2))) return NRL;
    Int_t sdir = ((MGNumc::Compare(len(2)) <= 0) ? 1 : -1); // (1 down-going) (-1 up-going)
    Int_t ssat = GetSimplicityEntry(sat(2), (sdir==1?0:1));
    Int_t send = GetSimplicityEntry(end(2), (sdir==1?1:0));
    if (ssat == -1 || send == -1) return NRL;

    for (Int_t isec = ssat; ((sdir == 1) ? (isec <= send) : (isec >= send)); isec+=sdir) {
        if (isec == 0 || isec == MgntMatPhyCal::NumOfSec) continue;
        if (MgntMatPhyCal::VacOfSec[isec-1]) continue;
        Bool_t satInSec = (isec == ssat);
        Bool_t endInSec = (isec == send);
        Int_t idx = isec - 1;

        Double_t reg[2] = { 
            (satInSec ? sat(2) : MgntMatPhyCal::NrlOfSec[idx][(sdir==1?0:1)]), 
            (endInSec ? end(2) : MgntMatPhyCal::NrlOfSec[idx][(sdir==1?1:0)]) 
        };

        Double_t nrl = dsz * NrlOfSec[idx][2];
        if (satInSec || endInSec) {
            Double_t stp = std::fabs(reg[1] - reg[0]);
            Double_t rat = NrlOfSec[idx][3];
            nrl = dsz * stp * rat;
        }
        NRL += nrl;
    }
    if (!MGNumc::Valid(NRL)) NRL = 0.; 

    return NRL;
}

Int_t MgntMatPhyCal::GetSimplicityEntry(Double_t zcoo, Bool_t ulEqOpt) { // ulEqOpt := (0 up) (1 low)
    Int_t entry = -1;
    if      ((zcoo > MgntMatPhyCal::RegOfSec[0]) || 
            ( ulEqOpt && MGNumc::Equal(zcoo, MgntMatPhyCal::RegOfSec[0]))
            ) entry = 0;
    else if ((zcoo < MgntMatPhyCal::RegOfSec[MgntMatPhyCal::NumOfSec]) ||
            (!ulEqOpt && MGNumc::Equal(zcoo, MgntMatPhyCal::RegOfSec[MgntMatPhyCal::NumOfSec]))
            ) entry = MgntMatPhyCal::NumOfSec+1;
    else {
        Int_t range[2] = { 0, MgntMatPhyCal::NumOfSec };
        while ((range[1] - range[0]) != 1) {
            Int_t mid = (range[0] + range[1]) / 2;
            Bool_t eq = MGNumc::Equal(zcoo, MgntMatPhyCal::RegOfSec[mid]);
            if ((zcoo > MgntMatPhyCal::RegOfSec[mid] && !eq) || (ulEqOpt && eq)) range[1] = mid;
            else                                                                 range[0] = mid;
        }
        entry = range[0] + 1;
    }
    return entry;
}

#endif // __MgntMat_H__
