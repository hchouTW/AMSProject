#include <CPPLibs.h>
#include <ROOTLibs.h>
#include <TrSys.h>

//#include "inc/SMGFit.h"
//#include "inc/ALGFit.h"

double SMG2func(double *x, double *p) {
    return TrSys::MultiGausFunc(x[0], 0.0, std::vector<std::array<double, 2>>{ {p[0], p[1]}, {p[2], p[3]} });
}
double SMG4func(double *x, double *p) {
    return TrSys::MultiGausFunc(x[0], 0.0, std::vector<std::array<double, 2>>{ {p[0], p[1]}, {p[2], p[3]}, {p[4], p[5]}, {p[6], p[7]} });
}
double ALGfunc(double *x, double *p) {
    return TrSys::LandauGausFunc(x[0], p[0], p[1], p[2], p[3]);
}

int main() {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
  
    const long MIN_NUM = 3;
    const long MAX_NUM = 25;
    Axis AXNUM("Number of Hits (N)", MAX_NUM-MIN_NUM+1, MIN_NUM, MAX_NUM+1);
    Axis AXSMG("Parameter", 4000, -2.0e-2, 2.0e-2);
    Axis AXALG("Parameter", 2000,     0.5,    5.5);
    
    std::vector<double> vals(25, 0.0);
    TFile* ifle = TFile::Open("out/nogaus_data.root");
   

    std::cerr << Form("==== SMG : LX ====\n");
    TTree* treeSMGLX = (TTree*) ifle->Get("SMGLX");
    treeSMGLX->SetBranchAddress("vals", vals.data());

    Hist* hSMGLXnum = Hist::New("hSMGLXnum", HistAxis(AXSMG));
    Hist* hSMGLXstd = Hist::New("hSMGLXstd", HistAxis(AXNUM, AXSMG));
    Hist* hSMGLXsmg = Hist::New("hSMGLXsmg", HistAxis(AXNUM, AXSMG));
    
    const std::array<double, 6> SMGLXpars({
        9.99196e-01, 1.88412e-03,
        3.89323e-01, 3.32112e-03,
        2.71260e-02, 6.98928e-03
    });
    TrSys::MultiGaus PDF_SMGLXstd(TrSys::Robust(), 0.0, SMGLXpars[1]);
    TrSys::MultiGaus PDF_SMGLXsmg(TrSys::Robust(), 0.0, 
        SMGLXpars[0], SMGLXpars[1], 
        SMGLXpars[2], SMGLXpars[3],
        SMGLXpars[4], SMGLXpars[5]
    );

    for (long entry = 0; entry < treeSMGLX->GetEntries(); ++entry) {
        if (entry%10000 == 0) std::cerr << Form("%ld\n", entry/10000);
        if (entry == 10000) break;
        treeSMGLX->GetEntry(entry);
        
        for (auto&& val : vals) hSMGLXnum->fillH1D(val);
        for (long num = MIN_NUM; num <= MAX_NUM; ++num) {
            std::vector<double> subvals(vals.begin(), vals.begin() + num - 1);
            
            TrSys::SimpleSMGFit stdfit(subvals, PDF_SMGLXstd);
            TrSys::SimpleSMGFit smgfit(subvals, PDF_SMGLXsmg);
            if (!stdfit.status() || !smgfit.status()) continue;

            hSMGLXstd->fillH2D(num, stdfit.param()); 
            hSMGLXsmg->fillH2D(num, smgfit.param()); 
        }
    }


    std::cerr << Form("==== SMG : LY ====\n");
    TTree* treeSMGLY = (TTree*) ifle->Get("SMGLY");
    treeSMGLY->SetBranchAddress("vals", vals.data());

    Hist* hSMGLYnum = Hist::New("hSMGLYnum", HistAxis(AXSMG));
    Hist* hSMGLYstd = Hist::New("hSMGLYstd", HistAxis(AXNUM, AXSMG));
    Hist* hSMGLYsmg = Hist::New("hSMGLYsmg", HistAxis(AXNUM, AXSMG));
    
    const std::array<double, 8> SMGLYpars({
        9.26431e-01, 7.50895e-04, 
        4.78682e-01, 1.37797e-03,
        8.21875e-02, 2.72244e-03,
        2.85380e-03, 6.83937e-03
    });
    TrSys::MultiGaus PDF_SMGLYstd(TrSys::Robust(), 0.0, SMGLYpars[1]);
    TrSys::MultiGaus PDF_SMGLYsmg(TrSys::Robust(), 0.0, 
        SMGLYpars[0], SMGLYpars[1], 
        SMGLYpars[2], SMGLYpars[3],
        SMGLYpars[4], SMGLYpars[5],
        SMGLYpars[6], SMGLYpars[7]
    );

    for (long entry = 0; entry < treeSMGLY->GetEntries(); ++entry) {
        if (entry%10000 == 0) std::cerr << Form("%ld\n", entry/10000);
        if (entry == 10000) break;
        treeSMGLY->GetEntry(entry);
        
        for (auto&& val : vals) hSMGLYnum->fillH1D(val);
        for (long num = MIN_NUM; num <= MAX_NUM; ++num) {
            std::vector<double> subvals(vals.begin(), vals.begin() + num);
            
            TrSys::SimpleSMGFit stdfit(subvals, PDF_SMGLYstd);
            TrSys::SimpleSMGFit smgfit(subvals, PDF_SMGLYsmg);
            if (!stdfit.status() || !smgfit.status()) continue;

            hSMGLYstd->fillH2D(num, stdfit.param()); 
            hSMGLYsmg->fillH2D(num, smgfit.param()); 
        }
    }


    std::cerr << Form("==== ALG : EL ====\n");
    TTree* treeALGEL = (TTree*) ifle->Get("ALGEL");
    treeALGEL->SetBranchAddress("vals", vals.data());

    Hist* hALGELnum = Hist::New("hALGELnum", HistAxis(AXALG));
    Hist* hALGELstd = Hist::New("hALGELstd", HistAxis(AXNUM, AXALG));
    Hist* hALGELalg = Hist::New("hALGELalg", HistAxis(AXNUM, AXALG));

    const std::array<double, 4> ALGELpars({ 
        4.07524e-01, 2.14736e+00, 2.40172e-01, 1.18389e-01
    });
    TrSys::MultiGaus  PDF_ALGELstd(TrSys::Robust(), 0.0, std::hypot(ALGELpars[2], ALGELpars[3]));
    TrSys::LandauGaus PDF_ALGELalg(TrSys::Robust(), ALGELpars[0], 0.0, ALGELpars[2], ALGELpars[3]);

    for (long entry = 0; entry < treeALGEL->GetEntries(); ++entry) {
        if (entry%10000 == 0) std::cerr << Form("%ld\n", entry/10000);
        if (entry == 10000) break;
        treeALGEL->GetEntry(entry);
        
        for (auto&& val : vals) hALGELnum->fillH1D(val);
        for (long num = MIN_NUM; num <= MAX_NUM; ++num) {
            std::vector<double> subvals(vals.begin(), vals.begin() + num - 1);
            
            TrSys::SimpleSMGFit stdfit(subvals, PDF_ALGELstd);
            TrSys::SimpleALGFit algfit(subvals, PDF_ALGELalg);
            if (!stdfit.status() || !algfit.status()) continue;

            hALGELstd->fillH2D(num, stdfit.param()); 
            hALGELalg->fillH2D(num, algfit.param()); 
        }
    }


    std::cerr << Form("==== ALG : EM ====\n");
    TTree* treeALGEM = (TTree*) ifle->Get("ALGEM");
    treeALGEM->SetBranchAddress("vals", vals.data());

    Hist* hALGEMnum = Hist::New("hALGEMnum", HistAxis(AXALG));
    Hist* hALGEMstd = Hist::New("hALGEMstd", HistAxis(AXNUM, AXALG));
    Hist* hALGEMalg = Hist::New("hALGEMalg", HistAxis(AXNUM, AXALG));

    const std::array<double, 4> ALGEMpars({ 
        2.19427e-02, 1.39837e+00, 1.01428e-01, 1.13978e-01
    });
    TrSys::MultiGaus  PDF_ALGEMstd(TrSys::Robust(), 0.0, std::hypot(ALGEMpars[2], ALGEMpars[3]));
    TrSys::LandauGaus PDF_ALGEMalg(TrSys::Robust(), ALGEMpars[0], 0.0, ALGEMpars[2], ALGEMpars[3]);

    for (long entry = 0; entry < treeALGEM->GetEntries(); ++entry) {
        if (entry%10000 == 0) std::cerr << Form("%ld\n", entry/10000);
        if (entry == 10000) break;
        treeALGEM->GetEntry(entry);
        
        for (auto&& val : vals) hALGEMnum->fillH1D(val);
        for (long num = MIN_NUM; num <= MAX_NUM; ++num) {
            std::vector<double> subvals(vals.begin(), vals.begin() + num - 1);
            
            TrSys::SimpleSMGFit stdfit(subvals, PDF_ALGEMstd);
            TrSys::SimpleALGFit algfit(subvals, PDF_ALGEMalg);
            if (!stdfit.status() || !algfit.status()) continue;

            hALGEMstd->fillH2D(num, stdfit.param()); 
            hALGEMalg->fillH2D(num, algfit.param()); 
        }
    }


    std::cerr << Form("==== ALG : EH ====\n");
    TTree* treeALGEH = (TTree*) ifle->Get("ALGEH");
    treeALGEH->SetBranchAddress("vals", vals.data());

    Hist* hALGEHnum = Hist::New("hALGEHnum", HistAxis(AXALG));
    Hist* hALGEHstd = Hist::New("hALGEHstd", HistAxis(AXNUM, AXALG));
    Hist* hALGEHalg = Hist::New("hALGEHalg", HistAxis(AXNUM, AXALG));

    const std::array<double, 4> ALGEHpars({ 
        4.82158e-03, 1.12452e+00, 8.38657e-02, 9.36490e-02  
    });
    TrSys::MultiGaus  PDF_ALGEHstd(TrSys::Robust(), 0.0, std::hypot(ALGEHpars[2], ALGEHpars[3]));
    TrSys::LandauGaus PDF_ALGEHalg(TrSys::Robust(), ALGEHpars[0], 0.0, ALGEHpars[2], ALGEHpars[3]);

    for (long entry = 0; entry < treeALGEH->GetEntries(); ++entry) {
        if (entry%10000 == 0) std::cerr << Form("%ld\n", entry/10000);
        if (entry == 10000) break;
        treeALGEH->GetEntry(entry);
        
        for (auto&& val : vals) hALGEHnum->fillH1D(val);
        for (long num = MIN_NUM; num <= MAX_NUM; ++num) {
            std::vector<double> subvals(vals.begin(), vals.begin() + num - 1);
            
            TrSys::SimpleSMGFit stdfit(subvals, PDF_ALGEHstd);
            TrSys::SimpleALGFit algfit(subvals, PDF_ALGEHalg);
            if (!stdfit.status() || !algfit.status()) continue;

            hALGEHstd->fillH2D(num, stdfit.param()); 
            hALGEHalg->fillH2D(num, algfit.param()); 
        }
    }

    ifle->Close();

    TFile * ofle = new TFile("out/nogaus_hist.root", "RECREATE");
    ofle->cd();

    (*hSMGLXnum)()->Write();
    (*hSMGLXstd)()->Write();
    (*hSMGLXsmg)()->Write();
    
    (*hSMGLYnum)()->Write();
    (*hSMGLYstd)()->Write();
    (*hSMGLYsmg)()->Write();
    
    (*hALGELnum)()->Write();
    (*hALGELstd)()->Write();
    (*hALGELalg)()->Write();
    
    (*hALGEMnum)()->Write();
    (*hALGEMstd)()->Write();
    (*hALGEMalg)()->Write();

    (*hALGEHnum)()->Write();
    (*hALGEHstd)()->Write();
    (*hALGEHalg)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}

/*
FCN=1113.36 FROM MIGRAD    STATUS=CONVERGED     111 CALLS         112 TOTAL
EDM=2.8215e-07    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           9.99196e-01   1.00365e-01   1.93522e-03   2.02760e-03
2  p1           1.88412e-03   1.03164e-05   2.53366e-06  -5.75100e-01
3  p2           3.89323e-01   1.78710e-02   6.03515e-05  -7.57883e-02
4  p3           3.32112e-03   3.25429e-05   3.07783e-06  -9.98096e-01
5  p4           2.71260e-02   1.73197e-03   1.96203e-05  -6.71387e-02
6  p5           6.98928e-03   8.04874e-05   1.05542e-05  -2.11808e-01
7  p6          -6.52651e-06   2.12891e-06   3.47000e-06  -1.87515e-01
8  p7           8.38696e-03   1.13485e-05   1.11949e-07   6.50503e+00


FCN=2549.16 FROM MIGRAD    STATUS=CONVERGED     164 CALLS         165 TOTAL
EDM=1.94554e-08    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           4.78682e-01   7.37566e-02   1.18522e-04   8.35804e-04
2  p1           1.37797e-03   1.83575e-05   3.11005e-06   1.47638e-01
3  p2           8.21875e-02   1.35245e-02   3.17199e-05   1.57474e-02
4  p3           2.72244e-03   3.56429e-05   6.16318e-06  -3.75867e-02
5  p4           2.85380e-03   4.95330e-04   1.72952e-05   2.66679e-03
6  p5           6.83937e-03   1.30370e-04   3.10273e-05  -4.23147e-03
7  p6           9.26431e-01   1.26574e-01   2.38865e-04  -1.59470e-03
8  p7           7.50895e-04   5.50307e-06   3.05995e-06   4.77938e-02
9  p8          -1.12822e-04   9.20708e-07   2.26516e-06  -3.81866e-02
10  p9           1.85680e-02   2.54554e-05   3.49628e-07   1.36540e+00


FCN=561.809 FROM MIGRAD    STATUS=CONVERGED      65 CALLS          66 TOTAL
EDM=1.43674e-07    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           4.07524e-01   2.21759e-02   6.69885e-05   1.56047e-03
2  p1           2.14736e+00   2.93272e-03   2.11130e-06   4.05238e-01
3  p2           2.40172e-01   5.45350e-03   3.15483e-06   2.77808e-01
4  p3           1.18389e-01   7.39729e-03   1.27997e-05   6.54317e-02
5  p4           1.10024e-02   1.65592e-04   3.33370e-07   2.61687e+00

FCN=559.187 FROM MIGRAD    STATUS=CONVERGED      59 CALLS          60 TOTAL
EDM=2.14883e-14    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           2.19427e-02   9.77297e-04  -0.00000e+00   2.56466e-05
2  p1           1.39837e+00   6.78995e-04  -0.00000e+00   1.18909e-05
3  p2           1.01428e-01   1.33337e-03  -0.00000e+00  -5.42759e-06
4  p3           1.13978e-01   1.36265e-03   0.00000e+00  -1.55312e-05
5  p4           1.83862e-02   1.91705e-04   0.00000e+00  -6.37058e-04


FCN=506.369 FROM MIGRAD    STATUS=CONVERGED      64 CALLS          65 TOTAL
EDM=5.00253e-10    STRATEGY= 1      ERROR MATRIX ACCURATE
EXT PARAMETER                                   STEP         FIRST
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
1  p0           4.82158e-03   1.85709e-04   1.48051e-05   9.32746e-04
2  p1           1.12452e+00   5.88167e-04   1.76892e-06   9.96169e-02
3  p2           8.38657e-02   7.49960e-04   2.40292e-06   6.49245e-02
4  p3           9.36490e-02   9.53060e-04   6.54752e-06   2.03025e-02
5  p4           2.04397e-02   1.67304e-04   6.19620e-07   2.52229e-01
*/
