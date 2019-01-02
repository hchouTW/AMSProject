{
	gInterpreter->AddIncludePath("${AMSSRC}/include");
	gSystem->Load("${AMSLIB}/libntuple_slc6_PG_dynamic.so");

    Printf("Current Version    = %s", gROOT->GetVersion());
	Printf("Build Architecture = %s", gSystem->GetBuildArch());
    
    TF1 * f1gs = new TF1("f1gs", "TMath::Abs([0])*TMath::Exp(-0.5*(x-[2])*(x-[2])/[1]/[1])");
    f1gs->SetParameters(10., 10., 0.);
    f1gs->SetNpx(100000);

    TF1 * f2gs = new TF1("f2gs", "TMath::Abs([0])*TMath::Exp(-0.5*(x-[4])*(x-[4])/[1]/[1]) + TMath::Abs([2])*TMath::Exp(-0.5*(x-[4])*(x-[4])/[3]/[3])");
    f2gs->SetParameters(10., 10., 10., 30., 0.);
    f2gs->SetNpx(100000);

    TF1 * f3gs = new TF1("f3gs", "TMath::Abs([0])*TMath::Exp(-0.5*(x-[6])*(x-[6])/[1]/[1]) + TMath::Abs([2])*TMath::Exp(-0.5*(x-[6])*(x-[6])/[3]/[3]) + TMath::Abs([4])*TMath::Exp(-0.5*(x-[6])*(x-[6])/[5]/[5])");
    f3gs->SetParameters(10., 10., 10., 30., 10., 50., 0.);
    f3gs->SetNpx(100000);
    
    TF1 * f4gs = new TF1("f4gs", "TMath::Abs([0])*TMath::Exp(-0.5*(x-[8])*(x-[8])/[1]/[1]) + TMath::Abs([2])*TMath::Exp(-0.5*(x-[8])*(x-[8])/[3]/[3]) + TMath::Abs([4])*TMath::Exp(-0.5*(x-[8])*(x-[8])/[5]/[5]) + TMath::Abs([6])*TMath::Exp(-0.5*(x-[8])*(x-[8])/[7]/[7])");
    f4gs->SetParameters(10., 10., 10., 30., 10., 50., 10., 70., 0.);
    f4gs->SetNpx(100000);
    
    TF1* flg = new TF1("flg", "[0] * TMath::Exp(TMath::Abs([1]) * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) + (1-TMath::Abs([1])) * TMath::Log(TMath::Landau(1.17741002*(x-[2])/TMath::Abs([3])-2.22782980e-01)/1.80655634e-01))");
    flg->SetParameters(1.0, 0.1, 0.0, 1.0);
    flg->SetParLimits(1, 0.0, 1.0);
    
    TF1* flgnum = new TF1("flgnum", "[0] * TMath::Exp(TMath::Abs([1]) * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) - (1-TMath::Abs([1])) * (4.90120e-01*(((x-[2])/[3])*1.16366e+00+TMath::Exp(-((x-[2])/[3])*1.16366e+00)-1) - 3.99384e+00*(((x-[2])/[3])*1.27500e-01+TMath::Exp(-((x-[2])/[3])*1.27500e-01)-1)))");
    flgnum->SetParameters(1.0, 0.1, 0.0, 1.0);
    flgnum->SetParLimits(1, 0.0, 1.0);
    
    TF1* fkpa = new TF1("fkpa", "0.5 * (1.0 + TMath::Erf([0] * TMath::Log(1+[1]*(1+x*x)^[2]) - [3]))");
    fkpa->SetParameters(1.0, 2.0, 2.0, 2.0);
    //
    TF1* fmpv = new TF1("fmpv", "[0] + [1] * (1+x*x)^[2]  - [3] * TMath::Log([4]+(x*x)^[5])");
    fmpv->SetParameters(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    
    TF1* fmpv2 = new TF1("fmpv2", "[0] + [1] * (1+x*x)  - [2] * TMath::Log([3]+(x*x))");
    fmpv2->SetParameters(1.0, 1.0, 1.0, 1.0);
    
    TF1* fmpv3 = new TF1("fmpv3", "[0] + [1] * (1+x*x) - [2] * TMath::Log([3]+(x*x)) + [4] * 0.5 * TMath::Erfc([5] * log(x*x) + [6])");
    fmpv3->SetParameters(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    
    TF1* fmpv4 = new TF1("fmpv4", "[0] + [1] * (1+x*x) - [2] * TMath::Log(1.0e-06+(x*x)) + [3] * 0.5 * TMath::Erfc([4] * log(x*x) + [5])");
    fmpv4->SetParameters(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    
    TF1* ftme = new TF1("ftme", "[0] + [1] * TMath::Erfc([2] * (1+x*x)^[3] - [4])");
    ftme->SetParameters(2.13914e+00, 2.40467e-01, 7.71764e+01, 1.70279e-02, 7.76559e+01);
    
    TF1* fIkpa = new TF1("fIkpa", "0.5*(1+TMath::Erf([0]*log(x*x)-[1])) + [2]*0.5*TMath::Erfc([3]*log(x*x)+[4])");
    fIkpa->SetParameters(3.03451e-01, 1.14913e+00, 1.19150e-01, 1.38026e-01, 7.91603e-01);
    
    TF1* fImpv = new TF1("fImpv", "[0]*(x*x+1)^[1] + [2]*0.5*TMath::Erfc([3]*log(x*x)+[4])");
    fImpv->SetParameters(1.70001e+00, 1.03690e+00, 2.12635e+00, 1.44258e-01, 9.07120e-01);
    
    TF1* fITkpa = new TF1("fITkpa", "0.5*(1+TMath::Erf([0]*log(x*x)-[1])) + [2]*0.5*TMath::Erfc([3]*log(x*x)+[4]) + [5]*0.5*TMath::Erfc([6]*log(x*x)+[7])");
    fITkpa->SetParameters(3.03451e-01, 1.14913e+00, 1.19150e-01, 1.38026e-01, 7.91603e-01, 4.31403e-01, 8.23065e-01, 1.15371e+01);
    
    TF1* fITmpv = new TF1("fITmpv", "[0]*(x*x+1)^[1] + [2]*0.5*TMath::Erfc([3]*log(x*x)+[4]) + [5]*0.5*TMath::Erfc([6]*log(x*x)+[7])");
    fITmpv->SetParameters(1.70001e+00, 1.03690e+00, 2.12635e+00, 1.44258e-01, 9.07120e-01, 6.70935e+00, 5.02516e-01, 7.28139e+00);

    //TF1* flgnum2 = new TF1("flgnum2", "[0] * TMath::Exp(TMath::Abs([1]) * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) - (1-TMath::Abs([1])) * (4.90120e-01*(((x-[2])/[3])*1.16366e+00+TMath::Exp(-((x-[2])/[3])*1.16366e+00)-1) - 3.99384e+00*(((x-[2])/[3])*1.27500e-01+TMath::Exp(-((x-[2])/[3])*1.27500e-01)-1)))     +     [4] * TMath::Power((x/[7]), [5]-1) * TMath::Exp(-[6]*(x/[7])) * (2.0-TMath::Erfc((x-[8])/[9])) ");
    TF1* flgnum2 = new TF1("flgnum2", "[0] * TMath::Exp(TMath::Abs([1]) * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) - (1-TMath::Abs([1])) * (4.90120e-01*(((x-[2])/[3])*1.16366e+00+TMath::Exp(-((x-[2])/[3])*1.16366e+00)-1) - 3.99384e+00*(((x-[2])/[3])*1.27500e-01+TMath::Exp(-((x-[2])/[3])*1.27500e-01)-1)))     +     [4] * TMath::Power(x, [5]-1) * TMath::Exp(-[6]*(x)) * 0.5 * (1+TMath::Erf( [7] * (x - [8]*([5]-1)/[6]))) ");
    flgnum2->SetParameters(2.07982e+04, 6.64481e-03, 2.26697e+00, 7.94455e-01, 1.51044e+02, 4.12811e+00, 3.84824e-01, 6.87232e-01, 7.66221e-01);
    
    
    TF1* flgnum3 = new TF1("flgnum3", "[0] * (1-[1]) * TMath::Exp(TMath::Abs([2]) * (-0.5)*((x-[3])*(x-[3])/[4]/[4]) - (1-TMath::Abs([2])) * (4.90120e-01*(((x-[3])/[4])*1.16366e+00+TMath::Exp(-((x-[3])/[4])*1.16366e+00)-1) - 3.99384e+00*(((x-[3])/[4])*1.27500e-01+TMath::Exp(-((x-[3])/[4])*1.27500e-01)-1)))     +     [0] * [1] * TMath::Exp(TMath::Abs([5]) * (-0.5)*((x-[6])/[7] - [8])*((x-[6])/[7] - [8]) - (1-TMath::Abs([5])) * (4.90120e-01*(((x-[6])/[7])*1.16366e+00+TMath::Exp(-((x-[6])/[7])*1.16366e+00)-1) - 3.99384e+00*(((x-[6])/[7])*1.27500e-01+TMath::Exp(-((x-[6])/[7])*1.27500e-01)-1)))");
    flgnum3->SetParameters(3.23646e+04, 3.57952e-01, 1.06771e-03, 2.26404e+00, 7.93093e-01, 5.56454e-02, 7.50501e+00, 1.41060e+00, 6.29173e+00);

    if (std::atof(gROOT->GetVersion()) < 6.00) return;

	//---------------//
	//  User Define  //
	//---------------//
	gROOT->LoadMacro("$AMSProjLibs/CPPLibs/CPPLibs.h");
	gROOT->LoadMacro("$AMSProjLibs/ROOTLibs/ROOTLibs.h");
	gROOT->LoadMacro("$AMSProjLibs/TRACKLibs/TRACKLibs.h");
	gROOT->ProcessLine("MGROOT::LoadDefaultEnvironment();");
	gROOT->ProcessLine("using namespace MGROOT;");
	gROOT->ProcessLine("using namespace TrackSys;");
}
