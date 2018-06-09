{
	gInterpreter->AddIncludePath("${AMSSRC}/include");
	gSystem->Load("${AMSLIB}/libntuple_slc6_PG_dynamic.so");

    Printf("Current Version    = %s", gROOT->GetVersion());
	Printf("Build Architecture = %s", gSystem->GetBuildArch());
    
    TF1 * f1gs = new TF1("f1gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1])");
    f1gs->SetParameters(10., 10.);
    f1gs->SetNpx(100000);

    TF1 * f2gs = new TF1("f2gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3])");
    f2gs->SetParameters(10., 10., 10., 30.);
    f2gs->SetNpx(100000);

    TF1 * f3gs = new TF1("f3gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5])");
    f3gs->SetParameters(10., 10., 10., 30., 10., 50.);
    f3gs->SetNpx(100000);
    
    TF1 * f4gs = new TF1("f4gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5]) + ([6]/[7])*TMath::Exp(-0.5*x*x/[7]/[7])");
    f4gs->SetParameters(10., 10., 10., 30., 10., 50., 10., 70.);
    f4gs->SetNpx(100000);
    
    TF1 * f5gs = new TF1("f5gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5]) + ([6]/[7])*TMath::Exp(-0.5*x*x/[7]/[7]) + ([8]/[9])*TMath::Exp(-0.5*x*x/[9]/[9])");
    f5gs->SetParameters(10., 10., 10., 30., 10., 50., 10., 70., 10., 100);
    f5gs->SetNpx(100000);
    
    TF1 * f1ags = new TF1("f1ags", "([0]/[1])*(TMath::Exp(-0.5*(x-[2])*(x-[2])/[1]/[1])+TMath::Exp(-0.5*(x+[2])*(x+[2])/[1]/[1]))");
    f1ags->SetParameters(10., 10., 0.);
    f1ags->SetNpx(100000);
    
    //TF1 * f1nags = new TF1("f1nags", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + (0.5*[2]/[3])*(TMath::Exp(-0.5*(x-[4])*(x-[4])/[3]/[3])+TMath::Exp(-0.5*(x+[4])*(x+[4])/[3]/[3]))");
    //f1nags->SetParameters(10., 10., 10., 30., 0.);
    //f1nags->SetNpx(100000);
    
    TF1 * f1nags = new TF1("f1nags", "([0]/[3])*([1]*TMath::Exp(-0.5*x*x/[3]/[3]) + TMath::Sqrt(2)*(TMath::Exp(-(x-[2])*(x-[2])/[3]/[3]) + TMath::Exp(-(x+[2])*(x+[2])/[3]/[3])))");
    f1nags->SetParameters(10., 1., 0., 10.);
    f1nags->SetNpx(100000);
    
    TF1 * feloss = new TF1("feloss", "[0] * TMath::Power( ([2]/x)/[1]/[1], ([2]/x)/[1]/[1] ) / TMath::Gamma( ([2]/x)/[1]/[1] ) * TMath::Exp(-(([2]/x)/[1]/[1]) * ((x-[2])/[3] + TMath::Exp(-(x-[2])/[3])) )");
    feloss->SetParameters(1000., 1.0, 0.0015, 0.0002); 
    
    TF1 * fseloss = new TF1("fseloss", "[0] * TMath::Power( (2*[2]/(x+[2]))/[4]/[4], (2*[2]/(x+[2]))/[4]/[4] ) / TMath::Gamma( (2*[2]/(x+[2]))/[4]/[4] ) * TMath::Exp(-( (([2]/x)/[1]/[1]) * ((x-[2])/[3] + TMath::Exp(-(x-[2])/[3])) ))");
    fseloss->SetParameters(1000., 1.0, 0.15, 0.02, 1.0);
    fseloss->FixParameter(4, 1.0);
   
    TF1 * fvavilov = new TF1("fvavilov", "[0] * TMath::Vavilov((x-[1])/[2], [3], [4])");
    fvavilov->SetParameters(1000, 0.15, 0.02, 1.0, 1.0);

    //const char* tt = "(1.0/[1]/[1]) * (([2]/[3]) / ((abs((x-[2])/[3])) + ([2]/[3])) )^[4]";
    //TF1 * feloss2 = new TF1("feloss2", Form("[0] * TMath::Power(%s, %s) / TMath::Gamma(%s) * TMath::Exp(-(%s) * ((x-[2])/[3] + TMath::Exp(-(x-[2])/[3])) )", tt, tt, tt, tt));
    //feloss2->SetParameters(1000., 1.0, 0.0015, 0.0002, 1.0); 
    
    TF1* fmpv = new TF1("fmpv", "[0] * (1+x*x)^[2] * ([1] - (1+x*x)^(-[2]) - TMath::Log([3]+x^[4]))");
    fmpv->SetParameters(10, 6.5, 1.0, 10.0, 1.0);
    
    TF1* fmpv2 = new TF1("fmpv2", "[0] * (1+x*x)^[3] * ([1] - [2]*(1+x*x)^(-[3]) - TMath::Log([4]+x^[5]))");
    fmpv2->SetParameters(10, 6.5, 1.0, 1.0, 10.0, 1.0);
    
    TF1* fkpa = new TF1("fkpa", "1.0 - 0.5*TMath::Erfc([0]*TMath::Log(1.0+[1]*x)+[2])");
    fkpa->SetParameters(1.0, 0.3, 0.0);
    
    TF1* fkpa2 = new TF1("fkpa2", "(1 - [0]) + [0] * 0.5 * (1.0 + TMath::Erf([1] * TMath::Log(1+[3]*(x*x)) + [2]))");
    fkpa2->SetParameters(1.0, 1.0, 0.3, 1.0);
    
    TF1* flg = new TF1("flg", "[0] * TMath::Exp( (1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/TMath::Landau(0)) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3]) )");
    flg->SetParameters(1.0, 0.1, 0.0, 1.0);
    flg->SetParLimits(1, 0.0, 1.0);
   
    TF1* fgm = new TF1("fgm", "[0] * (TMath::Erf((x - [1]) / [2]) + 1) * TMath::Power(x, [3]) * TMath::Exp(-[4] * x)");
    fgm->SetParameters(1.0, 6.0, 1.5, 3.0, 0.3);

    //TF1* flggm = new TF1("flggm", "[0] * TMath::Exp((1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/TMath::Landau(0)) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3])) * (TMath::Erf((x-[8])/[7])+1) + [4] * TMath::Power(x,[5]) * TMath::Exp(-[6]*x) * (TMath::Erf((x-[9])/[7])+1)");
    //flggm->SetParameters(2.94676e+04, 6.26581e-04, 1.89170e+00, 4.76106e-01, 1.20969e+02, 2.61678e+00, 2.87316e-01, 1.65832e+00, 2.37444e+00, 6.19537e+00);
    //flggm->SetParLimits(1, 0.0, 1.0);
    //flggm->SetNpx(10000);
    
    TF1* flggm = new TF1("flggm", "[0] * TMath::Exp((1-[1]) * TMath::Log(TMath::Landau((x-[2])/[3])/TMath::Landau(0)) + [1] * (-0.5)*((x-[2])*(x-[2])/[3]/[3])) * (TMath::Erf((x-[4])/[5])+1) + [6] * TMath::Power(x,[7]) * TMath::Exp(-[8]*x) * (TMath::Erf((x-[9])/[10])+1)");
    flggm->SetParameters(1.86244e+04, 8.01663e-03, 2.46953e+00, 6.79641e-01, 8.69772e-01, 2.90538e-01, 1.10869e+02, 3.18041e+00, 4.10684e-01, 6.50848e+00, 1.86977e+00);
    flggm->SetParLimits(1, 0.0, 1.0);
    flggm->SetNpx(10000);
    
    
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
