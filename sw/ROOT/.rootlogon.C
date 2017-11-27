{
	gInterpreter->AddIncludePath("${AMSSRC}/include");
	gSystem->Load("${AMSLIB}/libntuple_slc6_PG_dynamic.so");

	std::cout << Form("Current Version    = %s\n", gROOT->GetVersion());
	std::cout << Form("Build Architecture = %s\n", gSystem->GetBuildArch());

    TF1 * f3gs = new TF1("f3gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5])");
    f3gs->SetParameters(10., 10., 10., 30., 10., 50.);
    
    TF1 * f4gs = new TF1("f4gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5]) + ([6]/[7])*TMath::Exp(-0.5*x*x/[7]/[7])");
    f4gs->SetParameters(10., 10., 10., 30., 10., 50., 10., 70.);
    
    TF1 * f5gs = new TF1("f5gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5]) + ([6]/[7])*TMath::Exp(-0.5*x*x/[7]/[7]) + ([8]/[9])*TMath::Exp(-0.5*x*x/[9]/[9])");
    f5gs->SetParameters(10., 10., 10., 30., 10., 50., 10., 70., 10., 100);
    
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
