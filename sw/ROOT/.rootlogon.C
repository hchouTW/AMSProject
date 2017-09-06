{
	gInterpreter->AddIncludePath("${AMSSRC}/include");
	gSystem->Load("${AMSLIB}/libntuple_slc6_PG_dynamic.so");

	std::cout << Form("Current Version    = %s\n", gROOT->GetVersion());
	std::cout << Form("Build Architecture = %s\n", gSystem->GetBuildArch());

    TF1 * f3gs = new TF1("f3gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5])");

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
