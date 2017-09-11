{
	gInterpreter->AddIncludePath("${AMSSRC}/include");
	gSystem->Load("${AMSLIB}/libntuple_slc6_PG_dynamic.so");

	std::cout << Form("Current Version    = %s\n", gROOT->GetVersion());
	std::cout << Form("Build Architecture = %s\n", gSystem->GetBuildArch());

    TF1 * f3gs = new TF1("f3gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5])");
    f3gs->SetParameters(10., 0.001, 10., 0.003, 10., 0.01);
    TF1 * f4gs = new TF1("f4gs", "([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5]) + ([6]/[7])*TMath::Exp(-0.5*x*x/[7]/[7])");
    f4gs->SetParameters(10., 0.001, 10., 0.003, 10., 0.01, 10., 0.03);
    
    TF1 * fld = new TF1("fld", "[0]*TMath::Exp(-0.5 * ( (x-[1])/[2] + TMath::Exp(-(x-[1])/[2]) ) )");
    fld->SetParameters(10., 0.01, 0.005);

    TF1 * f3gsp = new TF1("f3gsp", "TMath::Power(([0]/[1])*TMath::Exp(-0.5*x*x/[1]/[1]) + ([2]/[3])*TMath::Exp(-0.5*x*x/[3]/[3]) + ([4]/[5])*TMath::Exp(-0.5*x*x/[5]/[5]), [6])");
    f3gsp->SetParameters(10., 0.001, 10., 0.003, 10., 0.01, 1.0);
    
    //TF1 * fldgs = new TF1("fldgs", "[0]*TMath::Exp(-0.5 * ( (x-[1])/[2] + TMath::Exp(-(x-[1])/[2]) ) ) * TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2]/[3]/[3])");
    TF1 * fldgs = new TF1("fldgs", "[0]*TMath::Exp(-0.5 * ( (x-[1])/[2] + TMath::Exp(-(x-[1])/[2]*[4]) + ((x-[1])*(x-[1])/[2]/[2])*[3]*[3]) )");
    fldgs->SetParameters(1000., 0.001, 0.0002, 1.0, 1.0);
    
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
