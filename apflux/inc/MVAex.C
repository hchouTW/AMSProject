#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int main() {
    TMVA::Tools::Instance();

    std::map<std::string,int> Use;

    Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)

#ifdef R__HAS_TMVAGPU
    Use["DNN_GPU"]         = 0; // CUDA-accelerated DNN training.
#else
    Use["DNN_GPU"]         = 0;
#endif

#ifdef R__HAS_TMVACPU
    Use["DNN_CPU"]         = 0; // Multi-core accelerated DNN.
#else
    Use["DNN_CPU"]         = 0;
#endif

    Use["BDTG"] = 1; // uses Gradient Boost
    Use["BDTD"] = 1; // decorrelation + Adaptive Boost

    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;
    

    //TFile* input = nullptr;
    //TString fname = "./tmva_class_example.root";
    //if (!gSystem->AccessPathName( fname )) {
    //    input = TFile::Open( fname ); // check if file in local directory exists
    //}
    //else {
    //    TFile::SetCacheFileDir(".");
    //    input = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD");
    //}
    //if (!input) {
    //    std::cout << "ERROR: could not open data file" << std::endl;
    //    exit(1);
    //}
   
    TString inputS = "/eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/dataset.HEXp.iss31.root";
    TString inputB = "/eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/dataset.HEXn.mcpr_l1o9flux31.root";

    std::cout << "--- TMVAClassification       : Using inputS files: " << inputS.Data() << std::endl;
    std::cout << "--- TMVAClassification       : Using inputB files: " << inputB.Data() << std::endl;
    
    TChain *signalTree = new TChain("varsHEXp");
    TChain *background = new TChain("varsHEXn");

    signalTree->Add(inputS.Data());
    background->Add(inputB.Data());

    //TString outfileName( "mva/MVAex/TMVA.root" );
    TString outfileName( "mva/test_MVAex/TMVA.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                                "!V:!Silent:Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

    //TMVA::DataLoader *dataloader=new TMVA::DataLoader("mva/MVAex");
    TMVA::DataLoader *dataloader=new TMVA::DataLoader("mva/test_MVAex");
    
    dataloader->AddVariable( "lxin" , "lxin" , "units", 'F' );
    dataloader->AddVariable( "lyin" , "lyin" , "units", 'F' );
    dataloader->AddVariable( "lxex" , "lxex" , "units", 'F' );
    dataloader->AddVariable( "lyex" , "lyex" , "units", 'F' );
    
    dataloader->AddVariable( "nlxin", "nlxin", "units", 'F' );
    dataloader->AddVariable( "nlyin", "nlyin", "units", 'F' );
    dataloader->AddVariable( "nlxex", "nlxex", "units", 'F' );
    dataloader->AddVariable( "nlyex", "nlyex", "units", 'F' );
    
    dataloader->AddVariable( "lxinC", "lxinC", "units", 'F' );
    dataloader->AddVariable( "lyinC", "lyinC", "units", 'F' );
    dataloader->AddVariable( "lxexC", "lxexC", "units", 'F' );
    dataloader->AddVariable( "lyexC", "lyexC", "units", 'F' );

    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
    dataloader->AddSignalTree    ( signalTree,     signalWeight );
    dataloader->AddBackgroundTree( background, backgroundWeight );
    //dataloader->SetSignalWeightExpression( "(wgt * mc_w27)" );
    //dataloader->SetBackgroundWeightExpression( "" );

    TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                            "nTrain_Signal=200000:"
                                            "nTrain_Background=50000:"
                                            "nTest_Signal=100000:"
                                            "nTest_Background=25000:"
                                            "SplitMode=Random:NormMode=NumEvents:!V" );

    if (Use["LikelihoodPCA"])
           factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
                                "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );

    if (Use["DNN_CPU"] or Use["DNN_GPU"]) {
        // General layout.
        TString layoutString ("Layout=TANH|16,TANH|16,TANH|16,LINEAR");

        // Training strategies.
        TString training0("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                          "ConvergenceSteps=30,BatchSize=256,TestRepetitions=10,"
                          "WeightDecay=1e-4,Regularization=None,"
                          "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
        TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                          "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                          "WeightDecay=1e-4,Regularization=L2,"
                          "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
        TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                          "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                          "WeightDecay=1e-4,Regularization=L2,"
                          "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
        TString trainingStrategyString ("TrainingStrategy=");
        trainingStrategyString += training0 + "|" + training1 + "|" + training2;

        // General Options.
        TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                            "WeightInitialization=XAVIERUNIFORM");
        dnnOptions.Append (":"); dnnOptions.Append (layoutString);
        dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

        // Cuda implementation.
        if (Use["DNN_GPU"]) {
            TString gpuOptions = dnnOptions + ":Architecture=GPU";
            factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_GPU", gpuOptions);
        }
        // Multi-core CPU implementation.
        if (Use["DNN_CPU"]) {
            TString cpuOptions = dnnOptions + ":Architecture=CPU";
            factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", cpuOptions);
        }
    }        

    if (Use["BDTD"]) // Decorrelation + Adaptive Boost
           factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
                               "V=True:"
                               "VerbosityLevel=Default:"
                               "VarTransform=Decorrelate:"
                               "H=True:"
                               "CreateMVAPdfs=False:"
                               "IgnoreNegWeightsInTraining=True:"
                               "NTrees=100:"
                               "MaxDepth=3:"
                               "MinNodeSize=5%:"
                               "nCuts=20:"
                               "BoostType=AdaBoost:"
                               "AdaBoostR2Loss=Quadratic:"
                               "UseBaggedGrad=False:"
                               "Shrinkage=0.1:"
                               "AdaBoostBeta=0.5:"
                               "UseRandomisedTrees=False:"
                               "UseNvars=2:"
                               "UsePoissonNvars=True:"
                               "BaggedSampleFraction=0.5:"
                               "UseYesNoLeaf=True:"
                               "NodePurityLimit=0.5:"
                               "SeparationType=GiniIndex:"
                               "DoBoostMonitor=False:"
                               "UseFisherCuts=False:"
                               "MinLinCorrForFisher=0.8:"
                               "UseExclusiveVars=False:"
                               "DoPreselection=False:"
                               "SigToBkgFraction=1:"
                               "PruneMethod=NoPruning:"
                               "PruneStrength=0:"
                               "PruningValFraction=0.5:"
                               "nEventsMin=0:"
                               "GradBaggingFraction=0.6:"
                               "UseNTrainEvents=0:"
                               "NNodesMax=0");

    if (Use["BDTG"]) // Gradient Boost
           factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                               "V=True:"
                               "VerbosityLevel=Default:"
                               "VarTransform=Decorrelate:"
                               "H=True:"
                               "CreateMVAPdfs=False:"
                               "IgnoreNegWeightsInTraining=True:"
                               "NTrees=100:"
                               "MaxDepth=3:"
                               "MinNodeSize=5%:"
                               "nCuts=20:"
                               "BoostType=Grad:"
                               "AdaBoostR2Loss=Quadratic:"
                               "UseBaggedGrad=False:"
                               "Shrinkage=0.1:"
                               "AdaBoostBeta=0.5:"
                               "UseRandomisedTrees=False:"
                               "UseNvars=2:"
                               "UsePoissonNvars=True:"
                               "BaggedSampleFraction=0.5:"
                               "UseYesNoLeaf=True:"
                               "NodePurityLimit=0.5:"
                               "SeparationType=GiniIndex:"
                               "DoBoostMonitor=False:"
                               "UseFisherCuts=False:"
                               "MinLinCorrForFisher=0.8:"
                               "UseExclusiveVars=False:"
                               "DoPreselection=False:"
                               "SigToBkgFraction=1:"
                               "PruneMethod=NoPruning:"
                               "PruneStrength=0:"
                               "PruningValFraction=0.5:"
                               "nEventsMin=0:"
                               "GradBaggingFraction=0.6:"
                               "UseNTrainEvents=0:"
                               "NNodesMax=0");

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
    delete dataloader;

    if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

    return 0;
}
