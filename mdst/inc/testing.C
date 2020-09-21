#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Size of argc is not 3." << std::endl;
        return 0;
    }

    std::string opt = argv[1];
    if (opt != "check" && opt != "kill") return 0;
    bool opt_check = (opt == "check");
    bool opt_kill  = (opt == "kill");
  
    bool status = false;
    std::string msg;

    std::string path = argv[2];
    TFile* file = TFile::Open(path.c_str());
    if (file != nullptr && !file->IsZombie()) {
        TTree* tree = (TTree*)file->Get("mdst");
        //TTree* treeZ1 = (TTree*)file->Get("mdstZ1");
        //TTree* treeZ2 = (TTree*)file->Get("mdstZ2");
        if (tree != nullptr && !tree->IsZombie()) {
        //if (treeZ1 != nullptr && treeZ2 != nullptr && !treeZ1->IsZombie() && !treeZ2->IsZombie()) {
            status = true;
            msg = "SUCCESS";
        }
        else {
            msg = "failure( TTree )";
        }
    }
    else {
        msg = "failure( TFile )";
    }
   
    if (opt_kill && file != nullptr && !status) {
        std::string command = Form("/bin/rm %s", path.c_str());
        std::cout << command << std::endl;
        system(command.c_str());
    }

    std::cout << msg << std::endl;
    return status;
}
