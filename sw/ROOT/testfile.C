// Usage :
//========================================================
// match=`root -b -q 'testfile.C("YiNtuple_ISS.0080457.root", "data")' | grep TESTFILE_IS_SUCCESS | wc -l`
// if (( $match == 1 )); then
//   echo "SUCCESS"
// else
//   echo "FAILURE"
// fi
//========================================================

#include <string>
#include "TFile.h"
#include "TTree.h"

bool testfile(const char* path = "", const char* name = "data") {
    std::string pathstr = path;
    std::string namestr = name;
    if (pathstr.size() == 0) return false;
    if (namestr.size() == 0) return false;

    TFile file(pathstr.c_str());
    if (file.IsZombie() || file.TestBit(TFile::kRecovered)) return false;

    TTree* tree = file.Get(namestr.c_str());
    if (tree == 0) return false;

    printf("TESTFILE_IS_SUCCESS\n");
    return true;
}
