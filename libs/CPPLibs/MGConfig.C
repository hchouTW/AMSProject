#ifndef __CPPLibs_MGConfig_C__
#define __CPPLibs_MGConfig_C__

#include "MGConfig.h"

#include <unistd.h>

namespace MGConfig {


JobOpt::JobOpt(int argc, char* argv[]) : succ_(false), type_("NONE"), ipath_("flist"), opath_("."), gi_(0), gs_(0) {
    int ch = 0;
    char* buf[5] = { nullptr };

    bool usht = false;
    const char* opt = "t:i:o:g:s:";
    while ((ch = getopt(argc, argv, opt)) != -1) {
        usht = true;
        bool def = true;
        switch (ch) {
            case 't':
                buf[0] = optarg;
                break;
            case 'i':
                buf[1] = optarg;
                break;
            case 'o':
                buf[2] = optarg;
                break;
            case 'g':
                buf[3] = optarg;
                break;
            case 's':
                buf[4] = optarg;
                break;
            default:
                def = false;
        }
        if (!def) MGSys::ShowErrorAndExit("JobOpt: invalid option! Exiting ...");
    }

    std::string type = (usht ? (buf[0] ? buf[0] : type_) : (argc>1 ? argv[1] : type_));
    std::use_facet<std::ctype<char> >(std::locale()).toupper(&type[0], &type[0] + type.size());
    if (type != "ISS" && type != "BT" && type != "MC" && type != "NONE")
        MGSys::ShowErrorAndExit("JobOpt: Can't find type(ISS, BT, MC, NONE)! Exiting ...");

    std::string opath = (usht ? (buf[2] ? buf[2] : opath_) : (argc>5 ? argv[5] : opath_));
	if (!MGSys::TestFile(opath, 'd'))
        MGSys::ShowErrorAndExit("JobOpt: Can't find opath! Exiting ...");

    std::string ipath = (usht ? (buf[1] ? buf[1] : ipath_) : (argc>2 ? argv[2] : ipath_));
    std::vector<std::string>&& list = MGIO::ReadFileContent(ipath);
    if (list.size() == 0)
        MGSys::ShowErrorAndExit("JobOpt: list empty ! Exiting ...");

    long gi = (usht ? (buf[3] ? std::atol(buf[3]) : gi_) : (argc>3 ? std::atol(argv[3]) : gi_));
    long gs = (usht ? (buf[4] ? std::atol(buf[4]) : gs_) : (argc>4 ? std::atol(argv[4]) : gs_));
    if (gs == 0) gs = list.size();
    if (gi < 0 || gs < 0 || gs > list.size() || (gs*gi) >= list.size())
        MGSys::ShowErrorAndExit("JobOpt: outside (gi, gs)! Exiting ...");

    long beg = gs * gi;
    long end = gs * (gi + 1);
    if (end > list.size()) end = list.size();
    std::vector<std::string> flist(list.begin()+beg, list.begin()+end);
    
    if (flist.size() == 0)
        MGSys::ShowErrorAndExit("JobOpt: flist empty ! Exiting ...");

    succ_  = true;
    type_  = type;
    ipath_ = ipath;
    opath_ = opath;
    gi_    = gi;
    gs_    = gs;
    flist_ = flist;
}


} // namesapce MGConfig


#endif // __CPPLibs_MGConfig_C__
