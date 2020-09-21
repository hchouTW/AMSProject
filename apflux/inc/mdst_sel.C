/**********************************************************
 * Author        : Hsin-Yi Chou
 * Email         : hchou@cern.ch
 * Last modified : 2015-07-22 16:29
 * Filename      : YiProdNtuple.C
 * Description   :
 * *******************************************************/
#ifndef __MDst_C__
#define __MDst_C__

#include "selector.h"
#include "selector.C"

// Args
DEFINE_string(type, "ISS", "Event Type: ISS, BT, MC");
DEFINE_string(inpath, "lst/flist.cern.iss.B1130.pass7", "Input Path");
DEFINE_string(outpath, "out", "Output Path");
DEFINE_uint64(gindex, 0, "Index of Group");
DEFINE_uint64(gsize, 1, "Size of Group");

static bool ValidateType(const char* flagname, const std::string& type) {
    std::string type_toupper = type;
    std::transform(type_toupper.begin(), type_toupper.end(), type_toupper.begin(), [](unsigned char c){ return std::toupper(c); });
    if (type_toupper == "ISS" || type_toupper == "BT" || type_toupper == "MC") return true;
    return false;
}

static bool ValidateInpath(const char* flagname, const std::string& inpath) {
    return (std::system((Format("test -f \"%s\"", inpath.c_str())).c_str()) == 0);
}

static bool ValidateOutpath(const char* flagname, const std::string& outpath) {
    return (std::system((Format("test -d \"%s\"", outpath.c_str())).c_str()) == 0);
}

//static const bool type_dummy = gflags::RegisterFlagValidator(&FLAGS_type, &ValidateType);
//static const bool type_inpath = gflags::RegisterFlagValidator(&FLAGS_inpath, &ValidateInpath);
//static const bool type_outpath = gflags::RegisterFlagValidator(&FLAGS_outpath, &ValidateOutpath);
static const bool type_dummy = google::RegisterFlagValidator(&FLAGS_type, &ValidateType);
static const bool type_inpath = google::RegisterFlagValidator(&FLAGS_inpath, &ValidateInpath);
static const bool type_outpath = google::RegisterFlagValidator(&FLAGS_outpath, &ValidateOutpath);
    

int main(int argc, char** argv) {
    std::string usage = Format("This program is used to analysis ntuple.");
    usage += Format("\nUsage: bin/mdst_sel -type=ISS -inpath=lst/flist.cern.iss.B1130.pass7 -outpath=out -gindex=0 -gsize=1\n");

    //gflags::SetVersionString("0.0.1");
    //gflags::SetUsageMessage(usage);
    //gflags::ParseCommandLineFlags(&argc, &argv, true);
    //google::SetVersionString("0.0.1");
    google::SetUsageMessage(usage);
    google::ParseCommandLineFlags(&argc, &argv, true);

    FLAGS_log_dir = "log";
    google::InitGoogleLogging(argv[0]);

    // Type
    std::string type = FLAGS_type;
    std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c){ return std::toupper(c); });

    // List
    std::vector<std::string>&& list = ReadListFile(FLAGS_inpath);
    LOG_IF(ERROR, FLAGS_gsize > list.size()) << Format("Size of Group (FLAGS_gsize %ld, LIST_size %ld)", FLAGS_gsize, list.size());
   
    unsigned long gindex = FLAGS_gindex;
    unsigned long gsize = (FLAGS_gsize == 0) ? list.size() : FLAGS_gsize;
    unsigned long gbeg  = gindex * gsize;
    unsigned long gend  = (gindex + 1) * gsize;
    if (gend > list.size()) gend = list.size();

    LOG_IF(ERROR, gbeg >= list.size()) << Format("Group (%ld, %ld) LIST_size (%ld)", gbeg, gend, list.size());

    std::vector<std::string> sublist = std::vector<std::string>(list.cbegin() + gbeg, list.cbegin() + gend);
    LOG_IF(ERROR, sublist.size() == 0) << Format("Size of List %ld", sublist.size());

    // Output
    std::string outpath = Format("%s/YiMdst.%07ld.root", FLAGS_outpath.c_str(), gindex);
   
    std::string args_summary;
    args_summary += Format("\n**----------------------    mdst args summary ----------------------**\n");
    args_summary += Format("Type  : %s\n", type.c_str());
    args_summary += Format("Input : %s\n", FLAGS_inpath.c_str());
    args_summary += Format("Output: %s\n", outpath.c_str());
    args_summary += Format("Index of Group: %ld\n", gindex);
    args_summary += Format("Size  of Group: %ld\n", gsize);
    args_summary += Format("List of ROOT Files: (size %ld)\n", (gend - gbeg));
    for (auto&& file : sublist) args_summary += Format("FILE: %s\n", file.c_str());
    args_summary += Format("**------------------------------------------------------------------**\n\n");
   
    std::cout << args_summary;
    LOG(INFO) << args_summary;

    std::string statement_start;
    statement_start += Format("\n**-------------------------**\n");
    statement_start += Format("**    ana ntuple START     **\n");
    statement_start += Format("**-------------------------**\n\n");
    std::cout << statement_start;
    LOG(INFO) << statement_start;

    TChain chain("mdstZ1");
    for (auto&& file : sublist) chain.Add(file.c_str());

	std::string chain_statement;
    chain_statement += Format("\nROOT files\n");
	chain_statement += Format("Totally Events : %ld\n\n", chain.GetEntries());

    std::cout << chain_statement;
    LOG(INFO) << chain_statement;

    if (chain.GetEntries() == 0) std::cerr << Format("Don't have event\n");
    LOG_IF(ERROR, chain.GetEntries() == 0) << Format("Don't have event");

    if      (type == "ISS") Selector::SetType(Selector::Type::ISS);
    else if (type == "BT" ) Selector::SetType(Selector::Type::BT);
    else if (type == "MC" ) Selector::SetType(Selector::Type::MC);

    Selector selector(&chain);
    selector.set_output(outpath);
    selector.set_environment();
    selector.process_events();
    selector.write();
    selector.close();
    
    std::string statement_end;
    statement_end += Format("\n**-------------------------**\n");
    statement_end += Format("**    ana ntuple END       **\n");
    statement_end += Format("**-------------------------**\n\n");
    std::cout << statement_end;
    LOG(INFO) << statement_end;
    
    google::ShutdownGoogleLogging();
    //gflags::ShutDownCommandLineFlags();
    //google::ShutDownCommandLineFlags();
	return 0;
}
#endif // ___MDst_C__
