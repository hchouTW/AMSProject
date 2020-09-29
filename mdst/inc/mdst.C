/**********************************************************
 * Author        : Hsin-Yi Chou
 * Email         : hchou@cern.ch
 * Last modified : 2015-07-22 16:29
 * Filename      : YiProdNtuple.C
 * Description   :
 * *******************************************************/
#ifndef _MDst_C__
#define _MDst_C__

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
    return (std::system((Fmt("test -f \"%s\"", inpath.c_str())).c_str()) == 0);
}

static bool ValidateOutpath(const char* flagname, const std::string& outpath) {
    return (std::system((Fmt("test -d \"%s\"", outpath.c_str())).c_str()) == 0);
}

//static const bool type_dummy = gflags::RegisterFlagValidator(&FLAGS_type, &ValidateType);
//static const bool type_inpath = gflags::RegisterFlagValidator(&FLAGS_inpath, &ValidateInpath);
//static const bool type_outpath = gflags::RegisterFlagValidator(&FLAGS_outpath, &ValidateOutpath);
static const bool type_dummy = google::RegisterFlagValidator(&FLAGS_type, &ValidateType);
static const bool type_inpath = google::RegisterFlagValidator(&FLAGS_inpath, &ValidateInpath);
static const bool type_outpath = google::RegisterFlagValidator(&FLAGS_outpath, &ValidateOutpath);
    

int main(int argc, char** argv) {
    std::string usage = Fmt("This program is used to produce mdst ntuple.");
    usage += Fmt("\nUsage: bin/mdst -type=ISS -inpath=lst/flist.cern.iss.B1130.pass7 -outpath=out -gindex=0 -gsize=1\n");

    //gflags::SetVersionString("0.0.1");
    //gflags::SetUsageMessage(usage);
    //gflags::ParseCommandLineFlags(&argc, &argv, true);
    //google::SetVersionString("0.0.1");
    google::SetUsageMessage(usage);
    google::ParseCommandLineFlags(&argc, &argv, true);

    FLAGS_log_dir = "log";
    FLAGS_minloglevel = 3;
    FLAGS_stop_logging_if_full_disk = true;
    google::InitGoogleLogging(argv[0]);

    // Type
    std::string type = FLAGS_type;
    std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c){ return std::toupper(c); });

    // List
    std::vector<std::string>&& list = ReadListFile(FLAGS_inpath);
    LOG_IF(ERROR, FLAGS_gsize > list.size()) << Fmt("Size of Group (FLAGS_gsize %ld, LIST_size %ld)", FLAGS_gsize, list.size());
   
    unsigned long gindex = FLAGS_gindex;
    unsigned long gsize = (FLAGS_gsize == 0) ? list.size() : FLAGS_gsize;
    unsigned long gbeg  = gindex * gsize;
    unsigned long gend  = (gindex + 1) * gsize;
    if (gend > list.size()) gend = list.size();

    LOG_IF(ERROR, gbeg >= list.size()) << Fmt("Group (%ld, %ld) LIST_size (%ld)", gbeg, gend, list.size());

    std::vector<std::string> sublist = std::vector<std::string>(list.cbegin() + gbeg, list.cbegin() + gend);
    LOG_IF(ERROR, sublist.size() == 0) << Fmt("Size of List %ld", sublist.size());

    // Output
    std::string outpath = Fmt("%s/YiMdst.%07ld.root", FLAGS_outpath.c_str(), gindex);
   
    std::string args_summary;
    args_summary += Fmt("\n**----------------------    mdst args summary ----------------------**\n");
    args_summary += Fmt("Type  : %s\n", type.c_str());
    args_summary += Fmt("Input : %s\n", FLAGS_inpath.c_str());
    args_summary += Fmt("Output: %s\n", outpath.c_str());
    args_summary += Fmt("Index of Group: %ld\n", gindex);
    args_summary += Fmt("Size  of Group: %ld\n", gsize);
    args_summary += Fmt("List of ROOT Files: (size %ld)\n", (gend - gbeg));
    for (auto&& file : sublist) args_summary += Fmt("FILE: %s\n", file.c_str());
    args_summary += Fmt("**------------------------------------------------------------------**\n\n");
   
    std::cout << args_summary;
    LOG(INFO) << args_summary;

    std::string statement_start;
    statement_start += Fmt("\n**--------------------------**\n");
    statement_start += Fmt("**    mdst ntuple START     **\n");
    statement_start += Fmt("**--------------------------**\n\n");
    std::cout << statement_start;
    LOG(INFO) << statement_start;

    // Read root files
    bool stagedonly = true;
	unsigned int timeout = 10;
	AMSChain ams_chain("AMSRoot");
	int chain_status = ams_chain.AddFromFile(FLAGS_inpath.c_str(), gbeg, gend, stagedonly, timeout);
	if (chain_status == -1) std::cerr << Fmt("ROOT files cannot be opend!");
    LOG_IF(ERROR, chain_status == -1) << Fmt("ROOT files cannot be opend!");

	std::string chain_statement;
    chain_statement += Fmt("\nAMS ROOT files\n");
    chain_statement += Fmt("Status : %d\n", chain_status);
	chain_statement += Fmt("Totally Events : %ld\n\n", ams_chain.GetEntries());

    std::cout << chain_statement;
    LOG(INFO) << chain_statement;

    if (ams_chain.GetEntries() == 0) std::cerr << Fmt("Don't have event\n");
    LOG_IF(ERROR, ams_chain.GetEntries() == 0) << Fmt("Don't have event");

    if      (type == "ISS") Selector::SetType(Selector::Type::ISS);
    else if (type == "BT" ) Selector::SetType(Selector::Type::BT);
    else if (type == "MC" ) Selector::SetType(Selector::Type::MC);

    Selector selector(&ams_chain);
    selector.set_output(outpath);
    selector.set_environment();
    selector.process_events();
    selector.write();
    selector.close();

    std::string statement_end;
    statement_end += Fmt("\n**--------------------------**\n");
    statement_end += Fmt("**    mdst ntuple END       **\n");
    statement_end += Fmt("**--------------------------**\n\n");
    std::cout << statement_end;
    LOG(INFO) << statement_end;
    
    google::ShutdownGoogleLogging();
    //gflags::ShutDownCommandLineFlags();
    //google::ShutDownCommandLineFlags();
	return 0;
}
#endif // __MDst_C__
