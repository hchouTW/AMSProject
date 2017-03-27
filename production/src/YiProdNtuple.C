/**********************************************************
 * Author        : Hsin-Yi Chou
 * Email         : hchou@cern.ch
 * Last modified : 2015-07-22 16:29
 * Filename      : YiProdNtuple.C
 * Description   :
 * *******************************************************/
#ifndef __YiProdNtuple_C__
#define __YiProdNtuple_C__
#include "YiProdNtuple.h"
#include "YiProdNtuple.tcc"
#include <string.h>

int main(int argc, const char ** argv) {
	std::cout << "\n**--------------------------**\n";
	std::cout << "\n**    YiProdNtuple START    **\n";
	std::cout << "\n**--------------------------**\n";

	YiNtuple::setSelectionMode(YiNtuple::NORM);
	//YiNtuple::setSelectionMode(YiNtuple::COPY);

	DataSelection::setOption(DataSelection::LIST, DataSelection::ON);
	DataSelection::setOption(DataSelection::RTI, DataSelection::ON);
	DataSelection::setOption(DataSelection::TRG, DataSelection::ON);
	DataSelection::setOption(DataSelection::TOF, DataSelection::ON);
	DataSelection::setOption(DataSelection::ACC, DataSelection::ON);
	DataSelection::setOption(DataSelection::TRK, DataSelection::ON);
	DataSelection::setOption(DataSelection::TRD, DataSelection::ON);
	DataSelection::setOption(DataSelection::RICH, DataSelection::ON);
	DataSelection::setOption(DataSelection::ECAL, DataSelection::ON);

	EventBase::setEventVersion(EventBase::B950);

	std::cout << std::endl << std::endl;
	std::cout << "Usage : YiProdNtuple event_mode file_list group_id group_size (path)\n";
	std::cout << "    Parameters : \n";
	std::cout << "    event_mode [ISS BT MC]\n";
	std::cout << "    file_list\n";
	std::cout << "    group_id\n";
	std::cout << "    group_size\n";
	std::cout << "    (path)\n";
	std::cout << std::endl << std::endl;

	if (argc != 5 && argc != 6) {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("Number of argument is not conform! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}

	std::string event_mode = argv[1];
	std::string file_list = argv[2];
	Long64_t group_id = atol(argv[3]);
	Long64_t group_size = atol(argv[4]);

	std::use_facet<std::ctype<char> >(std::locale()).toupper(&event_mode[0], &event_mode[0] + event_mode.size());
	if (event_mode == "ISS") EventBase::setEventMode(EventBase::ISS);
	else if (event_mode == "BT") EventBase::setEventMode(EventBase::BT);
	else if (event_mode == "MC") EventBase::setEventMode(EventBase::MC);
	else {
		MgntSys::Error(LocAddr(), MgntSys::MESSAGE("Can't find event mode (ISS, BT, MC)! Exiting ..."));
		MgntSys::Exit(EXIT_FAILURE);
	}

	std::string outputFile = "";
	if (YiNtuple::checkSelectionMode(YiNtuple::NORM))
		outputFile = CStrFmt("YiNtuple_%s.%07ld.root", event_mode.c_str(), group_id);
	else if (YiNtuple::checkSelectionMode(YiNtuple::COPY))
		outputFile = CStrFmt("YiMirror_%s.%07ld.root", event_mode.c_str(), group_id);

	std::string path = ".";
	if (argc == 6) path = argv[5];

	bool isMultiTree = false;
	YiNtuple * ntuple = new YiNtuple();
	ntuple->setOutputFile(outputFile, path, isMultiTree);
	ntuple->readDataFrom(file_list, group_id, group_size);
	ntuple->loopEventChain();
	if (ntuple != 0) delete ntuple;
	ntuple = 0;

	std::cout << "\n**------------------------**\n";
	std::cout << "\n**    YiProdNtuple END    **\n";
	std::cout << "\n**------------------------**\n";
	return 0;
}
#endif // __YiProdNtuple_C__
