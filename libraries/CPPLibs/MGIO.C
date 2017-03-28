#ifndef __MGIO_C__
#define __MGIO_C__

#include "MGIO.h"

std::vector<std::string> MGIO::ReadFileContent(const std::string& path) {
	std::vector<std::string> list;
	if (!MGSys::TestFile(path, 'f')) return list;
	MGIO::File file(path);
	if (!file.exist()) return list;
	file().seekg(0);
	for (std::string line; std::getline(file(), line); )
		list.push_back(line);
	return list;
}


bool MGIO::Exist(const std::string& filename) {
	if (!MGSys::TestFile(filename, 'f')) return false;
	std::fstream file(filename, File::ReadWrite);
	if (file.is_open()) { file.close(); return true; }
	else { return false; }
}


void MGIO::File::open(const std::string& filename, OpenMode mode) {
	close();
	bool isopt = ((mode&File::ReadWrite) != OpenMode(0));
	if (!isopt) { MGSys::ShowError(LocAddr(), StrFmt("File::OpenMode is wrong. (%d)", mode)); return; }
	bool isfine = (isopt && (MGIO::Exist(filename) || (mode&(File::Trunc|File::App)) != OpenMode(0)));
	if (isfine) fFstr.open(filename, mode);
	if (!fFstr.is_open()) { MGSys::ShowError(LocAddr(), StrFmt("Can't open file. (%s)", filename.c_str())); return; }
	fMode = mode;
}


#endif // __MGIO_C__
