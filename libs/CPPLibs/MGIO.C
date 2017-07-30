#ifndef __CPPLibs_MGIO_C__
#define __CPPLibs_MGIO_C__

#include "MGIO.h"

namespace MGIO {

std::vector<std::string> ReadFileContent(const std::string& path) {
	std::vector<std::string> list;
	if (!MGSys::TestFile(path, 'f')) return list;
	File file(path);
	if (!file.exist()) return list;
	file().seekg(0);
	for (std::string line; std::getline(file(), line); )
		list.push_back(line);
	return list;
}


bool Exist(const std::string& filename) {
	if (!MGSys::TestFile(filename, 'f')) return false;
	std::fstream file(filename, File::kRead);
	if (file.is_open()) { file.close(); return true; }
	else { return false; }
}


void File::open(const std::string& filename, OpenMode mode) {
	close();
	//bool isopt = ((mode&File::ReadWrite) != OpenMode(0));
	//if (!isopt) { MGSys::ShowError(LocAddr(), STR_FMT("File::OpenMode is wrong. (%d)", mode)); return; }
	//bool isfine = (isopt && (Exist(filename) || (mode&(File::Trunc|File::App)) != OpenMode(0)));
	//if (isfine) fFstr.open(filename, mode);
	fstr_.open(filename, mode);
	if (!fstr_.is_open()) { MGSys::ShowError(LOC_ADDR(), STR_FMT("Can't open file. (%s)", filename.c_str())); return; }
	mode_ = mode;
}

} // namespace MGIO

#endif // __CPPLibs_MGIO_C__
