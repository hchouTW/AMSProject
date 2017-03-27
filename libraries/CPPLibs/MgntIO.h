#ifndef __MgntIO_H__
#define __MgntIO_H__
#include <fstream>
#include <cstdio>
#include "MgntSys.h"

namespace MgntIO {
	class File;
	std::vector<std::string> ReadFileContent(const std::string& path);
}

class MgntIO::File {
	public :
		typedef std::ios_base::openmode OpenMode;
		static constexpr OpenMode ReadWrite = std::ios::in|std::ios::out; // open mode
		static constexpr OpenMode Read      = std::ios::in;               // open mode
		static constexpr OpenMode Write     = std::ios::out;              // open mode
		static constexpr OpenMode Trunc     = std::ios::trunc;            // destroy contents
		static constexpr OpenMode App       = std::ios::app;              // append to file
		static constexpr OpenMode Binary    = std::ios::binary;           // binary file
		
		inline static bool Exist(const std::string& filename) {
			if (!MgntSys::TestFile(filename, 'f')) return false;
			std::fstream file(filename, File::ReadWrite);
			if (file.is_open()) { file.close(); return true; }
			else { return false; }
		}

	public :
		File(const std::string& filename = std::string(), OpenMode mode = File::ReadWrite) { open(filename, mode); }
		~File() { close(); }

		void open(const std::string& filename, OpenMode mode = File::ReadWrite) {
			close();
			bool isopt = ((mode&File::ReadWrite) != OpenMode(0));
			if (MgntSys::Error(LocAddr(), !isopt, StrFmt("File::OpenMode is wrong. (%d)", mode))) return;
			bool isfine = (isopt && (File::Exist(filename) || (mode&(File::Trunc|File::App)) != OpenMode(0)));
			if (isfine) fFstr.open(filename, mode);
			if (MgntSys::Error(LocAddr(), !fFstr.is_open(), StrFmt("Can't open file. (%s)", filename.c_str()))) return;
			fMode = mode;
		}
	
		inline bool exist() { return fFstr.is_open(); }
		
		inline void close() { if (exist()) fFstr.close(); fMode = File::ReadWrite; }

		inline std::fstream& operator()() { return fFstr; }

	protected :
		std::fstream   fFstr;
		File::OpenMode fMode;
};

std::vector<std::string> MgntIO::ReadFileContent(const std::string& path) {
	std::vector<std::string> list;
	if (!MgntSys::TestFile(path, 'f')) return list;
	MgntIO::File file(path);
	if (!file.exist()) return list;
	file().seekg(0);
	for (std::string line; std::getline(file(), line); )
		list.push_back(line);
	return list;
}

#endif // __MgntIO_H__
