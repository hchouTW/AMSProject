#ifndef __MGIO_H__
#define __MGIO_H__



#include <fstream>
#include <cstdio>
#include "MGSys.h"


namespace MGIO {
	class File;
	std::vector<std::string> ReadFileContent(const std::string& path);
		
	bool Exist(const std::string& filename);
}


class MGIO::File {
	public :
		typedef std::ios_base::openmode OpenMode;
		static constexpr OpenMode ReadWrite = std::ios::in|std::ios::out; // open mode
		static constexpr OpenMode Read      = std::ios::in;               // open mode
		static constexpr OpenMode Write     = std::ios::out;              // open mode
		static constexpr OpenMode Trunc     = std::ios::trunc;            // destroy contents
		static constexpr OpenMode App       = std::ios::app;              // append to file
		static constexpr OpenMode Binary    = std::ios::binary;           // binary file

	public :
		File(const std::string& filename = std::string(), OpenMode mode = File::ReadWrite) { open(filename, mode); }
		~File() { close(); }

		void open(const std::string& filename, OpenMode mode = File::ReadWrite);
	
		inline bool exist() { return fFstr.is_open(); }
		
		inline void close() { if (exist()) fFstr.close(); fMode = File::ReadWrite; }

		inline std::fstream& operator()() { return fFstr; }

	protected :
		std::fstream   fFstr;
		File::OpenMode fMode;
};


#endif // __MGIO_H__
