#ifndef __CPPLibs_MGIO_H__
#define __CPPLibs_MGIO_H__

#include <fstream>
#include <cstdio>


namespace MGIO {

class File {
	public :
		typedef std::ios_base::openmode OpenMode;
		static constexpr OpenMode kReadWrite = std::ios::in|std::ios::out; // open mode
		static constexpr OpenMode kRead      = std::ios::in;               // open mode
		static constexpr OpenMode kWrite     = std::ios::out;              // open mode
		static constexpr OpenMode kTrunc     = std::ios::trunc;            // destroy contents
		static constexpr OpenMode kApp       = std::ios::app;              // append to file
		static constexpr OpenMode kBinary    = std::ios::binary;           // binary file

	public :
		File(const std::string& filename = std::string(), OpenMode mode = File::kRead) { open(filename, mode); }
		~File() { close(); }

		void open(const std::string& filename, OpenMode mode = File::kReadWrite);
	
		inline bool exist() { return fstr_.is_open(); }
		
		inline void close() { if (exist()) fstr_.close(); mode_ = File::kReadWrite; }

		inline std::fstream& operator()() { return fstr_; }

	protected :
		std::fstream   fstr_;
		File::OpenMode mode_;
};
    
std::vector<std::string> ReadFileContent(const std::string& path);
	
bool Exist(const std::string& filename);

} // namespace MGIO




#endif // __CPPLibs_MGIO_H__
