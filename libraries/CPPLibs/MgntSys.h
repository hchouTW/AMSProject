#ifndef __MgntSys_H__
#define __MgntSys_H__
#include <iostream>
#include <cstdlib>
#include <exception>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <dirent.h>


#define LocAddr() (MgntSys::StringFormat("LINE %d  {FILE  %s} [FUNC  %s] :  ", __LINE__, __FILE__, __func__))
#define StrFmt(fmt, ...)  (MgntSys::StringFormat(fmt, ##__VA_ARGS__))
#define CStrFmt(fmt, ...) ((MgntSys::StringFormat(fmt, ##__VA_ARGS__)).c_str())

namespace MgntSys {
	std::string StringFormat(const char * fmt, ...) {
		va_list args1; va_start(args1, fmt);
		va_list args2; va_copy(args2, args1);
		std::vector<char> buf(1+std::vsnprintf(nullptr, 0, fmt, args1));
		va_end(args1);
		std::vsnprintf(buf.data(), buf.size(), fmt, args2);
		va_end(args2);
		std::string str = std::string(buf.begin(), buf.end());
		str.erase(std::remove_if(str.begin(), str.end(), ([](const char& ch)->bool{return (ch==char('\0'));}) ), str.end());
		return str;
	}

	inline int                    System(const std::string& command) { return std::system(command.c_str()); }
	inline void                   Exit(int exit_code = EXIT_SUCCESS) { std::exit(exit_code); } // EXIT_SUCCESS, EXIT_FAILURE
	inline void                   AtExit( void(*func)() ) { std::atexit(func); }
	inline void                   Abort() { std::abort(); }
	inline void                   Terminate() { std::terminate(); }
	inline std::terminate_handler SetTerminate(std::terminate_handler func) { return std::set_terminate(func); }
	inline std::terminate_handler GetTerminate() { return std::get_terminate(); }

	inline bool TestFile(const std::string& file, char opt = 'e') { return (MgntSys::System(StrFmt("test -%c \"%s\"", opt, file.c_str())) == 0); }

	inline std::string GetEnv(const std::string& env_var) { return std::string(std::getenv(env_var.c_str())); }
	inline int         SetEnv(const std::string& envname, const std::string& envval, int overwrite = 0) { return setenv(envname.c_str(), envval.c_str(), overwrite); }
	inline int         UnsetEnv(const std::string& envname) { return unsetenv(envname.c_str()); }
	
	template< class Rep, class Period >
	inline void SleepFor(const std::chrono::duration<Rep, Period>& sleep_duration) { std::this_thread::sleep_for(sleep_duration); } 
	template< class Clock, class Duration >
	inline void SleepUntil(const std::chrono::time_point<Clock, Duration>& sleep_time) { std::this_thread::sleep_until(sleep_time); } 

	typedef std::string          MESSAGE;
	typedef std::vector<MESSAGE> MESSAGES;
	bool ShowMessage(const MESSAGE& mainInfo = MESSAGE(), bool expr = true, const MgntSys::MESSAGES& messages = MgntSys::MESSAGES(), std::ostream& out = std::cout, const MESSAGE& type = MESSAGE("MESSAGE")) {
		if (!expr) return expr;
		MESSAGE title = MESSAGE(StrFmt("==%s==", type.c_str()));
		out << CStrFmt("%-15s%s\n", title.c_str(), mainInfo.c_str());
		for (auto str : messages)
			if (str != "")
				out << CStrFmt("%-15s%s\n", "", str.c_str()); 
		return expr;
	}

	inline bool Message(const MESSAGE& mainInfo = MESSAGE(), bool expr = true, const MgntSys::MESSAGES& messages = MgntSys::MESSAGES(), std::ostream& out = std::cout) { 
		return ShowMessage(mainInfo, expr, messages, out, "MESSAGE");
	}
	inline bool Message(const MESSAGE& mainInfo, bool expr, const MESSAGE& message, std::ostream& out = std::cout) {
		return ShowMessage(mainInfo, expr, MgntSys::MESSAGES({message}), out, "MESSAGE");
	}
	inline bool Message(const MESSAGE& mainInfo, const MgntSys::MESSAGES& messages, std::ostream& out = std::cout) { 
		return ShowMessage(mainInfo, true, messages, out, "MESSAGE");
	}
	inline bool Message(const MESSAGE& mainInfo, const MESSAGE& message, std::ostream& out = std::cout) {
		return ShowMessage(mainInfo, true, MgntSys::MESSAGES({message}), out, "MESSAGE");
	}

	inline bool Warning(const MESSAGE& mainInfo = MESSAGE(), bool expr = true, const MgntSys::MESSAGES& messages = MgntSys::MESSAGES(), std::ostream& out = std::clog) { 
		return ShowMessage(mainInfo, expr, messages, out, "WARNING");
	}
	inline bool Warning(const MESSAGE& mainInfo, bool expr, const MESSAGE& message, std::ostream& out = std::clog) {
		return ShowMessage(mainInfo, expr, MgntSys::MESSAGES({message}), out, "WARNING");
	}
	inline bool Warning(const MESSAGE& mainInfo, const MgntSys::MESSAGES& messages, std::ostream& out = std::clog) { 
		return ShowMessage(mainInfo, true, messages, out, "WARNING");
	}
	inline bool Warning(const MESSAGE& mainInfo, const MESSAGE& message, std::ostream& out = std::clog) {
		return ShowMessage(mainInfo, true, MgntSys::MESSAGES({message}), out, "WARNING");
	}
	
	inline bool Error(const MESSAGE& mainInfo = MESSAGE(), bool expr = true, const MgntSys::MESSAGES& messages = MgntSys::MESSAGES(), std::ostream& out = std::cerr) { 
		return ShowMessage(mainInfo, expr, messages, out, "ERROR");
	}
	inline bool Error(const MESSAGE& mainInfo, bool expr, const MESSAGE& message, std::ostream& out = std::cerr) {
		return ShowMessage(mainInfo, expr, MgntSys::MESSAGES({message}), out, "ERROR");
	}
	inline bool Error(const MESSAGE& mainInfo, const MgntSys::MESSAGES& messages, std::ostream& out = std::cerr) { 
		return ShowMessage(mainInfo, true, messages, out, "ERROR");
	}
	inline bool Error(const MESSAGE& mainInfo, const MESSAGE& message, std::ostream& out = std::cerr) {
		return ShowMessage(mainInfo, true, MgntSys::MESSAGES({message}), out, "ERROR");
	}

	std::vector<std::string> ReadDirectory(const std::string& path = std::string("."), const std::string& patt = std::string("")) {
		bool searchOpt = (patt != "");
		std::vector<std::string> list;
		if (!MgntSys::TestFile(path, 'd')) return list;
		DIR * dp = opendir(path.c_str());
		if (dp == nullptr) return list;
		while (true)
		{
			dirent * de = readdir(dp);
			if (de == nullptr) break;
			std::string str(de->d_name);
			if (str.find_first_of('~') == str.size()-1) continue;
			if (str == std::string(".") || str == std::string("..")) continue;
			if (searchOpt && str.find(patt) == std::string::npos) continue;
			list.push_back(str);
		}
		closedir(dp);
		std::sort(list.begin(), list.end());
		return list;
	}

	// TODO : This is only for testing.
	//        1) std::cin timeout ?!  -> condition_variable::wait_for
	inline void Console() {
		long int count = 1;
		std::string cmdin = std::string();
		console :
		cmdin = "";
		while ("" == cmdin) { 
			std::cout << CStrFmt("C++ Bash Console [%ld] >  ", count);
			std::getline(std::cin, cmdin); 
		}
		if (cmdin == "exit" || cmdin == "EXIT") return;
		MgntSys::System(cmdin);
		count++;
		goto console;
	}
}

#endif // __MgntSys_H__
