#ifndef __MGSys_H__
#define __MGSys_H__
#include <iostream>
#include <cstdlib>
#include <exception>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <dirent.h>

#define StrFmt(fmt, ...)  (MGSys::StringFormat(fmt, ##__VA_ARGS__))
#define CStrFmt(fmt, ...) ((MGSys::StringFormat(fmt, ##__VA_ARGS__)).c_str())

#define COUT(fmt, ...) (std::cout << ((MGSys::StringFormat(fmt, ##__VA_ARGS__)).c_str()))
#define CLOG(fmt, ...) (std::clog << ((MGSys::StringFormat(fmt, ##__VA_ARGS__)).c_str()))
#define CERR(fmt, ...) (std::cerr << ((MGSys::StringFormat(fmt, ##__VA_ARGS__)).c_str()))

#define LocAddr() (MGSys::StringFormat("(LINE %d)  [FUNC  %s]  {FILE  %s} :  ", __LINE__, __func__, __FILE__))

namespace MGSys {
	inline int                    System(const std::string& command) { return std::system(command.c_str()); }
	inline void                   Exit(int exit_code = EXIT_SUCCESS) { std::exit(exit_code); } // EXIT_SUCCESS, EXIT_FAILURE
	inline void                   AtExit( void(*func)() ) { std::atexit(func); }
	inline void                   Abort() { std::abort(); }
	inline void                   Terminate() { std::terminate(); }
	inline std::terminate_handler SetTerminate(std::terminate_handler func) { return std::set_terminate(func); }
	inline std::terminate_handler GetTerminate() { return std::get_terminate(); }
	
	inline std::string GetEnv(const std::string& envvar) { return std::string(std::getenv(envvar.c_str())); }
	inline int         SetEnv(const std::string& envname, const std::string& envval, int overwrite = 0) { return setenv(envname.c_str(), envval.c_str(), overwrite); }
	inline int         UnsetEnv(const std::string& envname) { return unsetenv(envname.c_str()); }
	
	template< class Rep, class Period >
	inline void SleepFor(const std::chrono::duration<Rep, Period>& sleep_duration) { std::this_thread::sleep_for(sleep_duration); } 
	template< class Clock, class Duration >
	inline void SleepUntil(const std::chrono::time_point<Clock, Duration>& sleep_time) { std::this_thread::sleep_until(sleep_time); } 
	
	template <typename... Args>
	std::string StringFormat(const std::string& fmt, Args... args);

	typedef std::string          MESSAGE;
	typedef std::vector<MESSAGE> MESSAGES;
	void ShowInfo(const MESSAGE& mainInfo = MESSAGE(), const MGSys::MESSAGE& message = MGSys::MESSAGE(), const MESSAGE& type = MESSAGE("INFO"), std::ostream& out = std::cout);
	void ShowInfo(const MESSAGE& mainInfo, const MGSys::MESSAGES& messages, const MESSAGE& type = MESSAGE("MESSAGE"), std::ostream& out = std::cout);

	inline void ShowMessage(const MESSAGE& mainInfo = MESSAGE(), const MGSys::MESSAGE& message = MGSys::MESSAGE(), std::ostream& out = std::cout) { ShowInfo(mainInfo, message, "MESSAGE", out); }
	inline void ShowMessage(const MESSAGE& mainInfo, const MGSys::MESSAGES& messages, std::ostream& out = std::cout) { ShowInfo(mainInfo, messages, "MESSAGE", out); }
	
	inline void ShowMessageAndExit(const MESSAGE& mainInfo = MESSAGE(), const MGSys::MESSAGE& message = MGSys::MESSAGE(), std::ostream& out = std::cout) { ShowInfo(mainInfo, message, "MESSAGE", out); Exit(EXIT_SUCCESS); }
	inline void ShowMessageAndExit(const MESSAGE& mainInfo, const MGSys::MESSAGES& messages, std::ostream& out = std::cout) { ShowInfo(mainInfo, messages, "MESSAGE", out); Exit(EXIT_SUCCESS); }

	inline void ShowWarning(const MESSAGE& mainInfo = MESSAGE(), const MGSys::MESSAGE& message = MGSys::MESSAGE(), std::ostream& out = std::clog) { ShowInfo(mainInfo, message, "WARNING", out); }
	inline void ShowWarning(const MESSAGE& mainInfo, const MGSys::MESSAGES& messages, std::ostream& out = std::clog) { ShowInfo(mainInfo, messages, "WARNING", out); }
	
	inline void ShowWarningAndExit(const MESSAGE& mainInfo = MESSAGE(), const MGSys::MESSAGE& message = MGSys::MESSAGE(), std::ostream& out = std::clog) { ShowInfo(mainInfo, message, "WARNING", out); Exit(EXIT_FAILURE); }
	inline void ShowWarningAndExit(const MESSAGE& mainInfo, const MGSys::MESSAGES& messages, std::ostream& out = std::clog) { ShowInfo(mainInfo, messages, "WARNING", out); Exit(EXIT_FAILURE); }
	
	inline void ShowError(const MESSAGE& mainInfo = MESSAGE(), const MGSys::MESSAGE& message = MGSys::MESSAGE(), std::ostream& out = std::cerr) { ShowInfo(mainInfo, message, "ERROR", out); }
	inline void ShowError(const MESSAGE& mainInfo, const MGSys::MESSAGES& messages, std::ostream& out = std::cerr) { ShowInfo(mainInfo, messages, "ERROR", out); }
	
	inline void ShowErrorAndExit(const MESSAGE& mainInfo = MESSAGE(), const MGSys::MESSAGE& message = MGSys::MESSAGE(), std::ostream& out = std::cerr) { ShowInfo(mainInfo, message, "ERROR", out); Exit(EXIT_FAILURE); }
	inline void ShowErrorAndExit(const MESSAGE& mainInfo, const MGSys::MESSAGES& messages, std::ostream& out = std::cerr) { ShowInfo(mainInfo, messages, "ERROR", out); Exit(EXIT_FAILURE); }

	// TODO : update by C++17 (filesystem)
	inline bool TestFile(const std::string& file, char opt = 'e') { return (MGSys::System(StrFmt("test -%c \"%s\"", opt, file.c_str())) == 0); }
	
	// TODO : update by C++17 (filesystem)
	std::vector<std::string> ReadDirectory(const std::string& path = std::string("."), const std::string& patt = std::string("")) {
		bool searchOpt = (patt != "");
		std::vector<std::string> list;
		if (!MGSys::TestFile(path, 'd')) return list;
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
	void Console() {
		long int count = 1;
		std::string cmdin = std::string();
		console :
		cmdin = "";
		while ("" == cmdin) { 
			std::cout << CStrFmt("C++ Bash Console [%ld] >  ", count);
			std::getline(std::cin, cmdin); 
		}
		if (cmdin == "exit" || cmdin == "EXIT") return;
		MGSys::System(cmdin);
		count++;
		goto console;
	}
}

#endif // __MGSys_H__
