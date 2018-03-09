#ifndef __TRACKLibs_Sys_H__
#define __TRACKLibs_Sys_H__
#include <iostream>
#include <string>
#include <tuple>

#ifndef __HAS_LOGSYS__
#define __HAS_LOGSYS__

#define STR(fmt, ...)  (std::string(Form(fmt, ##__VA_ARGS__)))
#define CSTR(fmt, ...) (std::string(Form(fmt, ##__VA_ARGS__)).c_str())

#define COUT(fmt, ...) (std::cout << CSTR(fmt, ##__VA_ARGS__))
#define CLOG(fmt, ...) (std::clog << CSTR(fmt, ##__VA_ARGS__))
#define CERR(fmt, ...) (std::cerr << CSTR(fmt, ##__VA_ARGS__))

//#define LOC_ADDR() (std::string(Form("(LINE %d) [FUNC  %s] {FILE  %s} :  ", __LINE__, __func__, __FILE__)))

#endif // __HAS_LOGSYS__


namespace TrackSys {
namespace Sys {

inline std::string GetEnv(const std::string& var) { return std::string(std::getenv(var.c_str())); }

} // namespace Sys
} // namespace TrackSys


#endif // __TRACKLibs_Sys_H__
