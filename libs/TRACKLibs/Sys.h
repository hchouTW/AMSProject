#ifndef __HAS_LOGSYS__
#define __HAS_LOGSYS__

#define STR(fmt, ...)  (std::string(Form(fmt, ##__VA_ARGS__)))
#define CSTR(fmt, ...) (std::string(Form(fmt, ##__VA_ARGS__)).c_str())

#define COUT(fmt, ...) (std::cout << CSTR(fmt, ##__VA_ARGS__))
#define CLOG(fmt, ...) (std::clog << CSTR(fmt, ##__VA_ARGS__))
#define CERR(fmt, ...) (std::cerr << CSTR(fmt, ##__VA_ARGS__))

#define LOCADR() (std::string(Form("(LINE %d) [FUNC  %s] {FILE  %s} :  ", __LINE__, __func__, __FILE__)))

#endif // __HAS_LOGSYS__


#ifndef __TRACKLibs_Sys_H__
#define __TRACKLibs_Sys_H__
#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>


namespace TrackSys {
namespace Sys {

inline Bool_t       IsEnv(const std::string& var) { return (getenv(var.c_str()) != nullptr); }
inline std::string GetEnv(const std::string& var) { return std::string(getenv(var.c_str())); }
inline void        PutEnv(const std::string& var, const std::string& val = "") { putenv(const_cast<char*>( STR("%s=%s", var.c_str(), val.c_str()).c_str() )); }

using Msg  = std::string;
using Msgs = std::vector<Msg>;

inline void ShowMsg(const Msg& info = Msg(), const Msg& msg = Msg(), const Msg& type = Msg("MSG"), std::ostream& out = std::cout);
inline void ShowMsgExit(const Msg& info = Msg(), const Msg& msg = Msg(), const Msg& type = Msg("MSG"), std::ostream& out = std::cout) { ShowMsg(info, msg, type, out); std::exit(EXIT_FAILURE); }

inline void ShowMsg(const Msg& info, const Msgs& msgs, const Msg& type = Msg("MSG"), std::ostream& out = std::cout);
inline void ShowMsgExit(const Msg& info, const Msgs& msgs, const Msg& type = Msg("MSG"), std::ostream& out = std::cout) { ShowMsg(info, msgs, type, out); std::exit(EXIT_FAILURE); }

inline void ShowLog(const Msg& info = Msg(), const Msg& msg = Msg()) { ShowMsg(info, msg, "Log", std::cerr); }
inline void ShowLogExit(const Msg& info = Msg(), const Msg& msg = Msg()) { ShowMsgExit(info, msg, "Log", std::cerr); }

inline void ShowLog(const Msg& info, const Msgs& msgs) { ShowMsg(info, msgs, "Warning", std::cerr); }
inline void ShowLogExit(const Msg& info, const Msgs& msgs) { ShowMsgExit(info, msgs, "Warning", std::cerr); }

inline void ShowWarning(const Msg& info = Msg(), const Msg& msg = Msg()) { ShowMsg(info, msg, "Warning", std::cerr); }
inline void ShowWarningExit(const Msg& info = Msg(), const Msg& msg = Msg()) { ShowMsgExit(info, msg, "Warning", std::cerr); }

inline void ShowWarning(const Msg& info, const Msgs& msgs) { ShowMsg(info, msgs, "Warning", std::cerr); }
inline void ShowWarningExit(const Msg& info, const Msgs& msgs) { ShowMsgExit(info, msgs, "Warning", std::cerr); }

inline void ShowError(const Msg& info = Msg(), const Msg& msg = Msg()) { ShowMsg(info, msg, "ERROR", std::cerr); }
inline void ShowErrorExit(const Msg& info = Msg(), const Msg& msg = Msg()) { ShowMsgExit(info, msg, "ERROR", std::cerr); }

inline void ShowError(const Msg& info, const Msgs& msgs) { ShowMsg(info, msgs, "ERROR", std::cerr); }
inline void ShowErrorExit(const Msg& info, const Msgs& msgs) { ShowMsgExit(info, msgs, "ERROR", std::cerr); }

} // namespace Sys
} // namespace TrackSys


#endif // __TRACKLibs_Sys_H__
