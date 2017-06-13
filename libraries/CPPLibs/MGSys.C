#ifndef __CPPLibs_MGSys_C__
#define __CPPLibs_MGSys_C__

#include "MGSys.h"


namespace MGSys {

template<typename... Args>
std::string StringFormat(const std::string& fmt, Args... args) {
	std::vector<char> buf(1+std::snprintf(nullptr, 0, fmt.c_str(), args...));
	std::snprintf(buf.data(), buf.size(), fmt.c_str(), args...);
	std::string&& str = std::string(buf.begin(), buf.end());
	str.erase(std::remove_if(str.begin(), str.end(), ([](const char& ch)->bool{return (ch==char('\0'));}) ), str.end());
	return str;
}


void ShowInfo(const Message& mainInfo, const MGSys::Message& message, const Message& title, std::ostream& out) {
	std::string space(title.size()+6, ' ');
	out << CSTR_FMT("==%s==  %s\n", title.c_str(), mainInfo.c_str());
	if (message != "")
		out << CSTR_FMT("%s%s\n", space.c_str(), message.c_str());
}


void ShowInfo(const Message& mainInfo, const MGSys::Messages& messages, const Message& title, std::ostream& out) {
	std::string space(title.size()+6, ' ');
	out << CSTR_FMT("==%s==  %s\n", title.c_str(), mainInfo.c_str());
	for (auto&& str : messages)
		if (str != "")
			out << CSTR_FMT("%s%s\n", space.c_str(), str.c_str()); 
}

} // namespace MGSys	

#endif // __CPPLibs_MGSys_C__
