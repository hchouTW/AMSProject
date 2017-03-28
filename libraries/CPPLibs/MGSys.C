#ifndef __MGSys_C__
#define __MGSys_C__

#include "MGSys.h"


template<typename... Args>
std::string MGSys::StringFormat(const std::string& fmt, Args... args) {
	std::vector<char> buf(1+std::snprintf(nullptr, 0, fmt.c_str(), args...));
	std::snprintf(buf.data(), buf.size(), fmt.c_str(), args...);
	std::string&& str = std::string(buf.begin(), buf.end());
	str.erase(std::remove_if(str.begin(), str.end(), ([](const char& ch)->bool{return (ch==char('\0'));}) ), str.end());
	return str;
}


void MGSys::ShowInfo(const MESSAGE& mainInfo, const MGSys::MESSAGE& message, const MESSAGE& title, std::ostream& out) {
	std::string space(title.size()+6, ' ');
	out << CStrFmt("==%s==  %s\n", title.c_str(), mainInfo.c_str());
	if (message != "")
		out << CStrFmt("%s%s\n", space.c_str(), message.c_str());
}


void MGSys::ShowInfo(const MESSAGE& mainInfo, const MGSys::MESSAGES& messages, const MESSAGE& title, std::ostream& out) {
	std::string space(title.size()+6, ' ');
	out << CStrFmt("==%s==  %s\n", title.c_str(), mainInfo.c_str());
	for (auto&& str : messages)
		if (str != "")
			out << CStrFmt("%s%s\n", space.c_str(), str.c_str()); 
}
	

#endif // __MGSys_C__
