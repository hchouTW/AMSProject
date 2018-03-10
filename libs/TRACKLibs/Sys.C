#ifndef __TRACKLibs_Sys_C__
#define __TRACKLibs_Sys_C__


namespace TrackSys {
namespace Sys {


void ShowMsg(const Msg& info, const Msg& msg, const Msg& title, std::ostream& out) {
	Msg space(title.size()+6, ' ');
    Msg ostr = STR("==%s==  %s\n", title.c_str(), info.c_str());
	if (msg != "")
		ostr += STR("%s%s\n", space.c_str(), msg.c_str());
    out << ostr.c_str();
}

void ShowMsg(const Msg& info, const Msgs& msgs, const Msg& title, std::ostream& out) {
	Msg space(title.size()+6, ' ');
    Msg ostr = STR("==%s==  %s\n", title.c_str(), info.c_str());
	for (auto&& str : msgs)
		if (str != "")
			ostr += STR("%s%s\n", space.c_str(), str.c_str()); 
    out << ostr.c_str();
}


} // namespace Sys
} // namespace TrackSys


#endif // __TRACKLibs_Sys_C__
