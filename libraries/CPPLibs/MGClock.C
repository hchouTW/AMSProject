#ifndef __CPPLibs_MGClock_C__
#define __CPPLibs_MGClock_C__

#include "MGClock.h"

namespace MGClock {

MTime * ConvertFromUTimeToMTime(UTime utime, ClockType type) { 
	MTime * mtime = nullptr;
	switch(type) {
		case ClockType::UTC   : mtime = std::gmtime(&utime);    break;
		case ClockType::LOCAL : mtime = std::localtime(&utime); break;
		default               : break;
	}
	return mtime;
}


CTime ConvertFromMTimeToCTime(MTime * mtime, const std::string& fmt) {
	if (fmt == "") {
		CTime ctime = std::asctime(mtime);
		ctime.erase(ctime.begin()+ctime.find_last_of('\n'));
		return ctime;
	}
	else {
		// TODO : using std::put_time on c++11,  (strftime c++)
		char buffer[256];
		strftime(buffer, 256, fmt.c_str(), mtime);
		CTime ctime = buffer;
		std::cerr << ctime.size() << std::endl;
		return ctime;
	}
	return "";
}


template<class Time>
std::ostream& Print(Time timpnt, ClockType type, const std::string& fmt, std::ostream& out) {
	std::string timetype = "";
	switch(type) {
		case ClockType::UTC   : timetype = "  UTC"; break;
		case ClockType::LOCAL : timetype = "LOCAL"; break;
		default               : break;
	}
	UTime utime = ConvertToUTime(timpnt);
	CTime ctime = ConvertFromUTimeToCTime(utime, type, fmt);
	out << CSTR_FMT("UNIX (%ld)    %5s (%s)", utime, timetype.c_str(), ctime.c_str());
	return out;
}


template <class Clock, class Time, class Duration>
std::ostream& Stopwatch<Clock, Time, Duration>::print(ClockType type, std::ostream& out) const {
	Duration&& durt = duration();
	double time = std::chrono::duration<double>(durt).count();
	
	Hours        hours   = std::chrono::duration_cast<Hours>(durt);   durt -= hours;
	Minutes      minutes = std::chrono::duration_cast<Minutes>(durt); durt -= minutes;
	FloatSeconds seconds = std::chrono::duration_cast<FloatSeconds>(durt);
	int      hr = hours.count();
	int      mn = minutes.count();
	double   sc = seconds.count();
	CTime ctime = STR_FMT("%-3d HR %-2d MIN %12.9f SEC", hr, mn, sc);

    std::string outstr;
	outstr += "===========================  Stopwatch  ===========================\n";
	outstr += "==  START TIME : "; Print(times_.first , type, "", out) << "    ==\n";
	outstr += "==  STOP  TIME : "; Print(times_.second, type, "", out) << "    ==\n";
	outstr += STR_FMT("==  Duration   : %-30s     (%18.9f)  ==\n", ctime.c_str(), time);
	outstr += "============================================================================\n";
    out << outstr.c_str();
	return out;
}

} // namespace MGClock

#endif // __CPPLibs_MGClock_C__
