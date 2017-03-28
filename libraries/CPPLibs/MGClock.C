#ifndef __MGClock_C__
#define __MGClock_C__


#include "MGClock.h"


MGClock::TTime * MGClock::ConvertFromUTimeToTTime(UTime utime, ClockType type) { 
	TTime * ttime = nullptr;
	switch(type) {
		case ClockType::UTC   : ttime = std::gmtime(&utime);    break;
		case ClockType::LOCAL : ttime = std::localtime(&utime); break;
		default               : break;
	}
	return ttime;
}


MGClock::CTime MGClock::ConvertFromTTimeToCTime(TTime * ttime, const std::string& fmt) {
	if (fmt == "") {
		CTime ctime = std::asctime(ttime);
		ctime.erase(ctime.begin()+ctime.find_last_of('\n'));
		return ctime;
	}
	else {
		// TODO : using std::put_time on c++11,  (strftime c++)
		char buffer[256];
		strftime(buffer, 256, fmt.c_str(), ttime);
		CTime ctime = buffer;
		std::cerr << ctime.size() << std::endl;
		return ctime;
	}
	return "";
}


template<class Time>
std::ostream& MGClock::Print(Time timpnt, const std::string& fmt, ClockType type, std::ostream& out) {
	std::string timetype = "";
	switch(type) {
		case ClockType::UTC   : timetype = "  UTC"; break;
		case ClockType::LOCAL : timetype = "LOCAL"; break;
		default               : break;
	}
	UTime utime = ConvertToUTime(timpnt);
	CTime ctime = ConvertFromUTimeToCTime(utime, "", type);
	out << CStrFmt("UNIX (%ld)    %5s (%s)", utime, timetype.c_str(), ctime.c_str());
	return out;
}


template <class Clock, class Time, class Duration>
std::ostream& MGClock::Stopwatch<Clock, Time, Duration>::print(std::ostream& out, MGClock::ClockType type) {
	Duration&& durt = duration();
	double time = std::chrono::duration<double>(durt).count();
	
	MGClock::Hours        hours   = std::chrono::duration_cast<MGClock::Hours>(durt);   durt -= hours;
	MGClock::Minutes      minutes = std::chrono::duration_cast<MGClock::Minutes>(durt); durt -= minutes;
	MGClock::FloatSeconds seconds = std::chrono::duration_cast<MGClock::FloatSeconds>(durt);
	int      hr = hours.count();
	int      mn = minutes.count();
	double   sc = seconds.count();
	MGClock::CTime ctime = StrFmt("%-3d HR %-2d MIN %12.9f SEC", hr, mn, sc);

	out << "===========================  MGClock::Stopwatch  ===========================\n";
	out << "==  START TIME : "; MGClock::Print(fTime.first , "", type, out) << "    ==\n";
	out << "==  STOP  TIME : "; MGClock::Print(fTime.second, "", type, out) << "    ==\n";
	out << CStrFmt("==  Duration   : %-30s     (%18.9f)  ==\n", ctime.c_str(), time);
	out << "============================================================================\n";
	return out;
}


#endif // __MGClock_C__
