#ifndef __MgntClock_H__
#define __MgntClock_H__
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include <cstdint>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <string>
#include <iomanip>
#include <algorithm>
//#include <cxxabi.h>

#include "MgntSys.h"

namespace MgntClock {
	enum class ClockType { LOCAL = 0, UTC = 1 };
	
	typedef std::time_t                                    UTime;
	typedef struct std::tm                                 TTime;
	typedef std::string                                    CTime;
	typedef std::chrono::duration<double>                  FloatSeconds;
	typedef std::chrono::nanoseconds                       Nanoseconds;
	typedef std::chrono::microseconds                      Microseconds;
	typedef std::chrono::milliseconds                      Milliseconds;
	typedef std::chrono::seconds                           Seconds;
	typedef std::chrono::minutes                           Minutes;
	typedef std::chrono::hours                             Hours;

	typedef std::chrono::system_clock                      SysClock;
	typedef std::chrono::system_clock::time_point          SysTime;
	typedef std::chrono::system_clock::duration            SysDuration;
	typedef std::chrono::steady_clock                      StdClock;
	typedef std::chrono::steady_clock::time_point          StdTime;
	typedef std::chrono::steady_clock::duration            StdDuration;
	typedef std::chrono::high_resolution_clock             HrsClock;
	typedef std::chrono::high_resolution_clock::time_point HrsTime;
	typedef std::chrono::high_resolution_clock::duration   HrsDuration;

	inline MgntClock::UTime   ConvertFromTTimeToUTime(MgntClock::TTime * ttime) {
		MgntClock::UTime utime = std::mktime(ttime);
		return utime;
	}

	inline MgntClock::TTime * ConvertFromUTimeToTTime(MgntClock::UTime utime, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) { 
		MgntClock::TTime * ttime = nullptr;
		switch(type) {
			case ClockType::LOCAL : ttime = std::localtime(&utime); break;
			case ClockType::UTC   : ttime = std::gmtime(&utime);    break;
			default               : break;
		}
		return ttime;
	}

	inline MgntClock::CTime   ConvertFromTTimeToCTime(MgntClock::TTime * ttime, const std::string& fmt = "") {
		if (fmt == std::string("")) {
			MgntClock::CTime ctime = std::asctime(ttime);
			ctime.erase(ctime.begin()+ctime.find_last_of('\n'));
			return ctime;
		}
		else {
			// TODO : using std::put_time on c++11,  (strftime c++)
			char buffer[200];
			strftime(buffer, 200, fmt.c_str(), ttime);
			MgntClock::CTime ctime = buffer;
			return ctime;
		}
		return "";
	}
	
	inline MgntClock::CTime   ConvertFromUTimeToCTime(MgntClock::UTime utime, MgntClock::ClockType type = MgntClock::ClockType::LOCAL, const std::string& fmt = "") {
		MgntClock::TTime * ttime = MgntClock::ConvertFromUTimeToTTime(utime, type);
		MgntClock::CTime   ctime = MgntClock::ConvertFromTTimeToCTime(ttime, fmt);
		return ctime;
	}
	
	inline MgntClock::SysTime ConvertFromUTime(MgntClock::UTime timpnt) { return SysClock::from_time_t(timpnt); }
	inline MgntClock::SysTime ConvertFromTTime(MgntClock::TTime * ttime) { return MgntClock::ConvertFromUTime(MgntClock::ConvertFromTTimeToUTime(ttime)); }
	
	inline MgntClock::UTime   ConvertToUTime(MgntClock::SysTime timpnt) { return SysClock::to_time_t(timpnt); }
	inline MgntClock::UTime   ConvertToUTime(MgntClock::StdTime timpnt) {
		MgntClock::SysTime systimpnt = MgntClock::SysClock::now() + (timpnt - MgntClock::StdClock::now());
		return SysClock::to_time_t(systimpnt); 
	}

	inline MgntClock::TTime * ConvertToTTime(MgntClock::SysTime timpnt, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) {
		MgntClock::UTime   utime = MgntClock::ConvertToUTime(timpnt);
		MgntClock::TTime * ttime = MgntClock::ConvertFromUTimeToTTime(utime, type);
		return ttime;
	}
	inline MgntClock::TTime * ConvertToTTime(MgntClock::StdTime timpnt, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) {
		MgntClock::UTime   utime = MgntClock::ConvertToUTime(timpnt);
		MgntClock::TTime * ttime = MgntClock::ConvertFromUTimeToTTime(utime, type);
		return ttime;
	}

	inline MgntClock::CTime   ConvertToCTime(MgntClock::SysTime timpnt, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) {
		MgntClock::UTime utime = MgntClock::ConvertToUTime(timpnt);
		MgntClock::CTime ctime = MgntClock::ConvertFromUTimeToCTime(utime, type);
		return ctime;
	}
	static inline MgntClock::CTime   ConvertToCTime(MgntClock::StdTime timpnt, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) {
		MgntClock::UTime utime = MgntClock::ConvertToUTime(timpnt);
		MgntClock::CTime ctime = MgntClock::ConvertFromUTimeToCTime(utime, type);
		return ctime;
	}

	template<class Time>
	inline std::ostream& Print(Time timpnt, std::ostream& out = std::cout, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) {
		std::string timetype = "";
		switch(type) {
			case ClockType::LOCAL : timetype = "LOCAL"; break;
			case ClockType::UTC   : timetype = "  UTC"; break;
			default               : break;
		}
		MgntClock::UTime utime = MgntClock::ConvertToUTime(timpnt);
		MgntClock::CTime ctime = MgntClock::ConvertFromUTimeToCTime(utime, type);
		out << CStrFmt("UNIT( %ld )    %s{ %s }", utime, timetype.c_str(), ctime.c_str());
		return out;
	}
	
	template<class Clock, class Time>
	inline std::ostream& PrintNow(std::ostream& out = std::cout, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) {
		Time timpnt = Clock::now();
		return MgntClock::Print<Time>(timpnt, out, type);	
	}

	template <class Clock, class Time, class Duration> class Stopwatch;
	typedef MgntClock::Stopwatch<MgntClock::SysClock, MgntClock::SysTime, MgntClock::SysDuration> SysStopwatch;
	typedef MgntClock::Stopwatch<MgntClock::StdClock, MgntClock::StdTime, MgntClock::StdDuration> StdStopwatch;
	typedef MgntClock::Stopwatch<MgntClock::HrsClock, MgntClock::HrsTime, MgntClock::HrsDuration> HrsStopwatch;
}


template <class Clock, class Time, class Duration>
class MgntClock::Stopwatch {
	public :
		static inline Time Now() { return Clock::now(); }
		static inline std::ostream& Print(Time timpnt, std::ostream& out = std::cout, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) { return MgntClock::Print<Time>(timpnt, out, type); }
		static inline std::ostream& PrintNow(std::ostream& out = std::cout, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) { return MgntClock::PrintNow<Clock, Time>(out, type); }

	public :
		Stopwatch() { start(); stop(); }
		~Stopwatch() {}

		Time & start() { fTime.first  = Stopwatch::Now(); return fTime.first;  }
		Time & stop()  { fTime.second = Stopwatch::Now(); return fTime.second; }
		
		Duration duration()      { return (fTime.second - fTime.first); }
		double   time()          { return std::chrono::duration<double>(duration()).count(); }
		MgntClock::CTime ctime() { 
			Duration durt = duration();
			MgntClock::Hours        hours   = std::chrono::duration_cast<MgntClock::Hours>(durt);   durt -= hours;
			MgntClock::Minutes      minutes = std::chrono::duration_cast<MgntClock::Minutes>(durt); durt -= minutes;
			MgntClock::FloatSeconds seconds = std::chrono::duration_cast<MgntClock::FloatSeconds>(durt);
			int      h = hours.count();
			int      m = minutes.count();
			double   s = seconds.count();
			MgntClock::CTime ctime = StrFmt("%-3d HR %-2d MIN %12.9f SEC", h, m, s);
			return ctime;
		}

		std::ostream& print(std::ostream& out = std::cout, MgntClock::ClockType type = MgntClock::ClockType::LOCAL) {
			//int status = -4;
			//std::string className = abi::__cxa_demangle(typeid(decltype(this)).name(), nullptr, nullptr, &status);
			out << "============================= MgntClock::Stopwatch =============================\n";
			out << "==  START TIME : "; Print(fTime.first,  out, type) << "  ==\n";
			out << "==  STOP  TIME : "; Print(fTime.second, out, type) << "  ==\n";
			out << CStrFmt("==  Duration   : %-30s     (%18.9f)  ==\n", ctime().c_str(), time());
			out << "============================================================================\n";
			return out;
		}

	private :
		std::pair<Time, Time> fTime;
};


#endif // __MgntClock_H__
