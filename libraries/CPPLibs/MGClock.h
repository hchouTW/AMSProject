#ifndef __MGClock_H__
#define __MGClock_H__


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


#include "MGSys.h"


namespace MGClock {
	enum class ClockType { UTC = 0, LOCAL = 1 };
	
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
	
	template <class Clock, class Time, class Duration> 
	class Stopwatch;
	typedef Stopwatch<SysClock, SysTime, SysDuration> SysStopwatch;
	typedef Stopwatch<StdClock, StdTime, StdDuration> StdStopwatch;
	typedef Stopwatch<HrsClock, HrsTime, HrsDuration> HrsStopwatch;

	inline UTime   ConvertFromTTimeToUTime(TTime * ttime) { return std::mktime(ttime); }
	inline TTime * ConvertFromUTimeToTTime(UTime utime, ClockType type = ClockType::UTC); 
	inline CTime   ConvertFromTTimeToCTime(TTime * ttime, const std::string& fmt = "");
	inline CTime   ConvertFromUTimeToCTime(UTime utime, const std::string& fmt = "", ClockType type = ClockType::UTC) { return ConvertFromTTimeToCTime(ConvertFromUTimeToTTime(utime, type), fmt); }

	inline SysTime ConvertFromUTime(UTime timpnt)  { return SysClock::from_time_t(timpnt); }
	inline SysTime ConvertFromTTime(TTime * ttime) { return ConvertFromUTime(ConvertFromTTimeToUTime(ttime)); }
	inline UTime   ConvertToUTime(SysTime timpnt) { return SysClock::to_time_t(timpnt); }
	inline UTime   ConvertToUTime(StdTime timpnt) { return SysClock::to_time_t( (SysClock::now() + (timpnt - StdClock::now())) ); }

	inline TTime * ConvertToTTime(SysTime timpnt, ClockType type = ClockType::UTC) { return ConvertFromUTimeToTTime(ConvertToUTime(timpnt), type); }
	inline TTime * ConvertToTTime(StdTime timpnt, ClockType type = ClockType::UTC) { return ConvertFromUTimeToTTime(ConvertToUTime(timpnt), type); }

	inline CTime   ConvertToCTime(SysTime timpnt, const std::string& fmt = "", ClockType type = ClockType::UTC) { return ConvertFromUTimeToCTime(ConvertToUTime(timpnt), fmt, type); }
	inline CTime   ConvertToCTime(StdTime timpnt, const std::string& fmt = "", ClockType type = ClockType::UTC) { return ConvertFromUTimeToCTime(ConvertToUTime(timpnt), fmt, type); }

	template<class Time>
	std::ostream& Print(Time timpnt, const std::string& fmt = "", ClockType type = ClockType::UTC, std::ostream& out = std::cout);
}


template <class Clock, class Time, class Duration>
class MGClock::Stopwatch {
	public :
		Stopwatch() { start(); stop(); }
		~Stopwatch() {}

		Time & start() { fTime.first  = Clock::now(); return fTime.first;  }
		Time & stop()  { fTime.second = Clock::now(); return fTime.second; }
		
		Duration duration() { return (fTime.second - fTime.first); }
		double   time()     { return std::chrono::duration<double>( (fTime.second - fTime.first) ).count(); }

		std::ostream& print(std::ostream& out = std::cout, MGClock::ClockType type = MGClock::ClockType::UTC);

	protected :
		std::pair<Time, Time> fTime;
};


#endif // __MGClock_H__
