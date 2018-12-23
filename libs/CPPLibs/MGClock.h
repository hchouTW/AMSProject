#ifndef __CPPLibs_MGClock_H__
#define __CPPLibs_MGClock_H__

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


namespace MGClock {

enum class ClockType { UTC = 0, LOCAL = 1 };

typedef std::time_t                                    UTime;
typedef struct std::tm                                 MTime;
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


inline UTime   ConvertFromMTimeToUTime(MTime * mtime) { return std::mktime(mtime); }
inline MTime * ConvertFromUTimeToMTime(UTime utime, ClockType type = ClockType::UTC); 
inline CTime   ConvertFromMTimeToCTime(MTime * mtime, const std::string& fmt = "");
inline CTime   ConvertFromUTimeToCTime(UTime utime, ClockType type = ClockType::UTC, const std::string& fmt = "") { return ConvertFromMTimeToCTime(ConvertFromUTimeToMTime(utime, type), fmt); }

inline SysTime ConvertFromUTime(UTime timpnt)  { return SysClock::from_time_t(timpnt); }
inline SysTime ConvertFromMTime(MTime * mtime) { return ConvertFromUTime(ConvertFromMTimeToUTime(mtime)); }
inline UTime   ConvertToUTime(SysTime timpnt) { return SysClock::to_time_t(timpnt); }
inline UTime   ConvertToUTime(StdTime timpnt) { return SysClock::to_time_t( (SysClock::now() + (timpnt - StdClock::now())) ); }

inline MTime * ConvertToMTime(SysTime timpnt, ClockType type = ClockType::UTC) { return ConvertFromUTimeToMTime(ConvertToUTime(timpnt), type); }
inline MTime * ConvertToMTime(StdTime timpnt, ClockType type = ClockType::UTC) { return ConvertFromUTimeToMTime(ConvertToUTime(timpnt), type); }

inline CTime   ConvertToCTime(SysTime timpnt, ClockType type = ClockType::UTC, const std::string& fmt = "") { return ConvertFromUTimeToCTime(ConvertToUTime(timpnt), type, fmt); }
inline CTime   ConvertToCTime(StdTime timpnt, ClockType type = ClockType::UTC, const std::string& fmt = "") { return ConvertFromUTimeToCTime(ConvertToUTime(timpnt), type, fmt); }

template<class Time>
std::string PrintStr(const Time& timpnt, ClockType type = ClockType::UTC, const std::string& fmt = "");

template<class Time>
std::ostream& Print(const Time& timpnt, ClockType type = ClockType::UTC, const std::string& fmt = "", std::ostream& out = std::cout);

template <class Clock, class Time, class Duration> 
class Stopwatch {
	public :
		Stopwatch() { start(); stop(); }
		~Stopwatch() {}

		Time& start() { times_.first  = Clock::now(); return times_.first;  }
		Time& stop()  { times_.second = Clock::now(); return times_.second; }
		
		Duration duration() const { return (times_.second - times_.first); }
		double   time()     const { return std::chrono::duration<double>( (times_.second - times_.first) ).count(); }

		std::ostream& print(MGClock::ClockType type = MGClock::ClockType::UTC, std::ostream& out = std::cout) const;

	protected :
		std::pair<Time, Time> times_;
};

typedef Stopwatch<SysClock, SysTime, SysDuration> SysStopwatch;
typedef Stopwatch<StdClock, StdTime, StdDuration> StdStopwatch;
typedef Stopwatch<HrsClock, HrsTime, HrsDuration> HrsStopwatch;

} // namespace MGClock


#endif // __CPPLibs_MGClock_H__
