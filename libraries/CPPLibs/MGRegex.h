#ifndef __MGRegex_H__
#define __MGRegex_H__


#include <string>
#include <regex>

#include "MGSys.h"


namespace MGRegex {
	namespace Formula {
		static constexpr const char * Integral = "^([-+]?(0|[1-9][0-9]*))$";
		static constexpr const char * Float    = "^(([-+])?(\\d+(\\.)?\\d*|\\d*(\\.)?\\d+)(([eE]([-+])?)?\\d+)?)$";
		static constexpr const char * Space    = "(\\s+)";
		static constexpr const char * NotSpace = "(\\S+)";
		static constexpr const char * Period    = "(\\.+)";
		static constexpr const char * Comma     = "(,+)";
		static constexpr const char * Colon     = "(:+)";
		static constexpr const char * Semicolon = "(;+)";
		static constexpr const char * Slash     = "(/+)";
		static constexpr const char * Backslash = "(\\+)";
	}

	void TrimSelf(std::string& str);
	std::string Trim(const std::string& str);

	void EraseSelf(std::string& str, const char target = ' ');
	std::string Erase(const std::string& str, const char target = ' ');

	void ReplaceSelf(std::string& str, const std::string& expr = Formula::Space, const std::string& fmt = std::string());
	std::string Replace(const std::string& str, const std::string& expr = Formula::Space, const std::string& fmt = std::string());
	
	std::vector<std::string> Split(const std::string& str, const std::string& expr = Formula::Space);
	std::vector<std::string> Match(const std::string& str, const std::string& expr = Formula::NotSpace);

	inline std::string StringIntegral(const std::string& str); 
	inline bool IsIntegral(const std::string& str) { return !(MGRegex::StringIntegral(str).empty()); }
	
	template <class IntType = long long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	std::pair<bool, IntType> ConvertFromStringToIntegral(const std::string& str);
	
	inline std::string StringFloat(const std::string& str); 
	inline bool IsFloat(const std::string& str) { return !(MGRegex::StringFloat(str).empty()); }
	
	template <class RealType = long double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	std::pair<bool, RealType> ConvertFromStringToFloat(const std::string& str);
}


#endif // __MGRegex_H__
