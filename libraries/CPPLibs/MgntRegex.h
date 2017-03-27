#ifndef __MgntRegex_H__
#define __MgntRegex_H__
#include <string>
#include <regex>

#include "MgntSys.h"

namespace MgntRegex {
	namespace Formula {
		static constexpr const char * Integral = "^([-+]?(0|[1-9][0-9]*))$";
		static constexpr const char * Float    = "^(([-+])?(\\d+(\\.)?\\d*|\\d*(\\.)?\\d+)(([eE]([-+])?)?\\d+)?)$";
		static constexpr const char * Space    = "(\\s+)";
		static constexpr const char * NotSpace = "(\\S+)";
		static constexpr const char * Period    = "(\\.+)";
		static constexpr const char * Comma     = "(,+)";
		static constexpr const char * Colon     = "(:+)";
		static constexpr const char * Semicolon = "(;+)";
	}

	std::string Trim(const std::string& str) {
		const unsigned long first = str.find_first_not_of(' ');
		const unsigned long last  = str.find_last_not_of(' ');
		if (first == std::string::npos || last == std::string::npos) return std::string();
		return { (str.begin() + first), (str.begin() + last + 1) };
	}

	void TrimSelf(std::string& str) {
		const unsigned long limit = std::numeric_limits<unsigned long>::max();
		const unsigned long first = str.find_first_not_of(' ');
		const unsigned long last  = str.find_last_not_of(' ');
		if (first == std::string::npos || last == std::string::npos) str.clear();
		else str.assign((str.begin() + first), (str.begin() + last + 1));
	}

	inline std::string Erase(const std::string& str, const char target = ' ') {
		std::string newStr;
		std::copy_if(str.begin(), str.end(), newStr.begin(), ([&](const char& ch)->bool{return (target!=ch);}) );
		return newStr;
	}
	inline void EraseSelf(std::string& str, const char target = ' ') {
		str.erase(std::remove_if(str.begin(), str.end(), ([&](const char& ch)->bool{return (target==ch);}) ), str.end());
	}
	
	std::string Replace(const std::string& str, const std::string& expr = "(\\s+)", const std::string& fmt = std::string()) {
		try {
			return std::regex_replace(str, std::regex(expr), fmt);
		}
		catch (const std::regex_error& err) { 
			MgntSys::Error(StrFmt("<< MgntRegex::Replace >>  %s", err.what()),
			MgntSys::MESSAGES({
				StrFmt("SOURCE : \"%s\"", str.c_str()), 
				StrFmt("REGEX  : \"%s\"", expr.c_str())
			}));
			return str; 
		}
	}
	
	void ReplaceSelf(std::string& str, const std::string& expr = "(\\s+)", const std::string& fmt = std::string()) {
		try {
			str = std::regex_replace(str, std::regex(expr), fmt);
		}
		catch (const std::regex_error& err) { 
			MgntSys::Error(StrFmt("<< MgntRegex::ReplaceSelf >>  %s", err.what()),
			MgntSys::MESSAGES({
				StrFmt("SOURCE : \"%s\"", str.c_str()), 
				StrFmt("REGEX  : \"%s\"", expr.c_str())
			}));
		}
	}
	
	std::vector<std::string> Split(const std::string& str, const std::string& expr = "(\\s+)") {
		const std::string&& trimStr = MgntRegex::Trim(str);
		if (trimStr.empty()) return std::vector<std::string>();
		try {
			std::vector<std::string>&& strvec = { std::sregex_token_iterator(trimStr.cbegin(), trimStr.cend(), std::regex(expr), -1), std::sregex_token_iterator() };
			//strvec.erase(std::remove_if(strvec.begin(), strvec.end(), ([](std::string& str)->bool{return str.empty();}) ), strvec.end());
			return strvec;
		}
		catch (const std::regex_error& err) { 
			MgntSys::Error(StrFmt("<< MgntRegex::Split >>  %s", err.what()),
			MgntSys::MESSAGES({
				StrFmt("SOURCE : \"%s\"", str.c_str()), 
				StrFmt("REGEX  : \"%s\"", expr.c_str())
			}));
			return std::vector<std::string>(); 
		}
	}

	std::vector<std::string> Match(const std::string& str, const std::string& expr = "(\\S+)") {
		const std::string&& trimStr = MgntRegex::Trim(str);
		if (trimStr.empty()) return std::vector<std::string>();
		try { 
			return { std::sregex_token_iterator(trimStr.cbegin(), trimStr.cend(), std::regex(expr)), std::sregex_token_iterator() };
		}	
		catch (const std::regex_error& err) { 
			MgntSys::Error(StrFmt("<< MgntRegex::Match >>  %s", err.what()),
			MgntSys::MESSAGES({
				StrFmt("SOURCE : \"%s\"", str.c_str()), 
				StrFmt("REGEX  : \"%s\"", expr.c_str())
			}));
			return std::vector<std::string>(); 
		}
	}

	inline std::string StringIntegral(const std::string& str) { 
		std::vector<std::string>&& strvec = MgntRegex::Match(str, MgntRegex::Formula::Integral);
		if (strvec.empty()) return std::string();
		else return strvec.at(0);
	}
	inline bool IsIntegral(const std::string& str) { return !(MgntRegex::StringIntegral(str).empty()); }
	
	template <class IntType = long long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
	std::pair<bool, IntType> ConvertFromStringToIntegral(const std::string& str) {
		std::vector<std::string>&& strvec = MgntRegex::Match(str, MgntRegex::Formula::Integral);
		if (strvec.empty()) return std::make_pair(false, IntType(0));
		else {
			try {
				long long tmpval = std::stoll(strvec.at(0));
				IntType value = static_cast<IntType>(tmpval);
				if (value != tmpval) throw std::out_of_range(CStrFmt("static_cast<%s>(%s) failure", typeid(IntType).name(), typeid(long long).name()));
				return std::make_pair(true, value);
			}
			catch (const std::logic_error & err) {
				MgntSys::Error(StrFmt("<< MgntRegex::ConvertFromStringToIntegral >>  %s", err.what()));
				return std::make_pair(false, IntType(0));
			}
		}
	}
	
	inline std::string StringFloat(const std::string& str) { 
		std::vector<std::string>&& strvec = MgntRegex::Match(str, MgntRegex::Formula::Float);
		if (strvec.empty()) return std::string();
		else return strvec.at(0);
	}
	inline bool IsFloat(const std::string& str) { return !(MgntRegex::StringFloat(str).empty()); }
	
	template <class RealType = long double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
	std::pair<bool, RealType> ConvertFromStringToFloat(const std::string& str) {
		std::vector<std::string>&& strvec = MgntRegex::Match(str, MgntRegex::Formula::Float);
		if (strvec.empty()) return std::make_pair(false, RealType(0.0));
		else {
			try {
				long double tmpval = std::stold(strvec.at(0));
				RealType value = static_cast<RealType>(tmpval);
				if (!std::isfinite(tmpval) || !std::isfinite(value)) throw std::out_of_range("is not finite");
				return std::make_pair(true, value);
			}
			catch (const std::logic_error & err) {
				MgntSys::Error(StrFmt("<< MgntRegex::ConvertFromStringToFloat >>  %s", err.what()));
				return std::make_pair(false, RealType(0.0));
			}
		}
	}
}

#endif // __MgntRegex_H__
