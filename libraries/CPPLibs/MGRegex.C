#ifndef __MGRegex_C__
#define __MGRegex_C__

#include "MGRegex.h"


void MGRegex::TrimSelf(std::string& str) {
	const unsigned long first = str.find_first_not_of(' ');
	const unsigned long last  = str.find_last_not_of(' ');
	if (first == std::string::npos || last == std::string::npos) str.clear();
	else str.assign((str.begin() + first), (str.begin() + last + 1));
}


std::string MGRegex::Trim(const std::string& str) {
	std::string newStr = str;
	MGRegex::TrimSelf(newStr);
	return newStr;
}


void MGRegex::EraseSelf(std::string& str, const char target) {
	str.erase(std::remove_if(str.begin(), str.end(), ([&](const char& ch)->bool{return (target==ch);}) ), str.end());
}


std::string MGRegex::Erase(const std::string& str, const char target) {
	std::string newStr = str;
	EraseSelf(newStr, target);
	return newStr;
}

void MGRegex::ReplaceSelf(std::string& str, const std::string& expr, const std::string& fmt) {
	try {
		str = std::regex_replace(str, std::regex(expr), fmt);
	}
	catch (const std::regex_error& err) { 
		MGSys::ShowError(StrFmt("<< MGRegex::ReplaceSelf >>  %s", err.what()),
		MGSys::MESSAGES({
			StrFmt("SOURCE : \"%s\"", str.c_str()), 
			StrFmt("REGEX  : \"%s\"", expr.c_str())
		}));
	}
}

std::string MGRegex::Replace(const std::string& str, const std::string& expr, const std::string& fmt) {
	std::string newStr = str;
	MGRegex::ReplaceSelf(newStr, expr, fmt);
	return newStr;
}


std::vector<std::string> MGRegex::Split(const std::string& str, const std::string& expr) {
	const std::string&& trimStr = MGRegex::Trim(str);
	if (trimStr.empty()) return std::vector<std::string>();
	try {
		std::vector<std::string>&& strvec = { std::sregex_token_iterator(trimStr.cbegin(), trimStr.cend(), std::regex(expr), -1), std::sregex_token_iterator() };
		//strvec.erase(std::remove_if(strvec.begin(), strvec.end(), ([](std::string& str)->bool{return str.empty();}) ), strvec.end());
		return strvec;
	}
	catch (const std::regex_error& err) { 
		MGSys::ShowError(StrFmt("<< MGRegex::Split >>  %s", err.what()),
		MGSys::MESSAGES({
			StrFmt("SOURCE : \"%s\"", str.c_str()), 
			StrFmt("REGEX  : \"%s\"", expr.c_str())
		}));
		return std::vector<std::string>(); 
	}
}


std::vector<std::string> MGRegex::Match(const std::string& str, const std::string& expr) {
	const std::string&& trimStr = MGRegex::Trim(str);
	if (trimStr.empty()) return std::vector<std::string>();
	try { 
		return { std::sregex_token_iterator(trimStr.cbegin(), trimStr.cend(), std::regex(expr)), std::sregex_token_iterator() };
	}	
	catch (const std::regex_error& err) { 
		MGSys::ShowError(StrFmt("<< MGRegex::Match >>  %s", err.what()),
		MGSys::MESSAGES({
			StrFmt("SOURCE : \"%s\"", str.c_str()), 
			StrFmt("REGEX  : \"%s\"", expr.c_str())
		}));
		return std::vector<std::string>(); 
	}
}


inline std::string MGRegex::StringIntegral(const std::string& str) { 
	std::vector<std::string>&& strvec = MGRegex::Match(str, MGRegex::Formula::Integral);
	if (strvec.empty()) return std::string();
	else return strvec.at(0);
}


template <class IntType = long long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
std::pair<bool, IntType> ConvertFromStringToIntegral(const std::string& str) {
	std::vector<std::string>&& strvec = MGRegex::Match(str, MGRegex::Formula::Integral);
	if (strvec.empty()) return std::make_pair(false, IntType(0));
	else {
		try {
			long long tmpval = std::stoll(strvec.at(0));
			IntType value = static_cast<IntType>(tmpval);
			if (value != tmpval) throw std::out_of_range(CStrFmt("static_cast<%s>(%s) failure", typeid(IntType).name(), typeid(long long).name()));
			return std::make_pair(true, value);
		}
		catch (const std::logic_error & err) {
			MGSys::ShowError(StrFmt("<< MGRegex::ConvertFromStringToIntegral >>  %s", err.what()));
			return std::make_pair(false, IntType(0));
		}
	}
}

inline std::string MGRegex::StringFloat(const std::string& str) { 
	std::vector<std::string>&& strvec = MGRegex::Match(str, MGRegex::Formula::Float);
	if (strvec.empty()) return std::string();
	else return strvec.at(0);
}


template <class RealType = long double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
std::pair<bool, RealType> MGRegex::ConvertFromStringToFloat(const std::string& str) {
	std::vector<std::string>&& strvec = MGRegex::Match(str, MGRegex::Formula::Float);
	if (strvec.empty()) return std::make_pair(false, RealType(0.0));
	else {
		try {
			long double tmpval = std::stold(strvec.at(0));
			RealType value = static_cast<RealType>(tmpval);
			if (!std::isfinite(tmpval) || !std::isfinite(value)) throw std::out_of_range("is not finite");
			return std::make_pair(true, value);
		}
		catch (const std::logic_error & err) {
			MGSys::ShowError(StrFmt("<< MGRegex::ConvertFromStringToFloat >>  %s", err.what()));
			return std::make_pair(false, RealType(0.0));
		}
	}
}


#endif // __MGRegex_C__
