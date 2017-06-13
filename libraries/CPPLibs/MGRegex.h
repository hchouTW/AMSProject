#ifndef __CPPLibs_MGRegex_H__
#define __CPPLibs_MGRegex_H__

#include <string>
#include <regex>

namespace MGRegex {

namespace Formula {

constexpr const char * kIntegral = "^([-+]?(0|[1-9][0-9]*))$";
constexpr const char * kFloat    = "^(([-+])?(\\d+(\\.)?\\d*|\\d*(\\.)?\\d+)(([eE]([-+])?)?\\d+)?)$";
constexpr const char * kSpace    = "(\\s+)";
constexpr const char * kNotSpace = "(\\S+)";
constexpr const char * kPeriod    = "(\\.+)";
constexpr const char * kComma     = "(,+)";
constexpr const char * kColon     = "(:+)";
constexpr const char * kSemicolon = "(;+)";
constexpr const char * kSlash     = "(/+)";
constexpr const char * kBackslash = "(\\+)";

} // namespace Formula

inline void TrimSelf(std::string& str);
inline std::string Trim(const std::string& str);

inline void EraseSelf(std::string& str, const char target = ' ');
inline std::string Erase(const std::string& str, const char target = ' ');

inline void ReplaceSelf(std::string& str, const std::string& expr = Formula::kSpace, const std::string& fmt = std::string());
inline std::string Replace(const std::string& str, const std::string& expr = Formula::kSpace, const std::string& fmt = std::string());

inline std::vector<std::string> Split(const std::string& str, const std::string& expr = Formula::kSpace);
inline std::vector<std::string> Match(const std::string& str, const std::string& expr = Formula::kNotSpace);

inline std::string StringIntegral(const std::string& str); 
inline bool IsIntegral(const std::string& str) { return !(StringIntegral(str).empty()); }

template <class IntType = long long, typename std::enable_if<std::is_integral<IntType>::value, int>::type = 0>
inline std::pair<bool, IntType> ConvertFromStringToIntegral(const std::string& str);

inline std::string StringFloat(const std::string& str); 
inline bool IsFloat(const std::string& str) { return !(StringFloat(str).empty()); }

template <class RealType = long double, typename std::enable_if<std::is_floating_point<RealType>::value, int>::type = 0>
inline std::pair<bool, RealType> ConvertFromStringToFloat(const std::string& str);

} // namespace MGRegex


#endif // __CPPLibs_MGRegex_H__
