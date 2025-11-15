#ifndef __sto_h_
#define __sto_h_

#include "ConvertException.h"

#include <string>

template <typename T>
T strto(const std::string str)
{
  const char *cstr = str.c_str();
  char *end = 0;

  T v;
  std::string type_str;
  if (is_same<T, float>::value)
  {
    v = std::strtof(cstr, &end);
    type_str = "float";
  }
  else if (is_same<T, double>::value)
  {
    v = std::strtod(cstr, &end);
    type_str = "double";
  }
  else if (is_same<T, long double>::value)
  {
    v = std::strtold(cstr, &end);
    type_str = "long double";
  }
  else if (is_same<T, long>::value)
  {
    v = std::strtol(cstr, &end, 10);
    type_str = "long";
  }
  else if (is_same<T, long long>::value)
  {
    v = std::strtoll(cstr, &end, 10);
    type_str = "long long";
  }
  else if (is_same<T, unsigned long>::value)
  {
    v = std::strtoul(cstr, &end, 10);
    type_str = "unsigned long";
  }
  else if (is_same<T, unsigned long long>::value)
  {
    v = std::strtoull(cstr, &end, 10);
    type_str = "unsigned long long";
  }
  else
    throw ConvertException("unsupported type");

  if (*end != 0)
  {
    /*
      Because str is supplied by the user and can contain the % symbole it cannot
      be concatenated into the err_msg directly and then parsed on to the
      ConvertException's contructor with one parameter. If this is done, the user
      supplied % symbole would inside ConvertException be processed by vsnprintf
      as a fmt string, which uses % as control character, and hereby not be included
      into the ConvertException message. Instead the two argument ConvertException
      constructur is used to bypass this to report a message that includes a % sign.
    */
    std::string err_msg = "error converting %s";
    err_msg += " to a value of the type " + type_str;
    throw ConvertException(err_msg.c_str(), str.c_str());
  }
  return v;
};

#endif  // __sto_h_
