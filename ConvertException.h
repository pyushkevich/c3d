#ifndef __ConvertException_h_
#define __ConvertException_h_

#include <exception>
#include <string>
#include <cstdarg>

class ConvertException : public std::exception
{
public:
  ConvertException(const char *fmt, ...)
    {
    char buffer[1024];
    va_list parg;
    va_start(parg, fmt);
    vsprintf(buffer, fmt, parg);
    va_end(parg);
    message=buffer;
    }

  virtual ~ConvertException() throw() {}

  virtual const char *what() const throw()
    { return message.c_str(); }

private:
  std::string message;
};

#endif // __ConvertException_h_
