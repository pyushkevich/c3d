/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConvertException.h
  Language:  C++
  Website:   itksnap.org/c3d
  Copyright (c) 2014 Paul A. Yushkevich
  
  This file is part of C3D, a command-line companion tool to ITK-SNAP

  C3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

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

class StackAccessException : public ConvertException
{
public:
  StackAccessException()
  : ConvertException("Image Stack Access Exception") {}

};

#endif // __ConvertException_h_
