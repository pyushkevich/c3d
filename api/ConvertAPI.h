/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConvertImageND.cxx
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
#ifndef __ConvertAPI_h_
#define __ConvertAPI_h_

#include <iostream>
#include <string>

namespace itk {
  template <class TPixel, unsigned int VDim> class Image;
}

template <class TPixel, unsigned int VDim> class ImageConverter;

/** 
 * This header file provides the C++ interface to the c3d API. 
 *
 * It allows a third-party program to call C3D in its own process space, and
 * to pass images to C3D through memory (rather than through files). In the 
 * future, it will also support adding custom commands from the third party 
 * software.
 *
 * The third-party program should be able to interface with C3D by just including
 * this header file and linking with libraries cnd_api, cnd_driver and cnd_adapters
 */
template <class TPixel, unsigned int VDim>
class ConvertAPI
{
public:

  /** Standard ITK image - used to pass data in and out of the API */
  typedef itk::Image<TPixel, VDim> ImageType;

  ConvertAPI();
  ~ConvertAPI();

  /**
   * Pass an image as a variable to C3D. You can then access that variable using
   * the '-push' command during the call to Execute. This class will hold on to this
   * image until it is destroyed.
   */
 void AddImage(const char *varname, ImageType *image); 

  /** 
   * Execute a command in c3d. This is just like a regular command-line, and the 
   * output will be captured to the provided stream.
   *
   * Return Value: if an exception is caught during command execution, this will
   * return false, and you can access the error text using GetError()
   */ 
  bool Execute(const char *command, ...);

  /**
   * Get the exception string 
   */
  std::string GetError() const { return m_Error; }

  /**
   * Get a variable created by C3D during command execution, for example using
   * '-as' or '-popas' commands
   */
  ImageType *GetImage(const char *varname);

private:

  // Private instance of the image converter - not directly exposed through the API
  typedef ImageConverter<TPixel, VDim> ConverterType;
  ConverterType *m_Converter;

  // Error text from last execution
  std::string m_Error;
};



#endif // __ConvertAPI_h_
