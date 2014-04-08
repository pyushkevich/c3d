/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConvertAdapter.h
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

#ifndef __ConvertAdapter_h_
#define __ConvertAdapter_h_

#include "ConvertImageND.h"

// Common typedefs for all child classes
#define CONVERTER_STANDARD_TYPEDEFS \
  typedef ImageConverter<TPixel, VDim> Converter; \
  typedef typename Converter::ImageType ImageType; \
  typedef typename Converter::UnorientedImageType UnorientedImageType; \
  typedef typename Converter::ImagePointer ImagePointer; \
  typedef typename Converter::SizeType SizeType; \
  typedef typename Converter::IndexType IndexType; \
  typedef typename Converter::IntegerVector IntegerVector; \
  typedef typename Converter::RealVector RealVector; \
  typedef typename Converter::RegionType RegionType; \
  typedef typename Converter::ComplexPixel ComplexPixel; \
  typedef typename Converter::ComplexImageType ComplexImageType; \
  typedef typename Converter::Iterator Iterator; \
  typedef typename Converter::ConstIterator ConstIterator;

/**
 * Parent class for all adapters
 */
template<class TPixel, unsigned int VDim>
class ConvertAdapter
{
public:

  CONVERTER_STANDARD_TYPEDEFS

  ConvertAdapter() {}

protected:

};

#endif
