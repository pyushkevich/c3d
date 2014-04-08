/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    GeneralLinearModel.h
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

#ifndef __GeneralLinearModel_h_
#define __GeneralLinearModel_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class GeneralLinearModel : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  GeneralLinearModel(Converter *c) : c(c) {}

  void operator() (string fn_matrix, string fn_contrast);

private:
  Converter *c;

};

#endif

