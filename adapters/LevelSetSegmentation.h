/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    LevelSetSegmentation.h
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

#ifndef __LevelSetSegmentation_h_
#define __LevelSetSegmentation_h_

#include "ConvertAdapter.h"

struct LevelSetParameters
{
  double CurvatureWeight;
  double AdvectionWeight;

  LevelSetParameters()
    {
    CurvatureWeight = 0.2;
    AdvectionWeight = 0.0;
    }
};

template<class TPixel, unsigned int VDim>
class LevelSetSegmentation : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  LevelSetSegmentation(Converter *c) : c(c) {}

  void operator() (int nIter, const LevelSetParameters &param);

private:
  Converter *c;

};

#endif

