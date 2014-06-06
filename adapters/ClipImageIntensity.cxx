/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ClipImageIntensity.cxx
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

#include "ClipImageIntensity.h"

template <class TPixel, unsigned int VDim>
void
ClipImageIntensity<TPixel, VDim>
::operator() (double iMin, double iMax)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  *c->verbose << "Clipping out-of-range intensities in #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Intensity range: " << iMin << " to " << iMax << endl;

  // Simply replace values outside the clip range with new values
  for(Iterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(it.Value() < iMin)
      it.Set(iMin);
    else if(it.Value() > iMax)
      it.Set(iMax);
    }

  // Update the image
  img->Modified();
}

// Invocations
template class ClipImageIntensity<double, 2>;
template class ClipImageIntensity<double, 3>;
template class ClipImageIntensity<double, 4>;
