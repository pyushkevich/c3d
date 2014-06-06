/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ImageERF.cxx
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

#include "ImageERF.h"
#include "vnl/vnl_erf.h"

template <class TPixel, unsigned int VDim>
void
ImageERF<TPixel, VDim>
::operator() (double thresh, double scale)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Use iterator to apply erf
  itk::ImageRegionIteratorWithIndex<ImageType> 
    it(input, input->GetBufferedRegion());
  for(; !it.IsAtEnd(); ++it)
    {
    double x = it.Value();
    double y = vnl_erf((x - thresh) / scale);
    it.Set(y);
    }

  // Say what we are doing
  *c->verbose << "Taking ERF of #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  y = erf((x - " << thresh << ") / scale)" << endl;

  // Updated
  input->Modified();
}

// Invocations
template class ImageERF<double, 2>;
template class ImageERF<double, 3>;
template class ImageERF<double, 4>;
