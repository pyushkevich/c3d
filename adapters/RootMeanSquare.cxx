/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    RootMeanSquare.cxx
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

#include "RootMeanSquare.h"

template <class TPixel, unsigned int VDim>
void
RootMeanSquare<TPixel, VDim>
::operator() ()
{

  // number of images on stack
  size_t nimages = c->m_ImageStack.size();
  // Check dimensions
  if(!c->CheckStackSameDimensions(nimages))
    throw ConvertException("Images on the stack don't have same dimensions");

  *c->verbose << "Taking voxelwise RMS of images #"
    << 0 << " to #" << nimages << endl;

  ImagePointer iref = c->m_ImageStack.back();

  // Scale the last image
  size_t nvox = iref->GetBufferedRegion().GetNumberOfPixels();
  *c->verbose << "  Squaring #" << c->m_ImageStack.size() << endl;
  TPixel *pref = iref->GetBufferPointer();
  for(size_t i = 0; i < nvox; i++)
    pref[i] *= pref[i];

  // Square and add the rest of the images
  for(size_t i = 1; i < nimages; i++)
    {
    *c->verbose << "  Squaring #" << c->m_ImageStack.size()-i << endl;
    TPixel *p = c->m_ImageStack[nimages - (i+1)]->GetBufferPointer();
    for(size_t i = 0; i < nvox; i++)
      pref[i] += p[i] * p[i];
    }

    *c->verbose << "Computing RMS" << endl;

    // now divide by number of images
    double denom = 1.0 / nimages;
    for(size_t i = 0; i < nvox; i++)
      pref[i] = std::sqrt(pref[i] * denom);

  // Drop all but the last image
  for(size_t i = 0; i < nimages; i++)
    c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(iref);
}

// Invocations
template class RootMeanSquare<double, 2>;
template class RootMeanSquare<double, 3>;
template class RootMeanSquare<double, 4>;
