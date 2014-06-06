/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    WeightedSumVoxelwise.cxx
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

#include "WeightedSumVoxelwise.h"

template <class TPixel, unsigned int VDim>
void
WeightedSumVoxelwise<TPixel, VDim>
::operator() ()
{
  // Check input validity 
  if((c->m_ImageStack.size() % 2) == 1)
    throw ConvertException("Number of images on stack incompatible with"
      " weighted averaging (must be even");

  if(!c->CheckStackSameDimensions(0))
    throw ConvertException("Images on the stack don't have same dimensions");

  // Make a copy of the last image on the stack
  ImagePointer isum = ImageType::New();
  isum->SetRegions(c->m_ImageStack.back()->GetBufferedRegion());
  isum->CopyInformation(c->m_ImageStack.back());
  isum->Allocate();
  isum->FillBuffer(0);

  *c->verbose << "Taking weighted sum of " << c->m_ImageStack.size() / 2 << " images.";

  // Compute the weighted sum, two images at a time
  while(c->m_ImageStack.size())
    {
    // Get images from the stack
    ImagePointer img = c->m_ImageStack.back(); 
    c->m_ImageStack.pop_back();
    ImagePointer wgt = c->m_ImageStack.back(); 
    c->m_ImageStack.pop_back();

    // Compute the weighted sum
    Iterator is(img, img->GetBufferedRegion());
    Iterator iw(wgt, wgt->GetBufferedRegion());
    Iterator it(isum, isum->GetBufferedRegion());
    for(; !is.IsAtEnd(); ++is, ++it, ++iw)
      it.Set(it.Get() + is.Get() * iw.Get());
    }

  // Put result on stack
  c->m_ImageStack.push_back(isum);
}

// Invocations
template class WeightedSumVoxelwise<double, 2>;
template class WeightedSumVoxelwise<double, 3>;
template class WeightedSumVoxelwise<double, 4>;
