/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    OverlayLabelImage.cxx
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

#include "OverlayLabelImage.h"

template <class TPixel, unsigned int VDim>
void
OverlayLabelImage<TPixel, VDim>
::operator() (const char *fnLUT, double opacity)
{
  // Check consistency
  if(!c->CheckStackSameDimensions(2))
    throw ConvertException(
      "Images have inconsistent dimensions for -overlay-label-image");

  // Get images off the stack
  size_t n = c->m_ImageStack.size();
  ImagePointer imLabel = c->m_ImageStack[n-1];
  ImagePointer imGray = c->m_ImageStack[n-2];

  // Load the lookup table
  typename ImageConverter<TPixel, VDim>::LabelToRGBAMap lmap = 
    c->ReadLabelToRGBAMap(fnLUT);

  // Explain what we are doing
  *c->verbose << "Overlaying labels in #" << (n-1) << " over image #" << (n-2) << endl;

  // Create output RGB images
  ImagePointer rgb[3];
  for(size_t d = 0; d < 3; d++)
    {
    // Allocate image
    rgb[d] = ImageType::New();
    rgb[d]->SetRegions(imGray->GetBufferedRegion());
    rgb[d]->CopyInformation(imGray);
    rgb[d]->Allocate();

    // Iterate
    ConstIterator ig(imGray, imGray->GetBufferedRegion());
    ConstIterator il(imLabel, imLabel->GetBufferedRegion());
    Iterator ic(rgb[d], rgb[d]->GetBufferedRegion());
    for(; !ig.IsAtEnd(); ++ig, ++il, ++ic)
      {
      // Find the color entry
      typename ImageConverter<TPixel, VDim>::LabelToRGBAMap::iterator ilut = 
        lmap.find(il.Get());
      vnl_vector_fixed<double, 4> color(0.0);
      if(ilut != lmap.end())
        color = ilut->second;

      // Blend (like OpenGL)
      double alpha = color[3] * opacity;
      ic.Set((1-alpha) * ig.Get() + alpha * color[d]);
      }
    }

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(rgb[0]);
  c->m_ImageStack.push_back(rgb[1]);
  c->m_ImageStack.push_back(rgb[2]);
}

// Invocations
template class OverlayLabelImage<double, 2>;
template class OverlayLabelImage<double, 3>;
template class OverlayLabelImage<double, 4>;
