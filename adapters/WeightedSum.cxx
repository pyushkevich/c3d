/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    WeightedSum.cxx
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

#include "WeightedSum.h"

template <class TPixel, unsigned int VDim>
void
WeightedSum<TPixel, VDim>
::operator() (std::vector<double> weights)
{
  // The vector size must not exceed number of images
  size_t nw = weights.size();
  if(nw > c->m_ImageStack.size())
    throw ConvertException("Too many weights specified for -weighted-sum command");
  if(nw == 0)
    throw ConvertException("No weights specified for -weighted-sum command");

  ImagePointer iref = c->m_ImageStack.back();

  // Check dimensions
  if(!c->CheckStackSameDimensions(nw))
    throw ConvertException("Images on the stack don't have same dimensions");

  *c->verbose << "Taking weighted sum of images #" 
    << c->m_ImageStack.size() - nw << " to #" << c->m_ImageStack.size() << endl;

  // Scale the last image
  size_t nvox = iref->GetBufferedRegion().GetNumberOfPixels(); 
  double w = weights[nw-1];
  *c->verbose << "  Scaling #" << c->m_ImageStack.size() << " by " << w << endl;
  TPixel *ptrg = c->m_ImageStack[nw - 1]->GetBufferPointer();
  for(size_t i = 0; i < nvox; i++)
    ptrg[i] *= w;

  // Add the rest of the images
  for(size_t i = 1; i < nw; i++)
    {
    double w = weights[nw-(i+1)];
    *c->verbose << "  Scaling #" << c->m_ImageStack.size()-i << " by " << w << endl;
    TPixel *p = c->m_ImageStack[nw - (i+1)]->GetBufferPointer();
    for(size_t i = 0; i < nvox; i++)
      ptrg[i] += w * p[i];
    }

  // Drop all but the last image
  for(size_t i = 0; i < nw; i++)
    c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(iref);
}

// Invocations
template class WeightedSum<double, 2>;
template class WeightedSum<double, 3>;
template class WeightedSum<double, 4>;
