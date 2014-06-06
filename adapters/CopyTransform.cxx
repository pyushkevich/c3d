/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    CopyTransform.cxx
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

#include "CopyTransform.h"

template <class TPixel, unsigned int VDim>
void
CopyTransform<TPixel, VDim>
::operator() ()
{
  // Pop off the last image
  size_t n = c->m_ImageStack.size();
  if(n < 2)
    {
    throw string("Two images must be on the stack");
    }

  // Get data source and transform source
  ImagePointer dsrc = c->m_ImageStack[n-1];
  ImagePointer tsrc = c->m_ImageStack[n-2];

  // Make sure dimensions match
  if(dsrc->GetBufferedRegion().GetSize() != tsrc->GetBufferedRegion().GetSize())
    {
    throw string("Dimensions of images must match");
    }

  // Print out what is being done
  *c->verbose << "Copying transform from #" << n-2 << " to #" << n-1 << endl;

  // Update the data source
  dsrc->SetOrigin(tsrc->GetOrigin());
  dsrc->SetSpacing(tsrc->GetSpacing());
  dsrc->SetDirection(tsrc->GetDirection());

  // Push and pop
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(dsrc);
}

// Invocations
template class CopyTransform<double, 2>;
template class CopyTransform<double, 3>;
template class CopyTransform<double, 4>;
