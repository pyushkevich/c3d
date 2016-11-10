/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    MatchBoundingBoxes.cxx
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

#include "MatchBoundingBoxes.h"

template <class TPixel, unsigned int VDim>
void
MatchBoundingBoxes<TPixel, VDim>
::operator() ()
{
  // Pop off the last image
  size_t n = c->m_ImageStack.size();
  if(n < 2)
    {
    throw string("Two images must be on the stack");
    }

  // Get data source and transform source
  ImagePointer d = c->m_ImageStack[n-1];
  ImagePointer t = c->m_ImageStack[n-2];

  // The header of the second image (spacing and origin) so that the
  // bounding boxes of the two images are exactly the same
  // O1 - S1 * D * [0.5] = O2 - S2 * D * [0.5]

  vnl_matrix<double> dir = t->GetDirection().GetVnlMatrix();
  vnl_vector<double> v(VDim); v.fill(0.5);
  vnl_vector<double> Dv = dir * v;
  vnl_vector<double> orig(VDim), spc(VDim);

  // Compute the spacing of the second image
  for(int i = 0; i < VDim; i++)
    {
    double w = t->GetSpacing()[i] * t->GetLargestPossibleRegion().GetSize(i);
    spc[i] = w / d->GetLargestPossibleRegion().GetSize(i);
    orig[i] = t->GetOrigin()[i] + (spc[i] - t->GetSpacing()[i]) * Dv[i]; 
    }

  // Print out what is being done
  *c->verbose << "Matching bounding box of #" << n-1 << " to that of #" << n-2 << endl;

  // Update the data source
  d->SetOrigin(orig.data_block());
  d->SetSpacing(spc.data_block());
  d->SetDirection(t->GetDirection());

  // Push and pop
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(d);
}

// Invocations
template class MatchBoundingBoxes<double, 2>;
template class MatchBoundingBoxes<double, 3>;
template class MatchBoundingBoxes<double, 4>;
