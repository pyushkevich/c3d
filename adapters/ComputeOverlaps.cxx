/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ComputeOverlaps.cxx
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

#include "ComputeOverlaps.h"

template <class TPixel, unsigned int VDim>
void
ComputeOverlaps<TPixel, VDim>
::operator() (double v)
{
  // There must be two images available on the stack!
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Overlap requires two images on the stack!" << endl;
    throw -1;
    }

  // Get the last two images
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Report what the filter is doing
  *c->verbose << "Computing overlap #" << c->m_ImageStack.size() - 1 
    << " and #" << c->m_ImageStack.size() << endl;

  // The images must have the same size
  if(i1->GetBufferedRegion() != i2->GetBufferedRegion())
    {
    cerr << "Overlap requires the images to be of the same dimensions!" << endl;
    throw -1;
    }

  // Create iterators for the two images
  ConstIterator it1(i1, i1->GetBufferedRegion());
  ConstIterator it2(i2, i2->GetBufferedRegion());
  double eps = 0.000001;

  // Counters of the overlap scores
  long n1 = 0, n2 = 0, n12 = 0;

  // Iterate over all pixels
  for(; !it1.IsAtEnd(); ++it1, ++it2)
    {
    // Read the values
    double v1 = it1.Get();
    double v2 = it2.Get();

    // Compare the values to target (within machine error)
    bool m1 = (v1 == v) || (fabs(2 * (v1 - v) / (v1 + v)) < eps);
    bool m2 = (v2 == v) || (fabs(2 * (v2 - v) / (v2 + v)) < eps);
    if(m1) n1++;
    if(m2) n2++;
    if(m1 && m2) n12++;
    }

  // Report the overlaps on one line
  double xDice = n12 * 2.0 / (n1 + n2);
  double xRobust = n12 * 1.0 / (n1 + n2 - n12);

  cout << "OVL: " << v << ", " << n1 << ", " << n2 << ", " << n12;
  cout << ", " << xDice << ", " << xRobust << endl;
  
  // Print the overlap to c->verbose channel
  *c->verbose << "  Matching voxels in first image:  " << n1 << endl;
  *c->verbose << "  Matching voxels in second image: " << n2 << endl;
  *c->verbose << "  Size of overlap region:          " << n12 << endl;
  *c->verbose << "  Dice similarity coefficient:     " << xDice << endl;
  *c->verbose << "  Intersection / ratio:            " << xRobust << endl;

}

// Invocations
template class ComputeOverlaps<double, 2>;
template class ComputeOverlaps<double, 3>;
template class ComputeOverlaps<double, 4>;
