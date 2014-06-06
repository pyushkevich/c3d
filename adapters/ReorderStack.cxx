/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ReorderStack.cxx
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

#include "ReorderStack.h"

template <class TPixel, unsigned int VDim>
void
ReorderStack<TPixel, VDim>
::operator() (size_t k)
{
  // Get the number of images on the stack
  size_t n = c->m_ImageStack.size();

  // Make sure n is divisible by k
  if(n % k != 0)
    throw ConvertException("Can not reorder %d images using stride of %d;"
      " %d is not divisible by %d.", n, k, n, k);

  // Explain what we are doing
  *c->verbose << "Reordering " << n << " images with stride of " << k << endl;

  // Make a copy of the stack
  ImageStack<ImageType> temp_stack;
  for(size_t i=0; i < c->m_ImageStack.size(); i++)
  {
    temp_stack.push_back(c->m_ImageStack[i]);
  }

  c->m_ImageStack.clear();

  // Traverse the stack
  for(size_t i = 0; i < k; i++)
    for(size_t j = i; j < n; j += k)
      c->m_ImageStack.push_back(temp_stack[j]);
}

// Invocations
template class ReorderStack<double, 2>;
template class ReorderStack<double, 3>;
template class ReorderStack<double, 4>;
