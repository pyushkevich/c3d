/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    UnaryMathOperation.cxx
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

#include "UnaryMathOperation.h"

template <class TPixel, unsigned int VDim>
void
UnaryMathOperation<TPixel, VDim>
::operator() (double (*func) (double))
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Print debug info
  *c->verbose << "Applying unary math operation to #" << c->m_ImageStack.size() << endl;

  // Apply the appropriate math operation (in place)
  Iterator it(img, img->GetBufferedRegion());
  for(; !it.IsAtEnd(); ++it)
    it.Set(func(it.Get()));
}

// Invocations
template class UnaryMathOperation<double, 2>;
template class UnaryMathOperation<double, 3>;
template class UnaryMathOperation<double, 4>;
