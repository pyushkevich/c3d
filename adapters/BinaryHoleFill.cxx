/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    BinaryHoleFill.cxx
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

#include "BinaryHoleFill.h"
#include "itkBinaryFillholeImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BinaryHoleFill<TPixel, VDim>
::operator() (double foreground, bool full_conn)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Apply the fill hole algorithm
  typedef itk::BinaryFillholeImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(img);
  filter->SetForegroundValue(foreground);
  filter->SetFullyConnected(full_conn);

  // Some debug message
  *c->verbose << "Performing binary hole fill for intensity value " << foreground 
    << " in # " << c->m_ImageStack.size() << endl;

  filter->Update();
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class BinaryHoleFill<double, 2>;
template class BinaryHoleFill<double, 3>;
template class BinaryHoleFill<double, 4>;
