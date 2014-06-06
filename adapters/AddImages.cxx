/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    AddImages.cxx
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

#include "AddImages.h"
#include "itkAddImageFilter.h"

template <class TPixel, unsigned int VDim>
void
AddImages<TPixel, VDim>
::operator() ()
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Binary operations require two images on the stack" << endl;
    throw -1;
    }

  // Get the last two images
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Write something
  *c->verbose << "Adding #" << c->m_ImageStack.size() - 1 << " and "  
    << c->m_ImageStack.size() - 2 << endl;

  // Perform the multiplication
  typedef itk::AddImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer flt = FilterType::New();
  flt->SetInput1(i1);
  flt->SetInput2(i2);
  flt->Update();

  // Replace the images with the product
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(flt->GetOutput());
}

// Invocations
template class AddImages<double, 2>;
template class AddImages<double, 3>;
template class AddImages<double, 4>;
