/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    MeanFilter.cxx
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

#include "MeanFilter.h"
#include <itkMeanImageFilter.h>

template <class TPixel, unsigned int VDim>
void
MeanFilter<TPixel, VDim>
::operator() (const SizeType &radius)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Do some processing ...
  typedef itk::MeanImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::RadiusType rad;

  for(int i = 0; i < VDim; i++)
    rad[i] = radius[i];

  *c->verbose << "Applying mean filter with radius " << rad << " to #" << c->m_ImageStack.size() << std::endl;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetRadius(rad);
  filter->SetInput(img);
  filter->Update();

  // ImagePointer result = ...;
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class MeanFilter<double, 2>;
template class MeanFilter<double, 3>;
template class MeanFilter<double, 4>;
