/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ImageLaplacian.cxx
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

#include "ImageLaplacian.h"
#include "itkLaplacianImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ImageLaplacian<TPixel, VDim>
::operator() ()
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();
  
  // Describe what we are doing
  *c->verbose << "Taking Laplacian of #" << c->m_ImageStack.size() << endl;

  // Create a smoothing kernel and use it
  typedef itk::LaplacianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(input);
  filter->UseImageSpacingOn();
  filter->Update();

  // Save the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class ImageLaplacian<double, 2>;
template class ImageLaplacian<double, 3>;
template class ImageLaplacian<double, 4>;
