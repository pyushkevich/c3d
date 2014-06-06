/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    PeronaMalik.cxx
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

#include "PeronaMalik.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

template <class TPixel, unsigned int VDim>
void
PeronaMalik<TPixel, VDim>
::operator() (double conductance, size_t nIter)
{
  // Get an image off the stack
  ImagePointer image = c->m_ImageStack.back();

  // Create a filter
  typedef itk::GradientAnisotropicDiffusionImageFilter<
    ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  *c->verbose << "Performing anisotropic diffusion on #" << c->m_ImageStack.size() << endl;
  filter->SetInput(image);
  filter->SetConductanceParameter(conductance);
  filter->SetNumberOfIterations(nIter);
  filter->SetTimeStep(0.0125);
  filter->UseImageSpacingOn();
  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class PeronaMalik<double, 2>;
template class PeronaMalik<double, 3>;
template class PeronaMalik<double, 4>;
