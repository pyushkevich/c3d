/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    AntiAliasImage.cxx
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

#include "AntiAliasImage.h"
#include "itkAntiAliasBinaryImageFilter.h"

template <class TPixel, unsigned int VDim>
void
AntiAliasImage<TPixel, VDim>
::operator() (double xIsoSurface, double rms)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Report what the filter is doing
  *c->verbose << "Anti-aliasing #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Root Mean Square error: " << rms << endl;
  *c->verbose << "  Iterations: "; 
  if(c->m_Iterations == 0) 
    *c->verbose << "Unlimited" << endl; 
  else 
    *c->verbose << c->m_Iterations << endl;

  // Apply antialiasing to the image
  typedef itk::AntiAliasBinaryImageFilter<ImageType,ImageType> AntiFilterType;
  typename AntiFilterType::Pointer fltAnti = AntiFilterType::New();
  fltAnti->SetInput(input);
  fltAnti->SetMaximumRMSError(rms);
  if(c->m_Iterations > 0)
    fltAnti->SetNumberOfIterations(c->m_Iterations);
  fltAnti->SetIsoSurfaceValue(xIsoSurface);
  // fltAnti->AddObserver(itk::ProgressEvent(),command);
  fltAnti->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltAnti->GetOutput());
}

// Invocations
template class AntiAliasImage<double, 2>;
template class AntiAliasImage<double, 3>;
template class AntiAliasImage<double, 4>;
