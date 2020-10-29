/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    OtsuThreshold.cxx
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

#include "OtsuThreshold.h"
#include "itkOtsuThresholdImageFilter.h"

template <class TPixel, unsigned int VDim>
void
OtsuThreshold<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Set the inside/outside values to zero and one
  typedef itk::OtsuThresholdImageFilter<ImageType,ImageType> OtsuFilterType;
  typename OtsuFilterType::Pointer filter = OtsuFilterType::New();
  filter->SetInput(img);
  filter->SetOutsideValue(0.0);
  filter->SetInsideValue(1.0);

  *c->verbose << "Applying Otsu thresholding to #" << c->m_ImageStack.size() << endl;
  filter->Update();

  ImagePointer result = filter->GetOutput();

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class OtsuThreshold<double, 2>;
template class OtsuThreshold<double, 3>;
template class OtsuThreshold<double, 4>;
