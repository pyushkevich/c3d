/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ThresholdImage.cxx
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

#include "ThresholdImage.h"
#include "itkBinaryThresholdImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ThresholdImage<TPixel, VDim>
::operator() (double u1, double u2, double vIn, double vOut)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();
  
  // Say what we are doing
  *c->verbose << "Thresholding #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Mapping range [" << u1 << ", " << u2 << "] to " << vIn << endl;
  *c->verbose << "  Values outside are mapped to " << vOut << endl;

  // Do the thresholding
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(input);
  filter->SetLowerThreshold(u1);
  filter->SetUpperThreshold(u2);
  filter->SetInsideValue(vIn);
  filter->SetOutsideValue(vOut);
  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class ThresholdImage<double, 2>;
template class ThresholdImage<double, 3>;
template class ThresholdImage<double, 4>;
