/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    PadImage.cxx
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

#include "PadImage.h"
#include <itkConstantPadImageFilter.h>

template <class TPixel, unsigned int VDim>
void
PadImage<TPixel, VDim>
::operator() (IndexType padExtentLower, IndexType padExtentUpper, float padValue)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  typedef itk::ConstantPadImageFilter< ImageType, ImageType > ConstantPad;
  typename ConstantPad::Pointer padFilter = ConstantPad::New();

  // Pad first three dimensions only
  itk::SizeValueType lowerBound[VDim];
  itk::SizeValueType upperBound[VDim];

  for (int i = 0; i < 3 && i < VDim; i++) {
    lowerBound[i] = padExtentLower[i];
    upperBound[i] = padExtentUpper[i];
  }

  for (unsigned int i = 3; i < VDim; i++) {
    lowerBound[i] = 0;
    upperBound[i] = 0;
  }

  padFilter->SetPadLowerBound(lowerBound);
  padFilter->SetPadUpperBound(upperBound);

  padFilter->SetConstant(static_cast<typename ImageType::PixelType>(padValue));

  padFilter->SetInput(img);

  padFilter->Update();

  ImagePointer output = padFilter->GetOutput();

  // ITK 3.20 fixes the origin for you, don't mess with it
  //  cout << "INP_ORG: " << img->GetOrigin() << endl;
  //  cout << "OUT_ORG: " << output->GetOrigin() << endl;

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(output);

}

// Invocations
template class PadImage<double, 2>;
template class PadImage<double, 3>;
template class PadImage<double, 4>;
