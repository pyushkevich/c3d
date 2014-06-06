/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ScaleShiftImage.cxx
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

#include "ScaleShiftImage.h"
#include "itkShiftScaleImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ScaleShiftImage<TPixel, VDim>
::operator() (double a, double b)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();
  
  // Say what we are doing
  *c->verbose << "Scaling #" << c->m_ImageStack.size() << " by " << a << " and adding " << b << endl;

  // If a=0, this means setting the image to a constant
  if(a == 0.0)
    {
    c->CopyImage();
    c->m_ImageStack.back()->FillBuffer(b);
    return;
    }

  // Create and run filter
  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(input);
  filter->SetScale(a);
  filter->SetShift(b / a);
  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class ScaleShiftImage<double, 2>;
template class ScaleShiftImage<double, 3>;
template class ScaleShiftImage<double, 4>;
