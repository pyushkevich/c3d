/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    WrapDimension.cxx
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

#include "WrapDimension.h"
#include "itkCyclicShiftImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "vnl/vnl_finite.h"

template <class TPixel, unsigned int VDim>
void
WrapDimension<TPixel, VDim>
::operator() (const IndexType &xWrap)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  typedef itk::CyclicShiftImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename FilterType::OffsetType shift;

  for(int i = 0; i < VDim; i++)
    shift[i] = xWrap[i];

  filter->SetShift(shift);
  filter->SetInput(img);
  filter->UpdateLargestPossibleRegion();

  // Explain what we are doing
  *c->verbose << "Wrapping image #" << c->m_ImageStack.size() << " by " << xWrap << endl;

  // Save the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class WrapDimension<double, 2>;
template class WrapDimension<double, 3>;
template class WrapDimension<double, 4>;
