/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SmoothImage.cxx
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

#include "SmoothImage.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

template <class TPixel, unsigned int VDim>
void
SmoothImage<TPixel, VDim>
::operator() (RealVector &stdev, bool do_recursive)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Which filter to use 
  if(do_recursive)
    {
    // Describe what we are doing
    *c->verbose << "Fast recursive smoothing #" << c->m_ImageStack.size() << " with std.dev. " << stdev << endl;

    // Create a smoothing filter
    typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer fltSmooth = FilterType::New();
    typename FilterType::SigmaArrayType sarray;

    for(int i = 0; i < VDim; i++)
      sarray[i] = stdev[i];

    fltSmooth->SetInput(input);
    fltSmooth->SetSigmaArray(sarray);
    fltSmooth->Update();

    // Save the output
    c->m_ImageStack.pop_back();
    c->m_ImageStack.push_back(fltSmooth->GetOutput());
    }
  else
    {
    // Describe what we are doing
    *c->verbose << "Smoothing #" << c->m_ImageStack.size() << " with std.dev. " << stdev << endl;

    // Create a smoothing kernel and use it
    typedef itk::DiscreteGaussianImageFilter<ImageType,ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    typename FilterType::ArrayType variance;

    for(size_t i = 0; i < VDim; i++)
      variance[i] = stdev[i] * stdev[i];
    
    filter->SetInput(input);
    filter->SetVariance(variance);
    filter->SetUseImageSpacingOn();
    filter->Update();

    // Save the output
    c->m_ImageStack.pop_back();
    c->m_ImageStack.push_back(filter->GetOutput());
    }
}

// Invocations
template class SmoothImage<double, 2>;
template class SmoothImage<double, 3>;
template class SmoothImage<double, 4>;

