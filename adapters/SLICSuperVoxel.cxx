/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SLICSuperVoxel.cxx
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

#include "SLICSuperVoxel.h"
#include "itkSLICSuperVoxelImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkScalarConnectedComponentImageFilter.h"
#include "itkCastImageFilter.h"

template <class TPixel, unsigned int VDim>
void
SLICSuperVoxel<TPixel, VDim>
::operator() (int svPerDim, double m)
{
  // Get the image from the stack
  ImagePointer img = c->m_ImageStack.back();

  // Compute the gradient of the image
  typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradFilter;
  typename GradFilter::Pointer gradFilter = GradFilter::New();
  gradFilter->SetInput(img);
  gradFilter->Update();

  // Create the main filter
  typedef itk::SLICSuperVoxelImageFilter<ImageType, ImageType, ImageType> SLICFilter;
  typename SLICFilter::Pointer fltSlic = SLICFilter::New();
  fltSlic->SetInput(img);
  fltSlic->SetGradientImage(gradFilter->GetOutput());
  fltSlic->SetMParameter(m);
  fltSlic->SetSeedsPerDimension(svPerDim);
  fltSlic->Update();

  // Run a scalar connected components filter
  typedef itk::Image<short, VDim> IntegerImage;
  typedef itk::ScalarConnectedComponentImageFilter<ImageType, IntegerImage> CompFilter;
  typename CompFilter::Pointer fltComp = CompFilter::New();
  fltComp->SetInput(fltSlic->GetOutput());
  fltComp->SetDistanceThreshold(0.0);
  // fltComp->Update();

  typedef itk::CastImageFilter<IntegerImage, ImageType> CastFilter;
  typename CastFilter::Pointer fltCast = CastFilter::New();
  fltCast->SetInput(fltComp->GetOutput());
  // fltCast->Update();

  // Do some processing ...
  // ImagePointer result = ...;
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltSlic->GetOutput());
}

// Invocations
template class SLICSuperVoxel<double, 2>;
template class SLICSuperVoxel<double, 3>;
template class SLICSuperVoxel<double, 4>;
