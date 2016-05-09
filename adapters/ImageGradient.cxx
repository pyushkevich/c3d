/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ImageGradient.cxx
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

#include "ImageGradient.h"
#include "itkGradientImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkShiftScaleImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ImageGradient<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();

  // Define the gradient image filter
  typedef itk::GradientImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetUseImageSpacing(true);
  filter->SetUseImageDirection(true);
  filter->SetInput(img);
  filter->Update();

  // Report
  *c->verbose << "Taking the gradient of #" << c->m_ImageStack.size()
              << " (in physical space)" << std::endl;

  // Break into components 
  for(int i = 0; i < VDim; i++)
    {
    // Extract component
    typedef itk::VectorIndexSelectionCastImageFilter<
      typename FilterType::OutputImageType, ImageType> CompFilterType;

    typename CompFilterType::Pointer fltComp = CompFilterType::New();
    fltComp->SetInput(filter->GetOutput());
    fltComp->SetIndex(i);
    fltComp->Update();

    // Deal with RAS/LPS
    if(VDim >= 3 && i < 2)
      {
      typedef itk::ShiftScaleImageFilter<ImageType, ImageType> ScaleFilter;
      typename ScaleFilter::Pointer fltScale = ScaleFilter::New();
      fltScale->SetInput(fltComp->GetOutput());
      fltScale->SetScale(-1.0);
      fltScale->Update();
      c->m_ImageStack.push_back(fltScale->GetOutput());
      }
    else
      {
      c->m_ImageStack.push_back(fltComp->GetOutput());
      }
    }
}

// Invocations
template class ImageGradient<double, 2>;
template class ImageGradient<double, 3>;
template class ImageGradient<double, 4>;
