/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SignedDistanceTransform.cxx
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

#include "SignedDistanceTransform.h"
#include "ThresholdImage.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

template <class TPixel, unsigned int VDim>
void
SignedDistanceTransform<TPixel, VDim>
::operator() ()
{
  // The image is assumed to be binary. If background is non-zero, call binarize
  // to map the background to zero
  if(c->m_Background != 0.0)
    {
    ThresholdImage<TPixel, VDim> thresh(c);
    thresh(c->m_Background, c->m_Background, 0.0, 1.0);
    }

  // Get the last image on the stack
  ImagePointer image = c->m_ImageStack.back();

  // Construct the connected components filter
  typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> Filter;

  // Describe what we are doing
  *c->verbose << "Computing signed distance function of #" << c->m_ImageStack.size() << endl;

  // Plug in the filter's components
  typename Filter::Pointer flt = Filter::New();
  flt->SetInput(image);
  flt->SetUseImageSpacing(true);
  flt->SquaredDistanceOff();
  flt->Update();

  // Store the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(flt->GetOutput());

}

// Invocations
template class SignedDistanceTransform<double, 2>;
template class SignedDistanceTransform<double, 3>;
template class SignedDistanceTransform<double, 4>;
