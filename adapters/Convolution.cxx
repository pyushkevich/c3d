/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    Convolution.cxx
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

#include "Convolution.h"
#include "itkFFTConvolutionImageFilter.h"

template <class TPixel, unsigned int VDim>
void
Convolution<TPixel, VDim>
::operator() ()
{
  // Get kernel from stack
  ImagePointer kernel = c->PopImage();

  // Get image from stack
  ImagePointer img = c->PopImage();

  // Create convolution filter
  typedef itk::FFTConvolutionImageFilter<ImageType> FilterType;
  typename FilterType::Pointer f = FilterType::New();
  f->SetInput(img);
  f->SetKernelImage(kernel);
  f->SetOutputRegionModeToSame();
  f->NormalizeOff();

  // Report what we are doing
  *c->verbose << "Performing convolution of #" << c->GetStackSize()+1
   << " with kernel #" << c->GetStackSize() << endl;

  f->Update();

  // The output image will have same information as the input image
  ImagePointer result = f->GetOutput();
  result->CopyInformation(img);

  // Push the image
  c->PushImage(result);
}

// Invocations
template class Convolution<double, 2>;
template class Convolution<double, 3>;
template class Convolution<double, 4>;
