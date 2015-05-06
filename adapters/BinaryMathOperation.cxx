/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    BinaryMathOperation.cxx
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

#include "BinaryMathOperation.h"

#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMinimumImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkAtan2ImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BinaryMathOperation<TPixel, VDim>
::operator() (Operation op)
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Binary operations require two images on the stack" << endl;
    throw -1;
    }

  // Get the last two images
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Write something
  *c->verbose << "Adding #" << c->m_ImageStack.size() - 1 << " and "  
    << c->m_ImageStack.size() - 2 << endl;

  // Select the operation
  typedef typename itk::ImageToImageFilter<ImageType, ImageType> BaseFilterType;
  typename BaseFilterType::Pointer filter;

  switch(op) 
    {
    case ADD:
      filter = itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case ATAN2:
      filter = itk::Atan2ImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case SUBTRACT:
      filter = itk::SubtractImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case MULTIPLY:
      filter = itk::MultiplyImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case DIVIDE:
      filter = itk::DivideImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case MINIMUM:
      filter = itk::MinimumImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case MAXIMUM:
      filter = itk::MaximumImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    }

  // Run the filter
  filter->SetInput(0, i1);
  filter->SetInput(1, i2);
  filter->Update();

  // Replace the images with the product
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class BinaryMathOperation<double, 2>;
template class BinaryMathOperation<double, 3>;
template class BinaryMathOperation<double, 4>;
