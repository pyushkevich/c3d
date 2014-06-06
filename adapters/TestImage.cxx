/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    TestImage.cxx
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

#include "TestImage.h"
#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkMinimumMaximumImageFilter.h>

template <class TPixel, unsigned int VDim>
void
TestImage<TPixel, VDim>
::operator() (bool test_header, bool test_voxels, double tol) 
{
  // There must be two images on the stack
  if(c->m_ImageStack.size() < 2)
    throw ConvertException("Two images are requred for the test command");

  // Load the two images from the stack
  ImagePointer x = c->m_ImageStack.back(); c->m_ImageStack.pop_back();
  ImagePointer y = c->m_ImageStack.back(); c->m_ImageStack.pop_back();

  // Compare the headers
  if(test_header)
    {
    double max_abs_diff = 0;
    RegionType rx = x->GetLargestPossibleRegion(), ry = y->GetLargestPossibleRegion();
    for(int i = 0; i < VDim; i++)
      {
      max_abs_diff = vnl_math_max(max_abs_diff, vnl_math_abs((double) (rx.GetSize()[i] - ry.GetSize()[i])));
      max_abs_diff = vnl_math_max(max_abs_diff, vnl_math_abs((double) (rx.GetIndex()[i] - ry.GetIndex()[i])));
      max_abs_diff = vnl_math_max(max_abs_diff, vnl_math_abs((double) (x->GetOrigin()[i] - y->GetOrigin()[i])));
      max_abs_diff = vnl_math_max(max_abs_diff, vnl_math_abs((double) (x->GetSpacing()[i] - y->GetSpacing()[i])));
      for(int j = 0; j < VDim; j++)
        {
        max_abs_diff = vnl_math_max(max_abs_diff, vnl_math_abs((double) (x->GetDirection()[i][j] - y->GetDirection()[i][j])));
        }
      }

    if(max_abs_diff > tol)
      {
      std::cout << "Image header test failed. Max abs difference: " << max_abs_diff << endl;
      std::exit(1);
      }
    }

  if(test_voxels)
    {
    typedef typename itk::AbsoluteValueDifferenceImageFilter<ImageType, ImageType, ImageType> DiffFilter;
    typedef typename itk::MinimumMaximumImageFilter<ImageType> MinMaxFilter;
    typename DiffFilter::Pointer fltDiff = DiffFilter::New();
    fltDiff->SetInput(0, x);
    fltDiff->SetInput(1, y); 

    typename MinMaxFilter::Pointer fltMinMax = MinMaxFilter::New();
    fltMinMax->SetInput(fltDiff->GetOutput());
    fltMinMax->Update();

    double max_abs_diff = fltMinMax->GetMaximumOutput()->Get();
    if(max_abs_diff > tol)
      {
      std::cout << "Image voxel test failed. Max abs difference: " << max_abs_diff << endl;
      std::exit(1);
      }
    }

    // Test succeeded
    exit(0);
}

// Invocations
template class TestImage<double, 2>;
template class TestImage<double, 3>;
template class TestImage<double, 4>;
