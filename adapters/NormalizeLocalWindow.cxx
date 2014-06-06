/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    NormalizeLocalWindow.cxx
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

#include "NormalizeLocalWindow.h"
#include "itkImageKernelOperator.h"
#include "itkNeighborhoodOperatorImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkSqrtImageFilter.h"

template <class TPixel, unsigned int VDim>
void
NormalizeLocalWindow<TPixel, VDim>
::operator() (SizeType radius)
{
  // This filter replaces the intensity g(x) at each voxel by (g-mu)/sigma, where mu and sigma
  // are the mean and standard deviation of intensity in a neighborhood. To compute the filter,
  // we need to compute these sliding window statistics
  
  // Create the mean image kernel
  ImagePointer kernel = ImageType::New();
  RegionType kregion;
  for(int i = 0; i < VDim; i++)
    {
    kregion.SetIndex(i, 0);
    kregion.SetSize(i, radius[i]*2+1);
    }
  kernel->SetRegions(kregion);
  kernel->Allocate();
  kernel->FillBuffer(1.0 / kregion.GetNumberOfPixels());

  // Apply kernel
  typedef itk::ImageKernelOperator<TPixel, VDim> OperatorType;
  OperatorType op;
  op.SetImageKernel(kernel);
  op.CreateToRadius(radius);

  // Create a mean filter
  ImagePointer img = c->m_ImageStack.back();
  typedef itk::NeighborhoodOperatorImageFilter<ImageType, ImageType> MeanFilter;
  typename MeanFilter::Pointer fmean = MeanFilter::New();
  fmean->SetInput(img);
  fmean->SetOperator(op);

  // Square the input image
  typedef itk::MultiplyImageFilter<ImageType, ImageType> MultFilter;
  typename MultFilter::Pointer fltMult = MultFilter::New();
  fltMult->SetInput1(img);
  fltMult->SetInput2(img);

  // Create the mean filter for the squared image
  typename MeanFilter::Pointer fmeansq = MeanFilter::New();
  fmeansq->SetInput(fltMult->GetOutput());
  fmeansq->SetOperator(op);

  // Compute the standard deviation
  typedef itk::SubtractImageFilter<ImageType, ImageType> SubtractFilter;
  typedef itk::SqrtImageFilter<ImageType, ImageType> SqrtFilter;
  typename SubtractFilter::Pointer fltSub = SubtractFilter::New();
  typename SqrtFilter::Pointer fltSqrt = SqrtFilter::New();
  typename MultFilter::Pointer fltSqMean = MultFilter::New();

  // Square the mean
  fltSqMean->SetInput1(fmean->GetOutput());
  fltSqMean->SetInput2(fmean->GetOutput());

  // Subtract squared mean from sum of squares
  fltSub->SetInput1(fmeansq->GetOutput());
  fltSub->SetInput2(fltSqMean->GetOutput());

  // Take square root - for standard deviation
  fltSqrt->SetInput(fltSub->GetOutput());

  // Subtract the mean from the data
  typename SubtractFilter::Pointer fltMeanSub = SubtractFilter::New();
  fltMeanSub->SetInput1(img);
  fltMeanSub->SetInput2(fmean->GetOutput());

  // Divide by the standard deviation
  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DivideFilter;
  typename DivideFilter::Pointer fltDivide = DivideFilter::New();
  fltDivide->SetInput1(fltMeanSub->GetOutput());
  fltDivide->SetInput2(fltSqrt->GetOutput());
  fltDivide->Update();

  // Do some processing ...
  ImagePointer result = fltDivide->GetOutput();
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class NormalizeLocalWindow<double, 2>;
template class NormalizeLocalWindow<double, 3>;
template class NormalizeLocalWindow<double, 4>;
