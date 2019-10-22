/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    MomentsFeatures.cxx
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

#include "MomentsFeatures.h"

#include "CreateImage.h"
#include "CoordinateMap.h"
#include "ScaleShiftImage.h"
#include "MultiplyImages.h"
#include "Convolution.h"

#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIterator.h"

template <class TPixel, unsigned int VDim>
void
MomentsFeatures<TPixel, VDim>
::operator() (SizeType radius)
{
  // Get image from stack
  ImagePointer img = c->PopImage();

  // Generate a blank kernel image
  CreateImage<TPixel,VDim> a_create(c); 

  SizeType k_size;
  for(int j = 0; j < VDim; j++)
    k_size[j] = radius[j] * 2 + 1;
  a_create(k_size, img->GetSpacing().GetVnlVector());

  // Generate normalized coordinate maps
  CoordinateMap<TPixel,VDim> a_coord(c);
  a_coord(false);
  ImagePointer k_xyz[VDim];

  for(int j = VDim-1; j>=0; j--)
    {
    ScaleShiftImage<TPixel,VDim> a_scale(c);
    a_scale(1.0 / radius[j], -0.5);
    k_xyz[j] = c->PopImage();
    }

  // Allocate an image of symmetric tensors
  typedef typename itk::NumericTraits<TPixel>::RealType RealPixelType;
  typedef itk::SymmetricSecondRankTensor<RealPixelType, VDim> TensorPixelType;
  typedef itk::Image<TensorPixelType, VDim> TensorImageType;

  typename TensorImageType::Pointer tensor = TensorImageType::New();
  tensor->SetRegions(img->GetBufferedRegion());
  tensor->CopyInformation(img);
  tensor->Allocate();

  // Compute the second moments
  for(int i = 0; i < VDim; i++)
    {
    for(int j = i; j < VDim; j++)
      {
      // Multiply the i-th and j-th coordinate maps
      c->PushImage(k_xyz[i]);
      c->PushImage(k_xyz[j]);
      MultiplyImages<TPixel,VDim> a_mult(c);
      a_mult();

      // Get the kernel
      ImagePointer kernel = c->PopImage();
      c->PushImage(img);
      c->PushImage(kernel);

      // Perform convolution
      Convolution<TPixel,VDim> a_conv(c);
      a_conv();
      ImagePointer moment = c->PopImage();

      // Copy the moment into the tensor image
      typedef itk::ImageRegionIterator<TensorImageType> TensorIter;
      Iterator it(moment, moment->GetBufferedRegion());
      TensorIter itt(tensor, tensor->GetBufferedRegion());
      for(; !itt.IsAtEnd(); ++itt, ++it)
        {
        itt.Value()(j,i) = it.Value();
        }
      }
    }

  // Create the eigenvalue computer
  typedef itk::Image<itk::CovariantVector<RealPixelType, VDim>, VDim> EigenValueImageType;
  typedef itk::SymmetricEigenAnalysisImageFilter<TensorImageType, EigenValueImageType>
    EigenValueFilterType;
  typename EigenValueFilterType::Pointer eigen = EigenValueFilterType::New();
  eigen->SetInput(tensor);
  eigen->SetDimension(VDim);
  eigen->OrderEigenValuesBy(EigenValueFilterType::FunctorType::OrderByValue);
   
  // Update the eigen filter
  eigen->Update();

  // Get the component images
  typedef itk::VectorIndexSelectionCastImageFilter<EigenValueImageType, ImageType>
    ComponentFilter;

  for(int i = 0; i < VDim; i++)
    {
    typename ComponentFilter::Pointer comp = ComponentFilter::New();
    comp->SetInput(eigen->GetOutput());
    comp->SetIndex(i);
    comp->Update();

    c->PushImage(comp->GetOutput());
    }
}

// Invocations
template class MomentsFeatures<double, 2>;
template class MomentsFeatures<double, 3>;
template class MomentsFeatures<double, 4>;
