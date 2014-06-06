/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    CreateInterpolator.cxx
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

#include "CreateInterpolator.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"

template <class TPixel, unsigned int VDim>
void
CreateInterpolator<TPixel, VDim>
::CreateNN()
{
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double> NNInterpolatorType;
  c->SetInterpolator(NNInterpolatorType::New());
}

template <class TPixel, unsigned int VDim>
void
CreateInterpolator<TPixel, VDim>
::CreateLinear()
{
  typedef itk::LinearInterpolateImageFunction<ImageType,double> LinearInterpolatorType;
  c->SetInterpolator(LinearInterpolatorType::New());
}

template <class TPixel, unsigned int VDim>
void
CreateInterpolator<TPixel, VDim>
::CreateCubic()
{
  typedef itk::BSplineInterpolateImageFunction<ImageType,double> CubicInterpolatorType;
  c->SetInterpolator(CubicInterpolatorType::New());
}

template <class TPixel, unsigned int VDim>
void
CreateInterpolator<TPixel, VDim>
::CreateSinc()
{
  typedef itk::WindowedSincInterpolateImageFunction<ImageType, 4> SincInterpolatorType;
  c->SetInterpolator(SincInterpolatorType::New());
}

template <class TPixel, unsigned int VDim>
void
CreateInterpolator<TPixel, VDim>
::CreateGaussian(RealVector sigma)
{
  typedef itk::GaussianInterpolateImageFunction<ImageType, double> GaussianInterpolatorType;
  typename GaussianInterpolatorType::Pointer gi = GaussianInterpolatorType::New();
  gi->SetParameters(sigma.data_block(), 4.0);
  c->SetInterpolator(gi);
}

template <class TPixel, unsigned int VDim>
void
CreateInterpolator<TPixel, VDim>
::CreateMultiLabel(RealVector sigma)
{
  typedef itk::LabelImageGaussianInterpolateImageFunction<ImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer gi = InterpolatorType::New();
  gi->SetParameters(sigma.data_block(), 4.0);
  c->SetInterpolator(gi);
}

// Invocations
template class CreateInterpolator<double, 2>;
template class CreateInterpolator<double, 3>;
template class CreateInterpolator<double, 4>;
