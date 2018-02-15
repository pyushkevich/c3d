/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ApplyNoise.cxx
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

#include "ApplyNoise.h"
#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "itkSaltAndPepperNoiseImageFilter.h"
#include "itkShotNoiseImageFilter.h"
#include "itkSpeckleNoiseImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ApplyNoise<TPixel, VDim>
::operator() (Mode mode, double param)
{
  // Get image from stack
  ImagePointer img = c->PopImage();
  ImagePointer result;

  // Create the appropriate filter
  if(mode == GAUSSIAN)
    {
    typedef itk::AdditiveGaussianNoiseImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(img);
    filter->SetStandardDeviation(param);
    filter->Update();
    result = filter->GetOutput();
    }
  else if(mode == POISSON)
    {
    typedef itk::ShotNoiseImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(img);
    filter->SetScale(param);
    filter->Update();
    result = filter->GetOutput();
    }
  else if(mode == SALT_PEPPER)
    {
    typedef itk::SaltAndPepperNoiseImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(img);
    filter->SetProbability(param);
    filter->Update();
    result = filter->GetOutput();
    }
  else if(mode == SPECKLE)
    {
    typedef itk::SpeckleNoiseImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(img);
    filter->SetStandardDeviation(param);
    filter->Update();
    result = filter->GetOutput();
    }
  else throw ConvertException("Unknown noise type");

  c->PushImage(result);
}

// Invocations
template class ApplyNoise<double, 2>;
template class ApplyNoise<double, 3>;
template class ApplyNoise<double, 4>;
