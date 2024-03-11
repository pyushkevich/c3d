/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    TustisonWellComposedness.cxx
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

#include "TustisonWellComposedness.h"
#include "itkTopologyPreservingDigitalSurfaceEvolutionImageFilter.h"
#include "itkBinaryDiamondStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"

template <class TPixel, unsigned int VDim>
void
TustisonWellComposedness<TPixel, VDim>
::operator() (int label)
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    throw ConvertException("Binary operations require two images on the stack");

  // Get the last two images
  ImagePointer itrg = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer isrc = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Extract the label from both images
  typedef itk::Image<unsigned char, VDim> LabelImageType;
  typedef itk::BinaryThresholdImageFilter<ImageType,LabelImageType> ThreshType;
  typename ThreshType::Pointer tsrc = ThreshType::New();
  tsrc->SetInput(itrg);
  tsrc->SetLowerThreshold(label);
  tsrc->SetUpperThreshold(label);
  tsrc->SetInsideValue(1);
  tsrc->SetOutsideValue(0);
  tsrc->Update();

  typename ThreshType::Pointer ttrg = ThreshType::New();
  ttrg->SetInput(isrc);
  ttrg->SetLowerThreshold(label);
  ttrg->SetUpperThreshold(label);
  ttrg->SetInsideValue(1);
  ttrg->SetOutsideValue(0);
  ttrg->Update();

  // Create the evolution filter
  typedef itk::TopologyPreservingDigitalSurfaceEvolutionImageFilter<LabelImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(tsrc->GetOutput());
  filter->SetTargetImage(ttrg->GetOutput());
  filter->SetNumberOfIterations(1000);
  filter->SetForegroundValue(1);
  filter->SetUseInversionMode(false);

  typedef itk::CastImageFilter<LabelImageType, ImageType> Cast;
  typename Cast::Pointer caster = Cast::New();
  caster->SetInput(filter->GetOutput());
  caster->Update();

  ImagePointer result = caster->GetOutput();

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class TustisonWellComposedness<double, 2>;
template class TustisonWellComposedness<double, 3>;
template class TustisonWellComposedness<double, 4>;
