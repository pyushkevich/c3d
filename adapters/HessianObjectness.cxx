/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    HessianObjectness.cxx
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

#include "HessianObjectness.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"

template <class TPixel, unsigned int VDim>
void
HessianObjectness<TPixel, VDim>
::operator() (int dimension, double minscale, double maxscale)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Define the filter
  typedef typename itk::NumericTraits<TPixel>::RealType RealPixelType;
  typedef itk::SymmetricSecondRankTensor<RealPixelType, VDim> HessianPixelType;
  typedef itk::Image<HessianPixelType, VDim> HessianImageType;

  typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, ImageType> 
    ObjectnessFilterType;

  typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType, ImageType>
    MultiScaleEnhancementFilterType;

  // Create the objectness filter
  typename ObjectnessFilterType::Pointer of = ObjectnessFilterType::New();
  of->SetScaleObjectnessMeasure(true); // why?
  of->SetBrightObject(dimension > 0);
  of->SetObjectDimension(abs(dimension));
  of->SetAlpha(0.5);
  of->SetBeta(0.5);
  of->SetGamma(5.0);

  // Create the main filter
  typename MultiScaleEnhancementFilterType::Pointer filter = MultiScaleEnhancementFilterType::New();
  filter->SetInput(img);
  filter->SetHessianToMeasureFilter(of);
  filter->SetSigmaStepMethodToLogarithmic();
  filter->SetSigmaMaximum(maxscale);
  filter->SetSigmaMinimum(minscale);
  filter->SetNumberOfSigmaSteps(minscale == maxscale ? 1 : 10);
  
  // Run the filter
  *c->verbose << "Extracting Hessian Objectness from #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Object dimension: " << of->GetObjectDimension() << endl; 
  *c->verbose << "  Object type: " << (of->GetBrightObject() ? "bright" : "dark") << endl; 
  *c->verbose << "  Sigma range: " << filter->GetSigmaMinimum() << " " << filter->GetSigmaMaximum() << endl;

  filter->Update();

  // Do some processing ...
  ImagePointer result = filter->GetOutput();
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class HessianObjectness<double, 2>;
template class HessianObjectness<double, 3>;
template class HessianObjectness<double, 4>;
