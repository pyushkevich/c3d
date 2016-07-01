/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    MixtureModel.cxx
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

#include "MixtureModel.h"

#include "itkImageToListSampleAdaptor.h"
#include "itkGaussianMixtureModelComponent.h"
#include "itkExpectationMaximizationMixtureModelEstimator.h"
#include "itkComposeImageFilter.h"

using namespace itk;
using namespace itk::Statistics;
using namespace std;

template <class TPixel, unsigned int VDim>
void
MixtureModel<TPixel, VDim>
::operator() (std::vector<double> mu, std::vector<double> sigma)
{
  // Get image from stack and create a statistics sample
  ImagePointer img = c->m_ImageStack.back();

  // creat a new image with array pixel type from the source
  typedef itk::FixedArray< TPixel, 1 > ArrayPixelType ;
  typedef itk::Image< ArrayPixelType, VDim > ArrayPixelImageType ;
  typedef itk::ComposeImageFilter< ImageType, ArrayPixelImageType > ImageCastFilterType ;
  typename ImageCastFilterType::Pointer castFilter = ImageCastFilterType::New() ;
  castFilter->SetInput(img);
  castFilter->Update() ;
  
  typedef ImageToListSampleAdaptor<ArrayPixelImageType> DataSampleType;
  typename DataSampleType::Pointer sample = DataSampleType::New();
  sample->SetImage(castFilter->GetOutput());

  // Component definition
  typedef GaussianMixtureModelComponent<DataSampleType> ComponentType;
  typedef typename ComponentType::Pointer ComponentPointer;
  vector<ComponentPointer> comps;

  // Create components and initial proportions
  itk::Array<double> initprop(mu.size());
  for(size_t i = 0; i < mu.size(); i++)
    {
    ComponentPointer cmp = ComponentType::New();
    cmp->SetSample(sample);
   
    itk::Array<double> param(2);
    param[0] = mu[i];
    param[1] = sigma[i] * sigma[i];
    cmp->SetParameters(param);
    comps.push_back(cmp);

    initprop[i] = 1.0 / mu.size();
    }

  // Explain what we are doing
  *c->verbose << "Running Gaussian Mixture Modeling on #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Initial Parameters : " << endl;
  for(size_t i = 0; i < mu.size(); i++)
    {
    *c->verbose << "    " << 
      "Class " << i << ": " <<
      "mu = " << mu[i] << "; " <<
      "sigma = " << sigma[i] << "; " <<
      "alpha = " << initprop[i] << "; " << endl;

    }

  typedef ExpectationMaximizationMixtureModelEstimator<DataSampleType> EstimatorType;
  typename EstimatorType::Pointer estimator = EstimatorType::New();

  estimator->SetSample(sample);
  estimator->SetMaximumIteration(100);
  estimator->SetInitialProportions(initprop);
  for(size_t i = 0; i < mu.size(); i++)
    estimator->AddComponent(comps[i]);


  estimator->Update();

  *c->verbose << "  Estimated Parameters : " << endl;
  for(size_t i = 0; i < mu.size(); i++)
    {
    *c->verbose << "    " << 
      "Class " << i << ": " <<
      "mu = " << comps[i]->GetFullParameters()[0] << "; " <<
      "sigma = " << sqrt(comps[i]->GetFullParameters()[1]) << "; " <<
      "alpha = " << (estimator->GetProportions())[i] << "; " << endl;
    }

  // Do some processing ...
  // ImagePointer result = ...;
  
  // Put result on stack
  // c->m_ImageStack.pop_back();
  // c->m_ImageStack.push_back(result);
}

// Invocations
template class MixtureModel<double, 2>;
template class MixtureModel<double, 3>;
template class MixtureModel<double, 4>;
