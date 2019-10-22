/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    HessianEigenValues.cxx
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

#include "HessianEigenValues.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

template <class TPixel, unsigned int VDim>
void
HessianEigenValues<TPixel, VDim>
::operator() (double scale)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Define the filters
  typedef typename itk::NumericTraits<TPixel>::RealType RealPixelType;
  typedef itk::SymmetricSecondRankTensor<RealPixelType, VDim> HessianPixelType;
  typedef itk::Image<HessianPixelType, VDim> HessianImageType;
  typedef itk::Image<itk::CovariantVector<RealPixelType, VDim>, VDim > EigenValueImageType;

  // Filter to compute Hessian
  typedef itk::HessianRecursiveGaussianImageFilter<ImageType, HessianImageType> 
    HessianFilterType;

  typename HessianFilterType::Pointer hf = HessianFilterType::New();
  hf->SetInput(img);
  hf->SetSigma(scale);

  // Create the eigenvalue computer
  typedef itk::SymmetricEigenAnalysisImageFilter<HessianImageType, EigenValueImageType>
    EigenValueFilterType;
  typename EigenValueFilterType::Pointer eigen = EigenValueFilterType::New();
  eigen->SetInput(hf->GetOutput());
  eigen->SetDimension(VDim);
  eigen->OrderEigenValuesBy(EigenValueFilterType::FunctorType::OrderByValue);

  // Run the filter
  *c->verbose << "Computing Hessian Eigenvalues from #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Sigma: " << scale << endl;

  // Update the eigen filter
  eigen->Update();

  // Get the component images
  typedef itk::VectorIndexSelectionCastImageFilter<EigenValueImageType, ImageType>
    ComponentFilter;

  // Put result on stack
  c->m_ImageStack.pop_back();

  for(int i = 0; i < VDim; i++)
    {
    typename ComponentFilter::Pointer comp = ComponentFilter::New();
    comp->SetInput(eigen->GetOutput());
    comp->SetIndex(i);
    comp->Update();

    c->m_ImageStack.push_back(comp->GetOutput());
    }
}

// Invocations
template class HessianEigenValues<double, 2>;
template class HessianEigenValues<double, 3>;
template class HessianEigenValues<double, 4>;
