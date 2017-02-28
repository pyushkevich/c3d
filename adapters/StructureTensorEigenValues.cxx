/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    StructureTensorEigenValues.cxx
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

#include "StructureTensorEigenValues.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "OneDimensionalInPlaceAccumulateFilter.h"

template <class TPixel, unsigned int VDim>
class StructureTensorOuterProductFunctor
{
public:

  typedef itk::CovariantVector<TPixel, VDim> GradientType;
  typedef itk::SymmetricSecondRankTensor<TPixel, VDim> TensorType;

  TensorType operator()(const GradientType &grad)
    {
    TensorType T;
    for(int i = 0; i < VDim; i++)
      for(int j = i; j < VDim; j++)
        T(i,j) = grad[i] * grad[j];

    return T;
    }
};

template <class TPixel, unsigned int VDim>
void
StructureTensorEigenValues<TPixel, VDim>
::operator() (double sigma_grad, double sigma_window)
{
  // Get image from stack
  ImagePointer img = c->PopImage();

  // Define the filters
  typedef typename itk::NumericTraits<TPixel>::RealType RealPixelType;
  typedef itk::CovariantVector<RealPixelType, VDim> VectorPixelType;
  typedef itk::Image<VectorPixelType, VDim> VectorImageType;
  typedef itk::GradientRecursiveGaussianImageFilter<ImageType, VectorImageType> GradientFilterType;

  // Compute the gradient of the input
  typename GradientFilterType::Pointer grad = GradientFilterType::New();
  grad->SetInput(img);
  grad->SetSigma(sigma_grad);

  // Compute the outer product tensor image
  typedef itk::SymmetricSecondRankTensor<RealPixelType, VDim> TensorPixelType;
  typedef itk::Image<TensorPixelType, VDim> TensorImageType;
  typedef StructureTensorOuterProductFunctor<TPixel, VDim> OuterProduct;
  typedef itk::UnaryFunctorImageFilter<VectorImageType, TensorImageType, OuterProduct> OuterProductFilter;

  typename OuterProductFilter::Pointer op = OuterProductFilter::New();
  op->SetInput(grad->GetOutput());

  // Smooth the outer product iamge
  typedef itk::SmoothingRecursiveGaussianImageFilter<TensorImageType, TensorImageType> SmoothFilter;
  typename SmoothFilter::Pointer smooth = SmoothFilter::New();
  smooth->SetSigma(sigma_window);
  smooth->SetInput(op->GetOutput());

  // Compute the sorted eigenvalues
  typedef itk::SymmetricEigenAnalysisImageFilter<TensorImageType, VectorImageType>
    EigenValueFilterType;
  typename EigenValueFilterType::Pointer eigen = EigenValueFilterType::New();
  eigen->SetInput(smooth->GetOutput());
  eigen->SetDimension(VDim);
  eigen->OrderEigenValuesBy(EigenValueFilterType::FunctorType::OrderByValue);

  // Run the filter
  *c->verbose << "Computing Structure Tensor Eigenvalues from #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Image Gradient Sigma: " << sigma_grad << endl;
  *c->verbose << "  Window Sigma: " << sigma_window << endl;

  // Update the eigen filter
  eigen->Update();

  // Get the component images
  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType>
    ComponentFilter;

  // Put result on stack
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
template class StructureTensorEigenValues<double, 2>;
template class StructureTensorEigenValues<double, 3>;
template class StructureTensorEigenValues<double, 4>;
