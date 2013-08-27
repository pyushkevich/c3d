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
