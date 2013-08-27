#include "SmoothImage.h"
#include "itkDiscreteGaussianImageFilter.h"

template <class TPixel, unsigned int VDim>
void
SmoothImage<TPixel, VDim>
::operator() (RealVector &stdev)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();
  
  // Describe what we are doing
  *c->verbose << "Smoothing #" << c->m_ImageStack.size() << " with std.dev. " << stdev << endl;

  // Create a smoothing kernel and use it
  typedef itk::DiscreteGaussianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename FilterType::ArrayType variance;

  for(size_t i = 0; i < VDim; i++)
    variance[i] = stdev[i] * stdev[i];
  
  filter->SetInput(input);
  filter->SetVariance(variance);
  filter->SetUseImageSpacingOn();
  filter->Update();

  // Save the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class SmoothImage<double, 2>;
template class SmoothImage<double, 3>;

