#include "ImageLaplacian.h"
#include "itkLaplacianImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ImageLaplacian<TPixel, VDim>
::operator() ()
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();
  
  // Describe what we are doing
  *c->verbose << "Taking Laplacian of #" << c->m_ImageStack.size() << endl;

  // Create a smoothing kernel and use it
  typedef itk::LaplacianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(input);
  filter->UseImageSpacingOn();
  filter->Update();

  // Save the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class ImageLaplacian<double, 2>;
template class ImageLaplacian<double, 3>;
