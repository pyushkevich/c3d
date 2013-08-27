#include "SignedDistanceTransform.h"
#include "ThresholdImage.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

template <class TPixel, unsigned int VDim>
void
SignedDistanceTransform<TPixel, VDim>
::operator() ()
{
  // The image is assumed to be binary. If background is non-zero, call binarize
  // to map the background to zero
  if(c->m_Background != 0.0)
    {
    ThresholdImage<TPixel, VDim> thresh(c);
    thresh(c->m_Background, c->m_Background, 0.0, 1.0);
    }

  // Get the last image on the stack
  ImagePointer image = c->m_ImageStack.back();

  // Construct the connected components filter
  typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> Filter;

  // Describe what we are doing
  *c->verbose << "Computing signed distance function of #" << c->m_ImageStack.size() << endl;

  // Plug in the filter's components
  typename Filter::Pointer flt = Filter::New();
  flt->SetInput(image);
  flt->SetUseImageSpacing(true);
  flt->SquaredDistanceOff();
  flt->Update();

  // Store the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(flt->GetOutput());

}

// Invocations
template class SignedDistanceTransform<double, 2>;
template class SignedDistanceTransform<double, 3>;
