#include "ImageERF.h"
#include "vnl/vnl_erf.h"

template <class TPixel, unsigned int VDim>
void
ImageERF<TPixel, VDim>
::operator() (double thresh, double scale)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Use iterator to apply erf
  itk::ImageRegionIteratorWithIndex<ImageType> 
    it(input, input->GetBufferedRegion());
  for(; !it.IsAtEnd(); ++it)
    {
    double x = it.Value();
    double y = vnl_erf((x - thresh) / scale);
    it.Set(y);
    }

  // Say what we are doing
  *c->verbose << "Taking ERF of #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  y = erf((x - " << thresh << ") / scale)" << endl;

  // Updated
  input->Modified();
}

// Invocations
template class ImageERF<double, 2>;
template class ImageERF<double, 3>;
