#include "ClipImageIntensity.h"

template <class TPixel, unsigned int VDim>
void
ClipImageIntensity<TPixel, VDim>
::operator() (double iMin, double iMax)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  *c->verbose << "Clipping out-of-range intensities in #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Intensity range: " << iMin << " to " << iMax << endl;

  // Simply replace values outside the clip range with new values
  for(Iterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(it.Value() < iMin)
      it.Set(iMin);
    else if(it.Value() > iMax)
      it.Set(iMax);
    }

  // Update the image
  img->Modified();
}

// Invocations
template class ClipImageIntensity<double, 2>;
template class ClipImageIntensity<double, 3>;
