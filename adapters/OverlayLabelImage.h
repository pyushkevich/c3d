#ifndef __OverlayLabelImage_h_
#define __OverlayLabelImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class OverlayLabelImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  OverlayLabelImage(Converter *c) : c(c) {}

  void operator() (const char *fnLUT, double opacity);

private:
  Converter *c;

};

#endif

