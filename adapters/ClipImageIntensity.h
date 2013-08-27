#ifndef __ClipImageIntensity_h_
#define __ClipImageIntensity_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ClipImageIntensity : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ClipImageIntensity(Converter *c) : c(c) {}

  void operator() (double iMin, double iMax);

private:
  Converter *c;

};

#endif

