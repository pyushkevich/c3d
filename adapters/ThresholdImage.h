#ifndef __ThresholdImage_h_
#define __ThresholdImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ThresholdImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ThresholdImage(Converter *c) : c(c) {}

  void operator() (double u1, double u2, double vIn, double vOut);

private:
  Converter *c;

};

#endif

