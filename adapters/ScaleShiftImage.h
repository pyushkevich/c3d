#ifndef __ScaleShiftImage_h_
#define __ScaleShiftImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ScaleShiftImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ScaleShiftImage(Converter *c) : c(c) {}

  void operator() (double a, double b);

private:
  Converter *c;

};

#endif

