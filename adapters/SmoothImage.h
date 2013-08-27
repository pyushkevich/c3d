#ifndef __SmoothImage_h_
#define __SmoothImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class SmoothImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  SmoothImage(Converter *c) : c(c) {}

  void operator() (RealVector &stdev);

private:
  Converter *c;

};

#endif

