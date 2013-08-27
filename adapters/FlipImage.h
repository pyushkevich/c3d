#ifndef __FlipImage_h_
#define __FlipImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class FlipImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  FlipImage(Converter *c) : c(c) {}

  void operator() (string axis);

private:
  Converter *c;

};

#endif

