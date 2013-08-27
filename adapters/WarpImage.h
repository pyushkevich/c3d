#ifndef __WarpImage_h_
#define __WarpImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class WarpImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  WarpImage(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

