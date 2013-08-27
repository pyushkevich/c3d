#ifndef __CoordinateMap_h_
#define __CoordinateMap_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class CoordinateMap : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  CoordinateMap(Converter *c) : c(c) {}

  void operator() (bool physical);

private:
  Converter *c;

};

#endif

