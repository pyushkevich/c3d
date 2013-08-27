#ifndef __ExtractRegion_h_
#define __ExtractRegion_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ExtractRegion : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ExtractRegion(Converter *c) : c(c) {}

  void operator() (RegionType bbox);

private:
  Converter *c;

};

#endif

