#ifndef __ComputeOverlaps_h_
#define __ComputeOverlaps_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ComputeOverlaps : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ComputeOverlaps(Converter *c) : c(c) {}

  void operator() (double v);

private:
  Converter *c;

};

#endif

