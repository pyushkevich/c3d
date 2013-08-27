#ifndef __PadImage_h_
#define __PadImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class PadImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  PadImage(Converter *c) : c(c) {}

  void operator() (IndexType padExtentLower, IndexType padExtentUpper, float padValue);

private:
  Converter *c;

};

#endif

