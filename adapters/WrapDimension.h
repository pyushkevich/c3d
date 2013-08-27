#ifndef __WrapDimension_h_
#define __WrapDimension_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class WrapDimension : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  WrapDimension(Converter *c) : c(c) {}

  void operator() (const IndexType &xWrap);

private:
  Converter *c;

};

#endif

