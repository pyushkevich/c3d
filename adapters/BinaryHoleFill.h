#ifndef __BinaryHoleFill_h_
#define __BinaryHoleFill_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class BinaryHoleFill : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  BinaryHoleFill(Converter *c) : c(c) {}

  void operator() (double foreground, bool full_conn);

private:
  Converter *c;

};

#endif

