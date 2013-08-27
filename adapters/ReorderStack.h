#ifndef __ReorderStack_h_
#define __ReorderStack_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ReorderStack : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ReorderStack(Converter *c) : c(c) {}

  void operator() (size_t n);

private:
  Converter *c;

};

#endif

