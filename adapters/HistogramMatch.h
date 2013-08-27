#ifndef __HistogramMatch_h_
#define __HistogramMatch_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class HistogramMatch : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  HistogramMatch(Converter *c) : c(c) {}

  void operator()(int nmatch);

private:
  Converter *c;

};

#endif

