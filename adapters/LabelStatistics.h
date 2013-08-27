#ifndef __LabelStatistics_h_
#define __LabelStatistics_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class LabelStatistics : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  LabelStatistics(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

