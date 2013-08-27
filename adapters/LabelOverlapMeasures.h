#ifndef __LabelOverlapMeasures_h_
#define __LabelOverlapMeasures_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class LabelOverlapMeasures : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  LabelOverlapMeasures(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

