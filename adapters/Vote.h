#ifndef __Vote_h_
#define __Vote_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class Vote : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  Vote(Converter *c) : c(c) {}

  void operator() (bool flagUseSplitLabelSet);

private:
  Converter *c;

};

#endif

