#ifndef __StapleAlgorithm_h_
#define __StapleAlgorithm_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class StapleAlgorithm : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  StapleAlgorithm(Converter *c) : c(c) {}

  void operator() (double ival);

private:
  Converter *c;

};

#endif

