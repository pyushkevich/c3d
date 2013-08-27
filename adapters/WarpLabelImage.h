#ifndef __WarpLabelImage_h_
#define __WarpLabelImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class WarpLabelImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  WarpLabelImage(Converter *c) : c(c) {}

  void operator()(RealVector &stdev);

private:
  Converter *c;

};

#endif

