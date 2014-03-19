#ifndef __HessianObjectness_h_
#define __HessianObjectness_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class HessianObjectness : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  HessianObjectness(Converter *c) : c(c) {}

  void operator() (int codimension, double minscale, double maxscale);

private:
  Converter *c;

};

#endif

