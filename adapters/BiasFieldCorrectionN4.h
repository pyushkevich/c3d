#ifndef __BiasFieldCorrectionN4_h_
#define __BiasFieldCorrectionN4_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class BiasFieldCorrectionN4 : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  BiasFieldCorrectionN4(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

